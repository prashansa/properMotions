import numpy as np
from astropy.time import Time
from time import time
from esutil.coords import sphdist
import healpy as hp

# def calcResiduals(rBar, r):
#     residual = []*len(r)
#     #rBar is one value, r is in general an array
#     for i in r:
#         residual.append(rBar - i)
#     return residual

noOfYears = 5
pixelWidth = 3 #arcminutes
searchRadius = 10 #arcminutes

#load data into a variable
testFile= np.load('/home/gupta/projects/propermotion/test.npy')

#save unique galaxy ids in another variable
uniqueObjIDfile = np.unique(testFile['obj_id'])

i=0
j=0
k=0

medianArray=[0]*len(uniqueObjIDfile)
residualRaArray= [0]*len(uniqueObjIDfile)
residualDecArray= [0]len(uniqueObjIDfile)


t0 = time()

for row in uniqueObjIDfile:
    maskForID = testFile['obj_id'] == row
    raOfObj = testFile['ra'][maskForID]
    decOfObj = testFile['dec'][maskForID]
    medianRa = np.median(raOfObj)
    medianDec = np.median(decOfObj)
    medianArray[i] = medianRa, medianDec
    residualRa = raOfObj - medianRa
    residualDec = decOfObj - medianDec
    residualRaArray[i] = residualRa
    residualDecArray[i] = residualDec

print 'time taken by for loop',time()-t0
    
#sorts mjd values in ascending order
mjdSorted = np.sort(testFile['mjd'])

#takes differences of neighbors
deltaT = mjdSorted[0:-1] - mjdSorted[1:]

#tBreak = deltaT < (-100.0)
#tBreakAt=np.where(tBreak)

#returns indices in deltaT array where the condition holds, [0]is added to make the tuple an array.
tBreakAt = np.where(deltaT < (-100.0))[0]

#find the mjd value corresponding to tBreak value
mjdBreakAt = mjdSorted[tBreakAt+1] 

##inputArray will be mjds and the binArray will be mjd plus minus deltaT
#indexArray = np.digitize(inputArray, binArray)

#define epochs according to the mjd breakpoints

#epoch = np.zeros (len(mjdSorted))
epoch = [0]* noOfYears
#epoch = np.arange(noOfYears)
previousVar2 = np.min(mjdSorted)

for var1, var2 in enumerate(mjdBreakAt, start=0):
    mjdIndex = (mjdSorted < var2) & (previousVar2 < mjdSorted)
    #define your epochs 
    epoch[var1] = mjdSorted[mjdIndex]
    #print previousVar2, var2
    previousVar2 = var2

    

#residualRAyear = np.arange(noOfYears)
#for i in range(noOfYears):
 #   residualRAyear[i]= 




for i in range(noOfYears):
    for row in testFile:
        mask1 = epoch[i] == row['mjd']
        #print mask1
        epochObjID = testFile['obj_id'][mask]
        mask2 = epochObjID == uniqueObjIDfile
        raMedianForEpoch = np.median(residualRaArray[mask2])
        decMedianForEpoch = np.median(residualDecArray[mask2])



#pixelisation
#resolution
NSIDE = 2**10

#nside2npix give the number of pixels for a given NSIDE
m = np.arange(hp.nside2npix(NSIDE))

#pix2ang gives the angular coordinates
theta, phi = hp.pix2ang(NSIDE, m)

#pixelcenters
raCenter = 180*phi/np.pi
decCenter = 90 - theta*180/np.pi
    

#raCenter and decCenter are pixel centers
angularSeparation=sphdist(raCenter, decCenter, testFile['ra'], testFile['dec'])

temp = angularSeparation*60 < searchRadius

searchThisArrayIDs = np.unique(testFile['obj_id'][temp])

