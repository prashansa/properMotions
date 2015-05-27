import numpy as np
from astropy.time import Time
from esutil.coords import sphdist
import healpy as hp
# def calcResiduals(rBar, r):
#     residual = []*len(r)
#     #rBar is one value, r is in general an array
#     for i in r:
#         residual.append(rBar - i)
#     return residual

Noofyears = 4 
pixelWidth = 3 #arcminutes
searchRadius = 10 #arcminutes
#load data into a variable
testFile= np.load('/home/gupta/projects/propermotion/test.npy')

#save unique galaxy ids in another variable
uniqueObjIDfile = np.unique(testFile['obj_id'])

i=0
j=0
k=0

medianArray=[]*len(uniqueObjIDfile)
residualRaArray= []
residualDecArray= []


for row in testFile:
    maskForID = row['obj_id'] == uniqueObjIDfile
    raOfObj = testFile['ra'][maskForID]
    decOfObj = testFile['dec'][maskForID]
    medianRa = np.median(raOfObj)
    medianDec = np.median(decOfObj)
    medianArray.append([medianRa, medianDec])
    residualRa = raOfObj - medianRa
    residualDec = decOfObj - medianDec
    residualRaArray.append(residualRa)
    residualDecArray.append(residualDec)



mjdSorted=np.sort(testFile['mjd'])

deltaT=mjdSorted[0:-1]-mjdSorted[1:]

#where can you 
tBreakAt=np.where(condition)

#inputArray will be mjds and the binArray will be mjd plus minus deltaT
indexArray = np.digitize(inputArray, binArray)



#pixelisation

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

