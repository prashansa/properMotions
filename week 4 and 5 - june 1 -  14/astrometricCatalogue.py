import numpy as np
from astropy.time import Time
from time import time
from esutil.coords import sphdist
import healpy as hp

#load data into a variable
testFile= np.sort(np.load('/home/gupta/projects/propermotion/test.npy'))

#save unique galaxy ids in another variable
uniqueObjIDfileTemp = np.unique(testFile['obj_id'])

#sorting the unique object id file by object id -- also keeps all the objects with same object ids together! 
uniqueObjIDfile = np.sort(uniqueObjIDfileTemp)

#all the variables
#the maximum number of years for which the data is available
noOfYears = 5
#the width of the pixel -- if user defined
pixelWidth = 3 #in arcminutes
#the distance upto which the galaxies are included for updating ra/dec
searchRadius = 2#in arcminutes
#the number of observations actually available for each object
#noOfObs = np.zeros(len(uniqueObjIDfile)) 

#stores the median ra and dec values for each object
#medianRaArray=np.zeros(uniqueObjIDfile.size)
#medianDecArray=np.zeros(uniqueObjIDfile.size)

#stores the residual values for each object and for each epoch
#residualRaArray= np.zeros(len(testFile))
#residualDecArray= np.zeros(len(testFile))

#store the first object values in here
#currentObjID = testFile['obj_id'][0]
#raObj = [testFile['ra'][0]]
#decObj = [testFile['dec'][0]]
#mjdObj = [testFile['mjd'][0]]


#this code pixelizes the dataset
nside= 2**10

#nside2npix gives number of pixels for a given nside, create an array containing pixel indices  
pixelIndexArray= np.arange(hp.nside2npix(nside))
#obtain angular coordinates corresponding to nside and the pixel indices
theta, phi = hp.pix2ang(nside, pixelIndexArray)

#pixel centers in degrees
pixelRa = 180*phi/np.pi
pixelDec = 90 - theta*180/np.pi


#pick a pixel
pickPixelNo = pixelIndexArray[len(pixelIndexArray)/2]


#find objects within the searchRadius and whose noOfObs is atleast 3

#--------------------
#VERSION 1
#--------------------

angSepMask1 =  sphdist( pixelRa[pickPixelNo], pixelDec[pickPixelNo], testFile['ra'], testFile['dec'] ) <= (searchRadius*60)

searchFile = testFile[angSepMask1]

noOfObs = np.zeros(len(np.unique(searchFile))) 

medianRaArray=np.ones((np.unique(searchFile)).size)
medianDecArray=np.ones((np.unique(searchFile)).size)

#stores the residual values for each object and for each epoch
residualRaArray= np.ones(len(searchFile))
residualDecArray= np.ones(len(searchFile))

currentObjID = searchFile['obj_id'][0]
raObj = [searchFile['ra'][0]]
decObj = [searchFile['dec'][0]]
mjdObj = [searchFile['mjd'][0]]


t0 = time()
objIDcounter = 0
pos=0

mjdSorted = np.sort(searchFile['mjd'])
deltaT = mjdSorted[0:-1] - mjdSorted[1:]
tBreakAt = np.where(deltaT < (-100.0))[0]
mjdBreakAt = mjdSorted[tBreakAt+1]

    #it should start from the second row, so add one.
for i in  np.arange(searchFile.size-1)+1:#np.arange(100)+1:#
    print i
    objID = searchFile['obj_id'][i]
    ra = searchFile['ra'][i]
    dec = searchFile['dec'][i]
    mjd = searchFile['mjd'][i]
    print ra
    if (objID == currentObjID):
        mjdObj.append(mjd)
        raObj.append(ra)
        decObj.append(dec)
        print 'same object now'
        noOfObs[objIDcounter]+=1
    elif (noOfObs[objIDcounter] >= 3.0):
        #convert list into array
        raObj=np.array(raObj)
        decObj=np.array(decObj)
        #calculate and store medians
        medianRa = np.median(raObj)
        medianDec = np.median(decObj)
        print 'medianRa is', medianRa
        print 'medianDec is', medianDec
        #median value is for each object
        medianRaArray[objIDcounter] = medianRa
        medianDecArray[objIDcounter] = medianDec
        print 'we are on the', objIDcounter,'object'
        print 'no of obs', noOfObs[objIDcounter]
        #calculate residual
        residualRa = np.array(raObj) - medianRa
        residualDec = np.array(decObj) - medianDec
        #store in separate arrays
        #print 'residualRa is', residualRa
        #print 'residualDec is', residualDec
        #print 'pos is', pos
        #print 'i is',i
        residualRaArray[pos:i] = residualRa
        residualDecArray[pos:i] = residualDec
        objIDcounter +=1
        currentObjID = objID
        pos = i
        raObj = [searchFile['ra'][i]]
        decObj = [searchFile['dec'][i]]
        mjdObj = [searchFile['mjd'][i]]
        print 'different object now'
print time()-t0

#residualRaEpochWise = np.zeros( (len(searchFile),len(mjdBreakAt+1)) )
#residualDecEpochWise = np.zeros( (len(searchFile),len(mjdBreakAt+1)) )
medianResidualRaEpochWise = np.zeros(len(mjdBreakAt+1) )
medianResidualDecEpochWise = np.zeros(len(mjdBreakAt+1) )

#if there are atleast three epochs present
if (len(mjdBreakAt)>=3.0):
    previousVar2 = np.min(mjdSorted)
    for var1, var2 in enumerate(mjdBreakAt, start=0):
        mjdIndex = (searchFile['mjd']<= var2) & (previousVar2 <= searchFile['mjd'])
        print var1
        if any(mjdIndex):
            print residualRaArray[mjdIndex]
            #residualRaEpochWise[:, var1] = residualRaArray[mjdIndex]
            medianResidualRaEpochWise[var1] = np.median(residualRaArray[mjdIndex])
            medianResidualDecEpochWise[var1] = np.median(residualDecArray[mjdIndex])
        else:
            medianResidualRaEpochWise[var1+1] = np.median(residualRaArray[mjdIndex])
            medianResidualDecEpochWise[var1+1] = np.median(residualDecArray[mjdIndex])
        previousVar2 = var2



#JUNK 
#-----------------------------------------------------------------------------------------------------

objIDin1 =np.unique(testFile['obj_id'][angSepMask1])


arr= [1,2,2,3,3,3,4,4,4,4,5,5,5,5,5,6,6,6,6,6,6,7,7,7,7,7,7,7,8,8,8,8,8,8,8,8]
np.unique(ndarray, return_index = True, return_inverse = True, return_counts = True)



VERSION 2
angSepMask2 =  (sphdist( pixelRa[pickPixelNo], pixelDec[pickPixelNo], medianRaArray, medianDecArray ) <= (searchRadius*60))  &  (noOfObs>=3.0)
objIDin2 = uniqueObjIDfile[angSepMask2]

for var1, var2 in enumerate()



m1,m2 = match(objIDin2, testFile['obj_id'])
objIDin = objIDin2[m1]
testFileIn = testFile[m2]


residualRaArray[]


#angSepMask =  sphdist( pixelRa[pickPixelNo], pixelDec[pickPixelNo],pixelRa, pixelDec ) <= (searchRadius*60)
#raIN = pixelRa[angSepMask]
#decIN = pixelDec[angSepMask]




#take their residuals

#group by epoch
