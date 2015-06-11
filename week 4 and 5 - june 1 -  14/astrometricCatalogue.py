import numpy as np
from astropy.time import Time
from time import time
from esutil.coords import sphdist
import healpy as hp

#load data into a variable
t0 = time()
testFile= np.sort(np.load('/home/gupta/projects/propermotion/test.npy'))
print time() - t0

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
searchRadius = 0.28#in arcminutes
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
nside= 2**5

#nside2npix gives number of pixels for a given nside, create an array containing pixel indices  
pixelIndexArray= np.arange(hp.nside2npix(nside))
#obtain angular coordinates corresponding to nside and the pixel indices
theta, phi = hp.pix2ang(nside, pixelIndexArray)

#pixel centers in degrees
pixelRa = 180*phi/np.pi
pixelDec = 90 - theta*180/np.pi

#calculate phi and theta for ra/dec values of the testFile
phiForObj = (testFile['ra'] * np.pi)/180 
thetaForObj = (90 - testFile['dec'])* (np.pi/180)
pixelIndexForObj= hp.ang2pix(nside,thetaForObj, phiForObj)


#pick a pixel from the testFile
pickPixelNo = pixelIndexForObj[len(pixelIndexForObj)/2]

#find objects within the searchRadius and whose noOfObs is atleast 3

#--------------------
#VERSION 1
#--------------------

#take those objects whose atleast one value of ra/dec is inside the searchRadius
angSepMask1 =  sphdist( pixelRa[pickPixelNo], pixelDec[pickPixelNo], testFile['ra'], testFile['dec'] ) <= (searchRadius*60)

#separate these entire rows from the original testFile
searchFile = testFile[angSepMask1]

#putting unique object IDs in this file
uniqueSearchFile = np.unique(searchFile['obj_id'])

#file containing data for objects inside pixel
objInPixelFile = testFile[pixelIndexForObj == pickPixelNo]

#file containing unique objIDs of the pixel
uniqueObjInPixelFile = np.unique(objInPixelFile['obj_id'])

#for containing no of observations for each object 
noOfObs = np.ones(len(uniqueSearchFile))

#for containing median ra/dec values for each object
medianRaArray=np.zeros(uniqueSearchFile.size)
medianDecArray=np.zeros(uniqueSearchFile.size)

#stores the residual values for each object and for each epoch
residualRaArray= np.ones(len(searchFile))
residualDecArray= np.ones(len(searchFile))

mjdSorted = np.sort(searchFile['mjd'])
deltaT = mjdSorted[0:-1] - mjdSorted[1:]
tBreakAt = np.where(deltaT < (-100.0))[0]
mjdBreakAt = mjdSorted[tBreakAt+1]

currentObjID = searchFile['obj_id'][0]
raObj = [searchFile['ra'][0]]
decObj = [searchFile['dec'][0]]
mjdObj = [searchFile['mjd'][0]]

t1 = time()
objIDcounter = 0
pos=0

#it should start from the second row, so add one.
for i in  np.arange(searchFile.size-1)+1:#np.arange(10)+1:#
    #print i
    objID = searchFile['obj_id'][i]
    ra = searchFile['ra'][i]
    dec = searchFile['dec'][i]
    mjd = searchFile['mjd'][i]
    #print objID
    #print ra
    #print dec
    #print mjd
    if (objID == currentObjID):
        mjdObj.append(mjd)
        raObj.append(ra)
        decObj.append(dec)
        #print 'same object now'
        noOfObs[objIDcounter]+=1
    else:
        #print 'different object now'
        #print 'we are on the', objIDcounter,'object'
        if(noOfObs[objIDcounter] >= 3.0):
            #print 'no of obs >= 3.0', noOfObs[objIDcounter]
            #convert list into array
            raObj=np.array(raObj)
            decObj=np.array(decObj)
            #calculate and store medians
            medianRa = np.median(raObj)
            medianDec = np.median(decObj)
            #print 'medianRa is', medianRa
            #print 'medianDec is', medianDec
            #median value is for each object
            medianRaArray[objIDcounter] = medianRa
            medianDecArray[objIDcounter] = medianDec
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
print 'final time is',time() - t1

#residualRaEpochWise = np.zeros( (len(searchFile),len(mjdBreakAt)+1) )
#residualDecEpochWise = np.zeros( (len(searchFile),len(mjdBreakAt)+1) )
medianResidualRaEpochWise = np.zeros(len(mjdBreakAt)+1)
medianResidualDecEpochWise = np.zeros(len(mjdBreakAt)+1)

#if there are atleast three epochs present
if (len(mjdBreakAt)>=3.0):
    previousVar2 = np.min(mjdSorted)
    for var1, var2 in enumerate(mjdBreakAt, start=0):
        mjdIndex = (searchFile['mjd']< var2) & (previousVar2 < searchFile['mjd'])
        mjdIndexInPixel = (objInPixelFile['mjd']<var2) & (previousVar2 < objInPixelFile['mjd'])
        print var1
        if any(mjdIndex):
            #print 'hehe'
            #print residualRaArray[mjdIndex]
            #residualRaEpochWise[:, var1] = residualRaArray[mjdIndex]
            medianResidualRaEpochWise[var1] = np.median(residualRaArray[mjdIndex])
            medianResidualDecEpochWise[var1] = np.median(residualDecArray[mjdIndex])
            #replacing old values with new values
            objInPixelFile['ra'][mjdIndexInPixel] -=  medianResidualRaEpochWise[var1]
            objInPixelFile['dec'][mjdIndexInPixel] -=  medianResidualDecEpochWise[var1]
        else:
            #print 'haha'
            #print residualRaArray[mjdIndex]
            medianResidualRaEpochWise[var1+1] = np.median(residualRaArray[mjdIndex])
            medianResidualDecEpochWise[var1+1] = np.median(residualDecArray[mjdIndex])
            #replacing old values with new values
            objInPixelFile['ra'][mjdIndexInPixel] -=  medianResidualRaEpochWise[var1+1]
            objInPixelFile['dec'][mjdIndexInPixel] -=  medianResidualDecEpochWise[var1+1]
            #newRa = objInPixelFile['ra'][mjdIndexInPixel] - medianResidualRaEpochWise[var1+1]
            #newDec = objInPixelFile['dec'][mjdIndexInPixel] - medianResidualDecEpochWise[var1+1]
        previousVar2 = var2


#calculate final ra/dec as the median of newRa/newDec values
finalRaArray = np.zeros(uniqueObjInPixelFile.size)
finalDecArray = np.zeros(uniqueObjInPixelFile.size)

currObj = objInPixelFile['obj_id'][0]
raObj = [objInPixelFile['ra'][0]]
decObj = [objInPixelFile['dec'][0]]
mjdObj = [objInPixelFile['mjd'][0]]

counter = 0
t2 = time()
for i in np.arange(objInPixelFile.size-1)+1:#np.arange(10)+1: #
    objID = objInPixelFile['obj_id'][i]
    ra = objInPixelFile['ra'][i]
    dec = objInPixelFile['dec'][i]
    if(objID == currObj):
        raObj.append(ra)
        decObj.append(dec)
    else:
        raObj=np.array(raObj)
        decObj=np.array(decObj)
        finalRaArray[counter] = np.median(raObj)
        finalDecArray[counter] = np.median(decObj)
        # print finalRaArray[counter]
        #print i
        #print counter
        currObj = objID
        raObj = [objInPixelFile['ra'][i]]
        decObj = [objInPixelFile['dec'][i]]
        counter +=1
print time()-t2
    
    












        

#find objects in the pixel of interest -- give pixel indices to all the objects -- pixel indices act like bins -- use these pixel indices to match the index of the pixel in question and VOILA! you have all the ra/dec values of the objects inside the pixel


 

















        

#JUNK CODE
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


#calculate phi and theta given ra and dec
phiForObj = (searchFile['ra'] * np.pi)/180 
thetaForObj = (90 - searchFile['dec'])* (np.pi/180)
pixelIndexForObj= hp.ang2pix(nside,thetaForObj, phiForObj)

phiForObj = (testFile['ra'] * np.pi)/180 
thetaForObj = (90 - testFile['dec'])* (np.pi/180)
pixelIndexForObj= hp.ang2pix(nside,thetaForObj, phiForObj)



#angSepMask =  sphdist( pixelRa[pickPixelNo], pixelDec[pickPixelNo],pixelRa, pixelDec ) <= (searchRadius*60)
#raIN = pixelRa[angSepMask]
#decIN = pixelDec[angSepMask]




#take their residuals

#group by epoch
