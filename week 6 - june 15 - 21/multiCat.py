import numpy as np
from multiprocessing import Pool
from numpy.random import randn
from astropy.time import Time
from time import time
from esutil.coords import sphdist
import healpy as hp

''' GLOBAL TASKS'''

#variables
searchRadius = 10 #in arcminutes

#load data into a variable
#sorting is automatically being done by object ID
t0 = time()
testFile= np.sort(np.load('/home/gupta/projects/propermotion/test.npy'))
t1= time()
print t1 - t0

#save unique galaxy ids in another variable
uniqueObjIDfile = np.unique(testFile['obj_id'])

#the number of observations actually available for each object
noOfObs = np.zeros(len(uniqueObjIDfile)) 
  
#binning mjds for once
mjdSorted = np.sort(testFile['mjd'])
deltaT = mjdSorted[0:-1] - mjdSorted[1:]
tBreakAt = np.where(deltaT < (-100.0))[0]
mjdBreakAt = mjdSorted[tBreakAt+1]

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
pixelIndexForObjTest= hp.ang2pix(nside,thetaForObj, phiForObj)
#we want every pixel only once! 
pixelIndexForObj = np.unique(pixelIndexForObjTest)

#pick a pixel from the testFile
#pickPixelNo = pixelIndexForObj[len(pixelIndexForObj)/2]

#save ra and dec value of that pixel, so that we send only that value and not the whole array
#ra = pixelRa[pickPixelNo]
#dec = pixelDec[pickPixelNo]


#load residual ra and residual dec
#res_ra = ...
#res_dec = ...

#the list of pixels as well -- saved as pixelIndexForObj
#pixel_id_list = 

#pack parameters for workers
#Note: even if res_ra and res_dec are huge, they are passed as references, so "parameterListForPixel" does not use a lot of memory
parameterListForPixel = [(pickPixelNo, pixelRa[pickPixelNo], pixelDec[pickPixelNo], noOfObs, medianRaArray, medianDecArray, residualRaArray, residualDecArray) for pickPixelNo in pixelIndexForObj]

#start workers
pool = Pool(processes=8)

ti = time()
#chunk size is used to submit jobs in batches which reduces overhead 
iterator = pool.imap_unordered(pixelTasks, parameterListForPixel, chunksize=100)
ti2 = time()-ti
#def pixelTasks(residualRaArray, residualDecArray, medianRaArray, medianDecArray):
#def pixelTasks(pickPixelNo, ra, dec, residualRaArray, residualDecArray, noOfObs, medianRaArray, medianDecArray, testFile):
def pixelTasks(parameterListForPixel):
    """Workers will execute this function."""
    #unpack the input parameter list that is to say separate them into different arrays
    pixelNo, ra, dec, nObs, medRa, medDec, resRa, resDec = parameterListForPixel
    #define an index so that all the calculations can be put in place later
    #ind =
    global searchRadius
    #find objects within the searchRadius and whose noOfObs is atleast 3
    angSepMask =  (sphdist(ra, dec, medRa, medDec) <= (searchRadius*60)) & (nObs>=3.0)
    global testFile
    #separate these entire rows from the original testFile
    searchFile = testFile[angSepMask]
    #putting unique object IDs in this file
    uniqueSearchFile = np.unique(searchFile['obj_id'])
    #file containing data for objects inside pixel
    objInPixelFile = testFile[pixelIndexForObj == pixelNo]
    #file containing unique objIDs of the pixel
    uniqueObjInPixelFile = np.unique(objInPixelFile['obj_id'])
    global mjdBreakAt
    global mjdSorted
    medianResidualRaEpochWise = np.zeros(len(mjdBreakAt)+1)
    medianResidualDecEpochWise = np.zeros(len(mjdBreakAt)+1)
    #if there are atleast three epochs present
    if (len(mjdBreakAt)>=3.0):
        previousVar2 = np.min(mjdSorted)
        for var1, var2 in enumerate(mjdBreakAt, start=0):
            mjdIndex = (searchFile['mjd']< var2) & (previousVar2 < searchFile['mjd'])
            mjdIndexInPixel = (objInPixelFile['mjd']<var2) & (previousVar2 < objInPixelFile['mjd'])
            #print var1
            if any(mjdIndex):
                medianResidualRaEpochWise[var1] = np.median(resRa[mjdIndex])
                medianResidualDecEpochWise[var1] = np.median(resDec[mjdIndex])
                #replacing old values with new values
                objInPixelFile['ra'][mjdIndexInPixel] -=  medianResidualRaEpochWise[var1]
                objInPixelFile['dec'][mjdIndexInPixel] -=  medianResidualDecEpochWise[var1]
            else:
                medianResidualRaEpochWise[var1+1] = np.median(resRa[mjdIndex])
                medianResidualDecEpochWise[var1+1] = np.median(resDec[mjdIndex])
                #replacing old values with new values
                objInPixelFile['ra'][mjdIndexInPixel] -=  medianResidualRaEpochWise[var1+1]
                objInPixelFile['dec'][mjdIndexInPixel] -=  medianResidualDecEpochWise[var1+1]
            previousVar2 = var2
    #calculate final ra/dec as the median of newRa/newDec values
    finalRaArray = np.zeros(uniqueObjInPixelFile.size)
    finalDecArray = np.zeros(uniqueObjInPixelFile.size)
    #take the current values separately
    currObj = objInPixelFile['obj_id'][0]
    raObj = [objInPixelFile['ra'][0]]
    decObj = [objInPixelFile['dec'][0]]
    mjdObj = [objInPixelFile['mjd'][0]]
    #set variables
    counter = 0
    #t2 = time()
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
            currObj=objID
            raObj=[objInPixelFile['ra'][i]]
            decObj=[objInPixelFile['dec'][i]]
            counter +=1
    #t3 = time()
    #print 'time taken for calculating ra and dec values for a pixel is ', t3-t2
    return finalRaArray, finalDecArray 



# create storage arrays
raFinal = ...
decFinal = ...

# loop over parameter sets
# "it" is an iterator that returns the result in the "res" list
for result in iterator:
    raFinal[result[0]] = result[1]
    decFinal[result[0]] = result[2]

#in the above loop, a worker takes a set of (pixelID, res_ra, res_dec) parameters, executes the "pixelTasks" function using these parameters, and returns the results (res_ra_aux and res_dec_aux which are returned in result[1] and result[2]) *and* the index ind (which is returned in result[0]) that can be used to store these results into the correct rows in storage arrays

#terminate workers when done
pool.terminate()

#dump data
np.save('raFinal.npy', raFinal)
np.save('decFinal.npy', decFinal)





def calcMedianAndResiduals(testFile, uniqueObjIDfile):
    '''calculates median ra and dec along with residuals for a region of sky once'''
    #stores the median ra and dec values for each object
    medianRaArray=np.zeros(uniqueObjIDfile.size)
    medianDecArray=np.zeros(uniqueObjIDfile.size)
    #stores the residual values for each object and for each epoch
    residualRaArray= np.zeros(len(testFile))
    residualDecArray= np.zeros(len(testFile))
    #store the first object values in here
    currentObjID = testFile['obj_id'][0]
    raObj = [testFile['ra'][0]]
    decObj = [testFile['dec'][0]]
    mjdObj = [testFile['mjd'][0]]
    #set variables
    objIDcounter = 0
    pos=0
    t1 = time()
    #it should start from the second row, so add one.
    for i in  np.arange(testFile.size-1)+1:#np.arange(10)+1:#
        #print i
        objID = testFile['obj_id'][i]
        ra = testFile['ra'][i]
        dec = testFile['dec'][i]
        mjd = testFile['mjd'][i]
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
            raObj = [testFile['ra'][i]]
            decObj = [testFile['dec'][i]]
            mjdObj = [testFile['mjd'][i]]
    print 'total time taken is',time() - t1
    return medianRaArray, medianDecArray, residualRaArray, residualDecArray

