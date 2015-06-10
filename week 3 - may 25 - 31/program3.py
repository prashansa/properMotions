import numpy as np
from astropy.time import Time
from time import time
from esutil.coords import sphdist
import healpy as hp


noOfYears = 5
pixelWidth = 3 #arcminutes
searchRadius = 10 #arcminutes

#load data into a variable
testFile= np.load('/home/gupta/projects/propermotion/test.npy')

#save unique galaxy ids in another variable
uniqueObjIDfile = np.unique(testFile['obj_id'])


mjdSorted = np.sort(testFile['mjd'])
deltaT = mjdSorted[0:-1] - mjdSorted[1:]
tBreakAt = np.where(deltaT < (-100.0))[0]
mjdBreakAt = mjdSorted[tBreakAt+1]


medianRaArray=np.zeros(uniqueObjIDfile.size)
medianDecArray=np.zeros(uniqueObjIDfile.size)
residualRaArray= np.zeros(len(testFile))
residualDecArray= np.zeros(len(testFile))

currentObjID = testFile['obj_id'][0]
raObj = [testFile['ra'][0]]
decObj = [testFile['dec'][0]]
mjdObj = [testFile['mjd'][0]]

t0 = time()
objIDcounter = 0
pos=0
#it should start from the second row, so add one.
for i in np.arange(100)+1:#np.arange(testFile.size-1)+1:#
        print i
        objID = testFile['obj_id'][i]
        ra = testFile['ra'][i]
        dec = testFile['dec'][i]
        mjd = testFile['mjd'][i]
        if (objID == currentObjID):
                mjdObj.append(mjd)
                raObj.append(ra)
                decObj.append(dec)
                print 'same object now'
        else:
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
                print 'we are on the', objIDcounter,'object'
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
                print 'different object now'
print time()-t0



previousVar2 = np.min(mjdSorted)
for var1, var2 in enumerate(mjdBreakAt, start=0):
        mjdIndex = (mjdObj< var2) & (previousVar2 < mjdObj)
        if any(mjdIndex):
                residualRaEpochWise[var1].append(residualRa[var1])
                residualDecEpochWise[var1].append(residualRa[var1])
        else:
                #means it lies in the last epoch
                residualRaEpochWise[var1+1].append(residualRa[var1+1])
                residualDecEpochWise[var1+1].append(residualDec[var1+1])
























