import numpy as np
from numpy.random import randn
from astropy.time import Time
from time import time
from esutil.coords import sphdist
import healpy as hp
from esutil.numpy_util import match

from multiprocessing import Pool

import tables

import sys
import os

#http://research.majuric.org/trac/wiki/LargeSurveyDatabase
from lsd import DB
from lsd import bounds as lsdbounds

import sqlite3

#variables
searchRadius = 10 #in arcminutes


def dummyCalcMedianAndResiduals(triand, objID, objIDarray, medianRaArray, medianDecArray, noOfObsArray, objIDs, residualRaArray, residualDecArray ):
    #triand is the table
    #objID is the object we are looking at
    #objIDs are the size of residual arrays
    #the rest are arrays obtained from the calc function
    
    k   = objIDarray == objID
    ra  = medianRaArray[k]
    dec = medianDecArray[k]
    obs = noOfObsArray[k]
    
    l  = objIDs == objID
    resRa = residualRaArray[l]
    resDec = residualDecArray[l]
    
    
    print "ra, dec, noOfObs, residualRa/Dec values for object", objID, "are:", ra, dec, obs, resRa, resDec

    mask   = triand.get_where_list('obj_id == objID')
    objids = triand.col('obj_id')[mask]
    ra  =  np.array(triand.col('ra')[mask])
    dec =  np.array(triand.col('dec')[mask])

    medianRaVal  = np.median(ra)
    medianDecVal = np.median(dec)

    residualRa  = ra  - medianRaVal
    residualDec = dec - medianDecVal

    return medianRaVal, medianDecVal, residualRa, residualDec




def calcMedianAndResiduals(table, debug):
    '''calculates median ra and dec along with residuals for a region of sky once'''
    '''table: pytables table that contains positions of PS1 objects'''
    #select rows that contain galaxies
    maskForGal = table.get_where_list('gal == 1')

    #total number of *detections* of galaxies
    noOfDetOfGal = np.where(maskForGal)[0].size
    #unique number of galaxies
    noOfGal = np.unique(table.col('obj_id')[maskForGal]).size

    #store the median ra and dec values for each galaxy 
    medianRaArray  = np.zeros(noOfGal)
    medianDecArray = np.zeros(noOfGal)

    #store the residual values for each galaxy and for each epoch
    residualRaArray  = np.zeros(noOfDetOfGal)
    residualDecArray = np.zeros(noOfDetOfGal)
    
    #store the number of observations for each galaxy 
    noOfObsArray = np.zeros(noOfGal, dtype='u1')
    
    #store the objIDs of each galaxy -- we run it in the loop just for debugging purposes 
    #needs 64-bit unsigned integer datatype
    objIDarray = np.zeros(noOfGal, dtype='i8')
    
    # pluck out objIDs, RAs, DECs, and MJDs for each detection of galaxy
    # (works *much* faster than accessing individual pytable rows)
    objIDs  = table.col('obj_id')[maskForGal]
    raObjs  = table.col('ra')[maskForGal]
    decObjs = table.col('dec')[maskForGal]
    mjdObjs = table.col('mjd')[maskForGal]


    try:
        # store the first object values in here
        currObjID = objIDs[0]
        currRaObj  = [raObjs[0]]
        currDecObj = [decObjs[0]]
        currMjdObj = [mjdObjs[0]]
        noOfObsArray[0] = 1

        #set variables
        objIDcounter = 0
        pos = 0
        t1 = time()

        #it should start from the second row, so add one.
        for i in np.arange(noOfDetOfGal-1)+1:
        #for i in np.arange(30-1)+1:
            if debug: print 'on the row number',i
            objID = objIDs[i]
            ra    = raObjs[i]
            dec   = decObjs[i]
            mjd   = mjdObjs[i]
            if debug: print 'object ID',objID
            if debug: print 'ra of object',ra
            if debug: print 'dec of object',dec
            if debug: print 'mjd of object',mjd
            if (objID == currObjID):
                currRaObj.append(ra)
                currDecObj.append(dec)
                currMjdObj.append(mjd)
                if debug: print 'same object now'
                noOfObsArray[objIDcounter] += 1
            else:
                if debug: print 'different object now'
                if debug: print 'run functions for the', objIDcounter, 'object'
                if(noOfObsArray[objIDcounter] >= 2.0):
                    if debug: print 'no of obs ', noOfObsArray[objIDcounter]
                    # convert list into array
                    currRaObj = np.array(currRaObj)
                    currDecObj = np.array(currDecObj)
                    if debug: print 'current objID ra dec mjd :',currObjID, currRaObj, currDecObj, currMjdObj
                    # calculate and store medians
                    medianRa = np.median(currRaObj)
                    medianDec = np.median(currDecObj)
                    if debug: print 'medianRa is', medianRa
                    if debug: print 'medianDec is', medianDec
                    # store median value for each object
                    medianRaArray[objIDcounter] = medianRa
                    medianDecArray[objIDcounter] = medianDec
                    objIDarray[objIDcounter] = currObjID
                    # calculate residual
                    residualRa = currRaObj - medianRa
                    residualDec = currDecObj - medianDec
                    if debug: print 'residualRa is', residualRa
                    if debug: print 'residualDec is', residualDec
                    if debug: print 'pos is', pos
                    if debug: print 'i is ',i
                    # store residuals in separate arrays
                    residualRaArray[pos:i] = residualRa
                    residualDecArray[pos:i] = residualDec
                # store the objID no matter what
                objIDarray[objIDcounter] = currObjID
                objIDcounter += 1
                currObjID = objID
                pos = i
                # start new ra, dec, and mjd arrays and
                # don't forget to increment the noOfObsArray for this object
                currRaObj = [raObjs[i]]
                currDecObj = [decObjs[i]]
                currMjdObj = [mjdObjs[i]]
                noOfObsArray[objIDcounter] += 1

        # after going through all of the rows, make sure to process the last
        # object, if it has enough observations
        i +=1
        print 'processing the last object now'
        if(noOfObsArray[objIDcounter] >= 2.0):
            if debug: print 'obs greater than three, adding the last object too'
            currRaObj = np.array(currRaObj)
            currDecObj = np.array(currDecObj)
            # calculate and store medians
            medianRa = np.median(currRaObj)
            medianDec = np.median(currDecObj)
            if debug: print 'medianRa is', medianRa
            if debug: print 'medianDec is', medianDec
            medianRaArray[objIDcounter] = medianRa
            medianDecArray[objIDcounter] = medianDec
            # calculate residual
            residualRa = currRaObj - medianRa
            residualDec = currDecObj - medianDec
            if debug: print 'residualRa is', residualRa
            if debug: print 'residualDec is', residualDec
            if debug: print 'pos is', pos
            if debug: print 'i is ',i
            residualRaArray[pos : i] = residualRa
            residualDecArray[pos : i] = residualDec
            # also need to store the objID, no matter what
            objIDarray[objIDcounter] = currObjID
       
        print 'total time taken is', time() - t1


    except ValueError:
        print "currObjID medianRa medianDec residualRa residualDec, pos, i, objIDcounter", currObjID, medianRa, medianDec, residualRa, residualDec, pos, i, objIDcounter
        pass
        
    return objIDarray, medianRaArray, medianDecArray, noOfObsArray, residualRaArray, residualDecArray


def pixelTasks(parameterListForPixel):
    """Workers will execute this function."""
    #unpack the input parameter list (i.e., separate them into different arrays)
    #objID -- contains the good galaxy IDs -- unique no of times
    #objIDs -- contains all the galaxy IDs corresponding to all the detections of galaxies
    #medRa and medDec -- contain good median values 
    #resRa and resDec -- contain the residual values for each detection of galaxy-- not good
    #ra/decAll -- ra/dec values for all the galaxies -- not good
    #pixAll -- pixel indices for all the galaxies -- not good
    pixelNo, pixRa, pixDec, objID, medRa, medDec, objIDs, MJDs, resRa, resDec, raAll, decAll, pixAll,mjdBreakAt, mjdSorted  = parameterListForPixel

    #print 'pixelNo', pixelNo
    #print 'pixAll',pixAll
    #global searchRadius
    #global mjdBreakAt
    #global mjdSorted

    #find objects within the searchRadius
    angSepMask = sphdist(pixRa, pixDec, medRa, medDec) <= (searchRadius*60)
    objInRadius = objID[angSepMask]
    #print 'objInRadius',objInRadius

    #this function matches the given arrays, returns a boolean array of the size of first array.
    #this index points to rows that contain residuals (and MJDs) of galaxy detections within search radius
    resIndexInRadius = np.in1d(objIDs, objInRadius)

    #print 'objIDs[resIndexInRadius]',objIDs[resIndexInRadius]
    #print 'len(objIDs[resIndexInRadius])',len(objIDs[resIndexInRadius])

    #select unique obj ids in radius
    #uniqueObjInRadius = np.unique(objIDs[resIndexInRadius])

    #select the residual values corresponding to the objects in radius
    resRaValuesInRadius = resRa[resIndexInRadius]
    resDecValuesInRadius = resDec[resIndexInRadius]
    #print 'resRaValuesInRadius', resRaValuesInRadius
    #print 'resDecValuesInRadius', resDecValuesInRadius
    #print 'len(resRaValuesInRadius)', len(resRaValuesInRadius)
    #print 'len(resDecValuesInRadius)', len(resDecValuesInRadius)
    #print 'np.count_nonzero(resRaValuesInRadius)', np.count_nonzero(resRaValuesInRadius)
    #print 'np.count_nonzero(resDecValuesInRadius)', np.count_nonzero(resDecValuesInRadius)

    #select mjds in the radius
    mjdInRadius = MJDs[resIndexInRadius]
    #print 'mjdInRadius',mjdInRadius

    #for storing median values for each epoch
    medianResidualRaEpochWise  = np.zeros(len(mjdBreakAt)+1)
    medianResidualDecEpochWise = np.zeros(len(mjdBreakAt)+1)

    #phiForObj  = (medRa*np.pi)/180 # should be done outside
    #thetaForObj = (90 - medDec)* (np.pi/180) # outside
    #pixelIndexForObj = h.ang2pix(nside, thetaForObj, phiForObj) # calc outside, pass here
    
    indexInPixel = pixAll == pixelNo
    objIDinPixel = objID[indexInPixel]
    #medRaInPixel = medRa[indexInPixel]
    #medDecInPixel = medDec[indexInPixel]

    #note : objIDinPixel and ra/decInPixel are not the same size, so we need to create another array
    resIndexInPixel = np.in1d(objIDs, objIDinPixel)
    #print 'len(objIDs)', len(objIDs)
    #print 'len(objIDinPixel)',len(objIDinPixel)
    #raInPixel = triand.col('ra')[maskForGal][resIndexInPixel] # access RA and Dec from outside - processing table each time for each pixel is bad
    #decInPixel = triand.col('dec')[maskForGal][resIndexInPixel]
    raInPixel = raAll[resIndexInPixel]
    decInPixel = decAll[resIndexInPixel]
    objInPixel = objIDs[resIndexInPixel]
    # print 'raInPixel', raInPixel
    # print 'decInPixel', decInPixel
    # print 'objInPixel', objInPixel
    # print 'len(raInPixel)', len(raInPixel)
    # print 'len(decInPixel)', len(decInPixel)
    # print 'len(objInPixel)', len(objInPixel)
    
    # sort data in the pixel by objInPixel
    assert((objInPixel[1:] - objInPixel[0:-1]).all >= 0)
    mjdInPixel = MJDs[resIndexInPixel]
    #print mjdInPixel


    #phiForObj   = (triand.col('ra')*np.pi)/180
    #thetaForObj = (90 - triand.col('dec'))* (np.pi/180)
    #pixelIndexForObj = h.ang2pix(nside, thetaForObj, phiForObj)
    #raInRadius = triand.col('ra')[maskForGal][resIndexInRadius]
    #decInRadius = triand.col('dec')[maskForGal][resIndexInRadius]
    #pixelIndexInRadius = pixelIndexForObj[resIndexInRadius]
    #inPixelMask = pixelIndexInRadius == pixelNo
    #raInPixel = raInRadius[inPixelMask]
    #decInPixel = decInRadius[inPixelMask]

  
    #if there are atleast three epochs present
    if (len(mjdBreakAt)>=3.0):
        previousVar2 = np.min(mjdSorted)
        for var1, var2 in enumerate(mjdBreakAt, start=0):
            mjdIndex = (mjdInRadius< var2) & (previousVar2 < mjdInRadius)
            mjdIndexInPixel = (mjdInPixel<var2) & (previousVar2 < mjdInPixel)
            #print var1
            if any(mjdIndex):
                resRaValues  = resRaValuesInRadius[mjdIndex]
                resDecValues = resDecValuesInRadius[mjdIndex]
                #print 'resRaValues', resRaValues
                #print 'resDecValues', resDecValues
                #print 'len(resRaValues)',len(resRaValues)
                #print 'len(resDecValues)',len(resDecValues)
                #print 'np.count_nonzero(resRaValues)', np.count_nonzero(resRaValues)
                #print 'np.count_nonzero(resDecValues)', np.count_nonzero(resDecValues)
                #resRaValuesTrue = resRaValues[np.where(resRaValues)]
                medianResidualRaEpochWise[var1]  = np.median(resRaValues)
                medianResidualDecEpochWise[var1] = np.median(resDecValues)
                #print 'medianResidualRaEpochWise[',var1,']', medianResidualRaEpochWise[var1]
                #print 'medianResidualDecEpochWise[',var1,']', medianResidualDecEpochWise[var1] 
                #replacing old values with new values
                raInPixel[mjdIndexInPixel] -=  medianResidualRaEpochWise[var1]
                decInPixel[mjdIndexInPixel] -=  medianResidualDecEpochWise[var1]
            else:
                resRaValues  = resRaValuesInRadius[mjdIndex]
                resDecValues = resDecValuesInRadius[mjdIndex]
                medianResidualRaEpochWise[var1+1]  = np.median(resRaValues)
                medianResidualDecEpochWise[var1+1] = np.median(resDecValues)
                #print 'medianResidualRaEpochWise[',var1+1,']', medianResidualRaEpochWise[var1+1]
                #print 'medianResidualDecEpochWise[',var1+1,']', medianResidualDecEpochWise[var1+1] 
                #replacing old values with new values
                raInPixel[mjdIndexInPixel] -=  medianResidualRaEpochWise[var1+1]
                decInPixel[mjdIndexInPixel] -=  medianResidualDecEpochWise[var1+1]
            previousVar2 = var2

    #calculate final ra/dec as the median of newRa/newDec values
    finalRaArray  = np.zeros(objIDinPixel.size)
    finalDecArray = np.zeros(objIDinPixel.size)

    medianRaError  = np.zeros(objIDinPixel.size)
    medianDecError = np.zeros(objIDinPixel.size)
    
    #maskForStar = table.get_where_list('gal == 0')
    #total number of *detections* of stars
    #noOfDetOfStar = np.where(maskForStar)[0].size
    #arrays to store deltas for stars
    #deltaRaArray  = np.zeros(noOfDetOfStar)
    #deltaDecArray = np.zeros(noOfDetOfStar)

    #if(pixelNo==74):
    #    print "Its pixel no 74"
    #    print "objInPixel:", objInPixel
    #    print "raInPixel:", raInPixel
    #    print "decInPixel:", decInPixel

    try:
            #take the current values separately
        obj   = objInPixel[0]
        raObj = [raInPixel[0]]
        decObj = [decInPixel[0]]
        mjdObj = [mjdInPixel[0]]
        
        #set variables
        counter = 0
        pos2 = 0
        t2   = time()
        for i in np.arange(objInPixel.size-1)+1:#np.arange(10)+1: #
            #print 'on the row number',i
            objID = objInPixel[i]
            ra  = raInPixel[i]
            dec = decInPixel[i]
            if(objID == obj):
                raObj.append(ra)
                decObj.append(dec)
            else:
                raObj = np.array(raObj)
                decObj= np.array(decObj)
                #print 'obj', obj
                #print 'raObj', raObj
                #print 'decObj', decObj            
                #find the final ra and dec values
                finalRaArray[counter]  = np.median(raObj)
                finalDecArray[counter] = np.median(decObj)
                #calculate rms values for ra/dec
                rmsRa  = 0.741*(np.percentile(raObj, 75) - np.percentile(raObj, 25))
                rmsDec = 0.741*(np.percentile(decObj, 75) - np.percentile(decObj, 25))
                #calculate uncertainity in median coordinates
                #print 'rmsRa', rmsRa
                #print 'rmsDec', rmsDec
                #print 'rmsRa.size = ', rmsRa.size
                #print 'rmsDec.size=', rmsDec.size
                #print 'raObj.size', raObj.size
                #print 'decObj.size', decObj.size
                medianRaError[counter]  = np.sqrt((np.pi/2)/(raObj.size-1))*rmsRa
                medianDecError[counter] = np.sqrt((np.pi/2)/(decObj.size-1))*rmsDec
                # print finalRaArray[counter]
                #print i
                #print 'counter', counter
                #calculate how much the galaxy moved in each epoch
                #recall : the final ra and dec values were obtained after removing offsets epochwise and obj_id wise
                #thus implying that these values make the galaxies static as required
                #deltaRa  = raObj2 - finalRaArray[counter]
                #deltaDec = decObj2 - finalDecArray[counter]
                #deltaRaArray[pos2:i]  = deltaRa
                #deltaDecArray[pos2:i] = deltaDec
                #move to the next obj
                obj  = objID
                pos2 = i
                counter +=1
                #make them lists again
                raObj  =[raInPixel[i]]
                decObj =[decInPixel[i]]
                #raObj2=[objInPixelFile2['ra'][i]]
                #decObj2=[objInPixelFile2['dec'][i]]
                
        #processing the last object now
        raObj  = np.array(raObj)
        decObj = np.array(decObj)
        finalRaArray[counter]  = np.median(raObj)
        finalDecArray[counter] = np.median(decObj)
        #calculate rms values for ra/dec
        rmsRa  = 0.741*(np.percentile(raObj, 75) - np.percentile(raObj, 25))
        rmsDec = 0.741*(np.percentile(decObj, 75) - np.percentile(decObj, 25))
        medianRaError[counter]  = np.sqrt((np.pi/2)/(raObj.size-1))*rmsRa
        medianDecError[counter] = np.sqrt((np.pi/2)/(decObj.size-1))*rmsDec
        
        t3 = time()
        #print 'time taken for calculating ra and dec values for a pixel is ', t3-t2
        return objIDinPixel, finalRaArray, finalDecArray, medianRaError, medianDecError

    except IndexError:
        #pass
        print "jhahaa"
        print pixelNo
        #return objIDinPixel, finalRaArray, finalDecArray, medianRaError, medianDecError
        return np.array(objIDinPixel), np.array(finalRaArray), np.array(finalDecArray), np.array(medianRaError), np.array(medianDecError)

        
    #final updation for STARS now
    #if (len(mjdBreakAt)>=3.0):
    #    previousVar2 = np.min(mjdSorted)
    #    for var1, var2 in enumerate(mjdBreakAt, start=0):
    #        mjdIndexInPixel = (objInPixelFile['mjd']<var2) & (previousVar2 < objInPixelFile['mjd'])
    #        if any(mjdIndex):
    #            medianDeltaRaEpochWise[var1] = np.median(deltaRaArray[mjdIndexInPixel])
    #            medianDeltaDecEpochWise[var1] = np.median(deltaDecArray[mjdIndexInPixel])
    #            #replacing old values with new values
    #            starInPixelFile['ra'][mjdIndexInPixel] -=  medianDeltaRaEpochWise[var1]
    #            starInPixelFile['dec'][mjdIndexInPixel] -=  medianDeltaDecEpochWise[var1]
    #        else:
    #            medianDeltaRaEpochWise[var1+1] = np.median(deltaRaArray[mjdIndexInPixel])
    #            medianDeltaDecEpochWise[var1+1] = np.median(deltaDecArray[mjdIndexInPixel])
    #            #replacing old values with new values
    #            starInPixelFile['ra'][mjdIndexInPixel] -=  medianDeltaRaEpochWise[var1+1]
    #            starInPixelFile['dec'][mjdIndexInPixel] -=  medianDeltaDecEpochWise[var1+1]
    #        previousVar2 = var2


# def pushDataInTable(finalData):

#     db = sqlite3.connect('mydb')

#     cursor = db.cursor()
#     cursor.execute(''' DROP TABLE IF EXISTS finalData''')
#     db.commit()

#     cursor = db.cursor()
#     cursor.execute(''' CREATE TABLE finalData(objID INT NOT NULL, ra REAL NOT NULL, dec REAL NOT NULL, raError REAL NOT NULL, decError REAL NOT NULL)''')

#     db.commit()

#     dat = starCatCorrectedFunctions.pixelTasks(parameterListForPixel)
    
#     cursor.executemany('INSERT INTO finalData (objID, ra, dec, raError, decError) VALUES (?, ?, ?, ?, ?)', dat)
    





