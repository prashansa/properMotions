''' code for both stars and galaxy -- save h5 files as combinedFile%d.h5 and databases as stardb%d '''

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

#db = DB('/home/bsesar/projects/PS1/DVO:/a41233d1/LSD/from_cfa')

if os.getenv('HOSTNAME') == 'aida41147':
	db = DB('/home/bsesar/projects/PS1/DVO:/a41233d1/LSD/from_cfa')
else:
	db = DB('/a41147d1/bsesar/projects/PS1/DVO:/a41233d1/LSD/from_cfa')

# define the query
query = 's.obj_id as obj_id, t.ra as ra, t.dec as dec, \
        t.raErr as raErr, t.decErr as decErr, \
        t.nObs as nObs, t.mjd as mjd, \
        -2.5*np.log10(s.mean(1)/s.mean_ap(1)) as sg_r, \
        -2.5*np.log10(s.mean(2)/s.mean_ap(2)) as sg_i \
        FROM dvo as t, ucal_fluxqy(matchedto=t, nmax=1, dmax=1.5) as s'

# table definition
class Star(tables.IsDescription):
    # in LSD, obj_id is a 64-bit unsigned integer,
    # but pytables cannot index 64-bit unsigned integers
    # so they are saved as 64-bit signed integers
    obj_id = tables.Int64Col(pos=0)
    ra  = tables.Float64Col(pos=1)
    dec = tables.Float64Col(pos=2)
    raErr  = tables.Float32Col(pos=3)
    decErr = tables.Float32Col(pos=4)
    nObs = tables.UInt8Col(pos=5)
    mjd  = tables.Float64Col(pos=6)
    gal  = tables.UIntCol(pos=7)
    
searchRadius = 10 #in arcminutes


def executeChunkWise(raMin, decMin, raMax, decMax, chunkNo):
    #this function will take in the bounds and return a filename corresponding to one set of bounds
    #every file name and table name will be appended by a chunk no, which will identify uniquely with a set of bounds.
    #tableName  = "aTable"
    h5fileName = "combinedFile" + "%d" %chunkNo
    #tableName = tableName + "%d" %chunkNo
    # open a pytable file
    h5file = tables.open_file("/a41233d1/gupta/combinedFiles/%s.h5" %h5fileName, mode = "w", title = "try1")
    # define compression filters
    filters = tables.Filters(complib='blosc', complevel=5)
    # create the table
    #create_table(where, name, obj, title, expectedrows, filters)
    #"/" refers to h5file.root object
    table = h5file.create_table('/', 'table1', Star, "tryTable", expectedrows=40563159, filters=filters)
    star = table.row
    # define selection bounds
    #gal long lower bound, gal lat lower bound, gal long upper bound, gal lat upper bound
    bounds = lsdbounds.rectangle(raMin, decMin, raMax, decMax, coordsys="equ")#(ra,dec) bottomleft; (ra,dec) topright
    bounds = lsdbounds.make_canonical(bounds)
    # query LSD for rows and store them into the pytable
    counterForRows = 0 
    dtype = table.colnames
    for row in db.query(query).iterate(bounds=bounds):
        star['obj_id']  = (row[0]).astype('i8')
        counterForRows += 1
        for j in range(len(dtype)-2):
            star[dtype[j+1]] = row[j+1]
        if (row[7] > 0.3) & (row[7] < 1.0) & (row[8] > 0.3) & (row[8] < 1.0):
            star['gal'] = 1
        else:
            star['gal'] = 0
        star.append()
        if (counterForRows % 10000 == 0): 
            table.flush()
    table.flush()
    noOfRowsInTable = table.nrows
    # create a full index on the obj_id column
    indexrows = table.cols.obj_id.create_csindex(filters=tables.Filters(complib='blosc', complevel=5))
    # create a new table that is sorted by obj_id
    sI = table.cols.obj_id
    table2 = tables.Table.copy(table, newname='aTable', overwrite=True, sortby=sI, checkCSI=True, propindexes=True)
    #  noOfRowsInTable2 = table2.nrows
    # delete the unsorted table
    h5file.root.table1.remove()    
    h5file.close()
    #should return the file name
    return h5fileName,noOfRowsInTable, counterForRows


def pixelTasksCombinedData(parameterList):
    #pixelNo is one number
    #pixelRa and dec are the ra/dec values of the pixel in consideration
    #galID are the object IDs of galaxies -- unique -- from the database
    #ra/decFinalGal are final ra/dec values of galaxies obtained from the database
    #galIDs -- galaxy IDs corresponding to all detections of galaxies
    #ra/decGal -- original ra and dec values of galaxies
    #starMjds --  mjd values corresponding to all the detections of stars
    pixelNo, pixelRa, pixelDec, galIDfinal, raFinalGal, decFinalGal, galIDs, raGal, decGal, galMjd , starIDs, starMjds, raStar, decStar, pixAllStar, mjdSorted, mjdBreakAt = parameterList


    angSepMask = (sphdist(pixelRa, pixelDec, raFinalGal, decFinalGal)) <= (searchRadius/60.0)
    #select unique galaxy ids within searchRadius
    uniqueGalIDinRadius = galIDfinal[angSepMask]
    raFinalGalInRadius  = raFinalGal[angSepMask]
    #print "raFinalGalInRadius.size",raFinalGalInRadius.size
    decFinalGalInRadius = decFinalGal[angSepMask]

    a = []
    b = []
    c = []
    d = []
    e = []

    
    for i in range(galIDs.size):
        if galIDs[i] in uniqueGalIDinRadius:
            a.append(galIDs[i])
            b.append(i)
            c.append(raGal[i])
            d.append(decGal[i])
            e.append(galMjd[i])

    galIDinRadius = np.array(a)
    raGalInRadius = np.array(c)
    decGalInRadius = np.array(d)
    galMjdInRadius = np.array(e)
    #match all detections of galaxies with those in radius
    #indexInRadius  = np.in1d(galIDs,uniqueGalIDinRadius)
    #obtain original ra and dec values for galaxies within radius
    #raGalInRadius  = raGal[indexInRadius]
    #decGalInRadius = decGal[indexInRadius]
    #galIDinRadius  = galIDs[indexInRadius]
    #print "galIDinRadius[1171070:1171083]",galIDinRadius[1171070:1171083]
    #galMjdInRadius = galMjd[indexInRadius] 
         
    currentGalID  = galIDinRadius[0]
    currentRaGal  = [raGalInRadius[0]]
    currentDecGal = [decGalInRadius[0]]
    
    galCounter = 0 
    position = 0 
    noOfDetectionsOfGal = galIDinRadius.size
    print "noOfDetectionsOfGal",noOfDetectionsOfGal
    
    offsetRaArray  = np.zeros(noOfDetectionsOfGal)
    offsetDecArray = np.zeros(noOfDetectionsOfGal)
    
    #run over all the objects within the searchRadius --- i dont know why we did a -1 previously ! think/ask ! :( -- there is no need and then we donot need to think about the last object separately! 
    t =  time()
    for k in np.arange(noOfDetectionsOfGal-1)+1:
        #print "on galaxy no", k 
        ID  = galIDinRadius[k]
        ra  = raGalInRadius[k]
        dec = decGalInRadius[k]
	#print "ID", ID
        if (ID == currentGalID):
                currentRaGal.append(ra)
                currentDecGal.append(dec)
                #print "same object now"
        else:
                #make them numpy arrays
                currentRaGal  = np.array(currentRaGal)
                currentDecGal = np.array(currentDecGal)
                #calculate offsets
	        #print "k", k
	        #print "galCounter", galCounter
                offsetRa  = currentRaGal - raFinalGalInRadius[galCounter]
                offsetDec = currentDecGal - decFinalGalInRadius[galCounter]
                #print "offsetRa", offsetRa
                #print "offsetDec", offsetDec
                #store offsets
                offsetRaArray[position:k]  = offsetRa
                offsetDecArray[position:k] = offsetDec
                #update counters and variables
	        galCounter+= 1
	        position   = k
	        currentGalID  = ID
	        currentRaGal  = [raGalInRadius[k]]
	        currentDecGal = [decGalInRadius[k]]
    #print "time taken to calculate offsets of galaxies from their final ra/dec", time()-t
    #print "galCounter after loop", galCounter	
	
    # select stars inside the pixel
    resIndexInPixel =  pixAllStar == pixelNo
    #uniqueStarIDinPixel = starID[resIndexInPixel]
    #match all detections of stars with those in the pixel
    #indexInPixel = np.in1d(starIDs, uniqueStarIDinPixel)
    #obtain original ra and dec values for stars within pixel
    raStarInPixel  = raStar[resIndexInPixel]
    decStarInPixel = decStar[resIndexInPixel]
    starIDInPixel  = starIDs[resIndexInPixel]
    mjdStarInPixel = starMjds[resIndexInPixel]

    medianRaOffsetEpochWise = np.zeros(len(mjdBreakAt)+1)
    medianDecOffsetEpochWise = np.zeros(len(mjdBreakAt)+1)
    avgMjd = np.zeros(len(mjdBreakAt)+1)

    if (len(mjdBreakAt)>= 3.0):
        previousVar2 = np.min(mjdSorted)
        for var1, var2 in enumerate(mjdBreakAt, start=0):
            mjdIndex = (galMjdInRadius < var2) & (previousVar2 <galMjdInRadius)
            mjdIndexInPixel = (mjdStarInPixel<var2) & (previousVar2 < mjdStarInPixel)
            if any(mjdIndex):
                offsetRaValues  = offsetRaArray[mjdIndex] 
                offsetDecValues = offsetDecArray[mjdIndex]
                medianRaOffsetEpochWise[var1]  = np.median(offsetRaValues)
                medianDecOffsetEpochWise[var1] = np.median(offsetDecValues)
                #update the ra/dec - replace old values
                raStarInPixel[mjdIndexInPixel] -= medianRaOffsetEpochWise[var1]
                decStarInPixel[mjdIndexInPixel]-= medianDecOffsetEpochWise[var1]
                avgMjd[var1] = (mjdStarInPixel[mjdIndexInPixel].max() - mjdStarInPixel[mjdIndexInPixel].min() )/2
            else:
                offsetRaValues  = offsetRaArray[mjdIndex] 
                offsetDecValues = offsetDecArray[mjdIndex]
                medianRaOffsetEpochWise[var1+1]  = np.median(offsetRaValues)
                medianDecOffsetEpochWise[var1+1] = np.median(offsetDecValues)
                #update the ra/dec - replace old values
                raStarInPixel[mjdIndexInPixel] -= medianRaOffsetEpochWise[var1]
                decStarInPixel[mjdIndexInPixel]-= medianDecOffsetEpochWise[var1]
                avgMjd[var1+1] = (mjdStarInPixel[mjdIndexInPixel].max() - mjdStarInPixel[mjdIndexInPixel].min() )/2
        previousVar2 = var2

    finalObjIDstar    = starIDInPixel
    finalRaArrayStar  = raStarInPixel
    finalDecArrayStar = decStarInPixel
    
    return finalObjIDstar, finalRaArrayStar, finalDecArrayStar
