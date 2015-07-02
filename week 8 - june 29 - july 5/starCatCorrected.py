#processing a file with both stars and galaxies
#this code finds the residuals of galaxies for each epoch = ra/dec_epoch - ra/dec_final
#these residuals are used to update the ra/dec values for stars
#and then finally find ONE ra/dec value for each star for this catalogue

''' NOTE : OBJECT in general refers to a galaxy,  a star is referred explicitly as a STAR'''

'''MULTIPROCESSING + PYTABLES + SQLITE '''

import numpy as np
from time import time
import healpy as hp

from multiprocessing import Pool

import tables

#import sys
#import os

#http://research.majuric.org/trac/wiki/LargeSurveyDatabase
#from lsd import DB
#from lsd import bounds as lsdbounds

import sqlite3

# various functions are defined here
import starCatCorrectedFunctions

################################################################################


''' GLOBAL TASKS'''

#variables
searchRadius = 10 #in arcminutes

#load h5 file -- note the 's' in 'tableS.open'
#table is already sorted by object id
#this file contains stars and galaxies in the TriAnd region
#column called "gal" : gal = 1 for galaxies, gal = 0 for stars
t0 = time()
testH5file = tables.open_file('/home/bsesar/projects/PS1/proper_motion/prashansa/TriAnd_sample.h5')
t1= time()
print t1 - t0

# load the table
triand = testH5file.root.triand

# run calcMedianAndResiduals
objIDarray, medianRaArray, medianDecArray, noOfObsArray, residualRaArray, residualDecArray = starCatCorrectedFunctions.calcMedianAndResiduals(triand, False)


#select galaxies observed in at least 3 epochs, and junk the rest
#the for loop in the previous section is structured such that the median values are only stored when the no of obs is greater than 3, nevertheless we need this step to remove zero values from the rest of the places in the array.
goodGal = noOfObsArray >= 3.0
objIDarray = objIDarray[goodGal]
medianRaArray = medianRaArray[goodGal]
medianDecArray = medianDecArray[goodGal]

#binning mjds for once
maskForGal = triand.get_where_list('gal == 1') #already defined in calcmed..func
mjdValues = triand.col('mjd')[maskForGal]
mjdSorted = np.sort(mjdValues)
deltaT = mjdSorted[0:-1] - mjdSorted[1:]
tBreakAt = np.where(deltaT < (-100.0))[0]
mjdBreakAt = mjdSorted[tBreakAt+1]

# objIDs and MJDs that match the residual rows (need this for the pixelTasks function)
objIDs = triand.col('obj_id')[maskForGal] #already defined in calcmed..func
MJDs = triand.col('mjd')[maskForGal]

#this code pixelizes the sky
nside= 2**10

#nside2npix gives number of pixels for a given nside, create an array containing pixel indices
pixelIndexArray = np.arange(hp.nside2npix(nside))

#obtain angular coordinates corresponding to nside and the pixel indices
theta, phi = hp.pix2ang(nside, pixelIndexArray)

#pixel centers in degrees
pixelRa = 180*phi/np.pi
pixelDec = 90 - theta*180/np.pi

#store ra/dec values of all galaxies, to be passed to each pixel later
raGal = triand.col('ra')[maskForGal]
decGal = triand.col('dec')[maskForGal]

#store min and max values separately
minRaGal = triand.col('ra')[maskForGal].min()
maxRaGal = triand.col('ra')[maskForGal].max()
minDecGal = triand.col('dec')[maskForGal].min()
maxDecGal = triand.col('dec')[maskForGal].max()

# consider only pixels in the region where we have data while taking into account the 10 arcmin buffer
goodPixels = (pixelRa >= (minRaGal + 10./60)) & \
             (pixelRa <= (maxRaGal - 10./60)) & \
             (pixelDec >= (minDecGal + 10./60)) & \
             (pixelDec <= (maxDecGal - 10./60))

pixelRa = pixelRa[goodPixels]
pixelDec = pixelDec[goodPixels]
pixelIndexArray = pixelIndexArray[goodPixels]


#phiForObj   = (triand.col('ra')*np.pi)/180
#thetaForObj = (90 - triand.col('dec'))* (np.pi/180)
#pixelIndexForObj = h.ang2pix(nside, thetaForObj, phiForObj)

phiForObj  = (medianRaArray*np.pi)/180 
thetaForObj = (90 - medianDecArray)* (np.pi/180) 
pixelIndexForObj = hp.ang2pix(nside, thetaForObj, phiForObj) 

print "Going to process %d pixels." % pixelIndexArray.size

#pack parameters for workers
#Note: even if res_ra and res_dec are huge, they are passed as references, so "parameterListForPixel" does not use a lot of memory
parameterListForPixel = [(pickPixelNo, pixelRa[i], pixelDec[i], objIDarray, medianRaArray, medianDecArray, objIDs, MJDs, residualRaArray, residualDecArray, raGal, decGal, pixelIndexForObj, mjdBreakAt, mjdSorted ) for i, pickPixelNo in enumerate(pixelIndexArray)]

# for starters, test PixelTasks on only one pixel
print starCatCorrectedFunctions.pixelTasks(parameterListForPixel[10])

#start workers
pool = Pool(processes=8)
ti = time()
#chunk size is used to submit jobs in batches which reduces overhead
iterator = pool.imap_unordered(starCatCorrectedFunctions.pixelTasks, parameterListForPixel, chunksize=100)
print time()-ti


for param in parameterListForPixel:
    #res1, res2, res3, res4, res5 = starCatCorrectedFunctions.pixelTasks(param)
    res = starCatCorrectedFunctions.pixelTasks(param)
    print res[0]

res1, res2, res3, res4, res5 = starCatCorrectedFunctions.pixelTasks(parameterListForPixel)

    
