import numpy as np
from time import time
import healpy as hp
from multiprocessing import Pool
import tables
import sqlite3
import sys
import os
from lsd import DB
from lsd import bounds as lsdbounds

# various functions are defined here
import combinedFunctions




#####################################################################


#to see how floating point errors will be handled
#invalid = invalid floating point operations
#divide = division by zero
np.seterr(invalid='ignore', divide='ignore')

# connect to folders with LSD data
#keywords (case insensitive) select, from, where, into, as
#db = DB('/ssd-raid0/bsesar:/a41233d1/LSD/from_cfa')
#db = DB('/ssd-raid0/bsesar/projects/PS1/DVO:/a41233d1/LSD/from_cfa')
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
        FROM dvo as t, ucal_fluxqy(matchedto=t, nmax=1, dmax=1.5) as s \
        WHERE (sg_r > 0.3) & (sg_i > 0.3) & \
        (sg_r < 1.0) & (sg_i < 1.0)'

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
    
#specifying size of chunks, assume they are squares so both ra/dec have same breadth
chunkSize = 10.0 #degrees
#for future - once you have implemented this, you might also want to keep the no of obj same in each chunk of the sky -- so your chunksize would depend on that!
someBigNumber = 1000

raMin  = 0.0 #ra.min()#
decMin = 5.0 #dec.min()#
raMax  = 360.0 #ra.max()#
decMax = 90.0 #dec.max()#

#arrays to contain the bounds
#i donot know the size of the arrays, so i specify some big number! 
raMinArray  = np.zeros(someBigNumber)
raMaxArray  = np.zeros(someBigNumber)
decMinArray = np.zeros(someBigNumber)
decMaxArray = np.zeros(someBigNumber)


#assign the current min/max ra/dec values to the 
currentRaMin  = raMin
currentDecMin = decMin
currentDecMax = decMin + chunkSize
decAvg = (currentDecMin + currentDecMax)/2.0
currentRaMax  = raMin + chunkSize/np.cos(np.radians(decAvg))

raMinArray[0] = raMin
decMinArray[0]= decMin
raMaxArray[0] = currentRaMax
decMaxArray[0]= currentDecMax

chunkNo = 0
counter = 0

#loop to calculate bounds for the chunks of sky 
while(currentDecMax < decMax):
    decAvg = (currentDecMin + currentDecMax)/2.0
    print "decAvg", decAvg
    print "cosine of decAvg", np.cos(np.radians(decAvg))
    counter = 0 
    while(currentRaMax < raMax):
        chunkNo += 1
        raMinArray[chunkNo] = currentRaMax - ((40.0/60.0)/np.cos(np.radians(decAvg)))
        raMaxArray[chunkNo] = raMinArray[chunkNo] + chunkSize/np.cos(np.radians(decAvg))
        decMinArray[chunkNo]= currentDecMin
        decMaxArray[chunkNo]= currentDecMax
        counter +=1# to count after how many chunks in ra do we reach the end
        if (raMaxArray[chunkNo] > raMax):
            print "hi", counter
            raMaxArray[chunkNo] = raMaxArray[chunkNo] - raMax
            raMinArray[chunkNo] = raMinArray[chunkNo] - raMax
            break #break statement breaks out of the smallest enclosing of the while/for loop
        if (decMaxArray[chunkNo] > decMax):
            decMaxArray[chunkNo] = decMax
        currentRaMax = raMaxArray[chunkNo]
    currentRaMin  = raMin
    currentRaMax  = raMin + chunkSize/np.cos(np.radians(decAvg))
    currentDecMin = currentDecMax - (40.0/60.0)
    currentDecMax = currentDecMin + chunkSize
    print "chunk no", chunkNo  


#shorten the arrays to the required length
#add one so that it takes elements upto the chunkNo
raMinArray  = raMinArray[0:chunkNo+1]
raMaxArray  = raMaxArray[0:chunkNo+1]
decMinArray = decMinArray[0:chunkNo+1]
decMaxArray = decMaxArray[0:chunkNo+1]

#this code pixelizes the sky
nside = 2**10
#nside2npix gives number of pixels for a given nside, create an array containing pixel indices
#pixelIndexArray = np.arange(hp.nside2npix(nside))
#obtain angular coordinates corresponding to nside and the pixel indices
#theta, phi = hp.pix2ang(nside, pixelIndexArray)
#pixel centers in degrees
#pixelRa  = 180*phi/np.pi
#pixelDec = 90 - theta*180/np.pi

#we have the bounds, now loop over all of them
for chunkNo, value in enumerate(raMinArray):
    raMin  = value
    raMax  = raMaxArray[chunkNo]
    decMin = decMinArray[chunkNo]
    decMax = decMaxArray[chunkNo]
    print "sending chunkNo", chunkNo
    print "bounds = raMin, decMin, raMax, decMax = ",raMin, decMin, raMax, decMax 
    h5fileName,noOfRows, countForRows = combinedFunctions.executeChunkWise(raMin, decMin, raMax, decMax, chunkNo)
    print h5fileName
    if (noOfRows!= countForRows):
        "error : the no. of rows dont match for chunk no", chunkNo, "with bounds", raMin, decMin, raMax, decMax
        break
    t0 = time()
    testH5file = tables.open_file('/home/gupta/projects/propermotion/files/%s.h5' %h5fileName)
    # load the table
    table = testH5file.root.aTable
    t1= time()
    print "time taken to load h5file and table", t1 - t0
    # run calcMedianAndResiduals
    t2 = time()
    objIDarray, medianRaArray, medianDecArray, noOfObsArray, residualRaArray, residualDecArray = combinedFunctions.calcMedianAndResiduals(table, False)
    print "time taken to run calcMedianAndResiduals is ", time()-t2
    print "len(residualRaArray)", len(residualRaArray)
    print "len(residualDecArray)", len(residualDecArray)
    #select galaxies observed in at least 3 epochs, and junk the rest
    #the for loop in the previous section is structured such that the median values are only stored when the no of obs is greater than 3, nevertheless we need this step to remove zero values from the rest of the places in the array.
    goodGal = noOfObsArray >= 2.0
    objIDarray    = objIDarray[goodGal]
    medianRaArray = medianRaArray[goodGal]
    medianDecArray= medianDecArray[goodGal]
    noOfObsArray  = noOfObsArray[goodGal]
    #binning mjds for once
    maskForGal = table.get_where_list('gal == 1') #already defined in calcmed..func
    mjdValues  = table.col('mjd')[maskForGal]
    mjdSorted  = np.sort(mjdValues)
    deltaT = mjdSorted[0:-1] - mjdSorted[1:]
    tBreakAt   = np.where(deltaT < (-100.0))[0]
    mjdBreakAt = mjdSorted[tBreakAt+1]
    # objIDs and MJDs that match the residual rows (need this for the pixelTasks function)
    objIDs = table.col('obj_id')[maskForGal] #already defined in calcmed..func
    MJDs   = table.col('mjd')[maskForGal]
    #store ra/dec values of all galaxies, to be passed to each pixel later
    raGal  = table.col('ra')[maskForGal]
    decGal = table.col('dec')[maskForGal]

    #store min and max values separately
    minRaGal  = table.col('ra')[maskForGal].min()
    maxRaGal  = table.col('ra')[maskForGal].max()
    minDecGal = table.col('dec')[maskForGal].min()
    maxDecGal = table.col('dec')[maskForGal].max()

    # consider only pixels in the region where we have data while taking into account the 10 arcmin buffer
    #goodPixels = (pixelRa >= (minRaGal + 10./60)) & \
    #             (pixelRa <= (maxRaGal - 10./60)) & \
    #             (pixelDec >= (minDecGal + 10./60)) & \
    #             (pixelDec <= (maxDecGal - 10./60))

    #pixelRa  = pixelRa[goodPixels]
    #pixelDec = pixelDec[goodPixels]
    #pixelIndexArray = pixelIndexArray[goodPixels]

    phiForObj   = (medianRaArray*np.pi)/180 
    thetaForObj = (90 - medianDecArray)* (np.pi/180) 
    pixelIndexForObj = hp.ang2pix(nside, thetaForObj, phiForObj) 

    goodPixels = (medianRaArray >= (minRaGal + 10./60)) & \
                 (medianRaArray <= (maxRaGal - 10./60)) & \
                 (medianDecArray >= (minDecGal + 10./60)) & \
                 (medianDecArray <= (maxDecGal - 10./60))
    
    pixelIndexArray = np.unique(pixelIndexForObj[goodPixels])
    theta, phi = hp.pix2ang(nside, pixelIndexArray)
    pixelRa  = 180*phi/np.pi
    pixelDec = 90 - theta*180/np.pi


    print "Going to process %d pixels." % pixelIndexArray.size
    
    #pack parameters for workers
    #Note: even if res_ra and res_dec are huge, they are passed as references, so "parameterListForPixel" does not use a lot of memory
    t3 = time()
    
    parameterListForPixel = [(pickPixelNo, pixelRa[i], pixelDec[i], objIDarray, medianRaArray, medianDecArray, objIDs, MJDs, residualRaArray, residualDecArray, raGal, decGal, pixelIndexForObj, mjdBreakAt, mjdSorted ) for i, pickPixelNo in enumerate(pixelIndexArray)]

    db = sqlite3.connect('mydb%d'%chunkNo)

    cursor = db.cursor()
    cursor.execute('''DROP TABLE IF EXISTS finalData''')
    db.commit()
    
    cursor = db.cursor()
    cursor.execute( '''CREATE TABLE finalData(objID INT NOT NULL, ra REAL NOT NULL, dec REAL NOT NULL, raError REAL NOT NULL, decError REAL NOT NULL)''')
    db.commit()
        
    #start workers
    pool = Pool(processes=24)
    ti = time()
    #chunk size is used to submit jobs in batches which reduces overhead
    iterator = pool.imap_unordered(combinedFunctions.pixelTasks, parameterListForPixel, chunksize=100)

    counter = 0
    noOfPixelsIterated = 0
    dat1 = []
    dat2 = []
    dat  = []
    
    # for res in iterator:
    #     objID_, ra_, dec_, rErr_, dErr_ = res
    #     noOfPixelsIterated += 1
    #     #so that you donot encounter empty pixels! 
    #     if (objID_.size != 0):
    #         dat = [(int(objID),float(ra), float(dec), float(rErr), float(dErr)) for objID, ra, dec, rErr, dErr in zip(objID_, ra_, dec_, rErr_, dErr_)]
    #         counter = counter + len(dat)
    #         if (noOfPixelsIterated % 1000 == 0):
    #             cursor.executemany('INSERT INTO finalData (objID, ra, dec, raError, decError) VALUES (?, ?, ?, ?, ?)', dat)
    #             db.commit()
    #             print "committed now"
    #     if (objID_.size ==0):
    #         print "zero array"

    # print "counter", counter

    for res in iterator:
    for param in parameterListForPixel:
        res = combinedFunctions.pixelTasks(param)
        objID_, ra_, dec_, rErr_, dErr_ = res
        noOfPixelsIterated += 1
        #so that you donot encounter empty pixels! 
        if (objID_.size != 0):
            dat = dat + [(int(objID),float(ra), float(dec), float(rErr), float(dErr)) for objID, ra, dec, rErr, dErr in zip(objID_, ra_, dec_, rErr_, dErr_)]
            counter = counter + len(objID_)
            if (noOfPixelsIterated % 1000 == 0):
                cursor.executemany('INSERT INTO finalData (objID, ra, dec, raError, decError) VALUES (?, ?, ?, ?, ?)', dat)
                db.commit()
                print "committed now"
                dat = []
        if (objID_.size ==0):
            print "zero array"
    if (noOfPixelsIterated % 1000 != 0):
        cursor.executemany('INSERT INTO finalData (objID, ra, dec, raError, decError) VALUES (?, ?, ?, ?, ?)', dat)
        db.commit()
        print "committing last"    
    print "counter", counter


    db.commit()
    print "counter=total no of rows added to database", counter
    print "no of rows in table we began with", noOfRows
    #terminate the pool of multiprocessors
    pool.terminate()
    # close connection to db
    db.close()
    print "time taken to compute for pixels in this chunk and put in database =", time()-t3







    # #iterate over all pixels, process each one of them and write to database.
    # for res in iterator:
    #     #unpack parameters
    #     objID_, ra_, dec_, rErr_, dErr_ = res
    #     #increment loop counter
    #     noOfPixelsIterated += 1
    #     #so that you donot encounter empty pixels! 
    #     if (objID_.size != 0):
    #         dat1 = [(int(objID),float(ra), float(dec), float(rErr), float(dErr)) for objID, ra, dec, rErr, dErr in zip(objID_, ra_, dec_, rErr_, dErr_)]
    #         dat2.append(dat1)
    #         dat1 = []
    #         counter = counter + len(dat1)
    #         #print "len of dat=no of rows", len(dat1)
    #         if (noOfPixelsIterated % 1000 == 0):
    #             cursor.executemany('INSERT INTO finalData (objID, ra, dec, raError, decError) VALUES (?, ?, ?, ?, ?)', dat2)
    #             db.commit()
    #             print "committed now"
    #             #make dat empty again 
    #             dat2 = []
    #     if (objID_.size ==0):
    #         print "zero array"
    # #after the loop, if there are any more pixels left that are not executed because they aren't a  multiple of 1000
    # if (noOfPixelsIterated % 1000 != 0):
    #     cursor.executemany('INSERT INTO finalData (objID, ra, dec, raError, decError) VALUES (?, ?, ?, ?, ?)', dat2)
    #     db.commit()    

   
