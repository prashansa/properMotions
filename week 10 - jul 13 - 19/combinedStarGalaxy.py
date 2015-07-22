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
import combinedStarGalaxyFunctions

#####################################################################


np.seterr(invalid='ignore', divide='ignore')

if os.getenv('HOSTNAME') == 'aida41147':
    db = DB('/home/bsesar/projects/PS1/DVO:/a41233d1/LSD/from_cfa')
else:
    db = DB('/a41147d1/bsesar/projects/PS1/DVO:/a41233d1/LSD/from_cfa')


searchRadius = 10.0 #in arcminutes
nside = 2**10
#specifying size of chunks, assume they are squares so both ra/dec have same breadth
chunkSize = 10.0 #degrees
#for future - once you have implemented this, you might also want to keep the no of obj same in each chunk of the sky -- so your chunksize would depend on that!
someBigNumber = 1000


chunkNo = 9999 #arbitrary
raMin  = 137.483996347
decMin = 23.6666666667
raMax  = 148.880976675
decMax = 33.6666666667

#h5fileName,noOfRows, countForRows = combinedStarGalaxyFunctions.executeChunkWise(raMin, decMin, raMax, decMax, chunkNo)
#noOfRows = 3852132


t0 = time()

testH5file = tables.open_file('/home/gupta/projects/propermotion/combinedFiles/combinedFile9999.h5')#%s.h5' %h5fileName)

#if (noOfRows!= countForRows): 
#    "error : the no. of rows dont match for chunk no", chunkNo, "with bounds", raMin, decMin, raMax, decMax 
#    break

# load the table
table = testH5file.root.aTable
t1= time()
print "time taken to load h5file and table", t1 - t0
t2 = time()

maskForGal  = table.get_where_list('gal==1')
maskForStar = table.get_where_list('gal==0')

#get IDs corresponding to all the detections of galaxies
galaxyIDs = table.col('obj_id')[maskForGal]
raGalaxy   = table.col('ra')[maskForGal]
decGalaxy  = table.col('dec')[maskForGal]
mjdsGalaxy = table.col('mjd')[maskForGal]

#get IDs corresponding to all the detections of stars
starIDs  = table.col('obj_id')[maskForStar]
raStar   = table.col('ra')[maskForStar]
decStar  = table.col('dec')[maskForStar]
mjdsStar = table.col('mjd')[maskForStar]


# mjdValues  = table.col('mjd')[maskForGal]
# #len(mjdValues) = 1190749
# mjdSorted  = np.sort(mjdValues)
# deltaT     = mjdSorted[0:-1] - mjdSorted[1:]
# tBreakAt   = np.where(deltaT < (-100.0))[0]
# mjdBreakAt = mjdSorted[tBreakAt + 1]
# #array([ 55521.60506944,  55860.6159838 ,  56231.6315625 ])


mjdValuesStar  = table.col('mjd')[maskForStar]
mjdSortedStar  = np.sort(mjdValuesStar)
#len(mjdSortedStar) = 2661383
deltaTStar     = mjdSortedStar[0:-1] - mjdSortedStar[1:]
tBreakAtStar   = np.where(deltaTStar < (-100.0))[0]
mjdBreakAtStar = mjdSortedStar[tBreakAtStar + 1]
#array([ 55521.60238426,  55897.66508102,  56231.6315625 ])


#store min and max values separately
minRaGal  = table.col('ra')[maskForGal].min()
maxRaGal  = table.col('ra')[maskForGal].max()
minDecGal = table.col('dec')[maskForGal].min()
maxDecGal = table.col('dec')[maskForGal].max()



##################################################
#specify the database name
databaseName = 'mydb86'
#establish connection
connection = sqlite3.connect('%s'%databaseName)
#obtain cursor object
cursor = connection.cursor()
#obtain object_ids, finalRa and finalDec values from the database and put in separate arrays
rowCounter = 0 
#define arrays with a reasonable size
objIDarray = np.zeros(100000000, dtype = 'u8') # eight zeros
finalRaArray  = np.zeros(100000000)
finalDecArray = np.zeros(100000000)

#loop over all rows in the table and obtain values
for row in cursor.execute('SELECT * FROM finalData'):
    objIDarray[rowCounter]     = row[0]
    finalRaArray[rowCounter]  = row[1]
    finalDecArray[rowCounter]  = row[2]
    rowCounter +=1

print "rowCounter", rowCounter
#close the connection to database
connection.close()

#trim the array to the required size 
objIDarray = objIDarray[0:rowCounter]
finalRaArray = finalRaArray[0:rowCounter]
finalDecArray = finalDecArray[0:rowCounter]
######################################################


#pixelise using data from available galaxies. 

phiForObj   = (finalRaArray*np.pi)/180 
thetaForObj = (90 - finalDecArray)* (np.pi/180) 
pixelIndexForObj = hp.ang2pix(nside, thetaForObj, phiForObj) 

# goodPixels = (finalRaArray >= (minRaGal + 10./60)) & \
#     (finalRaArray <= (maxRaGal - 10./60)) & \
#     (finalDecArray >= (minDecGal + 10./60)) & \
#     (finalDecArray <= (maxDecGal - 10./60))
    
pixelIndexArray = np.unique(pixelIndexForObj)#[goodPixels])
theta, phi = hp.pix2ang(nside, pixelIndexArray)
pixelRa  = 180*phi/np.pi
pixelDec = 90 - theta*180/np.pi

print "Going to process %d pixels." % pixelIndexArray.size
#Going to process 28878 pixels -- the same as we obtained while running the previous galaxy code -- expected since the min ra/dec are the same


phiForStar   = (raStar*np.pi)/180 
thetaForStar = (90 - decStar)* (np.pi/180) 
pixelIndexForStar = hp.ang2pix(nside, thetaForStar, phiForStar) 


packParameterList = [ (pickPixelNo, pixelRa[index], pixelDec[index],objIDarray, finalRaArray, finalDecArray, galaxyIDs, raGalaxy, decGalaxy,mjdsGalaxy , starIDs, mjdsStar, raStar, decStar, pixelIndexForStar,mjdSortedStar, mjdBreakAtStar ) for index, pickPixelNo in enumerate(pixelIndexArray)]


# #start workers
# pool = Pool(processes=24)
# ti = time()
# #chunk size is used to submit jobs in batches which reduces overhead
# iterator = pool.imap_unordered(combinedStarGalaxyFunctions.pixelTasksCombinedData,packParameterList, chunksize=100)

counter = 0
noOfPixelsIterated = 0
dat1 = []
dat2 = []
dat  = []

#for res in iterator:
for param in packParameterList:
    noOfPixelsIterated +=1
    objIDstar, raStar, decStar =combinedStarGalaxyFunctions.pixelTasksCombinedData(param)
    #dat.append(objIDstar)
    #dat1.append(raStar)
    #dat2.append(decStar)
    #print 'obj id ',objIDstar#dat
    #print 'ra', raStar# dat1
    #print 'dec', decStar # dat2
    print 'noOfPixelsIterated',noOfPixelsIterated
    
