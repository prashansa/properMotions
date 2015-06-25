#processing a file with both stars and galaxies
#this code finds the residuals of galaxies for each epoch = ra/dec_epoch - ra/dec_final
#these residuals are used to update the ra/dec values for stars
#and then finally find ONE ra/dec value for each star for this catalogue

''' NOTE : OBJECT in general refers to a galaxy,  a star is referred explicitly as a STAR'''

'''MULTIPROCESSING + PYTABLES + SQLITE + LSD DATA'''

import numpy as np
from numpy.random import randn
from astropy.time import Time
from time import time
from esutil.coords import sphdist
import healpy as hp
from esutil.numpy_utils import match

from multiprocessing import Pool

import tables

import sys
import os

#http://research.majuric.org/trac/wiki/LargeSurveyDatabase
from lsd import DB
from lsd import bounds as lsdbounds

import sqlite3

#to see how floating point errors will be handled
#invalid = invalid floating point operations
#divide = division by zero
np.seterr(invalid='ignore', divide='ignore')

# connect to folders with LSD data
db = DB('/home/bsesar:/a41233d1/LSD/from_cfa')

# define the query
#keywords (case insensitive) select, from, where, into, as
query = 's.obj_id as obj_id, t.ra as ra, t.dec as dec, \
        t.raErr as raErr, t.decErr as decErr, \
        t.nObs as nObs, t.mjd as mjd, \
        -2.5*np.log10(s.mean(1)/s.mean_ap(1)) as sg_r, \
        -2.5*np.log10(s.mean(2)/s.mean_ap(2)) as sg_i \
        FROM dvo as t, ucal_fluxqy(matchedto=t, nmax=1, dmax=1.5) as s \
        WHERE (sg_r > 0.3) & (sg_i > 0.3) & \
        (sg_r < 1.0) & (sg_i < 1.0)'

#table definition
#TOCHECK - class Star(IsDescription): -- do we need "tables." everywhere? 
class Star(tables.IsDescription):
    #in LSD, obj_id is a 64-bit unsigned integer, but pytables cannot index 64-bit unsigned integers
    #so they are saved as 64-bit signed integers
    obj_id = tables.Int64Col(pos=0)
    ra = tables.Float64Col(pos=1)
    dec = tables.Float64Col(pos=2)
    raErr = tables.Float32Col(pos=3)
    decErr = tables.Float32Col(pos=4)
    nObs = tables.UInt8Col(pos=5)
    mjd = tables.Float64Col(pos=6)
    gal = tables.UIntCol(pos=7)

#open a pytable file
h5file = tables.open_file("TriAnd.h5", mode = "w", title = "PS1 data")

#define compression filters
filters = tables.Filters(complib='blosc', complevel=5)

#create the table
#create_table(where, name, obj, title, expectedrows, filters)
#"/" refers to h5file.root object
#
table = h5file.create_table('/', 'triand_unsorted', Star, "TriAnd region", expectedrows=40563159, filters=filters)

star = table.row

#define selection bounds
#gal long lower bound, gal lat lower bound, gal long upper bound, gal lat upper bound
bounds = lsdbounds.rectangle(95, -50, 165, -10, coordsys="gal") #(ra,dec) bottomleft; (ra,dec) topright
bounds = lsdbounds.make_canonical(bounds)

# query LSD for rows and store them into the pytable
dtype = table.colnames
for row in db.query(query).iterate(bounds=bounds):
    star['obj_id'] = (row[0]).astype('i8')
    for j in range(len(dtype)-2):
        star[dtype[j+1]] = row[j+1]
    if (row[7] > 0.3) & (row[7] < 1.0) & (row[8] > 0.3) & (row[8] < 1.0):
        star['gal'] = 1
    else:
        star['gal'] = 0
    star.append()

table.flush()

# create a full index on the obj_id column
indexrows = table.cols.obj_id.create_csindex(filters=tables.Filters(complib='blosc', complevel=5))

# create a new table that is sorted by obj_id
sI = table.cols.obj_id
table2 = tables.Table.copy(table, newname='triand', overwrite=True, sortby=sI, checkCSI=True, propindexes=True)

# delete the unsorted table
h5file.root.triand_unsorted.remove()

# close the pytable file
h5file.close()




################################################################################



def create_obs_table(ref, sc):
    """Creates and populate a table with observations for each SuperCosmos epoch"""

    # connect to the database
    db = sqlite3.connect('mydb')

    # DROP table OBS
    cursor = db.cursor()
    cursor.execute('''DROP TABLE IF EXISTS obs''')
    db.commit()

    # CREATE table obs
    cursor = db.cursor()
    cursor.execute('''
                   CREATE TABLE obs(obj_id INT NOT NULL,
                   ra REAL NOT NULL, decl REAL NOT NULL, raErr REAL NOT NULL,
                   declErr REAL NOT NULL, year REAL NOT NULL,
                   blend INT NOT NULL, quality INT NOT NULL)
                   ''')

    # commit changes
    db.commit()

    # INSERT data into table obs
    for year in np.unique(sc['epoch']):
        plate = np.abs(sc['epoch'] - year) < 0.00001
        sc_plate = sc[plate]
        m1, m2, d12 = htm_mesh.match(ref['ra'], ref['dec'], sc_plate['ra'], sc_plate['dec'], 2.0/3600., maxmatch=1)
        dat = [(int(obj_id), float(ra), float(dec), float(120.), float(120.), float(year), int(blend), int(quality)) for obj_id, ra, dec, year, blend, quality in zip(ref['obj_id'][m1], sc_plate['ra'][m2], sc_plate['dec'][m2], sc_plate['epoch'][m2], sc_plate['blend'][m2], sc_plate['quality'][m2])]
        cursor.executemany('INSERT INTO obs (obj_id, ra, decl, raErr, declErr, year, blend, quality) VALUES (?, ?, ?, ?, ?, ?, ?, ?)', dat)

    # commit changes
    db.commit()

    # create index on obj_id
    cursor.execute('''CREATE INDEX objid_idx on obs (obj_id)''')
    db.commit()

    # close connection to db
    db.close()

def create_ref_table(ref):
    """Create a table with reference PS1 positions"""

    # connect to the database
    db = sqlite3.connect('mydb')

    # DROP table OBS
    cursor = db.cursor()
    cursor.execute('''DROP TABLE IF EXISTS ref''')
    db.commit()

    # CREATE table ref
    cursor = db.cursor()
    cursor.execute('''
                   CREATE TABLE ref(obj_id INT PRIMARY KEY NOT NULL,
                   ra REAL NOT NULL, decl REAL NOT NULL, r REAL NOT NULL,
                   galaxy INT NOT NULL)
                   ''')

    # commit changes
    db.commit()

    # INSERT data into table ref
    galaxy = np.zeros(ref.size)

    galaxy[(ref['sg_r'] > 0.3) & (ref['sg_i'] > 0.3) & (ref['sg_r'] < 1) & (ref['sg_i'] < 1)] = 1
    dat = [(int(obj_id), float(ra), float(dec), float(r), int(gal)) for obj_id, ra, dec, r, gal in zip(ref['obj_id'], ref['ra'], ref['dec'], ref['r'], galaxy)]

    cursor.executemany('INSERT INTO ref (obj_id, ra, decl, r, galaxy) VALUES (?, ?, ?, ?, ?)', dat)
    db.commit()

    # close connection to db
    db.close()




################################################################################


''' GLOBAL TASKS'''

#variables
searchRadius = 10 #in arcminutes

#load h5 file -- note the 's' in 'tableS.open'
#table is already sorted by object id
#this file contains stars and galaxies in the TriAnd region
#column called "gal" : gal = 1 for galaxies, gal = 0 for stars
t0 = time()
testH5file= tables.open_file('/home/bsesar/projects/PS1/proper_motion/prashansa/TriAnd.h5')
t1= time()
print t1 - t0

#load the table
triand = testH5file.root.triand

#select galaxies
#mask = triand.col('gal')==1
#galaxy numpy array
#galaxyIDfile = triand[mask]

#where function works much faster than the masking option
#conditionForGal = 'gal == 1'
#conditionForStar= 'gal == 0'

maskForGal= triand.get_where_list('gal == 1')
maskForStar = triand.get_where_list('gal == 0')

#select galaxy ids 
passValuesGal = triand.col('obj_id')[maskForGal]
#save the number of galaxies
noOfGal = len(passValuesGal)
#save unique galaxy ids in another variable
uniqueObjIDfile = np.unique(passValuesGal)
#the number of observations actually available for each object
noOfObs = np.zeros(len(uniqueObjIDfile)) 

#select star ids 
passValuesStar = triand.col('obj_id')[maskForStar]
#save the number of stars
noOfStar = len(passValuesStar)
#save unique star ids in another variable
uniqueStarIDfile = np.unique(passValuesStar)


#binning mjds for once
mjdValues = triand.col('mjd')[maskForGal]
mjdSorted = np.sort(mjdValues)
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


# ind = triand.get_where_list(conditionForGal)
#raValues1 = []
#decValue1 = []
#for row in triand.where(conditionForGal):
#    raValues1.append(row['ra'])
#    decValue1.append(row['dec'])
#raValues1 = [row['ra'] for  in triand.where(conditionFoconditionForGalrGal)]
#decValues1 = [row['dec'] for row in triand.where(conditionForGal)]
#raValues = np.array(raValues1)
#decValues = np.array(raValues1)

#store ra/dec values of galaxies
raValues = triand.col('ra')[maskForGal]
decValues = triand.col('dec')[maskForGal]
#store ra/dec values of stars
raValuesStar = triand.col('ra')[maskForStar]
decValuesStar = triand.col('dec')[maskForStar]


#calculate phi and theta for ra/dec values of the testH5file
phiForObj = (raValues * np.pi)/180 
thetaForObj = (90 - decValues)* (np.pi/180)
#calculate pixel indices for each theta and phi in the dataset
pixelIndexForObjTest= hp.ang2pix(nside,thetaForObj, phiForObj)
#we want every pixel index only once! 
pixelIndexForObj = np.unique(pixelIndexForObjTest)

#pack parameters for workers
#Note: even if res_ra and res_dec are huge, they are passed as references, so "parameterListForPixel" does not use a lot of memory
parameterListForPixel = [(pickPixelNo, pixelRa[pickPixelNo], pixelDec[pickPixelNo], noOfObs, objIDarray, medianRaArray, medianDecArray, residualRaArray, residualDecArray) for pickPixelNo in pixelIndexForObj]

#start workers
pool = Pool(processes=8)
ti = time()
#chunk size is used to submit jobs in batches which reduces overhead 
iterator = pool.imap_unordered(pixelTasks, parameterListForPixel, chunksize=100)
ti2 = time()-ti

def pixelTasks(parameterListForPixel):
    """Workers will execute this function."""
    #unpack the input parameter list that is to say separate them into different arrays
    pixelNo, ra, dec, nObs, objID, medRa, medDec, resRa, resDec = parameterListForPixel
    #define an index so that all the calculations can be put in place later
    #ind =
    global searchRadius
    #find objects within the searchRadius and whose noOfObs is atleast 3
    angSepMask = sphdist(ra, dec, medRa, medDec) <= (searchRadius*60) #& (nObs>=3.0) -- not needed since this has been taken care of while calculating residuals and medians previously 
    global triand
    #separate these entire rows from the original testH5file
    #we donot need to check for galaxies now, since the parameters being referred to have already been checked

    #searchFile = triand[angSepMask]
    objInRadius = objID[angSepMask]
    
    k = triand[triand.get_where_list('obj_id == o') for o in objInRadius]
    #this might be stupid
    searchFile = []
    for o in objInRadius:
        searchFile.append(triand[triand.get_where_list('obj_id == o')])

    m1, m2 = match(triand.col('obj_id'), objInRadius)
    print triand.col('obj_id')[m1]
    print objInRadius[m2]


    

    #putting unique object IDs in this file
    uniqueSearchFile = np.unique(searchFile['obj_id'])
    #file containing data for objects inside pixel
    objInPixelFile =searchFile[pixelIndexForObj == pixelNo]
    #i make a copy of the file, since i need the original values later
    objInPixelFile2 = searchFile[pixelIndexForObj == pixelNo]
    #file containing unique objIDs of the pixel
    uniqueObjInPixelFile = np.unique(objInPixelFile['obj_id'])
    #obtain original ra and dec values for the galaxies in the pixel
    #originalRaGalInPixel  =  raValues[pixelIndexForObj == pixelNo]
    #originalDecGalInPixel =  decValues[pixelIndexForObj == pixelNo]
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
    #arrays to store deltas for stars
    deltaRaArray  = np.zeros()
    deltaDecArray = np.zeros()
    #take the current values separately
    currObj = objInPixelFile['obj_id'][0]
    raObj = [objInPixelFile['ra'][0]]
    decObj = [objInPixelFile['dec'][0]]
    mjdObj = [objInPixelFile['mjd'][0]]
    raObj2 = [objInPixelFile2['ra'][0]]
    decObj2 = [objInPixelFile2['dec'][0]]
    mjdObj2 = [objInPixelFile2['mjd'][0]]
    #set variables
    counter = 0
    pos2 = 0 
    #t2 = time()
    for i in np.arange(objInPixelFile.size-1)+1:#np.arange(10)+1: #
        objID = objInPixelFile['obj_id'][i]
        ra = objInPixelFile['ra'][i]
        dec = objInPixelFile['dec'][i]
        objID2 = objInPixelFile2['obj_id'][i]
        ra2 = objInPixelFile2['ra'][i]
        dec2 = objInPixelFile2['dec'][i]
        if(objID == currObj):
            raObj.append(ra)
            decObj.append(dec)
            raObj2.append(ra2)
            decObj2.append(dec2)
        else:
            #these contain the averaged/processed ra and dec values of galaxies
            raObj=np.array(raObj)
            decObj=np.array(decObj)
            #these contain the original ra and dec values of galaxies
            raObj2=np.array(raObj2)
            decObj2=np.array(decObj2)
            #find the final ra and dec values
            finalRaArray[counter] = np.median(raObj)
            finalDecArray[counter] = np.median(decObj)
            #calculate rms values for ra/dec
            rmsRa = 0.741*(np.percentile(raObj, 0.75) - np.percentile(raObj, 0.25))
            rmsDec = 0.741*(np.percentile(decObj, 0.75) - np.percentile(decObj, 0.25))
            #calculate uncertainity in median coordinates
            medianRaError  = np.sqrt((np.pi/2)/(rmsRa.size-1))*rmsRa
            medianDecError = np.sqrt((np.pi/2)/(rmsDec.size-1))*rmsDec
            # print finalRaArray[counter]
            #print i
            #print counter
            #calculate how much the galaxy moved in each epoch
            #recall : the final ra and dec values were obtained after removing offsets epochwise and obj_id wise
            #thus implying that these values make the galaxies static as required
            deltaRa  = raObj2 - finalRaArray[counter]
            deltaDec = decObj2 - finalDecArray[counter]
            deltaRaArray[pos2:i]  = deltaRa
            deltaDecArray[pos2:i] = deltaDec
            #move to the next obj
            currObj=objID
            pos2 = i 
            counter +=1
            #make them lists again
            raObj=[objInPixelFile['ra'][i]]
            decObj=[objInPixelFile['dec'][i]]
            raObj2=[objInPixelFile2['ra'][i]]
            decObj2=[objInPixelFile2['dec'][i]]
    #t3 = time()
    #print 'time taken for calculating ra and dec values for a pixel is ', t3-t2

    #final updation for STARS now
    if (len(mjdBreakAt)>=3.0):
        previousVar2 = np.min(mjdSorted)
        for var1, var2 in enumerate(mjdBreakAt, start=0):
            mjdIndexInPixel = (objInPixelFile['mjd']<var2) & (previousVar2 < objInPixelFile['mjd'])
            if any(mjdIndex):
                medianDeltaRaEpochWise[var1] = np.median(deltaRaArray[mjdIndexInPixel])
                medianDeltaDecEpochWise[var1] = np.median(deltaDecArray[mjdIndexInPixel])
                #replacing old values with new values
                starInPixelFile['ra'][mjdIndexInPixel] -=  medianDeltaRaEpochWise[var1]
                starInPixelFile['dec'][mjdIndexInPixel] -=  medianDeltaDecEpochWise[var1]
            else:
                medianDeltaRaEpochWise[var1+1] = np.median(deltaRaArray[mjdIndexInPixel])
                medianDeltaDecEpochWise[var1+1] = np.median(deltaDecArray[mjdIndexInPixel])
                #replacing old values with new values
                starInPixelFile['ra'][mjdIndexInPixel] -=  medianDeltaRaEpochWise[var1+1]
                starInPixelFile['dec'][mjdIndexInPixel] -=  medianDeltaDecEpochWise[var1+1]
            previousVar2 = var2


    
    return finalRaArray, finalDecArray, medianRaError, medianDecError 



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





def calcMedianAndResiduals(testH5file, uniqueObjIDfile):
    '''calculates median ra and dec along with residuals for a region of sky once'''
    #stores the median ra and dec values for each object
    medianRaArray  = np.zeros(uniqueObjIDfile.size)
    medianDecArray = np.zeros(uniqueObjIDfile.size)
    objIDarray = np.zeros(uniqueObjIDfile.size) #, dtype='u8')
    #stores the residual values for each object and for each epoch
    #residualRaArray= np.zeros(len(triand[maskForGal]))
    #residualDecArray= np.zeros(len(triand[maskForGal]))
    residualRaArray  = np.zeros(noOfGal)
    residualDecArray = np.zeros(noOfGal)
    #store the first object values in here
    currentObjID = triand.col('obj_id')[maskForGal][0]
    #instead of taking the entire thing, take only the zeroth for mask, and then the corresponding ra etc
    #raObj = [triand.col('ra') [maskForGal][0]]
    raObj  = [triand.col('ra')[maskForGal[0]]]
    decObj = [triand.col('dec')[maskForGal[0]]]
    mjdObj = [triand.col('mjd')[maskForGal[0]]]
    #set variables
    objIDcounter = 0
    pos=0
    #removed these, since they will occupy much space in memory
    #instead put in the i th component of the mask in the loop itself! CLEVER BRANI! 
    objIDgal = triand.col('obj_id')[maskForGal]
    raGal    = triand.col('ra')[maskForGal]
    decGal   = triand.col('dec')[maskForGal]
    mjdGal   = triand.col('mjd')[maskForGal]
    t1 = time()
    #it should start from the second row, so add one.
    for i in np.arange(30)+1:# np.arange(noOfGal-1)+1:#
        #print i
        #objID = triand.col('obj_id')[maskForGal[i]]
        #ra    = triand.col('ra')[maskForGal[i]]
        #dec   = triand.col('dec')[maskForGal[i]]
        #mjd   = triand.col('mjd')[maskForGal[i]]
        objID = objIDgal[i]
        ra    = raGal[i]
        dec   = decGal[i]
        mjd   = mjdGal[i]
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
                #also need to store the objID
                objIDarray[objIDcounter]=objID
                #calculate residual
                residualRa = raObj - medianRa
                residualDec = decObj - medianDec
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
            #raObj = [triand.col('ra')[maskForGal[i]]]
            #decObj = [triand.col('dec')[maskForGal[i]]]
            #mjdObj = [triand.col('mjd')[maskForGal[i]]]
            raObj  = [raGal[i]]
            decObj = [decGal[i]]
            mjdObj = [mjdGal[i]]
    print 'total time taken is',time() - t1
    return medianRaArray, medianDecArray, residualRaArray, residualDecArray, objIDarray





#calculate rms and median error
medianRaError  = np.zeros(uniqueObjIDfile.size)
medianDecError = np.zeros(uniqueObjIDfile.size)

rmsRa = 0.741*(np.percentile(raObj, 0.75) - np.percentile(raObj, 0.25))
rmsDec = 0.741*(np.percentile(decObj, 0.75) - np.percentile(decObj, 0.25))

medianRaError  = np.sqrt((np.pi/2)/(rmsRa.size-1))*rmsRa
medianDecError = np.sqrt((np.pi/2)/(rmsDec.size-1))*rmsDec
