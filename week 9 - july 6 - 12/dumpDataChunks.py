import numpy as np
import sys
import os
from lsd import DB
from lsd import bounds as lsdbounds
import tables

#to see how floating point errors will be handled
#invalid = invalid floating point operations
#divide = division by zero
np.seterr(invalid='ignore', divide='ignore')

# connect to folders with LSD data
#keywords (case insensitive) select, from, where, into, as
#db = DB('/ssd-raid0/bsesar:/a41233d1/LSD/from_cfa')
#db = DB('/ssd-raid0/bsesar/projects/PS1/DVO:/a41233d1/LSD/from_cfa')
db = DB('/home/bsesar/projects/PS1/DVO:/a41233d1/LSD/from_cfa')

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
#TOCHECK - class Star(IsDescription): -- do we need "tables." everywhere? yes
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

#specifying size of chunks, assume they are squares so both ra/dec have same breadth
chunkSize = 10.0 #degrees
#for future - once you have implemented this, you might also want to keep the no of obj same in each chunk of the sky -- so your chunksize would depend on that!

someBigNumber = 1000000

raMin  = 0.0 #ra.min()#
decMin = -90.0 #dec.min()#
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
currentRaMax  = raMin + chunkSize
currentDecMax = decMin + chunkSize

raMinArray[0] = raMin
decMinArray[0]= decMin
raMaxArray[0] = currentRaMax
decMaxArray[0]= currentDecMax

chunkNo = 1

while(currentDecMax < decMax):
    decAvg = (currentDecMin + currentDecMax)/2.0
    print "decAvg", decAvg
    print "cosine of decAvg", np.cos(np.radians(decAvg))
    while(currentRaMax < raMax):
        raMinArray[chunkNo] = (currentRaMax - (40.0/60.0))#*np.cos(np.radians(decAvg))
        raMaxArray[chunkNo] = raMinArray[chunkNo] + chunkSize
        decMinArray[chunkNo]= currentDecMin
        decMaxArray[chunkNo]= currentDecMax
        if (raMaxArray[chunkNo] > raMax):
            raMaxArray[chunkNo] = raMax
        if (decMaxArray[chunkNo] > decMax):
            decMaxArray[chunkNo] = decMax
        currentRaMax = raMaxArray[chunkNo]
        chunkNo += 1
    currentRaMin  = raMin
    currentRaMax  = raMin + chunkSize
    currentDecMin = currentDecMax - (40.0/60.0)
    currentDecMax = currentDecMin + chunkSize
    print "chunk no", chunkNo  


#reversing the loop, sweep over dec for each ra
# while(currentRaMax < raMax):
#      while(currentDecMax < decMax):
#         decAvg = (currentDecMin + currentDecMax)/2.0
#         decMinArray[chunkNo] = currentDecMax - (40.0/60.0)
#         decMaxArray[chunkNo] = decMinArray[chunkNo] + chunkSize
#         raMinArray[chunkNo]  = currentRaMin
#         raMaxArray[chunkNo]  = currentRaMax
#         currentDecMax = decMaxArray[chunkNo]
#         chunkNo += 1
#     currentDecMin = decMin
#     currentDecMax = decMax
#     currentRaMin  = (currentRaMax - (40.0/60.0))*np.cos(np.radians(decAvg))
       
    



    

def makeChunks(raMin, decMin, raMax, decMax, chunksize, fileName):
    #this function will make chunks of the sky area passed onto it, return ra/dec bounds of these chunks
    
      




   
# open a pytable file
h5file = tables.open_file("file1.h5", mode = "w", title = "try1")

# define compression filters
filters = tables.Filters(complib='blosc', complevel=5)

# create the table
#create_table(where, name, obj, title, expectedrows, filters)
#"/" refers to h5file.root object
table = h5file.create_table('/', 'table1', Star, "tryTable", expectedrows=40563159, filters=filters)

star = table.row

# define selection bounds
#gal long lower bound, gal lat lower bound, gal long upper bound, gal lat upper bound
bounds = lsdbounds.rectangle(raMin, decMin, raMax, decMax, coordsys="gal")#(ra,dec) bottomleft; (ra,dec) topright
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
table2 = tables.Table.copy(table, newname='haha', overwrite=True, sortby=sI, checkCSI=True, propindexes=True)

# delete the unsorted table
h5file.root.file1.remove()

# close the pytable file
h5file.close()












'''Fake data -- testing '''

#always write floats! it was taking 20/60 as zero!! :(

#specifying size of chunks, assume they are squares so both ra/dec have same breadth
chunkSize = 10.0 #degrees

#we take so many points, since at this value the uniform distribution seems nicely uniform!
ra  = np.random.rand(1000000)*360.0
dec = np.random.rand(1000000)*180.0-90.0

raMin  = ra.min()
decMin = dec.min()
raMax  = ra.max()
decMax = dec.max()

#write a loop that calls makeChunks

#define bounds of the first chunk
ra0  = raMin
ra1  = raMin + chunkSize
dec0 = decMin#+ 1.0 #40.0 
dec1 = dec0 + chunkSize

ra2  = ra1 - 20.0/60.0
ra3  = ra2 + chunkSize
dec2 = dec0
dec3 = dec1

ra0 = 10.2
ra1 = 10.7

dec0 = 31.8
dec1 = 32.0



raChunk1  = raWith2Obs[(raWith2Obs <= ra1) &  (raWith2Obs >= ra0)]
decChunk1 = decWith2Obs[(decWith2Obs <= dec1) & (decWith2Obs >= dec0)]

raChunk1_3obs  = raWith3Obs[(raWith3Obs <= ra1) &  (raWith3Obs >= ra0)]
decChunk1_3obs = decWith3Obs[(decWith3Obs <= dec1) & (decWith3Obs >= dec0)]


#STUPID ! PUT THE LOGICAL EXPRESSIONS TOGETHER! 
#define chunk arrays
raChunk1  = ra[(ra <= ra1) &  (ra >= ra0)]
decChunk1 = dec[(dec <= dec1) & (dec >= dec0)]

raChunk2  = ra[(ra <= ra3) &  (ra >= ra2)]
decChunk2 = dec[(dec <= dec3) & (dec >= dec2)]

#make the ra/dec arrays the same length -- YOU DONOT NEED TO DO THIS
if (len(raChunk1_3obs) < len(decChunk1_3obs)):
    decChunk1_3obs = decChunk1_3obs[0:len(raChunk1_3obs)]
if (len(decChunk1_3obs) < len(raChunk1_3obs)):
    raChunk1_3obs = raChunk1_3obs[0:len(decChunk1_3obs)]
if (len(raChunk2) < len(decChunk2)):
    decChunk2 = decChunk2[0:len(raChunk2)]
if (len(decChunk2) < len(raChunk2)):
    raChunk2 = raChunk2[0:len(decChunk2)]


plt.clf(), plt.plot(raChunk1_3obs, decChunk1_3obs, 'ro',  alpha=0.5)
plt.plot(raChunk2, decChunk2, 'bo', alpha=0.5)

plt.plot(raChunk1, decChunk1, 'ro', ms=0.5, alpha=0.5)
plt.plot(raChunk2, decChunk2, 'bo', ms=0.5, alpha=0.5)



plt.plot(raChunk1, decChunk1, raChunk2, decChunk2, alpha=0.6)
    

plt.plot(raChunk1, decChunk1),  plt.plot(raChunk2, decChunk2)

































    
