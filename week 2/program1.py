import numpy as np
from astropy.time import Time
from esutil.coords import sphdist

# Variables : i specifies year, k specifies galaxy
# r = RA , d = DEC
# r_ik |for all k :  rBar = med(r_i) | for a given k : deltaR_ik = rBar - r_ik
# deltaRbar_i = med(deltaR_i)
# rNew_ik = r_
#
# convert mjd to years and bin years into four periods.
# searchRadius updateRadius 


#FUNCTIONS

#----------------------------------------------------------------
def calcMedianForOneDarray(oneDarray):
    med = np.median(oneDarray)
    return med

def calcResiduals(rBar, r):
    residual = []*len(r)
    #rBar is one value, r is in general an array
    for i in r:
        residual.append(rBar - i)
    return residual
#----------------------------------------------------------------



noOfYears = 4 
pixelWidth = 3 #arcminutes
searchRadius = 10 #arcminutes
#load data into a variable
testFile= np.load('test.npy')

#save unique galaxy ids in another variable
uniqueObjIDfile = np.unique(testFile['obj_id'])

i=0
j=0
k=0


#obtain residuals for all galaxies and for a given year
#rBar_k = [] 
deltaR_ik = [] 
for k in uniqueObjIDfile:
    temp = []*noOfYears
    for j in testFile:
        if (j['obj_id'] == k):
        #if testFile['obj_id'] == uniqueObjIDfile[k]:
            temp.append(j['ra'])
            i=i+1
        if (i==noOfYears): 
            break
    #print temp
    #rBar_k[i] = [ uniqueObjIDfile[k], calcMedianForOneDarray(temp)]
    #rBar_k.append([k, calcMedianForOneDarray(temp)])
    #deltaR_ik[k]=[uniqueObjIDfile[k] ,calcResiduals(rBar_k[k,k], temp)]
    #deltaR_ik.append([k, calcResiduals(calcMedianForOneDarray(temp), temp)])
    deltaR_ik.append(calcResiduals(calcMedianForOneDarray(temp), temp))
    break

#choose a ra/dec value to be the center of a pixel - can be anything, I choose the middle value of data file

center = uniqueObjIDfile[len(uniqueObjIDfile)/2]

#obtain ra and dec values separately
for i in testFile:
    if (i['obj_id']==center):
        raCenter  = i['ra']
        print raCenter
        decCenter = i['dec']
        print decCenter
       # break    
#find values of ra and dec for which this is the center of the pixel
raInitial = raCenter - (pixelWidth/2)
raFinal = raCenter + (pixelWidth/2)
decInitial = decCenter - (pixelWidth/2)
decFinal = decCenter + (pixelWidth/2)

k=0

#........................................................
#obtain the galaxies within this pixel- the searchArray  
#.......................................................

# updateThisArray = []*
 

#for row in rBar_k:
#     for i in row:
#     #i dont think i can use 'ra' as index anymore
#         #if (raInitial<= rBar_k['ra'] <= raFinal) :
#         if (raInitial<= rBar_k[i] <= raFinal) :
#             #updateThisArray[k]= row
#             #k=k+1
#             updateThisArray.append(row)
#     #callingFunction is the function that would calcuate the spherical distance between the center and galaxies, given ra and dec values.(i will add dec later here.
#     if (callingFunction(rBar_k[1],raCenter) <= (searchRadius*60)):
#         searchThisArray.append(row)

# sphericalDistArray = []*len(rBar_k)
# for row in rBar_k:
#     #rBar_k[1] gives the ra values
#     sphericalDistArray.append(callingFunction(rBar_k[1], raCenter))



#......................................................................
#obtain galaxies within a searchRadius of three times the pixelWidth - the updateArray
#you would need to take distance between the galaxies on the sky, so convert ra dec values to spherical coordinates - obtain the spherical distance - check if this is less than the searchRadius
#.......................................................................

searchThisArrayTemp = []

for i in testFile:
    angSep = sphdist(raCenter, decCenter, i['ra'],i['dec'])
    #print 'angSept is ', angSep
    if ((angSep*60)<3000):
    #if((angSep*60)<10):
        #print 'i within 10 is', i
        #break
        searchThisArrayTemp.append(i)
        #print searchThisArrayTemp
        #break
    
searchThisArray =  np.unique(searchThisArrayTemp[0])       

#contains all the angular separations
angSep =  sphdist(raCenter, decCenter, testFile['ra'],testFile['dec'])

#dtype=boolean
gal = angSep*60 < searchRadius

#ask how is this unique statement working!
#contains IDs of galxies with angualar separation less than the searchRadius
searchThisArrayIDs = np.unique(testFile['obj_id'][gal])

#-------------------------------------------------
#from esutil.coords import sphdist
# ang_sep = sphdist(ra1, dec1, ra2, dec2) all in degrees
# https://code.google.com/p/esutil/source/browse/trunk/esutil/coords.py
# description : sphdist: Calculate the arc length between two sets of points on the sphere. Currently only takes ra,dec.
# ------------------------------------------------
# dat= [] all the data is in this
# ang_sep =sphdist(ra0, dec0, dat['ra'], dat['dec']
# select galaxies within 10 arcminutes
# ind_galaxies = ang_sep*60 < 10
# object id's of galaxies within 10 arcmins
# objid_gal = np.unique(dat['obj_id'][ind_galaxies])
# ------------------------------------------------



deltaRbar_i = []*noOfYears
#first loop for noOfYears
#searchThisArrayRes = []
searchThisArrayResiduals = []

for i in range(noOfYears):
    for j in deltaR_ik:
        for k in searchThisArrayIDs:
            if (j[0]== k):
                searchThisArrayRes.append(j)
    #searchThisArrayResiduals = np.array(searchThisArrayRes)
    deltaRbar_i=np.median(searchThisArrayResiduals[:,i])
    #column = []#*len(deltaR_ik)
    #second loop for going over all galaxies
    #for i in deltaR_ik:
    #column.append(deltaR_ik[:,k+1])#k+1 because the first one is 'obj_id'
    #deltaRbar_i.append(calcMedianForOneDarray(column))


rNew_k = raCenter - deltaRbar_i
rNewBar_k= np.median(rNew_k)
            
#obtain delta ra_i bar for galaxies in searchArray, but obtain raNew_k for galaxies in updateArray
            


        


        

