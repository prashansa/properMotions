import numpy as np
from astropy.time import Time
from time import time
from esutil.coords import sphdist
import healpy as hp




noOfYears = 5
pixelWidth = 3 #arcminutes
searchRadius = 10 #arcminutes

NSIDE = 2**10 # with this NSIDE, I oversample the correction map by a factor of 5
tempArray = np.arange(hp.nside2npix(NSIDE))
theta, phi = hp.pix2ang(NSIDE, tempArray)
pixelRa = 180*phi/np.pi
pixelDec = 90 - theta*180/np.pi


#take galaxies within searchRadius
index = sphdist(testFile['ra'], testFile['dec'] <= searchRadius)

useThisForUpdate = testFile[index]

#count the number of observations
#noOfObs = 0 #before the loop
#noOfObs +=1 # if it goes to the else part
