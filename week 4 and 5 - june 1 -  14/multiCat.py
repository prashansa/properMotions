import numpy as np
from multiprocessing import Pool
from numpy.random import randn


searchRadius = 10 #in arcminutes


def pixelTasks(pickPixelNo, ra, dec, residualRaArray, residualDecArray, noOfObs, medianRaArray, medianDecArray, testFile):
    """Workers will execute this function."""
    
    
