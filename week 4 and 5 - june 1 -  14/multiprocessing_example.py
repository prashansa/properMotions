#!/usr/bin/env python
import numpy as np
from multiprocessing import Pool
from numpy.random import randn

def work(param):
    """Workers will execute this function."""
    # unpack the input parameters
    pixel_id, ind, res_ra, res_dec = params
    # select all rows in this pixel
    ind = ...
    ## now do something
    # for example, add Gaussian noise to certain rows (indexed by ind array)
    res_ra_aux = res_ra[ind] + 0.5*randn(res_ra[ind].size)
    res_dec_aux = res_dec[ind] + 0.5*randn(res_dec[ind].size)
    # return the result *and* the index ind
    return ind, res_ra_aux, res_dec_aux

# load res_ra and res_dec
res_ra = ...
res_dec = ...

# need the list of pixels as well
pixel_id_list = ...

# pack parameters for workers
# Note: even if res_ra and res_dec are huge, they are passed as references,
## so "params" does not use a lot of memory
params = [(pixel_id, res_ra, res_dec) for pixel_id in pixel_id_list]

# start workers
pool = Pool(processes=8)

# chunk_size is used to submit jobs in batches --> this reduces overhead
it = pool.imap_unordered(work, params, chunksize=100)

# create storage arrays
ra_final = ...
dec_final = ...

# loop over parameter sets
# "it" is an iterator that returns the result in the "res" list
for res in it:
    ra_final[res[0]] = res[1]
    dec_final[res[0]] = res[2]

## in the above loop, a worker takes a set of (pixel_id, res_ra, res_dec)
## parameters, executes the "work" function using these parameters,
## and returns the results (res_ra_aux and res_dec_aux which are returned in
## res[1] and res[2]) *and* the index ind (which is returned in res[0]) that
## can be used to store these results into the correct rows in storage arrays

# terminate workers when done
pool.terminate()

# dump data
np.save('ra_final.npy', ra_final)
np.save('dec_final.npy', dec_final)
