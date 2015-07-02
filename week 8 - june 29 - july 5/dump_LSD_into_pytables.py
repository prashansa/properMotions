#!/usr/bin/env python
import numpy as np
import sys
import os
from lsd import DB
from lsd import bounds as lsdbounds
import tables
np.seterr(invalid='ignore', divide='ignore')

# connect to folders with LSD data
db = DB('/ssd-raid0/bsesar:/a41233d1/LSD/from_cfa')

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
    ra = tables.Float64Col(pos=1)
    dec = tables.Float64Col(pos=2)
    raErr = tables.Float32Col(pos=3)
    decErr = tables.Float32Col(pos=4)
    nObs = tables.UInt8Col(pos=5)
    mjd = tables.Float64Col(pos=6)
    gal = tables.UIntCol(pos=7)

# open a pytable file
h5file = tables.open_file("TriAnd.h5", mode = "w", title = "PS1 data")

# define compression filters
filters = tables.Filters(complib='blosc', complevel=5)

# create the table
table = h5file.create_table('/', 'triand_unsorted', Star, "TriAnd region", expectedrows=40563159, filters=filters)

star = table.row

# define selection bounds
bounds = lsdbounds.rectangle(95, -50, 165, -10, coordsys="gal")
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
