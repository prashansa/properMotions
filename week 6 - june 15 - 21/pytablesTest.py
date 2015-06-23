import tables
import numpy as np


hd5file = tables.open_file('/disk1/bsesar/projects/PS1/proper_motion/prashansa/TriAnd_galaxies.h5')

triand = hd5file.root.triand

#print out columns in the table
triand.cols

#find unique object ids from the table
uniqueObjID = np.unique(triand.col('obj_id'))

#create an index into rows that have objID
ind = triand.col('obj_id') == -2961750073693619764

#get all rows using the index
#obj_rows is a normal numpy array 
obj_rows = triand[ind]


#close pytable files
h5file.close()


