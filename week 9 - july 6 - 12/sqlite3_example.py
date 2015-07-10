import sqlite3

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
