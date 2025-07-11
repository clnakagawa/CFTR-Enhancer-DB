library(DBI)

# setup tables
db <- dbConnect(RSQLite::SQLite(), "enhancerDB.sqlite")

# drop tables
dbExecute(db, 'DROP TABLE IF EXISTS Activity')
dbExecute(db, 'DROP TABLE IF EXISTS Target')
dbExecute(db, 'DROP TABLE IF EXISTS Contact')
dbExecute(db, 'DROP TABLE IF EXISTS Source')
dbExecute(db, 'DROP TABLE IF EXISTS Factor')
dbExecute(db, 'DROP TABLE IF EXISTS Cell')
dbExecute(db, 'DROP TABLE IF EXISTS Region')

# make tables
# table of regions
dbExecute(db, 'CREATE TABLE Region (
          rid int PRIMARY KEY,
          chr varchar(5),
          start int,
          end int,
          type varchar(50),
          name varchar(50))')

# table of cell types
dbExecute(db, 'CREATE TABLE Cell (
          cid int PRIMARY KEY,
          name varchar(50),
          tissue varchar(50))')

# table of transcription factor (+ hist modifications?)
dbExecute(db, 'CREATE TABLE Factor (
          fid int PRIMARY KEY,
          name varchar(50))')

# table of sources for different associations
dbExecute(db, 'CREATE TABLE Source (
          sid int PRIMARY KEY,
          doi varchar(50))')

# table of region activity in cell type 
dbExecute(db, 'CREATE TABLE Activity (
          rid int,
          cid int,
          sid int,
          evidence varchar(50),
          activity varchar(20),
          PRIMARY KEY (rid, cid, sid),
          FOREIGN KEY (rid) REFERENCES Region(rid),
          FOREIGN KEY (cid) REFERENCES Cell(cid),
          FOREIGN KEY (sid) REFERENCES Source(sid))')

# table of region activity in cell type 
dbExecute(db, 'CREATE TABLE Target (
          rid int,
          cid int,
          fid int,
          sid int,
          evidence varchar(50),
          activity varchar(20),
          PRIMARY KEY (rid, cid, fid, sid),
          FOREIGN KEY (rid) REFERENCES Region(rid),
          FOREIGN KEY (cid) REFERENCES Cell(cid),
          FOREIGN KEY (fid) REFERENCES Factor(fid),
          FOREIGN KEY (sid) REFERENCES Source(sid))')

# table of region contacts
dbExecute(db, 'CREATE TABLE Contact (
          rid1 int,
          rid2 int,
          cid int,
          sid int,
          evidence varchar(50),
          activity varchar(20),
          PRIMARY KEY (rid1, rid2, cid, sid),
          FOREIGN KEY (rid1) REFERENCES Region(rid),
          FOREIGN KEY (rid2) REFERENCES Region(rid),
          FOREIGN KEY (cid) REFERENCES Cell(cid),
          FOREIGN KEY (sid) REFERENCES Source(sid))')

dbDisconnect(db)