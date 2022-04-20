import os
import os.path
import sys
import glob
import psycopg
import types
from nexoclom.utilities.database_connect import NexoclomConfig


DEFAULT_MESSENGER_DB = 'messengeruvvsdb'

def messengerdb_connect(configfile=None):
    config = NexoclomConfig(configfile=configfile)
    messengerdb = config.__dict__.get('mesdatabase', DEFAULT_MESSENGER_DB)
    
    if config.dbhost:
        con = psycopg.connect(host=config.dbhost, dbname=messengerdb,
                              port=config.port)
    else:
        con = psycopg.connect(dbname=messengerdb, port=config.port)
        
    con.autocommit = True
    return con


# def messenger_database_setup(configfile=None):
#     """Setup the database from SQL database dump files.
#     Repopulates the database using a SQL backup rather than the original
#     IDL save files. See :doc:`database_fields` for a description of the
#     tables and fields used by MESSENGERuvvs.
#
#     **Parameters**
#
#     force
#         If True, deletes old database tables and remakes them.
#         Default is False, which only creates the tables if necessary.
#
#     **Returns**
#
#     No output.
#
#     """
#     # Create MESSENGER database if necessary
#     config = NexoclomConfig()
#     messengerdb = config.__dict__.get('mesdatabase', DEFAULT_MESSENGER_DB)
#     if 'mesdatapath' in config.__dict__:
#         datapath = config.mesdatapath
#     else:
#         raise ConfigfileError(config.configfile, 'mesdatapath')
#
#     with psycopg.connect(port=config.port) as con:
#         con.autocommit = True
#         cur = con.cursor()
#         cur.execute('select datname from pg_database;')
#         dbs = [r[0] for r in cur.fetchall()]
#
#         if messengerdb not in dbs:
#             print(f'Creating database {messengerdb}')
#             cur.execute(f'create database {messengerdb};')
#         else:
#             pass
#
#     # Create the MESSENGER tables if necessary
#     with messengerdb_connect(configfile) as con:
#         cur = con.cursor()
#         cur.execute('select table_name from information_schema.tables')
#         tables = [r[0] for r in cur.fetchall()]
#
#         mestables = ['capointing', 'cauvvsdata', 'caspectra',
#                      'mgpointing', 'mguvvsdata', 'mgspectra',
#                      'napointing', 'nauvvsdata', 'naspectra',
#                      'mesmercyear']
#         there = [m in tables for m in mestables]
#
#         if not all(there):
#             # Delete any tables that may exist
#             for mestab in mestables:
#                 if mestab in tables:
#                     cur.execute(f'drop table {mestab}')
#                 else:
#                     pass
#
#             # Import the dumped tables
#             datafiles = glob.glob(os.path.join(datapath, 'UVVS*sql'))
#             for dfile in datafiles:
#                 print(f'Loading {os.path.basename(dfile)}')
#                 os.system(f'psql -d {config.mesdatabase} -p {config.port} -f {dfile}')
#         else:
#             pass
