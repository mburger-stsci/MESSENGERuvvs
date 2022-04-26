import psycopg2
from nexoclom.utilities.NexoclomConfig import NexoclomConfig


DEFAULT_MESSENGER_DB = 'messengeruvvsdb'

def messengerdb_connect():
    config = NexoclomConfig()
    messengerdb = config.__dict__.get('mesdatabase', DEFAULT_MESSENGER_DB)
    
    if config.dbhost:
        con = psycopg2.connect(host=config.dbhost, dbname=messengerdb,
                              port=config.port)
    else:
        con = psycopg2.connect(dbname=messengerdb, port=config.port)
        
    con.autocommit = True
    return con
