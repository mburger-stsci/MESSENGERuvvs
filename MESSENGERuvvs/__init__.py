from .databasebackups import databasebackups
from .database_setup import messenger_database_setup, database_connect
from .MESSENGERdata import MESSENGERdata
from .initialize_MESSENGERdata import initialize_MESSENGERdata


name = 'MESSENGERuvvs'
__author__ = 'Matthew Burger'
__email__ = 'mburger@stsci.edu'
__version__ = '1.2.10'


messenger_database_setup()
