from nexoclom import config
from MESSENGERuvvs.MESSENGERdata import MESSENGERdata
from .initialize_MESSENGERdata import initialize_MESSENGERdata


name = 'MESSENGERuvvs'
__author__ = 'Matthew Burger'
__email__ = 'mburger@stsci.edu'
__version__ = '1.9.13'
__date__ = '2022-08-24'


engine = config.create_engine(config.mesdatabase)
