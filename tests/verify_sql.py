import numpy as np
import os
import pandas as pd
import sqlalchemy as sqla
from astropy.visualization import quantity_support
import astropy.units as u
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib import rcParams
quantity_support()
rcParams.update({'lines.linewidth': 1})


query = sqla.text(
    'SELECT utc, x, y, z, loctimetan, longtan, lattan, alttan from data')
path = '/Users/mburger/Work/Data/MESSENGER/UVVS/Level2'

oldfile = os.path.join(path, '2025-10-30', 'MESSENGERuvvsCa.sqlite')
old_engine = sqla.create_engine(f'sqlite:///{oldfile}')
with old_engine.connect() as con:
    old = pd.read_sql(query, con)
old = old.iloc[:10000]
old.longtan = old.longtan.apply(np.degrees)
old.lattan = old.lattan.apply(np.degrees)

newfile = os.path.join(path, 'MESSENGERuvvsCa.sqlite')
new_engine = sqla.create_engine(f'sqlite:///{newfile}')
with new_engine.connect() as con:
    new = pd.read_sql(query, con)
new = new.iloc[:10000]
new.longtan = new.longtan.apply(np.degrees)
new.lattan = new.lattan.apply(np.degrees)

with  PdfPages('CoordComparison.pdf') as pdf:
    # x, y, z, alttan, longtan, lattan
    fig, ax = plt.subplots(3, 1, figsize=(16, 10))
    ax[0].plot(old.x, label='Summary File')
    ax[0].plot(new.x, label='SPICE')
    ax[0].set_title('Tanpoint x')
    
    ax[1].plot(old.y, label='Summary File')
    ax[1].plot(new.y, label='SPICE')
    ax[1].set_title('Tanpoint y')
    
    ax[2].plot(old.z, label='Summary File')
    ax[2].plot(new.z, label='SPICE')
    ax[2].set_title('Tanpoint z')
    pdf.savefig()
    plt.close()
    
    fig, ax = plt.subplots(3, 1, figsize=(16, 10))
    q = (old.loctimetan >= 0*u.deg) & (old.loctimetan < 360*u.deg)
    ax[0].plot(old.loctimetan[q], label='Summary File')
    ax[0].plot(new.loctimetan[q], label='SPICE')
    ax[0].set_title('Tanpoint Subsolar Local Time')
    
    q = (old.lattan >= -90*u.deg) & (old.lattan < 90*u.deg)
    ax[1].plot(old.lattan[q], label='Summary File')
    ax[1].plot(new.lattan[q], label='SPICE')
    ax[1].set_title('Tanpoint Subsolar Latitude')

    q = old.alttan <= 1e5*u.km
    ax[2].plot(old.alttan[q], label='Summary File')
    ax[2].plot(new.alttan[q], label='SPICE')
    ax[2].set_title('Tanpoint Altitude')

    ax[0].legend()
    pdf.savefig()
    plt.close()

    fig, ax = plt.subplots(3, 1, figsize=(16, 10))
    ax[0].plot(old.x-new.x, label='Summary File - SPICE')
    ax[0].set_title(r'Tanpoint $\Delta$x')

    ax[1].plot(old.y-new.y, label='Summary File - SPICE')
    ax[1].set_title(r'Tanpoint $\Delta$y')

    ax[2].plot(old.z-new.z, label='Summary File - SPICE')
    ax[2].set_title(r'Tanpoint $\Delta$z')
    
    ax[0].legend()
    pdf.savefig()
    plt.close()
    
    fig, ax = plt.subplots(3, 1, figsize=(16, 10))
    q = (old.longtan >= 0*u.deg) & (old.longtan < 360*u.deg)
    ax[0].plot(old.loctimetan[q]-new.loctimetan[q], label='Summary File - Spice')
    ax[0].set_title(r'Tanpoint $\Delta$ Subsolar Local Time')
    
    q = (old.lattan >= -90*u.deg) & (old.lattan < 90*u.deg)
    ax[1].plot(old.lattan[q]-new.lattan[q], label='Summary File - Spice')
    ax[1].set_title(r'Tanpoint $\Delta$ Subsolar Latitude')

    q = old.alttan <= 1e5*u.km
    ax[2].plot(old.alttan[q]-new.alttan[q], label='Summary File - Spice')
    ax[2].set_title(r'Tanpoint $\Delta$ Altitude')

    ax[0].legend()
    pdf.savefig()
    plt.close()

from inspect import currentframe, getframeinfo
frameinfo = getframeinfo(currentframe())
print(frameinfo.filename, frameinfo.lineno)
from IPython import embed; embed()
import sys; sys.exit()
