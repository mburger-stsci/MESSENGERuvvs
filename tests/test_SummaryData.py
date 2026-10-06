import os
import numpy as np
import astropy.units as u
from astropy.time import Time
import glob
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from nexoclom2 import SSObject

from inspect import currentframe, getframeinfo
from IPython import embed
import sys

def test_SummaryData():
    """ Compare old summary data to spice generated data
    """
    
    mercury = SSObject('Mercury')
    datadir = '/Users/mburger/Work/Data/MESSENGER/UVVS/Level1'
    species = 'ca', 'NA', 'MG'
    
    # cols = ('TAA', 'rmerc', 'drdt', 'subslong', 'x', 'y',
    #         'z', 'magbore', 'xbore', 'ybore', 'zbore', 'xcorn1', 'xcorn2', 'xcorn3', 'xcorn4',
    #         'ycorn1', 'ycorn2', 'ycorn3', 'ycorn4', 'zcorn1', 'zcorn2', 'zcorn3',
    #         'zcorn4', 'xtan', 'ytan', 'ztan', 'rtan', 'alttan', 'longtan', 'lattan',
    #         'loctimetan')
    cols = ('TAA', 'rmerc', 'drdt', 'subslong', 'x', 'y',
            'z', 'magbore', 'xbore', 'ybore', 'zbore',
            'xtan', 'ytan', 'ztan', 'rtan', 'alttan', 'longtan', 'lattan',
            'loctimetan')

    loadnew = False
    loadold = False
    for sp in species:
        pklfiles = glob.glob(os.path.join(datadir, f'{sp}_*.pkl'))
        
        if loadnew:
            print(sp, 'loading new')
            new_ = [pd.read_pickle(pklfile) for pklfile in pklfiles]
            new = pd.concat(new_, ignore_index=True)
            del new_
            new['mjd'] = new.UTC.apply(lambda x: x.mjd)
            new.sort_values('mjd', inplace=True)
            new.to_pickle(f'Summary_{sp}_new.pkl')
        else:
            new = pd.read_pickle(f'Summary_{sp}_new.pkl')
        
        if loadold:
            print(sp, 'loading old')
            old_ = [pd.read_pickle(os.path.join(datadir, '2026-09-17', os.path.basename(pklfile)))
                   for pklfile in pklfiles]
            old = pd.concat(old_, ignore_index=True)
            del old_
            old['UTC'] = old.UTC.apply(Time)
            old['mjd'] = old.UTC.apply(lambda x: x.mjd)
            old.sort_values('mjd', inplace=True)
            for col in ('TAA', 'subslong', 'longtan', 'lattan'):
                old[col] = old[col]*180/np.pi
            old.to_pickle(f'Summary_{sp}_old.pkl')
        else:
            old = pd.read_pickle(f'Summary_{sp}_old.pkl')
        
        new['magbore'] = np.sqrt(new.xbore**2 + new.ybore**2 + new.zbore**2)
        old['magbore'] = np.sqrt(old.xbore**2 + old.ybore**2 + old.zbore**2)
        
        q = np.abs(new['TAA'] - old['TAA']) < 10
        new, old = new[q], old[q]
        
        q = np.abs(new['subslong'] - old['subslong']) < 10
        new, old = new[q], old[q]
        
        q = np.abs(new['longtan'] - old['longtan']) < 10
        new, old = new[q], old[q]
        
        q = np.abs(new['loctimetan'] - old['loctimetan']) < 12
        new, old = new[q], old[q]
        
        q = np.isclose(old['magbore'], 1)
        new, old = new[q], old[q]
        
        # Remove where intersect the surface
        q = new.alttan > 0
        new, old = new[q], old[q]
        w = old.alttan > 0
        new, old = new[w], old[w]
        
        q = (new.rangetan == 0) == (old.rtan == 1e30)
        new, old = new[q], old[q]
        
        w = old.rtan == 1e30
        old.loc[w, ['longtan', 'lattan', 'lttan']] = -999
        
        # if set(new[q].index) != set(old[w].index):
        #     frameinfo = getframeinfo(currentframe())
        #     print(frameinfo.filename, frameinfo.lineno)
        #     embed()
        #     sys.exit()
        # new, old = new[np.logical_not(q)], old[np.logical_not(q)]
        
        
        # alttan = (new.rtan - 1)*mercury.radius.value
        
        pdf = PdfPages(f'Summary_{sp.capitalize()}_filtered.pdf')
        for col in cols:
            fig, ax = plt.subplots(2, 1, figsize=(12, 12), sharex=True)
            mx = np.max(np.abs(new[col] - old[col]))
            med = np.median(np.abs(new[col] - old[col]))
            
            print(f'* max($|\Delta {col}|$) = {mx:1.1e}, median($|\Delta {col}|$) = {med:1.1e}')
            
            ax[0].plot(new.mjd, new[col], label='With SPICE', color='black')
            ax[0].plot(new.mjd, old[col], label='From Summary', color='red', linewidth=1)
            ax[0].set_label('New vs. Old')
            ax[0].legend()
            
            ax[1].plot(new.mjd, new[col] - old[col])
            ax[1].set_ylabel('New - Old')
            
            fig.suptitle(col)
            pdf.savefig(fig)
            plt.close()
        pdf.close()
        
        embed()
        sys.exit()
        
        del new, old
    

if __name__ == '__main__':
    test_SummaryData()
