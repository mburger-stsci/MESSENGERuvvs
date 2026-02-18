import os
import numpy as np
import pandas as pd
from astropy.time import Time
import astropy.units as u
from nexoclom2 import SSObject
from nexoclom2.initial_state import GeometryTime
from MESSENGERuvvs import __path__


basepath = __path__[0]


def create_merc_year_table():
    """Insert/read start date for each Mercury year from database.

    This creates and reads from database table *MESmercyear*
    """
    tstart = Time('2011-03-18T00:00:00', format='isot', scale='utc')
    tend = Time('2015-05-02T23:59:59', format='isot', scale='utc')
    
    geometry = GeometryTime({'center': 'Mercury',
                             'start_point': 'Mercury',
                             'modeltime': tend.iso})
    mercury = SSObject('Mercury')

    n = 10
    dt = (tend-tstart).to(u.s)/n
    times, taa = np.zeros((0, )), np.zeros((0, ))*u.rad
    times_ = np.linspace(-dt, 0*u.s, 1000)
    for i in range(n):
        geometry.modeltime = tend - dt*i
        print(geometry.modeltime)
        mercury.get_geometry(geometry, dt)
        taa_ = mercury.taa(times_)
        times = np.concatenate([times, geometry.modeltime + times_])
        taa = np.concatenate([taa, taa_])
    
    s = np.argsort(times)
    times, taa = times[s], taa[s]
    
    diff = taa[:-1] - taa[1:]
    q = diff > 2*u.rad
    
    sttimes = [times[0]]
    sttimes.extend(list(times[:-1][q]))
    endtimes = list(times[:-1][q])
    endtimes.append(times[-1])
    
    mercyear = pd.DataFrame({'yrnum': np.arange(len(sttimes), dtype=int),
                             'yrstart': sttimes,
                             'yrend': endtimes})
    mercyear.to_pickle(os.path.join(basepath, 'data', 'MES_merc_year.pkl'))
    
    return mercyear

if __name__ == '__main__':
    create_merc_year_table()
