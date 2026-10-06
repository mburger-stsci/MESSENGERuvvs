""" Converts IDL savefiles to pickle and sqlite files.

Notes
-----
The data are converted from IAU_MERCURY frame to MercurySolarFixed (MSGR_MSO) using
SPICE
"""
import os
from sysconfig import expand_makefile_vars

import numpy as np
import pandas as pd
import pickle
import glob
from scipy import io
from astropy.time import Time
from astropy import units as u
import sqlalchemy as sqla
import spiceypy as spice
from nexoclom2.solarsystem import SSObject
from MESSENGERuvvs import path
from MESSENGERuvvs.get_datapath import get_datapath
from MESSENGERuvvs.load_MESSENGER_kernels import MESSENGER_kernels


def determine_mercyear(datatime, config):
    yrnum = pd.read_pickle(os.path.join(path, 'data', 'MES_merc_year.pkl'))
    
    # datatime_ = datatime.apply(pd.Timestamp)
    myear = np.zeros((len(datatime),), dtype=int) - 1
    for _, yr in yrnum.iterrows():
        q = (datatime >= yr.yrstart) & (datatime < yr.yrend)
        myear[q] = yr.yrnum
    
    assert np.all(myear > -1)
    
    return myear

    
def convert_IDL_to_pickle(path_to_out, path_for_pkl):
    """Store data from IDL summary files in a database.
    The IDL summary files were provided by Aimee Merkel.
    
    Two tables are created for each species (Ca, Na, and Mg): *xxuvvsdata* and
    *xxuvvspointing* where xx is the species. See :doc:`database_fields` for
    a description of these tables and fields.
    
    **Parameters**
    
    datapath
        Path to the IDL summary files
        
    **Returns**
    
    No output.
    """
    idlfiles = glob.glob(os.path.join(path_to_out, '*.sav'))
    picklefiles = []
    for idlfile in idlfiles:
        idldata = io.readsav(idlfile)
        pydict = {key: item for key, item in idldata.items()}
        newfile = os.path.join(path_for_pkl,
                               os.path.basename(idlfile).replace('.sav', '.L0.pkl'))
        picklefiles.append(newfile)
        with open(newfile, 'wb') as pfile:
            pickle.dump(pydict, pfile)
         
    return picklefiles

    
def process_L0_pickle(picklefiles):
    mercury = SSObject('Mercury')
    sun = SSObject('Sun')
    Rmerc = u.def_unit('R_Mercury', mercury.radius)
    l1files = []
    kernels = MESSENGER_kernels()
    
    merc_year = pd.read_pickle(os.path.join(path, 'data', 'MES_merc_year.pkl'))
    year_start = Time(merc_year.yrstart.values)
    year_end = Time(merc_year.yrend.values)
    
    for picklefile in picklefiles:
        print(picklefile)
        with open(picklefile, 'rb') as pfile:
            data = pickle.load(pfile)

        npts = len(data['orb_num'])
        species = os.path.basename(picklefile)[0:2].lower()
        
        # Determine UT for each spectrum
        t_iso = [f"20{time[0:2].decode('utf-8')}:{time[2:5].decode('utf-8')}:"
                 f"{time[6:].decode('utf-8')}"
                 for time in data['step_utc_time']]
        UTC = Time(t_iso, format='yday')
        times_et = spice.str2et(UTC.iso)
    
        # Orbit number for each data spectrum
        orbit = np.array([int(o) for o in data['orb_num']])
        
        # Mercury info
        rmerc, drdt_merc, taa = np.zeros(npts), np.zeros(npts), np.zeros(npts)
        subsolar_long, subsolar_lat = np.zeros(npts), np.zeros(npts)
       
        # Radiance
        radiance = data[f'{species.lower()}_tot_rad_kr']
        sigma = radiance/data[f'{species.lower()}_tot_rad_snr']
        
        # Spacecraft position and boresight in MSO
        xyz = np.ndarray((npts, 3))
        r = np.linalg.norm(xyz, axis=1)
        bore = np.ndarray((npts, 3))
        corn0 = np.ndarray((npts, 3))
        corn1 = np.ndarray((npts, 3))
        corn2 = np.ndarray((npts, 3))
        corn3 = np.ndarray((npts, 3))
        
        year = np.ndarray(npts)
        
        # tangent point
        tanpt, rtan, alttan = np.zeros((npts, 3)), np.zeros(npts), np.zeros(npts)
        longtan, lattan = np.zeros(npts), np.zeros(npts)
        rangetan = np.zeros(npts)
        longtanIAU, lattanIAU = np.zeros(npts), np.zeros(npts)
        minalt = np.zeros(npts)
        
        uvvs_atmos_naifid = -236600
        uvvs_surf_naifid = -236610
        msgr_naifid = '-236'
        _, inst_atmos_frame, boresight, _, corners_atmos = spice.getfov(uvvs_atmos_naifid,
                                                                        100, 100, 100)
        _, inst_surf_frame, boresight, _, corners_surf = spice.getfov(uvvs_surf_naifid,
                                                                      100, 100, 100)
        
        data['slit'] = data['slit'].astype(int)
        slit_frame = {0: inst_surf_frame, 1: inst_atmos_frame}
        slit_corners = {0: corners_surf, 1: corners_atmos}
        
        for i, et in enumerate(times_et):
            # Determine Mercury position
            state, _ = spice.spkezr('Mercury', et, 'J2000', 'CN+S', 'Sun')
            rmerc[i] = np.sqrt(np.sum(state[:3]**2))*u.km.to(u.au)
            drdt_merc[i] = np.sum(state[:3]*state[3:])/(rmerc[i]*u.au.to(u.km))
            taa[i] = np.degrees(spice.oscltx(state, et, -sun.GM.value)[8])
            
            # Determine mercury year
            yr = np.where((UTC[i] >= year_start) & (UTC[i] < year_end))[0]
            assert len(yr) == 1
            year[i] = yr[0]
            
            # Determine s/c location
            state, _ = spice.spkezr('MESSENGER', et, 'MercurySolarFixed', 'CN+S', 'Mercury')
            xyz[i,:] = state[:3]/mercury.radius.value
            
            # Determine s/c boresight
            s = data['slit'][i]
            R =  spice.pxform(slit_frame[s], 'MERCURYSOLARFIXED', et)
            bore[i,:] = np.matmul(R, boresight)
            # corn0[i,:] = np.matmul(R, slit_corners[s][0,:])
            # corn1[i,:] = np.matmul(R, slit_corners[s][1,:])
            # corn2[i,:] = np.matmul(R, slit_corners[s][2,:])
            # corn3[i,:] = np.matmul(R, slit_corners[s][3,:])
            tp = spice.tangpt('ELLIPSOID', 'Mercury', et, 'MercurySolarFixed', 'CN+S',
                              'TANGENT POINT', msgr_naifid, slit_frame[s], boresight)
            
            # Line of sight tangent point in MSO
            tanpt[i,:] = tp[0]/mercury.radius.value
            r, lon, lat = spice.reclat(tp[0])
            rtan[i] = r/mercury.radius.value
            longtan[i] = np.mod(np.degrees(lon), 360.)
            lattan[i] = np.degrees(lat)
            alttan[i] = tp[1]
            rangetan[i] = tp[2]/mercury.radius.value
            
            minalt[i] = 1e30
            for j in range(4):
                tp0 = spice.tangpt('ELLIPSOID', 'Mercury', et, 'MercurySolarFixed', 'CN+S',
                                   'TANGENT POINT', msgr_naifid, slit_frame[s],
                                   slit_corners[s][j,:])
                minalt[i] = min(minalt[i], tp0[1])
                
            # Tangent point lon/lat in IAU frame
            tp = spice.tangpt('ELLIPSOID', 'Mercury', et, 'IAU_Mercury', 'CN+S',
                              'TANGENT POINT', msgr_naifid, slit_frame[s], boresight)
            r, lon, lat = spice.reclat(tp[0])
            longtanIAU[i] = np.mod(np.degrees(lon), 360.)
            lattanIAU[i] = np.degrees(lat)
            
            # Subsolar point in IAU frame
            sublon, _, _ = spice.subslr('INTERCEPT/ELLIPSOID', 'Mercury', et, 'IAU_Mercury',
                                        'CN+S', 'Sun')
            r, lon, lat = spice.reclat(sublon)
            subsolar_long[i] = np.mod(np.degrees(lon), 360.)
            subsolar_lat[i] = np.degrees(lat)
            
        # xcorner = np.column_stack([corn0[:, 0], corn1[:, 0],
        #                            corn2[:, 0], corn3[:, 0]])
        # ycorner = np.column_stack([corn0[:, 1], corn1[:, 1],
        #                            corn2[:, 1], corn3[:, 1]])
        # zcorner = np.column_stack([corn0[:, 2], corn1[:, 2],
        #                            corn2[:, 2], corn3[:, 2]])
        
        # Determine tangent point
        lttan = np.mod(longtan*12/180 + 12, 24)
        
        point_away = rangetan == 0
        tanpt[point_away,:] = 1e30
        rtan[point_away] = 1e30
        alttan[point_away] = 1e30
        minalt[point_away] = 1e30
        longtan[point_away] = np.nan
        lattan[point_away] = np.nan
        lttan[point_away] = np.nan

        slit = np.array(['Surface'
                         if s == 0
                         else 'Atmospheric'
                         for s in data['slit']])
        obstype = np.array([str(ob).replace('b', '').replace("'", '').strip()
                            for ob in data['obs_typ']])
        
        # Add in the spectra
        spectra = data[species.lower()+'_rad_kr']
        wavelength = data['wavelength']
        raw = data['orig']
        try:
            corrected = data['fully_corr_cr']
        except KeyError:
            corrected = data['corr']
        dark = data['dark']
        solarfit = data['sol_fit']
        
        ndata = pd.DataFrame(
            {'species': species.capitalize(),
             'frame': 'MercurySolarFixed',
             'UTC': UTC,
             'orbit': orbit,
             'TAA': taa,
             'year': year,
             'rmerc': rmerc,
             'drdt': drdt_merc,
             'subslong': subsolar_long,
             'subslat': subsolar_lat,
             'radiance': radiance,
             'sigma': sigma,
             'x': xyz[:, 0]*Rmerc,
             'y': xyz[:, 1]*Rmerc,
             'z': xyz[:, 2]*Rmerc,
             'xbore': bore[:, 0], 'ybore': bore[:, 1], 'zbore': bore[:, 2],
             # 'xcorn1': xcorner[:, 0], 'xcorn2': xcorner[:, 1],
             # 'xcorn3': xcorner[:, 2], 'xcorn4': xcorner[:, 3],
             # 'ycorn1': ycorner[:, 0], 'ycorn2': ycorner[:, 1],
             # 'ycorn3': ycorner[:, 2], 'ycorn4': ycorner[:, 3],
             # 'zcorn1': zcorner[:, 0], 'zcorn2': zcorner[:, 1],
             # 'zcorn3': zcorner[:, 2], 'zcorn4': zcorner[:, 3],
             'obstype': obstype,
             # 'obstype_num': data['obs_typ_num'],
             'xtan': tanpt[:, 0],
             'ytan': tanpt[:, 1],
             'ztan': tanpt[:, 2],
             'rtan': rtan,
             'alttan': alttan,
             'minalt': minalt,
             'longtan': longtan,
             'lattan': lattan,
             'loctimetan': lttan,
             'rangetan': rangetan,
             'longtanIAU': longtanIAU,
             'lattanIAU': lattanIAU,
             'slit': slit})
        # ndata.fillna(-999, inplace=True)
        
        # spectra = [spectra[i,:] for i in range(spectra.shape[0])]
        # wavelength = [wavelength[i,:] for i in range(wavelength.shape[0])]
        # raw = [raw[i,:] for i in range(raw.shape[0])]
        # corrected = [corrected[i,:] for i in range(corrected.shape[0])]
        # dark = [dark[i,:] for i in range(dark.shape[0])]
        # solarfit = [solarfit[i,:] for i in range(solarfit.shape[0])]
        # spectra = pd.DataFrame(
        #     {'UTC': UTC,
        #      'spectra': spectra,
        #      'wavelength': wavelength,
        #      'raw': raw,
        #      'corrected': corrected,
        #      'dark': dark,
        #      'solarfit': solarfit})
        
        # save this for later
        newfile = picklefile.replace('L0.pkl', 'L1.pkl').replace('Level0', 'Level1')
        print(f'Saving {newfile}')
        l1files.append(newfile)
        ndata.to_pickle(newfile)
        # spectra.to_pickle(newfile.replace('.pkl', '_spectra.pkl'))
        
    kernels.unload()
    
    return l1files
    
 
def set_up_database(l1files):
    datapath = get_datapath()
    
    # print('creating MESmercyear table')
    # create_merc_year_table(datapath)
    
    print('creating UVVS tables')
    species = 'Ca', 'Na', 'Mg'
    
    for sp in species:
        files = [file for file in l1files
                 if ((os.path.basename(file).lower().startswith(sp.lower())) and
                     ('spectra' not in file))]
        data = pd.concat([pd.read_pickle(file) for file in files])
        data.sort_values('UTC', inplace=True)
        data.UTC = data.UTC.apply(lambda x: x.iso)
        
        # spectra = pd.concat([pd.read_pickle(file.replace('.pkl', '_spectra.pkl'))
        #                      for file in files])
        # spectra.sort_values('UTC', inplace=True)
        
        dbfile = os.path.join(datapath, 'Level2', f'MESSENGERuvvs{sp}.sqlite')
        engine = sqla.create_engine(f'sqlite:///{dbfile}',
                                    connect_args={"autocommit": True})
        with engine.connect() as con:
            print(sp, 'data')
            data.to_sql('data', con)
        
        # print(sp, 'spectra')
        # spectra.to_pickle(specfile)

# specfile = os.path.join(datapath, 'Level2', f'MESSENGERuvvs{sp}.pkl')


def initialize_MESSENGERdata(idl_convert=False, to_level1=False, to_sql=True):
    datapath = get_datapath()
    if idl_convert:
        pfiles = convert_IDL_to_pickle(
            os.path.join(datapath, 'SummaryFiles', 'V0001'),
            os.path.join(datapath, 'Level0'))
    else:
        pfiles = glob.glob(os.path.join(datapath, 'Level0', '*'))
        
    if to_level1:
        l1files = process_L0_pickle(pfiles)
    else:
        l1files = glob.glob(os.path.join(datapath, 'Level1', '*L1.pkl'))
        
    if to_sql:
        set_up_database(l1files)
    else:
        pass
    

if __name__ == '__main__':
    initialize_MESSENGERdata(to_level1=True, to_sql=True)
    # create_MESSSENGER_summary()
