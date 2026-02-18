import os
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
from nexoclom2.solarsystem.load_kernels import SpiceKernels
from MESSENGERuvvs import __path__
from MESSENGERuvvs.get_datapath import get_datapath

basepath = __path__[0]


def determine_mercyear(datatime, config):
    yrnum = pd.read_pickle(os.path.join(basepath, 'data', 'MES_merc_year.pkl'))
    
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
    Rmerc = u.def_unit('R_Mercury', mercury.radius)
    # kR = u.def_unit('kR', 1e3*u.R)
    # nm = u.def_unit('nm', 1e-9*u.m)
    l1files = []

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
        UTC = Time(t_iso, format='yday').iso
        
        # Orbit number for each data spectrum
        orbit = np.array([int(o) for o in data['orb_num']])
        
        rmerc = (np.sqrt(np.sum(data['planet_sun_vector_tg']**2,
                                axis=1))*u.km).to(u.AU)
        
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
        
        kernels = SpiceKernels('Mercury')
        for i in np.arange(npts):
            # xyz[i, :] = np.dot(data['mso_rotation_matrix'][i, :, :],
            #                    data['planet_sc_vector_tg'][i, :]
            #                    )/mercury.radius.value
            # bore[i, :] = np.dot(data['mso_rotation_matrix'][i, :, :],
            #                     data['boresight_unit_vector_center_tg'][i, :])
            # corn0[i, :] = np.dot(data['mso_rotation_matrix'][i, :, :],
            #                      data['boresight_unit_vector_c1_tg'][i, :])
            # corn1[i, :] = np.dot(data['mso_rotation_matrix'][i, :, :],
            #                      data['boresight_unit_vector_c2_tg'][i, :])
            # corn2[i, :] = np.dot(data['mso_rotation_matrix'][i, :, :],
            #                      data['boresight_unit_vector_c3_tg'][i, :])
            # corn3[i, :] = np.dot(data['mso_rotation_matrix'][i, :, :],
            #                      data['boresight_unit_vector_c4_tg'][i, :])
            
            et = spice.str2et(UTC[i])
            R = spice.pxform('IAU_MERCURY', 'MERCURYSOLAR', et)
            
            xyz[i,:] = np.matmul(R, data['planet_sc_vector_tg'][i,:])/mercury.radius.value
            bore[i,:] = np.matmul(R, data['boresight_unit_vector_center_tg'][i,:])
            corn0[i,:] = np.matmul(R, data['boresight_unit_vector_c1_tg'][i,:])
            corn1[i,:] = np.matmul(R, data['boresight_unit_vector_c2_tg'][i,:])
            corn2[i,:] = np.matmul(R, data['boresight_unit_vector_c3_tg'][i,:])
            corn3[i,:] = np.matmul(R, data['boresight_unit_vector_c4_tg'][i,:])
            
        xcorner = np.column_stack([corn0[:, 0], corn1[:, 0],
                                   corn2[:, 0], corn3[:, 0]])
        ycorner = np.column_stack([corn0[:, 1], corn1[:, 1],
                                   corn2[:, 1], corn3[:, 1]])
        zcorner = np.column_stack([corn0[:, 2], corn1[:, 2],
                                   corn2[:, 2], corn3[:, 2]])
        
        # Determine tangent point
        t = -np.sum(xyz*bore, axis=1)
        behind = t < 0
        tanpt = xyz + bore*t[:, np.newaxis]
        rtan = np.linalg.norm(tanpt, axis=1)
        alttan = (rtan - 1) * mercury.radius.value
        lattan = np.arcsin(tanpt[:,2]/rtan)
        longtan = np.mod(np.arctan2(tanpt[:, 1], tanpt[:, 0]), 2*np.pi)
        lttan = np.mod(longtan*12/np.pi + 12, 24)

        tanpt[behind,:] = 1e30
        rtan[behind] = 1e30
        alttan[behind] = 1e30
        lattan[behind] = np.nan
        lttan[behind] = np.nan
        longtan[behind] = np.nan
        
        below = alttan < 0
        alttan[below] = 0
        tanpt[below, :] /= rtan[below, np.newaxis]
        rtan[below] = 1

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
        except:
            corrected = data['corr']
        dark = data['dark']
        solarfit = data['sol_fit']
        
        ndata = pd.DataFrame(
            {'species': species,
             'frame': 'MSO',
             'UTC': UTC,
             'orbit': orbit,
             'TAA': data['true_anomaly']*np.pi/180.,
             'rmerc': rmerc.value,
             'drdt': data['rad_vel'],
             'subslong': data['subsolar_longitude']*np.pi/180.,
             'g': data['gvals']/u.s,
             'radiance': radiance,
             'sigma': sigma,
             'x': xyz[:, 0]*Rmerc,
             'y': xyz[:, 1]*Rmerc,
             'z': xyz[:, 2]*Rmerc,
             'xbore': bore[:, 0], 'ybore': bore[:, 1], 'zbore': bore[:, 2],
             'xcorn1': xcorner[:, 0], 'xcorn2': xcorner[:, 1],
             'xcorn3': xcorner[:, 2], 'xcorn4': xcorner[:, 3],
             'ycorn1': ycorner[:, 0], 'ycorn2': ycorner[:, 1],
             'ycorn3': ycorner[:, 2], 'ycorn4': ycorner[:, 3],
             'zcorn1': zcorner[:, 0], 'zcorn2': zcorner[:, 1],
             'zcorn3': zcorner[:, 2], 'zcorn4': zcorner[:, 3],
             'obstype': obstype,
             'obstype_num': data['obs_typ_num'],
             'xtan': tanpt[:, 0],
             'ytan': tanpt[:, 1],
             'ztan': tanpt[:, 2],
             'rtan': rtan,
             'alttan': alttan,
             'longtan': longtan,
             'lattan': lattan,
             'loctimetan': lttan,
             'minalt': data['minalt'],
             # 'longtan': data['target_longitude_set'][:, 0]*np.pi/180, #IAU coords
             # 'lattan': data['target_latitude_set'][:, 0]*np.pi/180,
             # 'alttan': data['target_altitude_set'][:, 0],
             # 'loctimetan': data['obs_solar_localtime'], ## This is s/c LT
             'slit': slit})
        ndata.fillna(-999, inplace=True)
        
        spectra = [spectra[i,:] for i in range(spectra.shape[0])]
        wavelength = [wavelength[i,:] for i in range(wavelength.shape[0])]
        raw = [raw[i,:] for i in range(raw.shape[0])]
        corrected = [corrected[i,:] for i in range(corrected.shape[0])]
        dark = [dark[i,:] for i in range(dark.shape[0])]
        solarfit = [solarfit[i,:] for i in range(solarfit.shape[0])]
        spectra = pd.DataFrame(
            {'UTC': UTC,
             'spectra': spectra,
             'wavelength': wavelength,
             'raw': raw,
             'corrected': corrected,
             'dark': dark,
             'solarfit': solarfit})
        
        # save this for later
        newfile = picklefile.replace('L0.pkl', 'L1.pkl').replace('Level0', 'Level1')
        print(f'Saving {newfile}')
        l1files.append(newfile)
        ndata.to_pickle(newfile)
        spectra.to_pickle(newfile.replace('.pkl', '_spectra.pkl'))
        
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
        
        spectra = pd.concat([pd.read_pickle(file.replace('.pkl', '_spectra.pkl'))
                             for file in files])
        spectra.sort_values('UTC', inplace=True)
        
        dbfile = os.path.join(datapath, 'Level2', f'MESSENGERuvvs{sp}.sqlite')
        engine = sqla.create_engine(f'sqlite:///{dbfile}',
                                    connect_args={"autocommit": True})
        with engine.connect() as con:
            print(sp, 'data')
            data.to_sql('data', con)
            
        print(sp, 'spectra')
        specfile = os.path.join(datapath, 'Level2', f'MESSENGERuvvs{sp}.pkl')
        spectra.to_pickle(specfile)


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
