import os
import numpy as np
import pandas as pd
import pickle
import glob
from scipy import io
from astropy.time import Time
from astropy import units as u
from sqlalchemy import text
from nexoclom.solarsystem import SSObject, planet_geometry
from nexoclom.utilities.NexoclomConfig import NexoclomConfig
from nexoclom.utilities.exceptions import ConfigfileError
from MESSENGERuvvs import __file__ as basefile


basepath = os.path.dirname(basefile)


def create_merc_year_table(config):
    """Insert/read start date for each Mercury year from database.

    This creates and reads from database table *MESmercyear*
    """
    tstart = Time('2011-03-18T00:00:00', format='isot', scale='utc')
    tend = Time('2015-05-02T23:59:59', format='isot', scale='utc')
    
    times = Time(np.arange(tstart.mjd, tend.mjd), format='mjd', scale='utc')
    taa = np.ndarray((len(times),))*u.rad
    for i, time in enumerate(times):
        geo = planet_geometry(time, 'Mercury')
        taa[i] = geo['taa']
        
    diff = taa[:-1] - taa[1:]
    q = np.where(diff > 0)[0]
    sttimes = np.append(times[0], times[q])
    endtimes = np.append(times[q], times[-1])
    
    mercyear = pd.DataFrame({'yrnum': np.arange(len(sttimes), dtype=int),
                             'yrstart': sttimes,
                             'yrend': endtimes})
    mercyear.to_pickle(os.path.join(basepath, 'data', 'MES_merc_year.pkl'))
    

def determine_mercyear(datatime, config):
    yrnum = pd.read_pickle(os.path.join(basepath, 'data', 'MES_merc_year.pkl'))
    
    datatime_ = datatime.apply(pd.Timestamp)
    myear = np.zeros((len(datatime),), dtype=int) - 1
    for _, yr in yrnum.iterrows():
        q = (datatime_ >= yr.yrstart) & (datatime_ < yr.yrend)
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
        bore = np.ndarray((npts, 3))
        corn0 = np.ndarray((npts, 3))
        corn1 = np.ndarray((npts, 3))
        corn2 = np.ndarray((npts, 3))
        corn3 = np.ndarray((npts, 3))
        for i in np.arange(npts):
            xyz[i, :] = np.dot(data['mso_rotation_matrix'][i, :, :],
                               data['planet_sc_vector_tg'][i, :]
                               )/mercury.radius.value
            bore[i, :] = np.dot(data['mso_rotation_matrix'][i, :, :],
                                data['boresight_unit_vector_center_tg'][i, :])
            corn0[i, :] = np.dot(data['mso_rotation_matrix'][i, :, :],
                                 data['boresight_unit_vector_c1_tg'][i, :])
            corn1[i, :] = np.dot(data['mso_rotation_matrix'][i, :, :],
                                 data['boresight_unit_vector_c2_tg'][i, :])
            corn2[i, :] = np.dot(data['mso_rotation_matrix'][i, :, :],
                                 data['boresight_unit_vector_c3_tg'][i, :])
            corn3[i, :] = np.dot(data['mso_rotation_matrix'][i, :, :],
                                 data['boresight_unit_vector_c4_tg'][i, :])
        
        xcorner = np.array([corn0[:, 0], corn1[:, 0],
                            corn2[:, 0], corn3[:, 0]]).transpose()
        ycorner = np.array([corn0[:, 1], corn1[:, 1],
                            corn2[:, 1], corn3[:, 1]]).transpose()
        zcorner = np.array([corn0[:, 2], corn1[:, 2],
                            corn2[:, 2], corn3[:, 2]]).transpose()
        
        # Determine tangent point
        t = -np.sum(xyz*bore, axis=1)
        tanpt = xyz+bore*t[:, np.newaxis]
        rtan = np.linalg.norm(tanpt, axis=1)
        
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
             'xtan': tanpt[:, 0], 'ytan': tanpt[:, 1],
             'ztan': tanpt[:, 2], 'rtan': rtan,
             'alttan': data['target_altitude_set'][:, 0],
             'minalt': data['minalt'],
             'longtan': data['target_longitude_set'][:, 0]*np.pi/180,
             'lattan': data['target_latitude_set'][:, 0]*np.pi/180,
             'loctimetan': data['obs_solar_localtime'],
             'slit': slit})
        ndata.fillna(-999, inplace=True)
        
        spectra = [spectra[i,:] for i in range(spectra.shape[0])]
        wavelength = [wavelength[i,:] for i in range(wavelength.shape[0])]
        raw = [raw[i,:] for i in range(raw.shape[0])]
        corrected = [corrected[i,:] for i in range(corrected.shape[0])]
        dark = [dark[i,:] for i in range(dark.shape[0])]
        solarfit = [solarfit[i,:] for i in range(solarfit.shape[0])]
        spectra = pd.DataFrame(
            {'spectra': spectra,
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
    config = NexoclomConfig()
    
    if 'mesdatabase' not in config.__dict__:
        raise ConfigfileError('mesdatabase', config.configfile)
    else:
        pass
    
    # Set up database
    try:
        os.system(f'dropdb {config.mesdatabase}')
    except:
        pass

    os.system(f'createdb {config.mesdatabase}')

    # print('creating MESmercyear table')
    # create_merc_year_table(config)

    print('creating UVVS tables')
    spec = ['Ca', 'Na', 'Mg']
    for sp in spec:
        with config.database_connect(config.mesdatabase) as con:
            cur = con.cursor()
            # Table with spectrum information
            print(f'Creating {sp}uvvsdata')
            cur.execute(f'''CREATE table {sp}uvvsdata (
                               unum SERIAL PRIMARY KEY,
                               species text,
                               frame text,
                               UTC timestamp,
                               orbit int,
                               merc_year int,
                               taa float,
                               rmerc float,
                               drdt float,
                               subslong float,
                               g float,
                               radiance float,
                               sigma float)''')

            # Table with MESSENGER geometry and UVVS pointing
            print(f'Creating {sp}pointing')
            cur.execute(f'''CREATE table {sp}pointing (
                               pnum SERIAL PRIMARY KEY,
                               x float,
                               y float,
                               z float,
                               xbore float,
                               ybore float,
                               zbore float,
                               obstype text,
                               obstype_num int,
                               xtan float,
                               ytan float,
                               ztan float,
                               rtan float,
                               alttan float,
                               longtan float,
                               lattan float,
                               loctimetan float,
                               slit text)''')  # Not including slit corners

            # Table with spectra
            print(f'Creating {sp}spectra')
            cur.execute(f'''CREATE table {sp}spectra (
                                snum SERIAL PRIMARY KEY,
                                wavelength float[],
                                calibrated float[],
                                raw float[],
                                dark float[],
                                solarfit float[])''')
 
    for l1file in l1files:
        print(f'Processing {l1file}')
        ndata = pd.read_pickle(l1file)
        species = ndata.species.unique()
        if len(species) == 1:
            species = species[0]
        else:
            assert False
            
        # add Mercury year
        ndata['merc_year'] = determine_mercyear(ndata.UTC, config)
        
        print('Inserting UVVS data')
        print(f'Saving {species} Data')
        
        data_query = text(
            f'''INSERT into {species}uvvsdata (species, frame, UTC, orbit,
                    merc_year, taa, rmerc, drdt, subslong, g, radiance,
                    sigma)
                VALUES (:species, :frame, :UTC, :orbit, :merc_year, :taa,
                        :rmerc, :drdt, :subslong, :g, :radiance, :sigma)''')
        
        point_query = text(
            f'''INSERT into {species}pointing (x, y, z, xbore, ybore, zbore,
                    obstype, obstype_num, xtan, ytan, ztan,
                    rtan, alttan, longtan, lattan,
                    loctimetan, slit)
                VALUES (:x, :y, :z, :xbore, :ybore, :zbore, :obstype,
                    :obstype_num, :xtan, :ytan, :ztan, :rtan, :alttan,
                    :longtan, :lattan, :loctimetan, :slit)''')
                    
        with config.create_engine(config.mesdatabase).begin() as con:
            for i, dpoint in ndata.iterrows():
                data_values = {'species': dpoint.species,
                                'frame': dpoint.frame,
                                'UTC': dpoint.UTC,
                                'orbit': dpoint.orbit,
                                'merc_year': dpoint.merc_year,
                                'taa': dpoint.TAA,
                                'rmerc': dpoint.rmerc,
                                'drdt': dpoint.drdt,
                                'subslong': dpoint.subslong,
                                'g': dpoint.g,
                                'radiance': dpoint.radiance,
                                'sigma': dpoint.sigma}
                con.execute(data_query, data_values)
                
                point_values = {'x': dpoint.x,
                                'y': dpoint.y,
                                'z': dpoint.z,
                                'xbore': dpoint.xbore,
                                'ybore': dpoint.ybore,
                                'zbore': dpoint.zbore,
                                'obstype': dpoint.obstype,
                                'obstype_num': dpoint.obstype_num,
                                'xtan': dpoint.xtan,
                                'ytan': dpoint.ytan,
                                'ztan': dpoint.ztan,
                                'rtan': dpoint.rtan,
                                'alttan': dpoint.alttan,
                                'longtan': dpoint.longtan,
                                'lattan': dpoint.lattan,
                                'loctimetan': dpoint.loctimetan,
                                'slit': dpoint.slit}
                con.execute(point_query, point_values)
                
        spectra = pd.read_pickle(l1file.replace('.pkl', '_spectra.pkl'))
        print(f'Saving {species} Spectra')
        spec_query = text(
            f'''INSERT into {species}spectra (wavelength,
                    calibrated, raw, dark, solarfit)
                VALUES (:wavelength, :calibrated, :raw, :dark, :solarfit)''')
        with config.create_engine(config.mesdatabase).begin() as con:
            for i, spectrum in spectra.iterrows():
                spec_values = {'wavelength': spectrum.wavelength.tolist(),
                               'calibrated': spectrum.corrected.tolist(),
                               'spectra': spectrum.spectra.tolist(),
                               'raw': spectrum.raw.tolist(),
                               'dark': spectrum.dark.tolist(),
                               'solarfit': spectrum.solarfit.tolist()}
                con.execute(spec_query, spec_values)


def initialize_MESSENGERdata(idl_convert=False, to_level1=False, to_sql=True):
    if idl_convert:
        pfiles = convert_IDL_to_pickle(
            '/Users/mburger/Work/Data/MESSENGER/ModelData/SummaryFiles/V0001',
            '/Users/mburger/Work/Data/MESSENGER/ModelData/Level0')
    else:
        pfiles = glob.glob(
            '/Users/mburger/Work/Data/MESSENGER/ModelData/Level0/*L0.pkl')
        
    if to_level1:
        l1files = process_L0_pickle(pfiles)
    else:
        config = NexoclomConfig()
        l1files = glob.glob(os.path.join(config.mesdatapath, 'Level1', '*L1.pkl'))
        
    if to_sql:
        set_up_database(l1files)
    else:
        pass


if __name__ == '__main__':
    initialize_MESSENGERdata(to_level1=False, to_sql=True)
