''' MESSENGER UVVS data class
Convert MESSENGER data into a useable class
*_temp.pickle are straight pickled versions of the IDL summmary files
Creates files that match my previous IDL summary files
This needs to be rerun if more methods are added to the class
'''
import numpy as np
import os
import sys
import glob
import pickle
import pandas as pd
import psycopg2
from scipy import io
from astropy.time import Time
from astropy import units as u
from astropy.visualization import PercentileInterval
from solarsystemMB import SSObject, planet_geometry


def get_database():
    configfile = os.path.join(os.environ['HOME'], '.nexoclom')
    config = {}
    if os.path.isfile(configfile):
        for line in open(configfile, 'r').readlines():
            key, value = line.split('=')
            if key.strip() == 'database':
                database = value.strip()
            else:
                pass
    else:
        database = None

    return database


def merc_year(datatime=None, initialize=False):
    '''Insert/read start date for each Mercury year from database.'''
    mercury = SSObject('Mercury')

    tstart = Time('2011-03-18T00:00:00', format='isot', scale='utc')
    tend = Time('2015-04-30T23:59:59', format='isot', scale='utc')

    database = get_database()
    con = psycopg2.connect(database=database)
    con.autocommit = True
    if initialize:
        times_ = np.arange(tstart.jd, tend.jd)
        times = [Time(t, format='jd', scale='utc') for t in times_]

        taa = np.ndarray(len(times))*u.rad
        for i,t in enumerate(times):
            time = Time(t, format='jd', scale='utc')
            geo = planet_geometry(time, 'Mercury')
            taa[i] = geo['taa']

        styear = [times[0]]
        for a,b,c in zip(taa[0:-1], taa[1:], times[1:]):
            if a > b:
                styear.append(c)
                print(c.iso)
        endyr = [*styear[1:], tend]

        cur = con.cursor()
        try:
            cur.execute('DROP table MESmercyear')
        except:
            pass
        cur.execute('''CREATE table MESmercyear
                        (yrnum int PRIMARY KEY,
                         yrstart timestamp,
                         yrend timestamp)''')
        for i,d in enumerate(zip(styear, endyr)):
            cur.execute(f'''INSERT into MESmercyear
                            values ({i}, '{d[0].iso}', '{d[1].iso}')''')
    else:
        pass

    if datatime is not None:
        yrnum = pd.read_sql('''SELECT * from MESmercyear''', con)
        myear = np.ndarray(len(datatime), dtype=int)
        for i,yr in yrnum.iterrows():
            q = (datatime > yr.yrstart) & (datatime < yr.yrend)
            myear[q] = yr.yrnum

        return myear
    else:
        return None

class MESSENGERdata():
    def __init__(self, species=None, comparisons=None):
        '''Retrieve MESSENGER data from database

        Species is required because each species is in different tables
        At least one SQL-formatted comparison is required. '''

        if (species is None):
            species = ''
        allspecies = ['Na', 'Ca', 'Mg']
        if species not in allspecies:
            print(f"Valid species are {', '.join(allspecies)}")
            return None

        database = get_database()
        con = psycopg2.connect(database=database)
        if comparisons is None:
            columns = pd.read_sql(
                f'''SELECT * from {species}uvvsdata, {species}pointing
                    WHERE 1=2''', con)
            print('Available fields are:')
            for col in columns.columns:
                print(f'\t{col}')
            return None
        else:
            query = f'''SELECT * from {species}uvvsdata, {species}pointing
                        WHERE unum=pnum and {comparisons}
                        ORDER BY unum'''
            try:
                data = pd.read_sql(query, con)
            except:
                print(query)
                assert 0, 'Problem with comparisons given.'

            if (len(data) > 0):
                self.species = species
                self.frame = data.frame[0]
                data.drop(['species', 'frame'], inplace=True, axis=1)
                self.data = data
            else:
                print(query)
                print('No data found')
                self.data = None

    @staticmethod
    def initialize(datapath, database='thesolarsystemmb'):
        '''Convert raw IDL sav files to pickles'''

        # Add to the database
        database = get_database()
        con = psycopg2.connect(database=database)
        con.autocommit = True
        cur = con.cursor()

        try:
            cur.execute('DROP table nauvvsdata')
            cur.execute('DROP table napointing')
            cur.execute('DROP table cauvvsdata')
            cur.execute('DROP table capointing')
            cur.execute('DROP table mguvvsdata')
            cur.execute('DROP table mgpointing')
        except:
            pass

        print('creating UVVS tables')
        spec = ['Ca', 'Na', 'Mg']
        for sp in spec:
            # Table with spectrum information
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
                               slit text)''')
        # Not including slit corners right now

        savfiles = glob.glob(datapath+'/*.sav')
        savfiles = sorted(savfiles)
        for oldfile in savfiles:
            realfile = oldfile.replace('.sav', '_temp.pkl')
            newfile = oldfile.replace('.sav', '.pkl')
            print('{}\n{}\n***'.format(oldfile, newfile))
            #data = io.readsav(oldfile, python_dict=True)
            data = pickle.load(open(realfile, 'rb'))

            kR = u.def_unit('kR', 1e3*u.R)
            Rmerc = u.def_unit('R_Mercury', mercury.radius)
            nm = u.def_unit('nm', 1e-9*u.m)

            npts = len(data['orb_num'])
            species = os.path.basename(oldfile).split('.')[0]

            # Determine UT for each spectrum
            t_iso = ['{}:{}:{}'.format('20' + time[0:2].decode('utf-8'),
                                       time[2:5].decode('utf-8'),
                                   time[6:].decode('utf-8'))
                     for time in data['step_utc_time']]
            UTC = Time(t_iso, format='yday')

            # Orbit number for each data spectrum
            orbit = np.array([int(o) for o in data['orb_num']])

            # determine Mercury year
            myear = merc_year(UTC)
            rmerc = (np.sqrt(np.sum(data['planet_sun_vector_tg']**2, axis=1)) *
                     u.km).to(u.AU)

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
                xyz[i,:] = np.dot(data['mso_rotation_matrix'][i,:,:],
                             data['planet_sc_vector_tg'][i,:])/mercury.radius.value
                bore[i,:] = np.dot(data['mso_rotation_matrix'][i,:,:],
                                 data['boresight_unit_vector_center_tg'][i,:])
                corn0[i,:] = np.dot(data['mso_rotation_matrix'][i,:,:],
                                    data['boresight_unit_vector_c1_tg'][i,:])
                corn1[i,:] = np.dot(data['mso_rotation_matrix'][i,:,:],
                                    data['boresight_unit_vector_c2_tg'][i,:])
                corn2[i,:] = np.dot(data['mso_rotation_matrix'][i,:,:],
                                    data['boresight_unit_vector_c3_tg'][i,:])
                corn3[i,:] = np.dot(data['mso_rotation_matrix'][i,:,:],
                                    data['boresight_unit_vector_c4_tg'][i,:])

            xcorner = np.array([corn0[:,0], corn1[:,0],
                                corn2[:,0], corn3[:,0]]).transpose()
            ycorner = np.array([corn0[:,1], corn1[:,1],
                                corn2[:,1], corn3[:,1]]).transpose()
            zcorner = np.array([corn0[:,2], corn1[:,2],
                                corn2[:,2], corn3[:,2]]).transpose()

            # Determine tangent point
            rr = np.linalg.norm(xyz, axis=1)
            t = -np.sum(xyz*bore, axis=1)
            tanpt = xyz + bore*t[:,np.newaxis]
            rtan = np.linalg.norm(tanpt, axis=1)

            slit = np.array(['Surface' if s == 0 else 'Atmospheric'
                             for s in data['slit']])
            obstype = np.array([str(ob).replace('b', '').replace("'", '').strip()
                       for ob in data['obs_typ']])

            # Add in the spectra
            spectra = data[species.lower()+'_rad_kr']*kR
            wavelength = data['wavelength']*nm
            raw = data['orig']*u.ct
            corrected = data['corr']*u.ct
            dark = data['dark']*u.ct
            solarfit = data['sol_fit']*u.ct

            ndata = pd.DataFrame({
                        'species': species,
                        'frame': 'MSO',
                        'UTC': UTC,
                        'orbit': orbit,
                        'merc_year': myear,
                        'TAA': data['true_anomaly']*np.pi/180.,
                        'rmerc': rmerc.value,
                        'drdt': data['rad_vel'],
                        'subslong': data['subsolar_longitude']*np.pi/180.,
                        'g': data['gvals']/u.s,
                        'radiance': radiance,
                        'sigma': sigma,
                        'x': xyz[:,0] * Rmerc,
                        'y': xyz[:,1] * Rmerc,
                        'z': xyz[:,2] * Rmerc,
                        'xbore': bore[:,0],
                        'ybore': bore[:,1],
                        'zbore': bore[:,2],
                        'xcorn1': xcorner[:,0],
                        'xcorn2': xcorner[:,1],
                        'xcorn3': xcorner[:,2],
                        'xcorn4': xcorner[:,3],
                        'ycorn1': ycorner[:,0],
                        'ycorn2': ycorner[:,1],
                        'ycorn3': ycorner[:,2],
                        'ycorn4': ycorner[:,3],
                        'zcorn1': zcorner[:,0],
                        'zcorn2': zcorner[:,1],
                        'zcorn3': zcorner[:,2],
                        'zcorn4': zcorner[:,3],
                        'obstype': obstype,
                        'obstype_num': data['obs_typ_num'],
                        'xtan': tanpt[:,0],
                        'ytan': tanpt[:,1],
                        'ztan': tanpt[:,2],
                        'rtan': rtan,
                        'alttan': data['target_altitude_set'][:,0],
                        'minalt': data['minalt'],
                        'longtan': data['target_longitude_set'][:,0]*np.pi/180,
                        'lattan': data['target_latitude_set'][:,0]*np.pi/180,
                        'loctimetan': data['obs_solar_localtime'],
                        'slit': slit})
            ndata.fillna(-999, inplace=True)

            spectra = {'spectra': spectra,
                       'wavelength': wavelength,
                       'corrected': corrected,
                       'dark': dark,
                       'solarfit': solarfit}

            # save this for later
            with open(newfile, 'wb') as f:
                pickle.dump(ndata, f, pickle.HIGHEST_PROTOCOL)
            with open(newfile.replace('.pkl', '_spectra.pkl'), 'wb') as f:
                pickle.dump(spectra, f, pickle.HIGHEST_PROTOCOL)

            print('Inserting UVVS data')
            for i,dpoint in ndata.iterrows():
                cur.execute(f'''INSERT into {species}uvvsdata (
                                    species, frame, UTC, orbit, merc_year,
                                    taa, rmerc, drdt, subslong, g, radiance,
                                    sigma) values (
                                    '{dpoint.species}',
                                    '{dpoint.frame}',
                                    '{dpoint.UTC.iso}',
                                    {dpoint.orbit},
                                    {dpoint.merc_year},
                                    {dpoint.TAA},
                                    {dpoint.rmerc},
                                    {dpoint.drdt},
                                    {dpoint.subslong},
                                    {dpoint.g},
                                    {dpoint.radiance},
                                    {dpoint.sigma})''')
                cur.execute(f'''INSERT into {species}pointing (
                                    x, y, z, xbore, ybore, zbore,
                                    obstype, obstype_num, xtan, ytan, ztan,
                                    rtan, alttan, longtan, lattan,
                                    loctimetan, slit) values (
                                    {dpoint.x},
                                    {dpoint.y},
                                    {dpoint.z},
                                    {dpoint.xbore},
                                    {dpoint.ybore},
                                    {dpoint.zbore},
                                    '{dpoint.obstype}',
                                    {dpoint.obstype_num},
                                    {dpoint.xtan},
                                    {dpoint.ytan},
                                    {dpoint.ztan},
                                    {dpoint.rtan},
                                    {dpoint.alttan},
                                    {dpoint.longtan},
                                    {dpoint.lattan},
                                    {dpoint.loctimetan},
                                    '{dpoint.slit}')''')

    def __str__(self):
        for a in self.__dict__:
            try:
                ll = self.__dict__[a].shape
            except:
                ll = ''
            print('{}\t{}\t{}'.format(a, type(self.__dict__[a]), ll))
        return ''

    def __len__(self):
        return len(self.data)

    def set_frame(self, frame=None):
        '''Convert between MSO and Model frames.

        More frames could be added if necessary.
        If Frame is not specified, flips between MSO and Model.
        '''

        if (frame is None) and (self.frame == 'MSO'):
            frame = 'Model'
        elif (frame is None) and (self.frame == 'Model'):
            frame = 'MSO'
        else:
            pass

        allframes = ['Model', 'MSO']
        if frame not in allframes:
            print('{} is not a valid frame.'.format(frame))
            return None
        elif frame == self.frame:
            pass
        elif (self.frame == 'MSO') and (frame == 'Model'):
            # Convert from MSO to Model
            self.data.x, self.data.y = self.data.y, -self.data.x
            self.data.xbore, self.data.ybore = (self.data.ybore,
                                                -self.data.xbore)
            # self.data.xcorner, self.data.ycorner = (self.data.ycorner,
            #                                         -self.data.xcorner)
            self.data.xtan, self.data.ytan = self.data.ytan, -self.data.xtan
            self.data.frame = 'Model'
        elif (self.data.frame == 'Model') and (frame == 'MSO'):
            self.data.x, self.data.y = -self.data.y, self.data.x
            self.data.xbore, self.data.ybore = (-self.data.ybore,
                                                self.data.xbore)
            # self.data.xcorner, self.data.ycorner = (-self.data.ycorner,
            #                                         self.data.xcorner)
            self.data.xtan, self.data.ytan = -self.data.ytan, self.data.xtan
            self.data.frame = 'MSO'
        else:
            assert 0, 'You somehow picked a bad combination.'

    def model(self, inputs_, npackets, quantity='radiance',
                                       dphi=3*u.deg,
                                       overwrite=False,
                                       filenames=None):
        from nexoclom import Input, modeldriver, LOSResult

        if isinstance(inputs_, str):
            inputs = Input(inputs_)
        elif isinstance(inputs_, Input):
            inputs = inputs_

        # TAA needs to match the data
        oldtaa = inputs.geometry.taa
        if len(set(self.data.orbit.tolist())) == 1:
            inputs.geometry.taa = np.median(self.data.taa)*u.rad
        elif np.max(self.data.taa)-np.min(self.data.taa) < 3*np.pi/180.:
            inputs.geometry.taa = np.median(self.data.taa)*u.rad
        else:
            assert 0, 'Too wide a range of taa'

        # Run the model
        modeldriver(inputs, npackets, overwrite=overwrite)

        # Simulate the data
        self.inputs = inputs
        self.set_frame('Model')
        self.modelresult = LOSResult(inputs, self.data, quantity,
                                     dphi=dphi, filenames=filenames)
        self.data['model'] = self.modelresult.radiance/1e3 # Convert to kR

        interval = PercentileInterval(50)
        lim = interval.get_limits(self.data.radiance)
        mask = ((self.data.radiance >= lim[0]) &
                (self.data.radiance <= lim[1]))

        strunit = u.def_unit('10**26 atoms/s', 1e26/u.s)
        m_data = np.mean(self.data.radiance[mask])
        m_model = np.mean(self.data.model[mask])
        self.modelstrength = m_data/m_model * strunit * strunit

        self.data.model = self.data.model * self.modelstrength.value

        # Put the old TAA back in.
        inputs.geometry.taa = oldtaa

    # def subset(self, q):
    #     ''' Reduces the object to a subset of the object.'''
    #
    #     self.UTC = self.UTC[q]
    #     self.orbit = self.orbit[q]
    #     self.mercyear = self.mercyear[q]
    #     self.taa = self.taa[q]
    #     self.rmerc = self.rmerc[q]
    #     self.drdt = self.drdt[q]
    #     self.subslong = self.subslong[q]
    #     self.g = self.g[q]
    #     self.radiance = self.radiance[q]
    #     self.sigma = self.sigma[q]
    #     self.x = self.x[q]
    #     self.y = self.y[q]
    #     self.z = self.z[q]
    #     self.xbore = self.xbore[q]
    #     self.ybore = self.ybore[q]
    #     self.zbore = self.zbore[q]
    #     self.xcorner = self.xcorner[q,:]
    #     self.ycorner = self.ycorner[q,:]
    #     self.zcorner = self.zcorner[q,:]
    #     self.obstype = self.obstype[q]
    #     self.obstype_num = self.obstype_num[q]
    #     self.xtan = self.xtan[q]
    #     self.ytan = self.ytan[q]
    #     self.ztan = self.ztan[q]
    #     self.rtan = self.rtan[q]
    #     self.alttan = self.alttan[q,:]
    #     self.minalt = self.minalt[q]
    #     self.longtan = self.longtan[q,:]
    #     self.lattan = self.lattan[q,:]
    #     self.loctimetan = self.loctimetan[q]
    #     self.slit = self.slit[q]
    #     self.wavelength = self.wavelength[q,:]
    #     self.spectra = self.spectra[q,:]
    #     self.raw = self.raw[q,:]
    #     self.corrected = self.corrected[q,:]
    #     self.dark = self.dark[q,:]
    #     self.solarfit = self.solarfit[q,:]
