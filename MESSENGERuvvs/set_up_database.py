import sqlalchemy as sqla
>>> from sqlalchemy.orm import DeclarativeBase

from MESSENGERuvvs.get_datapath import get_datapath


def set_up_database(l1files):
    datapath = get_datapath()
    
class Base(DeclarativeBase):
    pass
    
class
    
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
