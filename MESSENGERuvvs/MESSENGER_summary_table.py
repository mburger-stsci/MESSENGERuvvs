"""Create summary of MESSENGER orbits"""
import numpy as np
import sqlalchemy as sqla
import pandas as pd
from nexoclom.utilities import NexoclomConfig


config = NexoclomConfig(verbose=False)
config.verify_database_running()
engine = config.create_engine(config.mesdatabase)

def create_MESSSENGER_summary():
    # Create the database table
    metadata_obj = sqla.MetaData()
    summary = sqla.Table('uvvs_summary', metadata_obj,
                         sqla.Column('orbit', sqla.Integer),
                         sqla.Column('species', sqla.String(2)),
                         sqla.Column('ut_start', sqla.DateTime),
                         sqla.Column('ut_end', sqla.DateTime),
                         sqla.Column('year', sqla.Integer),
                         sqla.Column('taa', sqla.Float),
                         sqla.Column('r_sun', sqla.Float),
                         sqla.Column('drdt', sqla.Float),
                         sqla.Column('subslong', sqla.Float),
                         sqla.Column('g', sqla.Float),
                         sqla.Column('nspec', sqla.Integer))
    # summary.delete(engine, checkfirst=True)
    summary.create(engine, checkfirst=True)
    
    species = ['Na', 'Ca', 'Mg']
    for sp in species:
        table = sqla.Table(f'{sp.lower()}uvvsdata', metadata_obj, autoload_with=engine)
        
        orbit_query = sqla.select(table.columns.orbit)
        with engine.connect() as con:
            results = pd.DataFrame(con.execute(orbit_query))
        orbits = sorted(results.orbit.unique())
        
        for orbit in orbits:
            print(sp, orbit)
            info_query = sqla.select(table.columns.utc,
                                     table.columns.merc_year,
                                     table.columns.taa,
                                     table.columns.rmerc,
                                     table.columns.drdt,
                                     table.columns.subslong,
                                     table.columns.g).where(
                table.columns.orbit == int(orbit))
            with engine.connect() as con:
                info = pd.DataFrame(con.execute(info_query))
                
            insert = summary.insert().values(
                orbit = int(orbit),
                species = sp,
                ut_start = info.utc.min(),
                ut_end = info.utc.max(),
                year = int(np.ceil(np.median(info.merc_year))),
                taa = np.median(info.taa),
                r_sun = info.rmerc.mean(),
                drdt = info.drdt.mean(),
                subslong = info.subslong.mean(),
                g = info.g.mean(),
                nspec = len(info))
            with engine.connect() as con:
                con.execute(insert)
                con.commit()
                
def query_summary_table(query):
    with engine.connect() as con:
        result = pd.read_sql(sqla.text(query), con)
    
    return result
    

if __name__ == '__main__':
    create_MESSSENGER_summary()
