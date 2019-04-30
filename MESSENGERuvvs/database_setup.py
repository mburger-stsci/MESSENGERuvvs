import os
import glob


def messenger_database_setup():
    configfile = os.path.join(os.environ['HOME'], '.nexoclom')
    config = {}
    if os.path.isfile(configfile):
        for line in open(configfile, 'r').readlines():
            key, value = line.split('=')
            config[key.strip()] = value.strip()

    savepath = config['savepath']
    datapath = config['datapath']
    database = config['database']

    datafiles = glob.glob(datapath+'/UVVS*sql')
    for dfile in datafiles:
        os.system(f'psql -d {database} -f {dfile}')

if __name__ == '__main__':
    # Load MESSENGER data into database
    cfmes = input('Load the MESSENGER data? (y/n) ')
    if cfmes.lower() in ('y', 'yes'):
        messenger_database_setup()
