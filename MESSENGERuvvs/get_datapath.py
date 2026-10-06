import os


def get_datapath():
    configfile = os.environ['NEXOCLOMCONFIG']
    datapath = None
    for line in open(configfile):
        if 'messengerdir' in line.lower():
            datapath = line.split('=')[1].strip()
        else:
            pass

    assert datapath is not None, 'MESSENGERDIR not set in config file'
    return datapath
