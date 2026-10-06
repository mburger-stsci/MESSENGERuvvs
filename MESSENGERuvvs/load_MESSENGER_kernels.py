import os
import spiceypy as spice
import glob
import requests
from nexoclom2.solarsystem.load_kernels import SpiceKernels
from nexoclom2.utilities.NexoclomConfig import NexoclomConfig


class MESSENGER_kernels(SpiceKernels):
    def __init__(self):
        super().__init__('Mercury')
        
        naif_url = 'https://naif.jpl.nasa.gov/pub/naif/pds/data/mess-e_v_h-spice-6-v1.0/messsp_1000/data/'
        
        datapath = NexoclomConfig().savepath
        kernelpath = os.path.join(datapath, 'spice_kernels')
        
        # MESSENGER Spacecraft Clock
        sclk_kernel = os.path.join(kernelpath, 'messenger_2548.tsc')
        if not os.path.exists(sclk_kernel):
            print('Downloading MESSENGER clock kernel')
            sclk_url = os.path.join(naif_url, 'sclk', 'messenger_2548.tsc')
            page = requests.get(sclk_url)
            with open(sclk_kernel, 'w') as file:
                file.write(page.content.decode('utf-8'))
        else:
            pass
        self.kernels.append(sclk_kernel)
        spice.furnsh(sclk_kernel)
        
        # MESSENGER instrument frames
        fk_kernel = os.path.join(kernelpath, 'msgr_v231.tf')
        if not os.path.exists(fk_kernel):
            print('Downloading Mercury instrument frames kernel')
            fk_url = os.path.join(naif_url, 'fk', 'msgr_v231.tf')
            page = requests.get(fk_url)
            with open(fk_kernel, 'w') as file:
                file.write(page.content.decode('utf-8'))
        else:
            pass
        self.kernels.append(fk_kernel)
        spice.furnsh(fk_kernel)
        
        # Spacecraft pointing kernels
        ck_kernels = glob.glob(os.path.join(datapath, 'spice_kernels', '*.bc'))
        self.kernels.extend(ck_kernels)
        spice.furnsh(ck_kernels)

        # Instrument kernel
        ik_kernel = os.path.join(kernelpath, 'msgr_mascs_v100.ti')
        self.kernels.append(ik_kernel)
        spice.furnsh(ik_kernel)
