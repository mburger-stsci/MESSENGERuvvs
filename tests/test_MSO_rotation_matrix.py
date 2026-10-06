import os
import numpy as np
import pandas as pd
import astropy.units as u
from astropy.time import Time
from scipy.io import readsav
import spiceypy as spice
import matplotlib.pyplot as plt
from matplotlib import rcParams
from matplotlib.backends.backend_pdf import PdfPages
from nexoclom2 import SSObject
from nexoclom2.solarsystem.load_kernels import SpiceKernels
from astropy.visualization import quantity_support
quantity_support()
rcParams.update({'lines.linewidth': 1})

mercury = SSObject('Mercury')
datapath = '/Users/mburger/Work/Data/MESSENGER/UVVS/Level1'
data = pd.read_pickle(os.path.join(datapath,
    'NA_summary_20110329_20120304_filtered_NM_v4.L1.pkl'))

datapath = '/Users/mburger/Work/Data/MESSENGER/UVVS/SummaryFiles/V0001'
rawdata = readsav(os.path.join(datapath,
    'NA_summary_20110329_20120304_filtered_NM_v4.sav'))

num = len(data)
utc = data.UTC.values[:num]
mso_matrix = rawdata['mso_rotation_matrix'][:num,:,:]
x_iau = (rawdata['planet_sc_vector_tg'][:num,:]*u.km).to(mercury.unit)
bore_iau = rawdata['boresight_unit_vector_center_tg'][:num,:]
x_mes = np.column_stack([data.x.values[:num],
                         data.y.values[:num],
                         data.z.values[:num]])*mercury.unit
x_mso = np.zeros_like(x_mes)
bore_mso = np.zeros((len(utc), 3))

kernels = SpiceKernels('Mercury')
for i in range(len(utc)):
    if utc[i].startswith('201'):
        et = spice.str2et(utc[i])
        R = spice.pxform('IAU_MERCURY', 'MERCURYSOLAR', et)
        
        x_mso[i,:] = np.matmul(R, x_iau[i,:])
        bore_mso[i,:] = np.matmul(R, bore_iau[i,:])
    else:
        print(i, utc[i])
        
t = -np.sum(x_mso*bore_mso, axis=1)
behind = t < 0
tanpt_mso = x_mso + bore_mso*t[:, np.newaxis]
rtan_mso = np.linalg.norm(tanpt_mso, axis=1)
alttan_mso = (rtan_mso - 1*mercury.unit).to(u.km)
lattan_mso = (np.arcsin(tanpt_mso[:,2]/rtan_mso)).to(u.deg)
longtan_mso = np.mod(np.arctan2(tanpt_mso[:, 1], tanpt_mso[:, 0]),
                     2*np.pi*u.rad).to(u.deg)

tanpt_mso[behind,:] = 1e30*mercury.unit
rtan_mso[behind] = 1e30*mercury.unit
alttan_mso[behind] = 1e30*u.km
lattan_mso[behind] = np.nan
longtan_mso[behind] = np.nan

below = alttan_mso < 0
alttan_mso[below] = 0*u.km
tanpt_mso[below, :] /= rtan_mso[below, np.newaxis].value
rtan_mso[below] = 1*mercury.unit

alttan_mes = data.alttan.values[:num]*u.km
lattan_mes = (data.lattan.values[:num]*u.rad).to(u.deg)
longtan_mes = (data.longtan.values[:num]*u.rad).to(u.deg)

with  PdfPages('SummaryFile_vs_SPICE.pdf') as pdf:
    # x, y, z, alttan, longtan, lattan
    fig, ax = plt.subplots(3, 2, figsize=(16, 10))
    ax[0,0].plot(x_mes[:,0], label='Summary File')
    ax[0,0].plot(x_mso[:,0], label='SPICE')
    ax[0,0].set_title('Tanpoint x')

    ax[1,0].plot(x_mes[:,1], label='Summary File')
    ax[1,0].plot(x_mso[:,1], label='SPICE')
    ax[1,0].set_title('Tanpoint y')

    ax[2,0].plot(x_mes[:,2], label='Summary File')
    ax[2,0].plot(x_mso[:,2], label='SPICE')
    ax[2,0].set_title('Tanpoint z')

    q = (longtan_mes >= 0*u.deg) & (longtan_mes < 360*u.deg)
    ax[0,1].plot(longtan_mes[q], label='Summary File')
    ax[0,1].plot(longtan_mso[q], label='SPICE')
    ax[0,1].set_title('Tanpoint Subsolar Longitude')

    ax[1,1].plot(lattan_mes[q], label='Summary File')
    ax[1,1].plot(lattan_mso[q], label='SPICE')
    ax[1,1].set_title('Tanpoint Subsolar Latitude')

    q = alttan_mes <= 1e5*u.km
    ax[2,1].plot(alttan_mes[q], label='Summary File')
    ax[2,1].plot(alttan_mso[q], label='SPICE')
    ax[2,1].set_title('Tanpoint Altitude')

    ax[0,0].legend()
    pdf.savefig()
    plt.close()

    fig, ax2 = plt.subplots(3, 2, figsize=(16, 10))
    ax2[0,0].plot(x_mes[:,0]-x_mso[:,0], label='Summary File - SPICE')
    ax2[0,0].set_title(r'Tanpoint $\Delta$x')

    ax2[1,0].plot(x_mes[:,1]-x_mso[:,1], label='Summary File - SPICE')
    ax2[1,0].set_title(r'Tanpoint $\Delta$y')

    ax2[2,0].plot(x_mes[:,2]-x_mso[:,2], label='Summary File - SPICE')
    ax2[2,0].set_title(r'Tanpoint $\Delta$z')

    q = (longtan_mes >= 0*u.deg) & (longtan_mes < 360*u.deg)
    ax2[0,1].plot(longtan_mes[q]-longtan_mso[q], label='Summary File - Spice')
    ax2[0,1].set_title(r'Tanpoint $\Delta$ Subsolar Longitude')

    ax2[1,1].plot(lattan_mes[q]-lattan_mso[q], label='Summary File - Spice')
    ax2[1,1].set_title(f'Tanpoint $\Delta$ Subsolar Latitude')

    q = alttan_mes <= 1e5*u.km
    ax2[2,1].plot(alttan_mes[q]-alttan_mso[q], label='Summary File - Spice')
    ax2[2,1].set_title(f'Tanpoint $\Delta$ Altitude')

    ax2[0,0].legend()
    pdf.savefig()
    plt.close
