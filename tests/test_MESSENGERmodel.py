import os
import numpy as np
import matplotlib.pyplot as plt
import astropy.units as u
from nexoclom2 import Input, Output, path, SSObject
from MESSENGERuvvs import MESSENGERdata, MESSENGERModel
from astropy.visualization import quantity_support
quantity_support()


data = MESSENGERdata('Ca', 'orbit=206')

inputfile = os.path.join(os.path.dirname(path), 'tests',
                         'test_data', 'inputfiles',
                         'Mercury_Ca_spot_notime.input')
inputs = Input(inputfile)
inputs.geometry.taa = np.median(data.taa)
inputs.options.species = data.species

output = Output(inputs, 100000, overwrite=False)

from inspect import currentframe, getframeinfo
frameinfo = getframeinfo(currentframe())
print(frameinfo.filename, frameinfo.lineno)
from IPython import embed; embed()
import sys; sys.exit()

final = output.final_state()
model = MESSENGERModel(data, output)

plt.scatter(data.utc.to_datetime(), data.radiance, color='black')
plt.plot(data.utc.to_datetime(), model.radiance, color='red')
plt.pause(1)

fig, ax = plt.subplots(figsize=(10, 10))
ax.set_aspect('equal')
ax.scatter(final.x, final.y, s=1)
ax.plot(data.x, data.y, color='black')
for i in range(0, len(data), 10):
    ax.plot([data.x[i].value, data.x[i].value+10*data.xbore[i]],
            [data.y[i].value, data.y[i].value+10*data.ybore[i]], color='blue')

from inspect import currentframe, getframeinfo
frameinfo = getframeinfo(currentframe())
print(frameinfo.filename, frameinfo.lineno)
from IPython import embed; embed()
import sys; sys.exit()
