import os
import numpy as np
import matplotlib.pyplot as plt
import astropy.units as u
from nexoclom2 import Input, Output, path, SSObject
from MESSENGERuvvs import MESSENGERdata, MESSENGERModel
from astropy.visualization import quantity_support
from astropy.time import Time
quantity_support()


data = MESSENGERdata('Ca', 'orbit=206')

inputfile = os.path.join(os.path.dirname(path), 'tests',
                         'test_data', 'inputfiles',
                         'Mercury_Ca_spot_notime.input')
inputs = Input(inputfile)
inputs.geometry.taa = np.median(data.taa)
inputs.options.species = data.species

output = Output(inputs, 100000, overwrite=False)

t0 = Time.now()
model = MESSENGERModel(data, output)
t1 = Time.now()
print((t1-t0).to(u.s))


plt.scatter(data.utc.to_datetime(), data.radiance, color='black')
plt.plot(data.utc.to_datetime(), model.radiance, color='red')
plt.pause(1)


from shapely import Polygon

# far_pt = (data.x[0]
# cone = Polygon([(data.x[0], data.y[0], data.z[0]),


final = output.final_state(which=range(0, 10000))


from inspect import currentframe, getframeinfo
frameinfo = getframeinfo(currentframe())
print(frameinfo.filename, frameinfo.lineno)
from IPython import embed; embed()
import sys; sys.exit()
