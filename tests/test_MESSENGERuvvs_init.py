from MESSENGERuvvs import MESSENGERdata


data = MESSENGERdata('Ca', 'orbit = 36', load_spectra=True)

from inspect import currentframe, getframeinfo
frameinfo = getframeinfo(currentframe())
print(frameinfo.filename, frameinfo.lineno)
from IPython import embed; embed()
import sys; sys.exit()
