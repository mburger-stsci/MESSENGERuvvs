.. _database_fields_:

*************************************
Database Fields Used by MESSENGERuvvs
*************************************

Two tables are created for each species (Ca, Na, and Mg): *xxuvvsdata* and
*xxuvvspointing* where xx is the species. Fields in each table:

cauvvsdata, nauvvsdata, mguvvsdata
==================================

unum
    Unique spectrum identifier. This is used to join xxuvvsdata and
    xxuvvspointing (xxuvvsdata.unum = xxpoainting.pnum)

species
    Species in the table. This is uni-valued for each table.

frame
    Reference frame for the pointing information. In the database this value
    is 'MSO for all spectra, but in can be changed to 'model' after being
    loaded into a MESSENGERdata object.

UTC
    UTC start time of the observation.

orbit
    Orbit number for the observation.

merc_year
    Mercury year since begining of mission. Years are measured from
    perihelion, so Year 1 is incomplete.

taa
    Mercury true anomaly angle in radians.

rmerc
    Distance of Mercury from the Sun in AU.

drdt
    Mercury's radial distance from the Sun in km/s.

subslong
    Sub-solar longitude in radians.

g
    g-value for Mercury relative to the Sun in photons/s.

radiance
    Measured radiance in the spectrum in kR.

sigma
    Uncertainty in radiance measurement in kR.

capointing, napointing, mgpointing
===================================

pnum
    Unique spectrum identifier. This is used to join xxuvvsdata and
    xxuvvspointing (xxuvvsdata.unum = xxpoainting.pnum)

x, y, z
    Spacecraft coordinates in the current reference frame in Mercury radii
    (Measured relative to Mercury's center).

xbore, ybore, zbore
    Boresight direction (:math:`\sqrt{xbore^2 + ybore^2 + zbore^2} = 1`)

obstype
    Type of commanded UVVS observation (e.g., UVVSExoScan).

obstype_num
    Type number corresponding to obsype (e.g., UVVSExoScan=5).

xtan, ytan, ztan
    Location of the point at which the instrument line-of-sight is tangent
    to Mercury's surface in Mercury radii (Measured relative to Mercury's
    center).

rtan
    Distance of the instrument line-of-sight tangent point from Mercury's
    center (:math:`rtan = \sqrt{xtan^2 + ytan^2 + ztan^2}`)

alttan
    Height (altitude) of the line-of-sight tangent point from Mercury's
    surface in km (:math:`alttan = rtan - 1 \mathrm{R_M}`).

longtan
    Longitude of the tangent point in radians. I'm not sure what coordinate
    system this is measured in (suface fixed or MSO).

lattan
    Latitude of the tangent point in radians.

loctimetan
    Local time of the tangent point in hours.

slit
    UVVS slit used for the observation (Atmospheric or Surface)

