Welcome to MESSENGERuvvs's documentation
========================================

MESSENGERuvvs provides an interface to reduced MESSENGER ( MErcury Surface,
Space ENvironment, GEochemistry, and Ranging) MASCS (Mercury Atmospheric and
Surface Composition Spectrometer) UVVS (UltraViolet and Visible Spectrometer)
data and to the nexoclom Neutral Exospheres and Cloud Model. See McClintock
et al. Space Science Reviews, 131, 481-521, 2007, and Solomon et al., Mercury -
The View After MESSENGER, Cambridge University Press, 2019, Ch 14-15, for
more details on the MESSENGER observations and models.

.. toctree::
  :maxdepth: 2

  MESSENGERuvvs/MESSENGERdata
  MESSENGERuvvs/database_setup
  MESSENGERuvvs/databasebackups

Installation
==================
MESSENGERuvvs can be installed with pip:
::

    $ pip install MESSENGERuvvs

The MESSENGERuvvs data is not provided with this software package. I was not
responsible for reducing it, and I'm not sure who I'm allowed to distribute it
to. Please contact me if you are interested in this. MESSENGER MASCS/UVVS data
is available from the `Planetary Data System
<http://pds-geosciences.wustl.edu/missions/messenger/mascs.htm>`_, although
in a different format. We can probably work on a way to get that data in the
format used here if necessary.

Reporting Issues
================
This project is hosted on github at `MESSENGERuvvs
<https://github.com/mburger-stsci/MESSENGERuvvs>`_. Please report bugs or make
comments there.

Contributing
============
Please let me know if you would like to make contributions.

:Authors: Matthew Burger
:License: :doc:`LICENSE`
