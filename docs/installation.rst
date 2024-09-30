.. _installation:
.. index:: Installation

Installation
============

:program:`PyPWDFT` is available via both Anaconda as well as PyPi. If you are
using Windows and have little experience with Python, we recommend to use
Anaconda. Anaconda is a complete Python package containing many modules and
useful programs. Anaconda can be obtained `via this link
<https://www.anaconda.com/download>`_.

.. note::

	If you run into trouble installing :program:`PyPWDFT` in Anaconda, an easy
	solution can be to create a separate environment for Anaconda. Please
	consult the "Troubleshooting" section as seen below.

Anaconda
--------

Open a Anaconda command prompt and run the following command:

.. code:: bash

	conda install -c ifilot pypwdft mendeleev pyfftw

PyPi
----

Open a terminal and run the following command:

.. code:: bash

	pip install pypwdft mendeleev pyfftw

Testing
-------

To test that your installation is working, you can run the following snippet
of code

.. code:: python

	import pypwdft
	print(pypwdft.__version__)

The version number should be returned.