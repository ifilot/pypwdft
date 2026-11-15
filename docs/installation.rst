.. _installation:
.. index:: Installation

Installation
============

PyPWDFT requires Python 3.9 or newer and is available from the Python Package
Index (PyPI). You will need a terminal and a working Python installation.

Create a virtual environment
----------------------------

.. important::

   We strongly recommend installing PyPWDFT in a dedicated virtual
   environment. A virtual environment keeps PyPWDFT and its dependencies
   separate from your operating system and other Python projects. This avoids
   version conflicts, makes upgrades safer, and lets you remove the complete
   installation without affecting anything else.

Create an environment named ``.venv`` in your project directory:

.. code:: bash

   python -m venv .venv

Activate it on Linux or macOS:

.. code:: bash

   source .venv/bin/activate

On Windows PowerShell, activate it with:

.. code:: powershell

   .venv\Scripts\Activate.ps1

On the Windows Command Prompt, use:

.. code:: batch

   .venv\Scripts\activate.bat

Your terminal prompt will normally show ``(.venv)`` while the environment is
active. Run ``deactivate`` when you want to leave it. When you return to the
project later, activate the environment again before using PyPWDFT.

Install PyPWDFT
---------------

With the virtual environment active, upgrade ``pip`` and install PyPWDFT:

.. code:: bash

   python -m pip install --upgrade pip
   python -m pip install pypwdft

Using ``python -m pip`` ensures that the package is installed for the same
Python interpreter that created and runs the virtual environment. PyPWDFT's
required dependencies are installed automatically.

.. _gpu_installation:

GPU-compatible installation
---------------------------

PyPWDFT can perform calculations on an NVIDIA GPU through its optional CuPy
backend. A GPU installation requires:

* a CUDA-capable NVIDIA GPU;
* a compatible NVIDIA driver; and
* a CuPy package matching the CUDA major version.

Use Python 3.10 or newer for a new GPU environment because current CuPy
releases no longer support Python 3.9. First follow the virtual-environment and
PyPWDFT installation steps above. You can check that the NVIDIA driver detects
your GPU with:

.. code:: bash

   nvidia-smi

Then install exactly one of the following CuPy packages in the active virtual
environment. For CUDA 12.x, run:

.. code:: bash

   python -m pip install "cupy-cuda12x[ctk]"

For CUDA 13.x, run:

.. code:: bash

   python -m pip install "cupy-cuda13x[ctk]"

The ``ctk`` extra installs the required CUDA components into the virtual
environment; a compatible NVIDIA driver must still be installed on the host.
If a system-wide CUDA Toolkit is already installed, the extra can be omitted,
for example ``python -m pip install cupy-cuda12x``.

.. warning::

   Do not install more than one CuPy variant in the same environment. Packages
   such as ``cupy``, ``cupy-cuda12x``, and ``cupy-cuda13x`` conflict with one
   another. Consult the `official CuPy installation guide
   <https://docs.cupy.dev/en/stable/install.html>`_ if you are unsure which
   package is compatible with your driver and CUDA setup.

Verify that CuPy can access the GPU:

.. code:: bash

   python -c "import cupy; print(cupy.cuda.runtime.getDeviceCount())"

This command should print ``1`` or a larger number. To run a PyPWDFT
calculation on the GPU, select the CUDA device:

.. code:: python

   from pypwdft import PWDFT, Structure

   structure = Structure.from_name("H2", cell=10)
   calculation = PWDFT(
       structure,
       cutoff=40,
       xc="pbe",
       pseudopotential="gth",
       device="cuda",
   )
   result = calculation.run()

PyPWDFT automatically selects the CuPy FFT backend when ``device="cuda"``.
See :ref:`gpu_calculations` for more calculation options.

Verify the installation
-----------------------

Check that PyPWDFT can be imported and print its installed version:

.. code:: bash

   python -c "import pypwdft; print(pypwdft.__version__)"

If this command prints a version number without an error, the installation is
ready. Continue with :doc:`getting_started` to run your first calculation.

Troubleshooting
---------------

If Python reports that ``pypwdft`` cannot be found, confirm that the virtual
environment is active, then compare the interpreter and installer locations:

.. code:: bash

   python -c "import sys; print(sys.executable)"
   python -m pip --version

Both paths should refer to the same virtual environment. If creating a virtual
environment fails on a Debian- or Ubuntu-based system, install the
``python3-venv`` package supplied by the operating system and try again.

If CuPy cannot initialize CUDA, first run ``nvidia-smi`` to confirm that the
GPU and driver are available, then check the detected CUDA configuration:

.. code:: bash

   python -c "import cupy; cupy.show_config()"

An initialization error usually indicates an incompatible driver, CUDA
runtime, or CuPy package rather than a PyPWDFT installation problem.

To start over, deactivate the environment and delete its ``.venv`` directory,
then repeat the steps above. A dedicated environment makes this safe because
it contains no system-wide packages.
