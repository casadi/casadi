.. _installation:

Installation
============

Pip
---

The preferred way to install the alpaqa Python interface is using pip:

.. code-block:: sh

    python3 -m pip install --upgrade --pre 'alpaqa[casadi]'

(`PyPI <https://pypi.org/project/alpaqa>`_)

To compile problems using the Python interface, you will need a C compiler, such
as GCC or Clang on Linux, Xcode on macOS, or MSVC on Windows (see the
:ref:`Tips and tricks` page for details and alternatives).

From source
-----------

Building alpaqa from source requires the installation of some C++ dependencies, 
see `Installation (Doxygen) <../../Doxygen/installation.html>`_ for detailed
instructions.

C++ Library
-----------

The easiest way to install and use the alpaqa C++ libraries is by using Conan,
as explained in the
`CMake example <https://github.com/kul-optec/alpaqa/tree/develop/examples/CMake/Solver>`_.

Alternatively, pre-built binaries for Linux are available from the
`Releases page on GitHub <https://github.com/kul-optec/alpaqa/releases>`_.

For Debian-based systems, the .deb packages can be installed using

.. code-block:: sh

    sudo apt update
    sudo apt install ./libalpaqa*_1.1.0a1_amd64.deb

Different components are available:

* ``libalpaqa`` contains the shared libraries needed to run applications that
  use alpaqa.
* ``libalpaqa-debug`` contains the debugging symbols for those libraries.
* ``libalpaqa-dl_dev`` contains the C header files needed to compile problem
  specifications that can be dynamically loaded by alpaqa.
* ``libalpaqa-dev`` contains all development files such as headers and CMake
  configuration files needed to compile software that invokes alpaqa solvers.
* ``libalpaqa-extra`` contains additional solvers and problem loaders that fall
  outside of the core library.
* ``libalpaqa-extra-dev`` contains all development files for the extra
  libraries.
* ``libalpaqa-casadi`` contains classes for loading and building problem
  definitions using CasADi.
* ``libalpaqa-casadi-dev`` contains classes development files for those classes.
* ``libalpaqa-tools`` contains command line utilities such as alpaqa-driver,
  which can be used to invoke the solvers directly, without the need to write
  any C++ code.

More information about the different components of alpaqa can be found in the
:ref:`cmake_api_ref`.

The following distributions are tested:

* Debian: 11 (Bullseye), 12 (Bookworm), 13 (Trixie), Sid
* Ubuntu: 20.04 (Focal), 22.04 (Jammy), 24.04 (Noble), rolling

Alternatively, the .tar.gz file can be extracted and installed manually.

.. code-block:: sh

    sudo tar xzf alpaqa-1.1.0a1-Linux-x86_64.tar.gz -C /usr/local --strip-components=1

When using the development packages, it is important to use the correct version
of Eigen (currently 3.4.0) to avoid ABI incompatibilities.
You should also add the correct architecture-specific flags: ``-march=haswell``
for amd64 packages, and ``-mcpu=cortex-a53+crc+simd`` for aarch64 packages.
The reason is that Eigen's ABI depends on the selected SIMD ISA extensions.
ABI-related issues can be fully avoided by simply using Conan to install alpaqa
instead of using the pre-built packages.

The development packages of alpaqa depend on the guanaqo library, which can be
installed from source:

.. code-block:: sh

    git clone https://github.com/tttapa/guanaqo --branch=1.0.0-alpha.16
    cmake -B build-guanaqo -S guanaqo -G "Ninja Multi-Config" -DBUILD_TESTING=Off
    cmake --build build-guanaqo --config Debug
    sudo cmake --install build-guanaqo --config Debug --prefix=/usr/local
    cmake --build build-guanaqo --config Release
    sudo cmake --install build-guanaqo --config Release --prefix=/usr/local

MATLAB interface
----------------

The experimental MATLAB/MEX interface can be installed by running the following
command in the MATLAB command window:

.. code-block:: matlab

    unzip(['https://github.com/kul-optec/alpaqa/releases/download/1.1.0a1/alpaqa-matlab-' computer('arch') '.zip'], userpath)

You need CasADi to be installed as well: https://web.casadi.org/get

To uninstall the alpaqa MATLAB interface, simply remove the ``+alpaqa`` folder:

.. code-block:: matlab

    rmdir(fullfile(userpath, '+alpaqa'), 's');
