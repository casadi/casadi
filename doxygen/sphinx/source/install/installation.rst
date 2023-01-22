Installation 
==============

Pip
---

The preferred way to install the alpaqa Python interface is using pip:

.. code-block:: sh

    python3 -m pip install alpaqa

(`PyPI <https://pypi.org/project/alpaqa>`_)

To compile problems using the Python interface, you will need a C compiler, such
as GCC or Clang on Linux and MSVC on Windows.

From source
-----------

Building alpaqa from source requires the installation of some C++ dependencies, 
see `Installation (Doxygen) <../../Doxygen/installation.html>`_ for detailed
instructions.

C++ Library
-----------

Pre-built binaries for Linux are available from the
`Releases page on GitHub <https://github.com/kul-optec/alpaqa/releases>`_.

For Debian-based systems, the .deb packages can be installed using

.. code-block:: sh

    sudo dpkg -i libalpaqa*_1.0.0a7_amd64.deb

Different components are available:

* ``libalpaqa`` contains the shared libraries needed to run applications that
  use alpaqa.
* ``libalpaqa-debug`` contains the debugging symbols for those libraries.
* ``libalpaqa-dl_dev`` contains the header files needed to compile problem
  specifications that can be dynamically loaded by alpaqa.
* ``libalpaqa-dev`` contains all development files such as headers and CMake
  configuration files needed to compile software that invokes alpaqa solvers.

Alternatively, the .tar.gz file can be extracted and installed manually.

.. code-block:: sh

    sudo tar xzf alpaqa-1.0.0a7-Linux-x86_64.tar.gz -C /usr/local --strip-components=1
