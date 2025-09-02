.. _cmake_api_ref:

CMake API Reference
===================

To include alpaqa in your CMake project, use the ``find_package`` command:

.. code-block:: cmake

    find_package(alpaqa 1.1.0 [EXACT] [QUIET] [REQUIRED]
                 [COMPONENTS <components> ...]
                 [OPTIONAL_COMPONENTS <components> ...])


Components and targets
----------------------

alpaqa comes with multiple optional components, some of which can be packaged
and installed independently. The available components are:

+--------------+--------------------------------------------------------------------------------------------------------------+-----------------------------------------------------------------------------------------------------+
|  Component   |  Description                                                                                                 |  Targets                                                                                            |
+==============+==============================================================================================================+=====================================================================================================+
|  ``Core``    |  The main alpaqa solvers and other core functionality. If no components are specified, this is the default.  |  ``alpaqa``                                                                                         |
+--------------+--------------------------------------------------------------------------------------------------------------+-----------------------------------------------------------------------------------------------------+
|  ``Dl``      |  The dynamic loading C API headers for creating problems that can be loaded by the alpaqa solvers.           |  ``dl-api``                                                                                         |
+--------------+--------------------------------------------------------------------------------------------------------------+-----------------------------------------------------------------------------------------------------+
|  ``CasADi``  |  Classes for loading and building problem definitions using CasADi.                                          |  ``casadi-loader``, ``casadi-ocp-loader``                                                           |
+--------------+--------------------------------------------------------------------------------------------------------------+-----------------------------------------------------------------------------------------------------+
|  ``Extra``   |  Additional solvers and problem loaders that fall outside of the core library.                               |  ``dl-loader``, ``cutest-interface``, ``ipopt-adapter``, ``lbfgsb-adapter``, ``qpalm-adapter`` ...  |
+--------------+--------------------------------------------------------------------------------------------------------------+-----------------------------------------------------------------------------------------------------+
|  ``Tools``   |  Command-line tools for invoking the solvers.                                                                |  ``driver``, ``gradient-checker``                                                                   |
+--------------+--------------------------------------------------------------------------------------------------------------+-----------------------------------------------------------------------------------------------------+

Targets are prefixed with ``alpaqa::``. For example, to link with the main alpaqa library,
use:

.. code-block:: cmake

    target_link_libraries(<target> PUBLIC alpaqa::alpaqa)

To check whether a certain target is available, you can use:

.. code-block:: cmake

    if (TARGET alpaqa::qpalm-adapter)
        # ...
    endif()

Commands
--------

.. cmake-module::
    ../../../../src/cmake/dl-problem.cmake
