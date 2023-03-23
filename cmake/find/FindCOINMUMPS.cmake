# Distributed under the OSI-approved BSD 3-Clause License.  See accompanying
# file Copyright.txt or https://cmake.org/licensing for details.

#[=======================================================================[.rst:
FindCOINMUMPS
-------------

Finds the COIN-OR MUMPS library.

Imported Targets
^^^^^^^^^^^^^^^^

This module provides the following imported targets, if found:

``COIN::MUMPS``
  The COIN-OR MUMPS library

Result Variables
^^^^^^^^^^^^^^^^

This will define the following variables:

``COINMUMPS_FOUND``
  True if the system has the MUMPS library.
``COINMUMPS_VERSION``
  The version of the MUMPS library which was found.
``COINMUMPS_INCLUDE_DIRS``
  Include directories needed to use MUMPS.
``COINMUMPS_LIBRARIES``
  Libraries needed to link to MUMPS.

Cache Variables
^^^^^^^^^^^^^^^

The following cache variables may also be set:

``COINMUMPS_INCLUDE_DIR``
  The directory containing ``dmumps_c.h``.
``COINMUMPS_LIBRARY``
  The path to the MUMPS library ``libcoinmumps.a``.

#]=======================================================================]

find_package(PkgConfig)
pkg_check_modules(PC_MUMPS QUIET coinmumps)

find_path(COINMUMPS_INCLUDE_DIR
    NAMES dmumps_c.h
    PATHS ${PC_COINMUMPS_INCLUDE_DIRS}
    PATH_SUFFIXES coin-or/mumps
)
find_library(COINMUMPS_LIBRARY
    NAMES coinmumps
    PATHS ${PC_COINMUMPS_LIBRARY_DIRS}
)

if (COINMUMPS_INCLUDE_DIR)
    set(COINMUMPS_VERSION_REGEX " *#define +MUMPS_VERSION +\"([0-9.]*)\"")
    file(STRINGS ${COINMUMPS_INCLUDE_DIR}/dmumps_c.h COINMUMPS_VERSION_LINE
        REGEX ${COINMUMPS_VERSION_REGEX}
        LIMIT_COUNT 1)
    string(REGEX MATCH ${COINMUMPS_VERSION_REGEX}
        COINMUMPS_VERSION_REGEX_MATCHED ${COINMUMPS_VERSION_LINE})
    if (COINMUMPS_VERSION_REGEX_MATCHED)
        set(COINMUMPS_VERSION ${CMAKE_MATCH_1})
    else()
        set(COINMUMPS_VERSION ${PC_COINMUMPS_VERSION})
    endif()
endif()

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(COINMUMPS
    FOUND_VAR COINMUMPS_FOUND
    REQUIRED_VARS
        COINMUMPS_LIBRARY
        COINMUMPS_INCLUDE_DIR
    VERSION_VAR COINMUMPS_VERSION
)

if(COINMUMPS_FOUND)
    set(COINMUMPS_LIBRARIES ${COINMUMPS_LIBRARY})
    set(COINMUMPS_INCLUDE_DIRS ${COINMUMPS_INCLUDE_DIR})
    set(COINMUMPS_DEFINITIONS ${PC_COINMUMPS_CFLAGS_OTHER})
endif()

if(COINMUMPS_FOUND AND NOT TARGET COIN::MUMPS)
    add_library(COIN::MUMPS UNKNOWN IMPORTED)
    set_target_properties(COIN::MUMPS PROPERTIES
        IMPORTED_LOCATION "${COINMUMPS_LIBRARY}"
        INTERFACE_INCLUDE_DIRECTORIES "${COINMUMPS_INCLUDE_DIR}"
    )
    find_package(OpenBLAS REQUIRED)
    find_package(Threads REQUIRED)
    target_link_libraries(COIN::MUMPS INTERFACE Threads::Threads OpenBLAS::OpenBLAS)
    include(CheckLanguage)
    check_language(Fortran)
    if (CMAKE_Fortran_COMPILER)
        set(ALPAQA_HAVE_FORTRAN On)
        enable_language(Fortran)
    endif()
    set_target_properties(COIN::MUMPS PROPERTIES LINKER_LANGUAGE Fortran)
    if (CMAKE_Fortran_COMPILER_ID STREQUAL GNU)
        target_link_libraries(COIN::MUMPS INTERFACE gfortran)
    endif()
endif()

mark_as_advanced(
    COINMUMPS_INCLUDE_DIR
    COINMUMPS_LIBRARY
)

# compatibility variables
set(COINMUMPS_VERSION_STRING ${COINMUMPS_VERSION})
