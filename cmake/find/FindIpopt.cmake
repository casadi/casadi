# Distributed under the OSI-approved BSD 3-Clause License.  See accompanying
# file Copyright.txt or https://cmake.org/licensing for details.

#[=======================================================================[.rst:
FindIpopt
---------

Finds the COIN-OR Ipopt library.

Imported Targets
^^^^^^^^^^^^^^^^

This module provides the following imported targets, if found:

``COIN::Ipopt``
  The COIN-OR Ipopt library

Result Variables
^^^^^^^^^^^^^^^^

This will define the following variables:

``Ipopt_FOUND``
  True if the system has the Ipopt library.
``Ipopt_VERSION``
  The version of the Ipopt library which was found.
``Ipopt_INCLUDE_DIRS``
  Include directories needed to use Ipopt.
``Ipopt_LIBRARIES``
  Libraries needed to link to Ipopt.

Cache Variables
^^^^^^^^^^^^^^^

The following cache variables may also be set:

``Ipopt_INCLUDE_DIR``
  The directory containing ``IpoptConfig.h``.
``Ipopt_LIBRARY``
  The path to the Ipopt library ``libipopt.a``.

#]=======================================================================]

find_package(PkgConfig)
pkg_check_modules(PC_Ipopt QUIET ipopt)

find_path(Ipopt_INCLUDE_DIR
    NAMES IpoptConfig.h
    PATHS ${PC_Ipopt_INCLUDE_DIRS}
    PATH_SUFFIXES coin-or
)
find_library(Ipopt_LIBRARY
    NAMES ipopt
    PATHS ${PC_Ipopt_LIBRARY_DIRS}
)

if (Ipopt_INCLUDE_DIR)
    set(Ipopt_VERSION_REGEX " *#define +IPOPT_VERSION +\"([0-9.]*)\"")
    file(STRINGS ${Ipopt_INCLUDE_DIR}/IpoptConfig.h Ipopt_VERSION_LINE
        REGEX ${Ipopt_VERSION_REGEX}
        LIMIT_COUNT 1)
    string(REGEX MATCH ${Ipopt_VERSION_REGEX}
        Ipopt_VERSION_REGEX_MATCHED ${Ipopt_VERSION_LINE})
    if (Ipopt_VERSION_REGEX_MATCHED)
        set(Ipopt_VERSION ${CMAKE_MATCH_1})
    else()
        set(Ipopt_VERSION ${PC_Ipopt_VERSION})
    endif()
endif()

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(Ipopt
    FOUND_VAR Ipopt_FOUND
    REQUIRED_VARS
        Ipopt_LIBRARY
        Ipopt_INCLUDE_DIR
    VERSION_VAR Ipopt_VERSION
)

if(Ipopt_FOUND)
    set(Ipopt_LIBRARIES ${Ipopt_LIBRARY})
    set(Ipopt_INCLUDE_DIRS ${Ipopt_INCLUDE_DIR})
    set(Ipopt_DEFINITIONS ${PC_Ipopt_CFLAGS_OTHER})
endif()

if(Ipopt_FOUND AND NOT TARGET COIN::Ipopt)
    add_library(COIN::Ipopt UNKNOWN IMPORTED)
    set_target_properties(COIN::Ipopt PROPERTIES
        IMPORTED_LOCATION "${Ipopt_LIBRARY}"
        INTERFACE_INCLUDE_DIRECTORIES "${Ipopt_INCLUDE_DIR}"
    )
    find_package(COINMUMPS REQUIRED)
    target_link_libraries(COIN::Ipopt INTERFACE COIN::MUMPS ${CMAKE_DL_LIBS})
endif()

mark_as_advanced(
    Ipopt_INCLUDE_DIR
    Ipopt_LIBRARY
)

# compatibility variables
set(Ipopt_VERSION_STRING ${Ipopt_VERSION})
