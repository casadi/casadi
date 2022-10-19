#[=======================================================================[.rst:
FindqpOASES
--------

Find the qpOASES includes and library.

IMPORTED Targets
^^^^^^^^^^^^^^^^

This module defines IMPORTED target ``qpOASES::qpOASES``, if
qpOASES has been found.

#]=======================================================================]

# Try each search configuration.
find_path(qpOASES_INCLUDE_DIR NAMES qpOASES.h PATH_SUFFIXES include)

# Allow qpOASES_LIBRARY to be set manually, as the location of the qpOASES library
if(NOT qpOASES_LIBRARY)
  find_library(qpOASES_LIBRARY_RELEASE NAMES qpOASES)
  find_library(qpOASES_LIBRARY_DEBUG NAMES qpOASESd)

  include(${CMAKE_CURRENT_LIST_DIR}/SelectLibraryConfigurations.cmake)
  select_library_configurations(qpOASES)
endif()

unset(qpOASES_NAMES)
unset(qpOASES_NAMES_DEBUG)

mark_as_advanced(qpOASES_INCLUDE_DIR)

include(${CMAKE_CURRENT_LIST_DIR}/FindPackageHandleStandardArgs.cmake)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(qpOASES REQUIRED_VARS qpOASES_LIBRARY qpOASES_INCLUDE_DIR
                                          HANDLE_COMPONENTS)

if(qpOASES_FOUND)
    if(NOT TARGET qpOASES::qpOASES)
      add_library(qpOASES::qpOASES UNKNOWN IMPORTED)
      set_target_properties(qpOASES::qpOASES PROPERTIES
        INTERFACE_INCLUDE_DIRECTORIES "${qpOASES_INCLUDE_DIR}")

      if(qpOASES_LIBRARY_RELEASE)
        set_property(TARGET qpOASES::qpOASES APPEND PROPERTY
          IMPORTED_CONFIGURATIONS RELEASE)
        set_target_properties(qpOASES::qpOASES PROPERTIES
          IMPORTED_LOCATION_RELEASE "${qpOASES_LIBRARY_RELEASE}")
      endif()

      if(qpOASES_LIBRARY_DEBUG)
        set_property(TARGET qpOASES::qpOASES APPEND PROPERTY
          IMPORTED_CONFIGURATIONS DEBUG)
        set_target_properties(qpOASES::qpOASES PROPERTIES
          IMPORTED_LOCATION_DEBUG "${qpOASES_LIBRARY_DEBUG}")
      endif()

      if(NOT qpOASES_LIBRARY_RELEASE AND NOT qpOASES_LIBRARY_DEBUG)
        set_property(TARGET qpOASES::qpOASES APPEND PROPERTY
          IMPORTED_LOCATION "${qpOASES_LIBRARY}")
      endif()
    endif()
endif()
