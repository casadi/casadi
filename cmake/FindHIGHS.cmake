# Get package info using pkg-config
find_package(PkgConfig)
pkg_search_module(HIGHS highs)

include(canonicalize_paths)
canonicalize_paths(HIGHS_LIBRARY_DIRS)

# Set standard flags
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(HIGHS DEFAULT_MSG HIGHS_LIBRARIES HIGHS_INCLUDE_DIRS)

if(HIGHS_FOUND)
  add_library(highs INTERFACE)

  foreach(LIB ${HIGHS_LIBRARIES})
    find_library(LIB_FULL_${LIB} NAMES ${LIB} PATHS ${HIGHS_LIBRARY_DIRS})
    target_link_libraries(highs INTERFACE ${LIB_FULL_${LIB}})
  endforeach()
  target_include_directories(highs INTERFACE ${HIGHS_INCLUDE_DIRS})
endif()
