# Get package info using pkg-config
find_package(PkgConfig)
pkg_search_module(BONMIN bonmin)

include(canonicalize_paths)
canonicalize_paths(BONMIN_LIBRARY_DIRS)

# add osx frameworks to BONMIN_LIBRARIES
if(BONMIN_LIBRARIES)
  if(APPLE)
    # turn "-framework;foo;-framework;bar;other" into "-framework foo;-framework bar;other"
    string(REPLACE "-framework;" "-framework " BONMIN_LDFLAGS_OTHER "${BONMIN_LDFLAGS_OTHER}")
    # take "-framework foo;-framework bar;other" and add only frameworks to BONMIN_LIBRARIES
    foreach(arg ${BONMIN_LDFLAGS_OTHER})
      if("${arg}" MATCHES "-framework .+")
        set(BONMIN_LIBRARIES "${BONMIN_LIBRARIES};${arg}")
      endif("${arg}" MATCHES "-framework .+")
    endforeach(arg ${BONMIN_LDFLAGS_OTHER})
  endif(APPLE)
endif(BONMIN_LIBRARIES)

# Set standard flags
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(BONMIN DEFAULT_MSG BONMIN_LIBRARIES BONMIN_INCLUDE_DIRS)

if(BONMIN_FOUND)
  add_library(bonmin INTERFACE)

  foreach(LIB ${BONMIN_LIBRARIES})
    find_library(LIB_FULL_${LIB} NAMES ${LIB} PATHS ${BONMIN_LIBRARY_DIRS})
    target_link_libraries(bonmin INTERFACE ${LIB_FULL_${LIB}})
  endforeach()
  target_include_directories(bonmin INTERFACE ${IPOPT_INCLUDE_DIRS})
  target_include_directories(bonmin INTERFACE ${BONMIN_INCLUDE_DIRS})
  if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")
  else()
    target_compile_definitions(bonmin INTERFACE HAVE_CSTDDEF)
  endif()
  target_compile_options(bonmin INTERFACE ${BONMIN_CFLAGS_OTHER})
endif()
