# Get package info using pkg-config
find_package(PkgConfig)
pkg_search_module(CBC cbc)

include(canonicalize_paths)
canonicalize_paths(CBC_LIBRARY_DIRS)

# add osx frameworks to CBC_LIBRARIES
if(CBC_LIBRARIES)
  if(APPLE)
    # turn "-framework;foo;-framework;bar;other" into "-framework foo;-framework bar;other"
    string(REPLACE "-framework;" "-framework " CBC_LDFLAGS_OTHER "${CBC_LDFLAGS_OTHER}")
    # take "-framework foo;-framework bar;other" and add only frameworks to CBC_LIBRARIES
    foreach(arg ${CBC_LDFLAGS_OTHER})
      if("${arg}" MATCHES "-framework .+")
        set(CBC_LIBRARIES "${CBC_LIBRARIES};${arg}")
      endif("${arg}" MATCHES "-framework .+")
    endforeach(arg ${CBC_LDFLAGS_OTHER})
  endif(APPLE)
endif(CBC_LIBRARIES)

# Set standard flags
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(CBC DEFAULT_MSG CBC_LIBRARIES CBC_INCLUDE_DIRS)
