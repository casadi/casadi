# This file was taken from Sarcasm/irony-mode and modified
#
# Try to find libclang
#
# Once done this will define:
# - LIBCLANG_FOUND
#               System has libclang.
# - LIBCLANG_INCLUDE_DIRS
#               The libclang include directories.
# - LIBCLANG_LIBRARIES
#               The libraries needed to use libclang.
# - LIBCLANG_LIBRARY_DIR
#               The path to the directory containing libclang.
# - LIBCLANG_KNOWN_LLVM_VERSIONS
#               Known LLVM release numbers.

# most recent versions come first
# http://llvm.org/apt/
set(LIBCLANG_KNOWN_LLVM_VERSIONS 3.7.0 3.7
  3.6.2 3.6.1 3.6.0 3.6
  3.5.2 3.5.1 3.5.0 3.5
  3.4.2 3.4.1 3.4
  3.3
  3.2
  3.1)

set(libclang_llvm_header_search_paths)
set(libclang_llvm_lib_search_paths
  # LLVM Fedora
  /usr/lib/llvm
  )

foreach (version ${LIBCLANG_KNOWN_LLVM_VERSIONS})
  string(REPLACE "." "" undotted_version "${version}")
  list(APPEND libclang_llvm_header_search_paths
    # LLVM Debian/Ubuntu nightly packages: http://llvm.org/apt/
    "/usr/lib/llvm-${version}/include/"
    # LLVM MacPorts
    "/opt/local/libexec/llvm-${version}/include"
    # LLVM Homebrew
    "/usr/local/Cellar/llvm/${version}/include"
    # LLVM Homebrew/versions
    "/usr/local/lib/llvm-${version}/include"
    # FreeBSD ports versions
    "/usr/local/llvm${undotted_version}/include"
    )

  list(APPEND libclang_llvm_lib_search_paths
    # LLVM Debian/Ubuntu nightly packages: http://llvm.org/apt/
    "/usr/lib/llvm-${version}/lib/"
    # LLVM MacPorts
    "/opt/local/libexec/llvm-${version}/lib"
    # LLVM Homebrew
    "/usr/local/Cellar/llvm/${version}/lib"
    # LLVM Homebrew/versions
    "/usr/local/lib/llvm-${version}/lib"
    # FreeBSD ports versions
    "/usr/local/llvm${undotted_version}/lib"
    )
endforeach()

find_path(LIBCLANG_INCLUDE_DIR clang-c/Index.h
  PATHS ${libclang_llvm_header_search_paths}
  PATH_SUFFIXES LLVM/include #Windows package from http://llvm.org/releases/
  DOC "The path to the directory that contains clang-c/Index.h")

# On Windows with MSVC, the import library uses the ".imp" file extension
# instead of the comon ".lib"
if (MSVC)
  find_file(LIBCLANG_LIBRARY libclang.imp
    PATH_SUFFIXES LLVM/lib
    DOC "The file that corresponds to the libclang library.")
endif()

find_library(LIBCLANG_LIBRARY NAMES libclang.imp libclang clang
  PATHS ${libclang_llvm_lib_search_paths}
  PATH_SUFFIXES LLVM/lib #Windows package from http://llvm.org/releases/
  DOC "The file that corresponds to the libclang library.")

get_filename_component(LIBCLANG_LIBRARY_DIR ${LIBCLANG_LIBRARY} PATH)

set(LIBCLANG_LIBRARIES ${LIBCLANG_LIBRARY})
set(LIBCLANG_INCLUDE_DIRS ${LIBCLANG_INCLUDE_DIR})

include(FindPackageHandleStandardArgs)
# handle the QUIETLY and REQUIRED arguments and set LIBCLANG_FOUND to TRUE if
# all listed variables are TRUE
find_package_handle_standard_args(LibClang DEFAULT_MSG
  LIBCLANG_LIBRARY LIBCLANG_INCLUDE_DIR)

mark_as_advanced(LIBCLANG_INCLUDE_DIR LIBCLANG_LIBRARY)
