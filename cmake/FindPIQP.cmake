# Find PIQP through its CMake package (piqpConfig.cmake)
#
# Once done this will define
#
#  PIQP_FOUND - system has PIQP
#  piqp - the target to link against

find_package(piqp CONFIG QUIET)

if(piqp_FOUND)
  add_library(piqp INTERFACE IMPORTED)
  set_target_properties(piqp PROPERTIES INTERFACE_LINK_LIBRARIES piqp::piqp)
endif()

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(PIQP DEFAULT_MSG piqp_FOUND)
