find_path(
  SPRAL_INCLUDE_DIR
  NAMES spral.h
  PATHS $ENV{SPRAL})
find_library(
  SPRAL_LIBRARY
  NAMES spral
  PATHS ${_PREFIX})

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(SPRAL DEFAULT_MSG SPRAL_LIBRARY
                                  SPRAL_INCLUDE_DIR)
mark_as_advanced(SPRAL_INCLUDE_DIR SPRAL_LIBRARY)

if(SPRAL_FOUND AND NOT TARGET spral::spral)
  add_library(spral::spral SHARED IMPORTED)
  set_target_properties(
    spral::spral
    PROPERTIES INTERFACE_INCLUDE_DIRECTORIES "${SPRAL_INCLUDE_DIR}"
               IMPORTED_LOCATION_RELEASE "${SPRAL_LIBRARY}"
               IMPORTED_CONFIGURATIONS "RELEASE")
  message(
    STATUS "SPRAL found (include: ${SPRAL_INCLUDE_DIR}, lib: ${SPRAL_LIBRARY})")
endif()
