message("KNITRO: find via CONFIG")
find_package(KNITRO CONFIG NO_CMAKE_FIND_ROOT_PATH)
message("KNITRO: ${KNITRO_FOUND}")

if(NOT KNITRO_FOUND)

    find_path(KNITRO_INCLUDE_DIR
    knitro.h
    HINTS $ENV{KNITRO}/include
    )

    if(KNITRO_INCLUDE_DIR)
       set(KNITRO_FOUND_INCLUDE TRUE)
       message(STATUS "Found KNITRO include dir: ${KNITRO_INCLUDE_DIR}")
    else(KNITRO_INCLUDE_DIR)
       message(STATUS "Could not find KNITRO include dir")
    endif(KNITRO_INCLUDE_DIR)

    find_library(KNITRO_LIBRARY NAMES knitro HINTS $ENV{KNITRO}/lib)

    if(KNITRO_LIBRARY)
      set(KNITRO_LIBRARIES ${KNITRO_LIBRARY} $ENV{KNITRO_EXTRA_LIBS})
       message(STATUS "Found KNITRO libraries: ${KNITRO_LIBRARIES}")
    else(KNITRO_LIBRARY)
       message(STATUS "Could not find KNITRO library")
    endif(KNITRO_LIBRARY)

    if(KNITRO_FOUND_INCLUDE AND KNITRO_LIBRARIES)
      set(KNITRO_FOUND TRUE)
    endif(KNITRO_FOUND_INCLUDE AND KNITRO_LIBRARIES)


    add_library(knitro::knitro SHARED IMPORTED)
    set_target_properties(knitro::knitro PROPERTIES
      INTERFACE_INCLUDE_DIRECTORIES ${KNITRO_INCLUDE_DIR}
      IMPORTED_LOCATION "${KNITRO_LIBRARIES}"
    )
    
endif()
