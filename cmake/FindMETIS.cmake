message("METIS: find via CONFIG")
find_package(METIS CONFIG NO_CMAKE_FIND_ROOT_PATH)
message("METIS: ${METIS_FOUND}")

if(NOT METIS_FOUND)

    # TRY TO FIND THE METIS LIBRARY

    find_library(METIS_LIB
        names coinmetis metis
        HINTS $ENV{METIS})
        
    find_path(METIS_INCLUDE_DIR
        metis.h
        HINTS $ENV{METIS})

    if(METIS_LIB)
      message(STATUS "Found METIS: ${METIS_LIB}")
      set(METIS_LIBRARIES ${METIS_LIB})
    else()
      message(STATUS "Could not find METIS; looking in environmental variable METIS ($ENV{METIS})")
    endif()

    if(METIS_INCLUDE_DIR)
      message(STATUS "Found METIS include dir: ${METIS_INCLUDE_DIR}")
    else()
      message(STATUS "Could not find METIS include dir; looking in environmental variable METIS ($ENV{METIS})")
    endif()

    if(METIS_INCLUDE_DIR AND METIS_LIB)
      set(METIS_FOUND TRUE)
    else()
      set(METIS_FOUND FALSE)
    endif()



    if(METIS_FOUND)
        add_library(metis::metis INTERFACE IMPORTED)
        set_target_properties(metis::metis PROPERTIES
          INTERFACE_INCLUDE_DIRECTORIES ${METIS_INCLUDE_DIR}
        )
        if("${METIS_LIB}" MATCHES "coinmetis")
            set_target_properties(metis::metis PROPERTIES
              IMPORTED_LIBNAME coinmetis
            )
        else()
            set_target_properties(metis::metis PROPERTIES
              IMPORTED_LIBNAME metis
            )
        endif()
        target_link_libraries(metis::metis INTERFACE ${METIS_LIBRARIES})
    endif()
endif()
