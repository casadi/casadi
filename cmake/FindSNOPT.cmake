message("SNOPT: find via CONFIG")
find_package(SNOPT CONFIG NO_CMAKE_FIND_ROOT_PATH)
message("SNOPT: ${SNOPT_FOUND}")

if(NOT SNOPT_FOUND)


    find_path(SNOPT_INCLUDE_DIR
        snopt_cwrap.h
        HINTS $ENV{SNOPT}
    )

    if(SNOPT_INCLUDE_DIR)
       set(SNOPT_FOUND_INCLUDE TRUE)
       message(STATUS "Found SNOPT include dir: ${SNOPT_INCLUDE_DIR}")
    else(SNOPT_INCLUDE_DIR)
       message(STATUS "Could not find SNOPT include dir")
    endif(SNOPT_INCLUDE_DIR)

    # Library
    find_library(SNOPT_LIBRARY
    snopt7
    HINTS $ENV{SNOPT} $ENV{SNOPT}/lib
    )

    if(SNOPT_LIBRARY)
      set(SNOPT_LIBRARIES ${SNOPT_LIBRARY})
      message(STATUS "Found SNOPT library: ${SNOPT_LIBRARY}")
    else()
      message(STATUS "Could not find SNOPT libs")
    endif()

    if(SNOPT_LIBRARIES)
      set(SNOPT_FOUND TRUE)
    endif()
    
    if(SNOPT_FOUND)
        add_library(snopt::snopt SHARED IMPORTED)
        set_target_properties(snopt::snopt PROPERTIES
          INTERFACE_INCLUDE_DIRECTORIES ${SNOPT_INCLUDE_DIR}
          IMPORTED_LOCATION "${SNOPT_LIBRARIES}"
        )
    endif()
endif()
