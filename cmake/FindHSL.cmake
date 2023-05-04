

if(NOT HSL_FOUND)

    message("HSL: find via CONFIG")
    find_package(HSL CONFIG NO_CMAKE_FIND_ROOT_PATH)
    message("HSL: ${HSL_FOUND}")
    
    if(NOT HSL_FOUND)
        # TRY TO FIND THE HSL LIBRARY

        find_library(HSL_LIB
            names hsl coinhsl
            HINTS $ENV{HSL})

        if(HSL_LIB)
          message(STATUS "Found HSL: ${HSL_LIB}")
          set(HSL_LIBRARIES ${HSL_LIB})
        else(HSL_LIB)
          message(STATUS "Could not find HSL; looking in environmental variable HSL ($ENV{HSL})")
        endif(HSL_LIB)

        # Set standard flags
        include(FindPackageHandleStandardArgs)
        find_package_handle_standard_args(HSL HSL_LIB)

        if(HSL_FOUND)
          message("debug HSL_FOUND")
          add_library(hsl::hsl INTERFACE IMPORTED)
          target_link_libraries(hsl::hsl INTERFACE ${HSL_LIB})
        endif()
        
    endif()
endif()
