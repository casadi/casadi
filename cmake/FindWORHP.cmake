message("WORHP: find via CONFIG")
find_package(WORHP CONFIG NO_CMAKE_FIND_ROOT_PATH)
message("WORHP: ${WORHP_FOUND}")

if(NOT WORHP_FOUND)

    # TRY TO FIND THE INCLUDE DIRECTORY
    find_file(WORHP_LICENSE_KEY
    worhp.lic
    HINTS $ENV{WORHP} $ENV{WORHP_LICENSE_FILE}
    )

    if(WORHP_LICENSE_KEY)
      message(STATUS "Worhp license detected: ${WORHP_LICENSE_KEY}")
    else()
      message(STATUS "Could not find worhp license key, looking in $ENV{WORHP}")
    endif()

    find_path(WORHP_INCLUDE_DIR
    worhp.h
    HINTS $ENV{WORHP}/include/worhp/ $ENV{WORHP})

    if(WORHP_INCLUDE_DIR)
      message(STATUS "Worhp include dir detected: ${WORHP_INCLUDE_DIR}")
    else()
      message(STATUS "Could not find worhp include dir, looking in $ENV{WORHP}/include/worhp/")
    endif()

    find_library(WORHP_LIBRARY
      worhp
      HINTS $ENV{WORHP}/lib/ $ENV{WORHP})


    if(WORHP_LIBRARY)
      message(STATUS "Worhp library detected: ${WORHP_LIBRARY}")
    else()
      message(STATUS "Could not find worhp library, looking in $ENV{WORHP}/lib/")
    endif()

    set(WORHP_LIBRARIES ${WORHP_LIBRARY})


    if(WORHP_LIBRARY AND WORHP_INCLUDE_DIR)
      message(STATUS "Worhp interface ready with these libraries: ${WORHP_LIBRARIES} ")
      set(WORHP_FOUND TRUE)
    else()
      message(STATUS "Will not compile worhp interface")
    endif()
    
    add_library(worhp::worhp INTERFACE IMPORTED)
    set_target_properties(worhp::worhp PROPERTIES
      INTERFACE_INCLUDE_DIRECTORIES ${WORHP_INCLUDE_DIR}
    )
    target_link_libraries(worhp::worhp INTERFACE ${WORHP_LIBRARIES})
    
    if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang" OR "${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
        target_compile_options(worhp::worhp INTERFACE "-Wno-unused-function")
    endif()

endif()
