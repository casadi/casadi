#### Taken from http://www.openflipper.org/svnrepo/CoMISo/trunk/CoMISo/cmake/FindGUROBI.cmake


# - Try to find GUROBI
# Once done this will define
#  GUROBI_FOUND - System has Gurobi
#  GUROBI_INCLUDE_DIRS - The Gurobi include directories
#  GUROBI_LIBRARIES - The libraries needed to use Gurobi

message("GUROBI: find via CONFIG")
find_package(CPLEX CONFIG NO_CMAKE_FIND_ROOT_PATH)
message("GUROBI: ${GUROBI_FOUND}")

if(NOT GUROBI_FOUND)
    if (GUROBI_INCLUDE_DIR)
      # in cache already
      set(GUROBI_FOUND TRUE)
      set(GUROBI_INCLUDE_DIRS "${GUROBI_INCLUDE_DIR}" )
      set(GUROBI_LIBRARIES "${GUROBI_LIBRARY}" )
    else (GUROBI_INCLUDE_DIR)

    find_path(GUROBI_INCLUDE_DIR 
              NAMES gurobi_c.h
              PATHS "$ENV{GUROBI_HOME}/include"
                      "/Library/gurobi650/mac64/include"
                     "C:\\libs\\gurobi650\\include"
              )

    find_library( GUROBI_LIBRARY 
                  NAMES gurobi
                  gurobi65
                  gurobi70
                  gurobi75
                  gurobi80
                  gurobi81
                  gurobi90
                  PATHS "$ENV{GUROBI_HOME}/lib" 
                        "/Library/gurobi650/mac64/lib"
                        "C:\\libs\\gurobi650\\lib"
                  )


    set(GUROBI_INCLUDE_DIRS "${GUROBI_INCLUDE_DIR}" )
    set(GUROBI_LIBRARIES "${GUROBI_LIBRARY}" )

    # use c++ headers as default
    # set(GUROBI_COMPILER_FLAGS "-DIL_STD" CACHE STRING "Gurobi Compiler Flags")

    include(FindPackageHandleStandardArgs)
    # handle the QUIETLY and REQUIRED arguments and set LIBCPLEX_FOUND to TRUE
    # if all listed variables are TRUE
    find_package_handle_standard_args(GUROBI  DEFAULT_MSG
                                      GUROBI_LIBRARY GUROBI_INCLUDE_DIR)

    mark_as_advanced(GUROBI_INCLUDE_DIR GUROBI_LIBRARY)

    endif(GUROBI_INCLUDE_DIR)

    if(GUROBI_FOUND)
      add_library(gurobi::gurobi INTERFACE IMPORTED)
      set_target_properties(gurobi::gurobi PROPERTIES
        INTERFACE_INCLUDE_DIRECTORIES ${GUROBI_INCLUDE_DIR}
      )
      target_link_libraries(gurobi::gurobi INTERFACE ${GUROBI_LIBRARIES})
    endif()

endif()
