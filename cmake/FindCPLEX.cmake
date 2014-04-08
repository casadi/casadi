set(CPLEX_HEADER_LIST
   ilcplex/ilocplex.h
   ilconcert/ilomodel.h
)

set(CPLEX_LIBRARY_LIST
   ilocplex
   cplex
   concert
)

set(CPLEX_FOUND_HEADERS TRUE)
foreach(HEADER ${CPLEX_HEADER_LIST})
    find_path(CPLEX_PATH_${HEADER}
        NAMES ${HEADER}
    )
    if(CPLEX_PATH_${HEADER})
        set(CPLEX_INCLUDE_DIR ${CPLEX_INCLUDE_DIR} ${CPLEX_PATH_${HEADER}})
    else(CPLEX_PATH_${HEADER})
        set(CPLEX_FOUND_HEADERS FALSE)
    endif(CPLEX_PATH_${HEADER})
endforeach(HEADER)

if(CPLEX_FOUND_HEADERS)
    message(STATUS "Found CPLEX include dir: ${CPLEX_INCLUDE_DIR}")
else(CPLEX_FOUND_HEADERS)
    message(STATUS "Could not find CPLEX include dir")
endif(CPLEX_FOUND_HEADERS)

set(CPLEX_FOUND_LIBRARIES TRUE)
foreach(LIBRARY ${CPLEX_LIBRARY_LIST})
    find_library(CPLEX_LIB_${LIBRARY}
        NAMES ${LIBRARY}
    )
    if(CPLEX_LIB_${LIBRARY})
        set(CPLEX_LIBRARIES ${CPLEX_LIBRARIES} ${CPLEX_LIB_${LIBRARY}})
    else(CPLEX_LIB_${LIBRARY})
        set(CPLEX_FOUND_LIBRARIES FALSE)
    endif(CPLEX_LIB_${LIBRARY})
endforeach(LIBRARY)

if(CPLEX_FOUND_LIBRARIES)
    message(STATUS "Found CPLEX libraries: ${CPLEX_LIBRARIES}")
else(CPLEX_FOUND_LIBRARIES)
    message(STATUS "Could not find CPLEX libraries")
endif(CPLEX_FOUND_LIBRARIES)


if(CPLEX_FOUND_HEADERS AND CPLEX_FOUND_LIBRARIES)
  set(CPLEX_FOUND TRUE)
endif(CPLEX_FOUND_HEADERS AND CPLEX_FOUND_LIBRARIES)
