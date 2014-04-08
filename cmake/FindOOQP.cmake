# TRY TO FIND THE INCLUDE DIRECTORY
find_path(OOQP_INCLUDE_DIR 
QpGenData.h
HINTS /usr/local/include/ooqp/
)

set(OOQP_LIBRARIES ${BLAS_LIBRARIES} ${HSL_LIBRARIES})

if(OOQP_INCLUDE_DIR)
  set(OOQP_FOUND_INCLUDE TRUE)
  set(OOQP_INCLUDE_DIRS 
  ${OOQP_INCLUDE_DIR})
  message(STATUS "Found OOQP include dirs: ${OOQP_INCLUDE_DIRS}")
else(OOQP_INCLUDE_DIR)
  message(STATUS "Could not find OOQP include dir")
endif(OOQP_INCLUDE_DIR)

# TRY TO FIND THE LIBRARIES
set(OOQP_LIBS_LIST
  ooqpgensparse ooqpsparse ooqpgondzio ooqpbase
)

set(OOQP_FOUND_LIBS TRUE)
foreach(LIB ${OOQP_LIBS_LIST})
  find_library(OOQP_LIB_${LIB}
    NAMES ${LIB}
    HINTS /usr/local/libs/)
  if(OOQP_LIB_${LIB})
    set(OOQP_LIBRARIES ${OOQP_LIBRARIES} ${OOQP_LIB_${LIB}})
  else(OOQP_LIB_${LIB})
    set(OOQP_FOUND_LIBS FALSE)
  endif(OOQP_LIB_${LIB})
endforeach(LIB)

# print OOQP_LIBRARIES
if(OOQP_FOUND_LIBS)
  message(STATUS "Found OOQP libraries: ${OOQP_LIBRARIES}")  
elseif(OOQP_FOUND_INCLUDE AND OOQP_FOUND_LIBS)
  message(STATUS "Cound not find OOQP libraries")  
endif(OOQP_FOUND_LIBS)

# SUCCESS if BOTH THE LIBRARIES AND THE INCLUDE DIRECTORIES WERE FOUND
if(OOQP_FOUND_INCLUDE AND OOQP_FOUND_LIBS)
  set(OOQP_FOUND TRUE)
  message(STATUS "Found OOQP")  
elseif(OOQP_FOUND_INCLUDE AND OOQP_FOUND_LIBS)
  message(STATUS "Cound not find OOQP")  
endif(OOQP_FOUND_INCLUDE AND OOQP_FOUND_LIBS)


