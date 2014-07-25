find_path(MATLAB_INCLUDE_DIR 
extern/include/engine.h
HINTS $ENV{MATLAB}
)

# includes
if(MATLAB_INCLUDE_DIR)
   set(MATLAB_FOUND_INCLUDE TRUE)
   set(MATLAB_INCLUDE_DIRS ${MATLAB_INCLUDE_DIR}/extern/include)
   message(STATUS "Found MATLAB include dir: ${MATLAB_INCLUDE_DIRS}")
else(MATLAB_INCLUDE_DIR)
   message(STATUS "Could not find MATLAB include dir")
endif(MATLAB_INCLUDE_DIR)

find_path(MATLAB_LIBS_DIR
libeng.so
HINTS $ENV{MATLAB}/bin/*
)

# libraries
set(MATLAB_LIBS_LIST
	eng
	mx
)

set(MATLAB_LIBRARIES )
foreach(LIB in ${MATLAB_LIBS_LIST})
  find_library(MATLAB_LIB_${LIB}
    NAMES ${LIB}
    HINTS ${MATLAB_LIBS_DIR})
  if(MATLAB_LIB_${LIB})
    #message(STATUS "Found ${LIB}: ${MATLAB_LIB_${LIB}}")
    set(MATLAB_LIBRARIES ${MATLAB_LIBRARIES} ${MATLAB_LIB_${LIB}})
  else(MATLAB_LIB_${LIB})
    #message(STATUS "Could not find lib ${MATLAB_LIB_${LIB}}")
  endif(MATLAB_LIB_${LIB})
endforeach(LIB)

if(MATLAB_LIBRARIES)
   set(MATLAB_LIBRARIES ${MATLAB_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT} ${LAPACK_LIBRARIES} ${EXTRA_LIBRARIES} ${CMAKE_DL_LIBS})
   message(STATUS "Found Matlab libs: ${MATLAB_LIBRARIES}")
   set(MATLAB_FOUND_LIBS TRUE)
else(MATLAB_LIBRARIES)
   message(STATUS "Could not find Matlab libs")
endif(MATLAB_LIBRARIES)


if(MATLAB_FOUND_INCLUDE AND MATLAB_FOUND_LIBS)
  set(MATLAB_FOUND TRUE)
endif(MATLAB_FOUND_INCLUDE AND MATLAB_FOUND_LIBS)
