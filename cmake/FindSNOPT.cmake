# libraries
set(SNOPT_LIBS_LIST
  snopt7
  snprint7
  snblas)

set(SNOPT_LIBRARIES)
foreach(LIB in ${SNOPT_LIBS_LIST})
  find_library(SNOPT_LIB_${LIB}
    NAMES ${LIB}
    HINTS $ENV{SNOPT}/lib)
  if(SNOPT_LIB_${LIB})
    message(STATUS "Found ${LIB}: ${SNOPT_LIB_${LIB}}")
    set(SNOPT_LIBRARIES ${SNOPT_LIBRARIES} ${SNOPT_LIB_${LIB}})
  else()
    message(STATUS "Could not find lib${LIB}")
  endif()
endforeach()

if(SNOPT_LIBRARIES)
  set(SNOPT_LIBRARIES ${SNOPT_LIBRARIES})
  message(STATUS "Found Snopt libs: ${SNOPT_LIBRARIES}")
  set(SNOPT_FOUND_LIBS TRUE)
else()
  message(STATUS "Could not find Snopt libs")
endif()


if(SNOPT_FOUND_LIBS)
  set(SNOPT_FOUND TRUE)
endif()
