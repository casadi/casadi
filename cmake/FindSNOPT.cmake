# libraries
set(SNOPT_LIBS_LIST snopt7 snopt7_cpp)

set(SNOPT_LIBRARIES)
set(SNOPT_HAS_EVERYTHING TRUE)
foreach(LIB ${SNOPT_LIBS_LIST})
  find_library(SNOPT_LIB_${LIB}
    NAMES ${LIB}
    HINTS $ENV{SNOPT}/lib)
  if(SNOPT_LIB_${LIB})
    message(STATUS "Found ${LIB}: ${SNOPT_LIB_${LIB}}")
    set(SNOPT_LIBRARIES ${SNOPT_LIBRARIES} ${SNOPT_LIB_${LIB}})
  else()
    set(SNOPT_HAS_EVERYTHING FALSE)
    message(STATUS "Could not find lib${LIB}")
  endif()
endforeach()

if(SNOPT_LIBRARIES)
  set(SNOPT_LIBRARIES ${SNOPT_LIBRARIES})
  message(STATUS "Found Snopt libs")
  set(SNOPT_FOUND_LIBS TRUE)
else()
  message(STATUS "Could not find Snopt libs")
endif()

if(SNOPT_FOUND_LIBS)
  if(SNOPT_HAS_EVERYTHING)
    set(SNOPT_FOUND TRUE)
  endif()
endif()
