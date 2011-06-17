find_package(SNOPT QUIET)

IF (SNOPT_FOUND)
add_subdirectory(snopt)
ENDIF (SNOPT_FOUND)

find_package(MATLAB QUIET)
IF (MATLAB_FOUND)
add_subdirectory(matlab)
ENDIF (MATLAB_FOUND)

add_subdirectory(cppocp)
add_subdirectory(examples)

