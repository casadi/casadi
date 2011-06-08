find_package(SNOPT QUIET)

IF (SNOPT_FOUND)
add_subdirectory(snopt)
ENDIF (SNOPT_FOUND)

add_subdirectory(cppocp)
add_subdirectory(examples)
