find_package(SNOPT QUIET)

if(SNOPT_FOUND)
add_subdirectory(snopt)
endif()

find_package(MATLAB QUIET)
if(MATLAB_FOUND)
add_subdirectory(matlab)
endif()

add_subdirectory(cppocp)
add_subdirectory(examples)

