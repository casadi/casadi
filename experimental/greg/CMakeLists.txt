if(SNOPT_FOUND)
add_subdirectory(snopt)
endif()

if(MATLAB_FOUND)
add_subdirectory(matlab)
endif()

add_subdirectory(cppocp)
add_subdirectory(examples)

