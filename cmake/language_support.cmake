# cmake/modules/language_support.cmake
#
# Temporary additional general language support is contained within this
# file.

# This additional function definition is needed to provide a workaround for
# CMake bug 9220.

# On debian testing (cmake 2.6.2), I get return code zero when calling
# cmake the first time, but cmake crashes when running a second time
# as follows:
#
#  -- The Fortran compiler identification is unknown
#  CMake Error at /usr/share/cmake-2.6/Modules/CMakeFortranInformation.cmake:7 (GET_FILENAME_COMPONENT):
#    get_filename_component called with incorrect number of arguments
#  Call Stack (most recent call first):
#    CMakeLists.txt:3 (enable_language)
#
# My workaround is to invoke cmake twice.  If both return codes are zero,
# it is safe to invoke ENABLE_LANGUAGE(Fortran OPTIONAL)

function(workaround_9220 language language_works)
  #message("DEBUG: language = ${language}")
  set(text
    "project(test NONE)
cmake_minimum_required(VERSION 3.10.2)
enable_language(${language} OPTIONAL)
"
    )
  file(REMOVE_RECURSE ${CMAKE_BINARY_DIR}/language_tests/${language})
  file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/language_tests/${language})
  file(WRITE ${CMAKE_BINARY_DIR}/language_tests/${language}/CMakeLists.txt
    ${text})
  execute_process(
    COMMAND ${CMAKE_COMMAND} .
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/language_tests/${language}
    RESULT_VARIABLE return_code
    OUTPUT_VARIABLE variable
    ERROR_QUIET)

  string(FIND ${variable} "The Fortran compiler identification is unknown" find_loc)

  if(return_code EQUAL 0 AND find_loc EQUAL -1)
    # Second run
    execute_process(
      COMMAND ${CMAKE_COMMAND} .
      WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/language_tests/${language}
      RESULT_VARIABLE return_code
      OUTPUT_VARIABLE variable
      ERROR_QUIET)
    string(FIND ${variable} "The Fortran compiler identification is unknown" find_loc)
    if(return_code EQUAL 0 AND find_loc EQUAL -1)
      set(${language_works} ON PARENT_SCOPE)
    else()
      set(${language_works} OFF PARENT_SCOPE)
    endif()
  else()
    set(${language_works} OFF PARENT_SCOPE)
  endif()
endfunction()

# Temporary tests of the above function.
#workaround_9220(CXX CXX_language_works)
#message("CXX_language_works = ${CXX_language_works}")
#workaround_9220(CXXp CXXp_language_works)
#message("CXXp_language_works = ${CXXp_language_works}")
