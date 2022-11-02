cmake_minimum_required(VERSION 3.10)

set(BUILD_TESTING OFF)
set(BUILD_PYTHON_INTERFACE OFF)
set(SIMDE_EXCLUDE_FROM_ALL 1)

# init proxsuite and cmake-module
find_package(Git QUIET)
if(GIT_FOUND AND EXISTS "${PROJECT_SOURCE_DIR}/.git")
  set(proxqp_dir ${PROJECT_SOURCE_DIR}/external_packages/proxqp/proxqp)
  set(proxqp_cmake_module_dir ${proxqp_dir}/cmake-module)
  set(proxqp ${proxqp_dir}/README.md)
  set(cmake-module ${proxqp_cmake_module_dir}/README.md)
  if(NOT EXISTS ${proxqp})
    message("Initialize proxqp submodule.")
    execute_process(COMMAND ${GIT_EXECUTABLE} submodule update --init ${proxqp_dir}
        WORKING_DIRECTORY ${proxqp_dir}
        COMMAND_ERROR_IS_FATAL ANY)
    endif()
  if(NOT EXISTS ${cmake-module})
    message("Initialize cmake module.")
    execute_process(COMMAND ${GIT_EXECUTABLE} submodule update --init ${proxqp_cmake_module_dir}
        WORKING_DIRECTORY ${proxqp_cmake_module_dir}
        COMMAND_ERROR_IS_FATAL ANY)
  endif()
endif()

unset(CMAKE_MODULE_PATH)

add_subdirectory(proxqp)
