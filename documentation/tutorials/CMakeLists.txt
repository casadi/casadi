cmake_minimum_required (VERSION 2.0)
project(tutor)

find_package(PYREPORT)

if(FOUND_PYREPORT)
MESSAGE(STATUS "pyreport found")
add_subdirectory(python)
else()
MESSAGE(STATUS "pyreport not found -- try 'easy_install pyreport' ")
endif()