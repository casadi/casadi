cmake_minimum_required (VERSION 2.6)
project(tutor)

find_package(PYREPORT)

if(FOUND_PYREPORT)
  message(STATUS "pyreport found")
  add_subdirectory(python)
else()
  message(STATUS "pyreport not found -- try 'easy_install pyreport' ")
endif()