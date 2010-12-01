cmake_minimum_required(VERSION 2.0)

#SET(CMAKE_CXX_FLAGS "-O3 -std=c++0x -DWITH_CPP0X")
# SET(CMAKE_CXX_FLAGS "-O3 ")

# Test working with symbolic variables
add_executable(test_matrix test_matrix.cpp)
target_link_libraries(
  test_matrix
  casadi)

# Test working with symbolic variables
add_executable(test_expression test_expression.cpp)
target_link_libraries(
  test_expression
  casadi)

# Test working with the function trees
add_executable(test_nodes test_nodes.cpp)
target_link_libraries(
  test_nodes
  casadi)

# add_executable(code_generation code_generation.cpp)
# target_link_libraries(code_generation casadi dl)

# Test working with symbolic variables
add_executable(speed_test speed_test.cpp)
target_link_libraries(
  speed_test
  casadi)

# Test working with symbolic variables
add_executable(det_minor det_minor.cpp)
target_link_libraries(
  det_minor
  casadi)

# # Test working with symbolic variables
# add_executable(test_runtime test_runtime.cpp)
# target_link_libraries(
#   test_runtime
#   casadi)


