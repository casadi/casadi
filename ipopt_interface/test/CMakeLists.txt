find_package(IPOPT REQUIRED QUIET)

include_directories(../../)
include_directories(${IPOPT_INCLUDE_DIR})

find_package(MPI QUIET)

# My version of hs071
add_executable(my_hs071 my_hs071.cpp)
target_link_libraries(my_hs071 ipopt_interface casadi ${IPOPT_LIBRARIES})

# 16-th root of 2
add_executable(root_test root_test.cpp)
target_link_libraries(root_test ipopt_interface casadi ${IPOPT_LIBRARIES})

# Oscillator example
add_executable(oscillator oscillator.cpp)
target_link_libraries(oscillator ipopt_interface casadi ${IPOPT_LIBRARIES})

# Philippe Dreesen's structured least squares problem
add_executable(structured_least_squares structured_least_squares.cpp)
target_link_libraries(structured_least_squares ipopt_interface casadi ${IPOPT_LIBRARIES})

# Ipopt example
add_executable(test_ipopt test_ipopt.cpp)
target_link_libraries(test_ipopt ${IPOPT_LIBRARIES})

# Tutorial example hs071
add_executable(hs071_main hs071_main.cpp hs071_nlp.cpp hs071_nlp.hpp)
target_link_libraries(hs071_main ${IPOPT_LIBRARIES})

