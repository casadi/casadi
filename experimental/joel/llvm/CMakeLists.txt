# Tutorial example
add_executable(llvm_tutorial llvm_tutorial.cpp)
target_link_libraries(llvm_tutorial ${CASADI_DEPENDENCIES})

# Simple example
add_executable(test_llvm test_llvm.cpp)
target_link_libraries(test_llvm ${CASADI_DEPENDENCIES})

# Test JIT in CasADi
add_executable(test_llvm_casadi test_llvm_casadi.cpp)
target_link_libraries(test_llvm_casadi casadi ${CASADI_DEPENDENCIES})

