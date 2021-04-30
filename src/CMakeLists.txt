# Target options
# --------------

add_library(panoc-alm SHARED
    "src/inner/panoc.cpp"
    "src/inner/second-order-panoc.cpp"
    "src/inner/second-order-panoc-lbfgs.cpp"
    "src/util/box.cpp"
    "src/inner/directions/lbfgs.cpp"
    "src/alm.cpp"
    "src/util/problem.cpp"
    "src/util/solverstatus.cpp"
    "src/reference-problems/riskaverse-mpc.cpp"
    "src/reference-problems/himmelblau.cpp"
)
add_library(panoc-alm::panoc-alm ALIAS panoc-alm)
set_target_properties(panoc-alm PROPERTIES EXPORT_NAME panoc-alm::panoc-alm)
target_include_directories(panoc-alm
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:include>
    PRIVATE
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src>
)
target_link_libraries(panoc-alm PUBLIC Eigen3::Eigen)
if (LBFGSpp_FOUND)
    target_link_libraries(panoc-alm PUBLIC LBFGSpp::LBFGSpp)
endif()

# target_compile_definitions(panoc-alm PUBLIC PRINT_DEBUG_COUT=1)

# Interop
# -------

# CUTEst
if (TARGET CUTEst::headers)
    add_library(cutest-loader STATIC src/interop/cutest/CUTEstLoader.cpp)
    target_link_libraries(cutest-loader 
        PRIVATE
            CUTEst::headers
            ${CMAKE_DL_LIBS}
        PUBLIC
            panoc-alm::panoc-alm)
    add_library(panoc-alm::cutest-loader ALIAS cutest-loader)
endif()

# CasADi
if (TARGET casadi)
    add_library(casadi-loader STATIC src/interop/casadi/CasADiLoader.cpp)
    target_link_libraries(casadi-loader 
        PUBLIC
            casadi
            panoc-alm::panoc-alm)
    add_library(panoc-alm::casadi-loader ALIAS casadi-loader)
    set_target_properties(casadi-loader PROPERTIES POSITION_INDEPENDENT_CODE ON)
endif()

# Target installation
# -------------------

if (NOT SKBUILD)
    include(cmake/Install.cmake)
endif()

# Python bindings
# ---------------

add_subdirectory(panocpy)