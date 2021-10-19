# Target options
# --------------

add_library(panoc-alm-obj OBJECT
    "src/inner/panoc.cpp"
    "src/inner/second-order-panoc.cpp"
    "src/inner/structured-panoc-lbfgs.cpp"
    "src/util/box.cpp"
    "src/inner/directions/lbfgs.cpp"
    "src/alm.cpp"
    "src/util/problem.cpp"
    "src/util/solverstatus.cpp"
    "src/reference-problems/riskaverse-mpc.cpp"
    "src/reference-problems/himmelblau.cpp"

    "include/panoc-alm/inner/guarded-aa-pga.hpp"
    "include/panoc-alm/inner/second-order-panoc.hpp"
    "include/panoc-alm/inner/panoc.hpp"
    "include/panoc-alm/inner/lbfgspp.hpp"
    "include/panoc-alm/inner/directions/lbfgs.hpp"
    "include/panoc-alm/inner/directions/decl/lbfgs.hpp"
    "include/panoc-alm/inner/directions/decl/lbfgs-fwd.hpp"
    "include/panoc-alm/inner/directions/decl/panoc-direction-update.hpp"
    "include/panoc-alm/inner/directions/decl/specialized-lbfgs.hpp"
    "include/panoc-alm/inner/directions/anderson-acceleration.hpp"
    "include/panoc-alm/inner/directions/specialized-lbfgs.hpp"
    "include/panoc-alm/inner/decl/panoc-fwd.hpp"
    "include/panoc-alm/inner/decl/second-order-panoc.hpp"
    "include/panoc-alm/inner/decl/panoc.hpp"
    "include/panoc-alm/inner/decl/structured-panoc-lbfgs.hpp"
    "include/panoc-alm/inner/decl/panoc-stop-crit.hpp"
    "include/panoc-alm/inner/decl/lbfgs-stepsize.hpp"
    "include/panoc-alm/inner/detail/limited-memory-qr.hpp"
    "include/panoc-alm/inner/detail/anderson-helpers.hpp"
    "include/panoc-alm/inner/detail/panoc-helpers.hpp"
    "include/panoc-alm/inner/structured-panoc-lbfgs.hpp"
    "include/panoc-alm/inner/pga.hpp"
    "include/panoc-alm/inner/newton.hpp"
    "include/panoc-alm/alm.hpp"
    "include/panoc-alm/reference-problems/himmelblau.hpp"
    "include/panoc-alm/reference-problems/riskaverse-mpc.hpp"
    "include/panoc-alm/decl/alm.hpp"
    "include/panoc-alm/detail/alm-helpers.hpp"
    "include/panoc-alm/util/problem.hpp"
    "include/panoc-alm/util/solverstatus.hpp"
    "include/panoc-alm/util/box.hpp"
    "include/panoc-alm/util/atomic_stop_signal.hpp"
    "include/panoc-alm/util/vec.hpp"
    "include/panoc-alm/util/ringbuffer.hpp"
    "include/panoc-alm/util/lipschitz.hpp"
)
add_library(panoc-alm::panoc-alm-obj ALIAS panoc-alm-obj)
set_property(TARGET panoc-alm-obj PROPERTY POSITION_INDEPENDENT_CODE 1)
target_include_directories(panoc-alm-obj
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:include>
    PRIVATE
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src>
)
target_link_libraries(panoc-alm-obj PUBLIC Eigen3::Eigen)
if (LBFGSpp_FOUND)
    target_link_libraries(panoc-alm-obj PUBLIC LBFGSpp::LBFGSpp)
endif()

add_library(panoc-alm)
target_link_libraries(panoc-alm PUBLIC panoc-alm-obj)
add_library(panoc-alm::panoc-alm ALIAS panoc-alm)
set_target_properties(panoc-alm PROPERTIES EXPORT_NAME panoc-alm::panoc-alm)

# target_compile_definitions(panoc-alm PUBLIC PRINT_DEBUG_COUT=1)

# Interop
# -------

# CUTEst
if (TARGET CUTEst::headers)
    add_library(cutest-loader STATIC 
        "src/interop/cutest/CUTEstLoader.cpp"
        "include/panoc-alm/interop/cutest/CUTEstLoader.hpp")
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
    add_library(casadi-loader-obj OBJECT 
        "src/interop/casadi/CasADiLoader.cpp"
        "include/panoc-alm/interop/casadi/CasADiFunctionWrapper.hpp"
        "include/panoc-alm/interop/casadi/CasADiLoader.hpp")
    target_link_libraries(casadi-loader-obj 
        PUBLIC
            casadi
            panoc-alm::panoc-alm-obj)
    set_target_properties(casadi-loader-obj PROPERTIES POSITION_INDEPENDENT_CODE ON)
    add_library(panoc-alm::casadi-loader-obj ALIAS casadi-loader-obj)
    add_library(casadi-loader)
    target_link_libraries(casadi-loader PUBLIC casadi-loader-obj)
    add_library(panoc-alm::casadi-loader ALIAS casadi-loader)
endif()

# Target installation
# -------------------

if (NOT SKBUILD)
    include(cmake/Install.cmake)
endif()

# Python bindings
# ---------------

option(WITH_PANOCPY "Build the Python module" Off)

if (WITH_PANOCPY OR SKBUILD)
    add_subdirectory(panocpy)
endif()