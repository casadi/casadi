# Target options
# --------------

add_library(alpaqa-obj OBJECT
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

    "include/alpaqa/inner/guarded-aa-pga.hpp"
    "include/alpaqa/inner/second-order-panoc.hpp"
    "include/alpaqa/inner/panoc.hpp"
    "include/alpaqa/inner/lbfgspp.hpp"
    "include/alpaqa/inner/directions/lbfgs.hpp"
    "include/alpaqa/inner/directions/decl/lbfgs.hpp"
    "include/alpaqa/inner/directions/decl/lbfgs-fwd.hpp"
    "include/alpaqa/inner/directions/decl/panoc-direction-update.hpp"
    "include/alpaqa/inner/directions/decl/specialized-lbfgs.hpp"
    "include/alpaqa/inner/directions/anderson-acceleration.hpp"
    "include/alpaqa/inner/directions/specialized-lbfgs.hpp"
    "include/alpaqa/inner/decl/panoc-fwd.hpp"
    "include/alpaqa/inner/decl/second-order-panoc.hpp"
    "include/alpaqa/inner/decl/panoc.hpp"
    "include/alpaqa/inner/decl/structured-panoc-lbfgs.hpp"
    "include/alpaqa/inner/decl/panoc-stop-crit.hpp"
    "include/alpaqa/inner/decl/lbfgs-stepsize.hpp"
    "include/alpaqa/inner/detail/limited-memory-qr.hpp"
    "include/alpaqa/inner/detail/anderson-helpers.hpp"
    "include/alpaqa/inner/detail/panoc-helpers.hpp"
    "include/alpaqa/inner/structured-panoc-lbfgs.hpp"
    "include/alpaqa/inner/pga.hpp"
    "include/alpaqa/inner/newton.hpp"
    "include/alpaqa/alm.hpp"
    "include/alpaqa/reference-problems/himmelblau.hpp"
    "include/alpaqa/reference-problems/riskaverse-mpc.hpp"
    "include/alpaqa/decl/alm.hpp"
    "include/alpaqa/detail/alm-helpers.hpp"
    "include/alpaqa/util/problem.hpp"
    "include/alpaqa/util/solverstatus.hpp"
    "include/alpaqa/util/box.hpp"
    "include/alpaqa/util/atomic_stop_signal.hpp"
    "include/alpaqa/util/vec.hpp"
    "include/alpaqa/util/ringbuffer.hpp"
    "include/alpaqa/util/lipschitz.hpp"
)
add_library(alpaqa::alpaqa-obj ALIAS alpaqa-obj)
set_property(TARGET alpaqa-obj PROPERTY POSITION_INDEPENDENT_CODE 1)
target_include_directories(alpaqa-obj
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:include>
    PRIVATE
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src>
)
target_link_libraries(alpaqa-obj PUBLIC Eigen3::Eigen)
if (LBFGSpp_FOUND)
    target_link_libraries(alpaqa-obj PUBLIC LBFGSpp::LBFGSpp)
endif()

add_library(alpaqa)
target_link_libraries(alpaqa PUBLIC alpaqa-obj)
add_library(alpaqa::alpaqa ALIAS alpaqa)
set_target_properties(alpaqa PROPERTIES EXPORT_NAME alpaqa::alpaqa)

# target_compile_definitions(alpaqa PUBLIC PRINT_DEBUG_COUT=1)

# Interop
# -------

# CUTEst
if (TARGET CUTEst::headers)
    add_library(cutest-loader STATIC 
        "src/interop/cutest/CUTEstLoader.cpp"
        "include/alpaqa/interop/cutest/CUTEstLoader.hpp")
    target_link_libraries(cutest-loader 
        PRIVATE
            CUTEst::headers
            ${CMAKE_DL_LIBS}
        PUBLIC
            alpaqa::alpaqa)
    add_library(alpaqa::cutest-loader ALIAS cutest-loader)
endif()

# CasADi
if (TARGET casadi)
    add_library(casadi-loader-obj OBJECT 
        "src/interop/casadi/CasADiLoader.cpp"
        "include/alpaqa/interop/casadi/CasADiFunctionWrapper.hpp"
        "include/alpaqa/interop/casadi/CasADiLoader.hpp")
    target_link_libraries(casadi-loader-obj 
        PUBLIC
            casadi
            alpaqa::alpaqa-obj)
    set_target_properties(casadi-loader-obj PROPERTIES POSITION_INDEPENDENT_CODE ON)
    add_library(alpaqa::casadi-loader-obj ALIAS casadi-loader-obj)
    add_library(casadi-loader)
    target_link_libraries(casadi-loader PUBLIC casadi-loader-obj)
    add_library(alpaqa::casadi-loader ALIAS casadi-loader)
endif()

# Target installation
# -------------------

if (NOT SKBUILD)
    include(cmake/Install.cmake)
endif()

# Python bindings
# ---------------

option(WITH_PYTHON_BINDINGS "Build the Python module" Off)

if (WITH_PYTHON_BINDINGS OR SKBUILD)
    add_subdirectory(alpaqa)
endif()