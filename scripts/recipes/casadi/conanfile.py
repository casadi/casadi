import os

from conan import ConanFile
from conan.tools.cmake import CMakeToolchain, CMake, cmake_layout
from conan.tools.files import get


class CasADiRecipe(ConanFile):
    name = "casadi"
    user = "alpaqa"

    # Optional metadata
    license = "LGPL-3.0-or-later"
    author = "Pieter P <pieter.p.dev@outlook.com>"
    url = "https://github.com/kul-optec/alpaqa"
    description = (
        "CasADi is a symbolic framework for numeric optimization implementing "
        "automatic differentiation in forward and reverse modes on sparse matrix-"
        "valued computational graphs. It supports self-contained C-code generation "
        "and interfaces state-of-the-art codes such as SUNDIALS, IPOPT etc. It can "
        "be used from C++, Python or Matlab/Octave."
    )

    # Binary configuration
    settings = "os", "compiler", "build_type", "arch"
    casadi_cmake_options = {
        "with_python": [True, False],
        "with_matlab": [True, False],
        "with_octave": [True, False],
        "with_json": [True, False],
        "with_octave_import": [True, False],
        "with_python3": [True, False],
        "install_internal_headers": [True, False],
        "enable_export_all": [True, False],
        "swig_export": [True, False],
        "swig_import": [True, False],
        "with_openmp": [True, False],
        "with_thread": [True, False],
        "with_thread_mingw": [True, False],
        "with_opencl": [True, False],
        "with_deepbind": [True, False],
        "with_refcount_warnings": [True, False],
        "with_so_version": [True, False],
        "with_copysign_undef": [True, False],
        "with_selfcontained": [True, False],
        "with_extra_checks": [True, False],
        "with_dl": [True, False],
        "with_fmu": [True, False],
        "with_deprecated_features": [True, False],
        "with_build_required": [True, False],
        "with_build_sundials": [True, False],
        "with_sundials": [True, False],
        "with_build_csparse": [True, False],
        "with_csparse": [True, False],
        "with_blasfeo": [True, False],
        "with_build_blasfeo": [True, False],
        "with_hpipm": [True, False],
        "with_build_hpipm": [True, False],
        "with_superscs": [True, False],
        "with_build_superscs": [True, False],
        "with_osqp": [True, False],
        "with_build_osqp": [True, False],
        "with_build_eigen3": [True, False],
        "with_build_simde": [True, False],
        "with_proxqp": [True, False],
        "with_build_proxqp": [True, False],
        "with_build_tinyxml": [True, False],
        "with_tinyxml": [True, False],
        "with_build_dsdp": [True, False],
        "with_dsdp": [True, False],
        "old_llvm": [True, False],
        "with_clang": [True, False],
        "with_lapack": [True, False],
        "with_build_lapack": [True, False],
        "with_qpoases": [True, False],
        "with_no_qpoases_banner": [True, False],
        "with_blocksqp": [True, False],
        "with_ipopt": [True, False],
        "with_build_ipopt": [True, False],
        "with_mockup_required": [True, False],
        "with_knitro": [True, False],
        "with_mockup_knitro": [True, False],
        "with_snopt": [True, False],
        "with_mockup_snopt": [True, False],
        "with_worhp": [True, False],
        "with_mockup_worhp": [True, False],
        "with_cplex": [True, False],
        "with_cplex_shared": [True, False],
        "with_mockup_cplex": [True, False],
        "with_gurobi": [True, False],
        "with_mockup_gurobi": [True, False],
        "with_bonmin": [True, False],
        "with_build_bonmin": [True, False],
        "with_cbc": [True, False],
        "with_build_cbc": [True, False],
        "with_clp": [True, False],
        "with_build_clp": [True, False],
        "with_mumps": [True, False],
        "with_build_metis": [True, False],
        "with_build_mumps": [True, False],
        "with_spral": [True, False],
        "with_build_spral": [True, False],
        "with_hsl": [True, False],
        "with_mockup_hsl": [True, False],
        "with_build_hsl": [True, False],
        "with_highs": [True, False],
        "with_build_highs": [True, False],
        "allow_docker": [True, False],
        "with_ooqp": [True, False],
        "with_sqic": [True, False],
        "with_ampl": [True, False],
        "with_slicot": [True, False],
        "with_examples": [True, False],
        "with_doc": [True, False],
    }
    options = {
        "shared": [True, False],
        "fPIC": [True, False],
    } | casadi_cmake_options
    default_options = {
        "shared": False,
        "fPIC": True,
        "with_python": False,
        "with_matlab": False,
        "with_octave": False,
        "with_json": False,
        "with_octave_import": False,
        "with_python3": False,
        "install_internal_headers": False,
        "enable_export_all": False,
        "swig_export": False,
        "swig_import": False,
        "with_openmp": False,
        "with_thread": True,
        "with_thread_mingw": False,
        "with_opencl": False,
        "with_deepbind": True,
        "with_refcount_warnings": False,
        "with_so_version": True,
        "with_copysign_undef": False,
        "with_selfcontained": False,
        "with_extra_checks": False,
        "with_dl": True,
        "with_fmu": True,
        "with_deprecated_features": True,
        "with_build_required": False,
        "with_build_sundials": False,
        "with_sundials": False,
        "with_build_csparse": False,
        "with_csparse": False,
        "with_blasfeo": False,
        "with_build_blasfeo": False,
        "with_hpipm": False,
        "with_build_hpipm": False,
        "with_superscs": False,
        "with_build_superscs": False,
        "with_osqp": False,
        "with_build_osqp": False,
        "with_build_eigen3": False,
        "with_build_simde": False,
        "with_proxqp": False,
        "with_build_proxqp": False,
        "with_build_tinyxml": True,
        "with_tinyxml": True,
        "with_build_dsdp": False,
        "with_dsdp": False,
        "old_llvm": False,
        "with_clang": False,
        "with_lapack": False,
        "with_build_lapack": False,
        "with_qpoases": False,
        "with_no_qpoases_banner": False,
        "with_blocksqp": False,
        "with_ipopt": False,
        "with_build_ipopt": False,
        "with_mockup_required": False,
        "with_knitro": False,
        "with_mockup_knitro": False,
        "with_snopt": False,
        "with_mockup_snopt": False,
        "with_worhp": False,
        "with_mockup_worhp": False,
        "with_cplex": False,
        "with_cplex_shared": False,
        "with_mockup_cplex": False,
        "with_gurobi": False,
        "with_mockup_gurobi": False,
        "with_bonmin": False,
        "with_build_bonmin": False,
        "with_cbc": False,
        "with_build_cbc": False,
        "with_clp": False,
        "with_build_clp": False,
        "with_mumps": False,
        "with_build_metis": False,
        "with_build_mumps": False,
        "with_spral": False,
        "with_build_spral": False,
        "with_hsl": False,
        "with_mockup_hsl": False,
        "with_build_hsl": False,
        "with_highs": False,
        "with_build_highs": False,
        "allow_docker": False,
        "with_ooqp": False,
        "with_sqic": False,
        "with_ampl": False,
        "with_slicot": False,
        "with_examples": False,
        "with_doc": False,
    }

    def set_version(self):
        latest_version = next(reversed(self.conan_data["sources"].keys()))
        self.version = self.version or latest_version

    def source(self):
        get(self, **self.conan_data["sources"][self.version])

    def requirements(self):
        pass

    def config_options(self):
        if self.settings.os == "Windows":
            del self.options.fPIC

    def layout(self):
        cmake_layout(self)

    def generate(self):
        tc = CMakeToolchain(self)
        tc.variables["ENABLE_SHARED"] = self.options.shared
        tc.variables["ENABLE_STATIC"] = not self.options.shared
        for opt in self.casadi_cmake_options.keys():
            tc.variables[opt.upper()] = getattr(self.options, opt)
        tc.generate()

    def build(self):
        cmake = CMake(self)
        cmake.configure()
        cmake.build()

    def package(self):
        cmake = CMake(self)
        cmake.install()

    def package_info(self):
        self.cpp_info.set_property("cmake_find_mode", "none")
        self.cpp_info.builddirs.append(os.path.join("lib", "cmake", "casadi"))
