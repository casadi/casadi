import os

from conan import ConanFile
from conan.tools.cmake import CMakeToolchain, CMake, cmake_layout


class AlpaqaRecipe(ConanFile):
    name = "alpaqa"
    version = "1.0.0-alpha.8"

    # Optional metadata
    license = "LGPLv3"
    author = "Pieter P <pieter.p.dev@outlook.com>"
    url = "https://github.com/kul-optec/alpaqa"
    description = (
        "Augmented Lagrangian and PANOC solvers for nonconvex numerical optimization"
    )
    topics = ("optimization", "panoc", "alm", "mpc")

    # Binary configuration
    settings = "os", "compiler", "build_type", "arch"
    options = {
        "shared": [True, False],
        "fPIC": [True, False],
        "with_quad_precision": [True, False],
        "with_casadi": [True, False],
        "with_ocp": [True, False],
        "with_drivers": [True, False],
    }
    default_options = {
        "shared": False,
        "fPIC": True,
        "with_quad_precision": False,
        "with_casadi": True,
        "with_ocp": True,
        "with_drivers": True,
    }

    # Sources are located in the same place as this recipe, copy them to the recipe
    exports_sources = (
        "CMakeLists.txt",
        "src/*",
        "cmake/*",
        "interfaces/*",
        "python/*",
        "test/*",
        "LICENSE",
        "README.md",
    )

    generators = ("CMakeDeps", )

    def requirements(self):
        self.requires("eigen/3.4.0")
        self.test_requires("gtest/1.11.0")
        if self.options.with_casadi:
            self.requires("casadi/3.6.0@alpaqa")

    def config_options(self):
        if self.settings.os == "Windows":
            del self.options.fPIC

    def layout(self):
        cmake_layout(self)

    def generate(self):
        tc = CMakeToolchain(self)
        tc.variables["ALPAQA_WITH_QUAD_PRECISION"] = self.options.with_quad_precision
        tc.variables["ALPAQA_WITH_CASADI"] = self.options.with_casadi
        tc.variables["ALPAQA_WITH_EXAMPLES"] = False
        tc.variables["ALPAQA_WITH_DRIVERS"] = self.options.with_drivers
        tc.variables["ALPAQA_WITH_LBFGSB"] = True
        tc.generate()

    def build(self):
        cmake = CMake(self)
        cmake.configure()
        cmake.build()
        cmake.test()

    def package(self):
        cmake = CMake(self)
        cmake.install()

    def package_info(self):
        self.cpp_info.set_property("cmake_find_mode", "none")
        self.cpp_info.builddirs.append(os.path.join("lib", "cmake", "alpaqa"))
