import os

from conan import ConanFile
from conan.errors import ConanInvalidConfiguration
from conan.tools.cmake import CMakeDeps, CMakeToolchain, CMake, cmake_layout
from conan.tools.build import can_run


class AlpaqaRecipe(ConanFile):
    name = "alpaqa"
    version = "1.1.0-alpha.1"

    # Optional metadata
    license = "LGPL-3.0-or-later"
    author = "Pieter P <pieter.p.dev@outlook.com>"
    url = "https://github.com/kul-optec/alpaqa"
    description = (
        "Augmented Lagrangian and PANOC solvers for nonconvex numerical optimization"
    )
    topics = ("optimization", "panoc", "alm", "mpc")

    # Binary configuration
    settings = "os", "compiler", "build_type", "arch"
    bool_alpaqa_options = {
        "with_python": False,
        "with_matlab": False,
        "with_drivers": True,
        "with_examples": False,
        "with_python_problem_loader": False,
        "with_gradient_checker": False,
        "with_casadi": True,
        "with_external_casadi": False,
        "with_cutest": False,
        "with_qpalm": False,
        "with_json": True,
        "with_lbfgsb": False,
        "with_ipopt": False,
        "with_ocp": False,
        "with_casadi_ocp": False,
        "with_openmp": False,
        "with_quad_precision": False,
        "with_single_precision": False,
        "with_long_double": False,
        "debug_checks_eigen": False,
        "dont_parallelize_eigen": True,
        "no_dlclose": False,
        "with_blas": False,
        "with_coverage": False,
    }
    options = (
        {
            "shared": [True, False],
            "fPIC": [True, False],
        }
        | {k: [True, False] for k in bool_alpaqa_options}
        | {"with_conan_python": [True, False]}
    )
    default_options = (
        {
            "shared": False,
            "fPIC": True,
        }
        | bool_alpaqa_options
        | {"with_conan_python": False}
    )

    # Sources are located in the same place as this recipe, copy them to the recipe
    exports_sources = (
        "CMakeLists.txt",
        "src/*",
        "cmake/*",
        "examples/*",
        "interfaces/*",
        "python/*",
        "test/*",
        "LICENSE",
        "README.md",
    )

    def requirements(self):
        self.requires("eigen/tttapa.20250504", transitive_headers=True, force=True)
        self.requires("guanaqo/1.0.0-alpha.16", transitive_headers=True)
        self.test_requires("gtest/1.17.0")
        if self.options.with_external_casadi:
            self.requires("casadi/3.7.1", transitive_headers=True)
        if self.options.with_json:
            self.requires("nlohmann_json/3.12.0", transitive_headers=True)
        if self.options.with_ipopt:
            self.requires("ipopt/3.14.16", transitive_headers=True)
        if self.options.with_qpalm:
            self.requires("qpalm/1.2.6", transitive_headers=True)
        if self.options.with_python or self.options.with_python_problem_loader:
            self.requires("pybind11/2.13.6")
            if self.options.with_conan_python:
                self.requires("tttapa-python-dev/3.13.7")
        if self.options.with_matlab:
            self.requires("utfcpp/4.0.5")
        if self.options.with_blas:
            self.requires("openblas/0.3.30")

    def build_requirements(self):
        self.tool_requires("cmake/[>=4.1 <5]")

    def config_options(self):
        if self.settings.get_safe("os") == "Windows":
            self.options.rm_safe("fPIC")

    def validate(self):
        if self.options.with_matlab and not self.options.with_json:
            msg = "MATLAB MEX interface requires JSON. Set 'with_json=True'."
            raise ConanInvalidConfiguration(msg)
        if self.options.with_matlab and not self.options.with_external_casadi:
            msg = (
                "MATLAB MEX interface requires CasADi. Set 'with_external_casadi=True'."
            )
            raise ConanInvalidConfiguration(msg)

    def configure(self):
        if self.options.get_safe("with_quad_precision"):
            self.options["guanaqo/*"].with_quad_precision = True

    def layout(self):
        if self.folders.build_folder_vars is None:
            if self.options.with_python:
                self.folders.build_folder_vars = ["const.python"]
            if self.options.with_matlab:
                self.folders.build_folder_vars = ["const.matlab"]
        cmake_layout(self)

    def generate(self):
        deps = CMakeDeps(self)
        deps.set_property("utfcpp", "cmake_target_name", "utf8cpp::utf8cpp")
        deps.generate()
        tc = CMakeToolchain(self)
        tc.user_presets_path = "ConanPresets.json"
        for k in self.bool_alpaqa_options:
            value = getattr(self.options, k, None)
            if value is not None and value.value is not None:
                tc.variables["ALPAQA_" + k.upper()] = bool(value)
        if can_run(self):
            tc.variables["ALPAQA_FORCE_TEST_DISCOVERY"] = True
        tc.generate()

    def build(self):
        cmake = CMake(self)
        cmake.configure()
        cmake.build()
        cmake.test()

    def package(self):
        cmake = CMake(self)
        cmake.install()
        if self.options.with_python:
            cmake.install(component="python_source")
            cmake.install(component="python_modules")
            cmake.install(component="python_stubs")

    def package_info(self):
        self.cpp_info.set_property("cmake_find_mode", "none")
        self.cpp_info.builddirs.append(os.path.join("lib", "cmake", "alpaqa"))
        self.runenv_info.prepend_path("PATH", os.path.join(self.package_folder, "bin"))
