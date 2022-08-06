import sys

try:
    from skbuild import setup
except ImportError:
    print(
        "Please update pip, you need pip 10 or greater,\n"
        " or you need to install the PEP 518 requirements in pyproject.toml yourself",
        file=sys.stderr,
    )
    raise

from setuptools import find_packages

import os 
dir_path = os.path.dirname(os.path.realpath(__file__))

with open(os.path.join(dir_path, "README.md"), "r", encoding="utf-8") as fh:
    long_description = fh.read()

version = "0.0.1"

setup(
    name="alpaqa",
    version=version,
    description="Augmented Lagrangian and PANOC solvers for nonconvex numerical optimization",
    long_description=long_description,
    long_description_content_type="text/markdown",
    author="Pieter P",
    author_email="",
    url="https://github.com/kul-optec/alpaqa",
    project_urls={
        "Documentation": "https://kul-optec.github.io/alpaqa",
        "Source": "https://github.com/kul-optec/alpaqa",
        "Bug Tracker": "https://github.com/kul-optec/alpaqa/issues",
    },
    keywords=["optimization"],
    packages=find_packages(where="src"),
    package_dir={"": "src"},
    cmake_args=["-D", "VERIFY_VERSION=" + version],
    cmake_install_dir="src/alpaqa",
    include_package_data=False,
    install_requires=[
        "numpy",
        "casadi",
        "ninja",
        "cmake",
    ],
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering",
        "License :: OSI Approved :: GNU Lesser General Public License v3 (LGPLv3)",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Operating System :: POSIX :: Linux",
        "Operating System :: Microsoft :: Windows",
    ],
    python_requires=">=3.7",
)