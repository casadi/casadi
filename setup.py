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

setup(
    name="panocpy",
    version="0.0.2",
    description="PANOC+ALM solvers",
    long_description=long_description,
    long_description_content_type="text/markdown",
    author="Pieter P",
    url="https://github.com/tttapa/PANOC-ALM",
    keywords=["optimization"],
    packages=find_packages(where="src"),
    package_dir={"": "src"},
    cmake_args=[],
    cmake_install_dir="src/panocpy",
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
    ],
    python_requires=">=3.7",
)