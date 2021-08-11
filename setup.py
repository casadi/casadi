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


setup(
    name="panocpy",
    version="0.0.1",
    description="PANOC+ALM solvers",
    author="Pieter P",
    license="TBD",
    packages=find_packages(where = 'src'),
    package_dir={"": "src"},
    cmake_args=[],
    cmake_install_dir="src/panocpy",
    include_package_data = False,
    install_requires=[
          'numpy',
          'casadi',
    ],
)