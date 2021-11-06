# Documentation

View the online documentation at **https://tttapa.github.io/alpaqa**.

## Generate the documentation yourself

```sh
cd doxygen
doxygen
```

You'll need Doxygen and Dot (Graphviz).

## Generate the coverage information

```sh
mkdir build && cd build
cmake .. -DCMAKE_BUILD_TYPE=Coverage
make coverage
```

You'll need GCC or Clang, Google Test, LCOV and its dependencies.
