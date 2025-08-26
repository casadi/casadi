include(default)
[settings]
compiler.cppstd=23
[conf]
tools.build:skip_test=True
tools.build:cxxflags+=["/arch:AVX2"]
tools.build:cflags+=["/arch:AVX2"]
tools.cmake.cmaketoolchain:generator=Ninja Multi-Config
tools.build:compiler_executables*={"fortran": "FC-NOTFOUND" }
