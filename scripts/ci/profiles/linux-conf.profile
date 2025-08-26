[conf]
tools.build:skip_test=True
tools.build:cflags+=["-fdiagnostics-color"]
tools.build:cxxflags+=["-fdiagnostics-color"]
tools.build:exelinkflags+=["-flto=auto", "-static-libstdc++"]
tools.build:sharedlinkflags+=["-flto=auto", "-static-libstdc++"]
tools.cmake.cmaketoolchain:generator=Ninja Multi-Config
