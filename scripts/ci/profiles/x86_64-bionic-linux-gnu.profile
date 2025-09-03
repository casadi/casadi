include({{ os.path.join(profile_dir, "cross-linux.profile") }})

[settings]
arch=x86_64
os=Linux

[conf]
tools.build:cflags+=["-march=haswell"]
tools.build:cxxflags+=["-march=haswell"]
tools.cmake.cmaketoolchain:extra_variables*={"CPACK_DEBIAN_PACKAGE_ARCHITECTURE": "amd64"}

[options]
openblas/*:target=HASWELL
