include({{ os.path.join(profile_dir, "cross-linux.profile") }})

[settings]
arch=armv8
os=Linux

[conf]
tools.build:cflags+=["-mcpu=cortex-a53+crc+simd"]
tools.build:cxxflags+=["-mcpu=cortex-a53+crc+simd"]
tools.cmake.cmaketoolchain:extra_variables*={"CPACK_DEBIAN_PACKAGE_ARCHITECTURE": "arm64"}

[options]
openblas/*:target=CORTEXA53
