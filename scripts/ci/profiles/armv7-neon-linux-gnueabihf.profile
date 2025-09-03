include({{ os.path.join(profile_dir, "cross-linux.profile") }})

[settings]
arch=armv7hf
os=Linux

[conf]
tools.build:cflags+=["-mfpu=neon", "-mfloat-abi=hard"]
tools.build:cxxflags+=["-mfpu=neon", "-mfloat-abi=hard"]
tools.cmake.cmaketoolchain:extra_variables*={"CPACK_DEBIAN_PACKAGE_ARCHITECTURE": "armhf"}

[options]
openblas/*:target=ARMV7
