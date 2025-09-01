include({{ os.path.join(profile_dir, "linux-conf.profile") }})
include({{ os.path.join(profile_dir, "alpaqa-python-linux.profile") }})

[conf]
tools.build:cflags+=["-fstack-protector-all", "-mshstk", "-fcf-protection=full", "-fdiagnostics-color"]
tools.build:cxxflags+=["-fstack-protector-all", "-mshstk", "-fcf-protection=full", "-fdiagnostics-color"]
tools.build:exelinkflags+=["-static-libgcc"]
tools.build:sharedlinkflags+=["-static-libgcc"]
