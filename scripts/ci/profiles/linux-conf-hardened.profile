include({{ os.path.join(profile_dir, "linux-conf.profile") }})

[conf]
tools.build:cflags+=["-fstack-protector-all", "-mshstk", "-fcf-protection=full", "-fdiagnostics-color"]
tools.build:cxxflags+=["-fstack-protector-all", "-mshstk", "-fcf-protection=full", "-fdiagnostics-color"]
tools.build:exelinkflags+=["-static-libgcc"]
tools.build:sharedlinkflags+=["-static-libgcc"]

[options]
# Used in docs
alpaqa/*:with_cutest=True
