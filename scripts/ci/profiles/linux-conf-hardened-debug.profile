include({{ os.path.join(profile_dir, "linux-conf-hardened.profile") }})
[settings]
build_type=Debug
[conf]
tools.build:cflags+=["-U_FORTIFY_SOURCE", "-D_FORTIFY_SOURCE=0"]
tools.build:cxxflags+=["-U_FORTIFY_SOURCE", "-D_FORTIFY_SOURCE=0"]
[buildenv]
ALPAQA_PYTHON_DEBUG=1
