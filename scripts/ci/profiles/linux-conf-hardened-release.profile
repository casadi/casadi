include({{ os.path.join(profile_dir, "linux-conf-hardened.profile") }})
[settings]
build_type=Release
[conf]
tools.build:cflags+=["-U_FORTIFY_SOURCE", "-D_FORTIFY_SOURCE=3"]
tools.build:cxxflags+=["-U_FORTIFY_SOURCE", "-D_FORTIFY_SOURCE=3"]
