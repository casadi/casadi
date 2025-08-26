include({{ os.path.join(profile_dir, "cross-linux.profile") }})

[settings]
arch=x86_64
os=Linux
os.toolchain-vendor=centos7
compiler.version=13

[conf]
tools.build:cflags+=["-march=haswell"]
tools.build:cxxflags+=["-march=haswell"]

[options]
openblas/*:target=HASWELL
