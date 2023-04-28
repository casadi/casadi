cd "$( dirname "${BASH_SOURCE[0]:-${(%):-%x}}" )"/../..

# Determine architecture
case $(uname -m) in
    x86_64) triple="x86_64-centos7-linux-gnu" ;;
    armv6l) triple="armv6-rpi-linux-gnueabihf" ;;
    armv7l) triple="armv7-neon-linux-gnueabihf" ;;
    aarch64) triple="aarch64-rpi3-linux-gnu" ;;
esac
tools_dir="$PWD/toolchains"
pfx="$tools_dir/$triple"

export PATH="$PWD/staging/bin:$PATH"
export LD_LIBRARY_PATH="$pfx/x-tools/$triple/$triple/lib64:$PWD/staging/lib:$LD_LIBRARY_PATH"
export LD_PRELOAD="$pfx/x-tools/$triple/$triple/lib64/libstdc++.so.6.0.30:$pfx/x-tools/$triple/$triple/lib64/libgfortran.so.5.0.0:$LD_PRELOAD"
