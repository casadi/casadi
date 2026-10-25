# Build the CasADiNative package (what the registry ships) from a build tree, in a dev
# environment where `import CasADiNative` resolves to it, and precompile it strictly.
#
# Usage: julia make_package.jl <module dir> <env dir> [extra packages...]
# <module dir> holds the build's CasADiNative.jl + casadi.jl + libcasadi_wrap
# (build/swig/julia/CasADiNative). Creates the project <env dir> (run with
# --project=<env dir> or JULIA_PROJECT), developing <env dir>/CasADiNative, plus the
# extra registered packages (e.g. Aqua). src/ is the registered package's front matter
# (swig/julia/package) verbatim, except bindings.jl, which points at the build tree:
# the layer without its module wrapper (as CasADiSWIG_jll's casadi_native_body.jl) and
# the plugin dir (first entry of CASADI_PLUGIN_SEARCH_PATH).
using Pkg

src = abspath(ARGS[1])
out = abspath(ARGS[2])
pkg = joinpath(out, "CasADiNative")
front = joinpath(@__DIR__, "..", "..", "..", "swig", "julia", "package")
plugins = first(split(get(ENV, "CASADI_PLUGIN_SEARCH_PATH", ""), ':'))
isdir(plugins) || error("CASADI_PLUGIN_SEARCH_PATH must start with the build's plugin dir")
# the layer without its module wrapper, as CasADiSWIG_jll's casadi_native_body.jl
layer = read(joinpath(src, "CasADiNative.jl"), String)
body = replace(layer, r"^module CasADiNative\n"m => "", r"^end # module CasADiNative\n"m => "")
length(body) == length(layer) - length("module CasADiNative\n") - length("end # module CasADiNative\n") ||
    error("module wrapper of CasADiNative.jl not found")
bindings = """
    # Binary locations: the casadi build tree (written by test/julia/package/make_package.jl).
    _plugin_dir() = $(repr(plugins))
    const _BODY = joinpath(@__DIR__, "casadi_native_body.jl")
    """
function srcdir(dir)
    mkpath(dir)
    for f in readdir(src)
        (f == "casadi.jl" || startswith(f, "libcasadi_wrap")) && cp(joinpath(src, f), joinpath(dir, f); force=true)
    end
    write(joinpath(dir, "casadi_native_body.jl"), body)
    cp(joinpath(front, "CasADiNative.jl"), joinpath(dir, "CasADiNative.jl"); force=true)
    write(joinpath(dir, "bindings.jl"), bindings)
end
srcdir(joinpath(pkg, "src"))
write(joinpath(pkg, "Project.toml"), """
    name = "CasADiNative"
    uuid = "30fb3c3e-f728-433f-984c-e53dd5973698"
    version = "0.1.0"

    [deps]
    LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
    OpenBLAS32_jll = "656ef2d0-ae68-5445-9ca0-591084a874a2"

    [compat]
    LinearAlgebra = "1"
    OpenBLAS32_jll = "0.3.10"
    julia = "1.10"
    """)

Pkg.activate(out)
Pkg.develop(path=pkg)
Pkg.add(["OpenBLAS32_jll"; ARGS[3:end]])
Pkg.precompile(strict=true)
Base.isprecompiled(Base.identify_package("CasADiNative")) || error("CasADiNative failed to precompile")
