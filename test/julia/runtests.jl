# test/julia/runtests.jl -- thin runner for the SWIG-generated Julia bindings.
# Mirrors the wasm-js multi-file layout: this file just loads CasADiNative and
# includes the per-area test files. Needs CasADiNative importable (the environment
# test/julia/package/make_package.jl builds from a build tree) and casadi plugins on
# the search path.
using Test
using CasADiNative
import Serialization
const C = CasADiNative.C

const HAS_IPOPT = try
    C.has_nlpsol("ipopt")
catch
    false
end
HAS_IPOPT || @info "ipopt plugin not found: solver testsets will be skipped"

# Auto-discover per-area test files: every *.jl here except this runner and
# helpers (prefixed with _). Drop a new file in and it runs -- no edits here.
@testset "casadi julia bindings" begin
    _files = sort(filter(readdir(@__DIR__)) do f
        endswith(f, ".jl") && f != "runtests.jl" && !startswith(f, "_")
    end)
    for f in _files
        include(joinpath(@__DIR__, f))
    end
end
