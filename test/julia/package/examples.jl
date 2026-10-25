# docs/examples/julia, each in a fresh process of the active project; fails on a method overwrite.
using CasADiNative
using Test

# Examples needing a plugin or a host tool that may be absent (e.g. qpoases is not in the JLLs).
has_gcc() = !Sys.iswindows() && Sys.which("gcc") !== nothing
const REQUIRES = Dict(
    "chain_qp.jl" => () -> CasADiNative.C.has_conic("qpoases"),
    "simple_lp.jl" => () -> CasADiNative.C.has_conic("qpoases"),
    "c_code_generation.jl" => has_gcc,
    "nlp_codegen.jl" => has_gcc,
)

@testset "casadi examples" begin
    dir = joinpath(@__DIR__, "..", "..", "..", "docs", "examples", "julia")
    for f in sort(filter(endswith(".jl"), readdir(dir)))
        @testset "$f" begin
            if !get(REQUIRES, f, () -> true)()
                @test_skip f
            else
                out = IOBuffer()
                cmd = `$(Base.julia_cmd()) --startup-file=no --warn-overwrite=yes --project=$(Base.active_project()) $(joinpath(dir, f))`
                ok = success(pipeline(Cmd(cmd; dir = mktempdir()); stdout = out, stderr = out))
                log = String(take!(out))
                ok &= !occursin(r"WARNING: Method definition .* overwritten", log)
                ok || print(log)
                @test ok
            end
        end
    end
end
