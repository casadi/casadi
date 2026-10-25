# Glue of the package front matter (swig/julia/package): plugin path and BLAS set up by __init__.
using CasADiNative
using Test
import LinearAlgebra

@testset "CasADiNative package" begin
    @testset "plugin path from __init__" begin
        dir = CasADiNative._plugin_dir()
        @test first(split(ENV["CASADI_PLUGIN_SEARCH_PATH"], Sys.iswindows() ? ';' : ':')) == dir
        @test CasADiNative.C.has_nlpsol("sqpmethod")
    end

    @testset "LP64 BLAS behind libblastrampoline" begin
        @test any(lib -> lib.interface == :lp64, LinearAlgebra.BLAS.lbt_get_config().loaded_libs)
    end

    @testset "ipopt converges" begin
        x = SX.sym("x"); y = SX.sym("y")
        nlp = Dict("x" => [x; y], "f" => (1 - x)^2 + 100 * (y - x^2)^2)
        S = nlpsol("S", "ipopt", nlp, Dict{String,Any}("print_time" => false,
            "ipopt" => Dict{String,Any}("print_level" => 0)))
        @test Vector{Float64}(S(x0 = [0.0, 0.0])["x"]) ≈ [1.0, 1.0] atol = 1e-6
    end
end
