# Aqua.jl hygiene of the CasADiNative package; run in the env of make_package.jl ... Aqua.
using Aqua, CasADiNative, Test

@testset "Aqua" begin
    # [a 1 2] / [1; x]: the bracket forms take a casadi value anywhere, so their vararg union
    # also admits number/array-only calls; Base's more specific methods keep serving those.
    Aqua.test_all(CasADiNative; piracies = (treat_as_own = [Base.vcat, Base.hcat, Base.hvcat],))
end
