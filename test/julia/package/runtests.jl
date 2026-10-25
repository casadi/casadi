# What the CasADiNative package's Pkg.test runs (synced verbatim): its glue, casadi's tests and examples.
include(joinpath(@__DIR__, "package.jl"))
include(joinpath(@__DIR__, "..", "runtests.jl"))
include(joinpath(@__DIR__, "examples.jl"))
