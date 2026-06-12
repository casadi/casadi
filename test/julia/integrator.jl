# integrator.jl -- capability suite (wasm-js only had `rk`; this goes far beyond).
# Covers rk/collocation/cvodes/idas on known ODEs/DAEs, quadratures, grid
# output, parameters, and forward sensitivity via the integrator's jacobian.

const _IN_OD = Dict{String,Any}()

# x' = -x over [0, 1] -> x(1) = e^-1, for every ODE integrator plugin.
for plugin in ["rk", "collocation", "cvodes", "idas"]
    (try C.has_integrator(plugin) catch; false end) || continue
    @testset "integrator $plugin: xdot = -x" begin
        x = sym(SX, "x")
        F = integrator("F", plugin, Dict("x" => x, "ode" => -x), 0.0, 1.0)
        r = F(Dict("x0" => DM(1.0)))
        @test abs(Float64(r["xf"]) - exp(-1.0)) < 1e-4
    end
end

# x' = x over [0, 1] -> x(1) = e (exponential growth).
@testset "integrator cvodes: xdot = x -> e" begin
    x = sym(SX, "x")
    F = integrator("F", "cvodes", Dict("x" => x, "ode" => x), 0.0, 1.0)
    @test abs(Float64(F(Dict("x0" => DM(1.0)))["xf"]) - exp(1.0)) < 1e-4
end

# Quadrature: x' = -x, q' = x ; q(1) = 1 - e^-1.
@testset "integrator cvodes: quadrature" begin
    x = sym(SX, "x")
    F = integrator("F", "cvodes", Dict("x" => x, "ode" => -x, "quad" => x), 0.0, 1.0)
    r = F(Dict("x0" => DM(1.0)))
    @test abs(Float64(r["xf"]) - exp(-1.0)) < 1e-4
    @test abs(Float64(r["qf"]) - (1.0 - exp(-1.0))) < 1e-4
end

# Parameter: x' = -p*x ; x(1) = e^{-p}.  Vary p over two calls.
@testset "integrator cvodes: parametric" begin
    x = sym(SX, "x"); p = sym(SX, "p")
    F = integrator("F", "cvodes", Dict("x" => x, "p" => p, "ode" => -p * x), 0.0, 1.0)
    @test abs(Float64(F(Dict("x0" => DM(1.0), "p" => DM(2.0)))["xf"]) - exp(-2.0)) < 1e-4
    @test abs(Float64(F(Dict("x0" => DM(1.0), "p" => DM(0.5)))["xf"]) - exp(-0.5)) < 1e-4
end

# Grid output: collect the full trajectory of x' = -x at t = 0,0.25,...,1.
@testset "integrator cvodes: grid output" begin
    x = sym(SX, "x")
    grid = collect(0.0:0.25:1.0)
    F = C.integrator("F", "cvodes", Dict{String,Any}("x" => x, "ode" => -x), 0.0, grid, _IN_OD)
    traj = Vector{Float64}(F(Dict("x0" => DM(1.0)))["xf"])
    @test length(traj) == 5
    @test maximum(abs.(traj .- exp.(-grid))) < 1e-4
end

# 2-state linear ODE (harmonic oscillator): x'=v, v'=-x ; period 2*pi.
@testset "integrator cvodes: 2-state oscillator" begin
    s = sym(SX, "s", 2)        # s = [x, v]
    ode = vertcat(s[2], -s[1])
    F = integrator("F", "cvodes", Dict("x" => s, "ode" => ode), 0.0, pi / 2)
    r = Vector{Float64}(F(Dict("x0" => DM([1.0, 0.0])))["xf"])
    # quarter period: [cos(pi/2), -sin(pi/2)] = [0, -1]
    @test maximum(abs.(r .- [0.0, -1.0])) < 1e-4
end

# DAE via idas: x' = z, 0 = z - x  -> x(t) = e^t.
@testset "integrator idas: DAE x'=z, 0=z-x" begin
    x = sym(SX, "x"); z = sym(SX, "z")
    F = integrator("F", "idas",
        Dict("x" => x, "z" => z, "ode" => z, "alg" => z - x), 0.0, 1.0)
    r = F(Dict("x0" => DM(1.0), "z0" => DM(1.0)))
    @test abs(Float64(r["xf"]) - exp(1.0)) < 1e-3
    @test abs(Float64(r["zf"]) - exp(1.0)) < 1e-3
end

# Forward sensitivity via the integrator's derivative Function.
# x' = -p*x ; d(xf)/d(p) at p=1, x0=1 is -1 * e^{-1}.
@testset "integrator cvodes: forward sensitivity (jacobian)" begin
    x = sym(SX, "x"); p = sym(SX, "p")
    F = integrator("F", "cvodes", Dict("x" => x, "p" => p, "ode" => -p * x), 0.0, 1.0)
    JF = C.jacobian(F)          # derivative Function of the integrator
    @test JF isa C.Function
    # Locate the d(xf)/d(p) output by name.
    onames = C.name_out(JF)
    idx = findfirst(n -> occursin("xf", n) && occursin("p", n), onames)
    @test idx !== nothing
    # Build a call dict from input names, supplying x0, p (others default to 0).
    args = Dict{String,Any}()
    for nm in C.name_in(JF)
        args[nm] = nm == "p" ? DM(1.0) : (nm == "x0" ? DM(1.0) : DM(0.0))
    end
    out = C.call(JF, args)
    dxf_dp = Float64(out[onames[idx]])
    @test abs(dxf_dp - (-exp(-1.0))) < 1e-3
end
