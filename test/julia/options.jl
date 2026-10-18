# options.jl -- option values keep their Julia type on the way into GenericType (#4445).

@testset "options: whole numbers keep their type" begin
    x = sym(SX, "x")
    F = integrator("F", "rk", Dict("x" => x, "ode" => -x), 0.0, 1.0, Dict("number_of_finite_elements" => 7))
    h = 1 / 7
    @test Float64(F(Dict("x0" => DM(1.0)))["xf"]) ≈ (1 - h + h^2 / 2 - h^3 / 6 + h^4 / 24)^7 rtol = 1e-12
    F = integrator("F", "rk", Dict("x" => x, "ode" => -x), 0.0, 1.0, Dict("number_of_finite_elements" => 7.0, "expand" => true))
    @test Float64(F(Dict("x0" => DM(1.0)))["xf"]) ≈ (1 - h + h^2 / 2 - h^3 / 6 + h^4 / 24)^7 rtol = 1e-12
    A = sym(MX, "A", 2, 2)
    for margin in [2.0, 2.5]
        G = C.Function("G", [A], [C.convexify(A, Dict{String,Any}("strategy" => "regularize", "margin" => margin))])
        @test Matrix(G(DM(zeros(2, 2)))) == [margin 0.0; 0.0 margin]
    end
    if C.has_rootfinder("kinsol")
        z, p = sym(SX, "z", 2), sym(SX, "p")
        R = rootfinder("R", "kinsol", Dict("x" => z, "p" => p, "g" => vertcat(z[1]^3 - p, z[2] - 2p)), Dict("u_scale" => [2, 3], "print_level" => 0))
        @test vec(Matrix(R(Dict("x0" => DM([1.0, 1.0]), "p" => 8.0))["x"])) ≈ [2.0, 16.0]
    end
    if HAS_IPOPT
        y = sym(SX, "y", 2)
        S = nlpsol("S", "ipopt", Dict("x" => y, "f" => (y[1] - 1)^2 + 100 * (y[2] - y[1]^2)^2),
                   Dict("ipopt" => Dict("max_iter" => 5, "print_level" => 0), "print_time" => false))
        S(Dict("x0" => DM([-1.2, 1.0])))
        @test C.stats(S)["iter_count"] == 5
    end
end
