# nlpsol convenience wrapper (ipopt): rosenbrock-style and kwarg call forms.

HAS_IPOPT && @testset "nlpsol rosenbrock (ipopt)" begin
    x = sym(SX, "x"); y = sym(SX, "y"); z = sym(SX, "z")
    S = nlpsol("S", "ipopt",
               Dict("x" => vertcat(x, y, z), "f" => x^2 + 100 * z^2,
                    "g" => z + (1 - x)^2 - y),
               Dict("ipopt" => Dict{String,Any}("print_level" => 0), "print_time" => false))
    sol = S(Dict("x0" => DM([2.5, 3.0, 0.75]), "lbg" => DM(0.0), "ubg" => DM(0.0)))
    xs = sol["x"]
    @test Base.abs(Float64(xs[2]) - 1.0) < 1e-8
    @test Base.abs(Float64(xs[1])) < 1e-8
    @test Float64(sol["f"]) < 1e-10
end

HAS_IPOPT && @testset "kwarg solver call" begin
    x = sym(SX, "x"); y = sym(SX, "y"); z = sym(SX, "z")
    S = nlpsol("S", "ipopt",
               Dict("x" => vertcat(x, y, z), "f" => x^2 + 100 * z^2,
                    "g" => z + (1 - x)^2 - y),
               Dict("ipopt" => Dict{String,Any}("print_level" => 0), "print_time" => false))
    sol = S(x0 = [2.5, 3.0, 0.75], lbg = 0, ubg = 0)   # bare numbers/arrays
    @test Base.abs(Float64(sol["x"][2]) - 1.0) < 1e-8
end

HAS_IPOPT && @testset "stats() returns a Julia Dict" begin
    x = sym(SX, "x"); y = sym(SX, "y"); z = sym(SX, "z")
    S = nlpsol("S", "ipopt",
               Dict("x" => vertcat(x, y, z), "f" => x^2 + 100 * z^2,
                    "g" => z + (1 - x)^2 - y),
               Dict("ipopt" => Dict{String,Any}("print_level" => 0), "print_time" => false))
    S(Dict("x0" => DM([2.5, 3.0, 0.75]), "lbg" => DM(0.0), "ubg" => DM(0.0)))
    st = C.stats(S)                       # free-function form
    @test st isa AbstractDict
    @test st["return_status"] == "Solve_Succeeded"   # string GenericType value
    @test st["success"] === true                     # bool GenericType value
    @test st["iter_count"] isa Integer               # int GenericType value
    @test S.stats() == st                            # method-style accessor
end
