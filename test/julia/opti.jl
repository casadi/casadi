# Opti stack: variable/parameter/minimize/subject_to/bounded/solve/value.

HAS_IPOPT && @testset "Opti stack" begin
    opti = Opti()
    x = variable(opti)
    y = variable(opti)
    minimize(opti, (x - 1)^2 + y^2)
    subject_to(opti, bounded(0, x, 0.5))   # instance-style bounded, bare numbers
    solver!(opti, "ipopt", Dict("ipopt" => Dict{String,Any}("print_level" => 0), "print_time" => false))
    sol = solve!(opti)
    @test Base.abs(value(sol, x) - 0.5) < 1e-6
    @test Base.abs(value(sol, y)) < 1e-6
end

@testset "comparison constraints + 2-D Opti" begin
    opti = Opti()
    X = variable(opti, 2, 3)
    @test Base.size(X[1, :]) == (1, 3)
    @test Base.size(X[:, 2]) == (2, 1)
    con = X[1, 1] == 0
    @test con isa MX
    con2 = X[2, 1] <= 1
    @test con2 isa MX
end

HAS_IPOPT && @testset "double-bounded constraint (python-style chaining)" begin
    opti = Opti()
    v = variable(opti)
    minimize(opti, (v - 3)^2)
    subject_to(opti, (0 <= v) <= 1)     # (a<=b)<=c forms a double bound
    solver!(opti, "ipopt", Dict("ipopt" => Dict{String,Any}("print_level" => 0), "print_time" => false))
    sol = solve!(opti)
    @test Base.abs(value(sol, v) - 1.0) < 1e-6
end

HAS_IPOPT && @testset "OptiSol/Opti stats() return a Julia Dict" begin
    opti = Opti()
    v = variable(opti)
    minimize(opti, (v - 3)^2)
    subject_to(opti, (0 <= v) <= 1)
    solver!(opti, "ipopt", Dict("ipopt" => Dict{String,Any}("print_level" => 0), "print_time" => false))
    sol = solve!(opti)
    st = sol.stats()                      # method-style on OptiSol
    @test st isa AbstractDict
    @test st["return_status"] == "Solve_Succeeded"
    @test st["success"] === true
    @test opti.stats()["return_status"] == "Solve_Succeeded"   # method-style on Opti
end
