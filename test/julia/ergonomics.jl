# R3 ergonomics: typed cat/dot, 3-arg convenience defaults, method-style calls.

@testset "typed cat preserves DM (R3-1)" begin
    # function-form horzcat/vertcat/diagcat/blockcat stay DM when all args DM/Real
    h = horzcat(DM(1.0), DM(2.0))
    @test h isa DM
    @test Base.size(h) == (1, 2)
    @test Matrix{Float64}(h) == [1.0 2.0]
    v = vertcat(DM(1.0), DM(2.0), DM(3.0))
    @test v isa DM
    @test Base.size(v) == (3, 1)
    dc = diagcat(DM([1.0 2.0; 3.0 4.0]), DM(5.0))
    @test dc isa DM
    @test Base.size(dc) == (3, 3)
    bc = blockcat(DM(1.0), DM(2.0), DM(3.0), DM(4.0))
    @test bc isa DM
    @test Matrix{Float64}(bc) == [1.0 2.0; 3.0 4.0]
    # mixing bare reals stays DM
    @test horzcat(DM(1.0), 2.0) isa DM
    # single-vector call form still works and stays DM
    @test horzcat(DM[DM(1.0), DM(2.0)]) isa DM
    @test vertcat(DM[DM(1.0), DM(2.0)]) isa DM
end

@testset "typed cat preserves SX, promotes to MX (R3-1)" begin
    a = sym(SX, "a"); b = sym(SX, "b")
    @test horzcat(a, b) isa SX
    @test vertcat(a, b) isa SX
    @test diagcat(a, b) isa SX
    @test horzcat(a, 1.0) isa SX           # SX + Real -> SX
    m = sym(MX, "m")
    @test horzcat(m, DM(1.0)) isa MX        # any MX present -> MX
    @test vertcat(m, m) isa MX
end

@testset "dot handles mixed types (R3-1)" begin
    d = dot(DM([1.0, 2.0]), DM([3.0, 4.0]))
    @test d isa DM
    @test Float64(d) == 11.0
    # dot(DM, SX): previously no method; now promotes DM->SX
    s = sym(SX, "s", 2)
    ds = dot(DM([1.0, 2.0]), s)
    @test ds isa SX
    # dot(SX, DM) symmetric
    @test dot(s, DM([1.0, 2.0])) isa SX
end

@testset "3-arg convenience defaults (R3-2)" begin
    x = sym(SX, "x"); y = sym(SX, "y")
    z = vertcat(x, y)
    # Function ctor without trailing opts dict
    f = C.Function("f", [x], [sin(x)])
    @test f isa C.Function
    @test C.n_in(f) == 1
    # gradient / tangent / jacobian / hessian without opts dict
    g = gradient(x^2 + y, z)
    @test g isa SX
    @test Base.size(g) == (2, 1)
    t = tangent(x^2, x)
    @test t isa SX
    j = jacobian(vertcat(x^2, y), z)
    @test Base.size(j) == (2, 2)
    H = hessian(x^2 + y^2, z)
    @test Base.size(H) == (2, 2)
end

@testset "factory convenience: default + aux Dict (R3-2)" begin
    x = sym(SX, "x")
    f = C.Function("f", Dict("x" => x, "f" => x^2, "g" => sin(x)),
                       ["x"], ["f", "g"], Dict{String,Any}())
    # default aux/opts: just request a jacobian output
    fa = factory(f, "fa", ["x"], ["f", "jac:f:x"])
    @test fa isa C.Function
    @test C.n_out(fa) == 2
    # aux as a Julia Dict marshalled to dict-of-string-lists (positional, generic typemap)
    fb = factory(f, "fb", ["x"], ["jac:f:x"], Dict("combo" => ["f", "g"]))
    @test fb isa C.Function
end

@testset "method-style calls on Function (R3-3)" begin
    x = sym(SX, "x")
    f = C.Function("f", [x], [sin(x)])
    # introspection method-style mirrors C.<name>(f)
    @test f.name() == "f"
    @test f.n_in() == 1
    @test f.n_out() == 1
    @test f.name_in() == ["i0"]
    # AD / transforms method-style
    jf = f.jacobian()
    @test jf isa C.Function
    @test C.n_out(jf) == 1
    fwd = f.forward(1)
    @test fwd isa C.Function
    mp = f.map(3)
    @test mp isa C.Function
    # field access still works and is NOT intercepted
    @test f.ptr isa Ptr
    @test f.ptr != C.C_NULL
    # a non-allowlisted property falls through to getfield (errors, not a shim)
    @test_throws Exception f.this_is_not_a_method
end

@testset "method-style equals free-function form (R3-3)" begin
    x = sym(SX, "x")
    f = C.Function("f", [x], [x^2])
    # f.jacobian() and C.jacobian(f) build the same-shaped function
    @test C.n_in(f.jacobian()) == C.n_in(C.jacobian(f))
    @test C.n_out(f.jacobian()) == C.n_out(C.jacobian(f))
end

HAS_IPOPT && @testset "method-style on OptiSol (R3-3)" begin
    o = C.Opti()
    xx = variable(o)
    minimize(o, (xx - 3.0)^2)
    solver!(o, "ipopt", Dict("print_time" => false, "ipopt" => Dict("print_level" => 0, "sb" => "yes")))
    sol = solve!(o)
    @test sol isa C.OptiSol
    # sol.value(xx) method-style mirrors C.value(sol, xx, ...)
    @test Base.abs(Float64(sol.value(xx, MX[])) - 3.0) < 1e-6
    # field access preserved
    @test sol.ptr isa Ptr
end

@testset "generate default opts + external 2-arg (R4-2,3)" begin
    x = sym(SX, "x")
    f = C.Function("f", [x], [x^2 + 1.0], ["x"], ["y"])
    mktempdir() do dir
        cd(dir) do
            # generate(f, fname) with no spelled-out Dict
            cname = generate(f, "f.c")
            @test Base.isfile(cname)
            @test Base.endswith(cname, ".c")
            # compile + load back via the 2-arg external(name, path)
            run(`gcc -fPIC -shared -O0 $cname -o libf.so`)
            g = external("f", "./libf.so")     # 2-arg form (no trailing Dict)
            @test g isa C.Function
            @test Float64(g(3.0)) == 10.0       # 3^2 + 1
            # method-style f.generate also defaults opts
            cname2 = f.generate("f2.c")
            @test Base.isfile(cname2)
        end
    end
end

@testset "named-IO Function via positional names (R4-4)" begin
    x = sym(SX, "x"); y = sym(SX, "y")
    f = C.Function("f", [x, y], [x + y, x * y], ["a", "b"], ["sum", "prod"])
    @test f isa C.Function
    @test C.name_in(f) == ["a", "b"]
    @test C.name_out(f) == ["sum", "prod"]
    # call by kwargs uses the declared input names
    r = f(a=2.0, b=3.0)
    @test Float64(r["sum"]) == 5.0
    @test Float64(r["prod"]) == 6.0
    # plain 3-arg still defaults opts and auto-names
    g = C.Function("g", [x], [sin(x)])
    @test g.name_in() == ["i0"]
    h = C.Function("h", [x], [x], ["q"], ["w"])
    @test C.name_in(h) == ["q"]
end

@testset "save default opts + method-style (R4-5)" begin
    x = sym(SX, "x")
    f = C.Function("f", [x], [x^2])
    mktempdir() do dir
        fname = joinpath(dir, "f.casadi")
        save(f, fname)                          # no trailing Dict
        @test Base.isfile(fname)
        f2 = C.load(C.Function, fname)
        @test C.n_in(f2) == 1
        # method-style f.save also defaults opts
        fname2 = joinpath(dir, "f2.casadi")
        f.save(fname2)
        @test Base.isfile(fname2)
    end
end

@testset "no method ambiguities introduced (R3)" begin
    @test Base.isempty(Test.detect_ambiguities(CasADiNative; recursive=false))
end
