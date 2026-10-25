# overloads.jl -- C++ overloads sharing a Julia signature: resolution order and the
# "no matching overload" error.

@testset "overloads sharing a signature: no match lists the prototypes" begin
    x = sym(SX, "x"); y = sym(MX, "y")
    err = try C.Function("f", [x], [y]) catch e e end
    @test err isa C.SwigError
    @test occursin("Wrong number or type of arguments for overloaded function 'Function'", err.msg)
    @test occursin("Function(String,Vector{SX},Vector{SX},Dict)", err.msg)
    @test occursin("Function(String,Vector{MX},Vector{MX},Dict)", err.msg)
    @test occursin("You have: Function(::String, ::Vector{SX}, ::Vector{MX})", err.msg)
    err = try C.Function("f", [1.5, "a"], [x], Dict{String,Any}()) catch e e end
    @test occursin("You have: Function(::String, ::Vector{Any}[Float64, String], ::Vector{SX}, ::Dict{String, Any})", err.msg)
    err = try C.substitute(x, y, x) catch e e end
    @test occursin("overloaded function 'substitute'", err.msg) && occursin("substitute(MX,MX,MX)", err.msg)
    err = try C.nlpsol("s", "ipopt", 3) catch e e end
    @test occursin("nlpsol(String,String,Dict{String,SX},Dict)", err.msg)
    @test occursin("You have: nlpsol(::String, ::String, ::Int64)", err.msg)
    err = try C.call(C.Function("f", [x], [x]), "a") catch e e end
    @test occursin("You have: call(::Function, ::String)", err.msg)
    # a match still resolves: SX, MX and DM through the same method
    @test C.n_in(C.Function("f", [x], [x])) == 1
    @test C.n_in(C.Function("f", [y], [y])) == 1
    @test C.substitute(x, x, DM(3.0)) isa SX && C.substitute(y, y, DM(3.0)) isa MX
end

@testset "typed containers pick their overload by Julia dispatch" begin
    x = sym(SX, "x"); y = sym(MX, "y")
    sig(args...) = which(C.Function, typeof.(args)).sig.parameters
    # Vector{SX} / Vector{MX} reach their own overload, not the runtime-dispatched method
    @test Vector{SX} <: sig("f", [x], [x])[3] && !(Vector{MX} <: sig("f", [x], [x])[3])
    @test Vector{MX} <: sig("f", [y], [y])[3] && !(Vector{SX} <: sig("f", [y], [y])[3])
    @test sig("f", Any[x], Any[x])[3] === Any
    # DM and number vectors fit SX and MX alike: the SX overload, as SWIG's precedence orders them
    @test C.is_a(C.Function("f", DM[], [DM(1)]), "SXFunction")
    @test C.is_a(C.Function("f", [x], [1.0]), "SXFunction")
    @test C.is_a(C.Function("f", [y], [1.0]), "MXFunction")
    # Vector{Any} still resolves on its contents
    @test C.is_a(C.Function("f", Any[y], Any[y]), "MXFunction")
    @test C.is_a(C.Function("f", Any[], Any[]), "SXFunction")
end

@testset "overloads sharing a signature: the C++ overload order" begin
    # Slice(int, int) before Slice(casadi_int, casadi_int) only while the values fit an int
    @test (C.Slice(1, 5).start, C.Slice(1, 5).stop) == (1, 5)
    s = C.Slice(0, 2^40)
    @test (s.start, s.stop, s.step) == (0, 2^40, 1)
    s = C.Slice(2^40, 2^41, 2)
    @test (s.start, s.stop, s.step) == (2^40, 2^41, 2)
    @test_throws C.SwigError C.Slice(1, big(2)^70)
end
