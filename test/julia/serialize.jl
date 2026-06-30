# Serialization round-trips: high-level casadi_serialize/casadi_deserialize and
# the unpack() type-tag dispatcher (mirrors Python/JS/MATLAB unpack()).

@testset "scalar round-trips" begin
    x = sym(SX, "x")
    @test C.str(casadi_deserialize(casadi_serialize(x))) == C.str(x)

    m = sym(MX, "m")
    @test C.str(casadi_deserialize(casadi_serialize(m))) == C.str(m)

    # matrix-shaped DM (not just scalar)
    d = DM([1.0 2.0; 3.0 4.0])
    rd = casadi_deserialize(casadi_serialize(d))
    @test rd isa DM
    @test Base.Matrix(rd) == [1.0 2.0; 3.0 4.0]

    sp = C.Sparsity(3, 2)
    rsp = casadi_deserialize(casadi_serialize(sp))
    @test rsp isa Sparsity
    @test (C.size1(rsp), C.size2(rsp)) == (3, 2)
end

@testset "primitive round-trips" begin
    @test casadi_deserialize(casadi_serialize(7)) == 7
    @test casadi_deserialize(casadi_serialize(3.5)) == 3.5
    @test casadi_deserialize(casadi_serialize("hello")) == "hello"
end

@testset "Function round-trip then call" begin
    y = sym(SX, "y")
    f = C.Function("f", SX[y], SX[sin(y)], Dict{String,Any}())
    g = casadi_deserialize(casadi_serialize(f))
    @test g isa C.Function
    # reconstructed function must still compute correctly
    @test C.str(g(0.5)) == C.str(f(0.5))
    @test Float64(g(0.5)) ≈ sin(0.5)
end

@testset "vector round-trips" begin
    sv = casadi_deserialize(casadi_serialize(SX[sym(SX, "a"), sym(SX, "b")]))
    @test [C.str(e) for e in sv] == ["a", "b"]

    mv = casadi_deserialize(casadi_serialize(MX[sym(MX, "p"), sym(MX, "q")]))
    @test [C.str(e) for e in mv] == ["p", "q"]

    # DM vector exercises the slot-pinning fix (generated dispatch mis-routes it)
    dv = casadi_deserialize(casadi_serialize(DM[DM(1.0), DM(2.5)]))
    @test [Float64(e) for e in dv] == [1.0, 2.5]
end

@testset "explicit unpack via deserializer" begin
    # the low-level pair: pack -> encode -> StringDeserializer -> unpack
    s = StringSerializer(Dict{String,Any}())
    pack(s, sym(SX, "w"))
    pack(s, DM(9.0))
    d = StringDeserializer(encode(s))
    @test C.str(unpack(d)) == "w"          # first object
    @test Float64(unpack(d)) == 9.0        # second object
end

@testset "malformed input errors (no segfault)" begin
    @test_throws Exception casadi_deserialize("garbage-not-valid")
    @test_throws Exception casadi_deserialize("")
end
