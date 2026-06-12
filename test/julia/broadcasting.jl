# P0-2: every elementwise broadcast op routes to a CasADi op.

@testset "broadcast arithmetic (P0-2)" begin
    A = sym(SX, "A", 2, 2)
    @test Base.size(A .+ A) == (2, 2)
    @test Base.size(A .- A) == (2, 2)
    @test Base.size(A .* A) == (2, 2)
    @test Base.size(A ./ A) == (2, 2)
    @test Base.size(A .^ 2) == (2, 2)
    # .^2 must square the VALUES (not merely stringify to "sq"/"pow")
    N = DM([1.0 2.0; 3.0 4.0])
    @test C.str(N .^ 2) == C.str(DM([1.0 4.0; 9.0 16.0]))
end

@testset "broadcast math functions (P0-2)" begin
    A = sym(SX, "A", 2, 2)
    @test Base.size(sin.(A)) == (2, 2)
    @test Base.size(exp.(A)) == (2, 2)
    @test Base.size(sqrt.(A)) == (2, 2)
    x = sym(SX, "x")
    @test Base.string(sin.(x)) == Base.string(sin(x))   # scalar broadcast == scalar op
end

@testset "broadcast scalar/array mixing (P0-2)" begin
    x = sym(SX, "x")
    @test Base.string(x .+ 1) == "(x+1)"
    A = sym(SX, "A", 2, 2)
    @test Base.size(1 .+ A) == (2, 2)
    @test Base.size(A .+ 1) == (2, 2)
    @test Base.size(A .* DM(2.0)) == (2, 2)        # mixed SX/DM
    v = sym(SX, "v", 2)
    @test Base.string(v .+ [1.0, 2.0]) == "[(v_0+1), (v_1+2)]"   # SX .+ Julia array
    # numeric DM broadcasts evaluate
    d = DM([1.0, 2.0, 3.0])
    @test Vector{Float64}(d .+ 1) == [2.0, 3.0, 4.0]
    @test Vector{Float64}(2 .* d) == [2.0, 4.0, 6.0]
end

@testset "broadcast comparisons (P0-2)" begin
    x = sym(SX, "x")
    @test Base.string(x .< 1) == "(x<1)"
    @test Base.string(x .<= 1) == "(x<=1)"
end

@testset "broadcast .* is elementwise, not mtimes" begin
    A = DM([1.0 2.0; 3.0 4.0]); B = DM([1.0 0.0; 0.0 1.0]); D = DM([2.0 0.0; 0.0 2.0])
    # .* must be Hadamard (elementwise), distinct from * (matrix product)
    @test C.str(A .* B) == C.str(DM([1.0 0.0; 0.0 4.0]))
    @test C.str(A * B)  == C.str(A)               # mtimes by identity = A
    # fused chains must stay elementwise at every node (regression: flatten composed mtimes)
    @test C.str(A .* B .* D)  == C.str(DM([2.0 0.0; 0.0 8.0]))
    @test C.str(A .* B .+ D)  == C.str(DM([3.0 0.0; 0.0 6.0]))
    @test C.str(2 .* A .* B)  == C.str(DM([2.0 0.0; 0.0 8.0]))
end
