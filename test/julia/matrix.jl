# P0-1 2-D construction, P0-3 hcat/vcat/hvcat, P1-5 conversion to native Julia.

@testset "2-D construction from Julia matrix (P0-1)" begin
    for T in (DM, SX, MX)
        m = T([1.0 2.0; 3.0 4.0])
        @test m isa T
        @test Base.size(m) == (2, 2)
    end
    # values land in the right cells (row-major literal -> correct casadi matrix)
    d = DM([1.0 2.0; 3.0 4.0])
    @test Float64(d[1, 1]) == 1.0
    @test Float64(d[1, 2]) == 2.0
    @test Float64(d[2, 1]) == 3.0
    @test Float64(d[2, 2]) == 4.0
    # vectors -> a column
    @test Base.size(DM([1.0, 2.0, 3.0])) == (3, 1)
    @test Base.size(MX([1.0, 2.0, 3.0])) == (3, 1)
end

@testset "hcat / vcat / hvcat (P0-3)" begin
    a = DM(1.0); b = DM(2.0); c = DM(3.0); e = DM(4.0)
    @test Base.size([a b]) == (1, 2)        # hcat
    @test Base.size([a; b]) == (2, 1)        # vcat
    M = [a b; c e]                            # hvcat
    @test Base.size(M) == (2, 2)
    @test Matrix{Float64}(M) == [1.0 2.0; 3.0 4.0]
    @test M isa DM                            # numeric cat stays DM
    # mixing in bare reals
    @test Matrix{Float64}([a 5.0]) == [1.0 5.0]
    @test Vector{Float64}([a; 5.0]) == [1.0, 5.0]
    # symbolic concatenation keeps its type and shape
    q = sym(SX, "q")
    @test [q q] isa SX
    @test Base.size([q q]) == (1, 2)
    @test Base.size([q; q]) == (2, 1)
end

@testset "function-form cat preserves type (R3-1)" begin
    # the EXPORTED horzcat/vertcat/diagcat (not [a b] syntax) keep DM/SX
    @test horzcat(DM(1.0), DM(2.0)) isa DM
    @test Base.size(horzcat(DM(1.0), DM(2.0))) == (1, 2)
    @test vertcat(DM(1.0), DM(2.0)) isa DM
    @test Base.size(vertcat(DM(1.0), DM(2.0))) == (2, 1)
    @test diagcat(DM(1.0), DM(2.0)) isa DM
    # no evalf-wrapper needed anymore for numeric concatenation
    @test Matrix{Float64}(vertcat(DM([1.0, 2.0]), DM([3.0]))) == Base.reshape([1.0, 2.0, 3.0], 3, 1)
    # SX stays SX
    a = sym(SX, "a"); b = sym(SX, "b")
    @test horzcat(a, b) isa SX
    @test vertcat(a, b) isa SX
end

@testset "2-D setindex! (R4-1)" begin
    d = DM(zeros(3, 3))
    d[2, 3] = 5.0                              # scalar into one cell
    @test Float64(d[2, 3]) == 5.0
    d[1, :] = [1.0 2.0 3.0]                     # row assignment from native matrix
    @test Matrix{Float64}(d)[1, :] == [1.0, 2.0, 3.0]
    d[:, 1] = DM([7.0, 8.0, 9.0])              # column from a DM
    @test Matrix{Float64}(d)[:, 1] == [7.0, 8.0, 9.0]
    d[2:3, 2] = [10.0, 11.0]                    # range index
    @test Float64(d[2, 2]) == 10.0
    @test Float64(d[3, 2]) == 11.0
    # the linear (1-index) setindex! still works
    e = DM(zeros(2, 2)); e[1] = 4.0
    @test Float64(e[1, 1]) == 4.0
    # SX 2-D setindex! routes too
    s = sym(SX, "s", 2, 2); s[1, 1] = sym(SX, "v")
    @test s isa SX
end

@testset "conversion to native Julia (P1-5)" begin
    d = DM([1.0 2.0; 3.0 4.0])
    @test Matrix{Float64}(d) == [1.0 2.0; 3.0 4.0]
    @test Matrix(d) == [1.0 2.0; 3.0 4.0]
    @test Array(d) == [1.0 2.0; 3.0 4.0]
    @test collect(d) == [1.0 2.0; 3.0 4.0]
    @test eltype(DM) == Float64
    # Vector{Float64} stays column-major nonzeros (documented back-compat)
    @test Vector{Float64}(d) == [1.0, 3.0, 2.0, 4.0]
    @test Float64(DM(7.0)) == 7.0
end
