# MX symbolics path.

@testset "MX path" begin
    xm = sym(MX, "x", 2)
    @test xm isa MX
    fm = C.Function("fm", [xm], [xm + xm], Dict{String,Any}())
    r = fm(DM([1.0, 2.0]))
    @test Vector{Float64}(r) == [2.0, 4.0]
    @test Base.size(jacobian(xm + xm, xm)) == (2, 2)
    # 2*q must scale the VALUE (not merely stringify non-empty)
    q = sym(MX, "q")
    gq = C.Function("gq", [q], [2 * q], Dict{String,Any}())
    @test Vector{Float64}(gq(3.0)) == [6.0]
end
