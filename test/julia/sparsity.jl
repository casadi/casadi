# sparsity.jl -- capability port of test/javascript/sparsity.js (+ deeper).
# Pure structure API: factories, queries, set algebra, triplet round-trip.
# Raw API via C.* (static methods dispatch on the C.Sparsity type tag).

const _S = C.Sparsity

# Collect (row,col) nonzeros into a Set of "r,c" strings (column-major order).
function _nzset(s)
    cols = C.get_col(s); rows = C.row(s); n = C.nnz(s)
    Set("$(rows[k]),$(cols[k])" for k in 1:n)
end

@testset "factories: dense/lower/diag/triplet/rowcol" begin
    d = C.dense(_S, 3, 4)
    @test (C.size1(d), C.size2(d)) == (3, 4)
    @test C.nnz(d) == 12
    @test C.numel(d) == 12
    @test C.is_dense(d)
    l = C.lower(_S, 4)                       # 4+3+2+1 = 10
    @test C.nnz(l) == 10
    @test (C.size1(l), C.size2(l)) == (4, 4)
    dg = C.diag(_S, 5)
    @test C.nnz(dg) == 5
    @test C.is_dense(dg) == false
    t = C.triplet(_S, 4, 5, [0, 0, 2, 3], [0, 1, 0, 1])
    @test C.nnz(t) == 4
    @test (C.size1(t), C.size2(t)) == (4, 5)
    # rowcol: column j gets the listed rows -> 2 rows in 2 cols = 4 nnz
    rc = C.rowcol(_S, [2, 0], [0, 2], 3, 3)
    @test C.colind(rc) == [0, 2, 2, 4]
end

@testset "add_nz builds same pattern as triplet" begin
    nz = [(0, 0), (0, 1), (2, 0), (2, 3), (2, 4), (3, 1)]
    a = C.Sparsity(4, 5)
    for (r, c) in nz; C.add_nz(a, r, c); end
    b = C.triplet(_S, 4, 5, [r for (r, _) in nz], [c for (_, c) in nz])
    @test C.nnz(a) == C.nnz(b) == 6
    @test _nzset(a) == _nzset(b)
end

@testset "triplet round-trip via row/get_col" begin
    rows = [0, 0, 2, 3]; cols = [0, 1, 0, 1]
    s = C.triplet(_S, 4, 5, rows, cols)
    # Column-major storage: get_col is sorted by column.
    @test C.get_col(s) == [0, 0, 1, 1]
    @test C.row(s) == [0, 2, 0, 3]
    @test C.colind(s) == [0, 2, 4, 4, 4, 4]    # cumulative per-column counts
    @test _nzset(s) == Set(["0,0", "0,1", "2,0", "3,1"])
    # get_triplet returns both index vectors via std::vector& output-args (argout).
    gr, gc = C.get_triplet(s)
    @test gr == C.row(s)
    @test gc == C.get_col(s)
end

@testset "queries: size/nnz/numel/is_dense/is_symmetric" begin
    s = C.dense(_S, 3, 4)
    @test (C.size1(s), C.size2(s), C.numel(s), C.nnz(s)) == (3, 4, 12, 12)
    @test C.is_dense(s)
    @test C.is_symmetric(C.dense(_S, 3, 3))
    @test C.is_symmetric(C.diag(_S, 4))
    @test C.is_symmetric(C.lower(_S, 3)) == false
    t = C.dense(_S, 3, 5)
    tT = C.T(t)
    @test (C.size1(tT), C.size2(tT)) == (5, 3)
end

@testset "set algebra: unite/intersect/is_subset" begin
    a = C.lower(_S, 3); b = C.dense(_S, 3, 3)
    @test C.is_subset(a, b)
    @test C.is_subset(b, a) == false
    @test C.nnz(C.unite(a, b)) == C.nnz(b)       # union(lower, dense) == dense
    @test C.nnz(C.intersect(a, b)) == C.nnz(a)   # intersect == lower
    @test C.is_subset(C.diag(_S, 3), C.lower(_S, 3))
    # explicit union/intersection counts on disjoint patterns
    nza = [(0, 0), (0, 1), (2, 0), (3, 1)]
    nzb = [(0, 2), (0, 0), (2, 2)]
    A = C.Sparsity(4, 5); for (r, c) in nza; C.add_nz(A, r, c); end
    B = C.Sparsity(4, 5); for (r, c) in nzb; C.add_nz(B, r, c); end
    @test C.nnz(C.unite(A, B)) == length(union(Set(nza), Set(nzb)))
    @test C.nnz(C.intersect(A, B)) == length(intersect(Set(nza), Set(nzb)))
end

@testset "kron / pattern_inverse / reshape preserve nnz" begin
    a = C.dense(_S, 2, 3); b = C.dense(_S, 4, 3)
    k = C.kron(a, b)
    @test (C.size1(k), C.size2(k)) == (C.size1(a) * C.size1(b), C.size2(a) * C.size2(b))
    @test C.nnz(k) == C.nnz(a) * C.nnz(b)
    inv = C.pattern_inverse(C.lower(_S, 4))      # 16 cells - 10 lower = 6
    @test C.nnz(inv) == 6
    r = C.reshape(C.dense(_S, 4, 5), 2, 10)
    @test (C.size1(r), C.size2(r)) == (2, 10)
    @test C.nnz(r) == 20
end

@testset "sparsity of an SX / MX expression" begin
    x = sym(SX, "x", 3)
    # elementwise square -> dense 3x1; jacobian wrt x is diagonal.
    e = C.times(x, x)
    @test C.nnz(C.sparsity(e)) == 3
    J = jacobian(e, x)
    sJ = C.sparsity(J)
    @test (C.size1(sJ), C.size2(sJ)) == (3, 3)
    @test C.nnz(sJ) == 3                          # diagonal Jacobian
    m = sym(MX, "m", 4)
    @test C.nnz(C.sparsity(m + m)) == 4
end

@testset "DM built on a custom Sparsity has matching shape" begin
    sp = C.triplet(_S, 4, 3, [0, 2, 3], [1, 2, 1])
    d = DM(sp, DM([5.0, 6.0, 7.0]))
    @test (C.size1(d), C.size2(d)) == (4, 3)
    @test C.nnz(d) == 3
    @test C.get_col(C.sparsity(d)) == C.get_col(sp)
end

@testset "== / isequal / hash compare patterns" begin
    @test C.dense(_S, 2) == C.dense(_S, 2)
    @test C.dense(_S, 2) != C.diag(_S, 2)
    @test isequal(C.lower(_S, 3), C.lower(_S, 3))
    @test hash(C.lower(_S, 3)) == hash(C.lower(_S, 3))
    @test length(Set([C.dense(_S, 2), C.dense(_S, 2), C.diag(_S, 2)])) == 2
    @test Sparsity() == Sparsity()
end

@testset "size / length / isempty" begin
    s = C.dense(_S, 2, 3)
    @test size(s) == (2, 3) && size(s, 1) == 2 && size(s, 2) == 3 && size(s, 3) == 1
    @test_throws BoundsError size(s, 0)
    @test_throws BoundsError size(DM([1, 2]), 0)
    @test length(s) == 6
    @test !isempty(s) && isempty(Sparsity(2, 0))
    @test isempty(SX()) && isempty(DM(Sparsity(0, 3))) && isempty(MX())
    @test !isempty(sym(SX, "x")) && !isempty(DM([1, 2]))
end

@testset "output arguments are returned" begin
    sp, mapping = C.triplet(_S, 3, 3, [2, 0, 1], [1, 0, 0], false)
    @test C.nnz(sp) == 3 && mapping == [1, 2, 0]
    n, index, offset = C.scc(C.banded(_S, 3, 1))
    @test n == 1 && index == [0, 1, 2] && offset == [0, 3]
    T, m = C.transpose(C.lower(_S, 2), false)
    @test C.is_equal(T, C.T(C.lower(_S, 2))) && length(m) == 3
end

@testset "in/out arguments are returned" begin
    @test C.get_nz(C.dense(_S, 2, 2), [0, 3]) == [0, 3]
    @test C.get_nz(C.diag(_S, 2), [1, 3]) == [-1, 1]
end

@testset "pattern algebra: + * .* ' and concatenation" begin
    d = C.diag(_S, 2); D = C.dense(_S, 2, 2)
    @test d + D == D
    @test d .* D == d && d .+ D == D
    @test C.diag(_S, 2) * C.dense(_S, 2, 1) == C.dense(_S, 2, 1)
    @test C.lower(_S, 2)' == C.upper(_S, 2) && transpose(C.lower(_S, 2)) == C.upper(_S, 2)
    @test size([d d]) == (2, 4) && size([d; d]) == (4, 2) && size([d D; D d]) == (4, 4)
    @test C.nnz([d D; D d]) == 12
    @test diagcat(d, d) == C.diag(_S, 4) && diagcat([d, d]) == C.diag(_S, 4)
end

@testset "coloring defaults" begin
    @test C.size2(C.uni_coloring(C.diag(_S, 3))) == 1
    @test C.size2(C.star_coloring(C.diag(_S, 5))) == 1
    @test C.size2(C.star_coloring2(C.diag(_S, 5))) == 1
    @test C.uni_coloring(C.diag(_S, 3), Sparsity()) == C.uni_coloring(C.diag(_S, 3))
end

@testset "null Sparsity displays" begin
    @test repr("text/plain", Sparsity()) == "Sparsity(NULL)"
    @test repr(Sparsity()) == "NULL"
end
