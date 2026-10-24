# CasADiNative.jl -- handwritten ergonomics layer over the SWIG-generated module.
# Prototype of what ships as the package's src/CasADi.jl; the generated
# casadi.jl stays untouched underneath.
module CasADiNative

include(joinpath(@__DIR__, "casadi.jl"))
using .casadi

const C = casadi
import LinearAlgebra: dot, cross, norm, det, tr
export SX, MX, DM, Sparsity, sym, nlpsol, jacobian, hessian
# P1-3: high-value builders that live as C.* but were unreachable via `using`
# (`reshape` is deliberately NOT exported: we extend Base.reshape, so a bare
#  `reshape` resolves to Base with our methods; exporting would shadow it. P2-1)
export blockcat, diagcat, sum1, sum2, repmat, dot, cross, norm, det, tr,
    mtimes, gradient, jtimes, substitute, simplify,
    qpsol, conic, integrator, rootfinder, CodeGenerator, codegen,
    external, generate, save
# blockcat/diagcat/dot/gradient get typed/convenience wrappers below; concatenation
# itself is [a b] / [a; b] (no horzcat/vertcat in the Julia API)
for _n in (:sum1, :sum2, :repmat, :reshape, :mtimes, :jtimes,
           :substitute, :simplify)
    @eval const $_n = C.$_n
end
export tangent, factory   # convenience defaults added below
# `Function` collides with Base.Function (abstract supertype): it is deliberately
# NOT exported. Under a plain `using`, bare `Function` resolves to Base.Function
# ("no constructors have been defined for Function"); access it qualified
# (`ca.Function`) or pull it in explicitly (`using LibCasADi: Function`). (P2-1)
const Function = C.Function
const CodeGenerator = C.CodeGenerator

const CasadiFunction = casadi.Function

# ---- 2-D / vector construction from native Julia arrays (P0-1) ----
# The SWIG ctor wants a vector-of-rows ([[float]]); a 2-D Julia Matrix reaches a
# flattening typemap, so repack rows explicitly to get a correct 2-D matrix.
_rows(m::AbstractMatrix{<:Real}) = [collect(Float64, @view m[i, :]) for i in 1:size(m, 1)]
# A 0-row matrix has no rows to carry the column count, so a vector-of-rows would
# collapse to 0×1; build the empty 0×n pattern directly instead. (P1-1)
_emptycols(::Type{DM}, n) = C.DM(C.Sparsity(0, n))
_emptycols(::Type{SX}, n) = C.SX(C.Sparsity(0, n))
_emptycols(::Type{MX}, n) = C.MX(C.DM(C.Sparsity(0, n)))
C.DM(m::AbstractMatrix{<:Real}) = size(m, 1) == 0 ? _emptycols(DM, size(m, 2)) : C.DM(_rows(m))
C.SX(m::AbstractMatrix{<:Real}) = size(m, 1) == 0 ? _emptycols(SX, size(m, 2)) : C.SX(_rows(m))
C.MX(m::AbstractMatrix{<:Real}) = size(m, 1) == 0 ? _emptycols(MX, size(m, 2)) : C.MX(C.DM(_rows(m)))
# AbstractVector{<:Real} -> a column (generated ctor already does this for DM/SX)
C.MX(v::AbstractVector{<:Real}) = C.MX(C.DM(collect(Float64, v)))

# ---- conversion (promotion by DM < SX < MX rank); generated ctors handle reals ----
# Julia numbers and arrays promote exactly like DM (`_Num`).
const _Num = Union{Real,DM,AbstractArray{<:Real}}
const _Proxy = Union{(t for (_, t) in casadi._proxy_types)...}
Base.convert(::Type{SX}, v::_Num) = C.SX(v isa Real ? Float64(v) : v)
Base.convert(::Type{MX}, v::_Num) = C.MX(v isa Real ? Float64(v) : v)
Base.convert(::Type{DM}, v::Union{Real,AbstractArray{<:Real}}) = C.DM(v isa Real ? Float64(v) : v)

# ---- operators ----
for (op, fn) in ((:+, :plus), (:-, :minus), (:*, :mtimes))
    @eval begin
        Base.$op(a::SX, b::SX) = C.$fn(a, b)
        Base.$op(a::MX, b::MX) = C.$fn(a, b)
        Base.$op(a::DM, b::DM) = C.$fn(a, b)
        # promotion: numbers, Julia arrays and DM promote to the symbolic side
        Base.$op(a::SX, b::_Num) = C.$fn(a, convert(SX, b))
        Base.$op(a::_Num, b::SX) = C.$fn(convert(SX, a), b)
        Base.$op(a::MX, b::_Num) = C.$fn(a, convert(MX, b))
        Base.$op(a::_Num, b::MX) = C.$fn(convert(MX, a), b)
        Base.$op(a::DM, b::Union{Real,AbstractArray{<:Real}}) = C.$fn(a, convert(DM, b))
        Base.$op(a::Union{Real,AbstractArray{<:Real}}, b::DM) = C.$fn(convert(DM, a), b)
    end
end
Base.:-(a::Union{SX,MX,DM}) = C.times(a, typeof(a) === DM ? C.DM(-1.0) : convert(typeof(a), -1))

# Where Julia gives an array operation a meaning (A/B, A^n, matrix functions, max of
# arrays), a casadi matrix does the same or throws; never an elementwise surprise.
# Elementwise is spelled with dots (./ .^ sin.), as in Julia.
_scal(x::Union{SX,MX,DM}) = C.numel(x) == 1
_scal(x::Real) = true
_scal(x) = false
_shape(x) = "$(C.size1(x))×$(C.size2(x))"
Base.:/(a::Union{SX,MX,DM}, b::Union{SX,MX,DM,Real,AbstractArray{<:Real}}) = _div(a, b)
Base.:/(a::Union{Real,AbstractArray{<:Real}}, b::Union{SX,MX,DM}) = _div(a, b)
function _div(a, b)
    T = _catpromote((a, b))
    _scal(b) && return C.rdivide(_ascat(T, a), _ascat(T, b))
    _scal(a) && Base.throw(ArgumentError("x / A with a matrix A has no meaning in Julia; use x ./ A"))
    C.mrdivide(_ascat(T, a), _ascat(T, b))
end
Base.:^(a::Union{SX,MX,DM}, b::Union{SX,MX,DM,Real}) = _pow(a, b)
Base.:^(a::Real, b::Union{SX,MX,DM}) = _pow(a, b)
function _pow(a, b)
    T = _catpromote((a, b))
    _scal(a) && _scal(b) && return C.power(_ascat(T, a), _ascat(T, b))
    p = b isa DM && _scal(b) ? Float64(b) : b
    if a isa _Cat && p isa Real && Base.isinteger(p)                  # A^n: casadi's mpower
        C.size1(a) == C.size2(a) || Base.throw(DimensionMismatch("matrix is not square: dimensions are $(_shape(a))"))
        n = Int(p)
        return n >= 0 ? C.mpower(a, n) : C.mpower(C.inv(a), -n)
    end
    Base.throw(ArgumentError("A^p of a matrix needs an integer p (matrix power); use .^ for " *
        "elementwise, or a Julia Matrix(A) for Julia's A^p and p^A"))
end

# ---- broadcasting (P0-2): casadi ops are already elementwise on matrices, so
# every `.op` / `f.(x)` short-circuits to the matrix-wide elementwise casadi op ----
struct CasadiStyle <: Base.Broadcast.BroadcastStyle end
Base.broadcastable(x::Union{SX,MX,DM}) = x   # treat as a scalar leaf, not iterable
Base.BroadcastStyle(::Type{<:Union{SX,MX,DM}}) = CasadiStyle()
Base.BroadcastStyle(s::CasadiStyle, ::Base.Broadcast.BroadcastStyle) = s
Base.BroadcastStyle(s::CasadiStyle, ::CasadiStyle) = s
Base.BroadcastStyle(s::CasadiStyle, ::Base.Broadcast.Unknown) = s  # disambiguate

# rank of the casadi operands present, to coerce scalars/DM/arrays to that type
_rank(::Type{DM}) = 1; _rank(::Type{SX}) = 2; _rank(::Type{MX}) = 3
_castype() = DM
_castype(a, rest...) = _castype(rest...)
_castype(a::Union{SX,MX,DM}, rest...) =
    (T = _castype(rest...); _rank(typeof(a)) >= _rank(T) ? typeof(a) : T)
_coerce(::Type{T}, a::Union{SX,MX,DM}) where {T} = convert(T, a)
_coerce(::Type{T}, a::Real) where {T} = convert(T, a)
_coerce(::Type{T}, a::AbstractArray{<:Real}) where {T} = T(a)   # native array -> casadi matrix
_coerce(::Type{T}, a) where {T} = a   # ops, Refs, Val: pass through untouched

# the applied op: unwrap literal_pow (x.^2) into a plain casadi power
_apply(f, args...) = f(args...)
_apply(::typeof(Base.literal_pow), ::Base.RefValue{typeof(^)}, x, ::Base.RefValue{Val{n}}) where {n} =
    C.power(x, n)
# f.(args) with a known f: one whole-matrix casadi op (elementwise); `nothing`: unknown f
_bcop(f) = nothing
_bcop(::typeof(+)) = +
_bcop(::typeof(-)) = -
_bcop(::typeof(*)) = C.times
_bcop(::typeof(/)) = C.rdivide
_bcop(::typeof(^)) = C.power
_bcop(::typeof(<)) = C.lt
_bcop(::typeof(<=)) = C.le
_bcop(::typeof(>)) = C.gt
_bcop(::typeof(>=)) = C.ge
_bcop(::typeof(==)) = C.eq
_bcop(::typeof(!=)) = C.ne
_bcop(::typeof(Base.max)) = C.fmax
_bcop(::typeof(Base.min)) = C.fmin
_bcop(::typeof(Base.hypot)) = C.hypot
_bcop(::typeof(Base.atan)) = (a...) -> Base.length(a) == 2 ? C.atan2(a...) : C.atan(a...)
_bcop(::typeof(Base.ifelse)) = C.if_else
_bcop(::typeof(&)) = C.logic_and
_bcop(::typeof(|)) = C.logic_or
_bcop(::typeof(!)) = C.logic_not
for fn in (:sin, :cos, :tan, :exp, :log, :sqrt, :abs, :tanh, :sinh, :cosh, :asin, :acos,
           :floor, :ceil, :sign, :log10, :log1p, :expm1, :asinh, :acosh, :atanh)
    @eval _bcop(::typeof(Base.$fn)) = C.$fn
end

# materialize by recursively walking the (unflattened) tree, so nested `.* ` inside a
# fused chain still maps to elementwise times rather than Base.:* (mtimes).
_bceval(x) = x
function _bceval(bc::Base.Broadcast.Broadcasted{CasadiStyle})
    args = Base.map(_bceval, bc.args)
    bc.f === Base.literal_pow && return _apply(bc.f, args...)
    # Float64.(dm) etc.: a numeric conversion broadcast over a single casadi matrix
    # escapes to a Julia array of converted scalars (not a casadi op).
    if Base.length(args) == 1 && bc.f isa Type{<:Number} && args[1] isa _Cat
        return bc.f.(Base.collect(args[1]))
    end
    op = _bcop(bc.f)
    op === nothing && return _bcelem(bc.f, args)
    T = _castype(args...)
    op(map(a -> _coerce(T, a), args)...)
end
# any other f is applied to each element, with Julia's broadcast shape rules
function _bcelem(f, args)
    r = Base.broadcast(f, Base.map(a -> a isa _Cat ? Base.collect(a) : a, args)...)
    r isa AbstractArray || return r
    T = _catpromote(Base.vec(r))
    _restore(T, C.reshape(C._vertcat(_cvec(Base.vec(r))), Base.size(r, 1), Base.size(r, 2)))
end
Base.copy(bc::Base.Broadcast.Broadcasted{CasadiStyle}) = _bceval(bc)

# ---- concatenation: [a b], [a;b], [a b;c d] (P0-3) ----
const _Cat = Union{SX,MX,DM}
_catpromote(xs) = Base.any(x -> x isa MX, xs) ? MX : Base.any(x -> x isa SX, xs) ? SX : DM
_ascat(::Type{T}, x::T) where {T<:_Cat} = x                 # same type: identity
_ascat(::Type{T}, x::_Cat) where {T} = convert(T, x)        # promote DM->SX/MX etc.
_ascat(T, x::Number) = convert(T, x)
_ascat(T, x::AbstractArray) = T(x)
_cvec(xs) = (T = _catpromote(xs); [_ascat(T, x) for x in xs])
# generated _vertcat/_horzcat always return MX; restore the DM type for numeric cats
_restore(::Type{DM}, r) = C.evalf(r)
_restore(::Type{T}, r) where {T} = r
_cat(builder, xs) = _restore(_catpromote(xs), builder(_cvec(xs)))
# One vararg method per bracket form, over casadi values, numbers and Julia arrays:
# a casadi operand may sit anywhere ([0; 0; v], [1 2 x], [A x; y' 0], vcat(w...)).
# Calls without one keep reaching Base, whose Number/AbstractArray methods are
# strictly more specific (this is why the union says Number, not Real); the
# fallbacks below only guard odd mixes Base has no specific method for.  (P1-2)
const _CatArg = Union{_Cat,Number,AbstractArray}
_hascat(xs) = Base.any(x -> x isa _Cat, xs)
# The pieces become a Vector{Any} at once: passing the tuple on would compile
# every helper anew for each length of a splatted vcat(w...).
Base.vcat(@nospecialize(xs::_CatArg...)) = _vcat(Any[xs...])
Base.hcat(@nospecialize(xs::_CatArg...)) = _hcat(Any[xs...])
Base.hvcat(rows::Tuple{Vararg{Int}}, @nospecialize(xs::_CatArg...)) = _hvcat_any(rows, Any[xs...])
_vcat(v::Vector{Any}) = _hascat(v) ? _cat(C._vertcat, v) : Base.cat(v...; dims = 1)
_hcat(v::Vector{Any}) = _hascat(v) ? _cat(C._horzcat, v) : Base.cat(v...; dims = 2)
_hvcat_any(rows, v::Vector{Any}) = _hascat(v) ? _hvcat(rows, v) :
    Base.invoke(Base.hvcat, Tuple{Tuple{Vararg{Int}},Vararg{Any}}, rows, v...)
# reduce(vcat, w) over a vector of casadi values: one vertcat, not a pairwise fold
Base.reduce(::typeof(vcat), v::AbstractVector{<:_Cat}) = _cat(C._vertcat, v)
Base.reduce(::typeof(hcat), v::AbstractVector{<:_Cat}) = _cat(C._horzcat, v)
# [a b; c d] -> rows of horzcat, then vertcat
function _hvcat(rows, xs)
    T = _catpromote(xs)
    cv = [_ascat(T, x) for x in xs]
    blocks = Any[]
    k = 0
    for ncol in rows
        Base.push!(blocks, C._horzcat(cv[(k + 1):(k + ncol)]))
        k += ncol
    end
    _restore(T, C._vertcat(blocks))
end

# ---- typed function-form cat/dot (R3): preserve DM/SX, MX otherwise ----
# The generated builders promote everything to MX; restore DM (evalf-collapse)
# when all args are DM/Real and SX when all are SX/Real.
_catargs(xs::Tuple{AbstractVector}) = Base.collect(xs[1])  # single-vector call form
_catargs(xs::Tuple) = Base.collect(xs)                     # varargs call form
function _typedcat(builder, args...)
    xs = _catargs(args)
    Base.isempty(xs) && return builder(xs)
    T = _catpromote(xs)
    _restore(T, builder([_ascat(T, x) for x in xs]))
end
# Varargs may mix casadi values, numbers and Julia arrays; a single Julia vector
# is the list form (its elements are the pieces), as in Python.
diagcat(args::_CatArg...) = _typedcat(C.diagcat, args...)
diagcat(v::AbstractVector) = _typedcat(C.diagcat, v)
# blockcat takes a 2-D matrix-of-blocks; lower it to vertcat-of-horzcats
function blockcat(rows::AbstractVector{<:AbstractVector})
    xs = Base.collect(Base.Iterators.flatten(rows))
    T = _catpromote(xs)
    _restore(T, C._vertcat([C._horzcat([_ascat(T, x) for x in r]) for r in rows]))
end
blockcat(A::_Cat, B::_Cat, C_::_Cat, D::_Cat) =
    blockcat([[A, B], [C_, D]])
# dot: promote both sides to a common type, then the matching generated method
function dot(a::_Cat, b::_CatArg)
    T = _catpromote((a, b))
    _restore(T, C.dot(_ascat(T, a), _ascat(T, b)))
end
dot(a::Union{Real,AbstractArray{<:Real}}, b::_Cat) = dot(b, a)
cross(a::T, b::T, dim::Integer...) where {T<:_Cat} = C.cross(a, b, dim...)
# MATLAB-style backslash solve `A \ b` (linear solve), distinct from elementwise
function Base.:\(a::_Cat, b::_CatArg)
    T = _catpromote((a, b))
    _restore(T, C.solve(_ascat(T, a), _ascat(T, b)))
end
Base.:\(a::AbstractMatrix{<:Real}, b::_Cat) = convert(typeof(b), a) \ b

# ---- Base math extensions (no shadowing: extend Base methods) ----
# Julia's sin(A), exp(A), ... of a square matrix are matrix functions, which casadi does
# not compute: a non-scalar throws (use sin.(A)).  Scalars (1×1) are plain casadi ops.
for fn in (:sin, :cos, :tan, :exp, :log, :sqrt, :tanh, :sinh, :cosh,
           :asin, :acos, :atan, :asinh, :acosh, :atanh)
    @eval Base.$fn(x::Union{SX,MX,DM}) = _scal(x) ? C.$fn(x) : Base.throw(ArgumentError(
        "$($(QuoteNode(fn)))(A) of a $(_shape(x)) $(typeof(x)): in Julia a matrix function of a " *
        "square A, undefined otherwise; use $($(QuoteNode(fn))).(A) for elementwise, or a Julia Matrix(A)"))
end
# no array meaning in Julia for these: elementwise, as casadi
for fn in (:abs, :floor, :ceil, :sign, :log10, :log1p, :expm1)
    @eval Base.$fn(x::Union{SX,MX,DM}) = C.$fn(x)
end

# Reductions and binary math with Julia's own names and semantics (sum, diff,
# cumsum, maximum, max, norm, inv, det, tr, atan(y, x), hypot), so Julia code
# reads like Matlab's.  Results keep the casadi type.  `Base.` is spelled out:
# inside this module the generated casadi exports clash with Base's names.
_dims(x, dims, what) = dims !== nothing ? Int(dims) :
    C.size2(x) == 1 ? 1 : C.size1(x) == 1 ? 2 :
    Base.throw(ArgumentError("$what of a matrix needs dims, e.g. $what(A; dims=1)"))
function Base.sum(x::_Cat; dims = :)
    dims === Colon() && return C.sum2(C.sum1(x))       # all elements, as for Julia arrays
    dims == 1 ? C.sum1(x) : dims == 2 ? C.sum2(x) : x
end
Base.diff(x::_Cat; dims::Union{Integer,Nothing} = nothing) = C.diff(x, 1, _dims(x, dims, "diff") - 1)
Base.cumsum(x::_Cat; dims::Union{Integer,Nothing} = nothing) = C.cumsum(x, _dims(x, dims, "cumsum") - 1)
Base.maximum(x::_Cat) = C.mmax(x)
Base.minimum(x::_Cat) = C.mmin(x)
Base.inv(x::_Cat) = C.inv(x)
det(x::_Cat) = C.det(x)
tr(x::_Cat) = C.trace(x)
# Julia's norm: a vector's p-norm, or entrywise for a matrix (Frobenius for p = 2)
function norm(x::_Cat, p::Real = 2)
    isvec = C.size1(x) == 1 || C.size2(x) == 1
    p == 2 ? (isvec ? C.norm_2(x) : C.norm_fro(x)) :
    p == 1 ? (isvec ? C.norm_1(x) : Base.sum(Base.abs(x))) :
    p == Inf ? (isvec ? C.norm_inf(x) : Base.maximum(Base.abs(x))) :
    Base.throw(ArgumentError("norm(x, p) of a casadi value supports p = 1, 2 or Inf"))
end
# Elementwise two-argument functions; the generated overloads pick MX/SX/DM from
# the argument types, so numbers and Julia arrays mix in as anywhere else.
const _Native = Union{Real,AbstractArray{<:Real}}
# atan(y, x) and hypot have no array meaning in Julia: elementwise.  max/min of arrays
# do (a lexicographic pick of vectors), so they take scalars only.
for (fn, cfn) in ((:atan, :atan2), (:hypot, :hypot))
    @eval begin
        Base.$fn(a::_Cat, b::_Cat) = C.$cfn(a, b)
        Base.$fn(a::_Cat, b::_Native) = C.$cfn(a, b)
        Base.$fn(a::_Native, b::_Cat) = C.$cfn(a, b)
    end
end
for (fn, cfn) in ((:max, :fmax), (:min, :fmin))
    @eval begin
        Base.$fn(a::_Cat, b::Union{_Cat,_Native}) = _maxmin(Base.$fn, C.$cfn, a, b)
        Base.$fn(a::_Native, b::_Cat) = _maxmin(Base.$fn, C.$cfn, a, b)
    end
end
function _maxmin(f, cf, a, b)
    _scal(a) && _scal(b) && return cf(a, b)
    Base.throw(ArgumentError("$f(A, B) of arrays is a lexicographic pick in Julia; use $f.(A, B) for elementwise"))
end

# & | ! of scalars, as for Julia Bools; elementwise is .& .| .!
for (op, cfn) in ((:&, :logic_and), (:|, :logic_or))
    @eval begin
        Base.$op(a::_Cat, b::Union{_Cat,Real}) = _logic($(QuoteNode(op)), C.$cfn, a, b)
        Base.$op(a::Real, b::_Cat) = _logic($(QuoteNode(op)), C.$cfn, a, b)
    end
end
function _logic(op, cf, a, b)
    _scal(a) && _scal(b) || Base.throw(ArgumentError("a $op b of arrays has no meaning in Julia; use a .$op b"))
    T = _catpromote((a, b))
    cf(_ascat(T, a), _ascat(T, b))
end
Base.:!(a::_Cat) = _scal(a) ? C.logic_not(a) :
    Base.throw(ArgumentError("!A of a $(_shape(a)) $(typeof(a)) has no meaning in Julia; use .!A"))

# Integer-style division of scalars with Julia's rounding (rem/% truncates, mod floors);
# elementwise is rem.(A, B) etc., as for Julia arrays.
function _scalar2(f, a, b)
    _scal(a) && _scal(b) || Base.throw(ArgumentError("$f of arrays has no meaning in Julia; use $f.(A, B)"))
    T = _catpromote((a, b))
    _ascat(T, a), _ascat(T, b)
end
for f in (:rem, :mod, :fld, :cld, :div, :copysign)
    @eval begin
        Base.$f(a::_Cat, b::Union{_Cat,Real}) = _int_div(Base.$f, _scalar2($(QuoteNode(f)), a, b)...)
        Base.$f(a::Real, b::_Cat) = _int_div(Base.$f, _scalar2($(QuoteNode(f)), a, b)...)
    end
end
_int_div(::typeof(Base.rem), a, b) = C.fmod(a, b)
_int_div(::typeof(Base.mod), a, b) = C.fmod(C.fmod(a, b) + b, b)
_int_div(::typeof(Base.fld), a, b) = C.floor(a / b)
_int_div(::typeof(Base.cld), a, b) = C.ceil(a / b)
_int_div(::typeof(Base.div), a, b) = (q = a / b; C.sign(q) * C.floor(C.abs(q)))
_int_div(::typeof(Base.copysign), a, b) = C.copysign(a, b)
_bcop(::typeof(Base.rem)) = C.fmod
_bcop(::typeof(Base.copysign)) = C.copysign

# ifelse(c, a, b) on a casadi condition is casadi's if_else
Base.ifelse(c::_Cat, a, b) = C.if_else(c, a, b)

# ---- transpose: MATLAB-style `A'` and `transpose(A)` ----
Base.transpose(x::Union{SX,MX,DM}) = C.transpose(x)
Base.adjoint(x::Union{SX,MX,DM}) = C.transpose(x)   # real/symbolic: adjoint == transpose

# ---- show (compact 2-arg: used by string interpolation; KEEP) ----
Base.show(io::IO, x::SX) = print(io, C.str(x))
Base.show(io::IO, x::MX) = print(io, C.str(x))
Base.show(io::IO, x::DM) = print(io, C.str(x))
Base.show(io::IO, f::CasadiFunction) = print(io, C.str(f))
Base.show(io::IO, s::Sparsity) = print(io, C.str(s))
Base.show(io::IO, o::C.Opti) = print(io, C.str(o))

# ---- rich REPL display (3-arg MIME: header + body) ----
const _MIME = MIME"text/plain"

# Format a column-major Float64 matrix as right-justified aligned rows.
function _fmt_matrix(io::IO, M::AbstractMatrix{Float64})
    cells = map(x -> isinteger(x) ? string(Int(x)) : string(round(x; digits=4)), M)
    w = Base.maximum(Base.length, cells; init=0)
    for i in 1:size(cells, 1)
        print(io, " ")
        for j in 1:size(cells, 2)
            print(io, Base.lpad(cells[i, j], w), j == size(cells, 2) ? "" : "  ")
        end
        i == size(cells, 1) || println(io)
    end
end

function Base.show(io::IO, ::_MIME, d::DM)
    m, n = C.size1(d), C.size2(d)
    if m == 1 && n == 1
        return print(io, C.str(d))   # scalar: just the value
    end
    nz = C.nnz(d)
    print(io, m, "×", n, " DM", C.is_dense(d) ? "" : " ($(nz) nnz)", ":\n")
    _fmt_matrix(io, Base.Matrix{Float64}(d))
end

function Base.show(io::IO, ::_MIME, x::SX)
    print(io, C.size1(x), "×", C.size2(x), " SX:\n", C.str(x))
end

function Base.show(io::IO, ::_MIME, x::MX)
    desc = C.is_symbolic(x) ? "symbolic" : C.is_constant(x) ? "constant" : C.class_name(x)
    print(io, C.size1(x), "×", C.size2(x), " MX (", desc, "):\n", C.str(x))
end

function Base.show(io::IO, ::_MIME, s::Sparsity)
    m, n = C.size1(s), C.size2(s)
    print(io, m, "×", n, " Sparsity (", C.nnz(s), " nnz)")
    (m == 0 || n == 0) && return
    if m <= 20 && n <= 20
        grid = Base.fill('.', m, n)         # spy pattern: '*' stored, '.' zero
        rows, cols = C.get_triplet(s)
        for k in eachindex(rows)
            grid[rows[k] + 1, cols[k] + 1] = '*'
        end
        print(io, ":\n")
        for i in 1:m
            print(io, " ", Base.String(grid[i, :]))
            i == m || println(io)
        end
    end
end

# One in/out slot rendered as `name[rows]` (drop trailing `x1` for column vectors).
function _slot(name, sz)
    r, c = sz[1], sz[2]
    c == 1 ? "$(name)[$(r)]" : "$(name)[$(r)×$(c)]"
end
function Base.show(io::IO, ::_MIME, f::CasadiFunction)
    ins = join((_slot(C.name_in(f, i - 1), C.size_in(f, i - 1)) for i in 1:C.n_in(f)), ", ")
    outs = join((_slot(C.name_out(f, i - 1), C.size_out(f, i - 1)) for i in 1:C.n_out(f)), ", ")
    print(io, "Function ", C.name(f), "(", ins, ") -> (", outs, ")")
end

function Base.show(io::IO, ::_MIME, o::C.Opti)
    print(io, "Opti(", C.nx(o), " vars, ", C.ng(o), " constraints, ", C.np(o), " params)")
end
function Base.show(io::IO, ::_MIME, s::C.OptiSol)
    print(io, "OptiSol of ", Base.repr(_MIME(), C.opti(s)))
end

# ---- size / conversion ----
Base.size(x::Union{SX,MX,DM,Sparsity}) = (C.size1(x), C.size2(x))
Base.size(x::Union{SX,MX,DM,Sparsity}, d::Integer) =
    d < 1 ? Base.throw(BoundsError(Base.size(x), d)) : d <= 2 ? Base.size(x)[d] : 1
Base.length(x::Union{SX,MX,DM,Sparsity}) = C.numel(x)
Base.isempty(x::Union{SX,MX,DM,Sparsity}) = C.numel(x) == 0
Base.Float64(d::DM) = C.numel(d) == 1 ? Base.Matrix{Float64}(d)[1] :
    Base.throw(DimensionMismatch("Float64 of a $(_shape(d)) DM; only a 1×1 DM converts"))
# scalar conversions follow the Float64 value, as for Julia numbers (Bool(2.0) throws)
for T in (:Bool, :Int32, :Int64, :Float32)
    @eval Base.$T(x::_Cat) = $T(Float64(x))
end
for f in (:round, :trunc, :floor, :ceil)
    @eval Base.$f(::Type{T}, x::_Cat) where {T<:Integer} = Base.$f(T, Float64(x))
end
# P1-5: dense, shape-correct conversions (densify -> column-major nonzeros)
Base.Matrix{Float64}(d::DM) =
    Base.reshape(C.nonzeros(C.densify(d)), (C.size1(d), C.size2(d)))
Base.Matrix(d::DM) = Base.Matrix{Float64}(d)
Base.Array(d::DM) = Base.Matrix{Float64}(d)
Base.collect(d::DM) = Base.Matrix{Float64}(d)
Base.eltype(::Type{DM}) = Float64
# a row or column -> dense Vector, as Vector(x) for SX/MX
Base.Vector{T}(d::DM) where {T} = C.size1(d) == 1 || C.size2(d) == 1 ? Base.Vector{T}(Base.vec(Base.Matrix{Float64}(d))) :
    Base.throw(DimensionMismatch("Vector of a $(_shape(d)) matrix; use vec(Matrix(x))"))
Base.Vector(d::DM) = Base.Vector{Float64}(d)

# ---- zero/one: scalar for the type (so zeros(SX,m,n) works), shaped for a value ----
for T in (:SX, :MX, :DM)
    @eval Base.zero(::Type{$T}) = C.$T(0.0)
    @eval Base.one(::Type{$T})  = C.$T(1.0)
    @eval Base.zero(x::$T) = C.zeros($T, C.size1(x), C.size2(x))
    @eval Base.one(x::$T) = C.size1(x) == C.size2(x) ? C.eye($T, C.size1(x)) :
        Base.throw(DimensionMismatch("multiplicative identity defined only for square matrices"))
end

# ---- Float64 of a CONSTANT scalar SX/MX (DM already has Float64 above) ----
function _const_scalar_float(x::Union{SX,MX})
    (C.is_constant(x) && C.numel(x) == 1) ||
        Base.error("Float64 requires a constant scalar; got $(C.size1(x))×$(C.size2(x)) non-constant")
    Float64(C.evalf(x))   # evalf -> DM -> Float64
end
Base.Float64(x::SX) = _const_scalar_float(x)
Base.Float64(x::MX) = _const_scalar_float(x)

# ---- collect / Array: dense Matrix of scalar SX/MX elements (NOT iterable, so
# broadcasting's scalar-leaf treatment via broadcastable(x)=x is preserved) ----
_to_matrix(x::Union{SX,MX}) =
    [x[i, j] for i in 1:C.size1(x), j in 1:C.size2(x)]
Base.collect(x::SX) = _to_matrix(x)
Base.collect(x::MX) = _to_matrix(x)
Base.Array(x::SX) = _to_matrix(x)
Base.Array(x::MX) = _to_matrix(x)
# a row or column -> Vector of scalar elements
Base.Vector(x::Union{SX,MX}) = C.size1(x) == 1 || C.size2(x) == 1 ? Base.vec(_to_matrix(x)) :
    Base.throw(DimensionMismatch("Vector of a $(C.size1(x))×$(C.size2(x)) matrix; use vec(x) or collect(x)"))

# ---- Base bridges to casadi free functions (typed: preserve DM/SX/MX) ----
Base.vec(x::Union{SX,MX,DM}) = C.vec(x)
Base.kron(a::T, b::T) where {T<:Union{SX,MX,DM}} = C.kron(a, b)
Base.kron(a::_Cat, b::_Cat) = (T = _catpromote((a, b)); C.kron(_ascat(T, a), _ascat(T, b)))
Base.repeat(x::Union{SX,MX,DM}, m::Integer, n::Integer) = C.repmat(x, Int(m), Int(n))
Base.repeat(x::Union{SX,MX,DM}, m::Integer) = C.repmat(x, Int(m), 1)
# Base.reshape on our types -> C.reshape (coexists with the `reshape` const export,
# which is C.reshape itself, so no ambiguity; Base method wins on Union{SX,MX,DM}).
# one dimension may be `:` (inferred), as for Julia arrays
const _Dim = Union{Integer,Colon}
function Base.reshape(x::_Cat, nrow::_Dim, ncol::_Dim)
    n = C.numel(x)
    _inf(d, other) = other isa Colon ? Base.throw(DimensionMismatch("only one dimension can be inferred")) :
        (other != 0 && n % other == 0 ? n ÷ other : Base.throw(DimensionMismatch("cannot reshape $(n) elements into (:, $other)")))
    r = nrow isa Colon ? _inf(nrow, ncol) : Int(nrow)
    c = ncol isa Colon ? _inf(ncol, nrow) : Int(ncol)
    C.reshape(x, r, c)
end
Base.reshape(x::_Cat, dims::Tuple{_Dim,_Dim}) = Base.reshape(x, dims[1], dims[2])
Base.reshape(x::_Cat, sp::Sparsity) = C.reshape(x, sp)

# ---- functor: f(args...) like Python ----
# Julia arrays pass through as-is: the wrapper converts them like DM (`[]` -> empty).
function (f::CasadiFunction)(args::Vararg{Union{Real,DM,SX,MX,AbstractArray}})
    ins = Any[a isa Real ? C.DM(Float64(a)) : a for a in args]
    r = C.call(f, ins)
    Base.length(r) == 1 ? r[1] : Base.Tuple(r)
end
(f::CasadiFunction)(d::Dict{String,<:Any}) = C.call(f, Dict{String,Any}(d))
function (f::CasadiFunction)(; kwargs...)
    # f() (no args) is the positional zero-input call, not an empty dict-call (P2-2)
    if Base.isempty(kwargs)
        r = C.call(f, Any[])
        return Base.length(r) == 1 ? r[1] : Base.Tuple(r)
    end
    C.call(f, Dict{String,Any}(Base.string(k) => (v isa Real ? C.DM(Float64(v)) : v)
                               for (k, v) in kwargs))
end

# `sym` (and ca.MX.sym) come from the generated static method-style; no wrapper needed.
_optdict(opts) = Dict{String,Any}(Base.string(k) => v for (k, v) in opts)
# jacobian/hessian/gradient/tangent need no wrappers: the SWIG backend now emits
# their trailing `opts::Any=Dict()` default, so `jacobian(ex,arg)` comes from
# `using .casadi` directly. (Generic opts-defaulting moved into julia.cxx.)
# Function ctor needs no wrapper: backend defaults opts and the named-IO form is the
# generated positional `Function(name, ex_in, ex_out, name_in, name_out[, opts])`.
# factory needs no wrapper: the backend defaults aux/opts, and the generic Dict
# typemap already marshals any Julia Dict -> map<string,vector<string>>. Reachable
# as the free `factory` (end-loop binding) and method-style `f.factory(...)`.
# nlpsol/qpsol/conic/rootfinder/integrator/external/interpolant/generate/save/
# generate_dependencies need no wrappers: the backend defaults their trailing opts
# and the generic Dict typemap marshals a Julia problem dict directly. Options are
# passed casadi-style as a positional Dict (NOT Julia kwargs). They reach the user
# via `using .casadi` (exported) or the end-loop binding.
const codegen = C.CodeGenerator   # short alias; CodeGenerator(name[, opts]) opts-defaulted
export interpolant, generate_dependencies

# ---- comparisons: elementwise, in the casadi type (SX/MX constraints, DM 0/1 values) ----
const _Real = Union{Real,AbstractArray{<:Real}}
for (op, fn) in ((:<, :lt), (:<=, :le), (:>=, :ge), (:(==), :eq), (:(!=), :ne))
    @eval begin
        Base.$op(a::T, b::T) where {T<:Union{SX,MX}} = C.$fn(a, b)
        Base.$op(a::T, b::_Num) where {T<:Union{SX,MX}} = C.$fn(a, convert(T, b))
        Base.$op(a::_Num, b::T) where {T<:Union{SX,MX}} = C.$fn(convert(T, a), b)
        Base.$op(a::DM, b::DM) = C.$fn(a, b)
        Base.$op(a::DM, b::_Real) = C.$fn(a, convert(DM, b))
        Base.$op(a::_Real, b::DM) = C.$fn(convert(DM, a), b)
    end
end
# isequal/hash/iszero stay Bool (Dict keys): structural identity, casadi's is_equal at depth 0
Base.isequal(a::T, b::T) where {T<:Union{SX,MX}} = C.is_equal(a, b)
Base.hash(x::MX, h::UInt) = Base.hash(C.__hash__(x), h)
Base.hash(x::SX, h::UInt) = C.numel(x) == 1 && C.nnz(x) == 1 ?
    Base.hash(C.element_hash(x), h) : Base.hash(Base.size(x), h)
Base.iszero(x::Union{SX,MX}) = C.is_zero(x)
# a DM is a value: equal when shape and entries are
Base.isequal(a::DM, b::DM) = Base.isequal(Base.Matrix{Float64}(a), Base.Matrix{Float64}(b))
Base.hash(x::DM, h::UInt) = Base.hash(Base.Matrix{Float64}(x), h)
Base.iszero(x::DM) = Base.iszero(Base.Matrix{Float64}(x))

# a Sparsity is a pattern: equal when the patterns are
Base.:(==)(a::Sparsity, b::Sparsity) =
    C.is_null(a) || C.is_null(b) ? C.is_null(a) && C.is_null(b) : C.is_equal(a, b)
Base.isequal(a::Sparsity, b::Sparsity) = a == b
Base.hash(s::Sparsity, h::UInt) = C.is_null(s) ? Base.hash(nothing, h) :
    Base.hash((C.size1(s), C.size2(s), C.colind(s), C.row(s)), h)

# Pattern algebra with the meaning of the matrix operations: A + B, A * B (product),
# A .* B (elementwise), A', [A B], [A; B], diagcat
Base.:+(a::Sparsity, b::Sparsity) = C.unite(a, b)
Base.:*(a::Sparsity, b::Sparsity) = C.mtimes(a, b)
Base.broadcasted(::typeof(*), a::Sparsity, b::Sparsity) = C.intersect(a, b)
Base.broadcasted(::typeof(+), a::Sparsity, b::Sparsity) = C.unite(a, b)
Base.transpose(s::Sparsity) = C.T(s)
Base.adjoint(s::Sparsity) = C.T(s)
Base.hcat(x::Sparsity, xs::Sparsity...) = C._horzcat(Sparsity[x, xs...])
Base.vcat(x::Sparsity, xs::Sparsity...) = C._vertcat(Sparsity[x, xs...])
function Base.hvcat(rows::Tuple{Vararg{Int}}, x::Sparsity, xs::Sparsity...)
    v = Sparsity[x, xs...]; k = 0
    C._vertcat(Sparsity[(k += n; C._horzcat(v[k-n+1:k])) for n in rows])
end
diagcat(x::Sparsity, xs::Sparsity...) = C.diagcat(Sparsity[x, xs...])
diagcat(v::AbstractVector{Sparsity}) = C.diagcat(v)

# ---- Opti conveniences ----
# variable/parameter/subject_to need no wrappers: they reach the generated Opti
# methods directly (subject_to's opts is backend-defaulted) but stay EXPORTED so
# `using`-style bare calls resolve. Only the value-coercing ones keep a wrapper.
export Opti, variable, parameter, subject_to, minimize, bounded, solver!, solve!, value
const Opti = C.Opti
minimize(o::Opti, ex) = C.minimize(o, convert(MX, ex isa Union{SX,MX} ? ex : MX(ex)))
bounded(lb, ex::MX, ub) = C.bounded(C.Opti, convert(MX, lb), ex, convert(MX, ub))
# casadi-style positional option dicts (plugin_options, solver_options), NOT kwargs.
solver!(o::Opti, plugin::AbstractString, plugin_opts=Dict{String,Any}(), solver_opts=Dict{String,Any}()) =
    C.solver(o, plugin, plugin_opts, solver_opts)
solve!(o::Opti) = C.solve(o)
# value() comes back as native Julia from the wrapper itself (native_DM typemap):
# scalar Float64, row or column -> Vector, else Matrix -- also via `sol.value(x)`.
value(sol, x) = C.value(sol, x)
export set_initial
set_initial(o::Opti, x, v) = C.set_initial(o, x, v)

# ---- 1-based indexing: Julia ranges -> casadi Slice (0-based stop-exclusive) ----
# Reject i < 1 (Julia/MATLAB are 1-based; never Python-style negative wraparound).
_ckpos(i) = i >= 1 ? i : Base.throw(BoundsError())
_slice(i::Integer) = (_ckpos(i); C.Slice(Int(i) - 1, Int(i)))
_slice(r::AbstractUnitRange) = (_ckpos(Base.first(r)); C.Slice(Int(Base.first(r)) - 1, Int(Base.last(r))))
_slice(r::StepRange) = (Base.step(r) > 0 ?
    C.Slice(Int(Base.first(r)) - 1, Int(Base.last(r)), Int(Base.step(r))) :
    Base.error("negative-step ranges not supported yet"))
_slice(::Colon) = C.Slice()

# Upper-bound guard so OOR throws a Julia BoundsError, not casadi's 0-based assertion (P1-5).
_ckhi(x, i::Integer, n) = i <= n ? i : Base.throw(BoundsError(x, i))
_ckhi(x, r::AbstractRange, n) = (Base.isempty(r) || Base.last(r) <= n || Base.throw(BoundsError(x, r)); r)
_ckhi(x, ::Colon, n) = nothing
# Integers, ranges and `:` index through a casadi Slice; integer vectors and Bool or DM masks
# through explicit 0-based indices (on both axes, once either is a vector).
const _Ind = Union{Integer,AbstractRange,Colon}
const _AnyInd = Union{_Ind,AbstractVector{<:Integer},AbstractArray{Bool},DM}
_imask(x, i, n) = (Base.length(i) == n || Base.throw(BoundsError(x, i)); Base.findall(Base.vec(i)) .- 1)
_ivec(x, i::AbstractArray{Bool}, n) = _imask(x, i, n)
_ivec(x, i::AbstractVector{Bool}, n) = _imask(x, i, n)
_ivec(x, i::DM, n) = _imask(x, Base.Matrix{Float64}(i) .!= 0, n)   # e.g. x[x .> 0]: a mask
_ivec(x, i::AbstractVector{<:Integer}, n) = (Base.all(k -> 1 <= k <= n, i) || Base.throw(BoundsError(x, i)); Int.(i) .- 1)
_ivec(x, i::Integer, n) = _ivec(x, [i], n)
_ivec(x, ::Colon, n) = Base.collect(0:n-1)
# Shapes follow casadi (as in its Python/MATLAB and SciML CasADi.jl), not Julia's dropped
# dimensions: A[i, :] is a 1×n row and x[range] keeps a vector's orientation.
Base.getindex(x::Union{SX,MX,DM}, i::_AnyInd) = i isa _Ind ?
    (_ckhi(x, i, C.numel(x)); C.get(x, false, _slice(i))) :
    C.get(x, false, _ivec(x, i, C.numel(x)))
function Base.getindex(x::Union{SX,MX,DM}, i::_AnyInd, j::_AnyInd)
    i isa _Ind && j isa _Ind ?
        (_ckhi(x, i, C.size1(x)); _ckhi(x, j, C.size2(x)); C.get(x, false, _slice(i), _slice(j))) :
        C.get(x, false, _ivec(x, i, C.size1(x)), _ivec(x, j, C.size2(x)))
end
# setindex! mirrors getindex; coerce Real/array v to the matrix type, a column into a row slot
_fitrow(v) = v isa _Cat && C.size2(v) == 1 && C.size1(v) != 1 ? C.transpose(v) : v
function Base.setindex!(x::Union{SX,MX,DM}, v, i::_AnyInd)
    w = _asval(typeof(x), v)
    C.size1(x) == 1 && (w = _fitrow(w))
    i isa _Ind ? C.set(x, w, false, _slice(i)) : C.set(x, w, false, _ivec(x, i, C.numel(x)))
    x
end
function Base.setindex!(x::Union{SX,MX,DM}, v, i::_AnyInd, j::_AnyInd)
    w = _asval(typeof(x), v)
    i isa Integer && !(j isa Integer) && (w = _fitrow(w))
    i isa _Ind && j isa _Ind ? C.set(x, w, false, _slice(i), _slice(j)) :
        C.set(x, w, false, _ivec(x, i, C.size1(x)), _ivec(x, j, C.size2(x)))
    x
end
# value coercion for set(): scalar Real / native array -> the matrix type, else pass through
_asval(::Type{T}, v::Real) where {T} = convert(T, v)
_asval(::Type{T}, v::AbstractArray{<:Real}) where {T} = T(v)
_asval(::Type{T}, v) where {T} = v
# A symbolic MX index is 1-based like any other; a Sparsity index selects its pattern
Base.getindex(x::MX, k::MX) = C.get(x, true, k)
Base.getindex(x::MX, i::MX, j::MX) = C.get(x, true, i, j)
Base.getindex(x::_Cat, sp::Sparsity) = C.get(x, false, sp)
Base.setindex!(x::_Cat, v, sp::Sparsity) = (C.set(x, _asval(typeof(x), v), false, sp); x)

# x.nz[k]: the structural nonzeros as a 1-based vector (Python's x.nz[k])
struct NZ{T<:_Cat}
    x::T
end
C._swig_getproperty(x::_Cat, ::Val{:nz}) = NZ(x)
Base.length(n::NZ) = C.nnz(n.x)
Base.lastindex(n::NZ) = C.nnz(n.x)
Base.firstindex(n::NZ) = 1
Base.getindex(n::NZ, i::_AnyInd) = i isa _Ind ?
    (_ckhi(n.x, i, C.nnz(n.x)); C.get_nz(n.x, false, _slice(i))) :
    C.get_nz(n.x, false, _ivec(n.x, i, C.nnz(n.x)))
Base.getindex(n::NZ{MX}, k::MX) = C.get_nz(n.x, true, k)
Base.setindex!(n::NZ{MX}, v, k::MX) = (C.set_nz(n.x, _asval(MX, v), true, k); n)
function Base.setindex!(n::NZ, v, i::_AnyInd)
    w = _asval(typeof(n.x), v)
    i isa _Ind ? C.set_nz(n.x, w, false, _slice(i)) : C.set_nz(n.x, w, false, _ivec(n.x, i, C.nnz(n.x)))
    n
end

Base.lastindex(x::Union{SX,MX,DM}) = Base.length(x)
Base.lastindex(x::Union{SX,MX,DM}, d::Integer) = d == 1 ? C.size1(x) : C.size2(x)
Base.firstindex(x::Union{SX,MX,DM}) = 1

# method-style calls (f.jacobian(), o.solver(...), sol.value(x)) are now emitted
# generically by the SWIG backend: a per-class `Base.getproperty` routes every
# instance-method name to the free function, fields (.ptr) short-circuit first.

# ---- serialization (R3b) ----
# High-level names are casadi_*-prefixed to avoid Base/stdlib Serialization piracy.
export StringSerializer, StringDeserializer, pack, unpack, encode,
    casadi_serialize, casadi_deserialize
const StringSerializer = C.StringSerializer
const StringDeserializer = C.StringDeserializer
const encode = C.encode

# SWIG-Julia generates no struct inheritance, so methods are keyed on each
# concrete (de)serializer struct. Treat the base + string + file variants alike.
const _Serializer = Union{C.SerializerBase,C.StringSerializer,C.FileSerializer}
const _Deserializer = Union{C.DeserializerBase,C.StringDeserializer,C.FileDeserializer}

# pin the correct pack slot for casadi vector types: the generated `pack(::Any)`
# probes overloads in order and mis-routes Vector{DM} to the MX-vector slot (the
# can-pack typecheck is ambiguous across DM/SX/MX). Force the matching SWIG slot.
for (T, slot) in ((:Sparsity, 10), (:MX, 11), (:DM, 12), (:SX, 13),
                  (:Linsol, 14), (:Function, 15), (:GenericType, 16))
    sym = QuoteNode(Symbol("_swig_SerializerBase_SerializerBase_pack__SWIG_$slot"))
    @eval pack(s::_Serializer, e::AbstractVector{<:C.$T}) =
        GC.@preserve s begin
            ccall(($sym, C._lib), Cvoid, (Ptr{Cvoid}, Any), s.ptr, e); C._check(nothing)
        end
end
# scalars and primitive vectors already dispatch correctly in the generated module
pack(s::_Serializer, e) = C.pack(s, e)

# Map each SerializationType tag (0..23) to its blind_unpack_* reconstructor,
# mirroring Python/JS `blind_unpack_<type_to_string(tag)>`. Built once from the
# enum so it survives a regen; skips tags whose unpacker the module lacks.
function _build_unpack_table()
    t = Dict{Int,Base.Function}()
    for tag in 0:23
        name = C.type_to_string(C.SerializerBase, tag)   # e.g. "sx", "dm_vector"
        fn = Symbol("blind_unpack_", name)
        Base.isdefined(C, fn) && (t[tag] = Base.getfield(C, fn))
    end
    t
end
const _UNPACK_TABLE = _build_unpack_table()

"""
    unpack(d) -> reconstructed object

Read the next type tag from a deserializer and route to the matching
`blind_unpack_*`, returning the reconstructed casadi object (or Vector thereof).
"""
function unpack(d::_Deserializer)
    tag = C.pop_type(d)
    fn = Base.get(_UNPACK_TABLE, tag, nothing)
    fn === nothing && Base.error("unpack: no handler for serialization tag $tag")
    fn(d)
end

# One-shot string round-trip, mirroring Python's StringSerializer/StringDeserializer.
"""
    casadi_serialize(obj) -> String

Pack a single casadi object (or Vector thereof) into a portable string.
"""
function casadi_serialize(obj)
    s = StringSerializer(Dict{String,Any}())
    pack(s, obj)
    encode(s)
end

"""
    casadi_deserialize(str)

Reconstruct the object packed by [`casadi_serialize`](@ref).
"""
casadi_deserialize(str::AbstractString) = unpack(StringDeserializer(str))

# ---- Callback subclassing (R5) -------------------------------------------
# Subtype `CasadiCallback`, give your type a `cb_eval` (and optionally the
# `cb_n_in`/`cb_n_out`/`cb_sparsity_*` methods below), then wrap it with
# `callback(obj; name=...)` to get a normal casadi `Function` whose evaluation
# dispatches into your Julia `cb_eval`.
export CasadiCallback, callback, cb_eval, cb_n_in, cb_n_out,
    cb_sparsity_in, cb_sparsity_out, cb_name_in, cb_name_out,
    cb_init, cb_finalize, cb_has_jacobian, cb_get_jacobian,
    cb_has_forward, cb_get_forward, cb_has_reverse, cb_get_reverse,
    cb_has_jac_sparsity, cb_get_jac_sparsity
"""
    CasadiCallback

Supertype of a Julia-defined casadi Callback; `callback(obj; name, opts...)` makes it a Function.
Define `cb_eval(obj, args::Vector{DM})` returning a vector of outputs, and optionally
`cb_n_in`, `cb_n_out`, `cb_sparsity_in`, `cb_sparsity_out`, `cb_name_in`, `cb_name_out`
(1-based indices), `cb_init`, `cb_finalize`, and the derivative hooks
`cb_has_jacobian`/`cb_get_jacobian(obj, name, inames, onames, opts)`,
`cb_has_forward`/`cb_get_forward(obj, nfwd, name, inames, onames, opts)`,
`cb_has_reverse`/`cb_get_reverse(obj, nadj, name, inames, onames, opts)` and
`cb_has_jac_sparsity`/`cb_get_jac_sparsity(obj, oind, iind, symmetric)`.
"""
abstract type CasadiCallback end

# Defaults mirroring casadi's Callback base: scalar dense I/O, override as needed.
cb_n_in(::CasadiCallback)  = 1
cb_n_out(::CasadiCallback) = 1
cb_sparsity_in(::CasadiCallback, i)  = C.dense(C.Sparsity, 1, 1)
cb_sparsity_out(::CasadiCallback, i) = C.dense(C.Sparsity, 1, 1)
cb_name_in(::CasadiCallback, i)  = "i$(i)"
cb_name_out(::CasadiCallback, i) = "o$(i)"
# cb_eval has no default: a callback must compute something.
function cb_eval end
# Optional hooks, as in casadi's Callback; unimplemented ones fall back to the C++ base.
const _NO = C._swig_director_no_override
for f in (:cb_init, :cb_finalize, :cb_has_jacobian, :cb_get_jacobian, :cb_has_forward,
          :cb_get_forward, :cb_has_reverse, :cb_get_reverse, :cb_has_jac_sparsity, :cb_get_jac_sparsity)
    @eval $f(::CasadiCallback, args...) = _NO
end

# Route the SWIG director generic functions (invoked from C++) to the user's
# CasadiCallback methods.  Defined once on the abstract type; no ambiguity.
# casadi indices are 0-based; expose 1-based to Julia.
C.Callback_get_n_in(self::CasadiCallback, args...)   = Int64(cb_n_in(self))
C.Callback_get_n_out(self::CasadiCallback, args...)  = Int64(cb_n_out(self))
C.Callback_get_sparsity_in(self::CasadiCallback, i)  = cb_sparsity_in(self, Int(i) + 1)
C.Callback_get_sparsity_out(self::CasadiCallback, i) = cb_sparsity_out(self, Int(i) + 1)
C.Callback_get_name_in(self::CasadiCallback, i)      = cb_name_in(self, Int(i) + 1)
C.Callback_get_name_out(self::CasadiCallback, i)     = cb_name_out(self, Int(i) + 1)
C.Callback_eval(self::CasadiCallback, arg)           = collect(DM, cb_eval(self, collect(DM, arg)))
C.Callback_init(self::CasadiCallback)     = cb_init(self)
C.Callback_finalize(self::CasadiCallback) = cb_finalize(self)
C.Callback_has_jacobian(self::CasadiCallback) = cb_has_jacobian(self)
C.Callback_get_jacobian(self::CasadiCallback, name, inames, onames, opts) =
    _asfunction(cb_get_jacobian(self, name, inames, onames, opts))
C.Callback_has_forward(self::CasadiCallback, nfwd) = cb_has_forward(self, Int(nfwd))
C.Callback_get_forward(self::CasadiCallback, nfwd, name, inames, onames, opts) =
    _asfunction(cb_get_forward(self, Int(nfwd), name, inames, onames, opts))
C.Callback_has_reverse(self::CasadiCallback, nadj) = cb_has_reverse(self, Int(nadj))
C.Callback_get_reverse(self::CasadiCallback, nadj, name, inames, onames, opts) =
    _asfunction(cb_get_reverse(self, Int(nadj), name, inames, onames, opts))
C.Callback_has_jac_sparsity(self::CasadiCallback, oind, iind) =
    cb_has_jac_sparsity(self, Int(oind) + 1, Int(iind) + 1)
C.Callback_get_jac_sparsity(self::CasadiCallback, oind, iind, symmetric) =
    cb_get_jac_sparsity(self, Int(oind) + 1, Int(iind) + 1, symmetric)
# A Callback handed back to casadi (e.g. a Jacobian) must outlive its Julia proxy.
const _kept_callbacks = Base.IdSet{Any}()
_asfunction(r) = r
_asfunction(r::C.Callback) = (Base.push!(_kept_callbacks, r); CasadiFunction(r))

# opti.callback(f): the director calls back with the iteration index.
C.OptiCallback_call(f::Base.Function, i) = (f(Int(i)); nothing)
callback(o::Opti, args...) = C.callback(o, args...)

"""
    callback(obj::CasadiCallback; name="callback", opts...) -> Function

Wrap a `CasadiCallback` into a casadi `Function`.  The returned object calls
back into `cb_eval(obj, args)` (and the `cb_n_in`/`cb_sparsity_*` methods)
whenever the function is evaluated, differentiated or embedded in a graph.
"""
function callback(obj::CasadiCallback; name::AbstractString="callback", opts...)
    proxy = C.Callback(obj)                      # director proxy over `obj`
    C.construct(proxy, name, _optdict(opts))     # builds the internal function
    proxy
end

# Let a Callback be used wherever a casadi Function is expected (nlpsol, graphs).
# The Function copy-ctor accepts the Callback proxy via the SWIG cast chain.
CasadiFunction(cb::C.Callback) = C.Function(C._check(ccall((:_swig_Function_new_Function__SWIG_7,
        C._lib), Ptr{Cvoid}, (Any,), cb)))
Base.convert(::Type{CasadiFunction}, cb::C.Callback) = CasadiFunction(cb)
# A constructed Callback proxy is a casadi Function; call it like one (numeric or symbolic).
(cb::C.Callback)(args...; kwargs...) = CasadiFunction(cb)(args...; kwargs...)

det(a::Union{_Cat,C.Linsol}, args...) = C.det(a, args...)

# ---- full-surface qualified access ----
# `import CasADiNative as ca` -> `ca.<name>` for the ENTIRE generated surface, including
# member-style methods (ca.str, ca.nnz, ca.size1, ...) that casadi does not `export` and
# so aren't reachable by `using` alone. Ergonomic definitions above take precedence (skipped
# via isdefined). Mirrors Python's `import casadi as ca`; `using CasADiNative` is `import *`.
# `@__MODULE__` rather than a literal name: this file is also spliced, unwrapped, into
# the CasADiNative.jl package's own module.
# A casadi method named like a Base function (map, copy, get, min, ...) gets its own
# `ca.` function: casadi's for a casadi object first argument, Base's otherwise.
for _n in names(casadi; all=true)
    _s = Base.string(_n)
    (Base.isempty(_s) || _s[1] == '_' || _s[1] == '#' || _s[1] == '@') && continue
    Base.isdefined(casadi, _n) || continue
    if !Base.isdefined(@__MODULE__, _n)
        @eval const $_n = casadi.$_n
    elseif _n in (:eval, :include) || _n in Base.names(@__MODULE__; all=true)
        continue
    elseif (_f = Base.getfield(@__MODULE__, _n)) !== Base.getfield(casadi, _n) && _f isa Base.Function
        @eval begin
            $_n(a::_Proxy, args...; kw...) = casadi.$_n(a, args...; kw...)
            $_n(args...; kw...) = $_f(args...; kw...)
        end
    end
end

end # module CasADiNative
