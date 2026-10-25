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
export blockcat, diagcat, veccat, sum1, sum2, repmat, dot, cross, norm, det, tr,
    mtimes, gradient, jtimes, substitute, simplify,
    qpsol, conic, integrator, rootfinder, CodeGenerator, codegen,
    external, generate, save
# concatenation is [a b] / [a; b] (no horzcat/vertcat in the Julia API)
export tangent, factory
# `Function` collides with Base.Function (abstract supertype): it is deliberately
# NOT exported. Under a plain `using`, bare `Function` resolves to Base.Function
# ("no constructors have been defined for Function"); access it qualified
# (`ca.Function`) or pull it in explicitly (`using CasADiNative: Function`).
const Function = C.Function

# ---- conversion (promotion by DM < SX < MX rank); generated ctors handle reals ----
# Julia numbers and arrays promote exactly like DM (`_Num`).
const _Cat = Union{SX,MX,DM}
const _Native = Union{Real,AbstractArray{<:Real}}
const _Num = Union{_Native,DM}
const _Proxy = Union{(t for (_, t) in casadi._proxy_types)...}
Base.convert(::Type{SX}, v::_Num) = C.SX(v)
Base.convert(::Type{MX}, v::_Num) = C.MX(v)
Base.convert(::Type{DM}, v::_Native) = C.DM(v)

# ---- operators and comparisons: the generated overloads promote numbers, arrays and DM ----
# Comparisons are elementwise, in the casadi type (SX/MX constraints, DM 0/1 values).
for (op, fn) in ((:+, :plus), (:-, :minus), (:*, :mtimes),
                 (:<, :lt), (:<=, :le), (:>=, :ge), (:(==), :eq), (:(!=), :ne))
    @eval begin
        Base.$op(a::T, b::Union{T,_Num}) where {T<:Union{SX,MX}} = C.$fn(a, b)
        Base.$op(a::_Num, b::Union{SX,MX}) = C.$fn(a, b)
        Base.$op(a::DM, b::_Num) = C.$fn(a, b)
        Base.$op(a::_Native, b::DM) = C.$fn(a, b)
    end
end
Base.:-(a::_Cat) = C.times(a, -1)

# Where Julia gives an array operation a meaning (A/B, A^n, matrix functions, max of
# arrays), a casadi matrix does the same or throws; never an elementwise surprise.
# Elementwise is spelled with dots (./ .^ sin.), as in Julia.
_scal(x::_Cat) = C.numel(x) == 1
_scal(x::Real) = true
_scal(x) = false
_shape(x) = "$(C.size1(x))×$(C.size2(x))"
Base.:/(a::_Cat, b::Union{_Cat,_Native}) = _div(a, b)
Base.:/(a::_Native, b::_Cat) = _div(a, b)
function _div(a, b)
    _scal(b) && return C.rdivide(a, b)
    _scal(a) && Base.throw(ArgumentError("x / A with a matrix A has no meaning in Julia; use x ./ A"))
    C.mrdivide(a, b)
end
Base.:^(a::_Cat, b::Union{_Cat,Real}) = _pow(a, b)
Base.:^(a::Real, b::_Cat) = _pow(a, b)
function _pow(a, b)
    _scal(a) && _scal(b) && return C.power(a, b)
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
Base.broadcastable(x::_Cat) = x   # treat as a scalar leaf, not iterable
Base.BroadcastStyle(::Type{<:_Cat}) = CasadiStyle()
Base.BroadcastStyle(s::CasadiStyle, ::Base.Broadcast.BroadcastStyle) = s
Base.BroadcastStyle(s::CasadiStyle, ::CasadiStyle) = s
Base.BroadcastStyle(s::CasadiStyle, ::Base.Broadcast.Unknown) = s  # disambiguate

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
    op === nothing ? _bcelem(bc.f, args) : op(args...)
end
# any other f is applied to each element, with Julia's broadcast shape rules
function _bcelem(f, args)
    r = Base.broadcast(f, Base.map(a -> a isa _Cat ? Base.collect(a) : a, args)...)
    r isa AbstractArray || return r
    Base.isempty(r) && return DM(Base.size(r, 1), Base.size(r, 2))
    C.reshape(C._vertcat(Base.vec(r)), Base.size(r, 1), Base.size(r, 2))
end
Base.copy(bc::Base.Broadcast.Broadcasted{CasadiStyle}) = _bceval(bc)

# ---- concatenation: [a b], [a;b], [a b;c d] (P0-3) ----
# The generated builders take numbers, Julia arrays and casadi values and pick DM, SX or MX.
_cat(builder, xs) = builder(_catcheck(xs, builder === C._horzcat ? 1 : 2))
# Julia's DimensionMismatch for pieces that do not line up (empty pieces are skipped, as in casadi)
function _catcheck(cv, d)
    n = [Base.size(x, d) for x in cv if Base.length(x) > 0]
    Base.all(==(Base.first(n, 1)...), n) || Base.throw(DimensionMismatch(
        "number of $(d == 1 ? "rows" : "columns") of each array must match (got $(Base.Tuple(n)))"))
    cv
end
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
    blocks = Any[]
    k = 0
    for ncol in rows
        Base.push!(blocks, C._horzcat(_catcheck(xs[(k + 1):(k + ncol)], 1)))
        k += ncol
    end
    C._vertcat(_catcheck(blocks, 2))
end

# ---- function-form cat ----
_catargs(xs::Tuple{AbstractVector}) = Base.collect(xs[1])  # single-vector call form
_catargs(xs::Tuple) = Base.collect(xs)                     # varargs call form
# no pieces: an empty matrix of the vector's element type (DM by default)
_emptycat(xs, n) = (Base.eltype(xs) !== Union{} && Base.eltype(xs) <: _Cat ? Base.eltype(xs) : DM)(0, n)
function _typedcat(builder, args...)
    xs = _catargs(args)
    Base.isempty(xs) && return _emptycat(args isa Tuple{AbstractVector} ? args[1] : xs, builder === C.veccat ? 1 : 0)
    builder(xs)
end
# Varargs may mix casadi values, numbers and Julia arrays; a single Julia vector
# is the list form (its elements are the pieces), as in Python.
diagcat(args::_CatArg...) = _typedcat(C.diagcat, args...)
diagcat(v::AbstractVector) = _typedcat(C.diagcat, v)
veccat(args::_CatArg...) = _typedcat(C.veccat, args...)
veccat(v::AbstractVector) = _typedcat(C.veccat, v)
# blockcat takes a 2-D matrix-of-blocks; lower it to vertcat-of-horzcats
function blockcat(rows::AbstractVector{<:AbstractVector})
    C._vertcat([C._horzcat(Base.collect(Any, r)) for r in rows])
end
blockcat(A::_CatArg, B::_CatArg, C_::_CatArg, D::_CatArg) =
    blockcat([[A, B], [C_, D]])
dot(a::_Cat, b::Union{_Cat,_Native}) = C.dot(a, b)
dot(a::_Native, b::_Cat) = C.dot(a, b)
cross(a::T, b::T, dim::Integer...) where {T<:_Cat} = C.cross(a, b, dim...)
# MATLAB-style backslash solve `A \ b` (linear solve), distinct from elementwise
Base.:\(a::_Cat, b::Union{_Cat,_Native}) = C.solve(a, b)
Base.:\(a::AbstractMatrix{<:Real}, b::_Cat) = C.solve(a, b)

# ---- Base math extensions (no shadowing: extend Base methods) ----
# Julia's sin(A), exp(A), ... of a square matrix are matrix functions, which casadi does
# not compute: a non-scalar throws (use sin.(A)).  Scalars (1×1) are plain casadi ops.
for fn in (:sin, :cos, :tan, :exp, :log, :sqrt, :tanh, :sinh, :cosh,
           :asin, :acos, :atan, :asinh, :acosh, :atanh)
    @eval Base.$fn(x::_Cat) = _scal(x) ? C.$fn(x) : Base.throw(ArgumentError(
        "$($(QuoteNode(fn)))(A) of a $(_shape(x)) $(typeof(x)): in Julia a matrix function of a " *
        "square A, undefined otherwise; use $($(QuoteNode(fn))).(A) for elementwise, or a Julia Matrix(A)"))
end
# no array meaning in Julia for these: elementwise, as casadi
for fn in (:abs, :floor, :ceil, :sign, :log10, :log1p, :expm1)
    @eval Base.$fn(x::_Cat) = C.$fn(x)
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
    cf(a, b)
end
Base.:!(a::_Cat) = _scal(a) ? C.logic_not(a) :
    Base.throw(ArgumentError("!A of a $(_shape(a)) $(typeof(a)) has no meaning in Julia; use .!A"))

# Integer-style division of scalars with Julia's rounding (rem/% truncates, mod floors);
# elementwise is rem.(A, B) etc., as for Julia arrays.
function _scalar2(f, a, b)
    _scal(a) && _scal(b) || Base.throw(ArgumentError("$f of arrays has no meaning in Julia; use $f.(A, B)"))
    a, b
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
Base.transpose(x::_Cat) = C.transpose(x)
Base.adjoint(x::_Cat) = C.transpose(x)   # real/symbolic: adjoint == transpose

# ---- show: casadi's own str for print/interpolation, its repr TypeName(str) for display, as Python ----
const _MIME = MIME"text/plain"
_type_name(x::T) where {T} = Base.applicable(C.type_name, x) ? C.type_name(x) :
    Base.applicable(C.type_name, T) ? C.type_name(T) : Base.string(nameof(T))
for (_, T) in casadi._proxy_types
    Base.hasmethod(C.str, Tuple{T}) || continue
    @eval Base.show(io::IO, x::$T) = print(io, C.str(x))
    @eval Base.show(io::IO, ::_MIME, x::$T) = print(io, _type_name(x), "(", C.str(x), ")")
end

# casadi types print by their name, not their internal module path (Vector{MX}, MethodErrors)
for (_, T) in casadi._proxy_types
    @eval Base.show(io::IO, ::Type{$T}) = print(io, $(Base.string(nameof(T))))
    @eval function Base.show(io::IO, v::Vector{$T})
        print(io, $(Base.string(nameof(T))), "[")
        for (k, x) in enumerate(v)
            k > 1 && print(io, ", ")
            Base.show(io, x)
        end
        print(io, "]")
    end
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
Base.vec(x::_Cat) = C.vec(x)
Base.kron(a::_Cat, b::_Cat) = C.kron(a, b)
Base.repeat(x::_Cat, m::Integer, n::Integer) = C.repmat(x, Int(m), Int(n))
Base.repeat(x::_Cat, m::Integer) = C.repmat(x, Int(m), 1)
# Base.reshape on our types -> C.reshape (coexists with the `reshape` const export,
# which is C.reshape itself, so no ambiguity; Base method wins on _Cat).
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
_outs(r) = Base.length(r) == 1 ? r[1] : Base.Tuple(r)
(f::Function)(args::Vararg{Union{Real,DM,SX,MX,AbstractArray}}) = _outs(C.call(f, Any[args...]))
(f::Function)(d::Dict{String,<:Any}) = C.call(f, Dict{String,Any}(d))
# f() (no args) is the positional zero-input call, not an empty dict-call
(f::Function)(; kwargs...) = Base.isempty(kwargs) ? _outs(C.call(f, Any[])) :
    C.call(f, Dict{String,Any}(Base.string(k) => v for (k, v) in kwargs))

_optdict(opts) = Dict{String,Any}(Base.string(k) => v for (k, v) in opts)
const codegen = C.CodeGenerator   # short alias; CodeGenerator(name[, opts]) opts-defaulted
export interpolant, generate_dependencies

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
export Opti, variable, parameter, subject_to, minimize, bounded, solver, solve, value
# variable/parameter/subject_to record the user's call site ("stacktrace" metadata,
# shown by describe), as in Python and MATLAB
_ownfile(f) = (f = Base.string(f); f == (@__FILE__) || Base.endswith(f, "casadi.jl") || Base.startswith(f, "."))
function _stacktrace(depth = typemax(Int))
    st = Dict{String,Any}[]
    for fr in Base.stacktrace()
        (fr.from_c || _ownfile(fr.file)) && continue
        Base.push!(st, Dict{String,Any}("file" => Base.string(fr.file), "line" => fr.line, "name" => Base.string(fr.func)))
        Base.length(st) >= depth && break
    end
    st
end
for f in (:variable, :parameter)
    @eval function $f(o::Opti, args...)
        r = C.$(Symbol("_", f))(o, args...)
        st = _stacktrace(1)
        Base.isempty(st) || C.update_user_dict(o, r, Dict{String,Any}("stacktrace" => st))
        r
    end
end
function subject_to(o::Opti, args...)
    Base.isempty(args) && return C._subject_to(o)
    meta = Dict{String,Any}("stacktrace" => _stacktrace())
    a = Any[args...]
    if a[end] isa AbstractDict
        a[end] = Base.merge(meta, a[end])
    else
        Base.push!(a, meta)
    end
    C._subject_to(o, a...)
end
for f in (:variable, :parameter, :subject_to)
    @eval C._swig_getproperty(o::Opti, ::Val{$(QuoteNode(f))}) = C._SwigBound($f, o)
end
bounded(lb, ex::MX, ub) = C.bounded(C.Opti, convert(MX, lb), ex, convert(MX, ub))
# value() comes back as native Julia (native_DM typemap): scalar Float64, row or
# column -> Vector, else Matrix -- also via `sol.value(x)`.
export set_initial

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
Base.getindex(x::_Cat, i::_AnyInd) = i isa _Ind ?
    (_ckhi(x, i, C.numel(x)); C.get(x, false, _slice(i))) :
    C.get(x, false, _ivec(x, i, C.numel(x)))
function Base.getindex(x::_Cat, i::_AnyInd, j::_AnyInd)
    i isa _Ind && j isa _Ind ?
        (_ckhi(x, i, C.size1(x)); _ckhi(x, j, C.size2(x)); C.get(x, false, _slice(i), _slice(j))) :
        C.get(x, false, _ivec(x, i, C.size1(x)), _ivec(x, j, C.size2(x)))
end
# setindex! mirrors getindex; coerce Real/array v to the matrix type, a column into a row slot
_fitrow(v) = v isa _Cat && C.size2(v) == 1 && C.size1(v) != 1 ? C.transpose(v) : v
function Base.setindex!(x::_Cat, v, i::_AnyInd)
    w = _asval(typeof(x), v)
    C.size1(x) == 1 && (w = _fitrow(w))
    i isa _Ind ? C.set(x, w, false, _slice(i)) : C.set(x, w, false, _ivec(x, i, C.numel(x)))
    x
end
function Base.setindex!(x::_Cat, v, i::_AnyInd, j::_AnyInd)
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

Base.lastindex(x::_Cat) = Base.length(x)
Base.lastindex(x::_Cat, d::Integer) = d == 1 ? C.size1(x) : C.size2(x)
Base.firstindex(x::_Cat) = 1

# ---- serialization (R3b) ----
# High-level names are casadi_*-prefixed to avoid Base/stdlib Serialization piracy.
export StringSerializer, StringDeserializer, pack, unpack, encode,
    casadi_serialize, casadi_deserialize

"""
    unpack(d) -> reconstructed object

Read the next type tag from a deserializer and route to the matching
`blind_unpack_*`, returning the reconstructed casadi object (or Vector thereof).
"""
function unpack(d::Union{C.DeserializerBase,C.StringDeserializer,C.FileDeserializer})
    tag = C.pop_type(d)
    fn = Symbol("blind_unpack_", C.type_to_string(C.SerializerBase, tag))
    Base.isdefined(C, fn) || Base.error("unpack: no handler for serialization tag $tag")
    Base.getfield(C, fn)(d)
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

# ---- Callback subclassing ----
export Callback, callback
@doc """
    Callback

Abstract supertype of a Julia-defined casadi Callback; `callback(obj, name[, opts])` makes it a Function.
Define `ca.eval_(obj, arg::Vector{DM})` returning a vector of outputs, and optionally any other
virtual of casadi's Callback, with its C++ name and 0-based indices: `ca.get_n_in`, `ca.get_n_out`,
`ca.get_sparsity_in(obj, i)`, `ca.get_sparsity_out`, `ca.get_name_in`, `ca.get_name_out`,
`ca.init`, `ca.finalize`, `ca.uses_output`, `ca.has_eval_buffer`/`ca.eval_buffer(obj, arg, res)`,
`ca.get_factory`, and the derivative hooks
`ca.has_jacobian`/`ca.get_jacobian(obj, name, inames, onames, opts)`,
`ca.has_forward`/`ca.get_forward(obj, nfwd, name, inames, onames, opts)`,
`ca.has_reverse`/`ca.get_reverse(obj, nadj, name, inames, onames, opts)` and
`ca.has_jac_sparsity`/`ca.get_jac_sparsity(obj, oind, iind, symmetric)`.
A method you do not define runs casadi's default.
""" C.Callback
# casadi's finalize hook, not Base's
const finalize = C.finalize

# opti.callback(f): any callable f(i) runs every iteration; the Opti proxy keeps it alive
struct _OptiFunction <: C._OptiCallbackDirector
    f::Any
end
C._call(o::_OptiFunction, i) = o.f(i)
callback(o::Union{Opti,C.OptiAdvanced}, f) = C._callback_class(o, C._OptiCallback(_OptiFunction(f)))
callback(o::Union{Opti,C.OptiAdvanced}) = C._callback_class(o)
C._swig_getproperty(o::Union{Opti,C.OptiAdvanced}, ::Val{:callback}) = C._SwigBound(callback, o)

"""
    callback(obj::Callback, name[, opts::Dict]) -> CallbackFunction

Wrap `obj` into a casadi `CallbackFunction`, usable as a `Function`.  It calls
back into `eval_(obj, arg)` (and the other methods defined for `obj`) whenever the
function is evaluated, differentiated or embedded in a graph.  Keep it alive while
any Function uses it.
"""
function callback(obj::C.Callback, name::AbstractString, opts::AbstractDict = Dict{String,Any}())
    proxy = C.CallbackFunction(obj)              # director proxy over `obj`
    C.construct(proxy, name, _optdict(opts))     # builds the internal function
    proxy
end

Base.convert(::Type{Function}, cb::C.CallbackFunction) = Function(cb)
# A constructed CallbackFunction is a casadi Function; call it like one (numeric or symbolic).
(cb::C.CallbackFunction)(args...; kwargs...) = Function(cb)(args...; kwargs...)

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
