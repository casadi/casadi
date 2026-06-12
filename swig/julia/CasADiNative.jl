# CasADiNative.jl -- handwritten ergonomics layer over the SWIG-generated module.
# Prototype of what ships as the package's src/CasADi.jl; the generated
# casadi.jl stays untouched underneath.
module CasADiNative

include(joinpath(@__DIR__, "casadi.jl"))
using .casadi

const C = casadi
export SX, MX, DM, Sparsity, sym, nlpsol, vertcat, jacobian, hessian
# P1-3: high-value builders that live as C.* but were unreachable via `using`
# (`reshape` is deliberately NOT exported: we extend Base.reshape, so a bare
#  `reshape` resolves to Base with our methods; exporting would shadow it. P2-1)
export horzcat, blockcat, diagcat, sum1, sum2, repmat, dot, cross,
    mtimes, gradient, jtimes, substitute, simplify,
    qpsol, conic, integrator, rootfinder, CodeGenerator, codegen,
    external, generate, save
# horzcat/vertcat/blockcat/diagcat/dot/gradient get typed/convenience wrappers below
for _n in (:sum1, :sum2, :repmat, :reshape, :cross, :mtimes, :jtimes,
           :substitute, :simplify)
    @eval const $_n = C.$_n
end
export tangent, factory   # convenience defaults added below
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
Base.convert(::Type{SX}, v::Union{Real,DM}) = C.SX(v isa Real ? Float64(v) : v)
Base.convert(::Type{MX}, v::Union{Real,DM}) = C.MX(v isa Real ? Float64(v) : v)
Base.convert(::Type{DM}, v::Real) = C.DM(Float64(v))

# ---- operators ----
for (op, fn) in ((:+, :plus), (:-, :minus), (:*, :mtimes), (:/, :rdivide), (:^, :power))
    @eval begin
        Base.$op(a::SX, b::SX) = C.$fn(a, b)
        Base.$op(a::MX, b::MX) = C.$fn(a, b)
        Base.$op(a::DM, b::DM) = C.$fn(a, b)
        # promotion: scalars and DM promote to the symbolic side
        Base.$op(a::SX, b::Union{Real,DM}) = C.$fn(a, convert(SX, b))
        Base.$op(a::Union{Real,DM}, b::SX) = C.$fn(convert(SX, a), b)
        Base.$op(a::MX, b::Union{Real,DM}) = C.$fn(a, convert(MX, b))
        Base.$op(a::Union{Real,DM}, b::MX) = C.$fn(convert(MX, a), b)
        Base.$op(a::DM, b::Real) = C.$fn(a, convert(DM, b))
        Base.$op(a::Real, b::DM) = C.$fn(convert(DM, a), b)
    end
end
Base.:-(a::Union{SX,MX,DM}) = C.times(a, typeof(a) === DM ? C.DM(-1.0) : convert(typeof(a), -1))

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
    x ^ convert(typeof(x), n)
# broadcast op -> ELEMENTWISE casadi op: `*`/`/`/`^` must NOT be mtimes/solve/mpower
_bcop(f) = f                       # +,-,math funcs are already elementwise
_bcop(::typeof(*)) = C.times
_bcop(::typeof(/)) = C.rdivide
_bcop(::typeof(^)) = C.power

# materialize by recursively walking the (unflattened) tree, so nested `.* ` inside a
# fused chain still maps to elementwise times rather than Base.:* (mtimes).
_bceval(x) = x
function _bceval(bc::Base.Broadcast.Broadcasted{CasadiStyle})
    args = map(_bceval, bc.args)
    bc.f === Base.literal_pow && return _apply(bc.f, args...)
    # Float64.(dm) etc.: a numeric conversion broadcast over a single casadi matrix
    # escapes to a Julia array of converted scalars (not a casadi op).
    if Base.length(args) == 1 && bc.f isa Type{<:Number} && args[1] isa _Cat
        return bc.f.(Base.collect(args[1]))
    end
    T = _castype(args...)
    _bcop(bc.f)(map(a -> _coerce(T, a), args)...)
end
Base.copy(bc::Base.Broadcast.Broadcasted{CasadiStyle}) = _bceval(bc)

# ---- concatenation: [a b], [a;b], [a b;c d] (P0-3) ----
const _Cat = Union{SX,MX,DM}
_catpromote(xs) = Base.any(x -> x isa MX, xs) ? MX : Base.any(x -> x isa SX, xs) ? SX : DM
_ascat(::Type{T}, x::T) where {T<:_Cat} = x                 # same type: identity
_ascat(::Type{T}, x::_Cat) where {T} = convert(T, x)        # promote DM->SX/MX etc.
_ascat(T, x::Real) = convert(T, x)
_ascat(T, x::AbstractArray{<:Real}) = T(x)
_cvec(xs) = (T = _catpromote(xs); [_ascat(T, x) for x in xs])
# generated vertcat/horzcat always return MX; restore the DM type for numeric cats
_restore(::Type{DM}, r) = C.evalf(r)
_restore(::Type{T}, r) where {T} = r
# Anchor on the first arg being a casadi type so all-Real cats never reach here.
Base.vcat(x::_Cat, xs::Union{_Cat,Real}...) = _restore(_catpromote((x, xs...)), C.vertcat(_cvec((x, xs...))))
Base.hcat(x::_Cat, xs::Union{_Cat,Real}...) = _restore(_catpromote((x, xs...)), C.horzcat(_cvec((x, xs...))))
# Leading Julia scalar + a casadi type later: keeps `[0; x; 0]`, `[c; x]`, `[c x]` from
# silently building a Vector{Any}. (A purely-scalar prefix like `[1; 2; x]` still escapes
# to Base; reach for the free `vertcat(...)` there.) (P1-2)
Base.vcat(a::Real, x::_Cat, xs::Union{_Cat,Real}...) = _restore(_catpromote((a, x, xs...)), C.vertcat(_cvec((a, x, xs...))))
Base.hcat(a::Real, x::_Cat, xs::Union{_Cat,Real}...) = _restore(_catpromote((a, x, xs...)), C.horzcat(_cvec((a, x, xs...))))
# [a b; c d] -> rows of horzcat, then vertcat
function Base.hvcat(rows::Tuple{Vararg{Int}}, x::_Cat, ys::Union{_Cat,Real}...)
    xs = (x, ys...)
    T = _catpromote(xs)
    cv = [_ascat(T, x) for x in xs]
    blocks = Any[]
    k = 0
    for ncol in rows
        Base.push!(blocks, C.horzcat(cv[(k + 1):(k + ncol)]))
        k += ncol
    end
    _restore(T, C.vertcat(blocks))
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
horzcat(args::Union{_Cat,Real}...) = _typedcat(C.horzcat, args...)
horzcat(v::AbstractVector) = _typedcat(C.horzcat, v)
vertcat(args::Union{_Cat,Real}...) = _typedcat(C.vertcat, args...)
vertcat(v::AbstractVector) = _typedcat(C.vertcat, v)
diagcat(args::Union{_Cat,Real}...) = _typedcat(C.diagcat, args...)
diagcat(v::AbstractVector) = _typedcat(C.diagcat, v)
# blockcat takes a 2-D matrix-of-blocks; lower it to vertcat-of-horzcats
function blockcat(rows::AbstractVector{<:AbstractVector})
    xs = Base.collect(Base.Iterators.flatten(rows))
    T = _catpromote(xs)
    _restore(T, C.vertcat([C.horzcat([_ascat(T, x) for x in r]) for r in rows]))
end
blockcat(A::_Cat, B::_Cat, C_::_Cat, D::_Cat) =
    blockcat([[A, B], [C_, D]])
# dot: promote both sides to a common type, then the matching generated method
function dot(a::Union{_Cat,Real}, b::Union{_Cat,Real})
    T = _catpromote((a, b))
    _restore(T, C.dot(_ascat(T, a), _ascat(T, b)))
end
# MATLAB-style backslash solve `A \ b` (linear solve), distinct from elementwise
function Base.:\(a::_Cat, b::Union{_Cat,Real})
    T = _catpromote((a, b))
    _restore(T, C.solve(_ascat(T, a), _ascat(T, b)))
end

# ---- Base math extensions (no shadowing: extend Base methods) ----
for fn in (:sin, :cos, :tan, :exp, :log, :sqrt, :abs, :tanh, :sinh, :cosh,
           :asin, :acos, :atan, :floor, :ceil, :sign)
    @eval begin
        Base.$fn(x::SX) = C.$fn(x)
        Base.$fn(x::MX) = C.$fn(x)
        Base.$fn(x::DM) = C.$fn(x)
    end
end

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
Base.size(x::Union{SX,MX,DM}) = (C.size1(x), C.size2(x))
Base.length(x::Union{SX,MX,DM}) = C.numel(x)
Base.Float64(d::DM) = (n = C.nonzeros(d); Base.length(n) == 1 ? n[1] :
    Base.error("DM is not scalar"))
# P1-5: dense, shape-correct conversions (densify -> column-major nonzeros)
Base.Matrix{Float64}(d::DM) =
    Base.reshape(C.nonzeros(C.densify(d)), (C.size1(d), C.size2(d)))
Base.Matrix(d::DM) = Base.Matrix{Float64}(d)
Base.Array(d::DM) = Base.Matrix{Float64}(d)
Base.collect(d::DM) = Base.Matrix{Float64}(d)
Base.eltype(::Type{DM}) = Float64
# kept for back-compat: structural nonzeros (column-major), ignores shape/sparsity
Base.Vector{Float64}(d::DM) = C.nonzeros(d)

# ---- scalar zero/one (so zeros(SX,m,n)/ones(...) and zero(x)/one(x) work) ----
for T in (:SX, :MX, :DM)
    @eval Base.zero(::Type{$T}) = C.$T(0.0)
    @eval Base.one(::Type{$T})  = C.$T(1.0)
    @eval Base.zero(::$T) = C.$T(0.0)
    @eval Base.one(::$T)  = C.$T(1.0)
end

# ---- Float64 of a CONSTANT scalar SX/MX (DM already has Float64 above) ----
function _const_scalar_float(x::Union{SX,MX})
    (C.is_constant(x) && C.numel(x) == 1) ||
        Base.error("Float64 requires a constant scalar; got $(C.size1(x))×$(C.size2(x)) non-constant")
    Float64(C.evalf(x))   # evalf -> DM -> Float64
end
Base.Float64(x::SX) = _const_scalar_float(x)
Base.Float64(x::MX) = _const_scalar_float(x)

# ---- size(x, dim): 1->rows, 2->cols, else 1 (mirrors Base array semantics) ----
Base.size(x::Union{SX,MX,DM}, d::Integer) =
    d == 1 ? C.size1(x) : d == 2 ? C.size2(x) : 1

# ---- collect / Array: dense Matrix of scalar SX/MX elements (NOT iterable, so
# broadcasting's scalar-leaf treatment via broadcastable(x)=x is preserved) ----
_to_matrix(x::Union{SX,MX}) =
    [x[i, j] for i in 1:C.size1(x), j in 1:C.size2(x)]
Base.collect(x::SX) = _to_matrix(x)
Base.collect(x::MX) = _to_matrix(x)
Base.Array(x::SX) = _to_matrix(x)
Base.Array(x::MX) = _to_matrix(x)

# ---- Base bridges to casadi free functions (typed: preserve DM/SX/MX) ----
Base.vec(x::Union{SX,MX,DM}) = C.vec(x)
Base.kron(a::T, b::T) where {T<:Union{SX,MX,DM}} = C.kron(a, b)
Base.kron(a::_Cat, b::_Cat) = (T = _catpromote((a, b)); C.kron(_ascat(T, a), _ascat(T, b)))
Base.repeat(x::Union{SX,MX,DM}, m::Integer, n::Integer) = C.repmat(x, Int(m), Int(n))
Base.repeat(x::Union{SX,MX,DM}, m::Integer) = C.repmat(x, Int(m), 1)
# Base.reshape on our types -> C.reshape (coexists with the `reshape` const export,
# which is C.reshape itself, so no ambiguity; Base method wins on Union{SX,MX,DM}).
Base.reshape(x::Union{SX,MX,DM}, nrow::Integer, ncol::Integer) = C.reshape(x, Int(nrow), Int(ncol))
Base.reshape(x::Union{SX,MX,DM}, dims::Tuple{Integer,Integer}) = C.reshape(x, Int(dims[1]), Int(dims[2]))

# ---- functor: f(args...) like Python ----
function (f::CasadiFunction)(args::Vararg{Union{Real,DM,SX,MX}})
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
    C.call(f, Dict{String,Any}(Base.string(k) => (v isa Real ? C.DM(Float64(v)) :
                                                  v isa AbstractVector{<:Real} ? C.DM(collect(Float64, v)) : v)
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

# ---- comparisons -> symbolic constraints ----
for (op, fn) in ((:<, :lt), (:<=, :le), (:>=, :ge), (:(==), :eq))
    @eval begin
        Base.$op(a::T, b::T) where {T<:Union{SX,MX}} = C.$fn(a, b)
        Base.$op(a::T, b::Union{Real,DM}) where {T<:Union{SX,MX}} = C.$fn(a, convert(T, b))
        Base.$op(a::Union{Real,DM}, b::T) where {T<:Union{SX,MX}} = C.$fn(convert(T, a), b)
    end
end

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
# Return the solution as native Julia: scalar Float64, column -> Vector, else Matrix.
# (P1-3: the old hard Float64() cast made vector/matrix decision vars unreadable.)
function _dm_native(d::DM)
    (C.size1(d) == 1 && C.size2(d) == 1) && return Float64(d)
    M = Base.Matrix{Float64}(d)
    C.size2(d) == 1 ? M[:, 1] : M
end
value(sol, x) = _dm_native(C.value(sol, x, MX[]))
export set_initial
set_initial(o::Opti, x, v) =
    C.set_initial(o, x, v isa Union{MX,DM} ? v : C.DM(Float64(v)))

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
Base.getindex(x::Union{SX,MX,DM}, i::Union{Integer,AbstractRange,Colon}) =
    (_ckhi(x, i, C.numel(x)); C.get(x, false, _slice(i)))
Base.getindex(x::Union{SX,MX,DM}, i::Union{Integer,AbstractRange,Colon},
              j::Union{Integer,AbstractRange,Colon}) =
    (_ckhi(x, i, C.size1(x)); _ckhi(x, j, C.size2(x)); C.get(x, false, _slice(i), _slice(j)))
Base.setindex!(x::Union{SX,MX,DM}, v, i::Union{Integer,AbstractRange,Colon}) =
    (C.set(x, _asval(typeof(x), v), false, _slice(i)); x)
# 2-D setindex! mirrors the 2-index getindex; coerce Real/array v to the matrix type
Base.setindex!(x::Union{SX,MX,DM}, v, i::Union{Integer,AbstractRange,Colon},
               j::Union{Integer,AbstractRange,Colon}) =
    (C.set(x, _asval(typeof(x), v), false, _slice(i), _slice(j)); x)
# value coercion for set(): scalar Real / native array -> the matrix type, else pass through
_asval(::Type{T}, v::Real) where {T} = convert(T, v)
_asval(::Type{T}, v::AbstractArray{<:Real}) where {T} = T(v)
_asval(::Type{T}, v) where {T} = v
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
    cb_sparsity_in, cb_sparsity_out, cb_name_in, cb_name_out
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

# A constructed Callback proxy is a casadi Function; evaluate it like one.
# Pin the plain DM-vector call slot (the generated overloads collapse to one
# Julia signature, last-wins picks the wrong dict form).
function (cb::C.Callback)(args::Vararg{Union{Real,DM,SX,MX}})
    ins = Any[a isa Real ? C.DM(Float64(a)) : a for a in args]
    r = GC.@preserve cb C._check(ccall((:_swig_Function_Function_call__SWIG_0, C._lib),
            Any, (Any, Any, Any, Any), cb, ins, false, false))
    length(r) == 1 ? r[1] : Tuple(r)
end
# Let a Callback be used wherever a casadi Function is expected (nlpsol, graphs).
# The Function copy-ctor accepts the Callback proxy via the SWIG cast chain.
CasadiFunction(cb::C.Callback) = C.Function(C._check(ccall((:_swig_Function_new_Function__SWIG_8,
        C._lib), Ptr{Cvoid}, (Any,), cb)))
Base.convert(::Type{CasadiFunction}, cb::C.Callback) = CasadiFunction(cb)

# ---- full-surface qualified access ----
# `import CasADiNative as ca` -> `ca.<name>` for the ENTIRE generated surface, including
# member-style methods (ca.str, ca.nnz, ca.size1, ...) that casadi does not `export` and
# so aren't reachable by `using` alone. Ergonomic definitions above take precedence (skipped
# via isdefined). Mirrors Python's `import casadi as ca`; `using CasADiNative` is `import *`.
for _n in names(casadi; all=true)
    _s = Base.string(_n)
    (Base.isempty(_s) || _s[1] == '_' || _s[1] == '#' || _s[1] == '@') && continue
    (Base.isdefined(CasADiNative, _n) || !Base.isdefined(casadi, _n)) && continue
    @eval const $_n = casadi.$_n
end

end # module CasADiNative
