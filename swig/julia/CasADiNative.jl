# CasADiNative.jl -- handwritten ergonomics layer over the SWIG-generated module.
# Prototype of what ships as the package's src/CasADi.jl; the generated
# casadi.jl stays untouched underneath.
module CasADiNative

include(joinpath(@__DIR__, "casadi.jl"))
using .casadi

const C = casadi
export SX, MX, DM, Sparsity, CasadiFunction, sym, nlpsol, vertcat, jacobian, hessian

const CasadiFunction = casadi.Function

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
# broadcasting -> elementwise casadi ops (x .* y, sin.(x))
Base.Broadcast.broadcasted(::typeof(*), a::T, b::T) where {T<:Union{SX,MX,DM}} = C.times(a, b)
Base.Broadcast.broadcasted(::typeof(+), a::T, b::T) where {T<:Union{SX,MX,DM}} = C.plus(a, b)

# ---- Base math extensions (no shadowing: extend Base methods) ----
for fn in (:sin, :cos, :tan, :exp, :log, :sqrt, :abs, :tanh, :sinh, :cosh,
           :asin, :acos, :atan, :floor, :ceil, :sign)
    @eval begin
        Base.$fn(x::SX) = C.$fn(x)
        Base.$fn(x::MX) = C.$fn(x)
        Base.$fn(x::DM) = C.$fn(x)
    end
end

# ---- show ----
Base.show(io::IO, x::SX) = print(io, C.str(x))
Base.show(io::IO, x::MX) = print(io, C.str(x))
Base.show(io::IO, x::DM) = print(io, C.str(x))
Base.show(io::IO, f::CasadiFunction) = print(io, C.str(f))

# ---- size / conversion ----
Base.size(x::Union{SX,MX,DM}) = (C.size1(x), C.size2(x))
Base.length(x::Union{SX,MX,DM}) = C.numel(x)
Base.Float64(d::DM) = (n = C.nonzeros(d); Base.length(n) == 1 ? n[1] :
    Base.error("DM is not scalar"))
Base.Vector{Float64}(d::DM) = C.nonzeros(d)

# ---- functor: f(args...) like Python ----
function (f::CasadiFunction)(args::Vararg{Union{Real,DM,SX,MX}})
    ins = Any[a isa Real ? C.DM(Float64(a)) : a for a in args]
    r = C.call(f, ins)
    Base.length(r) == 1 ? r[1] : Base.Tuple(r)
end
(f::CasadiFunction)(d::Dict{String,<:Any}) = C.call(f, Dict{String,Any}(d))
(f::CasadiFunction)(; kwargs...) =
    C.call(f, Dict{String,Any}(Base.string(k) => (v isa Real ? C.DM(Float64(v)) :
                                                  v isa AbstractVector{<:Real} ? C.DM(collect(Float64, v)) : v)
                               for (k, v) in kwargs))

# ---- convenience re-exports with casadi names ----
sym(T, args...) = C.sym(T, args...)
vertcat(xs::AbstractVector) = C.vertcat(xs)
vertcat(xs...) = C.vertcat(Base.collect(xs))
jacobian(ex, arg) = C.jacobian(ex, arg, Dict{String,Any}())
hessian(ex, arg) = C.hessian(ex, arg, Dict{String,Any}())
nlpsol(name, plugin, nlp; opts...) =
    C.nlpsol(name, plugin, Dict{String,Any}(nlp), Dict{String,Any}(Base.string(k) => v for (k, v) in opts))

# ---- comparisons -> symbolic constraints ----
for (op, fn) in ((:<, :lt), (:<=, :le), (:>=, :ge), (:(==), :eq))
    @eval begin
        Base.$op(a::T, b::T) where {T<:Union{SX,MX}} = C.$fn(a, b)
        Base.$op(a::T, b::Union{Real,DM}) where {T<:Union{SX,MX}} = C.$fn(a, convert(T, b))
        Base.$op(a::Union{Real,DM}, b::T) where {T<:Union{SX,MX}} = C.$fn(convert(T, a), b)
    end
end

# ---- Opti conveniences ----
export Opti, variable, parameter, minimize, subject_to, bounded, solver!, solve!, value
const Opti = C.Opti
variable(o::Opti, args...) = C.variable(o, args...)
parameter(o::Opti, args...) = C.parameter(o, args...)
minimize(o::Opti, ex) = C.minimize(o, convert(MX, ex isa Union{SX,MX} ? ex : MX(ex)))
minimize(o::Opti, ex::MX) = C.minimize(o, ex)
subject_to(o::Opti, con::MX) = C.subject_to(o, con, Dict{String,Any}())
bounded(lb, ex::MX, ub) = C.bounded(C.Opti, convert(MX, lb), ex, convert(MX, ub))
solver!(o::Opti, plugin::AbstractString; opts...) =
    C.solver(o, plugin, Dict{String,Any}(Base.string(k) => v for (k, v) in opts), Dict{String,Any}())
solve!(o::Opti) = C.solve(o)
value(sol, x) = Float64(C.value(sol, x, MX[]))
export set_initial
set_initial(o::Opti, x, v) =
    C.set_initial(o, x, v isa Union{MX,DM} ? v : C.DM(Float64(v)))

# ---- 1-based indexing: Julia ranges -> casadi Slice (0-based stop-exclusive) ----
_slice(i::Integer) = C.Slice(Int(i) - 1, Int(i))
_slice(r::AbstractUnitRange) = C.Slice(Int(Base.first(r)) - 1, Int(Base.last(r)))
_slice(r::StepRange) = (Base.step(r) > 0 ?
    C.Slice(Int(Base.first(r)) - 1, Int(Base.last(r)), Int(Base.step(r))) :
    Base.error("negative-step ranges not supported yet"))
_slice(::Colon) = C.Slice()

Base.getindex(x::Union{SX,MX,DM}, i::Union{Integer,AbstractRange,Colon}) =
    C.get(x, false, _slice(i))
Base.getindex(x::Union{SX,MX,DM}, i::Union{Integer,AbstractRange,Colon},
              j::Union{Integer,AbstractRange,Colon}) =
    C.get(x, false, _slice(i), _slice(j))
Base.setindex!(x::Union{SX,MX,DM}, v, i::Union{Integer,AbstractRange,Colon}) =
    (C.set(x, v isa Real ? convert(typeof(x), v) : v, false, _slice(i)); x)
Base.lastindex(x::Union{SX,MX,DM}) = Base.length(x)
Base.lastindex(x::Union{SX,MX,DM}, d::Integer) = d == 1 ? C.size1(x) : C.size2(x)
Base.firstindex(x::Union{SX,MX,DM}) = 1

end # module CasADiNativeNative
