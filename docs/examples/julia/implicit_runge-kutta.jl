#
#     MIT No Attribution
#
#     Copyright (C) 2010-2023 Joel Andersson, Joris Gillis, Moritz Diehl, KU Leuven.
#
#     Permission is hereby granted, free of charge, to any person obtaining a copy of this
#     software and associated documentation files (the "Software"), to deal in the Software
#     without restriction, including without limitation the rights to use, copy, modify,
#     merge, publish, distribute, sublicense, and/or sell copies of the Software, and to
#     permit persons to whom the Software is furnished to do so.
#
#     THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED,
#     INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A
#     PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
#     HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
#     OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
#     SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
#
# Construct a fixed-step implicit Runge-Kutta (collocation) integrator by hand,
# using a newton rootfinder to solve the collocation equations, then exercise
# forward/adjoint sensitivities via factory.
#
# The python reference cvodes integrator is dropped (no sundials here); the
# block-extraction uses 1-based slicing instead of vertsplit (which the Julia
# binding does not yet marshal).
const _jl_dir = get(ENV, "CASADI_JL",
    joinpath(@__DIR__, "..", "..", "..", "julia-proto", "casadi-contact"))
include(joinpath(_jl_dir, "CasADiNative.jl"))
import .CasADiNative as ca

# End time
tf = 10.0

# Dimensions
nx = 3
nu = 1

# Declare variables
x = ca.SX.sym("x", nx)   # state
p = ca.SX.sym("u", nu)   # control

# ODE right hand side function
ode = ca.vertcat((1 - x[2] * x[2]) * x[1] - x[2] + p,
              x[1],
              x[1] * x[1] + x[2] * x[2] + p * p)
f = ca.Function("f", [x, p], [ode])

# Number of finite elements
n = 100
# Size of the finite elements
h = tf / n
# Degree of interpolating polynomial
d = 4

# Choose collocation points
tau_root = vcat(0.0, ca.collocation_points(d, "legendre"))

# Coefficients of the collocation (Coll) and continuity (Cont) equations
Coll = zeros(d + 1, d + 1)
Cont = zeros(d + 1)

# Dimensionless time inside one control interval
tau = ca.SX.sym("tau")

function build_F()
    # For all collocation points, construct the Lagrange basis
    for j in 1:d + 1
        L = ca.SX(1.0)
        for r in 1:d + 1
            if r != j
                L = L * (tau - tau_root[r]) / (tau_root[j] - tau_root[r])
            end
        end
        # Continuity coefficient: value of the polynomial at the final time
        lfcn = ca.Function("lfcn", [tau], [L])
        Cont[j] = Float64(lfcn(1.0))
        # Collocation coefficients: derivative at each collocation point
        tfcn = ca.Function("tfcn", [tau], [ca.tangent(L, tau)])
        for r in 1:d + 1
            Coll[j, r] = Float64(tfcn(tau_root[r]))
        end
    end

    # Variables for one finite element
    X0 = ca.MX.sym("X0", nx)
    P = ca.MX.sym("P", nu)
    V = ca.MX.sym("V", d * nx)

    # State at each collocation point: X[1]==X0, X[r] is block (r-1) of V
    block(r) = V[(r - 2) * nx + 1:(r - 1) * nx]   # r in 2:d+1
    Xcp = [r == 1 ? X0 : block(r) for r in 1:d + 1]

    # Collocation equations that implicitly define V
    V_eq = ca.MX[]
    for j in 2:d + 1
        # State derivative at the collocation point
        xp_j = Coll[1, j] * Xcp[1]
        for r in 2:d + 1
            xp_j = xp_j + Coll[r, j] * Xcp[r]
        end
        push!(V_eq, h * f(Xcp[j], P) - xp_j)
    end

    # Root-finding function, implicitly defines V as a function of X0 and P
    vfcn = ca.Function("vfcn", [V, X0, P], [ca.vertcat(V_eq)])
    vfcn_sx = ca.expand(vfcn)   # convert to SX to decrease overhead

    # Implicit function instance solving the collocation system
    ifcn = ca.rootfinder("ifcn", "newton", vfcn_sx)
    Vsol = ifcn(ca.MX(), X0, P)

    # Recover states and form the end-of-element state
    Xend = [r == 1 ? X0 : Vsol[(r - 2) * nx + 1:(r - 1) * nx] for r in 1:d + 1]
    XF = Cont[1] * Xend[1]
    for r in 2:d + 1
        XF = XF + Cont[r] * Xend[r]
    end

    # Discrete time dynamics for one finite element
    Fel = ca.Function("F", [X0, P], [XF])

    # Iterate over all finite elements
    Xacc = X0
    for i in 1:n
        Xacc = Fel(Xacc, P)
    end

    # Fixed-step integrator as a plain Function over x0, p -> xf
    ca.Function("irk_integrator", [X0, P], [Xacc], ["x0", "p"], ["xf"])
end

function main()
    irk = build_F()

    # Test values
    x0_val = ca.DM([0.0, 1.0, 0.0])
    p_val = 0.2

    println("-------")
    println("Testing ", ca.name(irk))
    println("-------")

    # Forward and reverse directional derivatives
    dF = ca.factory(irk, "dF",
                 ["x0", "p", "fwd:x0", "fwd:p", "adj:xf"],
                 ["xf", "fwd:xf", "adj:x0", "adj:p"])

    res = dF(x0 = x0_val, p = p_val,
             fwd_x0 = [1, 0, 0], fwd_p = 1,
             adj_xf = [0, 0, 1])

    println("xf = ", res["xf"])
    println("d(xf)/d(p)+d(xf)/d(x0[0]) = ", res["fwd_xf"])
    println("d(xf[2])/d(x0) = ", res["adj_x0"])
    println("d(xf[2])/d(p) = ", res["adj_p"])
end

main()
