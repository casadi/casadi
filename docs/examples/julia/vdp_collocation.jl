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
# Van der Pol OCP by direct collocation on Radau points, assembling the
# collocation NLP by hand from one MX decision vector V. numpy poly1d helpers
# become small array-based polynomial routines (coeffs highest-degree-first).
const _jl_dir = get(ENV, "CASADI_JL",
    joinpath(@__DIR__, "..", "..", "..", "julia-proto", "casadi-contact"))
include(joinpath(_jl_dir, "CasADiNative.jl"))
import .CasADiNative as ca

function polymul(a, b)
    r = zeros(length(a) + length(b) - 1)
    for i in 1:length(a), j in 1:length(b)
        r[i + j - 1] += a[i] * b[j]
    end
    r
end
polyval(p, x) = foldl((r, c) -> r * x + c, p; init = 0.0)
polyder(p) = (n = length(p) - 1; n == 0 ? [0.0] : [p[i] * (n - i + 1) for i in 1:n])
polyint(p) = vcat([p[i] / (length(p) - i + 1) for i in 1:length(p)], 0.0)

function main()
    # Degree of interpolating polynomial
    d = 3

    # Choose collocation points (Radau)
    tau_root = [0.0; ca.collocation_points(d, "radau")]

    # Coefficients of the collocation, continuity and quadrature equations
    Cc = zeros(d + 1, d + 1)
    D = zeros(d + 1)
    F = zeros(d + 1)

    # Construct polynomial basis
    for j in 1:(d + 1)
        p = [1.0]
        for r in 1:(d + 1)
            if r != j
                p = polymul(p, [1.0, -tau_root[r]]) ./ (tau_root[j] - tau_root[r])
            end
        end
        D[j] = polyval(p, 1.0)
        pder = polyder(p)
        for r in 1:(d + 1)
            Cc[j, r] = polyval(pder, tau_root[r])
        end
        F[j] = polyval(polyint(p), 1.0)
    end

    # Control discretization
    nk = 20
    tf = 10.0          # End time
    h = tf / nk        # Size of the finite elements

    # Declare variables (use scalar graph)
    t = ca.SX.sym("t")   # time
    u = ca.SX.sym("u")   # control
    x = ca.SX.sym("x", 2)  # state

    # ODE rhs function and quadratures
    xdot = ca.vertcat((1 - x[2] * x[2]) * x[1] - x[2] + u, x[1])
    qdot = x[1] * x[1] + x[2] * x[2] + u * u
    f = ca.Function("f", [t, x, u], [xdot, qdot])

    # Control bounds
    u_min, u_max, u_init = -0.75, 1.0, 0.0

    # State bounds and initial guess
    x_min = [-Inf, -Inf]; x_max = [Inf, Inf]
    xi_min = [0.0, 1.0];  xi_max = [0.0, 1.0]
    xf_min = [0.0, 0.0];  xf_max = [0.0, 0.0]
    x_init = [0.0, 0.0]

    # Dimensions
    nx, nu = 2, 1

    # Total number of variables
    NX = nk * (d + 1) * nx   # Collocated states
    NU = nk * nu             # Parametrized controls
    NXF = nx                 # Final state
    NV = NX + NU + NXF

    # NLP variable vector and bounds/guess
    V = ca.MX.sym("V", NV)
    vars_lb = zeros(NV); vars_ub = zeros(NV); vars_init = zeros(NV)
    offset = 0

    # Get collocated states (as 1-based [k][j]) and parametrized control [k]
    Xc = [Vector{Any}(undef, d + 1) for _ in 1:(nk + 1)]
    U = Vector{Any}(undef, nk)
    for k in 1:nk
        # Collocated states
        for j in 1:(d + 1)
            Xc[k][j] = V[(offset + 1):(offset + nx)]
            vars_init[(offset + 1):(offset + nx)] = x_init
            if k == 1 && j == 1
                vars_lb[(offset + 1):(offset + nx)] = xi_min
                vars_ub[(offset + 1):(offset + nx)] = xi_max
            else
                vars_lb[(offset + 1):(offset + nx)] = x_min
                vars_ub[(offset + 1):(offset + nx)] = x_max
            end
            offset += nx
        end
        # Parametrized controls
        U[k] = V[(offset + 1):(offset + nu)]
        vars_lb[offset + 1] = u_min; vars_ub[offset + 1] = u_max
        vars_init[offset + 1] = u_init
        offset += nu
    end
    # State at end time
    Xc[nk + 1][1] = V[(offset + 1):(offset + nx)]
    vars_lb[(offset + 1):(offset + nx)] = xf_min
    vars_ub[(offset + 1):(offset + nx)] = xf_max
    vars_init[(offset + 1):(offset + nx)] = x_init

    # Constraint function for the NLP
    g = ca.MX[]; J = ca.MX(ca.DM(0.0))

    # For all finite elements
    for k in 1:nk
        # For all collocation points
        for j in 2:(d + 1)
            # State derivative at the collocation point
            xp_jk = ca.MX(ca.DM(zeros(nx)))
            for r in 1:(d + 1)
                xp_jk = xp_jk + Cc[r, j] * Xc[k][r]
            end
            # Collocation equations
            Tkj = h * (k - 1 + tau_root[j])
            fk, qk = f(Tkj, Xc[k][j], U[k])
            push!(g, h * fk - xp_jk)
            # Contribution to objective
            J = J + F[j] * qk * h
        end
        # State at the end of the finite element
        xf_k = ca.MX(ca.DM(zeros(nx)))
        for r in 1:(d + 1)
            xf_k = xf_k + D[r] * Xc[k][r]
        end
        # Continuity equation
        push!(g, Xc[k + 1][1] - xf_k)
    end

    # NLP
    nlp = Dict("x" => V, "f" => J, "g" => ca.vertcat(g))

    # Allocate an NLP solver (mumps default; the python ma27 needs HSL)
    solver = ca.nlpsol("solver", "ipopt", nlp, Dict("expand" => true))

    # Solve (constraints are all equalities -> lbg = ubg = 0)
    res = solver(x0 = vars_init, lbx = vars_lb, ubx = vars_ub, lbg = 0, ubg = 0)

    println("optimal cost: ", Float64(res["f"]))

    # States at the start of each finite element (stride (d+1)*nx + nu)
    stride = (d + 1) * nx + nu
    x0_opt = res["x"][1:stride:end]
    x1_opt = res["x"][2:stride:end]
    println("x[0] at element starts = ", x0_opt)
    println("x[1] at element starts = ", x1_opt)
end

main()
