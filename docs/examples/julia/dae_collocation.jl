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
# Hand-built direct collocation of a crane/pendulum DAE optimal control problem
# (Mario Zanon & Sebastien Gross, KU Leuven 2012), solved as one big NLP with
# ipopt. No integrator plugin is used -- the collocation equations are assembled
# by hand -- so this ports directly. The original ipopt option
# linear_solver='ma27' is dropped (default mumps is used).
const _jl_dir = get(ENV, "CASADI_JL",
    joinpath(@__DIR__, "..", "..", "..", "julia-proto", "casadi-contact"))
include(joinpath(_jl_dir, "CasADiNative.jl"))
import .CasADiNative as ca

function main()
    # -------- Collocation setup --------
    nicp = 1
    xref = 0.1
    l = 1.0; mmass = 1.0; Mmass = 1.0; g = 9.81
    tf = 5.0
    nk = 50
    ndstate = 6; nastate = 1; ninput = 1
    deg = 4
    h = tf / nk / nicp

    tau = ca.SX.sym("tau")
    tau_root = [0.0; ca.collocation_points(deg, "radau")]

    # Lagrange-polynomial collocation coefficients
    Cc = zeros(deg + 1, deg + 1)
    D = zeros(deg + 1)
    for j in 1:(deg + 1)
        L = ca.SX(ca.DM(1.0))
        for j2 in 1:(deg + 1)
            if j2 != j
                L = L * (tau - tau_root[j2]) / (tau_root[j] - tau_root[j2])
            end
        end
        lfcn = ca.Function("lfcn", [tau], [L])
        D[j] = Float64(lfcn(1.0))
        tfcn = ca.Function("tfcn", [tau], [ca.tangent(L, tau)])
        for j2 in 1:(deg + 1)
            Cc[j, j2] = Float64(tfcn(tau_root[j2]))
        end
    end

    # -------- Model setup (implicit DAE) --------
    t = ca.SX.sym("t")
    u = ca.SX.sym("u")
    xd = ca.SX.sym("xd", ndstate)
    xa = ca.SX.sym("xa", nastate)
    xddot = ca.SX.sym("xdot", ndstate)
    p = ca.SX.sym("p", 0, 1)

    x, y, w, dx, dy, dw = (xd[i] for i in 1:6)
    xa0 = xa[1]

    res = ca.vertcat(
        xddot[1] - dx,
        xddot[2] - dy,
        xddot[3] - dw,
        mmass * xddot[4] + (x - w) * xa0,
        mmass * xddot[5] + y * xa0 - g * mmass,
        Mmass * xddot[6] + (w - x) * xa0 + u,
        (x - w) * (xddot[4] - xddot[6]) + y * xddot[5] + dy * dy + (dx - dw) * (dx - dw))

    ffcn = ca.Function("ffcn", [t, xddot, xd, xa, u, p], [res])
    MayerTerm = ca.Function("mayer", [t, xd, xa, u, p],
        [(x - xref)^2 + (w - xref)^2 + dx * dx + dy * dy])
    LagrangeTerm = ca.Function("lagrange", [t, xd, xa, u, p],
        [(x - xref)^2 + (w - xref)^2])

    # Bounds
    u_min = [-2.0]; u_max = [2.0]
    xD_min = fill(-Inf, 6); xD_max = fill(Inf, 6)
    xDi_min = [0.0, l, 0.0, 0.0, 0.0, 0.0]; xDi_max = copy(xDi_min)
    xD_init = [0.0, l, 0.0, 0.0, 0.0, 0.0]
    xA_min = [-Inf]; xA_max = [Inf]; xA_init = [sign(l) * 9.81]

    # -------- NLP variable layout --------
    nx = ndstate + nastate
    ndiff = ndstate; nalg = nastate; nu = ninput; NP = 0
    NXD = nicp * nk * (deg + 1) * ndiff
    NXA = nicp * nk * deg * nalg
    NU = nk * nu
    NXF = ndiff
    NV = NXD + NXA + NU + NXF + NP

    V = ca.MX.sym("V", NV)
    vars_lb = zeros(NV); vars_ub = zeros(NV); vars_init = zeros(NV)
    setrange!(arr, off, vals) = (arr[(off + 1):(off + length(vals))] = vals)

    # XD[k][i][j], XA[k][i][j-1] (1-based), U[k] (slices of V)
    XD = [[Vector{Any}(undef, deg + 1) for _ in 1:nicp] for _ in 1:(nk + 1)]
    XA = [[Vector{Any}(undef, deg) for _ in 1:nicp] for _ in 1:nk]
    U = Vector{Any}(undef, nk)

    offset = 0
    for k in 1:nk
        for i in 1:nicp
            for j in 1:(deg + 1)
                XD[k][i][j] = V[(offset + 1):(offset + ndiff)]
                if j != 1
                    XA[k][i][j - 1] = V[(offset + ndiff + 1):(offset + ndiff + nalg)]
                end
                if k == 1 && j == 1 && i == 1
                    setrange!(vars_init, offset, xD_init)
                    setrange!(vars_lb, offset, xDi_min)
                    setrange!(vars_ub, offset, xDi_max)
                    offset += ndiff
                elseif j != 1
                    setrange!(vars_init, offset, vcat(xD_init, xA_init))
                    setrange!(vars_lb, offset, vcat(xD_min, xA_min))
                    setrange!(vars_ub, offset, vcat(xD_max, xA_max))
                    offset += nx
                else
                    setrange!(vars_init, offset, xD_init)
                    setrange!(vars_lb, offset, xD_min)
                    setrange!(vars_ub, offset, xD_max)
                    offset += ndiff
                end
            end
            U[k] = V[(offset + 1):(offset + nu)]
            setrange!(vars_lb, offset, u_min)
            setrange!(vars_ub, offset, u_max)
            setrange!(vars_init, offset, [0.0])
            offset += nu
        end
    end
    XD[nk + 1][1][1] = V[(offset + 1):(offset + ndiff)]
    setrange!(vars_lb, offset, xD_min)
    setrange!(vars_ub, offset, xD_max)
    setrange!(vars_init, offset, xD_init)
    offset += ndiff
    @assert offset == NV

    P = ca.MX(0, 1)

    # -------- Constraints: collocation + continuity --------
    g_con = ca.MX[]
    for k in 1:nk
        for i in 1:nicp
            for j in 2:(deg + 1)
                xp_jk = ca.MX(ca.DM(zeros(ndiff)))
                for j2 in 1:(deg + 1)
                    xp_jk = xp_jk + Cc[j2, j] * XD[k][i][j2]
                end
                fk = ffcn(ca.DM(0.0), xp_jk / h, XD[k][i][j], XA[k][i][j - 1], U[k], P)
                push!(g_con, fk[1:ndiff])               # differential part == 0
                push!(g_con, fk[(ndiff + 1):(ndiff + nalg)])  # algebraic part == 0
            end
            xf_k = ca.MX(ca.DM(zeros(ndiff)))
            for j in 1:(deg + 1)
                xf_k = xf_k + D[j] * XD[k][i][j]
            end
            if i == nicp
                push!(g_con, XD[k + 1][1][1] - xf_k)
            else
                push!(g_con, XD[k][i + 1][1] - xf_k)
            end
        end
    end

    # -------- Objective --------
    # Mayer term at the final collocation node
    Obj = MayerTerm(ca.DM(0.0), XD[nk][1][deg + 1], XA[nk][1][deg], U[nk], P)

    # Lagrange term via the collocation quadrature weights:
    #   lDotAtTauRoot = ca.T ; ldInv = inv(ca.T[2:,2:]) ; lAtOne = D[2:]
    CT = collect(transpose(Cc))
    ldInv = ca.MX(ca.inv(ca.DM(CT[2:end, 2:end])))
    lAtOne1 = ca.MX(ca.DM(D[2:end]))
    for k in 1:nk
        for i in 1:nicp
            dqParts = ca.MX[]
            for j in 2:(deg + 1)
                push!(dqParts, LagrangeTerm(ca.DM(0.0), XD[k][i][j], XA[k][i][j - 1], U[k], P))
            end
            dQs = h * ca.vertcat(dqParts)               # deg x 1
            Qs = ca.mtimes(ldInv, dQs)                # deg x 1
            Obj = Obj + ca.mtimes(Qs', lAtOne1)  # scalar
        end
    end

    # -------- Solve --------
    nlp = Dict("x" => V, "f" => Obj, "g" => ca.vertcat(g_con))
    solver = ca.nlpsol("solver", "ipopt", nlp, Dict(
                    "expand" => true, "ipopt.tol" => 1e-4,
                    "ipopt.print_level" => 0, "print_time" => false))
    sol = solver(x0 = vars_init, lbx = vars_lb, ubx = vars_ub, lbg = 0, ubg = 0)

    println("optimal cost: ", Float64(sol["f"]))

    # Chariot position x (xd[1]) at the start of each finite element: fixed stride
    stride = nicp * ((deg + 1) * ndiff + deg * nalg) + nu
    xstarts = sol["x"][1:stride:end]   # nk element starts + final state
    println("chariot position x at element starts (every 10th):")
    println("  ", xstarts[1:10:end])
    println("final x = ", xstarts[nk + 1], " (target xref = ", xref, ")")
end

main()
