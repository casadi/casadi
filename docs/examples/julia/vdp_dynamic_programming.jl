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
# Dynamic programming over a discretized state/control grid for the Van der Pol
# problem. The Python original is pure numpy (no CasADi symbolics); this port
# keeps it as plain Julia array math. matplotlib output is dropped; we log the
# optimal cost and trajectory.
const _jl_dir = get(ENV, "CASADI_JL",
    joinpath(@__DIR__, "..", "..", "..", "julia-proto", "casadi-contact"))
include(joinpath(_jl_dir, "CasADiNative.jl"))
import .CasADiNative as ca

function main()
    T = 10.0     # End time
    N = 20       # Number of control intervals
    NK = 20      # RK4 steps per interval
    DT = T / (N * NK)
    NU = 101     # Number of discrete control values
    NX = 101     # Number of discrete state values per axis

    # System dynamics (x1_dot, x2_dot, q_dot)
    f(x1, x2, u) = ((1 - x2 * x2) * x1 - x2 + u, x1, x1 * x1 + x2 * x2 + u * u)

    U = range(-1, 1, length = NU)
    x1 = range(-1, 1, length = NX)
    x2 = range(-1, 1, length = NX)

    NG = NX * NX                         # total grid points (row = x2, col = x1)
    idx(i2, i1) = i2 * NX + i1 + 1       # 0-based (i2, i1) -> 1-based linear index

    # For each control action, precompute next-state index and stage cost
    stage_J = Vector{Vector{Float64}}(undef, NU)
    next_x1 = Vector{Vector{Int}}(undef, NU)
    next_x2 = Vector{Vector{Int}}(undef, NU)
    for uind in 1:NU
        u = U[uind]
        nx1 = zeros(Int, NG); nx2 = zeros(Int, NG); Q = zeros(NG)
        for i2 in 0:(NX - 1)
            for i1 in 0:(NX - 1)
                X1 = x1[i1 + 1]; X2 = x2[i2 + 1]; Qk = 0.0
                for _ in 1:NK
                    a1, a2, aq = f(X1, X2, u)
                    b1, b2, bq = f(X1 + DT / 2 * a1, X2 + DT / 2 * a2, u)
                    c1, c2, cq = f(X1 + DT / 2 * b1, X2 + DT / 2 * b2, u)
                    d1, d2, dq = f(X1 + DT * c1, X2 + DT * c2, u)
                    X1 += DT / 6 * (a1 + 2b1 + 2c1 + d1)
                    X2 += DT / 6 * (a2 + 2b2 + 2c2 + d2)
                    Qk += DT / 6 * (aq + 2bq + 2cq + dq)
                end
                # Round to nearest grid index (0-based)
                r1 = round(Int, (X1 + 1) / 2 * (NX - 1))
                r2 = round(Int, (X2 + 1) / 2 * (NX - 1))
                # Infinite cost if out-of-bounds
                if r1 < 0 || r1 >= NX || r2 < 0 || r2 >= NX
                    Qk = Inf; r1 = 0; r2 = 0
                end
                gg = idx(i2, i1)
                nx1[gg] = r1; nx2[gg] = r2; Q[gg] = Qk
            end
        end
        next_x1[uind] = nx1; next_x2[uind] = nx2; stage_J[uind] = Q
    end

    # Cost-to-go (no end cost) and optimal control
    J = zeros(NG)
    U_opt = Vector{Vector{Int}}()
    for _ in 1:N
        J_prev = fill(Inf, NG)
        u_prev = fill(-1, NG)
        for uind in 1:NU
            nx1 = next_x1[uind]; nx2 = next_x2[uind]; sj = stage_J[uind]
            for gg in 1:NG
                test = J[idx(nx2[gg], nx1[gg])] + sj[gg]
                if test < J_prev[gg]
                    J_prev[gg] = test; u_prev[gg] = uind
                end
            end
        end
        J = J_prev
        push!(U_opt, u_prev)
    end
    reverse!(U_opt)

    # Optimal control starting at x1 = 0, x2 = 1 (0-based grid indices)
    i1 = NX ÷ 2
    i2 = NX - 1
    u_opt = Float64[]; x1_opt = [x1[i1 + 1]]; x2_opt = [x2[i2 + 1]]
    cost = 0.0
    for k in 1:N
        u_ind = U_opt[k][idx(i2, i1)]
        cost += stage_J[u_ind][idx(i2, i1)]
        ni1 = next_x1[u_ind][idx(i2, i1)]; ni2 = next_x2[u_ind][idx(i2, i1)]
        i1 = ni1; i2 = ni2
        push!(u_opt, U[u_ind]); push!(x1_opt, x1[i1 + 1]); push!(x2_opt, x2[i2 + 1])
    end

    println("-----")
    println("Minimal cost: ", cost)
    # Consistency check (cf. python assert)
    @assert abs(cost - J[idx(NX - 1, NX ÷ 2)]) < 1e-8
    println("x1 trajectory = ", x1_opt)
    println("x2 trajectory = ", x2_opt)
    println("u  trajectory = ", u_opt)
end

main()
