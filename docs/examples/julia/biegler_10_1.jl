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
# Exercise 1, chapter 10 from Larry Biegler's book: collocation solution of the
# scalar ODE dz/dt = z^2 - 2z + 1 with z(0) = -3, swept over N = 1..10 finite
# elements. Lagrange-polynomial collocation coefficients are built symbolically.
const _jl_dir = get(ENV, "CASADI_JL",
    joinpath(@__DIR__, "..", "..", "..", "julia-proto", "casadi-contact"))
include(joinpath(_jl_dir, "CasADiNative.jl"))
import .CasADiNative as ca

function main()
    println("program started")

    # Test with different number of elements
    for N in 1:10
        println("N = ", N)

        # Degree of interpolating polynomial
        K = 2

        # Legendre roots
        tau_root = [0.0, 0.211325, 0.788675]

        # Differential equation
        z = ca.SX.sym("z")
        F = ca.Function("dz_dt", [z], [z * z - 2 * z + 1])

        z0 = -3

        # Collocation point and step size
        tau = ca.SX.sym("tau")
        h = 1.0 / N

        # Coefficients of the continuity and collocation equations
        D = zeros(K + 1)
        Cc = zeros(K + 1, K + 1)
        for j in 1:(K + 1)
            # Lagrange polynomial
            L = ca.SX(ca.DM(1.0))
            for k in 1:(K + 1)
                if k != j
                    L = L * (tau - tau_root[k]) / (tau_root[j] - tau_root[k])
                end
            end
            # Evaluate at end for continuity-equation coefficient
            lfcn = ca.Function("lfcn", [tau], [L])
            D[j] = Float64(lfcn(1.0))
            # Differentiate and evaluate at collocation points
            tfcn = ca.Function("tfcn", [tau], [ca.tangent(L, tau)])
            for k in 1:(K + 1)
                Cc[j, k] = Float64(tfcn(tau_root[k]))
            end
        end

        # Collocated states: an N x (K+1) symbolic matrix
        Z = ca.SX.sym("Z", N, K + 1)

        # Construct the NLP: x = vec(Z.T) (column-major over transposed rows)
        x = ca.vec(Z')
        g = ca.SX[]
        for i in 1:N
            for k in 2:(K + 1)
                # Collocation equations
                rhs = ca.SX(ca.DM(0.0))
                for j in 1:(K + 1)
                    rhs = rhs + Z[i, j] * Cc[j, k]
                end
                FF = F(Z[i, k])
                push!(g, h * FF - rhs)
            end
            # Continuity equation
            rhs = ca.SX(ca.DM(0.0))
            for j in 1:(K + 1)
                rhs = rhs + D[j] * Z[i, j]
            end
            if i < N
                push!(g, Z[i + 1, 1] - rhs)
            end
        end
        gg = ca.vertcat(g)

        # NLP
        nlp = Dict("x" => x, "f" => x[1]^2, "g" => gg)

        # Allocate an NLP solver
        solver = ca.nlpsol("solver", "ipopt", nlp, Dict("ipopt.tol" => 1e-10))

        # Bounds: fix x[1] = z0, free the rest in [-100, 100]
        nnz = ca.nnz(x)
        lbx = fill(-100.0, nnz); ubx = fill(100.0, nnz)
        lbx[1] = ubx[1] = z0

        # Solve the problem
        res = solver(x0 = zeros(nnz), lbx = lbx, ubx = ubx, lbg = 0, ubg = 0)

        println("optimal cost: ", Float64(res["f"]))
        println("optimal solution: ", ca.nonzeros(res["x"]))
    end
end

main()
