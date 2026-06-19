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
# A chain of N mass points connected by N-1 springs, hanging between two fixed
# supports. The equilibrium minimizes the potential energy subject to piecewise
# linear ground constraints -- a convex QP solved with qpoases.
#
#   minimize{y,z} Vchain(y, z)
#   subject to    z_i - 0.1*y_i >= 0.5
const _jl_dir = get(ENV, "CASADI_JL",
    joinpath(@__DIR__, "..", "..", "..", "julia-proto", "casadi-contact"))
include(joinpath(_jl_dir, "CasADiNative.jl"))
import .CasADiNative as ca

function main()
    # Constants
    N = 40
    m_i = 40.0 / N
    D_i = 70.0 * N
    g0 = 9.81
    zmin = 0.5   # ground

    # Objective, variables, bounds, constraints
    Vchain = ca.SX(ca.DM(0.0))
    x = ca.SX[]; lbx = Float64[]; ubx = Float64[]
    g = ca.SX[]; lbg = Float64[]; ubg = Float64[]

    y_prev = z_prev = nothing
    # Loop over all chain elements
    for i in 1:N
        # Create variables for the (y_i, z_i) coordinates
        y_i = ca.SX.sym("y_$i")
        z_i = ca.SX.sym("z_$i")
        push!(x, y_i, z_i)

        if i == 1
            append!(lbx, [-2.0, 1.0]); append!(ubx, [-2.0, 1.0])
        elseif i == N
            append!(lbx, [2.0, 1.0]); append!(ubx, [2.0, 1.0])
        else
            append!(lbx, [-Inf, zmin]); append!(ubx, [Inf, Inf])
        end

        # Spring potential
        if i > 1
            Vchain = Vchain + D_i / 2 * ((y_prev - y_i)^2 + (z_prev - z_i)^2)
        end

        # Gravitational potential
        Vchain = Vchain + g0 * m_i * z_i

        # Slanted ground constraint
        push!(g, z_i - 0.1 * y_i); push!(lbg, 0.5); push!(ubg, Inf)

        y_prev, z_prev = y_i, z_i
    end

    # Formulate QP
    qp = Dict("x" => ca.vertcat(x), "f" => Vchain, "g" => ca.vertcat(g))

    # Solve with qpoases
    solver = ca.qpsol("solver", "qpoases", qp, Dict("sparse" => true))
    sol = solver(lbx = lbx, ubx = ubx, lbg = lbg, ubg = ubg)

    x_opt = sol["x"]
    println("f_opt = ", sol["f"])

    # Retrieve the (y, z) coordinates
    Y0 = x_opt[1:2:end]
    Z0 = x_opt[2:2:end]
    println("y coordinates = ", Y0)
    println("z coordinates = ", Z0)
end

main()
