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
# Simulate a Van der Pol-like system over a grid of N intervals with piecewise
# constant controls, also accumulating the integral of a quadrature y.
# "rk" stands in for the python "cvodes" integrator.
const _jl_dir = get(ENV, "CASADI_JL",
    joinpath(@__DIR__, "..", "..", "..", "julia-proto", "casadi-contact"))
include(joinpath(_jl_dir, "CasADiNative.jl"))
import .CasADiNative as ca

# Nonlinear system with two states and one control input
x1 = ca.MX.sym("x1")
x2 = ca.MX.sym("x2")
x = ca.vertcat(x1, x2)
u = ca.MX.sym("u")
xdot = ca.vertcat((1 - x2^2) * x1 - x2 + u, x1)

# Sum-of-squares distance from the origin
y = x1^2 + x2^2

# Time horizon, discretization
N = 10
T = 10.0
tgrid = collect(range(0, T, length = N + 1))

# Integrate, also calculating the integral of y
dae = Dict("x" => x, "u" => u, "ode" => xdot, "quad" => y)
F = ca.integrator("F", "rk", dae, tgrid[1], tgrid[2:end], Dict("simplify" => true))

# Pointwise values of y
yfun = ca.Function("yfun", [x], [y], ["x"], ["y"])

function main()
    # Initial conditions and piecewise-constant controls (1 x N row)
    uvals = ca.DM(collect(range(-1, 1, length = N)))'
    x0 = ca.DM([0.0, 0.0])

    # Simulate
    Fk = F(x0 = x0, u = uvals)
    xf = Fk["xf"]
    qf = Fk["qf"]

    # Prepend the state at the initial time (typed horzcat keeps it numeric DM)
    xf = ca.horzcat(x0, xf)
    qf = ca.horzcat(0.0, qf)

    # y at each grid point
    yf = yfun(xf)

    println("-----")
    println("x1 = ", xf[1, :])
    println("x2 = ", xf[2, :])
    println("y  = ", yf)
    println("integral of y (quadrature) = ", qf)

    # Compare the integral with a Riemann sum of the pointwise y values
    println("integral of y (cumsum)     = ", ca.cumsum(T / N * yf))
end

main()
