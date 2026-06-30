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
# Race along a track with a speed limit: minimize lap time.
const _jl_dir = get(ENV, "CASADI_JL",
    joinpath(@__DIR__, "..", "..", "..", "julia-proto", "casadi-contact"))
include(joinpath(_jl_dir, "CasADiNative.jl"))
import .CasADiNative as ca

N = 100                        # number of control intervals
opti = ca.Opti()

X = ca.variable(opti, 2, N + 1)   # state trajectory: pos, speed
pos = X[1, :]
speed = X[2, :]
U = ca.variable(opti, 1, N)       # throttle
T = ca.variable(opti)             # final time

ca.minimize(opti, T)              # race in minimal time

f(x, u) = ca.vertcat(x[2], u - x[2])   # dx/dt = f(x, u)

dt = T / N
for k in 1:N                   # RK4 over each interval
    k1 = f(X[:, k], U[:, k])
    k2 = f(X[:, k] + dt / 2 * k1, U[:, k])
    k3 = f(X[:, k] + dt / 2 * k2, U[:, k])
    k4 = f(X[:, k] + dt * k3, U[:, k])
    x_next = X[:, k] + dt / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
    ca.subject_to(opti, X[:, k + 1] == x_next)
end

limit(p) = 1 - sin(2 * pi * p) / 2
ca.subject_to(opti, speed <= limit(pos))      # track speed limit
ca.subject_to(opti, ca.bounded(0, U, 1))         # control is limited
ca.subject_to(opti, pos[1] == 0)              # start at position 0
ca.subject_to(opti, speed[1] == 0)            # from stand-still
ca.subject_to(opti, pos[N + 1] == 1)          # finish line at position 1
ca.subject_to(opti, T >= 0)                   # time must be positive

ca.set_initial(opti, speed, 1)
ca.set_initial(opti, T, 1)

ca.solver!(opti, "ipopt", Dict("ipopt" => Dict{String,Any}("print_level" => 0), "print_time" => false))
sol = ca.solve!(opti)

println("minimal lap time: ", ca.value(sol, T))
