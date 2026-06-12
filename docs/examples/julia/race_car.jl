# Julia port of docs/examples/python/race_car.py (Opti showcase).
# Race along a track with a speed limit: minimize lap time.
const _jl_dir = get(ENV, "CASADI_JL",
    joinpath(@__DIR__, "..", "..", "..", "julia-proto", "casadi-contact"))
include(joinpath(_jl_dir, "CasADiNative.jl"))
using .CasADiNative
const C = CasADiNative.C

N = 100                        # number of control intervals
opti = Opti()

X = variable(opti, 2, N + 1)   # state trajectory: pos, speed
pos = X[1, :]
speed = X[2, :]
U = variable(opti, 1, N)       # throttle
T = variable(opti)             # final time

minimize(opti, T)              # race in minimal time

f(x, u) = vertcat(x[2], u - x[2])   # dx/dt = f(x, u)

dt = T / N
for k in 1:N                   # RK4 over each interval
    k1 = f(X[:, k], U[:, k])
    k2 = f(X[:, k] + dt / 2 * k1, U[:, k])
    k3 = f(X[:, k] + dt / 2 * k2, U[:, k])
    k4 = f(X[:, k] + dt * k3, U[:, k])
    x_next = X[:, k] + dt / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
    subject_to(opti, X[:, k + 1] == x_next)
end

limit(p) = 1 - C.sin(2 * pi * p) / 2
subject_to(opti, speed <= limit(pos))      # track speed limit
subject_to(opti, bounded(0, U, 1))         # control is limited
subject_to(opti, pos[1] == 0)              # start at position 0
subject_to(opti, speed[1] == 0)            # from stand-still
subject_to(opti, pos[N + 1] == 1)          # finish line at position 1
subject_to(opti, T >= 0)                   # time must be positive

set_initial(opti, speed, 1)
set_initial(opti, T, 1)

solver!(opti, "ipopt"; ipopt = Dict{String,Any}("print_level" => 0), print_time = false)
sol = solve!(opti)

println("minimal lap time: ", value(sol, T))
