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
# Indirect (Pontryagin) method for the Van der Pol OCP: form the Hamiltonian,
# eliminate the control by its analytic minimizer, integrate the augmented
# state/costate DAE, and solve the two-point boundary value problem with a
# multiple-shooting rootfinder (nlpsol/ipopt). "rk" stands in for "cvodes".
const _jl_dir = get(ENV, "CASADI_JL",
    joinpath(@__DIR__, "..", "..", "..", "julia-proto", "casadi-contact"))
include(joinpath(_jl_dir, "CasADiNative.jl"))
import .CasADiNative as ca

grad(e, a) = ca.gradient(e, a)

# Declare variables (use a simple, efficient DAG)
x0 = ca.SX.sym("x0"); x1 = ca.SX.sym("x1")
x = ca.vertcat(x0, x1)

# Control
u = ca.SX.sym("u")

# ODE right hand side
xdot = ca.vertcat((1 - x1 * x1) * x0 - x1 + u, x0)

# Lagrangian
L = x0 * x0 + x1 * x1 + u * u

# Costate
lam = ca.SX.sym("lam", 2)

# Hamiltonian
H = ca.dot(lam, xdot) + L

# Costate equations
ldot = -grad(H, x)

println("Hamiltonian: ", H)

# H is convex quadratic in u: H = u*u + p*u + q; extract p
p = grad(H, u)                     # gives 2*u + p
p = ca.substitute(p, u, ca.SX(0.0))  # replace u with zero: gives p

# Unconstrained minimizer u = -p/2, clipped to [-0.75, 1.0]
u_opt = -0.5 * p
u_opt = ca.fmin(u_opt, ca.SX(1.0))
u_opt = ca.fmax(u_opt, ca.SX(-0.75))
println("optimal control: ", u_opt)

# Augment f with ldot and substitute the optimal control
f = ca.vertcat(xdot, ldot)
f = ca.substitute(f, u, u_opt)

# Function for the optimal control given the augmented state
xlam = ca.vertcat(x, lam)
u_fcn = ca.Function("ufcn", [xlam], [u_opt])

# DAE for the augmented dynamics
dae = Dict("x" => xlam, "ode" => f)

nX = 4              # augmented state dimension
tf = 10.0           # end time
num_nodes = 20      # shooting nodes

# Integrator over one shooting interval
I = ca.integrator("I", "rk", dae, 0.0, tf / num_nodes,
                Dict("number_of_finite_elements" => 50))

function main()
    # States at each shooting node (4 x (num_nodes+1))
    X = ca.MX.sym("X", nX, num_nodes + 1)

    # Formulate the root finding problem
    G = ca.MX[]
    push!(G, X[1:2, 1] - ca.MX(ca.DM([0.0, 1.0])))   # states fixed, costates free at t=0
    for k in 1:num_nodes
        XF = I(x0 = X[:, k])["xf"]
        push!(G, XF - X[:, k + 1])
    end
    push!(G, X[3:4, num_nodes + 1] - ca.MX(ca.DM([0.0, 0.0])))  # costates fixed at t=tf

    # Root-finding problem (vectorize X), solved as an NLP feasibility problem
    rfp = ca.Function("rfp", [ca.vec(X)], [ca.vertcat(G)])
    solver = ca.rootfinder("solver", "nlpsol", rfp, Dict(
        "nlpsol" => "ipopt",
        "nlpsol_options" => Dict{String,Any}(
            "ipopt.hessian_approximation" => "limited-memory",
            "ipopt.print_level" => 0, "print_time" => false)))

    # Solve
    X_sol = solver(ca.zeros(ca.GenDM, nX * (num_nodes + 1), 1))
    println("-----")
    println("node-state solution (4 x ", num_nodes + 1, "):")
    println("X_sol = ", X_sol)

    # Re-integrate each shooting interval to recover the trajectories. A single
    # 10s re-integration would blow up on the unstable costate ODE.
    per = 5      # sub-samples per shooting interval
    tgrid_node = collect(range(0, tf / num_nodes, length = per + 1))
    vis = ca.integrator("vis", "rk", dae, 0.0, tgrid_node[2:end],
                     Dict("number_of_finite_elements" => 20))
    segs = ca.DM[]
    for k in 0:num_nodes - 1
        xk = X_sol[k * nX + 1:(k + 1) * nX]
        push!(segs, vis(x0 = xk)["xf"])
    end
    sol = ca.horzcat(segs)                                  # 4 x (num_nodes*per)
    u_traj = u_fcn(sol)                                  # columnwise control

    println("x trajectory = ", sol[1, :])
    println("y trajectory = ", sol[2, :])
    println("u trajectory = ", u_traj)
end

main()
