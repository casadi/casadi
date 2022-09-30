"""
This script implements the optimization problem of a P2P motion problem
for the amazing scara robot.
This implementation describes the end effector position in the cartesian space
with a double integrator model. The start and end position are given in
cartesian space, i.e. x- and  y-coordinates are given.
The joints are calculated by inverse kinematics
"""
import casadi as cs
import numpy as np
import copy
import matplotlib
from matplotlib import pyplot as plt
from matplotlib.patches import Rectangle

def latexify():
    params = {# 'text.latex.preamble': r"\usepackage{amsmath}",
              'axes.labelsize': 10,
              'axes.titlesize': 10,
              'legend.fontsize': 10,
              'xtick.labelsize': 10,
              'ytick.labelsize': 10,
              'text.usetex': True,
              'font.family': 'serif'
              }

    matplotlib.rcParams.update(params)


latexify()


class doubleintegrator:

    def __init__(self, N=20):
        """
        The constructor
        """
        self.N = N

        self.n_x = 4
        self.n_u = 2
        self.n_s0 = self.n_x
        self.n_sf = self.n_x

        self.T_end = 5

        self.ode_integrator = self.create_integrator_time_optimal()

    def create_integrator_time_optimal(self):
        """
        Returns the double integrator.
        """
        x = cs.MX.sym('x', self.n_x)  
        u = cs.MX.sym('u', self.n_u)

        xdot_temp = cs.vertcat(x[2], x[3], u[0], u[1])

        opts = {'tf': self.T_end/self.N}
        ode = {'x': x, 'p': u, 'ode': xdot_temp}
        ode_integrator = cs.integrator('ode_integrator', 'rk', ode, opts)

        return ode_integrator

    def create_problem(self):
        """
        We create the optimization problem here.
        """

        # ----- Define optimization problem -----
        self.opti = cs.Opti()

        # ----- Define optimization variables -----
        n_var = (self.N+1)*(self.n_x) + self.N * \
            (self.n_u) + self.n_s0 + self.n_sf
        X_tmp = []
        U_tmp = []

        S0 = self.opti.variable(self.n_s0, 1)
        for k in range(self.N+1):
            X_tmp.append(self.opti.variable(self.n_x, 1))
            if k == self.N+1:
                self.indeces_statef = list(
                    range(self.opti.debug.x.size[1]-self.n_x, self.opti.debug.x.size[1]))
            if k < self.N:
                U_tmp.append(self.opti.variable(self.n_u, 1))
        # Why minus 2???? because the torques do not matter, just angle at the end
        Sf = self.opti.variable(self.n_sf, 1)

        self.indeces_S0 = list(range(self.n_s0))
        self.indeces_state0 = list(range(self.n_s0, self.n_s0+self.n_x))
        self.indeces_Sf = list(range(n_var-self.n_sf, n_var))
        self.indeces_statef = list(
            range(n_var - self.n_sf - self.n_x, n_var - self.n_sf))

        # Transform list into variables
        X = cs.horzcat(*X_tmp)
        U = cs.horzcat(*U_tmp)

        # ----- Define the initial state -----

        # Transform start and end positions from cartesian to joint space
        # Define parameters
        # Define x0, compare with matlab line 62
        self.start_state = cs.vertcat(0, 0, 0, 0)
        self.end_state = cs.vertcat(0, 10, 0, 0)

        objective = 0

        # Shooting constraints
        for k in range(self.N+1):
            x = X[0, k]
            y = X[1, k]
            x_dot = X[2, k]
            y_dot = X[3, k]

            if k < self.N:
                # Gap closing constraints
                self.opti.subject_to(self.ode_integrator(
                    x0=X[:, k], p=U[:, k])['xf'] == X[:, k+1])

            if k == 0:
                # Slacked Initial Condition
                self.opti.subject_to(self.opti.bounded(0, S0, cs.inf))
                self.opti.subject_to(self.start_state <= X[:, 0] + S0)
                self.opti.subject_to(X[:, 0] - S0 <= self.start_state)

            if k < self.N:
                # Constraints on controls
                self.control_max = 5
                self.opti.subject_to(self.opti.bounded(-self.control_max, U[0, k], self.control_max))
                self.opti.subject_to(self.opti.bounded(-self.control_max, U[1, k], self.control_max))
                objective += 1e-2*U[0, k]**2 + 1e-2*U[1, k]**2

        # Slacked Constraints on terminal state
        self.opti.subject_to(self.end_state <= X[:, -1] + Sf)
        self.opti.subject_to(X[:, -1] - Sf <= self.end_state)
        self.opti.subject_to(self.opti.bounded(0, Sf, cs.inf))

        objective += 1e8*cs.sum1(S0) + 1e5*cs.sum1(Sf)

        self.opti.minimize(objective)

        # ----- Create feasible initialization -----
        # q_start = self.inverse_kin_fun(p_start, cs.DM.zeros(2, 1))[0]
        # qu_start = q_start[iqu]
        # qv_start = q_start[iqv]
        # x0_init = cs.vertcat(qu_start, cs.DM.zeros(nqu, 1))
        x0_init = cs.vertcat(-0.5, -1 , 0, 0)

        init = []

        # ----- Create feasible initialization -----

        # Define initialisation for S0
        S0_plus_init = cs.fmax(self.start_state - x0_init, 0)
        S0_minus_init = cs.fmax(x0_init - self.start_state, 0)
        S0_init = cs.fmax(S0_plus_init, S0_minus_init)
        init.append(S0_init)

        u_const = cs.DM([0.0, 0.5])

        self.X_init = x0_init
        U_init = []
        x_curr = x0_init
        for k in range(self.N):
            init.append(x_curr)
            x_curr = self.ode_integrator(
                x0=x_curr, p=u_const)['xf']
            self.X_init = cs.horzcat(self.X_init, x_curr)
            U_init = cs.horzcat(U_init, u_const)
            # Initialize separating hyperplanes
            init.append(u_const)
        init.append(x_curr)
        Sf_plus_init = cs.fmax(self.end_state - self.X_init[:, -1], 0)
        Sf_minus_init = cs.fmax(self.X_init[:, -1] - self.end_state, 0)
        Sf_init = cs.fmax(Sf_plus_init, Sf_minus_init)
        init.append(Sf_init)

        # self.plot_trajectory([X_init], ['init'])

        self.opti.set_initial(X, self.X_init)
        self.opti.set_initial(Sf, Sf_init)
        self.opti.set_initial(S0, S0_init)
        self.opti.set_initial(U, U_init)

        self.opti.solver('ipopt', {'error_on_fail': False, 'ipopt': {
            "max_iter": 2000, 'hessian_approximation': 'exact', 'limited_memory_max_history': 5, 'print_level': 5}})

        x0 = cs.vertcat(*init)
        x = self.opti.x
        f = self.opti.f
        g = self.opti.g

        self.n_vars = self.opti.x.shape[0]
        lbg = cs.evalf(self.opti.lbg)
        ubg = cs.evalf(self.opti.ubg)
        lbx = -cs.inf*cs.DM.ones(self.n_vars)
        ubx = cs.inf*cs.DM.ones(self.n_vars)

        return (x, f, g, lbg, ubg, lbx, ubx, x0)
        

    def get_state_sol(self, sol):
        """
        Given the solution of the Scara robot. Get all the states.
        """
        states_sol = []
        ind_count = self.n_s0
        for k in range(self.N+1):
            states_sol.append(sol[ind_count:ind_count+self.n_x])
            ind_count += self.n_x
            if k < self.N:
                ind_count += self.n_u

        return cs.horzcat(*states_sol)

    def get_control_sol(self, sol):
        """
        Given the solution of the Scara robot. Get all the constraints.
        """
        control_sol = []
        ind_count = self.n_s0
        for k in range(self.N+1):
            ind_count += self.n_x
            if k < self.N:
                control_sol.append(sol[ind_count:ind_count+self.n_u])
                ind_count += self.n_u

        return cs.horzcat(*control_sol)

    def get_slack0(self, sol):
        """
        Extracts the slack variable at the beginning from all decision 
        variables. (Use after optimization) 

        Args:
            sol (Casadi DM vector): solution vector of the OCP.

        Returns:
            Casadi DM vector: The optimal time of the OCP.
        """
        slack0 = sol[:self.n_s0]
        return slack0

    def get_slackf(self, sol):
        """
        Extracts the slack variable at the end from all decision 
        variables. (Use after optimization) 

        Args:
            sol (Casadi DM vector): solution vector of the OCP.

        Returns:
            Casadi DM vector: The optimal time of the OCP.
        """
        slackf = sol[-self.n_sf+self.n_t:-self.n_t]
        return slackf

    def plot_controls(self, sol):

        controls = self.get_control_sol(sol)
        time_grid = np.linspace(0, np.array(self.T_end).squeeze(), self.N)
        
        plt.figure(figsize=(5,5))
        plt.subplot(211)
        plt.step(time_grid, np.array(controls[0,:].T).squeeze(), label='torqueOn1', linestyle='solid')
        plt.step(time_grid, self.control_max*np.ones(self.N), linestyle='solid')
        plt.step(time_grid, -self.control_max*np.ones(self.N), linestyle='solid')

        plt.grid(alpha=0.5)
        plt.legend(loc='upper right')

        plt.subplot(212)
        plt.step(time_grid, controls[1,:].T, label='torqueOn3', linestyle='solid')
        plt.step(time_grid, self.control_max*np.ones(self.N), linestyle='solid')
        plt.step(time_grid, -self.control_max*np.ones(self.N), linestyle='solid')
        
        plt.legend(loc='upper right')
        plt.grid(alpha=0.5)
        plt.tight_layout()
        plt.show()

    def plot_states(self, sol):

        states = self.get_state_sol(sol)
        time_grid = np.linspace(0, np.array(self.T_end).squeeze(), self.N+1)
        
        plt.figure(figsize=(5,5))
        plt.subplot(211)
        plt.step(time_grid, np.array(states[0,:].T).squeeze(), label='angle1', linestyle='solid')        
        plt.grid(alpha=0.5)
        plt.legend(loc='upper right')
        plt.subplot(212)
        plt.step(time_grid, states[1,:].T, label='angle3', linestyle='solid')        
        plt.legend(loc='upper right')
        plt.grid(alpha=0.5)
        plt.tight_layout()
        plt.show()

    def plot_trajectory(self, list_X, list_labels):
        """
        Plot the trajectory of some given joint trajectories. 
        The given trajectories are given as angles and torques on the joints 
        by list_X.
        """
        # Transform the angles and torques into the end effector positions
        list_x1_sol = []
        list_x2_sol = []
        list_xee_dot_sol = []
        list_yee_dot_sol = []
        for i in range(len(list_X)):
            list_x1_sol.append(list_X[i][0,:])
            list_x2_sol.append(list_X[i][1,:])
            list_xee_dot_sol.append(list_X[i][2,:])
            list_yee_dot_sol.append(list_X[i][3,:])

        fig, ax = plt.subplots()
        ax.set_aspect('equal')
        for i in range(len(list_X)):
            ax.plot(list_x1_sol[i].T, list_x2_sol[i].T, label=list_labels[i])
        ax.set_xlabel('$x_1$-axis')
        ax.set_ylabel('$x_2$-axis')
        ax.set_title('Double integrator')
        plt.xlim((-5, 5))
        plt.ylim((-1, 11))
        plt.legend()
        plt.grid(alpha=0.5)
        plt.tight_layout()
        plt.show()



testproblem = doubleintegrator()
(x, f, g, lbg, ubg, lbx, ubx, x0) = testproblem.create_problem()

# Create an NLP solver
opts_sqpmethod = {  'qpsol': 'nlpsol',
                    'qpsol_options': {  "nlpsol": "ipopt", 
                                        "verbose": True, 
                                        "print_time": False, 
                                        "nlpsol_options": {"ipopt": {   "print_level": 0, 
                                                                        "sb": "yes", 
                                                                        "fixed_variable_treatment": "make_constraint", 
                                                                        "hessian_constant": "yes", 
                                                                        "jac_c_constant": "yes", 
                                                                        "jac_d_constant": "yes",
                                                                        "tol": 1e-12, 
                                                                        "tiny_step_tol": 1e-20, 
                                                                        "mumps_scaling": 0, 
                                                                        "honor_original_bounds": "no", 
                                                                        "bound_relax_factor": 0}, 
                                                                        "print_time": False}, 
                                        "error_on_fail": False},
                    'print_time': False,
                    'max_iter':300,
                    'max_inner_iter':50,
                    'tr_rad0': 10,
                    'feas_tol': 1e-7}
nlp = {'x':x, 'f':f, 'g':g}
# solver = cs.nlpsol("solver", "ipopt", nlp)
solver = cs.nlpsol("S", "feasiblesqpmethod", nlp, opts_sqpmethod)

# Solve the problem
res = solver(x0  = x0, ubg = ubg, lbg = lbg, lbx = lbx, ubx=ubx)

print("Optimal cost:", res["f"])
print("Primal solution:", res["x"])

state_sol = testproblem.get_state_sol(res["x"])
print(state_sol)
testproblem.plot_trajectory([testproblem.X_init, state_sol], ['init', 'FP-SQP'])