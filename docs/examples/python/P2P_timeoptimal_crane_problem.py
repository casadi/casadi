"""
This file implements the P2P motion problem of the overhead crane in casadi.

In rockit the slack reformulation did not work immediately as I wanted it to be
"""
import casadi as cs
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.patches import Rectangle
import copy


class crane_problem:

    def __init__(self, N=20):
        """
        The constructor for the optimization problem.

        Args:
            N (int, optional): The number of discretization intervals.
                               Defaults to 20.
        """
        self.N = N
        self.g = 9.81  # m/s**2

        # Set number of entries per decision variable
        self.n_x = 6
        self.n_u = 2
        self.n_sh = 3
        self.n_t = 1
        self.n_s0 = self.n_x
        self.n_sf = self.n_x

        # Some Parameters
        self.crane_angle_limit = 0.75  # in radians
        self.cart_pos_min = -0.1  # in metres
        self.cart_pos_max = 0.6  # in metres
        self.cart_max_vel = 0.4  # in m/s
        self.l_min = 1e-2
        self.l_max = 2
        self.l_dot_max = 0.25
        self.l_ddot_max = 5
        self.x_ddot_max = 5

        # Initial states
        self.x_0 = 0.0
        self.x_f = 0.5
        self.l_0 = 0.9
        self.l_f = 0.9
        self.theta_0 = 0.0
        self.theta_f = 0.0
        self.theta_dot_0 = 0
        self.theta_dot_f = 0
        self.x_dot_0 = 0
        self.x_dot_f = 0
        self.l_dot_0 = 0
        self.l_dot_f = 0

        self.rad = 0.08
        self.x_obstacle = [0.1, 0.2]
        self.y_obstacle = [-2.0, -0.7]

        self.lower_path_bounds = cs.vertcat(
            self.l_min,
            self.cart_pos_min,
            -self.crane_angle_limit,
            -self.l_dot_max,
            -self.cart_max_vel)
        self.upper_path_bounds = cs.vertcat(
            self.l_max,
            self.cart_pos_max,
            self.crane_angle_limit,
            self.l_dot_max,
            self.cart_max_vel)

        self.indeces_S0 = None
        self.indeces_Sf = None
        self.indeces_state0 = None
        self.indeces_statef = None

        self.ode_integrator = self.create_integrator_time_optimal()

    def create_integrator_time_optimal(self):
        """
        Returns the time optimal integrator.
        """
        # Define states
        length = cs.MX.sym('l', 1)  # length of the pendulum
        x_c = cs.MX.sym('x_c', 1)
        theta = cs.MX.sym('theta', 1)
        l_dot = cs.MX.sym('l_dot', 1)  # derivative of length of the pendulum
        x_c_dot = cs.MX.sym('x_c_dot', 1)
        theta_dot = cs.MX.sym('theta_dot', 1)
        T = cs.MX.sym('T')
        x_temp = cs.vertcat(length, x_c, theta, l_dot, x_c_dot, theta_dot)

        # Define controls
        x_ddot = cs.MX.sym('x_ddot', 1)
        l_ddot = cs.MX.sym('l_ddot', 1)
        u_temp = cs.vertcat(l_ddot, x_ddot)
        theta_ddot = (-x_ddot*cs.cos(theta) - 2*l_dot *
                      theta_dot - self.g*cs.sin(theta))/length
        xdot_temp = cs.vertcat(l_dot, x_c_dot, theta_dot,
                               l_ddot, x_ddot, theta_ddot)
        ode = {'x': x_temp, 'p': cs.vertcat(u_temp, T), 'ode': T*xdot_temp}
        opts = {'tf': 1}
        ode_integrator = cs.integrator('integrator', 'rk', ode, opts)
        return ode_integrator

    def create_problem(self, dev_lx_0=[0.0, 0.0], dev_lx_f=[0.0, 0.0]):
        """
        We create the optimization problem here.

        Args:
            dev_lx_0 (list, optional): Perturbation of starting point.
                                       Defaults to [0.0, 0.0].
            dev_lx_f (list, optional): Perturbation of end point. 
                                        Defaults to [0.0, 0.0].

        Returns:
            (Casadi MX.sym arrays, Casadi DM vectors): symbolic expression for 
            decision variable, objective function, constraint function,
            Casadi DM vectors for lower and upper bounds of variables,
            the initial guess for the optimization.
        """
        ## ----- DEFINE OPTIMIZATION PROBLEM -----
        # Define optimization problem
        self.opti = cs.Opti()

        # Define optimization variables
        n_var = (self.N+1)*(self.n_x) + self.N*(self.n_u+self.n_sh) +\
            self.n_s0 + self.n_sf + self.n_t
        X_tmp = []
        U_tmp = []
        sh_tmp = []
        
        S0 = self.opti.variable(self.n_s0, 1)
        for k in range(self.N+1):
            X_tmp.append(self.opti.variable(self.n_x, 1))
            if k == self.N+1:
                self.indeces_statef = list(
                    range(self.opti.debug.x.size[1]-self.n_x,
                          self.opti.debug.x.size[1]))
            if k < self.N:
                U_tmp.append(self.opti.variable(self.n_u, 1))
                sh_tmp.append(self.opti.variable(self.n_sh, 1))    
        Sf = self.opti.variable(self.n_sf, 1)
        T = self.opti.variable(self.n_t, 1)
        
        self.indeces_S0 = list(range(self.n_s0))
        self.indeces_state0 = list(range(self.n_s0, self.n_s0+self.n_x))
        self.indeces_Sf = list(range(n_var-self.n_sf-self.n_t, n_var-self.n_t))
        self.indeces_statef = list(range(n_var -self.n_sf -self.n_t -self.n_x, n_var -self.n_sf -self.n_t))

        # Transform list into variables
        X = cs.horzcat(*X_tmp)
        U = cs.horzcat(*U_tmp)
        sh = cs.horzcat(*sh_tmp)
        # p_load starts a the first index, since we do not need it for the 
        # zero index
        p_load = cs.vertcat(X[1, 1:] + X[0, 1:] *
                            cs.sin(X[2, 1:]), -X[0, 1:]*cs.cos(X[2, 1:]))

        # Define the constraints
        # Constraints on initial states with slacked initial state
        self.start_state = cs.vertcat(dev_lx_0[0] + self.l_0,
                        dev_lx_0[1] + self.x_0,
                        self.theta_0,
                        self.l_dot_0,
                        self.x_dot_0,
                        self.theta_dot_0)

        self.end_state = cs.vertcat(dev_lx_f[0] + self.l_f,
                               dev_lx_f[1] + self.x_f,
                               self.theta_f,
                               self.l_dot_f,
                               self.x_dot_f,
                               self.theta_dot_f)

        # Shooting constraints
        for k in range(self.N+1):
            # Gap closing constraints / Multiple Shooting
            if k < self.N:
                self.opti.subject_to(self.ode_integrator(
                    x0=X[:, k],
                    p=cs.vertcat(U[:, k], T/self.N))['xf'] == X[:, k+1])
                # self.opti.subject_to(T[k] == T[k+1])

            if k == 0:
                self.opti.subject_to(self.opti.bounded(0, S0, cs.inf))
                self.opti.subject_to(self.start_state <= X[:, 0] + S0)
                self.opti.subject_to(X[:, 0] - S0 <= self.start_state)

            if k > 0:
                # Path dynamic constraints
                self.opti.subject_to(X[0:5, k] <= self.upper_path_bounds)
                self.opti.subject_to(self.lower_path_bounds <= X[0:5, k])

                # Obstacle Avoiding Constraint
                # adding separating hyperplane constraint for crane
                self.opti.subject_to(sh[0, k-1]*p_load[0, k-1] +
                                sh[1, k-1]*p_load[1, k-1] +
                                sh[2, k-1] <= -self.rad)
                # adding separating hyperplane constraint for obstacle
                self.opti.subject_to(sh[0, k-1]*self.x_obstacle[0] +
                                sh[1, k-1]*self.y_obstacle[0] +
                                sh[2, k-1] >= 0)
                self.opti.subject_to(sh[0, k-1]*self.x_obstacle[0] +
                                sh[1, k-1]*self.y_obstacle[0] +
                                sh[2, k-1] >= 0)
                self.opti.subject_to(sh[0, k-1]*self.x_obstacle[1] +
                                sh[1, k-1]*self.y_obstacle[1] +
                                sh[2, k-1] >= 0)
                self.opti.subject_to(sh[0, k-1]*self.x_obstacle[1] +
                                sh[1, k-1]*self.y_obstacle[1] +
                                sh[2, k-1] >= 0)
                # Bound constraints on separating hyperplanes
                self.opti.subject_to(self.opti.bounded(-1, sh[:, k-1], 1))

            # Constraints on controls
            if k < self.N:
                self.opti.subject_to(self.opti.bounded(-self.l_ddot_max,
                                                U[0, k], self.l_ddot_max))
                self.opti.subject_to(self.opti.bounded(-self.x_ddot_max,
                                                U[1, k], self.x_ddot_max))
            # add bound constraints on time variable

        # Slacked Constraints on terminal state and temperature constraint
        self.opti.subject_to(self.end_state <= X[:, -1] + Sf)
        self.opti.subject_to(X[:, -1] - Sf <= self.end_state)
        self.opti.subject_to(self.opti.bounded(0, Sf, cs.inf))
        self.opti.subject_to(self.opti.bounded(0.5, T, 20))

        # For the Gauss-Newton Hessian
        self.V_mats = cs.vertcat((1/cs.sqrt(10))*1e-2*U[0, :].T,
            	                 1e-3*U[1,:].T,
                                 cs.sqrt(10)*1e2*S0,
                                 cs.sqrt(10)*1e2*Sf,
                                 1/cs.sqrt(self.N)*T)

        objective = 1e5*cs.sum1(Sf) + 1e5*cs.sum1(S0) + T

        self.opti.minimize(objective)

        ## ----- CREATE FEASIBLE INITIALIZATION -----
        x0_init = cs.vertcat(0.6,
                             self.x_0,
                             self.theta_0,
                             self.l_dot_0,
                             self.x_dot_0,
                             self.theta_dot_0)

        init = []

        # Define initialisation for S0
        S0_plus_init = cs.fmax(self.start_state - x0_init, 0)
        S0_minus_init = cs.fmax(x0_init - self.start_state, 0)
        S0_init = cs.fmax(S0_plus_init, S0_minus_init)
        init.append(S0_init)

        x_curr = x0_init
        sh_init0 = cs.DM([0, -1, -0.6 - self.rad - 0.01])
        sh_init = []
        self.p_load_init = cs.vertcat(x_curr[1] + x_curr[0]*cs.sin(x_curr[2]),
                                 -x_curr[0]*cs.cos(x_curr[2]))

        u_const = cs.DM([0.0, 0.1])
        t_init = cs.DM([0.125*self.N])

        X_init = x0_init

        for k in range(self.N):
            x_curr = self.ode_integrator(
                x0=x_curr, p=cs.vertcat(u_const, t_init/self.N))['xf']
            X_init = cs.horzcat(X_init, x_curr)
            sh_init = cs.horzcat(sh_init, sh_init0)
            p_load_curr = cs.vertcat(x_curr[1] + x_curr[0]*cs.sin(x_curr[2]),
                                     -x_curr[0]*cs.cos(x_curr[2]))
            self.p_load_init = cs.horzcat(self.p_load_init, p_load_curr)
            

        x_curr = x0_init
        for k in range(self.N+1):
            init.append(x_curr)

            if cs.sum1(cs.fmax(x_curr[0:5] - self.upper_path_bounds, 0)) + cs.sum1(cs.fmax(self.lower_path_bounds -x_curr[0:5], 0)) > 1e-8:
                print('here in iter', k)
            x_curr = self.ode_integrator(
                x0=x_curr, p=cs.vertcat(u_const, t_init/self.N))['xf']
            if k < self.N:
                init.append(u_const)
                init.append(sh_init0)

        Sf_plus_init = cs.fmax(self.end_state - X_init[:, -1], 0)
        Sf_minus_init = cs.fmax(X_init[:, -1] - self.end_state, 0)
        Sf_init = cs.fmax(Sf_plus_init, Sf_minus_init)
        init.append(Sf_init)
        init.append(t_init)

        self.opti.set_initial(sh, sh_init)
        self.opti.set_initial(X, X_init)
        self.opti.set_initial(Sf, Sf_init)
        self.opti.set_initial(S0, S0_init)
        self.opti.set_initial(U, cs.repmat(u_const, 1, self.N))
        self.opti.set_initial(T, t_init)
        self.opti.solver('ipopt', {'error_on_fail': False, 'ipopt': {
            'max_iter': 2000,
            'hessian_approximation': 'exact',
            'limited_memory_max_history': 5,
            'print_level': 5}})
        # sol = opti.solve_limited()
        # self.sol = self.opti.solve()

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

    def create_gn_hessian(self):
        """
        Creates the Gauss-Newton Hessian matrix and returns it as a function.
        """
        J = cs.jacobian(self.V_mats, self.opti.x)
        self.H_gn = 2*J.T @ J

        # lam_x = cs.MX.sym('lam_x', self.opti.x.shape[0])
        # lam_g = cs.MX.sym('lam_g', self.opti.g.shape[0])

        sigma = cs.MX.sym('sigma')

        H_gn_fun = cs.Function('gn_fun', [self.opti.x,
                                          self.opti.p,
                                          sigma,
                                          self.opti.lam_g], [sigma * self.H_gn])

        return H_gn_fun


    def create_scaling_matrices(self, states=True, controls=True,
            time=True, sep_hyp=False, slack0=False, slack_f=False):
        """
        Creates the scaling matrices and its inverse for the feasible SQP
        solver.
        As default we set a trust-region around the states, controls, and time.
        For additional trust-regions set the bools to TRUE.

        Args:
            states (bool, optional): Are states in trust-region? Defaults to True.
            controls (bool, optional): Are controls in trust-region?. Defaults to True.
            time (bool, optional): Is time in trust-region?. Defaults to True.
            sep_hyp (bool, optional): Are separating hyperplane variables in 
            trust-region? Defaults to False.
            slack0 (bool, optional): Is slack variable of initial condition 
            in trust-region? Defaults to False.
            slack_f (bool, optional): Is slack variable of terminal condition 
            in trust-region? Defaults to False.

        Returns:
            Casadi DM vectors: First entry of tuple is scaling matrix of the 
            decision variables. Second entry is the inverse of the scaling 
            matrix.
        """
        # Define the one or zero vectors, if variable is in trust-region or
        # not
        if states:
            vec_states = np.ones(self.n_x)
        else:
            vec_states = np.zeros(self.n_x)

        if controls:
            vec_controls = np.ones(self.n_u)
        else:
            vec_controls = np.zeros(self.n_u)

        if time:
            vec_time = np.ones(self.n_t)
        else:
            vec_time = np.zeros(self.n_t)

        if sep_hyp:
            vec_sep_hyp = np.ones(self.n_sh)
        else:
            vec_sep_hyp = np.zeros(self.n_sh)

        if slack0:
            vec_slack0 = np.ones(self.n_s0)
        else:
            vec_slack0 = np.zeros(self.n_s0)

        if slack_f:
            vec_slack_f = np.ones(self.n_sf)
        else:
            vec_slack_f = np.zeros(self.n_sf)

        # Create the trust-region scaling matrix
        list_tr_scale_mat_vec = []
        list_tr_scale_mat_vec.append(vec_slack0)
        for k in range(self.N+1):
            list_tr_scale_mat_vec.append(vec_states)
            if k < self.N:
                list_tr_scale_mat_vec.append(vec_controls)
                list_tr_scale_mat_vec.append(vec_sep_hyp)
        list_tr_scale_mat_vec.append(vec_slack_f)
        list_tr_scale_mat_vec.append(vec_time)
        tr_scale_mat_vec = np.hstack(list_tr_scale_mat_vec)
        n_nonzeros = np.count_nonzero(tr_scale_mat_vec)
        row_ind = np.arange(n_nonzeros)
        col_ind = np.where(tr_scale_mat_vec == 1)[0]
        tr_scale_mat = np.zeros((n_nonzeros, self.n_vars))
        tr_scale_mat[row_ind, col_ind] = 1
        tr_scale_mat = cs.DM(tr_scale_mat)

        # Create the inverse trust-region scaling matrix
        tr_scale_mat_inv_vec = copy.copy(tr_scale_mat_vec)
        tr_scale_mat_inv_vec[tr_scale_mat_inv_vec == 0] = np.inf
        tr_scale_mat_inv = cs.diag(tr_scale_mat_inv_vec)

        return tr_scale_mat_vec#tr_scale_mat, tr_scale_mat_inv

    def get_state_sol(self, sol):
        """
        Retrieves the states from all decision variables.
        (Use after optimization)

        Args:
            sol (Casadi DM vector): solution vector of the OCP.

        Returns:
            Casadi DM vector: vector just containing of the states.
        """
        states_sol = []
        ind_count = self.n_s0
        for k in range(self.N+1):
            states_sol.append(sol[ind_count:ind_count+self.n_x])
            ind_count += self.n_x
            if k < self.N:
                ind_count += self.n_u
                ind_count += self.n_sh
        ind_count += 1

        return cs.vertcat(*states_sol)

    def get_control_sol(self, sol):
        """
        Retrieves the controls from all decision variables.
        (Use after optimization)
        
        Args:
            sol (Casadi DM vector): solution vector of the OCP.

        Returns:
            Casadi DM vector: vector just containing the controls.
        """
        control_sol = []
        ind_count = self.n_s0
        for k in range(self.N+1):
            ind_count += self.n_x
            if k < self.N:
                control_sol.append(sol[ind_count:ind_count+self.n_u])
                ind_count += self.n_u
                ind_count += self.n_sh
        ind_count += 1

        control_sol = cs.reshape(cs.vertcat(*control_sol), self.n_u, self.N)

        return control_sol

    def get_optimal_time(self, sol):
        """
        Extracts the optimal time from all decision variables.
        (Use after optimization) 

        Args:
            sol (Casadi DM vector): solution vector of the OCP.

        Returns:
            Casadi DM vector: The optimal time of the OCP.
        """
        time = sol[-1]
        return time
        # return self.T0*cs.exp(time/(self.T0*self.beta))

    def get_initial_state(self, sol):
        """
        Extracts the iniial state from all decision variables.
        (Use after optimization) 

        Args:
            sol (Casadi DM vector): solution vector of the OCP.

        Returns:
            Casadi DM vector: The initial state.
        """
        x0 = sol[self.n_s0:self.n_s0+self.n_x]
        return x0

    def single_shooting_operator_and_new_sol(self, sol):
        """
        Performs the single shooting operator, i.e., starts at x0 and applies
        all controls to the system.
        This yields into a solution that is feasible w.r.t. the dynamics.
        
        Args:
            sol (Casadi DM vector): solution vector of the OCP.

        Returns:
            Casadi DM vector: the single shooting corrected solution.
        """
        time = self.get_optimal_time(sol)
        controls = self.get_control_sol(sol)
        x0 = self.get_initial_state(sol)
        # Get feasible state trajectory
        states = [x0]
        x_k = x0
        for i in range(self.N):
            x_k = self.ode_integrator(
                x0=x_k, p=cs.vertcat(controls[:,i], time/self.N))['xf']
            states.append(x_k)

        # Plug new states into original solution
        new_sol = sol
        ind_count = self.n_s0
        for k in range(self.N+1):
            new_sol[ind_count:ind_count+self.n_x] = states[k]
            ind_count += self.n_x
            if k < self.N:
                ind_count += self.n_u
                ind_count += self.n_sh
        ind_count += 1

        return new_sol

    def get_particular_states(self, sol):
        """
        Retrieves the particular states from the whole state vector. I.e., 
        this function gives the position, angle, length of cable, etc..

        Args:
            sol (Casadi DM vector): solution vector of the OCP.

        Returns:
            tuple of Casadi DM vector: returns the solution of the length of 
            the crane hoist, horizontal position of the crane, angle of the 
            hoist, and the optimal trajectory of the load.
        """
        states_sol = self.get_state_sol(sol)
        l_sol = states_sol[0::6]
        xc_sol = states_sol[1::6]
        theta_sol = states_sol[2::6]
        p_load = cs.horzcat(
            xc_sol + l_sol*cs.sin(theta_sol), -l_sol*cs.cos(theta_sol))

        return (l_sol, xc_sol, theta_sol, p_load)

    def plot(self, x_sol, x_sol_ipopt):
        """
        Plot the optimal trajectory of the crane load.
        
        Args:
            x_sol (Casadi DM vector): optimal solution calculated by FSLP
            x_sol_ipopt (Casadi DM vector): optimal solution calculated by IPOPT
        """
        (_, _, _, p_load) = self.get_particular_states(x_sol)
        (_, _, _, p_load_ipopt) = self.get_particular_states(x_sol_ipopt)

        _, ax = plt.subplots()
        ax.set_aspect('equal')
        rect = Rectangle((0.1, -2), 0.1, 1.3, linewidth=1,
                         edgecolor='r', facecolor='none')
        ax.add_patch(rect)
        ax.plot(p_load[:, 0], p_load[:, 1], label='FP-SQP opt')
        ax.plot(p_load_ipopt[:, 0], p_load_ipopt[:, 1], label='IPOPT')
        ax.set_ylim((-2, -0.02))
        plt.legend(loc='upper right')
        plt.xlabel('x1')
        plt.ylabel('x2')
        plt.title('Overhead Crane, P2P Motion with obstacle')
        plt.show()

testproblem = crane_problem()
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
                    'max_iter':20,
                    'max_inner_iter':50,
                    'tr_rad0': 1,
                    'feas_tol': 1e-8,
                    'hess_lag': testproblem.create_gn_hessian(),
                    'tr_scale_vector': testproblem.create_scaling_matrices()}
# opts_sqpmethod = {
#                     'print_time': False,
#                     'max_iter':1,
#                     'max_inner_iter':50,
#                     'tr_rad0': 1,
#                     'feas_tol': 1e-7,
#                     'hess_lag': testproblem.create_gn_hessian(),
#                     'tr_scale_vector': testproblem.create_scaling_matrices()}
nlp = {'x':x, 'f':f, 'g':g}
# solver = cs.nlpsol("solver", "ipopt", nlp)
solver = cs.nlpsol("S", "feasiblesqpmethod", nlp, opts_sqpmethod)

# Solve the problem
res = solver(x0  = x0, ubg = ubg, lbg = lbg, lbx = lbx, ubx=ubx)

print("Optimal cost:", res["f"])
print("Primal solution:", res["x"])

state_sol = testproblem.get_state_sol(res["x"])
print(state_sol)
testproblem.plot([testproblem.X_init, state_sol], ['init', 'FP-SQP'])
