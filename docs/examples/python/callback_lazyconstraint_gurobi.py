# This example shows how to implement a Gurobi callback to add LazyConstraint at MIPSOL using the Casadi/Gurobi interface
# The results obtained by this script are identical to the one in the official Gurobi documentation available at:
# https://colab.research.google.com/github/Gurobi/modeling-examples/blob/master/traveling_salesman/tsp.ipynb
# At the bottom of this file there is an implementation of the Colab notebook

import json
import math
from itertools import combinations
import casadi as ca

# Read capital data (same as your original example)
try:
    capitals_json = json.load(open('capitals.json'))
except:
    import urllib.request
    url = 'https://raw.githubusercontent.com/Gurobi/modeling-examples/master/traveling_salesman/capitals.json'
    data = urllib.request.urlopen(url).read()
    capitals_json = json.loads(data)

capitals = []
coordinates = {}
for state in capitals_json:
    if state not in ['AK', 'HI']:
        capital = capitals_json[state]['capital']
        capitals.append(capital)
        coordinates[capital] = (float(capitals_json[state]['lat']), float(capitals_json[state]['long']))

# Compute pairwise distance matrix
def distance(city1, city2):
    c1 = coordinates[city1]
    c2 = coordinates[city2]
    diff = (c1[0]-c2[0], c1[1]-c2[1])
    return math.sqrt(diff[0]*diff[0]+diff[1]*diff[1])

dist = {(c1, c2): distance(c1, c2) for c1, c2 in combinations(capitals, 2)}

class SubtourEliminationCallback(ca.Callback):
    def __init__(self, name, capitals_list, distances, opts={}):
        ca.Callback.__init__(self)
        self.capitals = capitals_list
        self.var_names = list(distances.keys())
        self.n_vars = len(distances.values())
        self.construct(name, opts)

    def get_n_in(self):
        return 5  # x_solution, obj_val, obj_best, obj_bound, sol_count

    def get_n_out(self):
        return 3  # add_lazy_flag, A_lazy, b_lazy

    def get_sparsity_in(self, i):
        if i == 0:  # x_solution - will be a dense vector of all variables
            return ca.Sparsity.dense(self.n_vars, 1)
        elif i == 1:  # obj_val
            return ca.Sparsity.dense(1, 1)
        elif i == 2:  # obj_best
            return ca.Sparsity.dense(1, 1)
        elif i == 3:  # obj_bound
            return ca.Sparsity.dense(1, 1)
        elif i == 4:  # sol_count
            return ca.Sparsity.dense(1, 1)
        else:
            return ca.Sparsity(0, 0)

    def get_sparsity_out(self, i):
        if i == 0:  # add_lazy_flag
            return ca.Sparsity.dense(1, 1)
        elif i == 1:  # A_lazy - sparse matrix for constraint coefficients
            return ca.Sparsity.dense(1, self.n_vars)
        elif i == 2:  # b_lazy
            return ca.Sparsity.dense(1, 1)
        else:
            return ca.Sparsity(0, 0)

    def eval(self, arg):
        # Extract solution data
        x_solution = arg[0]  # Vector of all variable values
        obj_val = float(arg[1])
        obj_best = float(arg[2])
        obj_bound = float(arg[3])
        sol_count = int(arg[4])

        # Find which edges are selected in the current solution
        selected_edges_1 = []
        selected_edges_2 = []
        for idx, (i, j) in enumerate(self.var_names):
            if x_solution[idx] > 0.5:
                selected_edges_1.append((i, j))
                selected_edges_2.append((j, i))  # need to handle the symmetry in the variable vector
        selected_edges = selected_edges_1 + selected_edges_2
        # Find subtour
        tour = self.find_subtour(selected_edges)
        cities_in_tour = [(i, j) for i, j in combinations(tour, 2)]
        result = [1 if city in cities_in_tour else 0 for city in self.var_names]
        var_names_symm = [(j, i) for (i, j) in self.var_names]
        result_symm = [1 if city in cities_in_tour else 0 for city in var_names_symm]  # need to handle the symmetry in the variable vector
        result = ca.DM(result) + ca.DM(result_symm)
        # Check if we found a valid subtour that needs constraint
        if len(tour) < len(self.capitals):
            # Create lazy constraint to eliminate this subtour
            # Constraint: sum of edges in subtour <= len(tour) - 1
            add_lazy_flag = 1.0
            # Create coefficient matrix A and right-hand side b
            A_lazy = result.T
            b_lazy = ca.DM([len(tour) - 1])
            return [add_lazy_flag, A_lazy, b_lazy]
        else:
            # No constraint to add
            add_lazy_flag = 0.0
            A_lazy = ca.DM(1, self.n_vars)
            b_lazy = ca.DM(1, 1)
            return [add_lazy_flag, A_lazy, b_lazy]

    def find_subtour(self, edges):
        """Find shortest subtour in the selected edges"""
        unvisited = self.capitals[:]
        cycle = self.capitals[:]  # Dummy

        while unvisited:
            thiscycle = []
            neighbors = unvisited.copy()

            while neighbors:
                current = neighbors[0]
                thiscycle.append(current)
                unvisited.remove(current)
                neighbors = [j for i, j in edges if i == current and j in unvisited]

            if len(thiscycle) <= len(cycle):
                cycle = thiscycle

        return cycle

# Create the callback
callback = SubtourEliminationCallback('subtour_callback', capitals, dist)

# Create CasADi optimization problem

# Variables: is city 'i' adjacent to city 'j' on the tour?
# Create binary variables for all pairs
vars = {}
x = []
for (i, j) in dist.keys():
    x.append(ca.SX.sym(f"x_{i}_{j}"))
    vars[(i, j)] = x[-1]
    vars[(j, i)] = vars[(i, j)]  # Make symmetric variables point to same symbol

# Add objective
objective = 0
for (i, j) in dist.keys():
    objective += dist[(i, j)] * vars[(i, j)]

# Constraints: two edges incident to each city
g = []
lbg = []
ubg = []
for c in capitals:
    constraint_sum = 0
    for (i, j) in dist.keys():
        if i == c or j == c:
            constraint_sum += vars[(i, j)]
    g.append(constraint_sum - 2)  # constraint_sum == 2, so constraint_sum - 2 == 0
    lbg.append(0)
    ubg.append(0)

# Set up the solver with callback
solver = ca.qpsol(
    "myqp",
    'gurobi',
    {"f": objective, "x": ca.vertcat(*x), "g": ca.vertcat(*g)},
    {
    'discrete': [1] * len(x),
    'enable_mipsol_callback': True,
    'mipsol_callback': callback,
})

# Solve
try:
    sol = solver(lbx=0, ubx=1, lbg=lbg, ubg=ubg)
    print("Optimization successful!")
    print("Objective value:", float(sol["f"]))
except Exception as e:
    print("Optimization failed:", str(e))


x_sol = sol["x"].full().squeeze()

# Find which edges are selected in the current solution
selected_edges_1 = []
selected_edges_2 = []
for idx, (i, j) in enumerate(list(dist.keys())):
    if x_sol[idx] > 0.5:
        selected_edges_1.append((i, j))
        selected_edges_2.append((j, i))
selected_edges = selected_edges_1 + selected_edges_2
tour = callback.find_subtour(selected_edges)

assert len(tour) == len(capitals)

# Map the solution

import folium

map = folium.Map(location=[40,-95], zoom_start = 4)

points = []
for city in tour:
  points.append(coordinates[city])
points.append(points[0])

folium.PolyLine(points).add_to(map)

map.show_in_browser()


# ------------------------------------------------------------------------------------------------------------------------------------
# Solving using Gurobipy
# ------------------------------------------------------------------------------------------------------------------------------------
# Need to run `pip install gurobipy`

import gurobipy as gp
from gurobipy import GRB

m = gp.Model()

# Variables: is city 'i' adjacent to city 'j' on the tour?
vars = m.addVars(dist.keys(), obj=dist, vtype=GRB.BINARY, name='x')

# Symmetric direction: use dict.update to alias variable with new key
vars.update({(j,i):vars[i,j] for i,j in vars.keys()})

# Constraints: two edges incident to each city
cons = m.addConstrs(vars.sum(c, '*') == 2 for c in capitals)

# Callback - use lazy constraints to eliminate sub-tours

def subtourelim(model, where):
    if where == GRB.Callback.MIPSOL:
        # make a list of edges selected in the solution
        vals = model.cbGetSolution(model._vars)
        selected = gp.tuplelist((i, j) for i, j in model._vars.keys()
                             if vals[i, j] > 0.5)
        # find the shortest cycle in the selected edge list
        tour = subtour(selected)
        if len(tour) < len(capitals):
            # add subtour elimination constr. for every pair of cities in subtour
            model.cbLazy(gp.quicksum(model._vars[i, j] for i, j in combinations(tour, 2)) <= len(tour)-1)

# Given a tuplelist of edges, find the shortest subtour

def subtour(edges):
    unvisited = capitals[:]
    cycle = capitals[:] # Dummy - guaranteed to be replaced
    while unvisited:  # true if list is non-empty
        thiscycle = []
        neighbors = unvisited
        while neighbors:
            current = neighbors[0]
            thiscycle.append(current)
            unvisited.remove(current)
            neighbors = [j for i, j in edges.select(current, '*')
                         if j in unvisited]
        if len(thiscycle) <= len(cycle):
            cycle = thiscycle # New shortest subtour
    return cycle


# Solve model
m._vars = vars
m.Params.lazyConstraints = 1
m.optimize(subtourelim)
# m.optimize()

# Retrieve solution

vals = m.getAttr('x', vars)
selected = gp.tuplelist((i, j) for i, j in vals.keys() if vals[i, j] > 0.5)

tour = subtour(selected)
assert len(tour) == len(capitals)

# Map the solution

import folium

map = folium.Map(location=[40,-95], zoom_start = 4)

points = []
for city in tour:
  points.append(coordinates[city])
points.append(points[0])

folium.PolyLine(points).add_to(map)

map.show_in_browser()

m.dispose()
gp.disposeDefaultEnv()