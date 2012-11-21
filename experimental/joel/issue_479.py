from casadi import *

H = DMatrix([[1,-1],[-1,2]])
G = DMatrix([-2,-6])
A =  DMatrix([[1, 1],[-1, 2],[2, 1]])
UBA = DMatrix([2, 2, 3])
LBA = DMatrix([-inf]*3)

LBX = DMatrix([0]*2)
UBX = DMatrix([inf]*2)

options = {"convex": True, "mutol": 1e-12, "artol": 1e-12, "tol":1e-12}


solver = NLPQPSolver(H.sparsity(),A.sparsity())
for key, val in options.iteritems():
  if solver.hasOption(key):
    solver.setOption(key,val)
solver.setOption({'nlp_solver': WorhpSolver})
solver.init()

solver.input(QP_H).set(H)
solver.input(QP_G).set(G)
solver.input(QP_A).set(A)
solver.input(QP_LBX).set(LBX)
solver.input(QP_UBX).set(UBX)
solver.input(QP_LBA).set(LBA)
solver.input(QP_UBA).set(UBA)

#solver.solve()

solver.input(QP_H).set(H*4)

solver.solve()

print solver.output()
