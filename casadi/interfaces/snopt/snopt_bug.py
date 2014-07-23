from instances import *
from read_instances import *

from lockfile import FileLock

import cPickle as pickle

  
data = readall()
  
alldata = []
  
for d in data[:1]:

  H = d["H"]
  nlp = createNLP(DMatrix([[1.2,0.3],[0.7,1.3]]),DMatrix([[0.2,0.4],[0.77,0.12]],lift=True,simple=False)
  
  log = []
  dists = []
  
  nlpsolver = SnoptSolver(nlp)
  #nlpsolver.setOption("tol",1e-12)
  nlpsolver.setOption("gather_stats",True)
  nlpsolver.setOption("_feasibility_tolerance",1e-12)
  nlpsolver.setOption("_optimality_tolerance",1e-12)
  nlpsolver.setOption("_major_iteration_limit",3000)
  nlpsolver.setOption("detect_linear",True)
  #nlpsolver.setOption("max_iter",3000)
  nlpsolver.init()
  
  nlpsolver.setInput(1e-5,"x0")
  
  bs_ = mul(d["problem"]["Bs"][0],1e-5*DMatrix.ones(2,2))
  nlpsolver.input("x0")[-bs_.size():] = vec(bs_)
  
  nlpsolver.setInput(0,"lbg")
  nlpsolver.setInput(0,"ubg")
  
  nlpsolver.evaluate()
  
  print nlpsolver.getStats()
  
  alldata.append({"f": nlpsolver.output("f"),"x":nlpsolver.output("x"), "stats": nlpsolver.getStats()})
  
