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
  
  nlpsol = SnoptSolver(nlp)
  #nlpsol.setOption("tol",1e-12)
  nlpsol.setOption("gather_stats",True)
  nlpsol.setOption("_feasibility_tolerance",1e-12)
  nlpsol.setOption("_optimality_tolerance",1e-12)
  nlpsol.setOption("_major_iteration_limit",3000)
  nlpsol.setOption("detect_linear",True)
  #nlpsol.setOption("max_iter",3000)
  nlpsol.init()
  
  nlpsol.setInput(1e-5,"x0")
  
  bs_ = mul(d["problem"]["Bs"][0],1e-5*DMatrix.ones(2,2))
  nlpsol.input("x0")[-bs_.size():] = vec(bs_)
  
  nlpsol.setInput(0,"lbg")
  nlpsol.setInput(0,"ubg")
  
  nlpsol.evaluate()
  
  print nlpsol.getStats()
  
  alldata.append({"f": nlpsol.output("f"),"x":nlpsol.output("x"), "stats": nlpsol.getStats()})
  
