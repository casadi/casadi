# -*- coding: utf-8 -*-


success = False
try:
  import pydot
  success = True
except:
  pass

if success:
  from graph import dotgraph, dotdraw, dotsave
