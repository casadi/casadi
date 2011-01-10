from casadi import *
matrix_double.getitem_old = matrix_double.__getitem__

def getitem_new(self,s):
  if isinstance(s,tuple) and len(s)==2 and isinstance(s[0],slice) and isinstance(s[1],slice):
    return (s[0].indices(self.size1()), s[1].indices(self.size2()))
  else:
    return self.getitem_old(s)
  
matrix_double.__getitem__ = getitem_new
  

a = matrix_double(2,3,1)
a.resize(10,10)
print a
print list(a[0:2])
print a[0,2]
print a[0:2,0:2]