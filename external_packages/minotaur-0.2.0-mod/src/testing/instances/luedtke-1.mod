#
# This is the simple "Luedtke Generalized Pooling Problem"
#  from our MINLP paper
#

var x1 >= -1, <= 4;
var x2 >= -2.4, <= 4;
var x3 >= -0.5, <= 6;
var x4 >= -0.5, <= 6;
var x5 >= 1, <= 6;
var x6 >= 1, <= 5;
var x7 >= 1, <= 3;
var x8 >= 1, <= 3;

minimize Obj: 4*(x1+x2+x3) + 3*(x4+x5) + 3.5*(x6+x7) + 2.5*x8 ;

subject to
  c1: x1*x2*x3*x4 - x1*x2 - x4*x5 + x5 + x6 >= 230 ;

  c2: x3*x4*x5*x6 - x1*x4 - x6*x7 + x2 + x8 = -2 ;




