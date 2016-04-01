##
##    MINOTAUR -- It's only 1/2 bull
##
##    (C)opyright 2009-- The MINOTAUR Team
##
##

# simple example of a minlp

var x0 integer;
var x1 binary;
var x2;
var x3 >=4, <=10;
var x4 >=0;

minimize obj: x0*x3 + x1*x2 + x4;

subject to cons0:
        x0*x0 + x1*x1 + x2*x2 = 1;

subject to cons1:
        x0*x0*x0 + x0*x0 <= 100;

subject to cons2:
        x0 + x1 - x2 >= 0;

subject to cons3:
        x0 + x1 + x2 <= 3;

subject to cons4:
        x3 + x4 <= 10;
