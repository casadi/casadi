##
##    MINOTAUR -- It's only 1/2 bull
##
##    (C)opyright 2009-- The MINOTAUR Team
##
##

# an example to test several cases of polynomial program.

var x0 integer;
var x1 binary;
var x2;
var x3 >=4, <=10;
var x4 >=0;
var x5 >=0 <= 5;

minimize obj: x0*x3 + x1*x2 + x4;

subject to cons0:
        x0*x0 + x1*x1 + x2*x2 = 1;

subject to cons1:
        x0*x0*x0 + x0*x0 <= 100;

subject to cons2:
        x0*x0*x0 - x5 = 0;

subject to cons3:
        x0^3 - x0^2 <= 101;

subject to cons4:
        3*x0^3 + 4*x0^2 <= 102;

subject to cons5:
        3*x0^3 + 7*x0*x0^2 + 4*x0^2 <= 103;

subject to cons6:
        3*x0^3 - 3*x0*x0^2 + 4*x0^2 <= 104;

subject to cons7:
        (3*x0 + x2)*x0 + (3*x0 + x2)*x1  - (3*x0 + x2)*x1 <= 105;

subject to cons8:
        (3*x0 + x2)*x0^3 + (3*x0 + 2*x2)*x1^2  <= 106;

subject to cons9:
        (3*x0 + x2)*(3*x0 + 2*x2)*x1^2  <= 107;

subject to cons10:
       (3*x0 + x2)*(3*x0 + 2*x2)*(3*x0 + x2)  <= 108;

subject to cons11:
        x0 + x1 - x2 >= 0;

subject to cons12:
        x0 + x1 + x2 <= 3;

subject to cons13:
        x3 + x4 <= 10;
