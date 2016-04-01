

#example of an infeasible lp 

var x0 >=0; 
var x1 >=0;
var x2 >=1;

minimize obj: x0 + x1 + x2;

subject to cons0:
        x0 + x1 <= 3; 

subject to cons1:
        x0 + x2 <= 0;
 

