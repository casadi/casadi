## This is instance on page 95 in Wolsey.
##  max 4x1 - x2, 7x1 - 2 x2 <= 14, x2 <= 3, 2 x1 - 2 x2 <= 3
##

var x_1;
var x_2;

maximize obj: 4*x_1 - x_2;

subject to cons1: 7*x_1 - 2*x_2 <= 14;
subject to cons2: x_2 <= 3;
subject to cons3: 2*x_1 - 2*x_2 <= 3;

