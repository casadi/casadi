% Example matrix from Davis p. 39
A = sparse(11,11);
A(3,2) = 1;
A(6,1) = 1;
A(6,4) = 1;
A(7,1) = 1;
A(8,2) = 1;
A(8,5) = 1;
A(9,6) = 1;
A(10,3) = 1;
A(10,4) = 1;
A(10,6) = 1;
A(10,8) = 1;
A(11,3) = 1;
A(11,5) = 1;
A(11,7) = 1;
A(11,8) = 1;
A(11,10) = 1;

% Make symmetric
A = A+A';

% Get the elimination tree
etr = etree(A,'col')
