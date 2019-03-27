A = [ 2  0  0 -1  2;
      0  4  0 -5 -1;
      0  0  6  0  4;
      1  5  0  0  0;
      2 -1  4  0  0];
 
b = [0; 0; 0; 12; 18];

x = A\b

f = x(1)^2 + 2*x(2)^2 + 3*x(3)^2



% The above shows that lambda2 is negative, indicating that it is not a
% binding constraint. Running again without g2(x):

A = [2 0 0 -1;
     0 4 0 -5;
     0 0 6  0;
     1 5 0  0];

b = [0; 0; 0; 12];

x = A\b

f = x(1)^2 + 2*x(2)^2 + 3*x(3)^2