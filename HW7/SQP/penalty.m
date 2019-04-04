function [P] = penalty(f,lamda,g)
    P = f + lamda*abs(g);
end