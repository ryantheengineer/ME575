%% Part a
x0 = [3; 2];

gradf0 = [6; 7];

N0 = [1 0;
      0 1];

s0 = -N0*gradf0;

mags0 = sqrt(s0(1)^2 + s0(2)^2);

s0norm = s0/mags0;

H = [4 -3;
     -3 8];

%% Part b
alpha_star0 = -(transpose(gradf0)*s0norm)/(transpose(s0norm)*H*s0norm);

x1 = x0 + alpha_star0*s0norm;

%% Part c
gradf1 = [(4*x1(1) - 3*x1(2)); (-3*x1(1) + 8*x1(2))];
gamma0 = gradf1 - gradf0;
deltax0 = x1 - x0;

N1 = N0 + (1 + (transpose(gamma0)*N0*gamma0)/(transpose(deltax0)*gamma0))*...
    ((deltax0*transpose(deltax0))/(transpose(deltax0)*gamma0)) - ...
    (deltax0*transpose(gamma0)*N0 + N0*gamma0*transpose(deltax0))/...
    (transpose(deltax0)*gamma0);

s1 = -N1*gradf1;

mags1 = sqrt(s1(1)^2 + s1(2)^2);

s1norm = s1/mags1;

%% Part d
alpha_star1 = -(transpose(gradf1)*s1norm)/(transpose(s1norm)*H*s1norm);

x2 = x1 + alpha_star1*s1norm;


%% Part e
eigenvaluesN1 = eig(N1);

%% Part f
quasicheck = N1*gamma0;

%% Part g
eig([4 -3;-3 8])