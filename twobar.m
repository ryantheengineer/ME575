x1 = 30.;
x2 = 3.;

tic
for i=1:100
[func, grad] = twobarAD(x1,x2);
end
el = toc;
func
grad
el



function [func, grad] = twobarAD(x1, x2)
%This function computes function values and derivatives
%for the two-bar truss using AD as given by Neidinger
%SIAM Review, Vol. 52, No. 3, pp. 543-563

wdth = 60.; thik = 0.15;
dens = 0.3; modu = 30000.; load = 66.;

%make height and diameter design variables of class valder
hght = valder(x1,[1,0]);
diam = valder(x2,[0,1]);
%hght = x1;
%diam = x2;

% compute intermediate functions
leng = ((wdth/2.)^2 + hght^2)^0.5;
area = pi * diam * thik;
iovera = (diam^2 + thik^2) / 8.;

% compute functions
wght = 2. * dens * area * leng;
strs = load * leng / (2. * area * hght);
buck = (pi^2 * modu * iovera / (leng^2));
buckcon = strs - buck;
defl = load * leng^3 / (2. * modu * area * hght^2);

%return functions and variables
func = [wght.val, strs.val, buckcon.val, defl.val];
grad = [wght.der; strs.der; buckcon.der; defl.der];
grad = grad';


end

