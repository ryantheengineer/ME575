% This program constructs a contour plot for HW2
clear; close all;
L = 15*5280;        % length of pipeline, feet
W = 12.67;          % flowrate of limestone, lbm/s
a = 0.01;           % avg lump size of limestone before grinding, ft
g = 32.17;          % acceleration due to gravity, ft/s^2
rho_w = 62.4;       % density of water, lbm/ft^3
gamma = 168.5;      % limestone density
S = gamma/rho_w;    % limestone specific gravity
mu = 7.392*10^-4;   % viscosity of water lbm/(ft*s)
gc = 32.17;         % conversion factor between lbf and lbm
        
V = 7.2398;                 % Value chosen for V from calculated optimums

 
% design variables at mesh points
[D,d] = meshgrid(0:0.01:0.5,0:0.0005:0.01);
 
% equations
Concentration = 4*W./(pi*gamma*V*(D.^2));
Area = (pi/4)*D.^2;
rho = rho_w + Concentration.*(gamma-rho_w);
Pg = (218*W*((1./sqrt(d)) - (1/sqrt(a))))/550;
CdRpsq_calculated = 4*g*rho_w*(d.^3)*((gamma-rho_w)/(3*mu^2));
Cd = dragReynolds(CdRpsq_calculated);
Rw = (rho_w.*V.*D)/mu;
fw = fw_function(Rw);
F = fw.*((rho_w./rho) + 150*Concentration.*(rho_w./rho).*((g.*D.*(S-1))./((V.^2).*sqrt(Cd))).^1.5);
delta_p = (F.*rho*L*V.^2)./(2.*D*gc);    % 
Qslurry = Area.*V;
Pf = (delta_p.*Qslurry)/550;     % pump power, hp
Vc = ((40*g*Concentration*(S-1).*D)./sqrt(Cd)).^0.5;
horsepower = Pf + Pg;
cost = NPV(horsepower,0.07,7); % develop cost function here

figure(1)
[C,h] = contour(D,d,cost,[5e+05,6.2501e+05,1e+06,1e+07,1.1e+07],'k-');
clabel(C,h,'Labelspacing',250);
title('Limestone Slurry Contour Plot');
xlabel('Pipe Diameter');
ylabel('Particle Diameter After Grinding');
hold on;

% solid lines to show constraint boundaries
contour(D,d,((1.1*Vc) - V),[0,0],'b','LineWidth',2);
contour(D,d,(Concentration - 0.4),[0,0],'r','LineWidth',2);
contour(D,d,(D - 0.5),[0,0],'g','LineWidth',2);
contour(D,d,(d - 0.01),[0,0],'y','LineWidth',2);
contour(D,d,(0.0005 - d),[0,0],'c','LineWidth',2);

% show a legend
legend('Cost','1.1V_c <= V','C <= 0.4','D <= 0.5','d <= 0.01','d >= 0.0005')
hold off;
