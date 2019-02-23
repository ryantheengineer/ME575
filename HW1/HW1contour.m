% This program constructs a contour plot for the two-bar truss
% constants
clear; close all;
h0 = 1.0;               % preload height, in
delta0 = 0.4;           % deflection, in
hdef = h0 - delta0;     % deflection height, in
G = 12*10^6;            % psi
Se = 45000;             % psi
w = 0.18;
Sf = 1.5;
Q = 150000;             % psi
delta_f = 0;            % in

n = 7.5928;             % Values chosen for n and hf from calculated optimums
hf = 1.3691;

 
% design variables at mesh points
[d,D] = meshgrid(0.01:0.01:0.12,0.2:0.05:1);
 
% equations
k = (G*d.^4)./(8*(D.^3)*n);                 % spring stiffness
delta_p = hf-h0;                            % deflection at preload
delta_def = delta0 + (hf-h0);               % greatest working deflection
hs = n.*d;                                  % solid height
delta_s = hf - hs;                          % deflection at solid height
F_f = k.*delta_f;                           % full height force (zero)
F_p = k.*delta_p;                           % preload force
F_def = k.*delta_def;                       % force at full deflection
F_s = k.*delta_s;                           % force at solid height
K = ((4.*D)-d)./(4.*(D-d))+0.62.*(d./D);    % Wahl factor
tau_f = ((8.*F_f.*D)./(pi.*d.^3)).*K;       % stress at full height (zero)
tau_p = ((8.*F_p.*D)./(pi.*d.^3)).*K;       % stress at preload height
tau_def = ((8.*F_def.*D)./(pi.*d.^3)).*K;   % stress at full deflection
tau_s = ((8.*F_s.*D)./(pi.*d.^3)).*K;       % stress at solid height
tau_max = tau_def;                          % max stress (assumed at max deflection)
tau_min = tau_p;                            % min stress (assumed at preload deflection)
tau_m = (tau_max + tau_min)./2;             % mean stress
tau_a = (tau_max - tau_min)./2;             % alternating stress
Sy = 0.44.*(Q./(d.^w));                     % yield stress
dratio = D./d;                              % ratio of diameters
dsum = D + d;                               % sum of diameters
clash = hdef - hs;                          % clash allowance
Sefratio = Se./Sf;                          % endurance limit to safety factor
Syfratio = Sy./Sf;                          % yield strength to safety factor
tau_amsum = tau_a + tau_m;                  % sum of alternating and mean stresses

figure(1)
[C,h] = contour(d,D,F_p,[0,2,4,6,8,10,20,50],'k-');
clabel(C,h,'Labelspacing',250);
title('Spring Contour Plot');
xlabel('Wire Diameter');
ylabel('Coil Diameter');
hold on;

% solid lines to show constraint boundaries
contour(d,D,(tau_a-Sefratio),[0.0,0.0],'b','LineWidth',2);
contour(d,D,(tau_amsum-Syfratio),[100,100],'r','LineWidth',2);
contour(d,D,(dratio-16),[0.25,0.25],'g','LineWidth',2);
contour(d,D,(4-dratio),[0.25,0.25],'y','LineWidth',2);
contour(d,D,(dsum-0.75),[0.25,0.25],'c','LineWidth',2);
contour(d,D,(0.05-clash),[0.25,0.25],'m','LineWidth',2);
contour(d,D,(tau_s-Sy),[0.25,0.25],'LineWidth',2);

% show a legend
legend('F_p','\tau_a <= S_e/S_f','\tau_a+\tau_m <= S_y/S_f','D/d <= 16',...
       'D/d >= 4','D + d <= 0.75','clash <= 0.05','\tau_s <= S_y')
hold off;
