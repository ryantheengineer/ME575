% This program constructs a contour plot for HW2
clear; close all;

% get search pattern data
load Rosenbrock_search.mat
x1search = xsearch(1,:);
x2search = xsearch(2,:);

% design variables at mesh points
[x1,x2] = meshgrid(-1.75:0.01:1.75,-1.5:0.01:3);
 
% equations
f = 100*(x2 - x1.^2).^2 + (1 - x1).^2;

figure(1)
[C,h] = contour(x1,x2,f,[1,10,50,100,200],'k-');
clabel(C,h,'Labelspacing',250);
title('Rosenbrock''s function');
xlabel('x_1');
ylabel('x_2');
hold on;

plot(x1search,x2search,'r','LineWidth',2)

% solid lines to show constraint boundaries
% contour(D,d,((1.1*Vc) - V),[0,0],'b','LineWidth',2);
% contour(D,d,(Concentration - 0.4),[0,0],'r','LineWidth',2);
% contour(D,d,(D - 0.5),[0,0],'g','LineWidth',2);
% contour(D,d,(d - 0.01),[0,0],'y','LineWidth',2);
% contour(D,d,(0.0005 - d),[0,0],'c','LineWidth',2);

% show a legend
legend('Rosenbrock''s function','Quasi-Newton search pattern')
hold off;
