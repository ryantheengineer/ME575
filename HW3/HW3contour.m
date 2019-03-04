% This program constructs a contour plot for HW2
clear; close all;

% get search pattern data
load Rosenbrock_search.mat
x1search = xsearch(1,:);
x2search = xsearch(2,:);

% design variables at mesh points
[x1,x2] = meshgrid(-3:0.01:3,-1.5:0.01:6);
 
% equations
f = 100*(x2 - x1.^2).^2 + (1 - x1).^2;

figure(1)
[C,h] = contour(x1,x2,f,[1,10,100,500,1000],'k-');
clabel(C,h,'Labelspacing',250);
title({'Rosenbrock''s function:';'Quasi-Newton BFGS'});
xlabel('x_1');
ylabel('x_2');
hold on;

plot(x1search,x2search,'r','LineWidth',2)

% show a legend
legend('Rosenbrock''s function','Quasi-Newton BFGS')
hold off;
