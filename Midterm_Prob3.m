%% Part a
h = 0.001;

x = 3;
fx = x^2.5 + log(x);
fxh = (x+h)^2.5 + log(x+h);

dfdx_forward = (fxh - fx)/h;

dfdx_analytic = 1/x + (5*x^(3/2))/2;

err_forward = abs(dfdx_analytic - dfdx_forward);

%% Part b

fxplush = fxh;
fxminush = (x-h)^2.5 + log(x-h);

dfdx_central = (fxplush - fxminush)/(2*h);

err_central = abs(dfdx_analytic - dfdx_central);

%% Part c

h = 1e-9;
dfdx_complex = imag((complex(x,h))^2.5 + log(complex(x,h)))/h;

err_complex = abs(dfdx_analytic - dfdx_complex);