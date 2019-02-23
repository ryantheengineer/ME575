function [totalcost] = NPV(horsepower,i,n)
%NPV.m: Calculates the total cost of operation using net present value
%methods, with interest calculated once per year (n = 7).

capital_cost = 300*horsepower + 200*horsepower;  % Capital cost of grinder and pump
power_cost = horsepower*8*300*0.07; % Power cost, $/year

totalcost = 0;
for t = 0:n
    if t == 0
        Rt = capital_cost + power_cost;
    else
        Rt = power_cost;
    end
    
    totalcost = totalcost + Rt/((1+i)^t);
end

end

