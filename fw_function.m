function [fw] = fw_function(Rw)
%fw_function: Calculates the value of fw based on a given Rw

if Rw <= 10^5
    fw = 0.3164/(Rw.^0.25);
else
    fw = 0.0032 + 0.221*Rw.^(-0.237);
end

end

