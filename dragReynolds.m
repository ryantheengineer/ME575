function [Cd_interpolated] = dragReynolds(CdRpsq_calculated)
%dragReynolds: Interpolates the empirical relationship between Cd and
%CdRp^2 using data as provided in HW2

load('Cd_CdRp_sq_correlation.mat','Cd_vector','CdRpsq');

Cd_interpolated = interp1(Cd_vector,CdRpsq,CdRpsq_calculated,'pchip');

end

