function output = exactPressure( xMatrix, yMatrix,miu,t,L,v0  )
%UNTITLED7 Summary of this function goes here
%   Detailed explanation goes here

output = exp(-16*pi^2*miu*t/L^2)*(cos(2*(xMatrix-v0*t*2*pi)/L)+cos(2*(yMatrix-v0*t*2*pi)/L));

end

