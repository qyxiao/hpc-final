function output= exactV( xMatrix, yMatrix,miu,t,L,v0  )
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here

output = v0 + 2*exp(-8*pi^2*miu*t/L^2)*sin((xMatrix-v0*t*2*pi)/L).*cos((yMatrix-v0*t*2*pi)/L);
end

