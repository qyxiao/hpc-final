function output = exactU(  xMatrix, yMatrix,miu,t,L,v0 )
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

 output =v0 - 2*exp(-8*pi^2*miu*t/L^2)*cos((xMatrix-v0*t*2*pi)/L).*sin((yMatrix-v0*t*2*pi)/L);

end

