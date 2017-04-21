function output = exactVorticity( xMatrix, yMatrix,miu,t,L,v0 )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
 
 output= 8*pi/L*exp(-8*pi^2*miu*t/L^2).*cos((xMatrix-v0*t*2*pi)/L).*cos((yMatrix-v0*t*2*pi)/L);

end

