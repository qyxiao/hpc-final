function output = CNfun( w,matrix,Source,dt)
%UNTITLED2 Summary of this function goes here
%   Crank-Nicolson
  
  output = ((1/dt - 0.5*matrix).*w+Source)./(1/dt + 0.5*matrix);

end

