
clear all,clc

L = 1;
miu = 0.05;
v0 = 1;


NX = 128;  NY = 128;


dt = 0.001;%0.1/NX;

KX=2*pi*ones(1,NY)'*(mod((1:NX)-ceil(NX/2+1),NX)-floor(NX/2)); % matrix of wavenumbers in x direction 
KY=2*pi*(mod((1:NY)'-ceil(NY/2+1),NY)-floor(NY/2))*ones(1,NX); % matrix of wavenumbers in y direction 

%dealias=KX<2/3*NX&KY<2/3*NY; % Cutting of frequencies using the 2/3 rule
dealias=abs(KX)<2*2*pi&abs(KY)<2*2*pi; 


Ksol = KX.^2+KY.^2;
Ksqure =Ksol; Ksqure(1,1)=4*pi^2;


[i,j]=meshgrid(0:NX-1,0:NY-1);
dx=2*pi/NX;
dy=2*pi/NY;


t=0;
u = exactU( i*dx, j*dy,miu,t,L,v0 );
v = exactV( i*dx, j*dy,miu,t,L,v0 );
Pressure = exactPressure( i*dx, j*dy,miu,t,L,v0  );


c = u; %initial 
chat = fft2(c);


chatX = 1i*KX.*chat;
chatY = 1i*KY.*chat;
cX = real(ifft2(chatX));
cY = real(ifft2(chatY));
uX = 4*pi/L*exp(-8*pi^2*miu*t/L^2)*sin((i*dx-v0*t*2*pi)/L).*sin((j*dy-v0*t*2*pi)/L);
uY = -4*pi/L*exp(-8*pi^2*miu*t/L^2)*cos((i*dx-v0*t*2*pi)/L).*cos((j*dy-v0*t*2*pi)/L);

N = u.*cX + v.*cY;

Phat = fft2(Pressure);
PhatX = 1i*KX.*Phat;


Nhat = fft2(N).*dealias;
chat = ((1/dt-0.5*miu*Ksol).*chat-(Nhat+PhatX))./(1/dt+0.5*miu*Ksol);

t=t+dt;
u = exactU(  i*dx, j*dy,miu,t,L,v0 );

norm(u-real(ifft2(chat)))

figure(1);clf
pcolor(u);shading flat;colorbar 
figure(2);clf
pcolor(real(ifft2(chat)));shading flat;colorbar 
Nhatold = Nhat;


while t<0.25
   
   
   u = exactU( i*dx, j*dy,miu,t,L,v0 );
   v = exactV( i*dx, j*dy,miu,t,L,v0 );
   Pressure = exactPressure( i*dx, j*dy,miu,t,L,v0  );
      
   
   Phat = fft2(Pressure);
   chatX = 1i*KX.*chat;
   chatY = 1i*KY.*chat;
  
   cX = real(ifft2(chatX));
   cY = real(ifft2(chatY));
   N = u.*cX + v.*cY;
   Nhat = fft2(N).*dealias;
   chat = CNfun( what,miu*Ksol,-1i*KX.*Phat-AdamsBashforth(Nhat,Nhatold),dt);
   
   Nhatold = Nhat;
   t=t+dt;
   figure(2);clf
   pcolor(real(ifft2(chat)));shading flat;colorbar 
   title(num2str(t));
   drawnow
    
end


