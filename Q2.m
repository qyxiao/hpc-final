
clear all

L = 1;
miu = 0.05;
v0 = 1;

NX = 128;  NY = 128;

dt = 0.1/NX;

KX=2*pi*ones(1,NY)'*(mod((1:NX)-ceil(NX/2+1),NX)-floor(NX/2)); % matrix of wavenumbers in x direction 
KY=2*pi*(mod((1:NY)'-ceil(NY/2+1),NY)-floor(NY/2))*ones(1,NX); % matrix of wavenumbers in y direction 

dealias=KX<2/3*NX&KY<2/3*NY; % Cutting of frequencies using the 2/3 rule
%dealias=abs(KX)<2*2*pi&abs(KY)<2*2*pi; 

Ksol = KX.^2+KY.^2;
Ksqure =Ksol; Ksqure(1,1)=4*pi^2;


[i,j]=meshgrid(0:NX-1,0:NY-1);
dx=2*pi/NX;
dy=2*pi/NY;


t=0;
u = exactU( i*dx, j*dy,miu,t,L,v0 );
v = exactV( i*dx, j*dy,miu,t,L,v0 );


%initial  
c = (sin(i*dx/2).*sin(j*dy/2)).^100;
chat = fft2(c);
figure(1);clf
pcolor(real(ifft2(chat)));shading flat;colorbar 

chatX = 1i*KX.*chat;
chatY = 1i*KY.*chat;
cX = real(ifft2(chatX));
cY = real(ifft2(chatY));


N = u.*cX + v.*cY;
Nhat = fft2(N).*dealias;

chat = chat - dt*Nhat;
Nhatold = Nhat;


while t<0.25
   
   t=t+dt;
   u = exactU( i*dx, j*dy,miu,t,L,v0 );
   v = exactV( i*dx, j*dy,miu,t,L,v0 );
   
   chatX = 1i*KX.*chat;
   chatY = 1i*KY.*chat;
  
   cX = real(ifft2(chatX));
   cY = real(ifft2(chatY));
   N = u.*cX + v.*cY;
   Nhat = fft2(N).*dealias;
   chat = chat - dt*(1.5*Nhat-0.5*Nhatold);
   Nhatold = Nhat;
   figure(1);clf
   pcolor(real(ifft2(chat)));shading flat;colorbar 
   title(num2str(t));
   drawnow
    
end


