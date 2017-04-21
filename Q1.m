

clear all,  clc

L = 1;
miu = 0.05;
v0 = 1;
% assume N to be even
timeNum = [ 32  64 128];

PSerror=zeros(3,3); % (i,j) i for norm and j for gridNum index

ratio = 5;
for index = 1:3   
NX = 128;  NY = 128;
dt = ratio/timeNum(index);
t=0;


KX=2*pi*ones(1,NY)'*(mod((1:NX)-ceil(NX/2+1),NX)-floor(NX/2)); % matrix of wavenumbers in x direction 
KY=2*pi*(mod((1:NY)'-ceil(NY/2+1),NY)-floor(NY/2))*ones(1,NX); % matrix of wavenumbers in y direction 

dealias=KX<2/3*NX&KY<2/3*NY; % Cutting of frequencies using the 2/3 rule

Ksol = KX.^2+KY.^2;
Ksqure =Ksol; Ksqure(1,1)=1;


[i,j]=meshgrid(0:NX-1,0:NY-1);
dx=2*pi/NX;
dy=2*pi/NY;

w = exactVorticity( i*dx, j*dy,miu,t,L,v0 );  % t = 0
what=fft2(w);

Stream = what./Ksqure;
uhat = 1i*KY.*Stream; vhat = -1i*KX.*Stream;
whatX = 1i*KX.*what; whatY = 1i*KY.*what;
u = real(ifft2(uhat))+v0; v = real(ifft2(vhat))+v0;  
wX = real(ifft2(whatX)); wY = real(ifft2(whatY));

N = u.*wX + v.*wY;
Nhat = (fft2(N)).*dealias;
what=CNfun( what,miu*Ksol,-Nhat,dt);
Nhatold = Nhat;   % t = 1*dt
t=t+dt;
vorticity =exactVorticity( i*dx, j*dy,miu,t,L,v0 );


while t<0.25

Stream = what./Ksqure;
uhat = 1i*KY.*Stream; vhat = -1i*KX.*Stream;
whatX = 1i*KX.*what; whatY = 1i*KY.*what;
u = real(ifft2(uhat))+v0; v = real(ifft2(vhat))+v0;  
wX = real(ifft2(whatX)); wY = real(ifft2(whatY));

N = u.*wX + v.*wY;
Nhat = (fft2(N)).*dealias;

what=CNfun( what,miu*Ksol,-AdamsBashforth(Nhat,Nhatold),dt);
Nhatold = Nhat;
t=t+dt;
vorticity =exactVorticity( i*dx, j*dy,miu,t,L,v0 );

end

  PSerror(1,index)=norm(vorticity-real(ifft2(what)),1)*1/(NX)^2;
  PSerror(2,index)=norm(vorticity-real(ifft2(what)),2)*1/NX;
  PSerror(3,index)=norm(vorticity-real(ifft2(what)),inf);

end


logPSerror = log(PSerror);

figure(1)
clf
plot(log([64 128 256 ]*1/ratio), logPSerror(1,:),'ro-',log([64 128 256 ]*1/ratio),logPSerror(2,:),'go-',log([64 128 256 ]*1/ratio), logPSerror(3,:),'bo-');hold on

comp = logPSerror(1,1)-2*(log([64 128 256 ]*1/ratio)-log(64*1/ratio));
plot(log([64 128 256 ]*1/ratio),comp,'k--')
legend('L1 norm', 'L2 norm', 'L inf norm','y=a-2*b')
ylabel('Log(error)')
xlabel('Log(time steps number)')



