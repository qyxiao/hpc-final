
clear all; clc

L = 1;
miu = 0.05;
v0 = 1;


error = zeros(5,3);
gridNum=[16 32 64 128 256];
for index = 1:5


NX = gridNum(index);  NY = gridNum(index);

KX=2*pi*ones(1,NY)'*(mod((1:NX)-ceil(NX/2+1),NX)-floor(NX/2)); % matrix of wavenumbers in x direction 
KY=2*pi*(mod((1:NY)'-ceil(NY/2+1),NY)-floor(NY/2))*ones(1,NX); % matrix of wavenumbers in y direction 

dealias=KX<2/3*NX&KY<2/3*NY; % Cutting of frequencies using the 2/3 rule

Ksol = KX.^2+KY.^2;
Ksqure =Ksol; Ksqure(1,1)=4*pi^2;


[i,j]=meshgrid(0:NX-1,0:NY-1);
dx=2*pi/NX;
dy=2*pi/NY;


t=0;

u = exactU(  i*dx, j*dy,miu,t,L,v0 );
v = exactV(  i*dx, j*dy,miu,t,L,v0  );


%initial  
c = (sin(i*dx).*sin(j*dy)).^100;
chat = fft2(c);

chatX = 1i*KX.*chat;
chatY = 1i*KY.*chat;
cX = real(ifft2(chatX));
cY = real(ifft2(chatY));


N = u.*cX + v.*cY;

Ntest = u.*(100*2*pi/L*(sin(i*dx).*sin(j*dy)).^99.*sin(j*dy).*cos(i*dx)) ...
    + v.*(100*2*pi/L*(sin(i*dx).*sin(j*dy)).^99.*cos(j*dy).*sin(i*dx));
error(index,1)=norm(Ntest-N,1)/NX^2;
error(index,2)=norm(Ntest-N,2)/NX;
error(index,3)=norm(Ntest-N,inf);
end

format long

logPSerror = log(error);

figure(1)
clf
plot(log([16 32 64 128 256]), logPSerror(:,1),'ro-',log([16 32 64 128 256]),logPSerror(:,2),'go-',log([16 32 64 128 256]), logPSerror(:,3),'bo-');hold on

legend('L1 norm', 'L2 norm', 'L inf norm')
ylabel('Log(error)')
xlabel('Log(grid point numbers)')
