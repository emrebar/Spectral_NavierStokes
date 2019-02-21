clc;clear;close;
%% u_t -  u_yy + k^2*u + b
% y Chebyshev, x Fourier

%% ---------------------EULER------------------------------------------ %%
clc;clear all;close all;
Ny=16;
kx = 3;
[D,y] = cheb(Ny);           %cos(pi*(0:Ny)/Ny);
yi = linspace(-1,1,100);
D2 = D*D;
ufun = @(y,t) exp(y+t);
% ufun_y = @(y,t) exp(y+t);
% ufun_t = @(y,t) exp(y+t);
% ufun_yy = @(y,t) exp(y+t);
u = ufun(y,0);
p=1e4;         %store
step=1e4;
dt= 1e-7;
t=(0:step)*dt;
I=eye(Ny+1);

A = dt*D2-kx^2*dt*I+I;
for j=1:step
    u=A*u;
    u(1)=ufun(y(1),t(j+1));    u(Ny+1)=ufun(y(Ny+1),t(j+1));
    if mod(j,p)==0;
        ue = ufun(y,t(j+1));
        uei = ufun(yi,t(j+1));
        subplot(121), plot(u,y,'b.-',uei,yi,'r-')
        subplot(122), plot(u(:)-ue(:),y,'g.-')
        title(sprintf('BD2   time=%-4.3e Nx %d Ny %d dt=%d.',t(j+1),Ny,Ny,dt)),grid on,drawnow,shg,pause(0.1)
        
    end
end
error=max(abs(ue(:)-u(:)))
%% -------------------------RK2-------------------------------------- %%
clc;clear all;close all;
Ny=16;
kx = 2;
[D,y] = cheb(Ny);           %cos(pi*(0:Ny)/Ny);
yi = linspace(-1,1,100);
D2 = D*D;
ky=6;
ufun = @(y,t) cos(ky*y+t);
ufun_y = @(y,t) -ky*sin(ky*y+t);
ufun_t = @(y,t) -sin(ky*y+t);
ufun_yy = @(y,t) -ky*ky*cos(ky*y+t);
bfun = @(y,t) ufun_t(y,t) - ufun_yy(y,t) + kx^2*ufun(y,t) ;
u = ufun(y,0);
p=5e2;         %store
step=5e2;
dt= 1e-4/1;
t=(0:step)*dt;
U=zeros(step,Ny+1);
I=eye(Ny+1);
A = D2-kx^2*I;

for j=1:step
    k1=(A*u+bfun(y,t(j)))*dt;
    k2=(A*(u+(k1/2))+bfun(y,t(j)+dt/2))*dt;
    
    u=u+k2;
    u(1)=ufun(y(1),t(j+1));    u(Ny+1)=ufun(y(Ny+1),t(j+1));
    if mod(j,p)==0;
        ue = ufun(y,t(j+1));
        uei = ufun(yi,t(j+1));
        subplot(121), plot(u,y,'b.-',uei,yi,'r-')
        subplot(122), plot(u(:)-ue(:),y,'g.-')
        title(sprintf('BD2   time=%-4.3e Nx %d Ny %d dt=%d.',t(j),Ny,Ny,dt)),grid on,drawnow,shg,pause(0.1)
        
    end
end
error=max(abs(ue(:)-u(:)))
%% -------------------------RK4-------------------------------------- %%
clc;clear all;close all;
Ny=16;
kx = 4;
[D,y] = cheb(Ny);           %cos(pi*(0:Ny)/Ny);
yi = linspace(-1,1,100);
D2 = D*D;
ky=4;
ufun = @(y,t) cos(ky*y+t);
ufun_y = @(y,t) -ky*sin(ky*y+t);
ufun_t = @(y,t) -sin(ky*y+t);
ufun_yy = @(y,t) -ky*ky*cos(ky*y+t);
bfun = @(y,t) ufun_t(y,t) - ufun_yy(y,t) + kx^2*ufun(y,t) ;
u = ufun(y,0);
p=1e6;         %store
time=0.5;
dt=1e-3/10;
step=floor(time/dt);
t=(0:step)*dt;


A = D2-kx^2*eye(Ny+1);

for j=1:step
    k1=A*u+bfun(y,t(j));
    k2=A*(u+(dt/2)*k1)+bfun(y,t(j)+dt/2);
    k3=A*(u+(dt/2)*k2)+bfun(y,t(j)+dt/2);
    k4=A*(u+(dt*k3))+bfun(y,t(j)+dt);
    u=u+dt/6*(k1+2*k2+2*k3+k4);
    u(1)=ufun(y(1),t(j+1));    u(Ny+1)=ufun(y(Ny+1),t(j+1));
    if mod(j,p)==0;
        ue = ufun(y,t(j+1));
        uei = ufun(yi,t(j+1));
        subplot(121), plot(u,y,'b.-',uei,yi,'r-')
        subplot(122), plot(u(:)-ue(:),y,'g.-')
        title(sprintf('BD2   time=%-4.3e Nx %d Ny %d dt=%d.',t(j+1),Ny,Ny,dt)),grid on,drawnow,shg,pause(0.1)
        
    end
end
uex=ufun(y,t(j+1));
error=max(abs(uex(:)-u(:)))


%     case 2
%% -------------------------BD2-------------------------------------- %%
clc;clear all;close all;
Ny=16;
kx = 31;
[D,y] = cheb(Ny);           %cos(pi*(0:Ny)/Ny);
yi = linspace(-1,1,100);
D2 = D*D;
ufun = @(y,t) cos(y+t);
ufun_y = @(y,t) -sin(y+t);
ufun_t = @(y,t) -sin(y+t);
ufun_yy = @(y,t) -cos(y+t);
bfun = @(y,t) ufun_t(y,t) - ufun_yy(y,t) + kx^2*ufun(y,t) ;
u = ufun(y,0);
I=eye(Ny+1);
p=1e5;   %store
time=4;
dt=1e-4/4;
step=floor(time/dt);
t=(0:step)*dt;
uold=ufun(y,-dt);

A=(3/(2*dt))*I-D2+kx^2*I;

for j=1:step;
    c=(2/dt)*u-(1/(2*dt))*uold+bfun(y,t(j)+dt);
    unew=A\c;
    unew(1)=ufun(y(1),t(j+1));    unew(end)=ufun(y(end),t(j+1));
    
    uold=u;u=unew;
    if mod(j,p)==0;
        %             U(j,:)=u;
        ue = ufun(y,t(j+1));
        uei = ufun(yi,t(j+1));
        subplot(121), plot(u,y,'b.-',uei,yi,'r-')
        subplot(122), plot(u(:)-ue(:),y,'g.-')
        title(sprintf('BD2   time=%-4.3e Nx %d Ny %d dt=%d.',t(j),Ny,Ny,dt)),grid on,drawnow,shg,pause(0.1)
        
    end
end
ue = ufun(y,t(j+1));
error=max(abs(ue(:)-u(:)))

%% -------------------------BD3-------------------------------------- %%
clc;clear all;close all;
Ny=16;
kx = 4;
[D,y] = cheb(Ny);           %cos(pi*(0:Ny)/Ny);
yi = linspace(-1,1,1000);
D2 = D*D;
ufun = @(y,t) cos(y+t);
ufun_y = @(y,t) -sin(y+t);
ufun_t = @(y,t) -sin(y+t);
ufun_yy = @(y,t) -cos(y+t);
bfun = @(y,t) ufun_t(y,t) - ufun_yy(y,t) + kx^2*ufun(y,t) ;
u = ufun(y,0);
I=eye(Ny+1);
p=1e5;     % store
step=5e4*2;
dt= 1e-4/2;
t=(0:step)*dt;
U=zeros(step,Ny+1);
uold=ufun(y,-dt);
uoldd=ufun(y,-2*dt);

A=(11/(6*dt))*I-D2+kx^2*I;

for j=1:step;
    c=(3/dt)*u-(3/(2*dt))*uold+1/(3*dt)*uoldd+bfun(y,t(j)+dt);
    unew=A\c;
    unew(1)=ufun(y(1),t(j+1));    unew(end)=ufun(y(end),t(j+1));
    uoldd=uold;uold=u;u=unew;
    
    if mod(j,p)==0;
        ue = ufun(y,t(j+1));
        uei = ufun(yi,t(j+1));
        subplot(121), plot(u,y,'b.-',uei,yi,'r-')
        subplot(122), plot(u(:)-ue(:),y,'g.-')
        title(sprintf('BD2   time=%-4.3e Nx %d Ny %d dt=%d.',t(j),Ny,Ny,dt)),grid on,drawnow,shg,pause(0.1)
        
    end
end
ue = ufun(y,t(j+1));
error=max(abs(ue(:)-u(:)))

