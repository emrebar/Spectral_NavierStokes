clc;clear;close;
%% u_t = u_yy+ u_xx   - y Chebyshev, x Fourier RK4
Nx=8;Ny=8;
x = (0:Nx-1)/Nx*2*pi;         % x coordinate
kx = fftshift(-Nx/2:Nx/2-1);   % wave vector
[D,y] = cheb(Ny);           %cos(pi*(0:Ny)/Ny);
D2 = D*D;
[xx,yy]= meshgrid (x,y);
p=1e4;
tfinal = 0.01;
dt   = 1e-4/2;
Nstep = ceil(tfinal/dt);
t    = (0:Nstep)*dt;

ufun    = @(x,y,t) cos(3*x+3*y+5*t);
yufun  = @(x,y,t) -3*sin(3*x+3*y+5*t);
dy2ufun = @(x,y,t) -9*cos(3*x+3*y+5*t);
dx2ufun = @(x,y,t) -9*cos(3*x+3*y+5*t);
ufun_t  = @(x,y,t) -5*sin(3*x+3*y+5*t);
bfun    = @(x,y,t) ufun_t(x,y,t) - dy2ufun(x,y,t) - dx2ufun(x,y,t);

u = ufun(xx,yy,0);


% K=cell(step/p,1);

for j=1:Nstep
    %     if mod(j,1)==0, waitbar(1*j/step,h),end
    
    uhat=(fft(u,[],2));
    b = bfun(xx,yy,t(j)); bhat1=(fft(b,[],2));   %b(t) for k1
    b2 =bfun(xx,yy,t(j)+dt/2); bhat2=(fft(b2,[],2));%b(t+dt/2) for k2&k3
    b4 =bfun(xx,yy,t(j)+dt); bhat4=(fft(b4,[],2));%b(t+dt) for k4
    
    for m=1:length(kx)
        k1=(D2-kx(m)^2*eye(Ny+1))*uhat(:,m)+bhat1(:,m);
        k2=(D2-kx(m)^2*eye(Ny+1))*(uhat(:,m)+(dt/2).*k1)+bhat2(:,m);
        k3=(D2-kx(m)^2*eye(Ny+1))*(uhat(:,m)+(dt/2).*k2)+bhat2(:,m);
        k4=(D2-kx(m)^2*eye(Ny+1))*(uhat(:,m)+dt.*k3)+bhat4(:,m);
        uhat(:,m)=uhat(:,m)+dt/6*(k1+2*k2+2*k3+k4);
    end
    u=real(ifft(uhat,[],2));
    u([1 Ny+1],:) = ufun(xx([1 Ny+1],:),yy([1 Ny+1],:),t(j+1));

    if mod(j,p)==0;
        %         K{j/p}=u;
        
        figure(1),clf(1);
        
        uexx=ufun(xx,yy,t(j+1));
        mesh(xx,yy,u,0.8*ones(size(xx)));hold on;
        mesh(xx,yy,uexx,0.2*ones(size(xx)));
        title(sprintf('time=%-4.3e',t(j))),drawnow,shg,pause(0.5)
    end
    
end


uex=ufun(xx,yy,t(j+1));
error=u-uex
max(abs(uex(:)-u(:)))


% %%
% close all;figure(1),clf(1);
% for ii=1:100:length(t);
%
%     mesh(xx,yy,ufun(xx,yy,t(ii)));
%     title(sprintf('time=%-4.3e',t(ii))),grid on,drawnow,shg
% end