clc;clear;close;
%% u_t = u_yy+ u_xx   - y Chebyshev, x Fourier RK4
Nx=12;Ny=20;
x = (0:Nx-1)/Nx*2*pi;         % x coordinate
kx = fftshift(-Nx/2:Nx/2-1);   % wave vector
y    = -cos(pi*(0:Ny)/Ny)';

% Dxi = zeros(Ny+1);
% for ii=1:Ny+1
%    Dxi(ii,:) = fdcoeffF(2,y(ii),y);
% end
% D2=Dxi;
[D]=cheb(Ny);D2=D^2:
[xx,yy]= meshgrid (x,y);
p=1e4;
tfinal = 0.01;
dt   = 1e-5;
Nstep = ceil(tfinal/dt);
t    = (0:Nstep)*dt;
time = t;
Lx = 17;
I=eye(Ny+1);

ufun     = @(x,y,t) sin(7*x+3*y+5*t);
dudyfun  = @(x,y,t) 3*cos(7*x+3*y+5*t);
d2udyfun = @(x,y,t) -9*sin(7*x+3*y+5*t);
dudxfun  = @(x,y,t) 7*cos(7*x+3*y+5*t);
d2udxfun = @(x,y,t) -49*sin(7*x+3*y+5*t);
ufun_t   = @(x,y,t) 5*cos(7*x+3*y+5*t);
ududxfun = @(x,y,t) ufun(x,y,t).*dudxfun(x,y,t);
bfun     = @(x,y,t) ufun_t(x,y,t)-d2udxfun(x,y,t)-d2udyfun(x,y,t);
u = ufun(xx,yy,0);
uold=ufun(xx,yy,-dt);

for j=1:Nstep
    
    uhat=fft(u,[],2);
    uoldhat=fft(uold,[],2);
    b =bfun(xx,yy,t(j)+dt); bh=fft(b,[],2);
    
    for m=1:length(kx);
        A=(3/(2*dt))*I-D2+kx(m)^2*I;
        c=(2/dt)*uhat(:,m)-(1/(2*dt))*uoldhat(:,m)+bh(:,m);
        unewhat(:,m)=A\c;
        
    end
    unew=real(ifft(unewhat,[],2));
    unew([1 Ny+1],:) = ufun(xx([1 Ny+1],:),yy([1 Ny+1],:),t(j+1));

    uold=real(ifft(uhat,[],2));
    u=unew;

    
end
%%
uex=ufun(xx,yy,t(end));
error=u-uex;
max(abs(uex(:)-u(:)))

%% Plot
%     if mod(j,p)==0;U{j/p}=u;
%            ue =  ufun(xx,yy,t(j+1));
% 
%         subplot(131),mesh(xx,yy,u);%set(s1,'FaceColor','none','EdgeColor','k');hold on;
%         subplot(132),mesh(xx,yy,ue);%set(s2,'FaceColor','interp','EdgeColor','interp')
%         subplot(133),mesh(xx,yy,u-ue);
%         title(sprintf('time=%-4.3e Nx %d Ny %d dt=%d.',t(j+1),Nx,Ny,dt));
%         grid on,drawnow,shg
%         title(sprintf('time=%-4.3e Nx %d Ny %d dt=%d.',t(j+1),Nx,Ny,dt));
%         grid on,drawnow;shg;
%     end
