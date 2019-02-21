clc;clear;close;
%% u_t = u_yy+ u_xx   - y Chebyshev, x Fourier RK4
Nx=20;Ny=16;
x = (0:Nx-1)/Nx*2*pi;         % x coordinate
kx = fftshift(-Nx/2:Nx/2-1);   % wave vector
eta    = -cos(pi*(0:Ny)/Ny)';
dVeta         = diag(1./sqrt(1-eta.^2))*sin(acos(eta)*(0:Ny))*diag(0:Ny);
dVeta(1,:)    = (-1).^(1:Ny+1).*(0:Ny).^2;
dVeta(Ny+1,:)  = (0:Ny).^2;
Veta          = cos(acos(eta)*(0:Ny));
D         = dVeta/Veta;
y=eta;
[D,y]=cheb(Ny);
D2 = D*D;
[xx,yy]= meshgrid (x,y);
p=1e4;
tfinal = 1;
dt   = 1e-4;
Nstep = ceil(tfinal/dt);
t    = (0:Nstep)*dt;
time = t;
Lx = 17;
I=eye(Ny+1);
BC = -D([1 Ny+1],[1 Ny+1])\D([1 Ny+1],2:Ny);

ufun     = @(x,y,t) sin(7*x+6*y+5*t);      % ufun =@(x,y,t) 2./(cos(5*(x+7*t)) + 2).* sinh(y);
dudyfun  = @(x,y,t) 6*cos(7*x+6*y+5*t);    % dudyfun  = @(x,y,t) (2.*cosh(y))./(cos(35*t + 5*x) + 2);
d2udyfun = @(x,y,t) -6*6*sin(7*x+6*y+5*t);  % d2udyfun = @(x,y,t) (2*sinh(y))./(cos(35*t + 5*x) + 2);
dudxfun  = @(x,y,t) 7*cos(7*x+6*y+5*t);     %dudxfun  = @(x,y,t) (10.*sin(35.*t + 5.*x).*sinh(y))./(cos(35.*t + 5.*x) + 2).^2;
d2udxfun = @(x,y,t) -49*sin(7*x+6*y+5*t);   %d2udxfun = @(x,y,t) (100.*sin(35.*t + 5.*x).^2.*sinh(y))./(cos(35.*t + 5.*x) + 2).^3 + (50*cos(35*t + 5*x).*sinh(y))./(cos(35*t + 5*x) + 2).^2; 
ufun_t   = @(x,y,t) 5*cos(7*x+6*y+5*t);    %ufun_t   = @(x,y,t) (70*sin(35*t + 5*x).*sinh(y))./(cos(35*t + 5*x) + 2).^2;

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
        c(1)=0;
        A(1,:)=D(1,:);
        unewhat(:,m)=A\c;
        
    end
    unew=real(ifft(unewhat,[],2));
%     unew([1 Ny+1],:) = ufun(xx([1 Ny+1],:),yy([1 Ny+1],:),t(j+1));

%     unew([1 Ny+1],:) = BC*unew(2:Ny,:);   % Neumann BCs for |y| = 1

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
