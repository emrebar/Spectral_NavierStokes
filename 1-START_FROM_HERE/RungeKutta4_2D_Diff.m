% This code solves the 2D time dependent equation
% u_t = u_yy+ u_xx + b
% using Runge Kutta Fourth order time marching scheme
% and Spatial discretization is carried out with :
% y Chebyshev Polynomial approximation, x Fourier approximation.
% b function can be used to study convergence;
% So-called modified functions method was implemented where the exact solution is known and
% with the help of an added right hand side the approximation was done. 
clc;clear;close;
%%
Nx=64;Ny=40; % number of modes in x and y directions
x = (0:Nx-1)/Nx*2*pi;        % x coordinate [0,2*pi) periodicity via Fourier
kx = fftshift(-Nx/2:Nx/2-1); % corresponding wave vector
[D,y] = cheb(Ny);            % Chebyshev polynomial first order differentiation matrix cos(pi*(0:Ny)/Ny);
D2 = D*D;                    % Chebyshev polynomial second order differentiation matrix 
[xx,yy]= meshgrid (x,y);     % Create the 2D grid.
tfinal = 2;                  % End time (non dimensional)
dt   = 1e-5;                 % time step
Nstep = ceil(tfinal/dt);     % number of time steps
t    = (0:Nstep)*dt;         % time array will be used to evaluate the exact funciton value
p=1e4/3;                     % every pth time step visualization will be shown
I=eye(Ny+1);                 % Identity matrix
% Function to be approximated and its derivatives
ufun    = @(x,y,t) cos(3*x+3*y+5*t);    % u
yufun   = @(x,y,t) -3*sin(3*x+3*y+5*t); % dudy
dy2ufun = @(x,y,t) -9*cos(3*x+3*y+5*t); % du2dy2
dx2ufun = @(x,y,t) -9*cos(3*x+3*y+5*t); % du2dx2
ufun_t  = @(x,y,t) -5*sin(3*x+3*y+5*t); % dudt
% right hand side function 
bfun    = @(x,y,t) ufun_t(x,y,t) - dy2ufun(x,y,t) - dx2ufun(x,y,t);
% initialize
u = ufun(xx,yy,0);
for j=1:Nstep % Start time marching
    uhat=(fft(u,[],2)); % Convert to Fourier Space in X dimension. (dim:2)
    % terms required for runge kutta
    b = bfun(xx,yy,t(j)); bhat1=(fft(b,[],2));   %b(t) for k1
    b2 =bfun(xx,yy,t(j)+dt/2); bhat2=(fft(b2,[],2));%b(t+dt/2) for k2&k3
    b4 =bfun(xx,yy,t(j)+dt); bhat4=(fft(b4,[],2));%b(t+dt) for k4
    % At each time step, for each Fourier Mode solve 
    % the y differentiation at 4 RK stages.
    for m=1:length(kx)
        k1=(D2-kx(m)^2*I)*uhat(:,m)+bhat1(:,m);
        k2=(D2-kx(m)^2*I)*(uhat(:,m)+(dt/2).*k1)+bhat2(:,m);
        k3=(D2-kx(m)^2*I)*(uhat(:,m)+(dt/2).*k2)+bhat2(:,m);
        k4=(D2-kx(m)^2*I)*(uhat(:,m)+dt.*k3)+bhat4(:,m);
        uhat(:,m)=uhat(:,m)+dt/6*(k1+2*k2+2*k3+k4);
    end
    u=real(ifft(uhat,[],2)); % Convert to 'real space' via inverse FFT
    % impose Dirichlet boundary counditions along y dimension
    u([1 Ny+1],:) = ufun(xx([1 Ny+1],:),yy([1 Ny+1],:),t(j+1));
    % Visualize
    if mod(j,p)==0
        clf;        
        uexx=ufun(xx,yy,t(j+1));
        subplot(131);surf(xx,yy,u);title('solved u')
        subplot(132);surf(xx,yy,uexx);title(sprintf('Exact U @ %-4.3e',t(j)))
        subplot(133);surf(xx,yy,u-uexx);title('error')
        for sp=1:3
            subplot(1,3,sp)
            shading interp;caxis([-1,1]);axis equal;view([-18,81])
            xlim([0 2*pi])
        end
%         saveas(gcf,[num2str(j),'.jpg'])
%         set(gcf,'visible','off')
    end
    
end


uex=ufun(xx,yy,t(j+1));
error=u-uex
max(abs(uex(:)-u(:)))
