%% Hybrid-Spectral Fourier-Galerkin/Chebyshev-Collocation Staggered NS 2d Solver

%% Clean
clc;clear all;close all;warning('off');

%% Input
Nx = 8;
Ny = 8;
Re = 10;
Lx= 30;Ly=4;
%% Preprocessing
% Parameters
vis = 1/Re;                % visc 1/Re

% Grids
eta_ygl  = 2/Ly;   etagl = -cos(pi*(0:Ny)/Ny)';   ygl = (etagl+1)/eta_ygl;
eta_yg  = 2/Ly;    etag  = -cos(pi*(0.5:Ny)/Ny)'; yg = (etag+1)/eta_yg;
xi  = (0:Nx-1)/Nx*2*pi; xi_x = 2*pi/Lx; x   = xi/xi_x;
[X,Y]   = meshgrid(x,ygl);
[Xg,Yg] = meshgrid(x,yg);
kx  = fftshift(-Nx/2:Nx/2-1);    % wave number vector

% Vandermonde Matrices for Gauss-Lobatto & Gauss
VGLeta  = cos(acos(etagl(:))*(0:Ny));
VGeta   = cos(acos(etag(:))*(0:Ny-1));
dVGLeta         = diag(1./sqrt(1-etagl.^2))*sin(acos(etagl)*(0:Ny))*diag(0:Ny);
dVGLeta(1,:)    = (-1).^(1:Ny+1).*(0:Ny).^2;
dVGLeta(Ny+1,:) = (0:Ny).^2;
dVGeta          = diag(1./sqrt(1-etag.^2))*sin(acos(etag)*(0:Ny-1))*diag(0:Ny-1);

% Differentiation Matrices for  Gauss & GaussLobatto
Dgeta  = dVGeta/VGeta;
Dgleta = dVGLeta/VGLeta;
Deta   = Dgleta;
D2eta  = Dgleta*Dgleta;

%% JUST FOR INITIAL
y  = -cos(pi*(0:Ny)/Ny)'; % y coordinate 
[XX,YY]   = meshgrid(x,y);
% Time
tfinal = 1e3;
dt     = 1e-3;
Nstep  = ceil(tfinal/dt);
t      = (0:Nstep)*dt;

% Other operators
ZNx = diag(ones(1,Ny+1));
ZNx(end,:) = Deta(end,:);
ZNy = diag([0 ones(1,Ny-1) 0]);

Igl = speye(Ny+1);
Ig  = speye(Ny);

% Divergence and gradient operators
VG_trunceta = [VGeta zeros(Ny,1)];
div_x  = VG_trunceta*inv(VGLeta)*Igl*ZNx;   % must be multiplied with 1i*kx(m) for each Fourier mode
div_y  = VG_trunceta*inv(VGLeta)*Dgleta*eta_ygl*ZNy;
grad_x = VGLeta(:,1:Ny)*inv(VGeta)*Ig;      % must be multiplied with 1i*kx(m) for each Fourier mode
grad_y = VGLeta(:,1:Ny)*inv(VGeta)*Dgeta;
div_x_act_on_grad_x = -Ig;            % must be multiplied with kx(m)^2 for each Fourier mode
div_y_act_on_grad_y = div_y*grad_y;

% Test Functions U,V,P
ufun = @(x,y,t) Ly-(y-Ly/2).^2;
ufun = @(x,y,t) x*0;
vfun = @(x,y,t) y*0;
ffun = @(x,y,t) (-2/Re)*ones(size((x)));

% Initial conditions, fft'ed
u  = ufun(X,Y,0);
v  = vfun(X,Y,0);
f  = ffun(X,Y,0);
fh = fft(f,[],2);
turb=zeros(size(X));
turb(ceil(3*Ny/7):ceil(5*Ny/7),ceil(4*Nx/16))=1;
th = fft(turb,[],2);
tol=1e-6;
%%  Time Marching.
tic
for tt=1:Nstep;
    uh = fft(u,[],2);
    vh = fft(v,[],2);
    % Convective terms
    [ududxh]=aapx(uh,uh*(diag(kx*1i*xi_x)));
    [vdudyh]=aapx(vh,Deta*eta_ygl*uh);
    [udvdxh]=aapx(uh,vh*(diag(kx*1i*xi_x)));
    [vdvdyh]=aapx(vh,Deta*eta_ygl*vh);
    % Loop over Fourier modes
    for m=1:length(kx);
        % Solving for intermediate velocities; ustar & vstar (eq.7.139-7.142)
        A = vis*D2eta*eta_ygl^2-(1/dt)*Igl-(vis*kx(m).^2*xi_x^2)*Igl;
        A(end,:)=Deta(end,:);
        convectiveu = ududxh(:,m)+vdudyh(:,m);
        rhsuint     = -(1/dt)*uh(:,m)+convectiveu+fh(:,m);%+th(:,m);
        rhsuint(end)=0;
        ustar       = A\rhsuint;
        ustar([1]) = 0;
        convectivev = udvdxh(:,m)+vdvdyh(:,m);
        rhsvint     = -(1/dt)*vh(:,m)+convectivev;
        A = vis*D2eta*eta_ygl^2-(1/dt)*Igl-(vis*kx(m).^2*xi_x^2)*Igl;
        vstar       = A\rhsvint;
        vstar([1 Ny+1])=0;
        % solve for pressure L*q=f;
        RHS = div_x*(kx(m)*1i*xi_x)*ustar + div_y*vstar;
        LHS = dt*( div_x_act_on_grad_x*(kx(m))^2*xi_x^2 + div_y_act_on_grad_y );
        pnewh(:,m) = LHS\RHS;
        % Update velocities
        unewh(:,m) = ustar-dt*grad_x*(kx(m)*1i*xi_x)*pnewh(:,m);
        vnewh(:,m) = vstar-dt*grad_y*pnewh(:,m);
        contcont(:,m)=div_x*(kx(m)*1i)*unewh(:,m)+div_y*vnewh(:,m);

    end
    unew = real(ifft(unewh,[],2));
    vnew = real(ifft(vnewh,[],2));
    % Enforce BC
%     unew=[zeros(1,Nx);  unew(2:Ny,:); zeros(1,Nx)];
%     vnew=[zeros(1,Nx);  vnew(2:Ny,:); zeros(1,Nx)];
   
    steady=(1/dt)*(unew-u);
    if max(abs(steady(:)))<tol; disp('solution stedy-tol 1e-6');break;end
    rescont=sum(contcont(:));
    if rescont>1e-12 ;disp('Cont not satisfied');break;end
	u=unew;v=vnew;
end
toc