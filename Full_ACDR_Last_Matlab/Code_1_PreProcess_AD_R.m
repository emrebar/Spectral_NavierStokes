clc;clear all;close all;
Nx =64;Ny=56;Nz=64; % number of harmonics
Lx= 30;Ly=4;Lz=6;   % mapping ''areas''
eta_ygl  = 2/Ly;   etagl = -cos(pi*(0:Ny)/Ny)';   ygl = (etagl+1)/eta_ygl;
eta_yg  = 2/Ly;    etag  = -cos(pi*(0.5:Ny)/Ny)'; yg = (etag+1)/eta_yg;
xi  = (0:Nx-1)/Nx*2*pi; xi_x = 2*pi/Lx; x   = xi/xi_x;dx=Lx/Nx;
zi  = (0:Nz-1)/Nz*2*pi; zi_z = 2*pi/Lz; z   = zi/zi_z;
dy=diff(ygl);dz=diff(z);vol=max(dy)*dx(1)*dz(1)

%% Modes
kx  = fftshift(-Nx/2:Nx/2-1);    % wave number vector
kz = fftshift(-Nz/2:Nz/2-1);
kzz=shiftdim(kz,-1);
%% Preprocessing
% Parametergausslast
Re=500;vis = 1/Re;                % visc 1/Re
[X,Y,Z]   = meshgrid(x,ygl,z);
[Xg,Yg,Zg] = meshgrid(x,yg,z);
% Vandermonde Matrices for Gauss-Lobatto & Gauss
Vetagl          = cos(acos(etagl)*(0:Ny));
dVetagl         = diag(1./sqrt(1-etagl.^2))*sin(acos(etagl)*(0:Ny))*diag(0:Ny);
dVetagl(1,:)    = (-1).^(1:Ny+1).*(0:Ny).^2;
dVetagl(Ny+1,:)  = (0:Ny).^2;
Detagl          = dVetagl/Vetagl;
D2etagl = Detagl^2;

Vetag  = cos(acos(etag)*(0:Ny-1));
Detag          = diag(1./sqrt(1-etag.^2))*sin(acos(etag)*(0:Ny-1))*diag(0:Ny-1);
Detag  = Detag/Vetag;
% Time
tfinal =200;  dt     = 1e-2/2.5;  Nstep  = ceil(tfinal/dt);  t      = (0:Nstep)*dt;
% Other operators
ZNx = diag(ones(1,Ny+1));
ZNz = ZNx;
ZNy = diag([0 ones(1,Ny-1) 0]);
Igl = speye(Ny+1);
Ig  = speye(Ny);
Igldt=(3/(2*dt))*Igl;
Iglvis=vis*Igl;
% Divergence and gradient operators
VG_trunc = [Vetag zeros(Ny,1)];
div_x  = VG_trunc*inv(Vetagl)*Igl*ZNx;   % must be multiplied with 1i*kx(oo) for each Fourier mode
div_y  = VG_trunc*inv(Vetagl)*Detagl*eta_ygl*ZNy;
grad_x = Vetagl(:,1:Ny)*inv(Vetag)*Ig;      % must be multiplied with 1i*kx(oo) for each Fourier mode
grad_y = Vetagl(:,1:Ny)*inv(Vetag)*Detag*eta_yg;
div_x_act_on_grad_x = -Ig;            % must be multiplied with kx(oo)^2 for each Fourier mode
div_y_act_on_grad_y = div_y*grad_y;
Detagleta_ygl=Detagl*eta_ygl;
operate=Igldt-D2etagl*eta_ygl^2*vis;
save preprocess