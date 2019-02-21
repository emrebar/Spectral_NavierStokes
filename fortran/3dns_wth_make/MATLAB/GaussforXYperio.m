clc;clear all;
load CENTERS.mat
Nx=Nx+1;Ny=Ny+1;Nz=Nz+1;
dx=Lx/Nx;
eta_z=2/Lz;eta_gl=-cos(pi*(0:Nz)/Nz)';z_gl=(eta_gl+1)/eta_z;
xi  = (0:Nx-1)/Nx*2*pi; xi_x = 2*pi/Lx; x   = xi/xi_x;
yi  = (0:Ny-1)/Ny*2*pi; yi_y = 2*pi/Ly; y   = yi/yi_y;
[X,Y,Z]   = meshgrid(x,y,z_gl);

% close all; plot(diff(z_gl));hold on;plot(diff(x),'r');plot(diff(y),'g');

% Gauss smearing
epi=1; ep=epi*dx; ep2=ep*ep;
%----- 3D -----%
gdim=3; pieps=1.0/(ep^gdim*pi^(gdim/2));
x=X;y=Y;z=Z;  ny=Ny;nz=Nz;nx=Nx;
for kkk=1:length(CENTERS)
% x=zeros(ny+1,nx,nz); y=x; z=x;  % vertex coords
x0=CENTERS(kkk,1);y0=CENTERS(kkk,2);z0=CENTERS(kkk,3);% arbitrary centre

%%
for k=2:nz; for j=2:nx-1;  for i=2:ny-1; 
            dvol(i,j,k)=(y(i+1,j-1,k+1)-y(i-1,j-1,k+1))...
                       *(x(i+1,j+1,k+1)-x(i+1,j-1,k+1))...
                       *(z(i+1,j-1,k+1)-z(i+1,j-1,k-1)); end; end; end; %use Jacobi in general case...
%%
sumvol=0;
for k=1:nz; for j=1:nx-1; for i=1:ny-1; sumvol=sumvol+dvol(i,j,k); end; end; end; sumvol;
%% center
sumf3d=0.0;
for k=2:nz; for j=2:nx-1;  for i=2:ny-1; 
    p=(x(i,j,k)-x0)^2+(y(i,j,k)-y0)^2+(z(i,j,k)-z0)^2;
    f(i,j,k)=exp(-p/ep2);
    sumf3d=sumf3d+f(i,j,k)*dvol(i,j,k)/8; end; end; end;
% smeared loading in 3d cfd, f_eps = f*dvol*pieps*g where g is disc/line loading at point(x0,y0,z0) 
sumf3d=sumf3d*pieps
err3d(kkk)=sumf3d-1.0
f_eps{kkk} = f.*dvol/8.*pieps;


end
save f_eps f_eps
