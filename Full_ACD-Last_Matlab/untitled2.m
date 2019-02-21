clc;clear all;
load CENTERS.mat
Nx=Nx+1;
Ny=Ny+1;
Nz=Nz+1;
dx=Lx/Nx;

eta_y=2/Ly;
eta_gl=-cos(pi*(0:Ny)/Ny)';
y_gl=(eta_gl+1)/eta_y;

xi  = (0:Nx-1)/Nx*2*pi; xi_x = 2*pi/Lx; x   = xi/xi_x;

zi  = (0:Nz-1)/Nz*2*pi; zi_z = 2*pi/Lz; z   = zi/zi_z;
[X,Y,Z]   = meshgrid(x,y_gl,z);

close all
figure
hold on
plot(diff(y_gl))
plot(diff(x),'r')
plot(diff(z),'g')


% Gauss smearing

epi=1.1; ep=epi*dx; ep2=ep*ep;
%-----
% 3D
%-----
% f_eps=zeros(length(CENTERS),Ny,Nx,Nz);
gdim=3; pieps=1.0/(ep^gdim*pi^(gdim/2));
for kkk=1:length(CENTERS)
% x=zeros(ny+1,nx,nz); y=x; z=x;  % vertex coords
x0=CENTERS(kkk,1);y0=CENTERS(kkk,2);z0=CENTERS(kkk,3);
% x0=0.5; y0=0.5; z0=0.5; % arbitrary centre 
x=X;y=Y;z=Z;
ny=Ny;nz=Nz;nx=Nx;
%%
for k=2:nz-1; for j=2:nx-1;  for i=2:ny; 
            dvol(i,j,k)=(y(i+1,j-1,k+1)-y(i-1,j-1,k+1))...
                       *(x(i+1,j+1,k+1)-x(i+1,j-1,k+1))...
                       *(z(i+1,j-1,k+1)-z(i+1,j-1,k-1)); end; end; end; %use Jacobi in general case...
%%
sumvol=0;
for k=1:nz-1; for j=1:nx-1; for i=1:ny; sumvol=sumvol+dvol(i,j,k); end; end; end; sumvol;
%% center
sumf3d=0.0;
for k=2:nz-1; for j=2:nx-1;  for i=2:ny; 
    p=(x(i,j,k)-x0)^2+(y(i,j,k)-y0)^2+(z(i,j,k)-z0)^2;
    f(i,j,k)=exp(-p/ep2);
    sumf3d=sumf3d+f(i,j,k)*dvol(i,j,k)/8; end; end; end;
% smeared loading in 3d cfd, f_eps = f*dvol*pieps*g where g is disc/line loading at point(x0,y0,z0) 
sumf3d=sumf3d*pieps
err3d(kkk)=sumf3d-1.0
f_eps{kkk} = f.*dvol/8.*pieps;


end
save f_eps f_eps
