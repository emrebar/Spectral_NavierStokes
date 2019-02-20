% function [thx]=fx3dsmear(u,v,w)

clear
load preprocess;load turbines
Lx=30;Ly=4;Lz=6;
dx=Lx/Nx;
ny=Ny+1;nz=Nz+1;nx=Nx+1;

eta_y=2/Ly;eta_gl=-cos(pi*(0:ny)/ny)';y_gl3d=(eta_gl+1)/eta_y;
xi3d  = (0:nx-1)/nx*2*pi; xi_x = 2*pi/Lx; x3d   = xi3d/xi_x;
zi3d = (0:nz-1)/nz*2*pi; zi_z = 2*pi/Lz; z3d   = zi3d/zi_z;
[X3d,Y3d,Z3d]   = meshgrid(x3d,y_gl3d,z3d);

%% Gauss smearing parameters
 ep2=ep*ep;
gdim=3; pieps=1.0/(ep^gdim*pi^(gdim/2));
%% Smear to each node
for k=2:nz-1; for j=2:nx-1;  for i=2:ny;
            dvol(i,j,k)=(Y3d(i+1,j-1,k+1)-Y3d(i-1,j-1,k+1))...
                *(X3d(i+1,j+1,k+1)-X3d(i+1,j-1,k+1))...
                *(Z3d(i+1,j-1,k+1)-Z3d(i+1,j-1,k-1)); end; end; end; %use Jacobi in general case...
%% Total volume;
sumvol=0;for k=1:nz-1; for j=1:nx-1; for i=1:ny; sumvol=sumvol+dvol(i,j,k); end; end; end; sumvol;
for kkk=1:length(CENTERS)
    x0=CENTERS(kkk,1);y0=CENTERS(kkk,2);z0=CENTERS(kkk,3);
    %% Volume per cell. 8 cell represented.
    %% Smear
    sumf3d=0.0;
    for k=2:nz-1; for j=2:nx-1;  for i=2:ny;
                p=(X3d(i,j,k)-x0)^2+(Y3d(i,j,k)-y0)^2+(Z3d(i,j,k)-z0)^2;
                f3d(i,j,k)=exp(-p/ep2);
                sumf3d=sumf3d+f3d(i,j,k)*dvol(i,j,k)/8; end; end; end;
    % smeared loading in 3d cfd, f_eps = f*dvol*pieps*g where g is disc/line loading at point(x0,y0,z0)
    sumf3d=sumf3d*pieps;
    err3d(kkk)=sumf3d-1.0
    if kkk<75;
    f_eps1{kkk} = f3d.*dvol/8.*pieps;
    else 
            f_eps2{kkk-74} = f3d.*dvol/8.*pieps;

    end
    
end

save('feps1','f_eps1')
save('feps2','f_eps2')
