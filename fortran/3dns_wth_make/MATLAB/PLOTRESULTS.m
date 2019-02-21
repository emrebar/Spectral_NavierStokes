% A program to create a plot of the computed results
% from the 3D Fortran Navier-Stokes solver
clear all; format compact; format short;
% Get coordinates
info=load('./data/info.dat');
Nx=info(1);Ny=info(2);Nz=info(3);Lx=info(4) ;Ly=info(5);Lz=info(6);step=info(7);Re=info(8);dt=info(9);
eta_zgl  = 2/Lz;   etagl = -cos(pi*(0:Nz)/Nz)';   zgl = (etagl+1)/eta_zgl;
eta_zg  = 2/Lz;    etag  = -cos(pi*(0.5:Nz)/Nz)'; zg = (etag+1)/eta_zg;
xi  = (0:Nx-1)/Nx*2*pi; xi_x = 2*pi/Lx; x   = xi/xi_x;
yi  = (0:Ny-1)/Ny*2*pi; yi_y = 2*pi/Ly; y   = yi/yi_y;
% reshape coordinates to allow easy plotting
[X,Y,Z]   = meshgrid(x,y,zgl);


%
% Open file and dataset using the default properties.
%
FILEX=['./data/u1.datbin'];
FILEY=['./data/v',num2str(9999999+i),'.datbin'];
% FILEZ=['./data/w',num2str(9999999+i),'.datbin'];
%     FILEPIC=['./data/pic',num2str(9999999+i),'.jpg'];
fid=fopen(FILEX,'r');
[fname,mode,mformat]=fopen(fid);
temp=fread(fid,Ny*Nx*Nz+1,'real*8');
u=reshape(temp,Ny,Nx,Nz+1);
fclose(fid);



