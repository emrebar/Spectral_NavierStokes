%% Mapping onto better resolution for visualization
close all;clc;
set(0,'DefaultSurfaceEdgeColor','none');%set(0,'DefaultFaceColor','interp');
u=(Ly-(Y-Ly/2).^2);v=u;w=u;
info=load('./data/info.dat');
Nx=info(1);Ny=info(2);Nz=info(3);Lx=info(4) ;Ly=info(5);Lz=info(6);step=info(7);Re=info(8);dt=info(9);

%% Finer Resolution in x
Nxi = 256; % Desired Resolution in x
Lx = 30; % Change,if the simulation was mapped
xi  = (0:Nx-1)/Nx*2*pi;
xii = (0:Nxi-1)/Nxi*2*pi;
xi_x = 2*pi/Lx;
x   = xi/xi_x;
xi  = xii/xi_x;
%% Finer Resolution in z
Nyi = 200; % Desired Resolution in z
Ly = 6; % Change,if the simulation was mapped
yi  = (0:Ny-1)/Ny*2*pi;
yii = (0:Nyi-1)/Nyi*2*pi;
yi_y = 2*pi/Ly;
y   = yi/yi_y;
yi  = yii/yi_y;
%% Finer Resolution in y
Nzi     = 100;
Lz=4;
eta_z  = 2/Lz;
eta    = -cos(pi*(0:Nz)/Nz)';
etai   = linspace(-1,1,Nzi)';
z      = (eta+1)/eta_z;
zi     = (etai+1)/eta_z;

Veta          = cos(acos(eta)*(0:Nz));
Vetai         = cos(acos(etai)*(0:Nz));
Heta          = Vetai/Veta;
[Xi Yi Zi]=meshgrid(xi,yi,zi);

[ny nx nz]=size(u);
uh = fft(fft(u,[],1),[],2); vh = fft(fft(v,[],1),[],2); wh = fft(fft(w,[],1),[],2);
uhp=[uh(:,1:nx/2,:) zeros(ny,(Nxi-nx),nz) uh(:,nx/2+1:nx,:)]*Nxi/Nx; % pad uhat with zeros
vhp=[vh(:,1:nx/2,:) zeros(ny,(Nxi-nx),nz) vh(:,nx/2+1:nx,:)]*Nxi/Nx; % pad vhat with zeros
whp=[wh(:,1:nx/2,:) zeros(ny,(Nxi-nx),nz) wh(:,nx/2+1:nx,:)]*Nxi/Nx; % pad vhat with zeros

uhpn=[uhp(1:ny/2,:,:) ;zeros((Nyi-ny),Nxi,nz); uhp(ny/2+1:ny,:,:)]*Nyi/Ny; % pad uhat with yeros
vhpn=[vhp(1:ny/2,:,:) ;zeros((Nyi-ny),Nxi,nz); vhp(ny/2+1:ny,:,:)]*Nyi/Ny; % pad uhat with yeros
whpn=[whp(1:ny/2,:,:) ;zeros((Nyi-ny),Nxi,nz); whp(ny/2+1:ny,:,:)]*Nyi/Ny; % pad uhat with yeros


ui=real(ifft(ifft(uhpn,[],1),[],2));
vi=real(ifft(ifft(vhpn,[],1),[],2));
wi=real(ifft(ifft(whpn,[],1),[],2));
for kkk = 1:Nx;
    for e=1:Ny;
        uie(e,kkk,:) = Heta*squeeze(ui(e,kkk,:)); vie(e,kkk,:) = Heta*squeeze(vi(e,kkk,:));wie(e,kkk,:) = Heta*squeeze(wi(e,kkk,:));
    end;
end;
% 
% figure,slice(Xi,Yi,Zi,uie,[],[],[3]);
% colorbar;caxis([0 1])
% 
% view([0 90])
% axis equal;box on;
% 
% xlim([0 Lx]);ylim([0 Ly]);zlim([0 Lz]);
% axis equal;box on;
