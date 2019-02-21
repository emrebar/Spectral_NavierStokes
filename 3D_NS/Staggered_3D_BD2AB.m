%% 3D Navier Stokes,x&z periodic B.C. y non periodic
%% Clean
tic
clc;clear all;close all;warning('off');
% profile on
%% Grid,Coordinates 
Nx =4;Ny =4;Nz=4;        % number of harmonics/modes
x  = 2*pi/Nx*(0:Nx-1);      % x coordinate
z  = 2*pi/Nz*(0:Nz-1);      % z coordinate
%% Parameters
Re=2;vis = 1/Re;            % Non-Dimensionalized Viscosity=1/Re
tol=1e-6;                   % tolerance-> Steady Check

% FFT Modes in X and Z
kx  = fftshift(-Nx/2:Nx/2-1);% wave number vector X
kz = fftshift(-Nz/2:Nz/2-1); % wave number vector Z
kzz=shiftdim(kz,-1);         % wave number vector Z adjustment
%% Grid Staggerin in y direction
% Gauss (Pressure)  & Gauss Lobatto (Velocity) Points 
% vertical y direction is staggered grid.
% The momentum equation is enforced at the Gauss�Lobatto points where the
% velocities are evaluated, and time stepped.
% whereas the continuity equations is enforced at the Gauss points/
% pressure nodes
ygl = -cos(pi*(0:Ny)/Ny)';
yg  = -cos(pi*(0.5:Ny)/Ny)';
%% Generate 3D modes
[X,Y,Z]   = meshgrid(x,ygl,z);
[Xg,Yg,Zg] = meshgrid(x,yg,z);
%% Preprocessing for Chebyshev differentiation matrices 
% they are used in y direction where the pressure nodes are on the
% staggered grid
%Vandermonde Matrices for Gauss-Lobatto & Gauss
VGL = cos(acos(ygl(:))*(0:Ny));
VG  = cos(acos(yg(:))*(0:Ny-1));
dVGL         = diag(1./sqrt(1-ygl.^2))*sin(acos(ygl)*(0:Ny))*diag(0:Ny);
dVGL(1,:)    = (-1).^(1:Ny+1).*(0:Ny).^2;
dVGL(Ny+1,:) = (0:Ny).^2;
dVG          = diag(1./sqrt(1-yg.^2))*sin(acos(yg)*(0:Ny-1))*diag(0:Ny-1);
% Differentiation Matrices for  Gauss & GaussLobatto
Dg  = dVG/VG;
Dgl = dVGL/VGL;
D   = Dgl;
D2  = Dgl*Dgl;
D2v=vis*D2;

%% Time Step etc.
tfinal = 0.001;  dt     = 1e-3;  Nstep  = ceil(tfinal/dt);  t      = (0:Nstep)*dt;

%% Other operators 
% Some operators for boundary conditions.
ZNx = diag(ones(1,Ny+1));
% ZNx([1 end],:) = D([1 end],:);
ZNy = diag([0 ones(1,Ny-1) 0]);
Igl = speye(Ny+1);
Ig  = speye(Ny);
Igldt=(3/(2*dt))*Igl;
Iglvis=vis*Igl;
% Divergence and gradient operators
VG_trunc = [VG zeros(Ny,1)];
div_x  = VG_trunc*inv(VGL)*Igl*ZNx;   % must be multiplied with 1i*kx(oo) for each Fourier mode
div_y  = VG_trunc*inv(VGL)*Dgl*ZNy;
grad_x = VGL(:,1:Ny)*inv(VG)*Ig;      % must be multiplied with 1i*kx(oo) for each Fourier mode
grad_y = VGL(:,1:Ny)*inv(VG)*Dg;
div_x_act_on_grad_x = -Ig;            % must be multiplied with kx(oo)^2 for each Fourier mode
div_y_act_on_grad_y = div_y*grad_y;

% Initial Velocity Functions U,V,W and Constant Pressure gradient -P=ffun
% everywhere -2/Re
ufun = @(x,y,z,t) ((1-y.^2));vfun = @(x,y,z,t) x*0;wfun = @(x,y,z,t) y*0;
ffun = @(x,y,z,t) (-2/Re)*ones(size((x)));

% Evaluate the initial conditions,
u=ufun(X,Y,Z,0); uold=ufun(X,Y,Z,-dt);
v=vfun(X,Y,Z,0); vold=vfun(X,Y,Z,-dt);
w=wfun(X,Y,Z,0); wold=wfun(X,Y,Z,-dt);
f=ffun(X,Y,Z,0);
%% Preallocate
[unewh,vnewh,wnewh]=deal(zeros(size(u)));[pnewh,contcont]=deal(zeros(size(Yg)));

%% Transform Constant Pressure gradient force in fourier space in the horizontal directions
fh = fft(fft(f,[],3),[],2);
%% Flat Plate representation
% a 3D empty matrix with constant forces 4 by 4 mounted in it with this configuration  
turb=zeros(size(X));
turb(ceil(7*Ny/16):ceil(10*Ny/16),ceil(40*Nx/256),ceil(7*Nz/16):ceil(10*Nz/16))=1;
th = fft(fft(turb,[],3),[],2);
%% Force distribution kernel
ep=2*pi/Nx;kern=(1/(ep^3*pi^1.5));
tic
for j=1:Nstep;
    %% 
    uoldh = fft(fft(uold,[],3),[],2);    voldh = fft(fft(vold,[],3),[],2);    woldh = fft(fft(wold,[],3),[],2);
    uh = fft(fft(u,[],3),[],2); vh = fft(fft(v,[],3),[],2); wh = fft(fft(w,[],3),[],2);
    %% Exact Convective Term Functions (Comment the ones to compare)
    %% Convective Terms Approximations
    %%%%%%%%%Convective Terms%%%%%%%%%
    %%%Pseudo spectral convective terms                                   %%%
    %%%After taking the derivatives in the fourier space                  %%%
    %%%transform back to real space and evaluated the non linear terms in %%%
    %%%real space.This configuration no dealisaing                        %%%
    
    % u*dudx & v*dudx & w*dudx (at time t & t-1) - non dealiased
    dudxh= bsxfun(@times, kx*1i, uh);    dudxoldh= bsxfun(@times, kx*1i, uoldh);
    dvdxh= bsxfun(@times, kx*1i, vh);    dvdxoldh= bsxfun(@times, kx*1i, voldh);
    dwdxh= bsxfun(@times, kx*1i, wh);    dwdxoldh= bsxfun(@times, kx*1i, woldh);
    
    ududx=ifft(ifft(dudxh,[],3),[],2).*u; ududxh = fft(fft(ududx,[],3),[],2); % Non Dealiased Products
    udvdx=ifft(ifft(dvdxh,[],3),[],2).*u; udvdxh = fft(fft(udvdx,[],3),[],2);
    udwdx=ifft(ifft(dwdxh,[],3),[],2).*u; udwdxh = fft(fft(udwdx,[],3),[],2);
    ududxold=ifft(ifft(dudxoldh,[],3),[],2).*uold; ududxoldh = fft(fft(ududxold,[],3),[],2); % Non Dealiased Products
    udvdxold=ifft(ifft(dvdxoldh,[],3),[],2).*uold; udvdxoldh  = fft(fft(udvdxold,[],3),[],2);
    udwdxold=ifft(ifft(dwdxoldh,[],3),[],2).*uold; udwdxoldh  = fft(fft(udwdxold,[],3),[],2);
    
    %% w*dudz & w*dvdz & w*dwdz (at time t & t-1) - dealiased
    % CLOSED LOOP NR2
    dudzh=bsxfun(@times,kzz*1i,uh);  dudzoldh=bsxfun(@times,kzz*1i,uoldh);
    dvdzh=bsxfun(@times,kzz*1i,vh);  dvdzoldh=bsxfun(@times,kzz*1i,voldh);
    dwdzh=bsxfun(@times,kzz*1i,wh);  dwdzoldh=bsxfun(@times,kzz*1i,woldh);
    
    wdudz=ifft(ifft(dudzh,[],3),[],2).*w; wdudzh = fft(fft(wdudz,[],3),[],2);
    wdvdz=ifft(ifft(dvdzh,[],3),[],2).*w; wdvdzh = fft(fft(wdvdz,[],3),[],2);
    wdwdz=ifft(ifft(dwdzh,[],3),[],2).*w; wdwdzh = fft(fft(wdwdz,[],3),[],2);
    wdudzold=ifft(ifft(dudzoldh,[],3),[],2).*wold; wdudzoldh = fft(fft(wdudzold,[],3),[],2);
    wdvdzold=ifft(ifft(dvdzoldh,[],3),[],2).*wold; wdvdzoldh = fft(fft(wdvdzold,[],3),[],2);
    wdwdzold=ifft(ifft(dwdzoldh,[],3),[],2).*wold; wdwdzoldh = fft(fft(wdwdzold,[],3),[],2);
    %% v*dudy & v*dvdy & v*dwdy (at time t & t-1) - dealiased
    for e = 1:Nz;
        dudy(:,:,e) = D*u(:,:,e); dvdy(:,:,e) = D*v(:,:,e); dwdy(:,:,e) = D*w(:,:,e);
        dudyold(:,:,e) = D*uold(:,:,e); dvdyold(:,:,e) = D*vold(:,:,e); dwdyold(:,:,e) = D*wold(:,:,e);
    end
    vdudy=dudy.*v; vdudyh = fft(fft(vdudy,[],3),[],2);
    vdvdy=dvdy.*v; vdvdyh = fft(fft(vdvdy,[],3),[],2);
    vdwdy=dwdy.*v; vdwdyh = fft(fft(vdwdy,[],3),[],2);
    vdudyold=dudyold.*vold; vdudyoldh = fft(fft(vdudyold,[],3),[],2);
    vdvdyold=dvdyold.*vold; vdvdyoldh = fft(fft(vdvdyold,[],3),[],2);
    vdwdyold=dwdyold.*vold; vdwdyoldh = fft(fft(vdwdyold,[],3),[],2);
    
    %% Loop over Fourier modes 
    % For each mode/k in X and Z directions solve a linear matrix system
    % with Chebyshev differentiation matrix.
    for m=1:length(kz);
        for oo=1:length(kx);
            %% Solving for intermediate velocities; ustar & vstar & wstar
            % Backward differencing and Adams Bashfort 2nd order time
            % stepping scheme...
            
            % Implicit Part, viscous terms and unew.LHS-->A matrix
            A=Igldt-D2v+(kx(oo)^2+kz(m)^2)*Iglvis;
            %% Construct RHSs for U,V,W 
            % Convective U where the forces for turbine (th) and constant
            % pressure gradient is fed in. 
            convectiveu=2*ududxh(:,oo,m)+2*vdudyh(:,oo,m)+2*wdudzh(:,oo,m)-ududxoldh(:,oo,m)-vdudyoldh(:,oo,m)-wdudzoldh(:,oo,m);
            rhsuint=(2/dt)*uh(:,oo,m)-(1/(2*dt))*uoldh(:,oo,m)+convectiveu-fh(:,oo,m)-th(:,oo,m);
            ustar       = A\rhsuint;
            ustar([1 Ny+1]) = 0;
            % Convective V
            convectivev=2*vdvdyh(:,oo,m)+2*udvdxh(:,oo,m)+2*wdvdzh(:,oo,m)-udvdxoldh(:,oo,m)-vdvdyoldh(:,oo,m)-wdvdzoldh(:,oo,m);
            rhsvint=(2/dt)*vh(:,oo,m)-(1/(2*dt))*voldh(:,oo,m)+convectivev;
            vstar       = A\rhsvint;
            vstar([1 Ny+1])=0;
            % Convective W
            convectivew=2*vdwdyh(:,oo,m)+2*udwdxh(:,oo,m)+2*wdwdzh(:,oo,m)-udwdxoldh(:,oo,m)-vdwdyoldh(:,oo,m)-wdwdzoldh(:,oo,m);
            rhswint=(2/dt)*wh(:,oo,m)-(1/(2*dt))*woldh(:,oo,m)+convectivew;
            wstar       = A\rhswint;
            wstar([1 Ny+1])=0;
            %% After solving for intermediate velocities 
            % solve for pressure with poisson LHS*q=RHS;   
            RHS = div_x*(kx(oo)*1i)*ustar + div_y*vstar +  div_x*(kz(m)*1i)*wstar;
            LHS = dt*( div_x_act_on_grad_x*(kx(oo)^2+kz(m)^2) + div_y_act_on_grad_y );
            
            pnewh(:,oo,m) = LHS\RHS;
            % Update velocities
            unewh(:,oo,m) = ustar-dt*grad_x*(kx(oo)*1i)*pnewh(:,oo,m);
            vnewh(:,oo,m) = vstar-dt*grad_y*pnewh(:,oo,m);
            wnewh(:,oo,m) = wstar-dt*grad_x*(kz(m)*1i)*pnewh(:,oo,m);
            
            %% Take the divergence of the new velocities and store 
            % later on in the code to make sure that it is satisfied
            contcont(:,oo,m)=div_x*(kx(oo)*1i)*unewh(:,oo,m)+div_y*vnewh(:,oo,m)+div_x*(kz(m)*1i)*wnewh(:,oo,m);
            
            
        end
    end
    % Transform to real space to check the steadiness and impose Boundary
    % conditions This can be done in fourier space as well.
    unew = real(ifft(ifft(unewh,[],3),[],2));
    vnew = real(ifft(ifft(vnewh,[],3),[],2));
    wnew = real(ifft(ifft(wnewh,[],3),[],2));
    % Enforce BC
    u([1 Ny+1],:) = ufun(X([1 Ny+1],:),Y([1 Ny+1],:),Z([1 Ny+1],:),t(j+1));    % B.C.
    v([1 Ny+1],:) = vfun(X([1 Ny+1],:),Y([1 Ny+1],:),Z([1 Ny+1],:),t(j+1));    % B.C.
    w([1 Ny+1],:) = wfun(X([1 Ny+1],:),Y([1 Ny+1],:),Z([1 Ny+1],:),t(j+1));    % B.C.
    
    %% Check if it has reached steady solution
    steady=(1/(2*dt))*(3*unew-4*u+uold);
    if max(abs(steady(:)))<tol; disp('solution stedy-tol 1e-6');break;end
    
    %% Check if the continuity is actually satisfied
    rescont=sum(contcont(:));
    if rescont>1e-10 ;disp('Cont not satisfied');break;end
    
    %% Transform to real space for the old time step velocity 
    % required for time scheme
    uold = real(ifft(ifft(uh,[],3),[],2));
    vold = real(ifft(ifft(vh,[],3),[],2));
    wold = real(ifft(ifft(wh,[],3),[],2));
    
    u=unew;v=vnew;w=wnew;
    
if mod(j,10000)==0;disp(num2str(j/Nstep)); filename=['Nx',num2str(Nx),'Ny',num2str(Ny),'Nz',num2str(Nz),'Re',num2str(Re),'dt',num2str(dt),'T',num2str(t(j)),'BDF2AB2nomap.mat'];	save(filename,'Nx','Ny','Nz','Re','t','j','Nstep','u','v','w','X','Y','Z','turb', 'dt','rescont');end

end
toc