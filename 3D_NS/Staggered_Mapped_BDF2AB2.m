%% Clean
clc;clear all;close all;warning('off');
% profile on
%% Grid,Coordinates
tic
Nx =192;Ny=42;Nz=64; % number of harmonics
Lx= 30;Ly=4;Lz=6;   % mapping ''areas''
eta_ygl  = 2/Ly;   etagl = -cos(pi*(0:Ny)/Ny)';   ygl = (etagl+1)/eta_ygl;
eta_yg  = 2/Ly;    etag  = -cos(pi*(0.5:Ny)/Ny)'; yg = (etag+1)/eta_yg;
xi  = (0:Nx-1)/Nx*2*pi; xi_x = 2*pi/Lx; x   = xi/xi_x;
zi  = (0:Nz-1)/Nz*2*pi; zi_z = 2*pi/Lz; z   = zi/zi_z;
%% JUST FOR INITIAL
% y  = -cos(pi*(0:Ny)/Ny)'; % y coordinate
% [XX,YY,ZZ]   = meshgrid(x,y,z);
%% Modes

kx  = fftshift(-Nx/2:Nx/2-1);    % wave number vector
kz = fftshift(-Nz/2:Nz/2-1);
kzz=shiftdim(kz,-1);
%% Preprocessing
% Parameters
Re=10;vis = 1/Re;                % visc 1/Re

[X,Y,Z]   = meshgrid(x,ygl,z);
[Xg,Yg,Zg] = meshgrid(x,yg,z);
% plot(ones(size(ygl)),ygl,'r*',ones(size(yg)),yg,'k*');legend('G.L.','G.')

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
tfinal = 0.001;  dt     = 1e-3;  Nstep  = ceil(tfinal/dt);  t      = (0:Nstep)*dt;

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

% Test Functions U,V,P
% ufun = @(x,y,z,t) Ly-(y-Ly/2).^2;
ufun = @(x,y,z,t) zeros(size(x));
vfun = @(x,y,z,t) x*0;
wfun = @(x,y,z,t) y*0;
ffun = @(x,y,z,t) (-2/Re)*ones(size((x)));
% Initial conditions,
u=ufun(X,Y,Z,0); uold=ufun(X,Y,Z,-dt);
v=vfun(X,Y,Z,0); vold=vfun(X,Y,Z,-dt);
w=wfun(X,Y,Z,0); wold=wfun(X,Y,Z,-dt);
f=ffun(X,Y,Z,0);
[unewh,vnewh,wnewh,dudy,dvdy,dwdy,dudyold,dvdyold,dwdyold]=deal(zeros(size(u)));[pnewh,contcont]=deal(zeros(size(Yg)));
% load stable_Nx8Ny6Nz8.mat;u=stable_Nx8Ny6Nz8;uold=u;
% Initial conditions, fft'ed
fh = fft(fft(f,[],3),[],2);
% %% Force Distribution
% ep=x(2)-x(1);kern=(1/(ep^3*pi^1.5));Ct=1;rad=0.5;
% Thrust=((max(u(:))^2)/2)*Ct*pi*rad^2;
% Totthrust=dot(Thrust,kern);
%% Force plug in
turb=zeros(size(X));
turb([Ny/2+1],ceil(5*Lx/30),ceil(15*Nz/30)+1)=10;


% %120 36 30
turb=zeros(size(X));
turb([Ny/2:Ny/2+2],ceil(5*Lx/30),(14*Nz/30:18*Nz/30))=2.28;
turb([Ny/2-1 Ny/2+3],ceil(5*Lx/30),(15*Nz/30:17*Nz/30))=2.28;
% %32^3
% turb=zeros(size(X));
% turb([Ny/2:Ny/2+2],ceil(7*Lx/30),(15*Nz/32:19*Nz/32))=10;
% turb([Ny/2-1 Ny/2+3],ceil(7*Lx/30),[16*Nz/32:18*Nz/32])=10;
% %64&32&64
% turb=zeros(size(X));
% turb(Ny/2+1,(31*Nx/64):(34*Nx/64),(28*Nz/64:38*Nz/64))=27;
% turb([Ny/2 Ny/2+2],(31*Nx/64):(34*Nx/64),(29*Nz/64:37*Nz/64))=27;
% turb([Ny/2-1 Ny/2+3],(31*Nx/64):(34*Nx/64),(30*Nz/64:36*Nz/64))=27;
% %16&32&16
% turb=zeros(size(X));
% turb(Ny/2+1,ceil(7*Nx/40),(8*Nz/16:10*Nz/16))=100;
% turb([Ny/2 Ny/2+2],ceil(7*Nx/40),(9*Nz/16:9*Nz/16))=100;
th = fft(fft(turb,[],3),[],2);
tol=1e-6;     % tolerance-> steady & continuity




% tic
for j=1:Nstep;
    %     while continuity>tol
    uh = fft(fft(u,[],3),[],2); vh = fft(fft(v,[],3),[],2); wh = fft(fft(w,[],3),[],2);
    uoldh = fft(fft(uold,[],3),[],2);    voldh = fft(fft(vold,[],3),[],2);    woldh = fft(fft(wold,[],3),[],2);
    %% Exact Convective Term Functions (Comment the ones to compare)
    %% Convective Terms Approximations
    %%%%%%%%%Convective Terms%%%%%%%%%
    % u*dudx & v*dudx & w*dudx (at time t & t-1) - non dealiased
    dudxh= bsxfun(@times, kx*1i*xi_x, uh);    dudxoldh= bsxfun(@times, kx*1i*xi_x, uoldh);
    dvdxh= bsxfun(@times, kx*1i*xi_x, vh);    dvdxoldh= bsxfun(@times, kx*1i*xi_x, voldh);
    dwdxh= bsxfun(@times, kx*1i*xi_x, wh);    dwdxoldh= bsxfun(@times, kx*1i*xi_x, woldh);
    ududxh=aap3dx(uh,dudxh);udvdxh=aap3dx(uh,dvdxh);udwdxh=aap3dx(uh,dwdxh);
    ududxoldh=aap3dx(uoldh,dudxoldh);udvdxoldh=aap3dx(uoldh,dvdxoldh);udwdxoldh=aap3dx(uoldh,dwdxoldh);
    
    %% w*dudz & w*dvdz & w*dwdz (at time t & t-1) - dealiased
    % CLOSED LOOP NR2
    dudzh=bsxfun(@times,kzz*1i*zi_z,uh);  dudzoldh=bsxfun(@times,kzz*1i*zi_z,uoldh);
    dvdzh=bsxfun(@times,kzz*1i*zi_z,vh);  dvdzoldh=bsxfun(@times,kzz*1i*zi_z,voldh);
    dwdzh=bsxfun(@times,kzz*1i*zi_z,wh);  dwdzoldh=bsxfun(@times,kzz*1i*zi_z,woldh);
    wdudzh=aap3dx(dudzh,wh);wdvdzh=aap3dx(dvdzh,wh);wdwdzh=aap3dx(wh,dwdzh);
    wdudzoldh=aap3dx(dudzoldh,woldh);wdvdzoldh=aap3dx(dvdzoldh,woldh);wdwdzoldh=aap3dx(woldh,dwdzoldh);
    %% v*dudy & v*dvdy & v*dwdy (at time t & t-1) - dealiased
    for e = 1:Nz;
        dudy(:,:,e) = Detagl *eta_ygl*u(:,:,e); dvdy(:,:,e) = Detagl *eta_ygl*v(:,:,e); dwdy(:,:,e) = Detagl *eta_ygl*w(:,:,e);
        dudyold(:,:,e) = Detagl *eta_ygl*uold(:,:,e); dvdyold(:,:,e) = Detagl *eta_ygl*vold(:,:,e); dwdyold(:,:,e) = Detagl *eta_ygl*wold(:,:,e);
    end
    vdudyh = aap3dx(vh,(fft(fft(dudy,[],3),[],2)));
    vdvdyh = aap3dx(vh,(fft(fft(dvdy,[],3),[],2)));
    vdwdyh = aap3dx(vh,(fft(fft(dwdy,[],3),[],2)));
    vdudyoldh = aap3dx(voldh,(fft(fft(dudyold,[],3),[],2)));
    vdvdyoldh = aap3dx(voldh,(fft(fft(dvdyold,[],3),[],2)));
    vdwdyoldh = aap3dx(voldh,(fft(fft(dwdyold,[],3),[],2)));
    
    % Loop over Fourier modes
    
    for m=1:length(kz);
        for oo=1:length(kx);
            % Solving for intermediate velocities; ustar & vsta
            A=Igldt-D2etagl*eta_ygl^2*vis+(xi_x^2*kx(oo)^2+zi_z^2*kz(m)^2)*Iglvis;
            
            convectiveu=2*ududxh(:,oo,m)+2*vdudyh(:,oo,m)+2*wdudzh(:,oo,m)-ududxoldh(:,oo,m)-vdudyoldh(:,oo,m)-wdudzoldh(:,oo,m);
            rhsuint=(2/dt)*uh(:,oo,m)-(1/(2*dt))*uoldh(:,oo,m)+convectiveu-fh(:,oo,m)-th(:,oo,m);
            ustar       = A\rhsuint;
            ustar([1 Ny+1]) = 0;
            %             ustar([1 Ny+1]) = dt*;
            convectivev=2*vdvdyh(:,oo,m)+2*udvdxh(:,oo,m)+2*wdvdzh(:,oo,m)-udvdxoldh(:,oo,m)-vdvdyoldh(:,oo,m)-wdvdzoldh(:,oo,m);
            rhsvint=(2/dt)*vh(:,oo,m)-(1/(2*dt))*voldh(:,oo,m)+convectivev;
            vstar       = A\rhsvint;
            vstar([1 Ny+1])=0;
            
            convectivew=2*vdwdyh(:,oo,m)+2*udwdxh(:,oo,m)+2*wdwdzh(:,oo,m)-udwdxoldh(:,oo,m)-vdwdyoldh(:,oo,m)-wdwdzoldh(:,oo,m);
            rhswint=(2/dt)*wh(:,oo,m)-(1/(2*dt))*woldh(:,oo,m)+convectivew;
            wstar       = A\rhswint;
            wstar([1 Ny+1])=0;
            
            % solve for pressure L*q=f;
            RHS = div_x*(kx(oo)*1i*xi_x)*ustar + div_y*vstar +  div_x*(kz(m)*1i*zi_z)*wstar;
            LHS = dt*( div_x_act_on_grad_x*(kx(oo)^2*xi_x^2+kz(m)^2*zi_z^2) + div_y_act_on_grad_y );
            
            pnewh(:,oo,m) = LHS\RHS;
            % Update velocities
            unewh(:,oo,m) = ustar-dt*grad_x*(kx(oo)*1i*xi_x)*pnewh(:,oo,m);
            vnewh(:,oo,m) = vstar-dt*grad_y*pnewh(:,oo,m);
            wnewh(:,oo,m) = wstar-dt*grad_x*(kz(m)*1i*zi_z)*pnewh(:,oo,m);
            
            contcont(:,oo,m)=div_x*(kx(oo)*1i*xi_x)*unewh(:,oo,m)+div_y*vnewh(:,oo,m)+div_x*(kz(m)*1i*zi_z)*wnewh(:,oo,m);
            
        end
    end
    unew = real(ifft(ifft(unewh,[],3),[],2));
    vnew = real(ifft(ifft(vnewh,[],3),[],2));
    wnew = real(ifft(ifft(wnewh,[],3),[],2));
    % Enforce BC
    unew([1 Ny+1],:) = ufun(X([1 Ny+1],:),Y([1 Ny+1],:),Z([1 Ny+1],:),t(j+1));    % B.C.
    vnew([1 Ny+1],:) = vfun(X([1 Ny+1],:),Y([1 Ny+1],:),Z([1 Ny+1],:),t(j+1));    % B.C.
    wnew([1 Ny+1],:) = wfun(X([1 Ny+1],:),Y([1 Ny+1],:),Z([1 Ny+1],:),t(j+1));    % B.C.
    
    steady=(1/(2*dt))*(3*unew-4*u+uold);s=max(abs(steady(:)));
    if s<tol; disp('solution stedy-tol 1e-6');
        save steadymapped;break;end
    
    %         uoldh = uh;        voldh = vh;        woldh = wh;
    uold = real(ifft(ifft(uh,[],3),[],2));
    vold = real(ifft(ifft(vh,[],3),[],2));
    wold = real(ifft(ifft(wh,[],3),[],2));
%      unew(:,14:16,:)=1;
    u=unew;v=vnew;w=wnew;
    rescont=sum(contcont(:));
    if rescont>1e-10 ;disp('Cont not satisfied');break;end
    if mod(j,10000)==0;disp(num2str(j/Nstep)); filename=['Nx',num2str(Nx),'Ny',num2str(Ny),'Nz',num2str(Nz),'Lx',num2str(Lx),'Ly',num2str(Ly),'Lz',num2str(Lz),'Re',num2str(Re),'dt',num2str(dt),'T',num2str(t(j)),'BDFAB2.mat'];	save(filename,'Nx','Ny','Nz','Re','j','Nstep','u','v','w','X','Y','Z','turb', 'dt','rescont','Lx','Ly','Lz','s');end
end

toc
sum(contcont(:))