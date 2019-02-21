%% Hybrid-Spectral Fourier-Galerkin/Chebyshev-Collocation Staggered NS 2d Solver

%% Clean
clc;clear all;close all;warning('off');

%% Input
Nx = 8;
Ny = 6;
Re = 1;

%% Preprocessing
% Parameters
vis = 1/Re;                % visc 1/Re

% Grids
x   = (0:Nx-1)/Nx*2*pi;          % x coordinate
kx  = fftshift(-Nx/2:Nx/2-1);    % wave number vector
ygl = -cos(pi*(0:Ny)/Ny)';
yg  = -cos(pi*(0.5:Ny)/Ny)';
[X,Y]   = meshgrid(x,ygl);
[Xg,Yg] = meshgrid(x,yg);
% plot(ones(size(ygl)),ygl,'r*',ones(size(yg)),yg,'k*');legend('G.L.','G.')

% Vandermonde Matrices for Gauss-Lobatto & Gauss
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

% Time
tfinal = 1e2;
dt     = 1e-3;
Nstep  = ceil(tfinal/dt);
t      = (0:Nstep)*dt;

% Other operators
ZNx = diag(ones(1,Ny+1));
% ZNx([1 end],:) = D([1 end],:);
ZNy = diag([0 ones(1,Ny-1) 0]);
Igl = speye(Ny+1);
Ig  = speye(Ny);

% Divergence and gradient operators
VG_trunc = [VG zeros(Ny,1)];
div_x  = VG_trunc*inv(VGL)*Igl*ZNx;   % must be multiplied with 1i*kx(m) for each Fourier mode
div_y  = VG_trunc*inv(VGL)*Dgl*ZNy;
grad_x = VGL(:,1:Ny)*inv(VG)*Ig;      % must be multiplied with 1i*kx(m) for each Fourier mode
grad_y = VGL(:,1:Ny)*inv(VG)*Dg;
div_x_act_on_grad_x = -Ig;            % must be multiplied with kx(m)^2 for each Fourier mode
div_y_act_on_grad_y = div_y*grad_y;

% Test Functions U,V,P
ufun = @(x,y,t) ones(size(X));
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

%%  Time Marching.
tic
for tt=1:Nstep;
    uh = fft(u,[],2);
    vh = fft(v,[],2);
    % Convective terms
    [ududxh]=aapx(uh,uh*(diag(kx*1i)));
    [vdudyh]=aapx(vh,D*uh);
    [udvdxh]=aapx(uh,vh*(diag(kx*1i)));
    [vdvdyh]=aapx(vh,D*vh);
    % Loop over Fourier modes
    for m=1:length(kx);
        % Solving for intermediate velocities; ustar & vstar (eq.7.139-7.142)
        A = vis*D2-(1/dt)*Igl-(vis*kx(m).^2)*Igl;
        convectiveu = ududxh(:,m)+vdudyh(:,m);
        rhsuint     = -(1/dt)*uh(:,m)+convectiveu+fh(:,m)+th(:,m);
        ustar       = A\rhsuint;
        ustar([1 Ny+1]) = 0;
        convectivev = udvdxh(:,m)+vdvdyh(:,m);
        rhsvint     = -(1/dt)*vh(:,m)+convectivev;
        vstar       = A\rhsvint;
        vstar([1 Ny+1])=0;
        % solve for pressure L*q=f;
        RHS = div_x*(kx(m)*1i)*ustar + div_y*vstar;
        LHS = dt*( div_x_act_on_grad_x*(kx(m))^2 + div_y_act_on_grad_y );
        pnewh(:,m) = LHS\RHS;
        % Update velocities
        unewh(:,m) = ustar-dt*grad_x*(kx(m)*1i)*pnewh(:,m);
        vnewh(:,m) = vstar-dt*grad_y*pnewh(:,m);
contcont(:,m)=div_x*(kx(m)*1i)*unewh(:,m)+div_y*vnewh(:,m);

    end
    unew = real(ifft(unewh,[],2));
    vnew = real(ifft(vnewh,[],2));
    u=unew;v=vnew;
    % Enforce BC
%     u=[zeros(1,Nx);  unew(2:Ny,:); zeros(1,Nx)];
%     v=[zeros(1,Nx);  vnew(2:Ny,:); zeros(1,Nx)];
     rescont=sum(contcont(:));
   if rescont>1e-12 ;disp('Cont not satisfied');break;end

end
toc