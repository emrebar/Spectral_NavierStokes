%% Clean
clc;clear all;close all;warning('off');
% profile on
tic
matlabpool
%% Grid,Coordinates
Nx =256;Ny=56;Nz=64; % number of harmonics
Lx= 30;Ly=4;Lz=6;   % mapping ''areas''
eta_ygl  = 2/Ly;   etagl = -cos(pi*(0:Ny)/Ny)';   ygl = (etagl+1)/eta_ygl;
eta_yg  = 2/Ly;    etag  = -cos(pi*(0.5:Ny)/Ny)'; yg = (etag+1)/eta_yg;
xi  = (0:Nx-1)/Nx*2*pi; xi_x = 2*pi/Lx; x   = xi/xi_x;
zi  = (0:Nz-1)/Nz*2*pi; zi_z = 2*pi/Lz; z   = zi/zi_z;
dy=diff(ygl);dx=diff(x);dz=diff(z);vol=max(dy)*dx(1)*dz(1)
%% Modes
kx  = fftshift(-Nx/2:Nx/2-1);    % wave number vector
kz = fftshift(-Nz/2:Nz/2-1);
kzz=shiftdim(kz,-1);
%% Preprocessing
% Parametergausslast
Re=1000;vis = 1/Re;                % visc 1/Re
[X,Y,Z]   = meshgrid(x,ygl,z);
[Xg,Yg,Zg] = meshgrid(x,yg,z);
% plot(ones(size(ygl)),ygl,'r*',ones(size(yg)),yg,'k*');legend('G.L.','G.')
        %%
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
% Timecomplex2real
tfinal =60;  dt     = 1e-2/2;  Nstep  = ceil(tfinal/dt);  t      = (0:Nstep)*dt;
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

% Test Functions U,V,P
% ufun = @(x,y,z,t) Ly-(y-Ly/2).^2;
ufun = @(x,y,z,t) -ones(size(x));
vfun = @(x,y,z,t) x*0;
wfun = @(x,y,z,t) y*0;
% ffun = @(x,y,z,t) (-2/Re)*ones(size((x)));
% Initial conditions,
u=ufun(X,Y,Z,0); uold=u;
v=vfun(X,Y,Z,0); vold=v;
w=wfun(X,Y,Z,0); wold=w;
% f=ffun(X,Y,Z,0);
[unewh,vnewh,wnewh,dudy,dvdy,dwdy,dudyold,dvdyold,dwdyold]=deal(zeros(size(u)));[pnewh,contcont]=deal(zeros(size(Yg)));
% Initial conditions, fft'ed
% fh = fft(fft(f,[],3),[],2);
% %% Force Distribution
%% CENTER etc
coord=[X(:) Y(:) Z(:)];
xcenter=find(coord(:,1)==x(ceil((Nx/5))));
[YY]=coord(xcenter,2);[ZZ]=coord(xcenter,3);
YY=YY-ygl(Ny/2);ZZ=ZZ-z(Nz/2);
hhh=zeros(length(YY));
 for i=1:size(YY(:));zz=ZZ(i);yy=YY(i);
  rad=norm([zz yy]);if rad<0.5;hold on;
     hhh(i)=1;end;
 end
circle=find(hhh);
CENTERS=[x(ceil(Nx/5))*ones(size(circle)) YY(circle)+ygl(Ny/2) ZZ(circle)+z(Nz/2)];
plot(YY(circle),ZZ(circle),'r*');
Ct=0.6;radius=0.5;
Thrust=((max(u(:))^2)/2)*Ct*pi*radius^2;
Thrusteach=(Thrust/vol*1.1)/length(CENTERS);

%%
save('CENTERS','CENTERS','Nx','Ny','Nz','Lx','Ly','Lz','x','ygl','z')
tol=1e-6;     % tolerance-> steady & continuity
%% Distribute FORCES
load f_eps;
% A=cell(1);
turb=0;for kkk=1:length(CENTERS);AAA{kkk}=Thrusteach*f_eps{kkk};turb=turb+AAA{kkk};end
th = fft(fft(turb,[],3),[],2);
%% Buffer logspace
buffer=floor(Nx/15);
alfa=logspace(0,1,buffer+1)/10;
Uinf=-ones(size(u(:,Nx-buffer:Nx,:)));

filen=['Nx',num2str(Nx),'Ny',num2str(Ny),'Nz',num2str(Nz),'Lx',num2str(Lx),'Ly',num2str(Ly),'Lz',num2str(Lz),'Re',num2str(Re),'dt',num2str(dt),'.mat'];	
        save(filen,'Nx','Ny','Nz','Re','X','Y','Z','turb', 'dt','Lx','Ly','Lz','Ct');


for j=1:Nstep;
   
    uh = fft(fft(u,[],3),[],2); vh = fft(fft(v,[],3),[],2); wh = fft(fft(w,[],3),[],2);
    uoldh = fft(fft(uold,[],3),[],2);    voldh = fft(fft(vold,[],3),[],2);    woldh = fft(fft(wold,[],3),[],2);
    %% Exact Convective Term Functions (Comment the ones to compare)
    %% Convective Terms Approximations
    %%%%%%%%%Convective Terms%%%%%%%%%
    % u*dudx & v*dudx & w*dudx (at time t & t-1) - non dealiased
    dudxh= bsxfun(@times, kx*1i*xi_x, uh);    dudxoldh= bsxfun(@times, kx*1i*xi_x, uoldh);
    dvdxh= bsxfun(@times, kx*1i*xi_x, vh);    dvdxoldh= bsxfun(@times, kx*1i*xi_x, voldh);
    dwdxh= bsxfun(@times, kx*1i*xi_x, wh);    dwdxoldh= bsxfun(@times, kx*1i*xi_x, woldh);
    ududx=ifft(ifft(dudxh,[],3),[],2).*u; ududxh = fft(fft(ududx,[],3),[],2); % Non Dealiased Products
    udvdx=ifft(ifft(dvdxh,[],3),[],2).*u; udvdxh = fft(fft(udvdx,[],3),[],2);
    udwdx=ifft(ifft(dwdxh,[],3),[],2).*u; udwdxh = fft(fft(udwdx,[],3),[],2);
    ududxold=ifft(ifft(dudxoldh,[],3),[],2).*uold; ududxoldh = fft(fft(ududxold,[],3),[],2); % Non Dealiased Products
    udvdxold=ifft(ifft(dvdxoldh,[],3),[],2).*uold; udvdxoldh  = fft(fft(udvdxold,[],3),[],2);
    udwdxold=ifft(ifft(dwdxoldh,[],3),[],2).*uold; udwdxoldh  = fft(fft(udwdxold,[],3),[],2);
    
    %% w*dudz & w*dvdz & w*dwdz (at time t & t-1) - dealiased
    % CLOSED LOOP NR2
    dudzh=bsxfun(@times,kzz*1i*zi_z,uh);  dudzoldh=bsxfun(@times,kzz*1i*zi_z,uoldh);
    dvdzh=bsxfun(@times,kzz*1i*zi_z,vh);  dvdzoldh=bsxfun(@times,kzz*1i*zi_z,voldh);
    dwdzh=bsxfun(@times,kzz*1i*zi_z,wh);  dwdzoldh=bsxfun(@times,kzz*1i*zi_z,woldh);
    wdudz=ifft(ifft(dudzh,[],3),[],2).*w; wdudzh = fft(fft(wdudz,[],3),[],2);
    wdvdz=ifft(ifft(dvdzh,[],3),[],2).*w; wdvdzh = fft(fft(wdvdz,[],3),[],2);
    wdwdz=ifft(ifft(dwdzh,[],3),[],2).*w; wdwdzh = fft(fft(wdwdz,[],3),[],2);
    wdudzold=ifft(ifft(dudzoldh,[],3),[],2).*wold; wdudzoldh = fft(fft(wdudzold,[],3),[],2);
    wdvdzold=ifft(ifft(dvdzoldh,[],3),[],2).*wold; wdvdzoldh = fft(fft(wdvdzold,[],3),[],2);
    wdwdzold=ifft(ifft(dwdzoldh,[],3),[],2).*wold; wdwdzoldh = fft(fft(wdwdzold,[],3),[],2);
	 %% v*dudy & v*dvdy & v*dwdy (at time t & t-1) - dealiased
    for e = 1:Nz;
        dudy(:,:,e) = Detagleta_ygl*u(:,:,e); dvdy(:,:,e) = Detagleta_ygl*v(:,:,e); dwdy(:,:,e) = Detagleta_ygl*w(:,:,e);
        dudyold(:,:,e) = Detagleta_ygl*uold(:,:,e); dvdyold(:,:,e) = Detagleta_ygl*vold(:,:,e); dwdyold(:,:,e) = Detagleta_ygl*wold(:,:,e);
    end
    vdudy=dudy.*v; vdudyh = fft(fft(vdudy,[],3),[],2);
    vdvdy=dvdy.*v; vdvdyh = fft(fft(vdvdy,[],3),[],2);
    vdwdy=dwdy.*v; vdwdyh = fft(fft(vdwdy,[],3),[],2);
    vdudyold=dudyold.*vold; vdudyoldh = fft(fft(vdudyold,[],3),[],2);
    vdvdyold=dvdyold.*vold; vdvdyoldh = fft(fft(vdvdyold,[],3),[],2);
    vdwdyold=dwdyold.*vold; vdwdyoldh = fft(fft(vdwdyold,[],3),[],2);
    
    % Loop over Fourier modes

    parfor m=1:Nz;
        for oo=1:Nx;
            warning('off');
            % Solving for intermediate velocities; ustar & vsta
            A=operate+(xi_x^2*kx(oo)^2+kz(m)^2*zi_z^2)*Iglvis;
            
            convectiveu=2*ududxh(:,oo,m)+2*vdudyh(:,oo,m)+2*wdudzh(:,oo,m)-ududxoldh(:,oo,m)-vdudyoldh(:,oo,m)-wdudzoldh(:,oo,m);
            rhsuint=(2/dt)*uh(:,oo,m)-(1/(2*dt))*uoldh(:,oo,m)+convectiveu+th(:,oo,m);%-fh(:,oo,m);%;
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
            if j>8000
            contcont(:,oo,m)=div_x*(kx(oo)*1i*xi_x)*unewh(:,oo,m)+div_y*vnewh(:,oo,m)+div_x*(kz(m)*1i*zi_z)*wnewh(:,oo,m);
            end
        end
    end
%     t_foruruer=toc
    unew = real(ifft(ifft(unewh,[],3),[],2));
    vnew = real(ifft(ifft(vnewh,[],3),[],2));
    wnew = real(ifft(ifft(wnewh,[],3),[],2));
    % Enforce BC
    unew([1 Ny+1],:) = -1;    % B.C.
    vnew([1 Ny+1],:) = 0;    % B.C.
    wnew([1 Ny+1],:) = 0;  % B.C.
    

    uold = real(ifft(ifft(uh,[],3),[],2));
    vold = real(ifft(ifft(vh,[],3),[],2));
    wold = real(ifft(ifft(wh,[],3),[],2));
%      for i=1:length(uout);unew(:,Nx-buffer+i-1,:)=unew(:,Nx-buffer+i-1,:)*uout(i)+Uinf(:,i,:)*uin(i);end
for i=1:buffer+1;unew(:,Nx-buffer+i-1,:)=unew(:,Nx-buffer+i-1,:)*(1-alfa(i))+Uinf(:,i,:)*alfa(i);end


    u=unew;v=vnew;w=wnew;
    rescont=sum(contcont(:));
    if j>8000; 
    if rescont>1e-10 ;disp('Cont not satisfied');save('asd');break;end
    if mod(j,1000)==0;disp(num2str(j/Nstep)); filename=['Nx',num2str(Nx),'Ny',num2str(Ny),'Nz',num2str(Nz),'Lx',num2str(Lx),'Ly',num2str(Ly),'Lz',num2str(Lz),'Re',num2str(Re),'dt',num2str(dt),'T',num2str(t(j)),'BDFAB2buffered.mat'];	
        save(filename,'j','u','v','w','rescont');end
    end
end

 twhole_codett=toc
