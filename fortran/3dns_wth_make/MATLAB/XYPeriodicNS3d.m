%% Clean
clc;clear all;close all;warning('off');
% profile on
tic

%% Grid,Coordinates
Nx =128;Ny=32;Nz=28; % number of harmonics
Lx= 30;Ly=6;Lz=4;   % mapping ''areas''
eta_zgl  = 2/Lz;   etagl = -cos(pi*(0:Nz)/Nz)';   zgl = (etagl+1)/eta_zgl;
eta_zg  = 2/Lz;    etag  = -cos(pi*(0.5:Nz)/Nz)'; zg = (etag+1)/eta_zg;
xi  = (0:Nx-1)/Nx*2*pi; xi_x = 2*pi/Lx; x   = xi/xi_x;
yi  = (0:Ny-1)/Ny*2*pi; yi_y = 2*pi/Ly; y   = yi/yi_y;
%% Modes
kx  = fftshift(-Nx/2:Nx/2-1);    % wave number vector
ky = fftshift(-Ny/2:Ny/2-1);
%% Preprocessing
% Parameters
Re=10;vis = 1/Re;                % visc 1/Re
[X,Y,Z]   = meshgrid(x,y,zgl);
[Xg,Yg,Zg] = meshgrid(x,y,zg);

% Vandermonde Matrices for Gauss-Lobatto & Gauss
Vetagl          = cos(acos(etagl)*(0:Nz));
dVetagl         = diag(1./sqrt(1-etagl.^2))*sin(acos(etagl)*(0:Nz))*diag(0:Nz);
dVetagl(1,:)    = (-1).^(1:Nz+1).*(0:Nz).^2;
dVetagl(Nz+1,:)  = (0:Nz).^2;
Detagl          = dVetagl/Vetagl;
D2etagl = Detagl^2;

Vetag  = cos(acos(etag)*(0:Nz-1));
Detag          = diag(1./sqrt(1-etag.^2))*sin(acos(etag)*(0:Nz-1))*diag(0:Nz-1);
Detag  = Detag/Vetag;

% Time
tfinal = 0.01;  dt     = 1e-3;  Nstep  = ceil(tfinal/dt);  t      = (0:Nstep)*dt;

% Other operators
ZNx = diag(ones(1,Nz+1));
ZNz = diag([0 ones(1,Nz-1) 0]);
Igl = speye(Nz+1);
Ig  = speye(Nz);
Igldt=(3/(2*dt))*Igl;
Iglvis=vis*Igl;
% Divergence and gradient operators
VG_trunc = [Vetag zeros(Nz,1)];
div_x  = VG_trunc*inv(Vetagl)*Igl*ZNx;   % must be multiplied with 1i*kx(oo) for each Fourier mode
div_z  = VG_trunc*inv(Vetagl)*Detagl*eta_zgl*ZNz;
grad_x = Vetagl(:,1:Nz)*inv(Vetag)*Ig;      % must be multiplied with 1i*kx(oo) for each Fourier mode
grad_z = Vetagl(:,1:Nz)*inv(Vetag)*Detag*eta_zg;
div_x_act_on_grad_x = -Ig;            % must be multiplied with kx(oo)^2 for each Fourier mode
div_z_act_on_grad_z = div_z*grad_z;

% Test Functions U,V,P
ufun = @(x,y,z,t) ones(size(x));
%  ufun = @(x,y,z,t) Lz-(z-Lz/2).^2;
vfun = @(x,y,z,t) x*0;
wfun = @(x,y,z,t) y*0;
ffun = @(x,y,z,t) (-2/Re)*ones(size((x)));
% Initial conditions,
u=ufun(X,Y,Z,0);uold=u;v=vfun(X,Y,Z,0);vold=v;w=wfun(X,Y,Z,0); wold=w;
% f=ffun(X,Y,Z,0);

% [unewh,vnewh,wnewh,dudy,dvdy,dwdy,dudyold,dvdyold,dwdyold]=deal(zeros(size(u)));[pnewh,contcont]=deal(zeros(size(Yg)));
% fh = fft(fft(f,[],1),[],2);
tol=1e-6;
%% CENTER etc
coord=[X(:) Y(:) Z(:)];
xcenter=find(coord(:,1)==x(ceil((Nx/5))));
[YY]=coord(xcenter,2);[ZZ]=coord(xcenter,3);
YY=YY-y(Ny/2);ZZ=ZZ-zgl(Nz/2);
hhh=zeros(length(YY));
 for i=1:size(YY(:));
     zz=ZZ(i);yy=YY(i);
     rad=norm([zz yy]);
     if rad<0.5; hhh(i)=1;end;
 end
circle=find(hhh);
CENTERS=[x(ceil(Nx/5))*ones(size(circle)) YY(circle)+y(Ny/2) ZZ(circle)+zgl(Nz/2)];
% plot(YY(circle),ZZ(circle),'r*');
Ct=1;radius=0.5;
Thrust=((max(u(:))^2)/2)*Ct*pi*radius^2;
Thrusteach=Thrust/length(CENTERS);

save('CENTERS','CENTERS','Nx','Ny','Nz','Lx','Ly','Lz','x','y','zgl')
%% Distribute FORCES
load f_eps;
turb=0;for kkk=1:length(CENTERS);AAA{kkk}=Thrusteach*f_eps{kkk};turb=turb+AAA{kkk};end
th = fft(fft(turb,[],1),[],2);
%% Save turbine forces for fortran 
nbytes = Ny*Nx*(Nz+1)*4;fid = fopen('turbine.dat','w');
fwrite(fid, nbytes, 'int32');fwrite(fid, turb, 'single');fwrite(fid, nbytes, 'int32');fclose(fid);
%% Bufferzone cos impl
% buffer=floor(Nx/15);
% alpha = linspace(0,2*pi,buffer+1);
% uout=zeros(1,buffer+1);uin=uout;
% uin(end:-1:buffer/2+1)=(cos(alpha(alpha<=pi))/2 + 0.5);
% uout(1:buffer/2+1)=(cos(alpha(alpha<=pi))/2 + 0.5);
% Uinf=(Ly-(Y(:,Nx-buffer:Nx,:)-Ly/2).^2).*ones(size(u(:,Nx-buffer:Nx,:)));
%% Buffer logspace
buffer=floor(Nx/15);
alfa=logspace(0,1,buffer+1)/10;
nbytes = (length(buffer)+1)*4;fid = fopen('alfaalfa.dat','w');
fwrite(fid, nbytes, 'int32');fwrite(fid, alfa, 'single');fwrite(fid, nbytes, 'int32');fclose(fid);

Uinf=ones(size(u(:,Nx-buffer:Nx,:)));

for j=1:Nstep;
    %     %     while continuity>tol
    tic
    uh = fft(fft(u,[],1),[],2); vh = fft(fft(v,[],1),[],2); wh = fft(fft(w,[],1),[],2);
    uoldh = fft(fft(uold,[],1),[],2);    voldh = fft(fft(vold,[],1),[],2);    woldh = fft(fft(wold,[],1),[],2);
    %     %% Exact Convective Term Functions (Comment the ones to compare)
    %     %% Convective Terms Approximations
    %     %%%%%%%%%Convective Terms%%%%%%%%%
    %     % u*dudx & v*dudx & w*dudx (at time t & t-1) - non dealiased
    dudxh= bsxfun(@times, kx*1i*xi_x, uh);    dudxoldh= bsxfun(@times, kx*1i*xi_x, uoldh);
    dvdxh= bsxfun(@times, kx*1i*xi_x, vh);    dvdxoldh= bsxfun(@times, kx*1i*xi_x, voldh);
    dwdxh= bsxfun(@times, kx*1i*xi_x, wh);    dwdxoldh= bsxfun(@times, kx*1i*xi_x, woldh);
    ududx=ifft(ifft(dudxh,[],1),[],2).*u; ududxh = fft(fft(ududx,[],1),[],2); % Non Dealiased Products
    udvdx=ifft(ifft(dvdxh,[],1),[],2).*u; udvdxh = fft(fft(udvdx,[],1),[],2);
    udwdx=ifft(ifft(dwdxh,[],1),[],2).*u; udwdxh = fft(fft(udwdx,[],1),[],2);
    ududxold=ifft(ifft(dudxoldh,[],1),[],2).*uold; ududxoldh = fft(fft(ududxold,[],1),[],2); % Non Dealiased Products
    udvdxold=ifft(ifft(dvdxoldh,[],1),[],2).*uold; udvdxoldh  = fft(fft(udvdxold,[],1),[],2);
    udwdxold=ifft(ifft(dwdxoldh,[],1),[],2).*uold; udwdxoldh  = fft(fft(udwdxold,[],1),[],2);
    % 	 %% v*dudy & v*dvdy & v*dwdy (at time t & t-1) - dealiased
    dudyh=bsxfun(@times,ky'*1i*yi_y,uh);  dudyoldh=bsxfun(@times,ky'*1i*yi_y,uoldh);
    dvdyh=bsxfun(@times,ky'*1i*yi_y,vh);  dvdyoldh=bsxfun(@times,ky'*1i*yi_y,voldh);
    dwdyh=bsxfun(@times,ky'*1i*yi_y,wh);  dwdyoldh=bsxfun(@times,ky'*1i*yi_y,woldh);
    vdudy=ifft(ifft(dudyh,[],1),[],2).*v; vdudyh = fft(fft(vdudy,[],1),[],2); % Non Dealiased Products
    vdvdy=ifft(ifft(dvdyh,[],1),[],2).*v; vdvdyh = fft(fft(vdvdy,[],1),[],2);
    vdwdy=ifft(ifft(dwdyh,[],1),[],2).*v; vdwdyh = fft(fft(vdwdy,[],1),[],2);
    vdudyold=ifft(ifft(dudyoldh,[],1),[],2).*vold; vdudyoldh = fft(fft(vdudyold,[],1),[],2); % Non Dealiased Products
    vdvdyold=ifft(ifft(dvdyoldh,[],1),[],2).*vold; vdvdyoldh  = fft(fft(vdvdyold,[],1),[],2);
    vdwdyold=ifft(ifft(dwdyoldh,[],1),[],2).*vold; vdwdyoldh  = fft(fft(vdwdyold,[],1),[],2);
    f1=toc
    
    
    
    
    %     %% w*dudz & w*dvdz & w*dwdz (at time t & t-1) - dealiased
    %     tic
    for kkk = 1:Nx;for e=1:Ny;
            dudz(e,kkk,:) = Detagl *eta_zgl*squeeze(u(e,kkk,:)); dvdz(e,kkk,:) = Detagl *eta_zgl*squeeze(v(e,kkk,:));
            dwdz(e,kkk,:) = Detagl *eta_zgl*squeeze(w(e,kkk,:));
            dudzold(e,kkk,:) = Detagl *eta_zgl*squeeze(uold(e,kkk,:)); dvdzold(e,kkk,:) = Detagl *eta_zgl*squeeze(vold(e,kkk,:));
            dwdzold(e,kkk,:) = Detagl *eta_zgl*squeeze(wold(e,kkk,:));
        end;end
    %     toc
    wdudz=dudz.*w; wdudzh = fft(fft(wdudz,[],1),[],2);
    wdvdz=dvdz.*w; wdvdzh = fft(fft(wdvdz,[],1),[],2);
    wdwdz=dwdz.*w; wdwdzh = fft(fft(wdwdz,[],1),[],2);
    wdudzold=dudzold.*wold; wdudzoldh = fft(fft(wdudzold,[],1),[],2);
    wdvdzold=dvdzold.*wold; wdvdzoldh = fft(fft(wdvdzold,[],1),[],2);
    wdwdzold=dwdzold.*wold; wdwdzoldh = fft(fft(wdwdzold,[],1),[],2);
    %     % Loop over Fourier modes
    tic
    for m=1:length(ky);
        for oo=1:length(kx);
            % Solving for intermediate velocities; ustar & vsta
            A=Igldt-D2etagl*eta_zgl^2*vis+(xi_x^2*kx(oo)^2+yi_y^2*ky(m)^2)*Iglvis;
            
            convectiveu=2*ududxh(m,oo,:)+2*vdudyh(m,oo,:)+2*wdudzh(m,oo,:)-ududxoldh(m,oo,:)-vdudyoldh(m,oo,:)-wdudzoldh(m,oo,:);
            rhsuint=(2/dt)*uh(m,oo,:)-(1/(2*dt))*uoldh(m,oo,:)+convectiveu-th(m,oo,:);%-fh(m,oo,:);%
            rhsuint=squeeze(rhsuint);
            ustar       = A\rhsuint;
            ustar([1 Nz+1]) = 0;
            convectivev=2*vdvdyh(m,oo,:)+2*udvdxh(m,oo,:)+2*wdvdzh(m,oo,:)-udvdxoldh(m,oo,:)-vdvdyoldh(m,oo,:)-wdvdzoldh(m,oo,:);
            rhsvint=(2/dt)*vh(m,oo,:)-(1/(2*dt))*voldh(m,oo,:)+convectivev;
            rhsvint=squeeze(rhsvint);
            vstar       = A\rhsvint;
            vstar([1 Nz+1])=0;
            
            convectivew=2*vdwdyh(m,oo,:)+2*udwdxh(m,oo,:)+2*wdwdzh(m,oo,:)-udwdxoldh(m,oo,:)-vdwdyoldh(m,oo,:)-wdwdzoldh(m,oo,:);
            rhswint=(2/dt)*wh(m,oo,:)-(1/(2*dt))*woldh(m,oo,:)+convectivew;
            rhswint=squeeze(rhswint);
            wstar       = A\rhswint;
            wstar([1 Nz+1])=0;
            %
            %             % solve for pressure L*q=f;
            RHS = div_x*(kx(oo)*1i*xi_x)*ustar + div_x*(ky(m)*1i*yi_y)*vstar +  div_z*wstar;
            LHS = dt*( div_x_act_on_grad_x*(kx(oo)^2*xi_x^2+ky(m)^2*yi_y^2) + div_z_act_on_grad_z );
            %
            pnewh= LHS\RHS;
            %             % Update velocities
            unewh(m,oo,:) = ustar-dt*grad_x*(kx(oo)*1i*xi_x)*pnewh;
            vnewh(m,oo,:) = vstar-dt*grad_x*(ky(m)*1i*yi_y)*pnewh;
            wnewh(m,oo,:) = wstar-dt*grad_z*pnewh;
            %
            contcont(m,oo,:)=div_x*(kx(oo)*1i*xi_x)*squeeze(unewh(m,oo,:))+div_x*(ky(m)*1i*yi_y)*squeeze(vnewh(m,oo,:))+div_z*squeeze(wnewh(m,oo,:));
            %
        end
    end
    t_foruruer=toc
    unew = real(ifft(ifft(unewh,[],1),[],2));
    vnew = real(ifft(ifft(vnewh,[],1),[],2));
    wnew = real(ifft(ifft(wnewh,[],1),[],2));
    %     % Enforce BC
    unew(:,:,Nz+1)=1;unew(:,:,1)=1;
    vnew(:,:,Nz+1)=0;vnew(:,:,1)=0;
    wnew(:,:,Nz+1)=0;wnew(:,:,1)=0;
    
    %     steady=(1/(2*dt))*(3*unew-4*u+uold);s=max(abs(steady(:)));
    %     if s<tol; disp('solution stedy-tol 1e-6');break;end;
    %         save steadymapped;break;end
    %
    uold = real(ifft(ifft(uh,[],1),[],2));
    vold = real(ifft(ifft(vh,[],1),[],2));
    wold = real(ifft(ifft(wh,[],1),[],2));
    
%      for i=1:length(uout);unew(:,Nx-buffer+i-1,:)=unew(:,Nx-buffer+i-1,:)*uout(i)+Uinf(:,i,:)*uin(i);end
for i=1:buffer;unew(:,Nx-buffer+i-1,:)=unew(:,Nx-buffer+i-1,:)*(1-alfa(i))+Uinf(:,i,:)*alfa(i);end

    u=unew;v=vnew;w=wnew;
    rescont=sum(contcont(:));
    if rescont>1e-10 ;disp('Cont not satisfied');break;end
    %     if mod(j,10000)==0;disp(num2str(j/Nstep)); filename=['Nx',num2str(Nx),'Ny',num2str(Ny),'Nz',num2str(Nz),'Lx',num2str(Lx),'Ly',num2str(Ly),'Lz',num2str(Lz),'Re',num2str(Re),'dt',num2str(dt),'T',num2str(t(j)),'BDFAB2.mat'];	save(filename,'Nx','Ny','Nz','Re','j','Nstep','u','v','w','X','Y','Z','turb', 'dt','rescont','Lx','Ly','Lz','s');end
end
%
% twhole_codett=toc