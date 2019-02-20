%% Clean
clc;clear all;close all;warning('off');
%  matlabpool
%% global variables for force distribution functions
global CENTERS x ygl z X Y Z Nx Ny j Nz turbindex omega dx ep f_eps  vol
global alfclcd12 alfclcd15 alfclcd18 alfclcd21 alfclcd24 bladerad chordtjae elementarea tempbladerad thicktja twisttja toverc
global aind apind radinside SampleTurbineR profcatalog yncent zncent forcewithinarea f_eps1 f_eps2
%% Distribute FORCES
load feps1;load feps2;
%% Load necessry
load preprocess.mat;load turbines.mat;load tjaturbine
%% Turbulence Inflow Initialize
% load turbinflowacd
% uold=u; vold=v; wold=w;
[wold,vold,uold,u,v,w,unewh,vnewh,wnewh,dudy,dvdy,dwdy,dudyold,dvdyold,dwdyold]=deal(zeros(size(X)));



%% Operation condition of turbnine
omega=0.3;%lam=7
%%
thx = zeros(size(u));thy = thx;thz = thx;
%% Starting
aind=0.3;apind=(1-3*aind)/(4*aind-1);%tangential induction
%% Buffer Zone log-increment
buffer=floor(Nx/15);
alfa=logspace(0,1,buffer+1)/10;
Uinf=-ones(size(u(:,Nx-buffer:Nx,:)));

%% Save for information - keep the track of simulation
filen=['Nx',num2str(Nx),'Ny',num2str(Ny),'Nz',num2str(Nz),'Lx',num2str(Lx),'Ly',num2str(Ly),'Lz',num2str(Lz),'Re',num2str(Re),'dt',num2str(dt),'.mat'];
save(filen,'Nx','Ny','Nz','Re','X','Y','Z', 'dt','Lx','Ly','Lz');
disp('Saved nicely')
%% Transform initial cond.
uoldh = fft(fft(uold,[],3),[],2);    voldh = fft(fft(vold,[],3),[],2);    woldh = fft(fft(wold,[],3),[],2);
%% old X derivatives initialize
dudxoldh= bsxfun(@times, kx*1i*xi_x, uoldh); dwdxoldh= bsxfun(@times, kx*1i*xi_x, woldh);    dvdxoldh= bsxfun(@times, kx*1i*xi_x, voldh);
ududxold=ifft(ifft(dudxoldh,[],3),[],2).*uold; ududxoldh = fft(fft(ududxold,[],3),[],2); % Non Dealiased Products
udvdxold=ifft(ifft(dvdxoldh,[],3),[],2).*uold; udvdxoldh  = fft(fft(udvdxold,[],3),[],2);
udwdxold=ifft(ifft(dwdxoldh,[],3),[],2).*uold; udwdxoldh  = fft(fft(udwdxold,[],3),[],2);
%% old Z derivatives initialize
dudzoldh=bsxfun(@times,kzz*1i*zi_z,uoldh);   dvdzoldh=bsxfun(@times,kzz*1i*zi_z,voldh);  dwdzoldh=bsxfun(@times,kzz*1i*zi_z,woldh);
wdudzold=ifft(ifft(dudzoldh,[],3),[],2).*wold;wdudzoldh = fft(fft(wdudzold,[],3),[],2);
wdvdzold=ifft(ifft(dvdzoldh,[],3),[],2).*wold; wdvdzoldh = fft(fft(wdvdzold,[],3),[],2);
wdwdzold=ifft(ifft(dwdzoldh,[],3),[],2).*wold; wdwdzoldh = fft(fft(wdwdzold,[],3),[],2);
%% old Y derivatives initialize
for e = 1:Nz;
    dudyold(:,:,e) = Detagleta_ygl*uold(:,:,e); dvdyold(:,:,e) = Detagleta_ygl*vold(:,:,e); dwdyold(:,:,e) = Detagleta_ygl*wold(:,:,e);
end
vdudyold=dudyold.*vold; vdudyoldh = fft(fft(vdudyold,[],3),[],2);
vdvdyold=dvdyold.*vold; vdvdyoldh = fft(fft(vdvdyold,[],3),[],2);
vdwdyold=dwdyold.*vold; vdwdyoldh = fft(fft(vdwdyold,[],3),[],2);
disp('Initials Transformed')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%START TIME LOOP %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Time loop Started')
disp(['Re',num2str(Re),'dt',num2str(dt),'Tfinal',num2str(tfinal) ])
for j=1:1;
    %%   TRANSFORM THE INITIAL CONDITIONS TO FOURIER SPACE; %   WILL BE UPDATED   EVERY STEP
    uh = fft(fft(u,[],3),[],2); vh = fft(fft(v,[],3),[],2); wh = fft(fft(w,[],3),[],2);
    
    %%  PSEUDEO SPECTRAL CONVECTIVE TERMS at time T
    % u*dudx & v*dudx & w*dudx (at time t) - not dealiased
    dudxh= bsxfun(@times, kx*1i*xi_x, uh);
    dvdxh= bsxfun(@times, kx*1i*xi_x, vh);
    dwdxh= bsxfun(@times, kx*1i*xi_x, wh);
    
    ududx=ifft(ifft(dudxh,[],3),[],2).*u; ududxh = fft(fft(ududx,[],3),[],2); % Non Dealiased Products
    udvdx=ifft(ifft(dvdxh,[],3),[],2).*u; udvdxh = fft(fft(udvdx,[],3),[],2);
    udwdx=ifft(ifft(dwdxh,[],3),[],2).*u; udwdxh = fft(fft(udwdx,[],3),[],2);
    
    %   w*dudz & w*dvdz & w*dwdz (at time t) - dealiased
    dudzh=bsxfun(@times,kzz*1i*zi_z,uh);
    dvdzh=bsxfun(@times,kzz*1i*zi_z,vh);
    dwdzh=bsxfun(@times,kzz*1i*zi_z,wh);
    
    wdudz=ifft(ifft(dudzh,[],3),[],2).*w; wdudzh = fft(fft(wdudz,[],3),[],2);
    wdvdz=ifft(ifft(dvdzh,[],3),[],2).*w; wdvdzh = fft(fft(wdvdz,[],3),[],2);
    wdwdz=ifft(ifft(dwdzh,[],3),[],2).*w; wdwdzh = fft(fft(wdwdz,[],3),[],2);
    
    
    %   v*dudy & v*dvdy & v*dwdy (at time t ) -Chebyshev
    for e = 1:Nz;
        dudy(:,:,e) = Detagleta_ygl*u(:,:,e); dvdy(:,:,e) = Detagleta_ygl*v(:,:,e); dwdy(:,:,e) = Detagleta_ygl*w(:,:,e);
    end
    vdudy=dudy.*v; vdudyh = fft(fft(vdudy,[],3),[],2);
    vdvdy=dvdy.*v; vdvdyh = fft(fft(vdvdy,[],3),[],2);
    vdwdy=dwdy.*v; vdwdyh = fft(fft(vdwdy,[],3),[],2);
    
    % Loop over Fourier modes
%     parfor m=1:Nz;
    for m=1:Nz;
        for oo=1:Nx;
            warning('OFF');
            
            % Solving for intermediate velocities; ustar & vstar & wstar
            A=operate+(xi_x^2*kx(oo)^2+kz(m)^2*zi_z^2)*Iglvis;
            
            convectiveu=2*ududxh(:,oo,m)+2*vdudyh(:,oo,m)+2*wdudzh(:,oo,m)-ududxoldh(:,oo,m)-vdudyoldh(:,oo,m)-wdudzoldh(:,oo,m);
            rhsuint=(2/dt)*uh(:,oo,m)-(1/(2*dt))*uoldh(:,oo,m)+convectiveu+thx(:,oo,m);%-fh(:,oo,m);%;
                       
            convectivev=2*vdvdyh(:,oo,m)+2*udvdxh(:,oo,m)+2*wdvdzh(:,oo,m)-udvdxoldh(:,oo,m)-vdvdyoldh(:,oo,m)-wdvdzoldh(:,oo,m);
            rhsvint=(2/dt)*vh(:,oo,m)-(1/(2*dt))*voldh(:,oo,m)+convectivev+thy(:,oo,m);
           
            convectivew=2*vdwdyh(:,oo,m)+2*udwdxh(:,oo,m)+2*wdwdzh(:,oo,m)-udwdxoldh(:,oo,m)-vdwdyoldh(:,oo,m)-wdwdzoldh(:,oo,m);
            rhswint=(2/dt)*wh(:,oo,m)-(1/(2*dt))*woldh(:,oo,m)+convectivew+thz(:,oo,m);
            
            total=A\[rhsuint rhsvint rhswint];total([1 Ny+1],:)=0;
            ustar=total(:,1);vstar=total(:,2);wstar=total(:,3);
            
            
            % solve for pressure L*q=f;
            RHS = div_x*(kx(oo)*1i*xi_x)*ustar + div_y*vstar +  div_x*(kz(m)*1i*zi_z)*wstar;
            LHS = dt*( div_x_act_on_grad_x*(kx(oo)^2*xi_x^2+kz(m)^2*zi_z^2) + div_y_act_on_grad_y );
            
            pnewh = LHS\RHS;
            % Update velocities
            unewh(:,oo,m) = ustar-dt*grad_x*(kx(oo)*1i*xi_x)*pnewh;
            vnewh(:,oo,m) = vstar-dt*grad_y*pnewh;
            wnewh(:,oo,m) = wstar-dt*grad_x*(kz(m)*1i*zi_z)*pnewh;
            
            % Check the continuity
%             if mod(j,150)==0; contcont(:,oo,m)=div_x*(kx(oo)*1i*xi_x)*unewh(:,oo,m)+div_y*vnewh(:,oo,m)+div_x*(kz(m)*1i*zi_z)*wnewh(:,oo,m);
%             end
            
        end
    end
    
    unew = real(ifft(ifft(unewh,[],3),[],2));    vnew = real(ifft(ifft(vnewh,[],3),[],2));    wnew = real(ifft(ifft(wnewh,[],3),[],2));
    % Enforce BC
    unew([1 Ny+1],:) = -1;        vnew([1 Ny+1],:) = 0;    wnew([1 Ny+1],:) = 0;  % B.C.
    % Update 'uold'
    uoldh =uh;  voldh = vh;    woldh = wh;
    ududxoldh= ududxh;    udvdxoldh= udvdxh;    udwdxoldh= udwdxh;
    wdudzoldh= wdudzh;    wdvdzoldh= wdvdzh;    wdwdzoldh= wdwdzh;
    vdudyoldh=vdudyh;     vdvdyoldh=vdvdyh;     vdwdyoldh=vdwdyh;

    
    % Impose buffer zone
    for i=1:buffer+1;        unew(:,Nx-buffer+i-1,:)=unew(:,Nx-buffer+i-1,:)*(1-alfa(i))+Uinf(:,i,:)*alfa(i);    end
    
    % Update 'u'
    u=unew;v=vnew;w=wnew;
    % Update the forces for ACD-R
    if j<300
        [turbxx,turbyy,turbzz]=fxftheta(u,v,w);
    else
        if mod(j,20)==0;
            [turbxx,turbyy,turbzz]=fxftheta(u,v,w);
        end
    end
    % Transform the calculated forces to fourier space
    thy = fft(fft(turbyy,[],3),[],2);    thz = fft(fft(turbzz,[],3),[],2);    thx = fft(fft(turbxx,[],3),[],2);
    
    if mod(j,50)==0;
        uvelside=u(:,:,Nz/2);filen=['U_sideACDR500/umiddomain',num2str(j*dt*10)];
        save(filen,'uvelside');
        vvelside=v(:,:,Nz/2);filen=['V_sideACDR500/vmiddomain',num2str(j*dt*10)];
        save(filen,'vvelside');
        wvelside=w(:,:,Nz/2);filen=['W_sideACDR500/wmiddomain',num2str(j*dt*10)];
        save(filen,'wvelside');
    end
    
    
    % Save every ? steps
    if mod(j,125)==0;
%         rescont=sum(contcont(:));
%         if rescont>1e-10 ;disp('Cont not satisfied');save('asd');break;end
        sumy=sum(abs(turbyy(:)));        sumz=sum(abs(turbzz(:)));        sumx=sum(abs(turbxx(:)));
        disp(num2str(j/Nstep)); filename=['Nx',num2str(Nx),'Ny',num2str(Ny),'Nz',num2str(Nz),'Lx',num2str(Lx),'Ly',num2str(Ly),'Lz',num2str(Lz),'Re',num2str(Re),'dt',num2str(dt),'T',num2str(t(j)),'BDFAB2buffered.mat'];
        save(filename,'j','u','v','w','sumx','sumy','sumz');
    end
end

