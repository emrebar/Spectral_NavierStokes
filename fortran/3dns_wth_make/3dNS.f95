program halfmerged

USE DISPMODULE
!USE omp_lib

IMPLICIT NONE

include '/usr/local/include/fftw3.f'

INTEGER, PARAMETER :: Nx=156,Ny=32,Nz=28,Lx= 30,Ly=6,Nstep=1
REAL, PARAMETER:: pi=3.1415926535897932384626433832795028841971693993,Re=10,vis=1/Re,dt=1e-3,&
xi_x=2*pi/Lx,yi_y=2*pi/Ly,Lz=4,eta_zgl=2/Lz,eta_zg=2/Lz
INTEGER(KIND=4):: l,m,oo, pivot(Nz+1), pivott(Nz),ok,vandrow,vandcol,j,k,h,n,tt,planf,planb,count,iol
COMPLEX,PARAMETER::i = (0.0,1.0)
COMPLEX, DIMENSION(Nz+1,3):: BB
REAL ,DIMENSION(Nz+1):: zgl,etagl
REAL ,DIMENSION(Nz):: zg,etag
REAL ,DIMENSION(Nx):: x,kx
REAL, DIMENSION(Ny):: y,ky
REAL ,DIMENSION(Nz,Nz)::VetaG,dVetaG,ARAVG,DIAG1VG,DIAG2VG,Ig,invVG,Detag,&
VGforinv,div_x_act_on_grad_x,&
div_z_act_on_grad_z
REAL ,DIMENSION(Nz+1,Nz+1)::Vetagl,dVetagl,DIAG1,ARAVGLeta,DIAG2,Dgl,&
ZNz,ZNx,Igldt,Iglvis,invVGL,&
Igl,VGLforinv,D2etagl,Detagl
COMPLEX ,DIMENSION(Nz+1,Nz+1)::A
COMPLEX, DIMENSION(Nz+1):: convectiveu,rhsuint,convectivev,rhsvint,&
convectivew,rhswint,ustar,vstar,wstar
real ( kind = 8 ), DIMENSION(Ny,Nx,Nz+1):: temp
complex ( kind = 8 ), DIMENSION(Ny,Nx,Nz+1):: ududxh,vdudyh,wdudzh,ududxoldh,vdudyoldh,wdudzoldh,&
uh,uoldh,vdvdyh,udvdxh,wdvdzh,udvdxoldh,vdvdyoldh,wdvdzoldh,&
vdwdyh,udwdxh,wdwdzh,udwdxoldh,vdwdyoldh,wdwdzoldh,unewh,vnewh,wnewh,vh,voldh,wh,woldh,f,&
fh,u,v,w,uold,vold,wold,dudxh,dvdxh,dwdxh,dudxoldh,dvdxoldh,dwdxoldh,&
dudyh,dvdyh,dwdyh,dudyoldh,dvdyoldh,dwdyoldh,&
dudx,dudz,dvdx,dvdz,dwdx,dwdz,dudxold,dudzold,dvdxold,dvdzold,dwdxold,dwdzold,&
ududx,udvdx,udwdx,ududxold,udvdxold,udwdxold,&
wdvdz,wdudz,wdwdz,wdudzold,wdvdzold,wdwdzold,dudy,dudyold,dvdy,dvdyold,dwdy,dwdyold,&
vdudy,vdudyold,vdvdy,vdvdyold,vdwdy,vdwdyold,unew,vnew,wnew,contcont,savearray
COMPLEX, DIMENSION(Nz,Nz):: LHS
COMPLEX, DIMENSION(Nz):: RHS
COMPLEX, DIMENSION(Ny,Nx,Nz):: pnewh
REAL ,DIMENSION(Nz,Nz+1)::VG_trunc,div_x,div_z
REAL ,DIMENSION(Nz+1,Nz)::grad_x,grad_z
COMPLEX(kind=8), DIMENSION(:,:,:), ALLOCATABLE :: inf,outf
 CHARACTER*100 :: name_config,action,file
real :: start1, finish1,startloop,finishloop,start2, finish2,start3, finish3,start,finish,start4, finish4,start5, finish5
real :: startstart,finishfinish
INTEGER(kind=4)      :: rank,idist,odist,istride,ostride,howmany
integer,dimension(2)::ppp,inembed,onembed
INTEGER   :: numthreads,ierr

call cpu_time(start1)
!------------------------------------------------------------------------------------------------------------------------------------------!
!------------------------------------------------------------------------------------------------------------------------------------------!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! REQUIRED OPERATORS MODES,GRID!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!------------------------------------------------------------------------------------------------------------------------------------------!
!------------------------------------------------------------------------------------------------------------------------------------------!
! X,Y Coordinates
do j=1,Nx;x(j)=(((2*pi)/Nx)*(j-1))/xi_x;end do
do j=1,Ny;y(j)=(((2*pi)/Ny)*(j-1))/yi_y;end do
! Z COORDINATES- GAUSS & GL zg/zgl coord.
DO vandrow=1,Nz+1;etagl(vandrow)=-cos(pi*(vandrow-1)/Nz);END DO; zgl=(etagl+1)/eta_zgl;
DO vandrow=1,Nz;  etag(vandrow)=-cos(pi*(vandrow-0.5)/Nz);END DO; zg=(etag+1)/eta_zg;
! X,Y modes-k's
DO l=1,Nx/2;kx(l)=l-1;ENDDO;DO l=1,Nx/2;kx(l+Nx/2)=-Nx/2+l-1;ENDDO
DO l=1,Ny/2;ky(l)=l-1;ENDDO;DO l=1,Ny/2;ky(l+Ny/2)=-Ny/2+l-1;ENDDO
! VANDERMONDE MATRICES REQUIRED TO CONSTRUCT CHEBYSHEV DIFFERENTIATION MATRICES
!GAUSS LOBATTO dVGL
do vandrow=1,Nz+1;  do vandcol=1,Nz+1;
 Vetagl(vandrow,vandcol)=cos(acos(etagl(vandrow))*(vandcol-1));
 ARAVGLeta(vandrow,vandcol)=sin(acos(etagl(vandrow))*(vandcol-1));
end do ;end do
 DIAG1=0.0;DIAG2=0.0;
do vandcol=1,Nz+1;DIAG2(vandcol,vandcol)=vandcol-1;  end do 
do vandcol=2,Nz;DIAG1(vandcol,vandcol)=1/SQRT(1-etagl(vandcol)**2);  end do
 dVetagl=MATMUL(DIAG1,ARAVGLeta)
 dVetagl=MATMUL(dVetagl,DIAG2)
do j=1,Nz+1;  dVetagl(1,j)=((-1)**j)*((j-1)**2);  dVetagl(Nz+1,j)=(j-1)**2;end do
! VANDERMONDE GAUSS dVG
do vandrow=1,Nz;  do vandcol=1,Nz;  
Vetag(vandrow,vandcol)=cos(acos(etag(vandrow))*(vandcol-1));ARAVG(vandrow,vandcol)=sin(acos(etag(vandrow))*(vandcol-1));
end do ;end do;
 DIAG1VG=0.0;DIAG2VG=0.0;
do j=1,Nz;DIAG1VG(j,j)=1/SQRT(1-etag(j)**2);end do
do j=1,Nz;DIAG2VG(j,j)=j-1;  end do 
 dVetaG=MATMUL(DIAG1VG,ARAVG);dVetaG=MATMUL(dVetaG,DIAG2VG)
!DIFFERENTIATION MATRIX  GAUSS-LOBATTO
 VGLforinv=Vetagl;CALL inverse(VGLforinv,invVGL,Nz+1);Detagl=MATMUL(dVetagl,invVGL)
D2etagl=MATMUL(Detagl,Detagl);
!DIFFERENTIATION MATRIX GAUSS 
 VGforinv=vetag; CALL inverse(VGforinv,invVG,Nz);Detag=MATMUL(dVetaG,invVG)

! OPERATORs
 ZNx=0.0;do l = 1,Nz+1;ZNx(l, l) = 1; end do;
 ZNz=0.0;do l = 2,Nz;ZNz(l, l) = 1; end do;
 Igl=0.0;do l = 1,Nz+1;Igl(l, l) = 1; end do; 
 Ig=0.0;do l = 1,Nz;Ig(l, l) = 1; end do; 
 Igldt=0.0;do l = 1,Nz+1;Igldt(l, l) = (3/(2*dt)); end do; 
 Iglvis=0.0;do l = 1,Nz+1;Iglvis(l, l) = vis; end do;

 VG_trunc=0;do l=1,Nz; VG_trunc(:,l)=Vetag(:,l); end do;
 div_x=MATMUL(VG_trunc,invVGL);div_x=MATMUL(div_x,Igl);div_x=MATMUL(div_x,ZNx)
 div_z=MATMUL(VG_trunc,invVGL);div_z=MATMUL(div_z,Detagl*eta_zgl);div_z=MATMUL(div_z,ZNz)
 grad_x=MATMUL(Vetagl(:,1:Nz),invVG);grad_x=MATMUL(grad_x,Ig);
 grad_z=MATMUL(Vetagl(:,1:Nz),invVG);grad_z=MATMUL(grad_z,Detag);grad_z=grad_z*eta_zg
 div_x_act_on_grad_x=-Ig;
 div_z_act_on_grad_z=MATMUL(div_z,grad_z)
call cpu_time(finish1)
!------------------------------------------------------------------------------------------------------------------------------------------!
!------------------------------------------------------------------------------------------------------------------------------------------!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! GENERAL FFT PLAN BACKWARD & FORWARD!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!------------------------------------------------------------------------------------------------------------------------------------------!
!------------------------------------------------------------------------------------------------------------------------------------------!
ALLOCATE(inf(1:Ny,1:Nx,1:Nz+1),outf(1:Ny,1:Nx,1:Nz+1))
  rank = 2;
  ppp(1) =Ny;  ppp(2) =Nx;
  howmany =Nz+1
  odist = ppp(1)*ppp(2);    idist =  ppp(1)*ppp(2); 
  ostride = 1;  istride = 1;
  inembed(:) = ppp(:); onembed(:) = ppp(:);

call dfftw_plan_many_dft(planf,rank,ppp,howmany,inf,inembed,istride,idist,outf, onembed,ostride,odist, FFTW_FORWARD, FFTW_MEASURE);
call dfftw_plan_many_dft(planb,rank,ppp,howmany,outf,inembed,istride,idist,inf, onembed,ostride,odist, FFTW_BACKWARD, FFTW_MEASURE);


!------------------------------------------------------------------------------------------------------------------------------------------!
!------------------------------------------------------------------------------------------------------------------------------------------!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Constant Pressure Gradient Real&Fourier Space  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!------------------------------------------------------------------------------------------------------------------------------------------!
!------------------------------------------------------------------------------------------------------------------------------------------!
!f=-2/Re;
!call dfftw_execute_dft(planf,f,fh);
!------------------------------------------------------------------------------------------------------------------------------------------!
!------------------------------------------------------------------------------------------------------------------------------------------!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! INITIAL  CONDITIONS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!------------------------------------------------------------------------------------------------------------------------------------------!
!------------------------------------------------------------------------------------------------------------------------------------------!
!call random_number(temp)
u=1;uold=1;v=0;vold=0;w=0;wold=0;
!u=temp;uold=temp;v=temp;vold=temp;w=temp;wold=temp;



do 100 tt=1,Nstep
!------------------------------------------------------------------------------------------------------------------------------------------!
!------------------------------------------------------------------------------------------------------------------------------------------!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! VELOCITY COMPONENTS  TRANSFROM TO FOURIER SPACE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!------------------------------------------------------------------------------------------------------------------------------------------!
!------------------------------------------------------------------------------------------------------------------------------------------!
call cpu_time(start2)
! U & UOLD
call dfftw_execute_dft_(planf,u,uh);call dfftw_execute_dft_(planf,uold,uoldh);
! V & Vold
call dfftw_execute_dft_(planf,v,vh);call dfftw_execute_dft_(planf,vold,voldh);
! w & Wold
call dfftw_execute_dft_(planf,w,wh);call dfftw_execute_dft_(planf,wold,woldh);
!------------------------------------------------------------------------------------------------------------------------------------------!
!------------------------------------------------------------------------------------------------------------------------------------------!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!PSEUDO SPECTRAL CONVECTIVE TERMS -NOT DEALIASED!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!------------------------------------------------------------------------------------------------------------------------------------------!
!------------------------------------------------------------------------------------------------------------------------------------------!
! X DERIVATIVES IN FOURIER SPACE du/dx&dv/dx&dw/dx at time t & t-1

do l=1,Nx;
dudxh(:,l,:)= uh(:,l,:)*kx(l)*i*xi_x;      dvdxh(:,l,:)= vh(:,l,:)*kx(l)*i*xi_x;      dwdxh(:,l,:)= wh(:,l,:)*kx(l)*i*xi_x;
dudxoldh(:,l,:)= uoldh(:,l,:)*kx(l)*i*xi_x;dvdxoldh(:,l,:)= voldh(:,l,:)*kx(l)*i*xi_x;dwdxoldh(:,l,:)= woldh(:,l,:)*kx(l)*i*xi_x;
end do

! Y DERIVATIVES IN FOURIER SPACE du/dz&dv/dz&dw/dz at time t & t-1
call cpu_time(startstart)


do l=1,Ny; 
dudyh(l,:,:)= uh(l,:,:)*ky(l)*i*yi_y;      dvdyh(l,:,:)= vh(l,:,:)*ky(l)*i*yi_y;      dwdyh(l,:,:)= wh(l,:,:)*ky(l)*i*yi_y;
dudyoldh(l,:,:)= uoldh(l,:,:)*ky(l)*i*yi_y;dvdyoldh(l,:,:)= voldh(l,:,:)*ky(l)*i*yi_y;dwdyoldh(l,:,:)= woldh(l,:,:)*ky(l)*i*yi_y;
end do



call cpu_time(finishfinish)
! PSEUDO SPECTRAL BACKWARD TRANSFORM
call dfftw_execute_dft_(planb,dudxh,dudx);call dfftw_execute_dft_(planb,dudxoldh,dudxold);
dudx=dudx/(Nx*Ny);  dudxold=dudxold/(Nx*Ny);
call dfftw_execute_dft_(planb,dvdxh,dvdx);call dfftw_execute_dft_(planb,dvdxoldh,dvdxold);
dvdx=dvdx/(Nx*Ny);  dvdxold=dvdxold/(Nx*Ny);
call dfftw_execute_dft_(planb,dwdxh,dwdx);call dfftw_execute_dft_(planb,dwdxoldh,dwdxold);
dwdx=dwdx/(Nx*Ny);  dwdxold=dwdxold/(Nx*Ny);
call dfftw_execute_dft_(planb,dudyh,dudy);call dfftw_execute_dft_(planb,dudyoldh,dudyold);
dudy=dudy/(Nx*Ny);  dudyold=dudyold/(Nx*Ny);
call dfftw_execute_dft_(planb,dvdyh,dvdy);call dfftw_execute_dft_(planb,dvdyoldh,dvdyold);
dvdy=dvdy/(Nx*Ny); dvdyold=dvdyold/(Nx*Ny);  
call dfftw_execute_dft_(planb,dwdyh,dwdy);call dfftw_execute_dft_(planb,dwdyoldh,dwdyold);
dwdy=dwdy/(Nx*Ny);  dwdyold=dwdyold/(Nx*Ny);
!!!!!NON LINEAR TERM MULTIPLICATION IN REAL SPACE
ududx=u*dudx;udvdx=u*dvdx;udwdx=u*dwdx;ududxold=uold*dudxold;udvdxold=uold*dvdxold;udwdxold=uold*dwdxold;
vdudy=v*dudy;vdvdy=v*dvdy;vdwdy=v*dwdy;vdudyold=vold*dudyold;vdvdyold=vold*dvdyold;vdwdyold=vold*dwdyold;

!! FORWARD TRANSFORM OF NON LINEAR TERMS
call dfftw_execute_dft_(planf,ududx,ududxh);call dfftw_execute_dft_(planf,ududxold,ududxoldh);
call dfftw_execute_dft_(planf,udvdx,udvdxh);call dfftw_execute_dft_(planf,udvdxold,udvdxoldh);
call dfftw_execute_dft_(planf,udwdx,udwdxh);call dfftw_execute_dft_(planf,udwdxold,udwdxoldh);

call dfftw_execute_dft_(planf,vdudy,vdudyh);call dfftw_execute_dft_(planf,vdudyold,vdudyoldh);
call dfftw_execute_dft_(planf,vdvdy,vdvdyh);call dfftw_execute_dft_(planf,vdvdyold,vdvdyoldh);
call dfftw_execute_dft_(planf,vdwdy,vdwdyh);call dfftw_execute_dft_(planf,vdwdyold,vdwdyoldh);
call cpu_time(finish2)



!------------------------------------------------------------------------------------------------------------------------------------------!
!------------------------------------------------------------------------------------------------------------------------------------------!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! NON PERIODIC BOUNDARIES IN Z DIRECTION  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!------------------------------------------------------------------------------------------------------------------------------------------!
!------------------------------------------------------------------------------------------------------------------------------------------!
!Z DERIVATIVES IN REAL SPACE FOR NON LINEAR TERMS 
call cpu_time(start3)

do l=1,Nx;
do m=1,Ny;
dudz(m,l,:)=MATMUL(Detagl*eta_zgl,u(m,l,:));dvdz(m,l,:)=MATMUL(Detagl*eta_zgl,v(m,l,:)); dwdz(m,l,:)=MATMUL(Detagl*eta_zgl,w(m,l,:));
dudzold(m,l,:)=MATMUL(Detagl*eta_zgl,uold(m,l,:));
dvdzold(m,l,:)=MATMUL(Detagl*eta_zgl,vold(m,l,:)); dwdzold(m,l,:)=MATMUL(Detagl*eta_zgl,wold(m,l,:));
end do;
end do;


call cpu_time(finish3)
!  MULTIPLICATION IN REAL SPACE FOR NON LINEAR TERMS
wdudz=w*dudz;wdvdz=w*dvdz;wdwdz=w*dwdz;wdudzold=wold*dudzold;wdvdzold=wold*dvdzold;wdwdzold=wold*dwdzold;

! FORWARD FFT FOR NONLINEAR Z DERIVATIVE TERMS
call cpu_time(start4)
call dfftw_execute_dft_(planf,wdudz,wdudzh);call dfftw_execute_dft_(planf,wdudzold,wdudzoldh);
call dfftw_execute_dft_(planf,wdvdz,wdvdzh);call dfftw_execute_dft_(planf,wdvdzold,wdvdzoldh);
call dfftw_execute_dft_(planf,wdwdz,wdwdzh);call dfftw_execute_dft_(planf,wdwdzold,wdwdzoldh);
call cpu_time(finish4)
!------------------------------------------------------------------------------------------------------------------------------------------!
!------------------------------------------------------------------------------------------------------------------------------------------!
!!!!!!!!!!!!!!!!!!!!!! SOLVING THE LINEAR SYSTEMS FOR INT. VELOCITIES AND PRESSURE  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!---------------------------loop nr.23 for Y modes & loop nr.36 for X modes  --------------------------------------------------------------!
!------------------------------------------------------------------------------------------------------------------------------------------!
!------------------------------------------------------------------------------------------------------------------------------------------!

!call cpu_time(start5)
!!!!!!!!!!doprivate(A,convectiveu,rhsuint,convectivev,rhsvint,convectivew,rhswint,BB,ustar,vstar,wstar,RHS,LHS)
DO 23 m=1,Nx           !size(kz)
        DO 36 oo=1,Ny  !size(kx);
A=Igldt-D2etagl*(eta_zgl**2)*vis+((xi_x**2)*kx(m)**2)*Iglvis+((yi_y**2)*ky(m)**2)*Iglvis;
! Convective U
convectiveu=2*ududxh(oo,m,:)+2*vdudyh(oo,m,:)+2*wdudzh(oo,m,:)-ududxoldh(oo,m,:)-vdudyoldh(oo,m,:)-wdudzoldh(oo,m,:);
rhsuint=(2/dt)*uh(oo,m,:)-(1/(2*dt))*uoldh(oo,m,:)+convectiveu;!-fh(oo,m,:);
! Convective V
convectivev=2*vdvdyh(oo,m,:)+2*udvdxh(oo,m,:)+2*wdvdzh(oo,m,:)-udvdxoldh(oo,m,:)-vdvdyoldh(oo,m,:)-wdvdzoldh(oo,m,:);
rhsvint=(2/dt)*vh(oo,m,:)-(1/(2*dt))*voldh(oo,m,:)+convectivev;
! Convective W
convectivew=2*vdwdyh(oo,m,:)+2*udwdxh(oo,m,:)+2*wdwdzh(oo,m,:)-udwdxoldh(oo,m,:)-vdwdyoldh(oo,m,:)-wdwdzoldh(oo,m,:);
rhswint=(2/dt)*wh(oo,m,:)-(1/(2*dt))*woldh(oo,m,:)+convectivew;
! Assign B for simultaneous G.E.
BB(:,1)=rhsuint;BB(:,2)=rhsvint;BB(:,3)=rhswint;

! DO the G.E.
call CGESV(Nz+1, 3, A, Nz+1, pivot, BB, Nz+1, ok)
! Assign to inter. velocities
ustar=BB(:,1);vstar=BB(:,2);wstar=BB(:,3);
! Impose the boundary conditions.
ustar(1)=0;ustar(Nz+1)=0;vstar(1)=0;vstar(Nz+1)=0;wstar(1)=0;wstar(Nz+1)=0;
RHS = MATMUL(div_x,((kx(m)*i*xi_x)*ustar)) +  MATMUL(div_x,((ky(m)*i*yi_y)*vstar))+MATMUL(div_z,wstar);
LHS = (div_x_act_on_grad_x*(kx(m)**2*xi_x**2+ky(m)**2*yi_y**2) + div_z_act_on_grad_z );
LHS=dt*LHS;

call CGESV(Nz, 1, LHS, Nz, pivott, RHS, Nz, ok)
pnewh(oo,m,:)=RHS;
! Solve for velocities at new time step
unewh(oo,m,:) = ustar-MATMUL(dt*grad_x*(kx(m)*i*xi_x),pnewh(oo,m,:));
vnewh(oo,m,:) = vstar-MATMUL(dt*grad_x*(ky(m)*i*yi_y),pnewh(oo,m,:));
wnewh(oo,m,:) = wstar-MATMUL(dt*grad_z,pnewh(oo,m,:));



36 end do
23 end do

!call cpu_time(finish5)



!------------------------------------------------------------------------------------------------------------------------------------------!
!------------------------------------------------------------------------------------------------------------------------------------------!
!!!!!!!!!!!!!!!!!!!!!!  UPDATE THE VELOCITIES & IMPOSE BOUNDARY CONDITIONS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!------------------------------------------------------------------------------------------------------------------------------------------!
!------------------------------------------------------------------------------------------------------------------------------------------!
! assign u to uold 
call dfftw_execute_dft_(planb,unewh,u);u=REAL(u/(Nx*Ny));
! assign v to wold
call dfftw_execute_dft_(planb,vnewh,v);v=REAL(v/(Nx*Ny));
! assign w to wold
call dfftw_execute_dft_(planb,wnewh,w);w=REAL(w/(Nx*Ny));

!impose boundaries.
u(:,:,Nz+1)=1;u(:,:,1)=1;
v(:,:,Nz+1)=0;v(:,:,1)=0;
w(:,:,Nz+1)=0;w(:,:,1)=0;

! assign u to uold 
call dfftw_execute_dft_(planb,uh,uold);uold=REAL(uold/(Nx*Ny));
! assign v to wold
call dfftw_execute_dft_(planb,vh,vold);vold=REAL(vold/(Nx*Ny));
! assign w to wold
call dfftw_execute_dft_(planb,wh,wold);wold=REAL(wold/(Nx*Ny));


100 end do

call cpu_time(finish)

write(*,*) finishfinish-startstart
!write(*,*) finish1-start1
!write(*,*) finish2-start2
!write(*,*) finish3-start3
!write(*,*) finish4-start4
!write(*,*) finish5-start5

!------------------------------------------------------------------------------------------------------------------------------------------!
!------------------------------------------------------------------------------------------------------------------------------------------!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  SAVE DATA !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!------------------------------------------------------------------------------------------------------------------------------------------!
!------------------------------------------------------------------------------------------------------------------------------------------!

PRINT *, 'Nx = ', Nx, ', Ny = ', Ny,',Nz=',Nz,',timestep=',tt

name_config = './data/info.dat'
OPEN ( unit =11 , FILE = name_config  );REWIND (11);WRITE (11 ,*) Nx,Ny,Nz,Lx,Ly,Lz,tt,Re,dt;CLOSE (11)

name_config = './data/u';savearray = REAL(u);CALL savedata(Ny,Nx,Nz+1,1,name_config,savearray)
name_config = './data/v';savearray = REAL(v);CALL savedata(Ny,Nx,Nz+1,1,name_config,savearray)
name_config = './data/w';savearray = REAL(w);CALL savedata(Ny,Nx,Nz+1,1,name_config,savearray)

!do l=1,Nz+1
!write(*,*) 'u',l
!CALL DISP((u(:,:,l)))
!end do

end program halfmerged





