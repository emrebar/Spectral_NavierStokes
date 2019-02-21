program test
use DISPMODULE
INTEGER, PARAMETER :: Nx=128,Ny=32,Nz=28
integer :: l,buffer,k,oo
REAL, DIMENSION(:), ALLOCATABLE :: alfa
Real(kind=8), DIMENSION(:,:,:), ALLOCATABLE ::Uinf
  real, dimension(Ny,Nx,Nz+1) :: turb

buffer=floor((REAL(Nx)/15));
ALLOCATE (alfa(1:buffer+1))
!open(12, file = 'alfaalfa.dat', form = 'unformatted');  read(12) alfa;  close(12)
  open(unit=1, file='alfaalfa.dat', status='OLD', action='READ');
read(1) alfa;  close(1);

ALLOCATE (Uinf(1:Ny,(Nx-buffer):Nx,1:Nz+1))
do l=1,Nz+1; do k=Nx-buffer,Nx;do oo=1,Ny;Uinf(oo,k,l)=1;end do;end do;end do;


  open(11, file = 'turbine.dat', form = 'unformatted');  read(11) turb;  close(11)

!do l=1,Nz+1;write(*,*) 'u'; call disp(turb(:,:,l));end do
call DISP(alfa)
CALL DISP(sum(turb))
end program test


