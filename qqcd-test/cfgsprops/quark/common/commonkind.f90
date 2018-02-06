!
! common.f90
!

!----------------changes made to the module---------------
!July 5th, 2007 changes made by Abdou
!In order to include the 4 components of the axial current
!in the loop calculations we change nsav from 12 to 20 .
!nsav=2*nop since it counts both real and imaginary parts.
!----------------------------------------------------------


  module commonkind
  use kinds
!   integer, parameter                           :: KI=8,KR=8,KR2=8,KKC=8
    integer(kind=KI), parameter                  :: nd=4,nc=3,nri=2
!   integer(kind=KI), parameter                  :: nx=8,ny=8,nz=8,nt=8
!   integer(kind=KI), parameter                  :: nxyzt=nx*ny*nz*nt
!   integer(kind=KI), parameter                  :: nxyz=nx*ny*nz
    !integer(kind=KI), parameter                  :: nsav=12
    integer(kind=KI), parameter                  :: nsav=20

!   integer(kind=KI), public, parameter          :: npx=1,npy=1,npz=1,npt=1
!   integer(kind=KI), public, parameter          :: npx=2,npy=2,npz=2,npt=2
!   integer(kind=KI), public, parameter          :: nps=npx*npy*npz*npt
!   integer(kind=KI), public, parameter          :: nvhalf=nxyzt/32/nps
!   integer(kind=KI), public, parameter          :: nbmax=max(2*nvhalf*npx/nx, &
!                                                             2*nvhalf*npy/ny, &
!                                                             2*nvhalf*npz/nz, &
!                                                             2*nvhalf*npt/nt)
!   integer(kind=KI), public, parameter          :: ntotal=nvhalf+nbmax

! For MPI_REAL select mpitype==1, else for MPI_DOUBLE_PRECISION select mpitype==2

!!! ...REMEMBER TO CHANGE KI,KR,KR2,KKC... !!!

!   integer(kind=KI), public, parameter          :: mpitype=2



    
  end module commonkind
