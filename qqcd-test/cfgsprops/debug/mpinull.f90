! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
! mpinull.f90, Randy Lewis, randy.lewis@uregina.ca
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
! This module mimics the true MPI module for f90.
! It is intended for testing compilation, but not for execution.
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 module MPI

    use kinds
    implicit none
    private

    integer, public, parameter :: MPI_STATUS_SIZE=4
    integer, public, parameter :: MPI_REAL=26
    integer, public, parameter :: MPI_DOUBLE_PRECISION=27
    integer, public, parameter :: MPI_COMM_WORLD=91
    integer, public, parameter :: MPI_SUM=102

    public :: MPI_WTIME, MPI_BCAST, MPI_REDUCE, MPI_INIT, MPI_FINALIZE, &
              MPI_COMM_RANK, MPI_COMM_SIZE, MPI_SENDRECV_REPLACE

contains

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 function MPI_WTIME() result (x)
    real(kind=KR2) :: x
    x = 0.0_KR2
 end function MPI_WTIME

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 subroutine MPI_BCAST(plaq,a,b,c,d,ierr)
    real(kind=KR), intent(inout) :: plaq
    integer(kind=KI), intent(in) :: a,b,c,d,ierr
    plaq = 0.0_KR*real(a+b+c+d+ierr,KR)
 end subroutine MPI_BCAST

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 subroutine MPI_REDUCE(plaq,plaqtot,a,b,c,d,e,ierr)
    real(kind=KR), intent(in)    :: plaq
    real(kind=KR), intent(out)   :: plaqtot
    integer(kind=KI), intent(in) :: a,b,c,d,e,ierr
    plaqtot = 0.0_KR*plaq*real(a+b+c+d+e+ierr)
 end subroutine MPI_REDUCE

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 subroutine MPI_INIT(ierr)
    integer(kind=KI), intent(out) :: ierr
    ierr = 0
 end subroutine MPI_INIT

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 subroutine MPI_FINALIZE(ierr)
    integer(kind=KI), intent(out) :: ierr
    ierr = 0
 end subroutine MPI_FINALIZE

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 subroutine MPI_COMM_RANK(a,b,ierr)
    integer(kind=KI), intent(in)  :: a
    integer(kind=KI), intent(out) :: b,ierr
    b = 0
    ierr = 0*a
 end subroutine MPI_COMM_RANK

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 subroutine MPI_COMM_SIZE(a,b,ierr)
    integer(kind=KI), intent(in)  :: a
    integer(kind=KI), intent(out) :: b,ierr
    b = 1
    ierr = 0*a
 end subroutine MPI_COMM_SIZE

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 subroutine MPI_SENDRECV_REPLACE(a,b,c,d,e,f,g,h,istatus,ierr)
    real(kind=KR),    intent(inout)             :: a
    integer(kind=KI), intent(in)                :: b, c, d, e, f, g, h
    integer(kind=KI), intent(out), dimension(:) :: istatus
    integer(kind=KI), intent(out)               :: ierr
    ierr = 0*int(a)*b*c*d*e*f*g*h
    istatus = 0
 end subroutine MPI_SENDRECV_REPLACE

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 end module MPI

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
