!
! multstorage.f90
!
  module multstorage

  use commonkind
  use kinds
  use basics

  public :: allocatemults, deallocatemults
 
  real(kind=KR), public, allocatable, dimension(:,:,:)    :: multsb1r,multsb1i
 
  contains

! These subroutine allocate/dealloacte memory for a large work vector
! in subroutine MULTIPLY in vevbleft.f90. When the lattice is large 
! these values cause seg fault on a 32-bit Linux machine. 
! Allocation/dealloaction corrects this.

  subroutine allocatemults
 
     if (.not. allocated(multsb1r)) then
        allocate(multsb1r(nxyzt,nc,nd))
        allocate(multsb1i(nxyzt,nc,nd))
     end if
 
  end subroutine allocatemults
 
  subroutine deallocatemults
 
    if (allocated(multsb1r)) then
       deallocate(multsb1r)
       deallocate(multsb1i)
    end if
 
  end subroutine deallocatemults

  end module multstorage
