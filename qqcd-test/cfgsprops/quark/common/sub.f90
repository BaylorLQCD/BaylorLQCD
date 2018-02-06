!
! sub.f90
!
  module sub

  use commonkind
  use kinds
  use basics

  public :: allocatesubs, deallocatesubs, allocateus, deallocateus

  real(kind=KR), public, allocatable, dimension(:,:,:,:,:) :: uss
  real(kind=KR), public, allocatable, dimension(:,:,:,:)   :: usp
  real(kind=KR), public, allocatable, dimension(:,:,:)     :: sb0r,sb0i,sb1r,sb1i,&
                                                              sb2r,sb2i,sb3r,sb3i,&
                                                              sb4r,sb4i,sb5r,sb5i,&
                                                              sb6r,sb6i

  contains

! These subroutine allocate/dealloacte memory for the subraction
! variables. When the lattice is large these values cause seg fault
! on a 32-bit Linux machine. Allocation/dealloaction corrects this. 



  subroutine allocateus

  if (.not. allocated(uss)) then
     allocate(uss(nxyzt,3,2,nc,nc))
  endif

  if (.not. allocated(usp)) then
     allocate(usp(nxyzt,2,nc,nc))
  endif ! allocate

  end subroutine allocateus

  subroutine deallocateus
  if (allocated(uss)) then
     deallocate(uss)
  endif

  if (allocated(usp)) then
     deallocate(usp)
  endif ! allocate
  end subroutine deallocateus

  subroutine allocatesubs
 
     if (.not. allocated(sb1r)) then
        allocate(sb0r(nxyzt,nc,nd))
        allocate(sb0i(nxyzt,nc,nd))
        allocate(sb1r(nxyzt,nc,nd))
        allocate(sb1i(nxyzt,nc,nd))
        allocate(sb2r(nxyzt,nc,nd))
        allocate(sb2i(nxyzt,nc,nd))
        allocate(sb3r(nxyzt,nc,nd))
        allocate(sb3i(nxyzt,nc,nd))
        allocate(sb4r(nxyzt,nc,nd))
        allocate(sb4i(nxyzt,nc,nd))
        allocate(sb5r(nxyzt,nc,nd))
        allocate(sb5i(nxyzt,nc,nd))
        allocate(sb6r(nxyzt,nc,nd))
        allocate(sb6i(nxyzt,nc,nd))
     endif ! allocate

  end subroutine allocatesubs
 
  subroutine deallocatesubs
 
    if (allocated(sb1r)) then
       deallocate(sb0r)
       deallocate(sb0i)
       deallocate(sb1r)
       deallocate(sb1i)
       deallocate(sb2r)
       deallocate(sb2i)
       deallocate(sb3r)
       deallocate(sb3i)
       deallocate(sb4r)
       deallocate(sb4i)
       deallocate(sb5r)
       deallocate(sb5i)
       deallocate(sb6r)
       deallocate(sb6i)
    endif ! deallocate

  end subroutine deallocatesubs
 
  end module sub
