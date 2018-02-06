!
! operator.f90
!
  module operator
  
  use kinds
  use commonkind
  use basics
    
    real(kind=KR), dimension(nxyzt)          :: psir,psii,rhor,rhoi,&
                                                j1pr,j1pi,j2pr,j2pi,&
                                                j3pr,j3pi,psur,psui,&
                                                arhor,arhoi,&
                                                aj1pr,aj1pi,aj2pr,aj2pi,&
                                                aj3pr,aj3pi
   real(kind=KR), allocatable, dimension(:,:,:) :: z2,z3

  contains

! These subroutine allocate/dealloacte memory for the exact unit vector
! variable in subroutine twvev. When the lattice is large these values
! cause seg fault on a 32-bit Linux machine. Allocation/dealloaction 
! corrects this.

  subroutine allocatez2
 
     if (.not. allocated(z2)) then
        allocate(z2(nxyzt,nc,nd))
     end if
 
     if (.not. allocated(z3)) then
        allocate(z3(nxyzt,nc,nd))
     end if

  end subroutine allocatez2
 
  subroutine deallocatez2
 
    if (allocated(z2)) then
       deallocate(z2)
    end if
 
    if (allocated(z3)) then
       deallocate(z3)
    end if
  end subroutine deallocatez2
   
  end module operator
