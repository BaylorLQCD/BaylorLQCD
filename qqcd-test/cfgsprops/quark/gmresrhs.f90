! modlue gmresrhs.f90

  module gmresrhs

  use basics
  use kinds
  use latdims
   private

  integer, public   :: isignal
  public :: setrhs, setrhsnoise
  contains

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


  subroutine setrhs(a,b)
  integer(kind=KI), intent(in)   :: a, b 

  if ((a == 1) .and. (b == 1)) then
     isignal = 1
  else
     isignal = 2
  endif
  
  end subroutine setrhs

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  subroutine setrhsnoise(a)
  
  integer(kind=KI), intent(in)  :: a
 
  if (a==1) then
    isignal = 1 
  else
    isignal = 2
  endif ! a==1
  
  end subroutine setrhsnoise
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  end module gmresrhs
