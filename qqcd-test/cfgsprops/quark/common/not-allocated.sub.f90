!
! sub.f90
!
  module sub

  use commonkind
  use kinds
  use basics


  real(kind=KR), public,  dimension(nxyzt,nc,nd)    :: sb1r,sb1i,sb2r,sb2i,&
                                                       sb3r,sb3i,sb4r,sb4i,&
                                                       sb5r,sb5i,sb6r,sb6i

  end module sub
