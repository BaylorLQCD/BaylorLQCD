!
! operator.f90
!
  module operator
  
  use kinds
  use commonkind
  use basics
    
    real(kind=KR), dimension(nxyzt)          :: psir,psii,rhor,rhoi,&
                                                j1pr,j1pi,j2pr,j2pi,&
                                                j3pr,j3pi
     real(kind=KR), dimension(nxyzt,nc,nd)    :: z2

   
  end module operator
