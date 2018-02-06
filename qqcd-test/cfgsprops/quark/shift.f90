! 
! shift.f90
! 
! The parameters public to gmresdrshfit and gmresproject are large
! so they are shared via this module. Equivalent to a common block
! in F77.

  module shift

    use basics
    use latdims
    use kinds


! HEY~ need to make kmax,nshift public to gmres routine......(future clean up)

! Need the next line when doing mulit-twisted mass
!   real(kind=KR2),   dimension(6,ntotal,4,2,8,nshifts)     :: xt,xbt
!   real(kind=KR2),   dimension(6,ntotal,4,2,8,2,nshifts)     :: xvshift
    real(kind=KR2),   dimension(6,ntotal,4,2,8)             :: xb
    real(kind=KR2),   dimension(6,ntotal,4,2,8,nmaxGMRES+1) :: v
    real(kind=KR2),   dimension(6,ntotal,4,2,8,nmaxGMRES+1) :: vprime
    real(kind=KR2),   dimension(6,ntotal,4,2,8,nmaxGMRES+1) :: vprime1
!    real(kind=KR2),   dimension(6,ntotal,4,2,8,nmaxGMRES+1) :: vprime3
!    real(kind=KR2),   dimension(6,ntotal,4,2,8,nmaxGMRES+1) :: vprime4
!    real(kind=KR2),   dimension(6,ntotal,4,2,8,nmaxGMRES+1) :: vprime5
!    real(kind=KR2),   dimension(6,ntotal,4,2,8,nmaxGMRES+1) :: vprime6
    real(kind=KR),   dimension(6,ntotal,4,2,8,nmaxGMRES+1) :: vprime2
    real(kind=KR2),   dimension(6,ntotal,4,2,8,1)           :: y
    real(kind=KR2),   dimension(6,ntotal,4,2,8,nmaxGMRES+1) :: ve, vo
    real(kind=KR2),   dimension(6,ntotal,4,2,8,nmaxGMRES+1) :: evector
    real(kind=KR2),   dimension(6,ntotal,4,2,8,nmaxGMRES+1) :: evectore,evectoro
    real(kind=KR2),   dimension(6,ntotal,4,2,8,nmaxGMRES+1) :: ovo,ove


    real(kind=KR2),   dimension(2,kmaxGMRES)                :: evalue,oneval

    real(kind=KR2),   dimension(6,nvhalf,4,2,8,nshifts)     :: rtemp
    real(kind=KR2),   dimension(6,ntotal,4,2,8,nmaxGMRES+1) :: vtemp
    real(kind=KR2),   dimension(6,nvhalf,4,2,8)           :: beg
    real(kind=KR2),   dimension(2,nmaxGMRES+1,nmaxGMRES)    :: hcnew
    real(kind=KR2),   dimension(nmaxGMRES+1,nmaxGMRES)    :: tcnew
    real(kind=KR2),   dimension(nmaxGMRES)                  :: heigen
    integer(kind=KI), parameter                             :: gdrsize = nxyzt
    integer(kind=KI)                                        :: psignal
    !real(kind=KR2),   dimension(2,nmaxGMRES+1)            :: wtemp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!parameters for p.p!!!!!!!!!!!!!!!!!!!!!
!    integer(kind=KI)                                        :: p!order of p
    integer(kind=KI), dimension(11)                          :: ipiv2!order of p
    real(kind=KR2),   dimension(2,11,11)                      :: lsmat
    real(kind=KR2),   dimension(2,11,1)                      :: cls
    real(kind=KR2),   dimension(2,11)                        :: co,co1!coefficients
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer(kind=KI), dimension(8)                          :: ipiv3!order of p
    real(kind=KR2),   dimension(2,8,8)                      :: lsmat1
    real(kind=KR2),   dimension(2,8,1)                      :: cls1
    real(kind=KR2),   dimension(2,8)                        :: co2,co3!coefficien
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer(kind=KI), dimension(5)                          :: ipiv4!order of p
    real(kind=KR2),   dimension(2,5,5)                      :: lsmat2
    real(kind=KR2),   dimension(2,5,1)                      :: cls2
    real(kind=KR2),   dimension(2,5)                        :: co4,co5!coefficien
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=KR2),   dimension(6,nvhalf,4,2,8)             :: w1, z1 ,z2 
    real(kind=KR2), dimension(6,ntotal,4,2,8,8) :: try!Last entry is the degree of P(A)*A
    real(kind=KR2), dimension(6,ntotal,4,2,8,1) :: test   
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  end module shift
