! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
! debug.f90, Randy Lewis, randy.lewis@uregina.ca
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
! Routines for debugging lattice code.
! I.e. perform a local gauge transformation!
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 module debug

    use kinds
    use latdims
    use basics
    use lagfib
    use gaugetools
    use diracops
    implicit none
    private

! Define access to subroutines.
    public  :: qqcdrot, gaugerot
    private :: ranlinks

 contains

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 
 subroutine qqcdrot(u,be,bo,myid,nn,ldiv,nms,lv,ib,lbd,iblv,vecbl,MRT)
! Perform a local gauge transformation at every lattice site, for both
! gauge fields and propagator sources.
! OUTPUT:
!   The links, u(), and the sources, be() and bo(), are updated.

    real(kind=KR),    intent(inout), dimension(:,:,:,:,:) :: u, be, bo
    integer(kind=KI), intent(in)                          :: myid, MRT
    integer(kind=KI), intent(in),    dimension(:,:)       :: nn, iblv, vecbl
    logical,          intent(in),    dimension(:)         :: ldiv
    integer(kind=KI), intent(in),    dimension(:)         :: nms
    integer(kind=KI), intent(in),    dimension(:,:,:)     :: lv
    integer(kind=KI), intent(in),    dimension(:,:,:,:)   :: ib
    logical,          intent(in),    dimension(:,:)       :: lbd

    real(kind=KR),   dimension(18,ntotal,2,16) :: rot
    real(kind=KR),   dimension(6,nvhalf)       :: tmp1, tmp2
    integer(kind=KI)                           :: ieo, jeo, ibl, jbl, ibleo, &
                                                  mu, isite

! Create a random SU(3) matrix for each lattice site.
    call ranlinks(rot)
! Share those matrices along process boundaries with the neighbouring process.
    do mu = 1,4
     if (ldiv(mu)) then
      do ieo = 1,2
       do ibl = 1,16
        call setlink(rot(:,:,ieo,ibl),ib(:,mu,2,ieo),mu,nms,myid,nn,MRT)
       enddo ! ibl
      enddo ! ieo
     endif
    enddo ! mu
! Multiply each link by the 2 rotation matrices (one at each end of the link).
    call rotlinks(u,rot,myid,nn,ldiv,nms,lv,ib,lbd,iblv,MRT)
! Multiply each source by the appropriate rotation matrix.
    do ieo = 1,2
     jeo = 3 - ieo
     do ibleo = 1,8
      ibl = vecbl(1,ibleo)
      jbl = vecbl(2,ibleo)
      do mu = 1,4
       call mv(nvhalf,rot(:,:,ieo,ibl),be(:,:,mu,ieo,ibleo),tmp1)
       call mv(nvhalf,rot(:,:,ieo,jbl),bo(:,:,mu,ieo,ibleo),tmp2)
       do isite = 1,nvhalf
        be(:,isite,mu,ieo,ibleo) = tmp1(:,isite)
        bo(:,isite,mu,ieo,ibleo) = tmp2(:,isite)
       enddo ! isite
      enddo ! mu
     enddo ! ibl
    enddo ! ieo

 end subroutine qqcdrot

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 
 subroutine gaugerot(u,myid,nn,ldiv,nms,lv,ib,lbd,iblv,MRT)
! Perform a local gauge transformation at every lattice site (no fermions).
! OUTPUT:
!   The links, u(), are updated.

    real(kind=KR),    intent(inout), dimension(:,:,:,:,:) :: u
    integer(kind=KI), intent(in)                          :: myid, MRT
    integer(kind=KI), intent(in),    dimension(:,:)       :: nn
    logical,          intent(in),    dimension(:)         :: ldiv
    integer(kind=KI), intent(in),    dimension(:)         :: nms
    integer(kind=KI), intent(in),    dimension(:,:,:)     :: lv
    integer(kind=KI), intent(in),    dimension(:,:,:,:)   :: ib
    logical,          intent(in),    dimension(:,:)       :: lbd
    integer(kind=KI), intent(in),    dimension(:,:)       :: iblv

    real(kind=KR),   dimension(18,ntotal,2,16) :: rot
    integer(kind=KI)                           :: ieo, ibl, mu

! Create a random SU(3) matrix for each lattice site.
    call ranlinks(rot)
! Share those matrices along process boundaries with the neighbouring process.
    do mu = 1,4
     if (ldiv(mu)) then
      do ieo = 1,2
       do ibl = 1,16
        call setlink(rot(:,:,ieo,ibl),ib(:,mu,2,ieo),mu,nms,myid,nn,MRT)
       enddo ! ibl
      enddo ! ieo
     endif
    enddo ! mu
! Multiply each link by the 2 rotation matrices (one at each end of the link).
    call rotlinks(u,rot,myid,nn,ldiv,nms,lv,ib,lbd,iblv,MRT)

 end subroutine gaugerot

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 
 subroutine ranlinks(rot)
! Build a random SU(3) matrix for each lattice site.
! OUTPUT:
!   rot() contains the SU(3) gauge transformation matrices.
!         expected size: rot(18,ntotal,neo,nblk)

    real(kind=KR),    intent(out), dimension(:,:,:,:) :: rot

    real(kind=KR), dimension(256*nvhalf) :: rn
    real(kind=KR), dimension(18)         :: one
    real(kind=KR)            :: ang1, ang2, ang3, cos4, sin4, ang5, &
                                ang6, cos7, sin7, ang8, rmax, rtry, &
                                tempa, tempb, tempc, tempd, fac1, fac2
    integer(kind=KI)         :: icount, ieo, ibl, isite, icol, jcol, kcol
    real(kind=KR), parameter :: Pi=3.141592653589793_KR

    icount = 256*nvhalf
    call rnds(icount,rn)
    icount = -8
    do ieo = 1,2
     do ibl = 1,16
      do isite = 1,nvhalf
       icount = icount + 8

! Choose the 8 free parameters for this SU(3) matrix.
       ang1 = Pi*rn(icount+1)
       ang2 = Pi*rn(icount+2)
       ang3 = Pi*rn(icount+3)
       cos4 = -1.0_KR + 2.0_KR*rn(icount+4)
       sin4 = sqrt(1.0_KR-cos4**2)
       ang5 = 2.0_KR*Pi*rn(icount+5)
       ang6 = Pi*rn(icount+6)
       cos7 = -1.0_KR + 2.0_KR*rn(icount+7)
       sin7 = sqrt(1.0_KR-cos7**2)
       ang8 = 2.0_KR*Pi*rn(icount+8)
! Calculate the first row of the SU(3) matrix.
       one(1) = cos(ang1)
       one(2) = sin(ang1)*cos(ang2)
       one(3) = sin(ang1)*sin(ang2)*cos(ang3)
       one(4) = sin(ang1)*sin(ang2)*sin(ang3)*cos4
       one(5) = sin(ang1)*sin(ang2)*sin(ang3)*sin4*cos(ang5)
       one(6) = sin(ang1)*sin(ang2)*sin(ang3)*sin4*sin(ang5)
! Choose the entry with the largest magnitude to avoid division by zero
! when calculating the second and third rows.
       rmax = one(1)**2 + one(2)**2
       icol = 3
       jcol = 5
       kcol = 1
       rtry = one(3)**2 + one(4)**2
       if (rmax<rtry) then
        rmax = rtry
        icol = 5
        jcol = 1
        kcol = 3
       endif
       rtry = one(5)**2 + one(6)**2
       if (rmax<rtry) then
        rmax = rtry
        icol = 1
        jcol = 3
        kcol = 5
       endif
! Calculate the second row of the SU(3) matrix.
       tempa = cos(ang6)
       tempb = sin(ang6)*cos7
       tempc = sin(ang6)*sin7*cos(ang8)
       tempd = sin(ang6)*sin7*sin(ang8)
       fac1 = sqrt(1.0_KR-one(jcol)**2-one(jcol+1)**2)
       fac2 = sqrt(one(kcol)**2+one(kcol+1)**2)
       one(icol+6) = ((tempa+tempb)*fac2                                   &
                  -tempc*(one(icol  )*one(jcol  )+one(icol+1)*one(jcol+1)  &
                         +one(icol  )*one(jcol+1)-one(icol+1)*one(jcol  )) &
                  -tempd*(one(icol+1)*one(jcol  )+one(icol  )*one(jcol  )  &
                         +one(icol+1)*one(jcol+1)-one(icol  )*one(jcol+1)) &
                  )/fac1/sqrt(2.0_KR)
       one(icol+7) = ((tempa-tempb)*fac2                                   &
                  -tempc*(one(icol  )*one(jcol  )+one(icol+1)*one(jcol+1)  &
                         -one(icol  )*one(jcol+1)+one(icol+1)*one(jcol  )) &
                  -tempd*(one(icol+1)*one(jcol  )-one(icol  )*one(jcol  )  &
                         -one(icol+1)*one(jcol+1)-one(icol  )*one(jcol+1)) &
                  )/fac1/sqrt(2.0_KR)
       one(jcol+6) = (tempc+tempd)*fac1/sqrt(2.0_KR)
       one(jcol+7) = (tempc-tempd)*fac1/sqrt(2.0_KR)
       one(kcol+6) = -1.0_KR/fac2**2*(                                      &
                (one(icol)*one(kcol  )+one(icol+1)*one(kcol+1))*one(icol+6) &
               -(one(icol)*one(kcol+1)-one(icol+1)*one(kcol  ))*one(icol+7) &
               +(one(jcol)*one(kcol  )+one(jcol+1)*one(kcol+1))*one(jcol+6) &
               -(one(jcol)*one(kcol+1)-one(jcol+1)*one(kcol  ))*one(jcol+7) )
       one(kcol+7) = -1.0_KR/fac2**2*(                                      &
                (one(icol)*one(kcol  )+one(icol+1)*one(kcol+1))*one(icol+7) &
               +(one(icol)*one(kcol+1)-one(icol+1)*one(kcol  ))*one(icol+6) &
               +(one(jcol)*one(kcol  )+one(jcol+1)*one(kcol+1))*one(jcol+7) &
               +(one(jcol)*one(kcol+1)-one(jcol+1)*one(kcol  ))*one(jcol+6) )
! Calculate the third row of the SU(3) matrix.
       one(icol+12) = one(jcol)*one(kcol+6)-one(jcol+1)*one(kcol+7) &
                     -one(kcol)*one(jcol+6)+one(kcol+1)*one(jcol+7)
       one(icol+13) = one(kcol)*one(jcol+7)+one(kcol+1)*one(jcol+6) &
                     -one(jcol)*one(kcol+7)-one(jcol+1)*one(kcol+6)
       one(jcol+12) = one(kcol)*one(icol+6)-one(kcol+1)*one(icol+7) &
                     -one(icol)*one(kcol+6)+one(icol+1)*one(kcol+7)
       one(jcol+13) = one(icol)*one(kcol+7)+one(icol+1)*one(kcol+6) &
                     -one(kcol)*one(icol+7)-one(kcol+1)*one(icol+6)
       one(kcol+12) = one(icol)*one(jcol+6)-one(icol+1)*one(jcol+7) &
                     -one(jcol)*one(icol+6)+one(jcol+1)*one(icol+7)
       one(kcol+13) = one(jcol)*one(icol+7)+one(jcol+1)*one(icol+6) &
                     -one(icol)*one(jcol+7)-one(icol+1)*one(jcol+6)
! Include this SU(3) matrix into the main collection.
       rot(:,isite,ieo,ibl) = one(:)

      enddo ! isite
     enddo ! ibl
    enddo ! ieo

 end subroutine ranlinks

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 end module debug

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
