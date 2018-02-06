! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
! gaugetools.f90, Randy Lewis, randy.lewis@uregina.ca
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
! Routines for manipulation of gauge field links.
!
! parts of this code are based on:
! *QCDMPI (version 5, 1996) by Shinji Hioki, Parallel Computing 22 (1997) 1335.
! and
! *QCDimMPI (version 0.9, 1998) by Astushi Nakamura and Shinji Hioki,
!           Nucl Phys B (Proc Suppl) 73 (1999) 895.
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 module gaugetools

 !   use MPI
    use kinds
    use latdims
    use basics
    implicit none
    private

! Use the following line if the MPI module is not available.
    include 'mpif.h'

! Define access to subroutines.
    public  :: staple, loop, makeclover, setmmf, mm, mmd, mmdg, setlink, &
               GFwrite, GFread, aveplaq, tadPcalc, tadLcalc, fuzz, &
               wilsonloop, rotlinks
    private :: seteo, mmg, mdm, slidematrix, trace, fixsu3, gfixL

 contains

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 subroutine staple(mu,ieo,ibl,coact,gaction,c,u,myid,nn,ldiv,nms,lv,ib,lbd, &
                   iblv,MRT)
! Calculate the staples which accompany a mu-direction link.
! INPUT:
!   mu is the spacetime direction of the links for which staples are required.
!   ieo tells whether even(ieo=1) or odd(ieo=2) links need staples computed.
!   ibl is the block of the sublattice containing the links of interest.
!   coact(2,4,4) contains the coefficients from the action.
!   gaction(i) = 1 for nearest neighbour action in the i'th spacetime direction.
!   gaction(i) = 2 for next-nearest neighbour action in the i'th direction.
! OUTPUT:
!   c(18,nvhalf) contains the staples.
!
! NOTATION:
!    +------+                +------+------+
!    |      |                |             |
!    |      |                |             |
!    o      +    o      +    o      +------+   o       ------+
!                |      |                      |             |
!                |      |                      |             |
!                +------+                      +------+------+
!
!    staple 1    staple 2       staple 3          staple 4
!
!                                          +------+
!                                          |      |
!                                          |      |
!    +------+------+                       +      +
!    |             |                       |      |
!    |             |                       |      |
!    +------o      +   +------o      +     o      +     o      +
!                      |             |                  |      |
!                      |             |                  |      |
!                      +------+------+                  +      +
!                                                       |      |
!                                                       |      |
!                                                       +------+
!       staple 5          staple 6         staple 7     staple 8
!
!   where the coordinate system is...
!
!   nu
!   /\
!    |
!    |
!    |____\
!         /mu

    integer(kind=KI), intent(in)                          :: mu, ieo, ibl, MRT
    integer(kind=KI), intent(in),    dimension(:)         :: gaction
    real(kind=KR),    intent(in),    dimension(:,:,:)     :: coact
    real(kind=KR),    intent(out),   dimension(:,:)       :: c
    real(kind=KR),    intent(inout), dimension(:,:,:,:,:) :: u
    integer(kind=KI), intent(in)                          :: myid
    integer(kind=KI), intent(in),    dimension(:,:)       :: nn, iblv
    logical,          intent(in),    dimension(:)         :: ldiv
    integer(kind=KI), intent(in),    dimension(:)         :: nms
    integer(kind=KI), intent(in),    dimension(:,:,:)     :: lv
    integer(kind=KI), intent(in),    dimension(:,:,:,:)   :: ib
    logical,          intent(in),    dimension(:,:)       :: lbd

    integer(kind=KI) :: jeo, ieo0, jeo0, ieo1, jeo1, nu, mu0, nu0, i, ibl0, &
                        ibl1, icount
    real(kind=KR), dimension(18,nvhalf) :: t
    real(kind=KR), dimension(18,ntotal) :: w, x

    icount = nvhalf
    jeo = 3 - ieo
    c = 0.0_KR
 
    do nu = 1,4
     if (nu/=mu) then
 
!*staple 1

! Multiply first two links (from initial site go in +nu direction, then +mu).
      call setmmf(nu,lbd(ibl,nu),u(:,:,mu,ieo,iblv(ibl,nu))           &
                 ,ldiv(nu),ib(:,nu,1,jeo),u(:,:,nu,ieo,ibl)                  &
                 ,u(:,:,mu,jeo,iblv(ibl,nu)),t,lv(:,nu,ieo),1,nms,myid,nn,MRT)
! Multiply in the third link (dagger of link in +nu direction).
      call setmmf(mu,lbd(ibl,mu),u(:,:,nu,ieo,iblv(ibl,mu))           &
                 ,ldiv(mu),ib(:,mu,1,jeo),t                                  &
                 ,u(:,:,nu,jeo,iblv(ibl,mu)),w,lv(:,mu,ieo),2,nms,myid,nn,MRT)
! Add this staple to the total.
      do i = 1,nvhalf
       c(:,i) = c(:,i) + coact(1,mu,nu)*w(:,i)
      enddo ! i

!*staple 2
! Note: in this case each process will actually compute its neighbour's
!       staple, and then they will be moved to the neighbour.

      ibl0 = iblv(ibl,nu)
      call seteo(lbd(ibl,nu),jeo,ieo,ieo0,jeo0)
! Multiply first two links (from initial site go in -nu direction, then +mu).
      call mdm(icount,u(:,:,nu,ieo0,ibl0),u(:,:,mu,ieo0,ibl0),w)
! Multiply in the third link (link in +nu direction).
      call setmmf(mu,lbd(ibl0,mu),u(:,:,nu,ieo0,iblv(ibl0,mu))               &
                 ,ldiv(mu),ib(:,mu,1,jeo0),w                                 &
               ,u(:,:,nu,jeo0,iblv(ibl0,mu)),t,lv(:,mu,ieo0),1,nms,myid,nn,MRT)

! If this is the first of two blocks in this direction, then move this
! staple to its true owner (neighbouring process in the +nu direction) 
! using slidematrix.  In any case, add the staple to the total.
      if (lbd(ibl,nu)) then
       do i = 1,nvhalf
        w(:,lv(i,nu,jeo)) = t(:,i)
       enddo ! i
       if (ldiv(nu)) call slidematrix(w,ib(:,nu,1,ieo),nu,nms,myid,nn,MRT)
       do i = 1,nvhalf
        c(:,i) = c(:,i) + coact(1,mu,nu)*w(:,i)
       enddo ! i
      else
       do i = 1,nvhalf
        c(:,i) = c(:,i) + coact(1,mu,nu)*t(:,i)
       enddo ! i
      endif

!*staples 3, 4, 5 and 6
      if (gaction(mu)/=1) then
! staple 3
       ibl0 = iblv(ibl,mu)
       call seteo(lbd(ibl,mu),ieo,jeo,ieo0,jeo0)
       mu0 = nu
       nu0 = mu
       call setmmf(nu0,lbd(ibl0,nu0),u(:,:,mu0,ieo0,iblv(ibl0,nu0))          &
                  ,ldiv(nu0),ib(:,nu0,1,jeo0),u(:,:,nu0,ieo0,ibl0)           &
           ,u(:,:,mu0,jeo0,iblv(ibl0,nu0)),t,lv(:,nu0,ieo0),1,nms,myid,nn,MRT)
       call setmmf(mu0,lbd(ibl0,mu0),u(:,:,nu0,ieo0,iblv(ibl0,mu0))          &
                  ,ldiv(mu0),ib(:,mu0,1,jeo0),t                              &
           ,u(:,:,nu0,jeo0,iblv(ibl0,mu0)),w,lv(:,mu0,ieo0),2,nms,myid,nn,MRT)
       call setmmf(nu,lbd(ibl,nu),u(:,:,mu,ieo,iblv(ibl,nu))                 &
                  ,ldiv(nu),ib(:,nu,1,jeo),u(:,:,nu,ieo,ibl)                 &
                 ,u(:,:,mu,jeo,iblv(ibl,nu)),t,lv(:,nu,ieo),1,nms,myid,nn,MRT)
       call setmmf(mu,lbd(ibl,mu),w                                          &
                  ,ldiv(mu),ib(:,mu,1,jeo),t                                 &
                  ,w,x,lv(:,mu,ieo),2,nms,myid,nn,MRT)
       do i = 1,nvhalf
        c(:,i) = c(:,i) + coact(2,mu,nu)*x(:,i)
       enddo ! i
! staple 4
       ibl0 = iblv(iblv(ibl,nu),mu)
       call seteo(lbd(ibl,nu),jeo,ieo,ieo1,jeo1)
       call seteo(lbd(ibl,mu),ieo1,jeo1,ieo0,jeo0)
       mu0 = nu
       nu0 = mu
       call setmmf(nu0,lbd(ibl0,nu0),u(:,:,mu0,ieo0,iblv(ibl0,nu0))          &
                  ,ldiv(nu0),ib(:,nu0,1,jeo0),u(:,:,nu0,ieo0,ibl0)           &
           ,u(:,:,mu0,jeo0,iblv(ibl0,nu0)),t,lv(:,nu0,ieo0),1,nms,myid,nn,MRT)
       call setmmf(mu0,lbd(ibl0,mu0),u(:,:,nu0,ieo0,iblv(ibl0,mu0))          &
                  ,ldiv(mu0),ib(:,mu0,1,jeo0),t                              &
           ,u(:,:,nu0,jeo0,iblv(ibl0,mu0)),w,lv(:,mu0,ieo0),2,nms,myid,nn,MRT)
       ibl0 = iblv(ibl,nu)
       call seteo(lbd(ibl,nu),jeo,ieo,ieo0,jeo0)
       call mdm(icount,u(:,:,nu,ieo0,ibl0),u(:,:,mu,ieo0,ibl0),x)
       call setmmf(mu,lbd(ibl0,mu),w                                         &
                  ,ldiv(mu),ib(:,mu,1,jeo0),x                                &
                  ,w,t,lv(:,mu,ieo0),1,nms,myid,nn,MRT)
       if (lbd(ibl,nu)) then
        do i = 1,nvhalf
         w(:,lv(i,nu,jeo)) = t(:,i)
        enddo ! i
        if (ldiv(nu)) call slidematrix(w,ib(:,nu,1,ieo),nu,nms,myid,nn,MRT)
        do i = 1,nvhalf
         c(:,i) = c(:,i) + coact(2,mu,nu)*w(:,i)
        enddo ! i
       else
        do i = 1,nvhalf
         c(:,i) = c(:,i) + coact(2,mu,nu)*t(:,i)
        enddo ! i
       endif
 
! staple 5
       ibl0 = iblv(ibl,mu)
       call seteo(lbd(ibl,mu),jeo,ieo,ieo0,jeo0)
       mu0 = nu
       nu0 = mu
       call mdm(icount,u(:,:,nu0,ieo0,ibl0),u(:,:,mu0,ieo0,ibl0),w)
       call setmmf(mu0,lbd(ibl0,mu0),u(:,:,nu0,ieo0,iblv(ibl0,mu0))          &
                  ,ldiv(mu0),ib(:,mu0,1,jeo0),w                              &
           ,u(:,:,nu0,jeo0,iblv(ibl0,mu0)),t,lv(:,mu0,ieo0),1,nms,myid,nn,MRT)
       if (lbd(ibl,nu0)) then
        do i = 1,nvhalf
         w(:,lv(i,nu0,jeo)) = t(:,i)
        enddo ! i
        if (ldiv(nu0)) call slidematrix(w,ib(:,nu0,1,ieo),nu0,nms,myid,nn,MRT)
        call setmmf(nu,lbd(ibl,nu),u(:,:,mu,ieo,iblv(ibl,nu))                &
                   ,ldiv(nu),ib(:,nu,1,jeo),w                                &
                 ,u(:,:,mu,jeo,iblv(ibl,nu)),x,lv(:,nu,ieo),1,nms,myid,nn,MRT)
       else
        call setmmf(nu,lbd(ibl,nu),u(:,:,mu,ieo,iblv(ibl,nu))                &
                   ,ldiv(nu),ib(:,nu,1,jeo),t                                &
                 ,u(:,:,mu,jeo,iblv(ibl,nu)),x,lv(:,nu,ieo),1,nms,myid,nn,MRT)
       endif
       call setmmf(mu,lbd(ibl,mu),u(:,:,nu,ieo,iblv(ibl,mu))                 &
                  ,ldiv(mu),ib(:,mu,1,jeo),x                                 &
                 ,u(:,:,nu,jeo,iblv(ibl,mu)),w,lv(:,mu,ieo),2,nms,myid,nn,MRT)
       do i = 1,nvhalf
        c(:,i) = c(:,i) + coact(2,mu,nu)*w(:,i)
       enddo ! i
 
! staple 6
       ibl0 = iblv(iblv(ibl,mu),nu)
       call seteo(lbd(ibl,nu),jeo,ieo,ieo1,jeo1)
       call seteo(lbd(ibl,mu),jeo1,ieo1,ieo0,jeo0)
       mu0 = nu
       nu0 = mu
       call mdm(icount,u(:,:,nu0,ieo0,ibl0),u(:,:,mu0,ieo0,ibl0),w)
       call setmmf(mu0,lbd(ibl0,mu0),u(:,:,nu0,ieo0,iblv(ibl0,mu0))          &
                  ,ldiv(mu0),ib(:,mu0,1,jeo0),w                              &
           ,u(:,:,nu0,jeo0,iblv(ibl0,mu0)),t,lv(:,mu0,ieo0),1,nms,myid,nn,MRT)
       if (lbd(ibl,nu0)) then
        do i = 1,nvhalf
         w(:,lv(i,nu0,jeo1)) = t(:,i)
        enddo ! i
        if (ldiv(nu0)) call slidematrix(w,ib(:,nu0,1,ieo1),nu0,nms,myid,nn,MRT)
       endif
       ibl0 = iblv(ibl,nu)
       call seteo(lbd(ibl,nu),jeo,ieo,ieo0,jeo0)
       if (lbd(ibl,nu0)) then
        call mdm(icount,w,u(:,:,mu,ieo0,ibl0),x)
       else
        call mdm(icount,t,u(:,:,mu,ieo0,ibl0),x)
       endif
       call setmmf(mu,lbd(ibl0,mu),u(:,:,nu,ieo0,iblv(ibl0,mu))              &
                  ,ldiv(mu),ib(:,mu,1,jeo0),x                                &
              ,u(:,:,nu,jeo0,iblv(ibl0,mu)),t,lv(:,mu,ieo0),1,nms,myid,nn,MRT)
       if (lbd(ibl,nu)) then
        do i = 1,nvhalf
         w(:,lv(i,nu,jeo)) = t(:,i)
        enddo ! i
        if (ldiv(nu)) call slidematrix(w,ib(:,nu,1,ieo),nu,nms,myid,nn,MRT)
         do i = 1,nvhalf
          c(:,i) = c(:,i) + coact(2,mu,nu)*w(:,i)
         enddo ! i
       else
        do i = 1,nvhalf
         c(:,i) = c(:,i) + coact(2,mu,nu)*t(:,i)
        enddo ! i
       endif
      endif
 
!*staples 7 and 8
      if (gaction(nu)/=1) then
! staple 7
       ibl0 = iblv(ibl,nu)
       call seteo(lbd(ibl,nu),ieo,jeo,ieo0,jeo0)
       call setmmf(nu,lbd(ibl0,nu),u(:,:,mu,ieo0,iblv(ibl0,nu))              &
                  ,ldiv(nu),ib(:,nu,1,jeo0),u(:,:,nu,ieo0,ibl0)              &
              ,u(:,:,mu,jeo0,iblv(ibl0,nu)),t,lv(:,nu,ieo0),1,nms,myid,nn,MRT)
       call setmmf(mu,lbd(ibl0,mu),u(:,:,nu,ieo0,iblv(ibl0,mu))              &
                  ,ldiv(mu),ib(:,mu,1,jeo0),t                                &
              ,u(:,:,nu,jeo0,iblv(ibl0,mu)),w,lv(:,mu,ieo0),2,nms,myid,nn,MRT)
       call setmmf(nu,lbd(ibl,nu),w                                          &
                  ,ldiv(nu),ib(:,nu,1,jeo),u(:,:,nu,ieo,ibl)                 &
                  ,w,t,lv(:,nu,ieo),1,nms,myid,nn,MRT)
       call setmmf(mu,lbd(ibl,mu),u(:,:,nu,ieo,iblv(ibl,mu))                 &
                  ,ldiv(mu),ib(:,mu,1,jeo),t                                 &
                 ,u(:,:,nu,jeo,iblv(ibl,mu)),w,lv(:,mu,ieo),2,nms,myid,nn,MRT)
       do i = 1,nvhalf
        c(:,i) = c(:,i) + coact(2,nu,mu)*w(:,i)
       enddo ! i
 
! staple 8
       ibl1 = iblv(ibl,nu)
       call seteo(lbd(ibl,nu),jeo,ieo,ieo1,jeo1)
       ibl0 = iblv(ibl1,nu)
       call seteo(lbd(ibl1,nu),jeo1,ieo1,ieo0,jeo0)
       call mdm(icount,u(:,:,nu,ieo0,ibl0),u(:,:,mu,ieo0,ibl0),w)
       call setmmf(mu,lbd(ibl0,mu),u(:,:,nu,ieo0,iblv(ibl0,mu))              &
                  ,ldiv(mu),ib(:,mu,1,jeo0),w                                &
              ,u(:,:,nu,jeo0,iblv(ibl0,mu)),t,lv(:,mu,ieo0),1,nms,myid,nn,MRT)
       if (.not.lbd(ibl,nu)) then
        do i = 1,nvhalf
         w(:,lv(i,nu,jeo1)) = t(:,i)
        enddo ! i
        if (ldiv(nu)) call slidematrix(w,ib(:,nu,1,ieo1),nu,nms,myid,nn,MRT)
       endif
       ibl0 = iblv(ibl,nu)
       call seteo(lbd(ibl,nu),jeo,ieo,ieo0,jeo0)
       if (.not.lbd(ibl,nu)) then
        call mdm(icount,u(:,:,nu,ieo0,ibl0),w,x)
       else
        call mdm(icount,u(:,:,nu,ieo0,ibl0),t,x)
       endif
       call setmmf(mu,lbd(ibl0,mu),u(:,:,nu,ieo0,iblv(ibl0,mu))              &
                  ,ldiv(mu),ib(:,mu,1,jeo0),x                                &
              ,u(:,:,nu,jeo0,iblv(ibl0,mu)),t,lv(:,mu,ieo0),1,nms,myid,nn,MRT)
       if (lbd(ibl,nu)) then
        do i = 1,nvhalf
         w(:,lv(i,nu,jeo)) = t(:,i)
        enddo ! i
         if (ldiv(nu)) call slidematrix(w,ib(:,nu,1,ieo),nu,nms,myid,nn,MRT)
         do i = 1,nvhalf
          c(:,i) = c(:,i) + coact(2,nu,mu)*w(:,i)
         enddo ! i
       else
         do i = 1,nvhalf
          c(:,i) = c(:,i) + coact(2,nu,mu)*t(:,i)
         enddo ! i
       endif
      endif
 
     endif
    enddo ! nu
 
 end subroutine staple

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 
 subroutine loop(mu,nu,ieo,ibl,c,u,myid,nn,ldiv,nms,lv,ib,lbd,iblv,MRT)
! Calculate the 2 nearest-neighbour staples which accompany a mu-direction link.
! INPUT:
!   mu is the spacetime direction of the links for which staples are required.
!   nu: if nu=mu then sum all 3 pairs of staples in the mu-anything planes.
!       if nu<>mu then only sum the two staples in the mu-nu plane.
!   ieo tells whether even(ieo=1) or odd(ieo=2) links need staples computed.
!   ibl is the block of the sublattice containing the links of interest.
!   MRT is MPIREALTYPE.
! OUTPUT:
!   c contains the staples.  expected size: c(18,nvhalf)
 
    integer(kind=KI), intent(in)                 :: mu, nu, ieo, ibl, myid, MRT
    real(kind=KR),    intent(out),   dimension(:,:)       :: c
    real(kind=KR),    intent(inout), dimension(:,:,:,:,:) :: u
    integer(kind=KI), intent(in),    dimension(:,:)       :: nn
    logical,          intent(in),    dimension(:)         :: ldiv
    integer(kind=KI), intent(in),    dimension(:)         :: nms
    integer(kind=KI), intent(in),    dimension(:,:,:)     :: lv
    integer(kind=KI), intent(in),    dimension(:,:,:,:)   :: ib
    logical,          intent(in),    dimension(:,:)       :: lbd
    integer(kind=KI), intent(in),    dimension(:,:)       :: iblv
 
    integer(kind=KI)                :: jeo, ieo0, jeo0, i, ibl0, lambda, icount
    real(kind=KR),   dimension(18,nvhalf) :: t
    real(kind=KR),   dimension(18,ntotal) :: w
 
    icount = nvhalf
    jeo = 3 - ieo
    c = 0.0_KR
 
    do lambda = 1,4
     if (lambda/=mu.and.(mu==nu.or.lambda==nu)) then
 
!*staple 1
 
! Multiply first two links (from initial site go in +lambda dir, then +mu).
      call setmmf(lambda,lbd(ibl,lambda),u(:,:,mu,ieo,iblv(ibl,lambda))   &
                 ,ldiv(lambda),ib(:,lambda,1,jeo),u(:,:,lambda,ieo,ibl)   &
                 ,u(:,:,mu,jeo,iblv(ibl,lambda)),t,lv(:,lambda,ieo),1,nms &
                 ,myid,nn,MRT)
! Multiply in the third link (dagger of link in +lambda direction).
      call setmmf(mu,lbd(ibl,mu),u(:,:,lambda,ieo,iblv(ibl,mu))        &
                 ,ldiv(mu),ib(:,mu,1,jeo),t                            &
                 ,u(:,:,lambda,jeo,iblv(ibl,mu)),w,lv(:,mu,ieo),2,nms  &
                 ,myid,nn,MRT)
! Add this staple to the total.
      do i = 1,nvhalf
       c(:,i) = c(:,i) + w(:,i)
      enddo ! i
 
!*staple 2
! Note: in this case each process will actually compute its neighbour's
!       staple, and then they will be moved to the neighbour.
      ibl0 = iblv(ibl,lambda)
      call seteo(lbd(ibl,lambda),jeo,ieo,ieo0,jeo0)
! Multiply first two links (from initial site go in -lambda dir, then +mu).
      call mdm(icount,u(:,:,lambda,ieo0,ibl0),u(:,:,mu,ieo0,ibl0),w)
! Multiply in the third link (link in +lambda direction).
      call setmmf(mu,lbd(ibl0,mu),u(:,:,lambda,ieo0,iblv(ibl0,mu))        &
                 ,ldiv(mu),ib(:,mu,1,jeo0),w                              &
                 ,u(:,:,lambda,jeo0,iblv(ibl0,mu)),t,lv(:,mu,ieo0),1,nms  &
                 ,myid,nn,MRT)
! If this is the first of two blocks in this direction, then move this
! staple to its true owner (neighbouring process in the -lambda direction)
! using slidematrix.  In any case, add the staple to the total.
      if (lbd(ibl,lambda)) then
       do i = 1,nvhalf
        w(:,lv(i,lambda,jeo)) = t(:,i)
       enddo ! i
       if (ldiv(lambda)) call slidematrix(w,ib(:,lambda,1,ieo),lambda,nms, &
                                          myid,nn,MRT)
       do i = 1,nvhalf
        c(:,i) = c(:,i) + w(:,i)
       enddo ! i
      else
       do i = 1,nvhalf
        c(:,i) = c(:,i) + t(:,i)
       enddo ! i
      endif
 
     endif
    enddo ! lambda
 
 end subroutine loop

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 subroutine makeclover(sigF,coact,u,myid,nn,ldiv,nms,lv,ib,lbd,iblv,MRT)
! Construct the clover term in the fermion matrix.
! INPUT:
!   u() contains the gauge fields for this sublattice.
! OUTPUT:
!   sigF(x) = sum_(mu,nu) coact(4,mu,nu)*(-i)*sigma_(mu,nu)*F_(mu,nu)(x)
!   where -i*sum_(mu<nu) sigma_(mu,nu)*F_(mu,nu) = /   V  i*W \
!                                                  \ -i*W  V  /
!         with V = i*sig3*F_(1,2) - i*sig2*F_(1,3) + i*sig1*F_(2,3)
!              W = i*sig1*F_(1,4) + i*sig2*F_(2,4) + i*sig3*F_(3,4)
!              i = sqrt(-1)
!              sig1 = / 0 1 \   sig2 = / 0 -i \   sig3 = / 1  0 \
!                     \ 1 0 /          \ i  0 /          \ 0 -1 /
!   expected size: sigF(18,nvhalf,8,2,16)
!      a=1,18 runs over Re and Im parts of the 3x3 colour matrix as follows:
!                / sF(1,i)+i*sF(2,i)   sF(3,i)+i*sF(4,i)   sF(5,i)+i*sF(6,i) \
!          sF = |  sF(7,i)+i*sF(8,i)   sF(9,i)+i*sF(10,i) sF(11,i)+i*sF(12,i) |
!                \sF(13,i)+i*sF(14,i) sF(15,i)+i*sF(16,i) sF(17,i)+i*sF(18,i)/
!      b=1,nvhalf runs over even or odd sites on one block of a sublattice.
!      c=1,8 is the 4x4 Dirac structure (4 entries in each of V and W above,
!            in this order: V11, V12, V21, V22, W11, W12, W21, W22.)
!      d=1,2 runs over the even(1) and odd(2) lattice sites.
!      e=1,16 runs over the blocks of a sublattice.
! References: Luscher, Sint, Sommer & Weisz, hep-lat/9605038.
!             Gockeler et al., NPB(PS)53, 312(1997) and PRD57, 5562 (1998).
!             Lewis & Woloshyn, PRD58, 074506 (1998).
!             Lepage, Magnea, Nakhleh, Magnea & Hornbostel, PRD46, 4052 (1992).
!
! NOTATION:
!    +------+ +------+
!    |  2   | |  1   |
!    |      | v      |     nu
!    +--->  o o------+     /\
!    +------o o  <---+      |
!    |  3   ^ |  4   |      |
!    |      | |      |      |____\
!    +------+ +------+           /mu
!
! Leaf 1 is computed locally.
! Leaf 2 is computed at a neighbouring site, then moved +mu.
! Leaf 3 is computed at a neighbouring site, then moved +mu+nu.
! Leaf 4 is computed at a neighbouring site, then moved +nu.

    real(kind=KR),    intent(out),   dimension(:,:,:,:,:) :: sigF
    real(kind=KR),    intent(in),    dimension(:,:,:)     :: coact
    real(kind=KR),    intent(inout), dimension(:,:,:,:,:) :: u
    integer(kind=KI), intent(in)                          :: myid, MRT
    integer(kind=KI), intent(in),    dimension(:,:)       :: nn
    logical,          intent(in),    dimension(:)         :: ldiv
    integer(kind=KI), intent(in),    dimension(:)         :: nms
    integer(kind=KI), intent(in),    dimension(:,:,:)     :: lv
    integer(kind=KI), intent(in),    dimension(:,:,:,:)   :: ib
    logical,          intent(in),    dimension(:,:)       :: lbd
    integer(kind=KI), intent(in),    dimension(:,:)       :: iblv

    real(kind=KR), dimension(18,ntotal) :: leaftot, temp1, temp2
    real(kind=KR)                       :: sgnVW
    integer(kind=KI)                    :: i, j, icc, jcc, icount, ieo, jeo, &
                                           ieo0, jeo0, ieo1, jeo1, ibl, ibl0, &
                                           ibl1, iplane, mu, nu
    integer(kind=KI),  dimension(2)     :: VorW
    integer(kind=KI),  dimension(2,9)   :: kcc

! Initialize matrices that will only get defined via a subroutine call.
    temp1 = 0.0_KR
    temp2 = 0.0_KR
    leaftot = 0.0_KR

! Set some parameters and initiate the outer iterations.
    icount = nvhalf
    kcc(1,1:9) = (/ 1, 7, 13, 3, 9, 15, 5, 11, 17 /)
    kcc(2,1:9) = (/ 2, 8, 14, 4, 10, 16, 6, 12, 18 /)
    sigF = 0.0_KR
    do ieo = 1,2
     jeo = 3 - ieo
     do ibl = 1,16
      iplane = 0
      do mu = 1,3
       do nu = mu+1,4
        iplane = iplane + 1

! Leaf 1.
        call setmmf(mu,lbd(ibl,mu),u(:,:,nu,ieo,iblv(ibl,mu)),ldiv(mu), &
                    ib(:,mu,1,jeo),u(:,:,mu,ieo,ibl), &
                    u(:,:,nu,jeo,iblv(ibl,mu)),temp1,lv(:,mu,ieo),1,nms,myid, &
                    nn,MRT)
        call setmmf(nu,lbd(ibl,nu),u(:,:,mu,ieo,iblv(ibl,nu)),ldiv(nu), &
                    ib(:,nu,1,jeo),temp1,u(:,:,mu,jeo,iblv(ibl,nu)),temp2, &
                    lv(:,nu,ieo),2,nms,myid,nn,MRT)
        call mmd(icount,temp2,u(:,:,nu,ieo,ibl),leaftot)

! Leaf 2.
        ibl0 = iblv(ibl,mu)
        call seteo(lbd(ibl,mu),jeo,ieo,ieo0,jeo0)
        call setmmf(nu,lbd(ibl0,nu),u(:,:,mu,ieo0,iblv(ibl0,nu)),ldiv(nu), &
                    ib(:,nu,1,jeo0),u(:,:,nu,ieo0,ibl0), &
                    u(:,:,mu,jeo0,iblv(ibl0,nu)),temp1,lv(:,nu,ieo0),1,nms, &
                    myid,nn,MRT)
        call setmmf(mu,lbd(ibl0,mu),u(:,:,nu,ieo0,iblv(ibl0,mu)),ldiv(mu), &
                    ib(:,mu,1,jeo0),temp1,u(:,:,nu,jeo0,iblv(ibl0,mu)),temp2, &
                    lv(:,mu,ieo0),2,nms,myid,nn,MRT)
        call mdm(icount,temp2,u(:,:,mu,ieo0,ibl0),temp1)
        if (lbd(ibl,mu)) then
         do i = 1,nvhalf
          temp2(:,lv(i,mu,jeo)) = temp1(:,i)
         enddo ! i
         if (ldiv(mu)) call slidematrix(temp2,ib(:,mu,1,ieo),mu,nms,myid,nn,MRT)
         leaftot = leaftot + temp2
        else
         leaftot = leaftot + temp1
        endif

! Leaf 3.
        ibl1 = iblv(ibl,nu)
        ibl0 = iblv(ibl1,mu)
        call seteo(lbd(ibl,nu),jeo,ieo,ieo1,jeo1)
        call seteo(lbd(ibl,mu),jeo1,ieo1,ieo0,jeo0)
        call mdm(icount,u(:,:,mu,ieo0,ibl0),u(:,:,nu,ieo0,ibl0),temp1)
        call setmmf(nu,lbd(ibl0,nu),u(:,:,mu,ieo0,iblv(ibl0,nu)),ldiv(nu), &
                    ib(:,nu,1,jeo0),temp1,u(:,:,mu,jeo0,iblv(ibl0,nu)),temp2, &
                    lv(:,nu,ieo0),1,nms,myid,nn,MRT)
        if (lbd(ibl,mu)) then
         do i = 1,nvhalf
          temp1(:,lv(i,mu,jeo1)) = temp2(:,i)
         enddo ! i
         if (ldiv(mu)) call slidematrix(temp1,ib(:,mu,1,ieo1),mu,nms,myid, &
                                        nn,MRT)
        else
         temp1 = temp2
        endif
        call mdm(icount,temp1,u(:,:,nu,ieo1,ibl1),temp2)
        if (lbd(ibl,nu)) then
         do i = 1,nvhalf
          temp1(:,lv(i,nu,jeo)) = temp2(:,i)
         enddo ! i
         if (ldiv(nu)) call slidematrix(temp1,ib(:,nu,1,ieo),nu,nms,myid,nn,MRT)
         leaftot = leaftot + temp1
        else
         leaftot = leaftot + temp2
        endif

! Leaf 4.
        ibl0 = iblv(ibl,nu)
        call seteo(lbd(ibl,nu),jeo,ieo,ieo0,jeo0)
        call mdm(icount,u(:,:,nu,ieo0,ibl0),u(:,:,mu,ieo0,ibl0),temp1)
        call setmmf(mu,lbd(ibl0,mu),u(:,:,nu,ieo0,iblv(ibl0,mu)),ldiv(mu), &
                    ib(:,mu,1,jeo0),temp1,u(:,:,nu,jeo0,iblv(ibl0,mu)),temp2, &
                    lv(:,mu,ieo0),1,nms,myid,nn,MRT)
        call setmmf(nu,lbd(ibl0,nu),u(:,:,mu,ieo0,iblv(ibl0,nu)),ldiv(nu), &
                    ib(:,nu,1,jeo0),temp2,u(:,:,mu,jeo0,iblv(ibl0,nu)),temp1, &
                    lv(:,nu,ieo0),2,nms,myid,nn,MRT)
        if (lbd(ibl,nu)) then
         do i = 1,nvhalf
          temp2(:,lv(i,nu,jeo)) = temp1(:,i)
         enddo ! i
         if (ldiv(nu)) call slidematrix(temp2,ib(:,nu,1,ieo),nu,nms,myid,nn,MRT)
         leaftot = leaftot + temp2
        else
         leaftot = leaftot + temp1
        endif

! Put in the tadpole and anisotropy coefficients of the clover term.
        leaftot = leaftot*coact(4,mu,nu)

! Build sigF from the four clover leaves.
        select case(iplane)
         case(1,6) ! the mu=1,nu=2 plane or the mu=3,nu=4 plane
          if (iplane==1) VorW(1:2)=(/1,4/)
          if (iplane==6) VorW(1:2)=(/5,8/)
          do j = 1,9
           jcc = 2*j
           icc = jcc - 1
           do i = 1,nvhalf
            sigF(icc,i,VorW(1),ieo,ibl) = &
            sigF(icc,i,VorW(1),ieo,ibl) - leaftot(jcc,i) - leaftot(kcc(2,j),i)
            sigF(jcc,i,VorW(1),ieo,ibl) = &
            sigF(jcc,i,VorW(1),ieo,ibl) + leaftot(icc,i) - leaftot(kcc(1,j),i)
            sigF(icc,i,VorW(2),ieo,ibl) = &
            sigF(icc,i,VorW(2),ieo,ibl) + leaftot(jcc,i) + leaftot(kcc(2,j),i)
            sigF(jcc,i,VorW(2),ieo,ibl) = &
            sigF(jcc,i,VorW(2),ieo,ibl) - leaftot(icc,i) + leaftot(kcc(1,j),i)
           enddo ! i
          enddo ! j
         case(2,5) ! the mu=1,nu=3 plane or the mu=2,nu=4 plane
          if (iplane==2) then
           VorW(1:2)=(/2,3/)
           sgnVW = -1.0_KR
          endif
          if (iplane==5) then
           VorW(1:2)=(/6,7/)
           sgnVW = 1.0_KR
          endif
          do j = 1,9
           jcc = 2*j
           icc = jcc - 1
           do i = 1,nvhalf
            sigF(icc,i,VorW(1),ieo,ibl) = sigF(icc,i,VorW(1),ieo,ibl) &
                            + sgnVW*( leaftot(icc,i) - leaftot(kcc(1,j),i) )
            sigF(jcc,i,VorW(1),ieo,ibl) = sigF(jcc,i,VorW(1),ieo,ibl) &
                            + sgnVW*( leaftot(jcc,i) + leaftot(kcc(2,j),i) )
            sigF(icc,i,VorW(2),ieo,ibl) = sigF(icc,i,VorW(2),ieo,ibl) &
                            - sgnVW*( leaftot(icc,i) - leaftot(kcc(1,j),i) )
            sigF(jcc,i,VorW(2),ieo,ibl) = sigF(jcc,i,VorW(2),ieo,ibl) &
                            - sgnVW*( leaftot(jcc,i) + leaftot(kcc(2,j),i) )
           enddo ! i
          enddo ! j
         case(3,4) ! the mu=1,nu=4 plane or the mu=2,nu=3 plane
          if (iplane==3) VorW(1:2)=(/6,7/)
          if (iplane==4) VorW(1:2)=(/2,3/)
          do j = 1,9
           jcc = 2*j
           icc = jcc - 1
           do i = 1,nvhalf
            sigF(icc,i,VorW(1),ieo,ibl) = &
            sigF(icc,i,VorW(1),ieo,ibl) - leaftot(jcc,i) - leaftot(kcc(2,j),i)
            sigF(jcc,i,VorW(1),ieo,ibl) = &
            sigF(jcc,i,VorW(1),ieo,ibl) + leaftot(icc,i) - leaftot(kcc(1,j),i)
            sigF(icc,i,VorW(2),ieo,ibl) = &
            sigF(icc,i,VorW(2),ieo,ibl) - leaftot(jcc,i) - leaftot(kcc(2,j),i)
            sigF(jcc,i,VorW(2),ieo,ibl) = &
            sigF(jcc,i,VorW(2),ieo,ibl) + leaftot(icc,i) - leaftot(kcc(1,j),i)
           enddo ! i
          enddo ! j
         case default
          open(unit=8,file="GAUGETOOLS.ERROR",action="write", &
               status="replace",form="formatted")
           write(unit=8,fmt=*) "subroutine makeclover: iplane =", iplane
          close(unit=8,status="keep")
          stop
        end select

       enddo ! nu
      enddo ! mu
     enddo ! ibl
    enddo ! ieo

! A final overall normalization of F_(mu,nu).
    sigF = 0.125_KR*sigF

 end subroutine makeclover

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 subroutine setmmf(mu,l1,b1,l,ibbit,a,b,c,lvbit,iflg,nms,myid,nn,MRT)
! Get links from a neighbouring process, and multiply two links to form part
! of a staple.

    integer(kind=KI), intent(in)                    :: mu, iflg, myid, MRT
    logical,          intent(in)                    :: l1, l
    real(kind=KR),    intent(in),    dimension(:,:) :: b1, a
    integer(kind=KI), intent(in),    dimension(:)   :: ibbit
    real(kind=KR),    intent(inout), dimension(:,:) :: b
    real(kind=KR),    intent(out),   dimension(:,:) :: c
    integer(kind=KI), intent(in),    dimension(:)   :: lvbit
    integer(kind=KI), intent(in),    dimension(:)   :: nms
    integer(kind=KI), intent(in),    dimension(:,:) :: nn

    integer(kind=KI) :: icount
    icount = nvhalf

!*If this is the first of the two blocks in this direction, then no links are 
! needed from a neighbouring process.
    if (l1) then
     if (iflg==1) then
      call mm(icount,a,b1,c)
     elseif (iflg==2) then
      call mmd(icount,a,b1,c)
     endif
!*If this is the second of the two blocks in this direction, then get the 
! second link from a neighbouring process if there are multiple processes 
! in the mu direction.
    else
     if (l) call setlink(b,ibbit,mu,nms,myid,nn,MRT)
     if (iflg==1) then
      call mmg(icount,a,b,c,lvbit)
     elseif (iflg==2) then
      call mmdg(icount,a,b,c,lvbit)
     endif
    endif

 end subroutine setmmf

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 subroutine seteo(l,ieo,jeo,ieo0,jeo0)
! Decide if the neighbouring site is even or odd.
! (This is nontrivial if there are next-nearest neighbour interactions.)

    logical,          intent(in)  :: l
    integer(kind=KI), intent(in)  :: ieo, jeo
    integer(kind=KI), intent(out) :: ieo0, jeo0

    if (l) then
     ieo0 = ieo
     jeo0 = jeo
    else
     ieo0 = jeo
     jeo0 = ieo
    endif

 end subroutine seteo

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 subroutine mm(n,y,v,w)
! Matrix multiplier: calculate w(i) = y(i)*v(i)
! INPUT:
!   y() is a set of 3x3 matrices 
!       (e.g. all links in mu direction on even sites),
!       expected size: y(18,n)
!   v() is a set of 3x3 matrices,
!       expected size: v(18,n)
!   n is the number of y matrices = number of v matrices,
! OUTPUT:
!   w(i) = y(i)*v(i) is the output.
!          expected size: w(18,n)

    integer(kind=KI), intent(in)                  :: n
    real(kind=KR),    intent(in),  dimension(:,:) :: y, v
    real(kind=KR),    intent(out), dimension(:,:) :: w

    integer(kind=KI) :: i
 
    do i = 1,n
     w(1 ,i) =   y(1 ,i) * v(1 ,i) - y(2 ,i) * v(2 ,i) &
               + y(3 ,i) * v(7 ,i) - y(4 ,i) * v(8 ,i) &
               + y(5 ,i) * v(13,i) - y(6 ,i) * v(14,i)
     w(2 ,i) =   y(1 ,i) * v(2 ,i) + y(2 ,i) * v(1 ,i) &
               + y(3 ,i) * v(8 ,i) + y(4 ,i) * v(7 ,i) &
               + y(5 ,i) * v(14,i) + y(6 ,i) * v(13,i)
     w(3 ,i) =   y(1 ,i) * v(3 ,i) - y(2 ,i) * v(4 ,i) &
               + y(3 ,i) * v(9 ,i) - y(4 ,i) * v(10,i) &
               + y(5 ,i) * v(15,i) - y(6 ,i) * v(16,i)
     w(4 ,i) =   y(1 ,i) * v(4 ,i) + y(2 ,i) * v(3 ,i) &
               + y(3 ,i) * v(10,i) + y(4 ,i) * v(9 ,i) &
               + y(5 ,i) * v(16,i) + y(6 ,i) * v(15,i)
     w(5 ,i) =   y(1 ,i) * v(5 ,i) - y(2 ,i) * v(6 ,i) &
               + y(3 ,i) * v(11,i) - y(4 ,i) * v(12,i) &
               + y(5 ,i) * v(17,i) - y(6 ,i) * v(18,i)
     w(6 ,i) =   y(1 ,i) * v(6 ,i) + y(2 ,i) * v(5 ,i) &
               + y(3 ,i) * v(12,i) + y(4 ,i) * v(11,i) &
               + y(5 ,i) * v(18,i) + y(6 ,i) * v(17,i)
     w(7 ,i) =   y(7 ,i) * v(1 ,i) - y(8 ,i) * v(2 ,i) &
               + y(9 ,i) * v(7 ,i) - y(10,i) * v(8 ,i) &
               + y(11,i) * v(13,i) - y(12,i) * v(14,i)
     w(8 ,i) =   y(7 ,i) * v(2 ,i) + y(8 ,i) * v(1 ,i) &
               + y(9 ,i) * v(8 ,i) + y(10,i) * v(7 ,i) &
               + y(11,i) * v(14,i) + y(12,i) * v(13,i)
     w(9 ,i) =   y(7 ,i) * v(3 ,i) - y(8 ,i) * v(4 ,i) &
               + y(9 ,i) * v(9 ,i) - y(10,i) * v(10,i) &
               + y(11,i) * v(15,i) - y(12,i) * v(16,i)
     w(10,i) =   y(7 ,i) * v(4 ,i) + y(8 ,i) * v(3 ,i) &
               + y(9 ,i) * v(10,i) + y(10,i) * v(9 ,i) &
               + y(11,i) * v(16,i) + y(12,i) * v(15,i)
     w(11,i) =   y(7 ,i) * v(5 ,i) - y(8 ,i) * v(6 ,i) &
               + y(9 ,i) * v(11,i) - y(10,i) * v(12,i) &
               + y(11,i) * v(17,i) - y(12,i) * v(18,i)
     w(12,i) =   y(7 ,i) * v(6 ,i) + y(8 ,i) * v(5 ,i) &
               + y(9 ,i) * v(12,i) + y(10,i) * v(11,i) &
               + y(11,i) * v(18,i) + y(12,i) * v(17,i)
     w(13,i) =   y(13,i) * v(1 ,i) - y(14,i) * v(2 ,i) &
               + y(15,i) * v(7 ,i) - y(16,i) * v(8 ,i) &
               + y(17,i) * v(13,i) - y(18,i) * v(14,i)
     w(14,i) =   y(13,i) * v(2 ,i) + y(14,i) * v(1 ,i) &
               + y(15,i) * v(8 ,i) + y(16,i) * v(7 ,i) &
               + y(17,i) * v(14,i) + y(18,i) * v(13,i)
     w(15,i) =   y(13,i) * v(3 ,i) - y(14,i) * v(4 ,i) &
               + y(15,i) * v(9 ,i) - y(16,i) * v(10,i) &
               + y(17,i) * v(15,i) - y(18,i) * v(16,i)
     w(16,i) =   y(13,i) * v(4 ,i) + y(14,i) * v(3 ,i) &
               + y(15,i) * v(10,i) + y(16,i) * v(9 ,i) &
               + y(17,i) * v(16,i) + y(18,i) * v(15,i)
     w(17,i) =   y(13,i) * v(5 ,i) - y(14,i) * v(6 ,i) &
               + y(15,i) * v(11,i) - y(16,i) * v(12,i) &
               + y(17,i) * v(17,i) - y(18,i) * v(18,i)
     w(18,i) =   y(13,i) * v(6 ,i) + y(14,i) * v(5 ,i) &
               + y(15,i) * v(12,i) + y(16,i) * v(11,i) &
               + y(17,i) * v(18,i) + y(18,i) * v(17,i)
    enddo ! i

 end subroutine mm

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 subroutine mmg(n,y,v,w,j)
! Matrix multiplier: calculate w(i) = y(i)*v(j(i))
! INPUT:
!   y() is a set of 3x3 matrices 
!       (e.g. all links in mu direction on even sites),
!       expected size: y(18,n)
!   v() is a set of 3x3 matrices,
!       expected size: v(18,m)
!   m is number of v matrices (does not appear in this subroutine),
!   j(i) is typically the neighbouring site (to i) in the appropriate
!        direction so w(i) will represent part of a staple,
!        expected size: j(n)
!   n is the number of entries in j(),
! OUTPUT:
!   w(i) = y(i)*v(j(i)) is the output.
!        expected size: w(18,n)

! QCDMPI&QCDimMPI dimension y(18,m) which seems incorrect (but harmless)
    integer(kind=KI), intent(in)                  :: n
    integer(kind=KI), intent(in),  dimension(:)   :: j
    real(kind=KR),    intent(in),  dimension(:,:) :: y, v
    real(kind=KR),    intent(out), dimension(:,:) :: w

    integer(kind=KI) :: i, k
 
    do i = 1,n
     k = j(i)
     w(1 ,i) =   y(1 ,i) * v(1 ,k) - y(2 ,i) * v(2 ,k) &
               + y(3 ,i) * v(7 ,k) - y(4 ,i) * v(8 ,k) &
               + y(5 ,i) * v(13,k) - y(6 ,i) * v(14,k)  
     w(2 ,i) =   y(1 ,i) * v(2 ,k) + y(2 ,i) * v(1 ,k) &
               + y(3 ,i) * v(8 ,k) + y(4 ,i) * v(7 ,k) &
               + y(5 ,i) * v(14,k) + y(6 ,i) * v(13,k)  
     w(3 ,i) =   y(1 ,i) * v(3 ,k) - y(2 ,i) * v(4 ,k) &
               + y(3 ,i) * v(9 ,k) - y(4 ,i) * v(10,k) &
               + y(5 ,i) * v(15,k) - y(6 ,i) * v(16,k)  
     w(4 ,i) =   y(1 ,i) * v(4 ,k) + y(2 ,i) * v(3 ,k) &
               + y(3 ,i) * v(10,k) + y(4 ,i) * v(9 ,k) &
               + y(5 ,i) * v(16,k) + y(6 ,i) * v(15,k)  
     w(5 ,i) =   y(1 ,i) * v(5 ,k) - y(2 ,i) * v(6 ,k) &
               + y(3 ,i) * v(11,k) - y(4 ,i) * v(12,k) &
               + y(5 ,i) * v(17,k) - y(6 ,i) * v(18,k)  
     w(6 ,i) =   y(1 ,i) * v(6 ,k) + y(2 ,i) * v(5 ,k) &
               + y(3 ,i) * v(12,k) + y(4 ,i) * v(11,k) &
               + y(5 ,i) * v(18,k) + y(6 ,i) * v(17,k)  
     w(7 ,i) =   y(7 ,i) * v(1 ,k) - y(8 ,i) * v(2 ,k) &
               + y(9 ,i) * v(7 ,k) - y(10,i) * v(8 ,k) &
               + y(11,i) * v(13,k) - y(12,i) * v(14,k)  
     w(8 ,i) =   y(7 ,i) * v(2 ,k) + y(8 ,i) * v(1 ,k) &
               + y(9 ,i) * v(8 ,k) + y(10,i) * v(7 ,k) &
               + y(11,i) * v(14,k) + y(12,i) * v(13,k)  
     w(9 ,i) =   y(7 ,i) * v(3 ,k) - y(8 ,i) * v(4 ,k) &
               + y(9 ,i) * v(9 ,k) - y(10,i) * v(10,k) &
               + y(11,i) * v(15,k) - y(12,i) * v(16,k)  
     w(10,i) =   y(7 ,i) * v(4 ,k) + y(8 ,i) * v(3 ,k) &
               + y(9 ,i) * v(10,k) + y(10,i) * v(9 ,k) &
               + y(11,i) * v(16,k) + y(12,i) * v(15,k)  
     w(11,i) =   y(7 ,i) * v(5 ,k) - y(8 ,i) * v(6 ,k) &
               + y(9 ,i) * v(11,k) - y(10,i) * v(12,k) &
               + y(11,i) * v(17,k) - y(12,i) * v(18,k)  
     w(12,i) =   y(7 ,i) * v(6 ,k) + y(8 ,i) * v(5 ,k) &
               + y(9 ,i) * v(12,k) + y(10,i) * v(11,k) &
               + y(11,i) * v(18,k) + y(12,i) * v(17,k)  
     w(13,i) =   y(13,i) * v(1 ,k) - y(14,i) * v(2 ,k) &
               + y(15,i) * v(7 ,k) - y(16,i) * v(8 ,k) &
               + y(17,i) * v(13,k) - y(18,i) * v(14,k)  
     w(14,i) =   y(13,i) * v(2 ,k) + y(14,i) * v(1 ,k) &
               + y(15,i) * v(8 ,k) + y(16,i) * v(7 ,k) &
               + y(17,i) * v(14,k) + y(18,i) * v(13,k)  
     w(15,i) =   y(13,i) * v(3 ,k) - y(14,i) * v(4 ,k) &
               + y(15,i) * v(9 ,k) - y(16,i) * v(10,k) &
               + y(17,i) * v(15,k) - y(18,i) * v(16,k)  
     w(16,i) =   y(13,i) * v(4 ,k) + y(14,i) * v(3 ,k) &
               + y(15,i) * v(10,k) + y(16,i) * v(9 ,k) &
               + y(17,i) * v(16,k) + y(18,i) * v(15,k)  
     w(17,i) =   y(13,i) * v(5 ,k) - y(14,i) * v(6 ,k) &
               + y(15,i) * v(11,k) - y(16,i) * v(12,k) &
               + y(17,i) * v(17,k) - y(18,i) * v(18,k)  
     w(18,i) =   y(13,i) * v(6 ,k) + y(14,i) * v(5 ,k) &
               + y(15,i) * v(12,k) + y(16,i) * v(11,k) &
               + y(17,i) * v(18,k) + y(18,i) * v(17,k)  
    enddo ! i

 end subroutine mmg

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 subroutine mmd(n,y,v,w)
! Matrix multiplier: calculate w(i) = y(i)*v^dagger(i)
! INPUT:
!   y() is a set of 3x3 matrices 
!       (e.g. all links in mu direction on even sites),
!       expected size: y(18,n)
!   v() is a set of 3x3 matrices,
!       expected size: v(18,n)
!   n is the number of y matrices = number of v matrices,
! OUTPUT:
!   w(i) = y(i)*v^dagger(i) is the output.
!       expected size: w(18,n)

    integer(kind=KI), intent(in)                  :: n
    real(kind=KR),    intent(in),  dimension(:,:) :: y, v
    real(kind=KR),    intent(out), dimension(:,:) :: w

    integer(kind=KI) :: i
 
    do i = 1,n
     w(1 ,i) =   y(1 ,i) * v(1 ,i) + y(2 ,i) * v(2 ,i) &
               + y(3 ,i) * v(3 ,i) + y(4 ,i) * v(4 ,i) &
               + y(5 ,i) * v(5 ,i) + y(6 ,i) * v(6 ,i) 
     w(2 ,i) = - y(1 ,i) * v(2 ,i) + y(2 ,i) * v(1 ,i) &
               - y(3 ,i) * v(4 ,i) + y(4 ,i) * v(3 ,i) &
               - y(5 ,i) * v(6 ,i) + y(6 ,i) * v(5 ,i) 
     w(3 ,i) =   y(1 ,i) * v(7 ,i) + y(2 ,i) * v(8 ,i) &
               + y(3 ,i) * v(9 ,i) + y(4 ,i) * v(10,i) &
               + y(5 ,i) * v(11,i) + y(6 ,i) * v(12,i) 
     w(4 ,i) = - y(1 ,i) * v(8 ,i) + y(2 ,i) * v(7 ,i) &
               - y(3 ,i) * v(10,i) + y(4 ,i) * v(9 ,i) &
               - y(5 ,i) * v(12,i) + y(6 ,i) * v(11,i) 
     w(5 ,i) =   y(1 ,i) * v(13,i) + y(2 ,i) * v(14,i) &
               + y(3 ,i) * v(15,i) + y(4 ,i) * v(16,i) &
               + y(5 ,i) * v(17,i) + y(6 ,i) * v(18,i) 
     w(6 ,i) = - y(1 ,i) * v(14,i) + y(2 ,i) * v(13,i) &
               - y(3 ,i) * v(16,i) + y(4 ,i) * v(15,i) &
               - y(5 ,i) * v(18,i) + y(6 ,i) * v(17,i) 
     w(7 ,i) =   y(7 ,i) * v(1 ,i) + y(8 ,i) * v(2 ,i) &
               + y(9 ,i) * v(3 ,i) + y(10,i) * v(4 ,i) &
               + y(11,i) * v(5 ,i) + y(12,i) * v(6 ,i) 
     w(8 ,i) = - y(7 ,i) * v(2 ,i) + y(8 ,i) * v(1 ,i) &
               - y(9 ,i) * v(4 ,i) + y(10,i) * v(3 ,i) &
               - y(11,i) * v(6 ,i) + y(12,i) * v(5 ,i) 
     w(9 ,i) =   y(7 ,i) * v(7 ,i) + y(8 ,i) * v(8 ,i) &
               + y(9 ,i) * v(9 ,i) + y(10,i) * v(10,i) &
               + y(11,i) * v(11,i) + y(12,i) * v(12,i) 
     w(10,i) = - y(7 ,i) * v(8 ,i) + y(8 ,i) * v(7 ,i) &
               - y(9 ,i) * v(10,i) + y(10,i) * v(9 ,i) &
               - y(11,i) * v(12,i) + y(12,i) * v(11,i) 
     w(11,i) =   y(7 ,i) * v(13,i) + y(8 ,i) * v(14,i) &
               + y(9 ,i) * v(15,i) + y(10,i) * v(16,i) &
               + y(11,i) * v(17,i) + y(12,i) * v(18,i) 
     w(12,i) = - y(7 ,i) * v(14,i) + y(8 ,i) * v(13,i) &
               - y(9 ,i) * v(16,i) + y(10,i) * v(15,i) &
               - y(11,i) * v(18,i) + y(12,i) * v(17,i) 
     w(13,i) =   y(13,i) * v(1 ,i) + y(14,i) * v(2 ,i) &
               + y(15,i) * v(3 ,i) + y(16,i) * v(4 ,i) &
               + y(17,i) * v(5 ,i) + y(18,i) * v(6 ,i) 
     w(14,i) = - y(13,i) * v(2 ,i) + y(14,i) * v(1 ,i) &
               - y(15,i) * v(4 ,i) + y(16,i) * v(3 ,i) &
               - y(17,i) * v(6 ,i) + y(18,i) * v(5 ,i) 
     w(15,i) =   y(13,i) * v(7 ,i) + y(14,i) * v(8 ,i) &
               + y(15,i) * v(9 ,i) + y(16,i) * v(10,i) &
               + y(17,i) * v(11,i) + y(18,i) * v(12,i) 
     w(16,i) = - y(13,i) * v(8 ,i) + y(14,i) * v(7 ,i) &
               - y(15,i) * v(10,i) + y(16,i) * v(9 ,i) &
               - y(17,i) * v(12,i) + y(18,i) * v(11,i) 
     w(17,i) =   y(13,i) * v(13,i) + y(14,i) * v(14,i) &
               + y(15,i) * v(15,i) + y(16,i) * v(16,i) &
               + y(17,i) * v(17,i) + y(18,i) * v(18,i) 
     w(18,i) = - y(13,i) * v(14,i) + y(14,i) * v(13,i) &
               - y(15,i) * v(16,i) + y(16,i) * v(15,i) &
               - y(17,i) * v(18,i) + y(18,i) * v(17,i) 
    enddo ! i

 end subroutine mmd

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 subroutine mmdg(n,y,v,w,j)
! Matrix multiplier: calculate w(i) = y(i)*v^dagger(j(i))
! INPUT:
!   y() is a set of 3x3 matrices 
!       (e.g. all links in mu direction on even sites),
!       expected size: y(18,n)
!   v() is a set of 3x3 matrices,
!       expected size: v(18,m)
!   m is the number of v matrices (not needed in this subroutine),
!   j(i) is typically the neighbouring site (to i) in the appropriate
!        direction so w(i) will represent part of a staple,
!        expected size: j(n)
!   n is the number of y matrices = number of entries in j(),
! OUTPUT:
!   w(i) = y(i)*v^dagger(j(i)) is the output.
!        expected size: w(18,n)

    integer(kind=KI), intent(in)                  :: n
    integer(kind=KI), intent(in),  dimension(:)   :: j
    real(kind=KR),    intent(in),  dimension(:,:) :: y, v
    real(kind=KR),    intent(out), dimension(:,:) :: w

    integer(kind=KI) :: i, k
 
    do i = 1,n
     k = j(i)
     w(1 ,i) =   y(1 ,i) * v(1 ,k) + y(2 ,i) * v(2 ,k) &
               + y(3 ,i) * v(3 ,k) + y(4 ,i) * v(4 ,k) &
               + y(5 ,i) * v(5 ,k) + y(6 ,i) * v(6 ,k) 
     w(2 ,i) = - y(1 ,i) * v(2 ,k) + y(2 ,i) * v(1 ,k) &
               - y(3 ,i) * v(4 ,k) + y(4 ,i) * v(3 ,k) &
               - y(5 ,i) * v(6 ,k) + y(6 ,i) * v(5 ,k) 
     w(3 ,i) =   y(1 ,i) * v(7 ,k) + y(2 ,i) * v(8 ,k) &
               + y(3 ,i) * v(9 ,k) + y(4 ,i) * v(10,k) &
               + y(5 ,i) * v(11,k) + y(6 ,i) * v(12,k) 
     w(4 ,i) = - y(1 ,i) * v(8 ,k) + y(2 ,i) * v(7 ,k) &
               - y(3 ,i) * v(10,k) + y(4 ,i) * v(9 ,k) &
               - y(5 ,i) * v(12,k) + y(6 ,i) * v(11,k) 
     w(5 ,i) =   y(1 ,i) * v(13,k) + y(2 ,i) * v(14,k) &
               + y(3 ,i) * v(15,k) + y(4 ,i) * v(16,k) &
               + y(5 ,i) * v(17,k) + y(6 ,i) * v(18,k) 
     w(6 ,i) = - y(1 ,i) * v(14,k) + y(2 ,i) * v(13,k) &
               - y(3 ,i) * v(16,k) + y(4 ,i) * v(15,k) &
               - y(5 ,i) * v(18,k) + y(6 ,i) * v(17,k) 
     w(7 ,i) =   y(7 ,i) * v(1 ,k) + y(8 ,i) * v(2 ,k) &
               + y(9 ,i) * v(3 ,k) + y(10,i) * v(4 ,k) &
               + y(11,i) * v(5 ,k) + y(12,i) * v(6 ,k) 
     w(8 ,i) = - y(7 ,i) * v(2 ,k) + y(8 ,i) * v(1 ,k) &
               - y(9 ,i) * v(4 ,k) + y(10,i) * v(3 ,k) &
               - y(11,i) * v(6 ,k) + y(12,i) * v(5 ,k) 
     w(9 ,i) =   y(7 ,i) * v(7 ,k) + y(8 ,i) * v(8 ,k) &
               + y(9 ,i) * v(9 ,k) + y(10,i) * v(10,k) &
               + y(11,i) * v(11,k) + y(12,i) * v(12,k) 
     w(10,i) = - y(7 ,i) * v(8 ,k) + y(8 ,i) * v(7 ,k) &
               - y(9 ,i) * v(10,k) + y(10,i) * v(9 ,k) &
               - y(11,i) * v(12,k) + y(12,i) * v(11,k) 
     w(11,i) =   y(7 ,i) * v(13,k) + y(8 ,i) * v(14,k) &
               + y(9 ,i) * v(15,k) + y(10,i) * v(16,k) &
               + y(11,i) * v(17,k) + y(12,i) * v(18,k) 
     w(12,i) = - y(7 ,i) * v(14,k) + y(8 ,i) * v(13,k) &
               - y(9 ,i) * v(16,k) + y(10,i) * v(15,k) &
               - y(11,i) * v(18,k) + y(12,i) * v(17,k) 
     w(13,i) =   y(13,i) * v(1 ,k) + y(14,i) * v(2 ,k) &
               + y(15,i) * v(3 ,k) + y(16,i) * v(4 ,k) &
               + y(17,i) * v(5 ,k) + y(18,i) * v(6 ,k) 
     w(14,i) = - y(13,i) * v(2 ,k) + y(14,i) * v(1 ,k) &
               - y(15,i) * v(4 ,k) + y(16,i) * v(3 ,k) &
               - y(17,i) * v(6 ,k) + y(18,i) * v(5 ,k) 
     w(15,i) =   y(13,i) * v(7 ,k) + y(14,i) * v(8 ,k) &
               + y(15,i) * v(9 ,k) + y(16,i) * v(10,k) &
               + y(17,i) * v(11,k) + y(18,i) * v(12,k) 
     w(16,i) = - y(13,i) * v(8 ,k) + y(14,i) * v(7 ,k) &
               - y(15,i) * v(10,k) + y(16,i) * v(9 ,k) &
               - y(17,i) * v(12,k) + y(18,i) * v(11,k) 
     w(17,i) =   y(13,i) * v(13,k) + y(14,i) * v(14,k) &
               + y(15,i) * v(15,k) + y(16,i) * v(16,k) &
               + y(17,i) * v(17,k) + y(18,i) * v(18,k) 
     w(18,i) = - y(13,i) * v(14,k) + y(14,i) * v(13,k) &
               - y(15,i) * v(16,k) + y(16,i) * v(15,k) &
               - y(17,i) * v(18,k) + y(18,i) * v(17,k) 
    enddo ! i

 end subroutine mmdg

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 subroutine mdm(n,y,v,w)
! Matrix multiplier: calculate w(i) = y^dagger(i)*v(i)
! INPUT:
!   y() is a set of 3x3 matrices 
!       (e.g. all links in nu direction on even sites),
!       expected size: y(18,n)
!   v() is a set of 3x3 matrices,
!       (e.g. all links in mu direction on even sites),
!       expected size: v(18,n)
!   n is the number of y matrices = number of v matrices,
! OUTPUT:
!   w(i) = y^dagger(i)*v(i) is the output.
!       expected size: w(18,n)

    integer(kind=KI), intent(in)                  :: n
    real(kind=KR),    intent(in),  dimension(:,:) :: y, v
    real(kind=KR),    intent(out), dimension(:,:) :: w

    integer(kind=KI) :: i
 
    do i = 1,n
     w(1 ,i) =   y(1 ,i) * v(1 ,i) + y(2 ,i) * v(2 ,i) &
               + y(7 ,i) * v(7 ,i) + y(8 ,i) * v(8 ,i) &
               + y(13,i) * v(13,i) + y(14,i) * v(14,i) 
     w(2 ,i) =   y(1 ,i) * v(2 ,i) - y(2 ,i) * v(1 ,i) &
               + y(7 ,i) * v(8 ,i) - y(8 ,i) * v(7 ,i) &
               + y(13,i) * v(14,i) - y(14,i) * v(13,i) 
     w(3 ,i) =   y(1 ,i) * v(3 ,i) + y(2 ,i) * v(4 ,i) &
               + y(7 ,i) * v(9 ,i) + y(8 ,i) * v(10,i) &
               + y(13,i) * v(15,i) + y(14,i) * v(16,i) 
     w(4 ,i) =   y(1 ,i) * v(4 ,i) - y(2 ,i) * v(3 ,i) &
               + y(7 ,i) * v(10,i) - y(8 ,i) * v(9 ,i) &
               + y(13,i) * v(16,i) - y(14,i) * v(15,i) 
     w(5 ,i) =   y(1 ,i) * v(5 ,i) + y(2 ,i) * v(6 ,i) &
               + y(7 ,i) * v(11,i) + y(8 ,i) * v(12,i) &
               + y(13,i) * v(17,i) + y(14,i) * v(18,i) 
     w(6 ,i) =   y(1 ,i) * v(6 ,i) - y(2 ,i) * v(5 ,i) &
               + y(7 ,i) * v(12,i) - y(8 ,i) * v(11,i) &
               + y(13,i) * v(18,i) - y(14,i) * v(17,i) 
     w(7 ,i) =   y(3 ,i) * v(1 ,i) + y(4 ,i) * v(2 ,i) &
               + y(9 ,i) * v(7 ,i) + y(10,i) * v(8 ,i) &
               + y(15,i) * v(13,i) + y(16,i) * v(14,i) 
     w(8 ,i) =   y(3 ,i) * v(2 ,i) - y(4 ,i) * v(1 ,i) &
               + y(9 ,i) * v(8 ,i) - y(10,i) * v(7 ,i) &
               + y(15,i) * v(14,i) - y(16,i) * v(13,i) 
     w(9 ,i) =   y(3 ,i) * v(3 ,i) + y(4 ,i) * v(4 ,i) &
               + y(9 ,i) * v(9 ,i) + y(10,i) * v(10,i) &
               + y(15,i) * v(15,i) + y(16,i) * v(16,i) 
     w(10,i) =   y(3 ,i) * v(4 ,i) - y(4 ,i) * v(3 ,i) &
               + y(9 ,i) * v(10,i) - y(10,i) * v(9 ,i) &
               + y(15,i) * v(16,i) - y(16,i) * v(15,i) 
     w(11,i) =   y(3 ,i) * v(5 ,i) + y(4 ,i) * v(6 ,i) &
               + y(9 ,i) * v(11,i) + y(10,i) * v(12,i) &
               + y(15,i) * v(17,i) + y(16,i) * v(18,i) 
     w(12,i) =   y(3 ,i) * v(6 ,i) - y(4 ,i) * v(5 ,i) &
               + y(9 ,i) * v(12,i) - y(10,i) * v(11,i) &
               + y(15,i) * v(18,i) - y(16,i) * v(17,i) 
     w(13,i) =   y(5 ,i) * v(1 ,i) + y(6 ,i) * v(2 ,i) &
               + y(11,i) * v(7 ,i) + y(12,i) * v(8 ,i) &
               + y(17,i) * v(13,i) + y(18,i) * v(14,i) 
     w(14,i) =   y(5 ,i) * v(2 ,i) - y(6 ,i) * v(1 ,i) &
               + y(11,i) * v(8 ,i) - y(12,i) * v(7 ,i) &
               + y(17,i) * v(14,i) - y(18,i) * v(13,i) 
     w(15,i) =   y(5 ,i) * v(3 ,i) + y(6 ,i) * v(4 ,i) &
               + y(11,i) * v(9 ,i) + y(12,i) * v(10,i) &
               + y(17,i) * v(15,i) + y(18,i) * v(16,i) 
     w(16,i) =   y(5 ,i) * v(4 ,i) - y(6 ,i) * v(3 ,i) &
               + y(11,i) * v(10,i) - y(12,i) * v(9 ,i) &
               + y(17,i) * v(16,i) - y(18,i) * v(15,i) 
     w(17,i) =   y(5 ,i) * v(5 ,i) + y(6 ,i) * v(6 ,i) &
               + y(11,i) * v(11,i) + y(12,i) * v(12,i) &
               + y(17,i) * v(17,i) + y(18,i) * v(18,i) 
     w(18,i) =   y(5 ,i) * v(6 ,i) - y(6 ,i) * v(5 ,i) &
               + y(11,i) * v(12,i) - y(12,i) * v(11,i) &
               + y(17,i) * v(18,i) - y(18,i) * v(17,i)
    enddo ! i

 end subroutine mdm

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 subroutine setlink(amatrix,ibbit,mu,nms,myid,nn,MRT)
! Send/receive links to/from a neighbouring process which are required
! for completion of staples (e.g. staples for all links in mu direction
! on even sites).
! INPUT:
!   amatrix() is the matrix containing the links to be transferred.
!             expected size: amatrix(18,ntotal)
!   ibbit() is a numbering of boundary sites for this block of the sublattice
!           taken from EITHER the +mu edge OR the -mu edge of the
!           neighbouring sublattice.
!           expected size: ibbit(nbmax)
!   mu is the direction of the links to be transferred.
!   MRT is MPIREALTYPE.
! OUTPUT:
!   amatrix() is updated.

    real(kind=KR),    intent(inout), dimension(:,:) :: amatrix
    integer(kind=KI), intent(in),    dimension(:)   :: ibbit, nms
    integer(kind=KI), intent(in)                    :: mu, myid, MRT
    integer(kind=KI), intent(in),    dimension(:,:) :: nn

    integer(kind=KI)                             :: i, ierr
    integer(kind=KI), dimension(MPI_STATUS_SIZE) :: istatus
    integer(kind=KI), dimension(4)               :: npct

! Locate links needed by neighbouring process and put them into the buffer.
    do i = 1,nms(mu)
     amatrix(:,nvhalf+i) = amatrix(:,ibbit(i))
    enddo ! i
    npct = 18*nms
! Send buffer links to neighbour and receive buffer links from other neighbour.
    call MPI_SENDRECV_REPLACE(amatrix(1,nvhalf+1),npct(mu), &
                              MRT,nn(mu,2),myid,nn(mu,1),   &
                              nn(mu,1),MPI_COMM_WORLD,istatus,ierr)

 end subroutine setlink

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 
 subroutine slidematrix(w,ibbit,mu,nms,myid,nn,MRT)
! Move the matrix w from each process to a neighbouring process.
! INPUT:
!   w() is the set of matrices to be transferred.
!   ibbit() is a numbering of boundary sites for this block of the sublattice
!           in EITHER the +mu direction OR the -mu direction.
!   mu is the direction of the links to be transferred.
!   MRT is MPIREALTYPE.
! OUTPUT:
!   w() is updated.
 
    real(kind=KR),    intent(inout), dimension(:,:) :: w
    integer(kind=KI), intent(in)                    :: mu, myid, MRT
    integer(kind=KI), intent(in),    dimension(:)   :: ibbit
    integer(kind=KI), intent(in),    dimension(:)   :: nms
    integer(kind=KI), intent(in),    dimension(:,:) :: nn
 
    integer(kind=KI)                             :: i, ierr
    integer(kind=KI), dimension(MPI_STATUS_SIZE) :: istatus
    integer(kind=KI), dimension(4)               :: npct
 
! Send buffer links to neighbour and receive buffer links from other neighbour.
    npct = 18*nms
    call MPI_SENDRECV_REPLACE(w(1,nvhalf+1),npct(mu),MRT,nn(mu,1),myid, &
                              nn(mu,2),nn(mu,2),MPI_COMM_WORLD,istatus,ierr)
 
! Move links from the buffer to their true position on the lattice.
    do i = 1,nms(mu)
     w(:,ibbit(i)) = w(:,nvhalf+i)
    enddo ! i
 
 end subroutine slidematrix
 
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 
 subroutine trace(mu,ieo,ibl,u,c,plaq)
! Compute ReTr(u * c^dagger) and add it to the existing value of plaq.
! INPUT:
!   mu is the spacetime direction of the links of interest.
!   ieo is the part of the checkerboard containing the links of interest.
!   ibl is the block of the sublattice containing the links of interest.
!   c() is a sum over staples.  expected size: c(18,nvhalf)
!   plaq is a partial result for ReTr(u * c^dagger).
! OUTPUT:
!   plaq is updated.
 
    integer(kind=KI), intent(in)                       :: mu, ieo, ibl
    real(kind=KR),    intent(in), dimension(:,:,:,:,:) :: u
    real(kind=KR),    intent(in), dimension(:,:)       :: c
    real(kind=KR),    intent(inout)                    :: plaq

    integer(kind=KI) :: i, j
    real(kind=KR2)   :: plaqbit
 
    plaqbit = 0.0_KR2
    do i = 1,nvhalf
     do j = 1,18
      plaqbit = plaqbit + real( u(j,i,mu,ieo,ibl)*c(j,i) ,KR2)
     enddo ! j
    enddo ! i
    plaq = plaq + real(plaqbit,KR)
 
 end subroutine trace

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 subroutine GFwrite(u,filename,gaction,beta,alat,tad,icfgsave,myid)
! Writes truncated gauge field configuration to disk.
! INPUT:
!   filename is the file to be created.
!   gaction(i) = 1 for nearest neighbour action in the i'th spacetime direction.
!   gaction(i) = 2 for next-nearest neighbour action in the i'th direction.
!   beta is the plaquette coupling.
!   alat() is the ratios of lattice spacings in each spacetime direction.
!   tad() is the tadpole factors in each spacetime direction.
!   icfgsave is the integer which will be used to identify this configuration.
 
    integer(kind=KI), intent(in)                          :: icfgsave, myid
    real(kind=KR),    intent(inout), dimension(:,:,:,:,:) :: u
    character(len=*), intent(in)                          :: filename
    real(kind=KR),    intent(in)                          :: beta
    real(kind=KR),    intent(in),    dimension(:)         :: alat, tad
    integer(kind=KI), intent(in),    dimension(:)         :: gaction
 
    integer(kind=KI) :: icri, isite, ibl, mu
 
    call fixsu3(u)
 
    open(unit=9,file=trim(filename),action="write",status="new", &
         form="unformatted")
     do icri = 1,12
      do isite = 1,nvhalf
       do ibl = 1,16
        write(unit=9) (u(icri,isite,mu,1,ibl),u(icri,isite,mu,2,ibl),mu=1,4)
       enddo ! ibl
      enddo ! isite
     enddo ! icri
     write(unit=9) gaction, beta, alat, tad, icfgsave, myid
    close(unit=9,status="keep")
 
 end subroutine GFwrite

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 subroutine GFread(u,filename,gaction,beta,alat,tad,icfgsave,myidin)
! Read in truncated gauge field configuration, and complete the entries.
! INPUT:
!   filename is the file to be read in.
! OUTPUT:
!   gaction(i) = 1 for nearest neighbour action in the i'th spacetime direction.
!   gaction(i) = 2 for next-nearest neighbour action in the i'th direction.
!   beta is the plaquette coupling.
!   alat(ndim) is the ratios of lattice spacings in each spacetime direction.
!   tad(ndim) is the tadpole factors in each spacetime direction.
!   icfgsave is the integer which has been used to identify this configuration.
!   myidin is "myid" for the process that created this partial configuration.
 
    real(kind=KR),    intent(inout), dimension(:,:,:,:,:) :: u
    character(len=*), intent(in)                          :: filename
    real(kind=KR),    intent(out)                         :: beta
    real(kind=KR),    intent(out),   dimension(:)         :: alat, tad
    integer(kind=KI), intent(out),   dimension(:)         :: gaction
    integer(kind=KI), intent(out)                         :: icfgsave, myidin
 
    integer(kind=KI) :: isite, icri, ibl, mu
 
    open(unit=9,file=trim(filename),action="read",status="old", &
         position="rewind",form="unformatted")
     do icri = 1,12
      do isite = 1,nvhalf
       do ibl = 1,16
        read(unit=9) (u(icri,isite,mu,1,ibl),u(icri,isite,mu,2,ibl),mu=1,4)
       enddo ! ibl
      enddo ! isite
     enddo ! icri
     read(unit=9) gaction, beta, alat, tad, icfgsave, myidin
    close(unit=9,status="keep")
 
    call fixsu3(u)
 
 end subroutine GFread

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 subroutine fixsu3(u)
! Unitarizes a lattice array of 3x3 complex matrices.
 
    real(kind=KR),    intent(inout), dimension(:,:,:,:,:) :: u

    integer(kind=KI) :: i, ieo, mu, isite, ibl
    real(kind=KR)    :: c1r,c1i
 
    do ieo = 1,2
     do mu = 1,4
      do isite = 1,nvhalf
       do ibl = 1,16
 
! Normalize the first row.
        c1r = sqrt(u(1,isite,mu,ieo,ibl)**2 &
                 + u(2,isite,mu,ieo,ibl)**2 &
                 + u(3,isite,mu,ieo,ibl)**2 &
                 + u(4,isite,mu,ieo,ibl)**2 &
                 + u(5,isite,mu,ieo,ibl)**2 &
                 + u(6,isite,mu,ieo,ibl)**2)
        do i = 1,6
         u(i,isite,mu,ieo,ibl) = u(i,isite,mu,ieo,ibl)/c1r
        enddo ! i
! Create the second row: row2 - (row2 dot row1)row1.
        c1r = u(7 ,isite,mu,ieo,ibl)*u(1,isite,mu,ieo,ibl) &
            + u(8 ,isite,mu,ieo,ibl)*u(2,isite,mu,ieo,ibl) &
            + u(9 ,isite,mu,ieo,ibl)*u(3,isite,mu,ieo,ibl) &
            + u(10,isite,mu,ieo,ibl)*u(4,isite,mu,ieo,ibl) &
            + u(11,isite,mu,ieo,ibl)*u(5,isite,mu,ieo,ibl) &
            + u(12,isite,mu,ieo,ibl)*u(6,isite,mu,ieo,ibl)
        c1i = u(8 ,isite,mu,ieo,ibl)*u(1,isite,mu,ieo,ibl) &
            - u(7 ,isite,mu,ieo,ibl)*u(2,isite,mu,ieo,ibl) &
            + u(10,isite,mu,ieo,ibl)*u(3,isite,mu,ieo,ibl) &
            - u(9 ,isite,mu,ieo,ibl)*u(4,isite,mu,ieo,ibl) &
            + u(12,isite,mu,ieo,ibl)*u(5,isite,mu,ieo,ibl) &
            - u(11,isite,mu,ieo,ibl)*u(6,isite,mu,ieo,ibl)
        u(7 ,isite,mu,ieo,ibl) = u(7 ,isite,mu,ieo,ibl) &
               - c1r*u(1,isite,mu,ieo,ibl) + c1i*u(2,isite,mu,ieo,ibl)
        u(8 ,isite,mu,ieo,ibl) = u(8 ,isite,mu,ieo,ibl) &
               - c1r*u(2,isite,mu,ieo,ibl) - c1i*u(1,isite,mu,ieo,ibl)
        u(9 ,isite,mu,ieo,ibl) = u(9 ,isite,mu,ieo,ibl) &
               - c1r*u(3,isite,mu,ieo,ibl) + c1i*u(4,isite,mu,ieo,ibl)
        u(10,isite,mu,ieo,ibl) = u(10,isite,mu,ieo,ibl) &
               - c1r*u(4,isite,mu,ieo,ibl) - c1i*u(3,isite,mu,ieo,ibl)
        u(11,isite,mu,ieo,ibl) = u(11,isite,mu,ieo,ibl) &
               - c1r*u(5,isite,mu,ieo,ibl) + c1i*u(6,isite,mu,ieo,ibl)
        u(12,isite,mu,ieo,ibl) = u(12,isite,mu,ieo,ibl) &
               - c1r*u(6,isite,mu,ieo,ibl) - c1i*u(5,isite,mu,ieo,ibl)
! Normalize the second row.
        c1r = sqrt(u(7 ,isite,mu,ieo,ibl)**2 &
                 + u(8 ,isite,mu,ieo,ibl)**2 &
                 + u(9 ,isite,mu,ieo,ibl)**2 &
                 + u(10,isite,mu,ieo,ibl)**2 &
                 + u(11,isite,mu,ieo,ibl)**2 &
                 + u(12,isite,mu,ieo,ibl)**2)
        do i = 7,12
         u(i,isite,mu,ieo,ibl) = u(i,isite,mu,ieo,ibl)/c1r
        enddo ! i
! Generate the third row: row3 = row1 cross row2.
        u(17,isite,mu,ieo,ibl) = &
            u(1,isite,mu,ieo,ibl)*u(9 ,isite,mu,ieo,ibl) &
           -u(2,isite,mu,ieo,ibl)*u(10,isite,mu,ieo,ibl) &
           -u(3,isite,mu,ieo,ibl)*u(7 ,isite,mu,ieo,ibl) &
           +u(4,isite,mu,ieo,ibl)*u(8 ,isite,mu,ieo,ibl)
        u(15,isite,mu,ieo,ibl) = &
            u(5,isite,mu,ieo,ibl)*u(7 ,isite,mu,ieo,ibl) &
           -u(6,isite,mu,ieo,ibl)*u(8 ,isite,mu,ieo,ibl) &
           -u(1,isite,mu,ieo,ibl)*u(11,isite,mu,ieo,ibl) &
           +u(2,isite,mu,ieo,ibl)*u(12,isite,mu,ieo,ibl)
        u(13,isite,mu,ieo,ibl) = &
            u(3,isite,mu,ieo,ibl)*u(11,isite,mu,ieo,ibl) &
           -u(4,isite,mu,ieo,ibl)*u(12,isite,mu,ieo,ibl) &
           -u(5,isite,mu,ieo,ibl)*u(9 ,isite,mu,ieo,ibl) &
           +u(6,isite,mu,ieo,ibl)*u(10,isite,mu,ieo,ibl)
        u(18,isite,mu,ieo,ibl) = &
           -u(1,isite,mu,ieo,ibl)*u(10,isite,mu,ieo,ibl) &
           -u(2,isite,mu,ieo,ibl)*u(9 ,isite,mu,ieo,ibl) &
           +u(3,isite,mu,ieo,ibl)*u(8 ,isite,mu,ieo,ibl) &
           +u(4,isite,mu,ieo,ibl)*u(7 ,isite,mu,ieo,ibl)
        u(16,isite,mu,ieo,ibl) = &
           -u(5,isite,mu,ieo,ibl)*u(8 ,isite,mu,ieo,ibl) &
           -u(6,isite,mu,ieo,ibl)*u(7 ,isite,mu,ieo,ibl) &
           +u(1,isite,mu,ieo,ibl)*u(12,isite,mu,ieo,ibl) &
           +u(2,isite,mu,ieo,ibl)*u(11,isite,mu,ieo,ibl)
        u(14,isite,mu,ieo,ibl) = &
           -u(3,isite,mu,ieo,ibl)*u(12,isite,mu,ieo,ibl) &
           -u(4,isite,mu,ieo,ibl)*u(11,isite,mu,ieo,ibl) &
           +u(5,isite,mu,ieo,ibl)*u(10,isite,mu,ieo,ibl) &
           +u(6,isite,mu,ieo,ibl)*u(9 ,isite,mu,ieo,ibl)
       enddo ! ibl
      enddo ! isite
     enddo ! mu
    enddo ! ieo
 
 end subroutine fixsu3

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 subroutine aveplaq(plaq,u,myid,nn,ldiv,nms,lv,ib,lbd,iblv,MRT)
! Compute the average plaquette.
! OUTPUT:
!   plaq is ReTr(sum over elementary plaquettes).

    real(kind=KR),    intent(out)                         :: plaq
    real(kind=KR),    intent(inout), dimension(:,:,:,:,:) :: u
    integer(kind=KI), intent(in)                          :: myid, MRT
    integer(kind=KI), intent(in),    dimension(:,:)       :: nn, iblv
    logical,          intent(in),    dimension(:)         :: ldiv
    integer(kind=KI), intent(in),    dimension(:)         :: nms
    integer(kind=KI), intent(in),    dimension(:,:,:)     :: lv
    integer(kind=KI), intent(in),    dimension(:,:,:,:)   :: ib
    logical,          intent(in),    dimension(:,:)       :: lbd

    integer(kind=KI)                    :: mu, ieo, ibl, ierr
    real(kind=KR), dimension(18,nvhalf) :: c
    real(kind=KR)                       :: plaqtot

    plaq = 0.0_KR
    do ieo = 1,2
     do ibl = 1,16
      do mu = 1,4
       call loop(mu,mu,ieo,ibl,c,u,myid,nn,ldiv,nms,lv,ib,lbd,iblv,MRT)
       call trace(mu,ieo,ibl,u,c,plaq)
      enddo ! mu
     enddo ! ibl
    enddo ! ieo
!-combine contributions to the average plaquette from all processes.
! and correctly normalize the average plaquette.
    if (nps==1) then
     plaq = plaq/real(nxyzt*72,KR)
    else
     call MPI_REDUCE(plaq,plaqtot,1,MRT,MPI_SUM,0,MPI_COMM_WORLD,ierr)
     call MPI_BCAST(plaqtot,1,MRT,0,MPI_COMM_WORLD,ierr)
     plaq = plaqtot/real(nxyzt*72,KR)
    endif

 end subroutine aveplaq

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 subroutine tadPcalc(tadmu,u,myid,nn,ldiv,nms,lv,ib,lbd,iblv,MRT)
! Compute the tadpole factor (mean link) in each spacetime direction, using
! the 4th root of an elementary plaquette.
! OUTPUT:
!   plaq is ReTr(sum over elementary plaquettes).

    real(kind=KR),    intent(out),   dimension(:)         :: tadmu
    real(kind=KR),    intent(inout), dimension(:,:,:,:,:) :: u
    integer(kind=KI), intent(in)                          :: myid, MRT
    integer(kind=KI), intent(in),    dimension(:,:)       :: nn, iblv
    logical,          intent(in),    dimension(:)         :: ldiv
    integer(kind=KI), intent(in),    dimension(:)         :: nms
    integer(kind=KI), intent(in),    dimension(:,:,:)     :: lv
    integer(kind=KI), intent(in),    dimension(:,:,:,:)   :: ib
    logical,          intent(in),    dimension(:,:)       :: lbd

    integer(kind=KI)                    :: mu, nu, lambda, ieo, ibl, ierr
    real(kind=KR), dimension(18,nvhalf) :: c
    real(kind=KR)                       :: plaqmu, plaqnu, plaqla, &
                                           plaqmutot, plaqnutot, plaqlatot

    do mu = 1,4
     nu = 1 + modulo(mu,4)
     lambda = 1 + modulo(nu,4)
     plaqmu = 0.0_KR
     plaqnu = 0.0_KR
     plaqla = 0.0_KR
     do ieo = 1,2
      do ibl = 1,16
       call loop(mu,nu,ieo,ibl,c,u,myid,nn,ldiv,nms,lv,ib,lbd,iblv,MRT)
       call trace(mu,ieo,ibl,u,c,plaqmu)
       call loop(nu,lambda,ieo,ibl,c,u,myid,nn,ldiv,nms,lv,ib,lbd,iblv,MRT)
       call trace(nu,ieo,ibl,u,c,plaqnu)
       call loop(lambda,mu,ieo,ibl,c,u,myid,nn,ldiv,nms,lv,ib,lbd,iblv,MRT)
       call trace(lambda,ieo,ibl,u,c,plaqla)
      enddo ! ibl
     enddo ! ieo
!-combine contributions to the average plaquette from all processes.
! and correctly normalize the average plaquette.
     if (nps==1) then
      plaqmutot = plaqmu
      plaqnutot = plaqnu
      plaqlatot = plaqla
     else
      call MPI_REDUCE(plaqmu,plaqmutot,1,MRT,MPI_SUM,0,MPI_COMM_WORLD,ierr)
      call MPI_REDUCE(plaqnu,plaqnutot,1,MRT,MPI_SUM,0,MPI_COMM_WORLD,ierr)
      call MPI_REDUCE(plaqla,plaqlatot,1,MRT,MPI_SUM,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(plaqmutot,1,MRT,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(plaqnutot,1,MRT,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(plaqlatot,1,MRT,0,MPI_COMM_WORLD,ierr)
     endif
!..the double sqrt is because the F compiler can't handle x**(real number)
     tadmu(mu) = sqrt(sqrt(plaqmutot*plaqlatot/plaqnutot/real(6*nxyzt,KR)))
    enddo ! mu

 end subroutine tadPcalc

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 subroutine tadLcalc(tadmu,u,alpha,omega,thetamax,landimp,alat,tad,myid, &
                     nn,ldiv,nms,lv,ib,lbd,iblv,MRT)
! Compute the tadpole factor (mean link) in each spacetime direction, using
! the mean link in Landau gauge.

    real(kind=KR),    intent(out),   dimension(:)         :: tadmu
    real(kind=KR),    intent(inout), dimension(:,:,:,:,:) :: u
    real(kind=KR),    intent(in)                      :: alpha, omega, thetamax
    integer(kind=KI), intent(in),    dimension(:)         :: landimp, nms
    real(kind=KR),    intent(in),    dimension(:)         :: alat, tad
    integer(kind=KI), intent(in)                          :: myid, MRT
    integer(kind=KI), intent(in),    dimension(:,:)       :: nn, iblv
    logical,          intent(in),    dimension(:)         :: ldiv
    integer(kind=KI), intent(in),    dimension(:,:,:)     :: lv
    integer(kind=KI), intent(in),    dimension(:,:,:,:)   :: ib
    logical,          intent(in),    dimension(:,:)       :: lbd

    integer(kind=KI)                             :: ieo, ibl, isite
    real(kind=KR),   dimension(18,ntotal,4,2,16) :: ulandau
    real(kind=KR2),  dimension(4)                :: tadmubit

! Create a working copy of the gauge field configuration.
    ulandau = u

! Gauge fix the working copy of the configuration into Landau gauge.
    call gfixL(ulandau,alpha,omega,thetamax,landimp,alat,tad,myid,nn,ldiv, &
               nms,lv,ib,lbd,iblv,MRT)

! Compute the tadpole factor in each spacetime direction.
    tadmubit = 0.0_KR2
    do ieo = 1,2
     do ibl = 1,16
      do isite = 1,nvhalf
       tadmubit(:) = tadmubit(:) + ulandau( 1,isite,:,ieo,ibl) &
                                 + ulandau( 9,isite,:,ieo,ibl) &
                                 + ulandau(17,isite,:,ieo,ibl)
      enddo ! isite
     enddo ! ibl
    enddo ! ieo
    tadmu = real(tadmubit,KR)/real(3*nxyzt,KR)

 end subroutine tadLcalc

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 subroutine gfixL(u,alpha,omega,thetamax,landimp,alat,tad,myid,nn,ldiv,nms, &
                  lv,ib,lbd,iblv,MRT)
! Apply Landau gauge fixing to the configuration u().
! The gauge fixing condition is improved only in spacetime directions where
! the landimp(mu)=1.
! The overrelaxation method of Mandula and Ogilvie, Phys Lett B248, 156 (1990)
! is used if omega>1.0.
! The Landau gauge fixing is discussed more explicitly by
! Davies et al, Phys Rev D37, 1581 (1988), although their FFT method is not
! used here.  A discussion of improvement can be found in Bonnet et al,
! Austral J Phys 52, 939 (1999).
! INPUT:
!   u() contains the links before gauge fixing.
!   alpha is a tunable parameter: Davies et al suggest alpha=0.1 .
!   omega is a tunable overrelaxation parameter: Mandula and Ogilvie suggest
!         omega=1.7 .  If omega<1.0, then overrelaxation is not used here.
!   thetamax is the limit for an acceptable accuracy in the gauge fixing.
!   landimp(mu)=1 for improvement in the mu direction, else it is 0.
!   alat,tad,myid,nn,ldiv,nms,lv,ib,lbd,iblv as defined in module gauge.
! OUTPUT:
!   u() contains the gauge-fixed links.

    real(kind=KR),    intent(inout), dimension(:,:,:,:,:) :: u
    real(kind=KR),    intent(in)                      :: alpha, omega, thetamax
    integer(kind=KI), intent(in),    dimension(:)         :: landimp, nms
    real(kind=KR),    intent(in),    dimension(:)         :: alat, tad
    integer(kind=KI), intent(in)                          :: myid, MRT
    integer(kind=KI), intent(in),    dimension(:,:)       :: nn, iblv
    logical,          intent(in),    dimension(:)         :: ldiv
    integer(kind=KI), intent(in),    dimension(:,:,:)     :: lv
    integer(kind=KI), intent(in),    dimension(:,:,:,:)   :: ib
    logical,          intent(in),    dimension(:,:)       :: lbd

    integer(kind=KI) :: nimp, mu, ieo, jeo, ibl, jbl, isite, ieo1, jeo1, &
                        ieo2, icount
    real(kind=KR), dimension(18,ntotal,1,2,16) :: G, Gtmp
    real(kind=KR), dimension(18,ntotal,2,16)   :: Dmu
    real(kind=KR), dimension(18,ntotal)        :: Dbit1, Dbit2
    real(kind=KR)                              :: mufac
    real(kind=KR2)                             :: theta

    icount = nvhalf
    maindo: do

! Construct the transformation which moves closer to Landau gauge
! via the steepest-descents method.
! (Use the notation of Davies et al and Bonnet et al, but generalize to
! the anisotropic case and implicitly choose alat(4) for normalization.)
     G = 0.0_KR
     nimp = 0
     do mu = 1,4
      Dmu = 0.0_KR
!-Delta1 is the basic term.
      Dbit1 = 0.0_KR
      do ieo = 1,2
       jeo = 3 - ieo
       do ibl = 1,16
        jbl = iblv(ibl,mu)
        do isite = 1,nvhalf
         Dmu(:,isite,ieo,ibl) = Dmu(:,isite,ieo,ibl) - u(:,isite,mu,ieo,ibl)
        enddo ! isite
        if (lbd(ibl,mu)) then
         do isite = 1,nvhalf
          Dmu(:,isite,ieo,jbl) = Dmu(:,isite,ieo,jbl) + u(:,isite,mu,ieo,ibl)
         enddo ! isite
        else
         do isite = 1,nvhalf
          Dbit1(:,lv(isite,mu,ieo)) = u(:,isite,mu,ieo,ibl)
         enddo ! isite
         if (ldiv(mu)) call slidematrix(Dbit1,ib(:,mu,1,jeo), &
                                        mu,nms,myid,nn,MRT)
         Dmu(:,isite,jeo,jbl) = Dmu(:,isite,jeo,jbl) + Dbit1(:,isite)
        endif
       enddo ! ibl
      enddo ! ieo
!-Delta2 is the improvement term.
      if (landimp(mu)==1) then
       Dbit1 = 0.0_KR
       Dbit2 = 0.0_KR
       nimp = nimp + 1
       do ieo = 1,2
        jeo = 3 - ieo
        do ibl = 1,16
         jbl = iblv(ibl,mu)
         call setmmf(mu,lbd(ibl,mu),u(:,:,mu,ieo,iblv(ibl,mu))        &
                    ,ldiv(mu),ib(:,mu,1,jeo),u(:,:,mu,ieo,ibl)        &
                    ,u(:,:,mu,jeo,iblv(ibl,mu)),Dbit1                 &
                    ,lv(:,mu,ieo),1,nms,myid,nn,MRT)
         Dmu(:,:,ieo,ibl) = Dmu(:,:,ieo,ibl) + Dbit1(:,:)/(16.0_KR*tad(mu))
         if (lbd(ibl,mu)) then
          ieo1 = ieo
          Dbit2 = Dbit1
         else
          ieo1 = jeo
          do isite = 1,nvhalf
           Dbit2(:,lv(isite,mu,ieo)) = Dbit1(:,isite)
          enddo ! isite
          if (ldiv(mu)) call slidematrix(Dbit2,ib(:,mu,1,jeo), &
                                         mu,nms,myid,nn,MRT)
         endif
         jeo1 = 3 - ieo1
         if (lbd(jbl,mu)) then
          ieo2 = ieo1
          Dbit1 = Dbit2
         else
          ieo2 = jeo1
          do isite = 1,nvhalf
           Dbit1(:,lv(isite,mu,ieo1)) = Dbit2(:,isite)
          enddo ! isite
          if (ldiv(mu)) call slidematrix(Dbit1,ib(:,mu,1,jeo1), &
                                         mu,nms,myid,nn,MRT)
         endif
         Dmu(:,:,ieo2,ibl) = Dmu(:,:,ieo2,ibl) - Dbit1(:,:)/(16.0_KR*tad(mu))
        enddo ! ibl
       enddo ! ieo
      endif
!-include the hermitian conjugate and subtract the trace.
      mufac = 1.0_KR/(tad(mu)*alat(mu)**2)
      G( 2,:,1,:,:) = G( 2,:,1,:,:) + 2.0_KR/3.0_KR*mufac &
                           *(2.0_KR*Dmu(2,:,:,:)-Dmu(10,:,:,:)-Dmu(18,:,:,:))
      G( 3,:,1,:,:) = G( 3,:,1,:,:) + mufac*(Dmu( 3,:,:,:)-Dmu( 7,:,:,:))
      G( 4,:,1,:,:) = G( 4,:,1,:,:) + mufac*(Dmu( 4,:,:,:)+Dmu( 8,:,:,:))
      G( 5,:,1,:,:) = G( 5,:,1,:,:) + mufac*(Dmu( 5,:,:,:)-Dmu(13,:,:,:))
      G( 6,:,1,:,:) = G( 6,:,1,:,:) + mufac*(Dmu( 6,:,:,:)+Dmu(14,:,:,:))
      G( 7,:,1,:,:) = G( 7,:,1,:,:) + mufac*(Dmu( 7,:,:,:)-Dmu( 3,:,:,:))
      G( 8,:,1,:,:) = G( 8,:,1,:,:) + mufac*(Dmu( 8,:,:,:)+Dmu( 4,:,:,:))
      G(10,:,1,:,:) = G(10,:,1,:,:) + 2.0_KR/3.0_KR*mufac &
                           *(2.0_KR*Dmu(10,:,:,:)-Dmu(2,:,:,:)-Dmu(18,:,:,:))
      G(11,:,1,:,:) = G(11,:,1,:,:) + mufac*(Dmu(11,:,:,:)-Dmu(15,:,:,:))
      G(12,:,1,:,:) = G(12,:,1,:,:) + mufac*(Dmu(12,:,:,:)+Dmu(16,:,:,:))
      G(13,:,1,:,:) = G(13,:,1,:,:) + mufac*(Dmu(13,:,:,:)-Dmu( 5,:,:,:))
      G(14,:,1,:,:) = G(14,:,1,:,:) + mufac*(Dmu(14,:,:,:)+Dmu( 6,:,:,:))
      G(15,:,1,:,:) = G(15,:,1,:,:) + mufac*(Dmu(15,:,:,:)-Dmu(11,:,:,:))
      G(16,:,1,:,:) = G(16,:,1,:,:) + mufac*(Dmu(16,:,:,:)+Dmu(12,:,:,:))
      G(18,:,1,:,:) = G(18,:,1,:,:) + 2.0_KR/3.0_KR*mufac &
                           *(2.0_KR*Dmu(18,:,:,:)-Dmu(2,:,:,:)-Dmu(10,:,:,:))
     enddo ! mu
     G = 16.0_KR*G/(16.0_KR-real(nimp,KR))

! Compute the criterion for convergence (Davies et al).
     theta = 0.0_KR2
     do ieo = 1,2
      do ibl = 1,16
       call mmd(icount,G(:,:,1,ieo,ibl),G(:,:,1,ieo,ibl),Dbit1)
       do isite = 1,nvhalf
        theta = theta + Dbit1(1,isite) + Dbit1(9,isite) + Dbit1(17,isite)
       enddo ! isite
      enddo ! ibl
     enddo ! ieo
     theta = theta/real(3*nxyzt,KR2)

! Exit the subroutine only when the configuration is sufficiently close
! to Landau gauge.
     if (abs(theta)<thetamax) exit maindo

! Normalize the transformation, add the identity, then reunitarize.
     G = 0.5_KR*alpha*G
     G( 1,:,1,:,:) = 1.0_KR
     G( 9,:,1,:,:) = 1.0_KR
     G(17,:,1,:,:) = 1.0_KR
     call fixsu3(G)

! Perform the optional overrelaxation step (Mandula and Ogilvie).
     if (omega>1.0_KR) then
      Gtmp = G
      do ieo = 1,2
       do ibl = 1,16
        do isite = 1,nvhalf
         Gtmp( 1,isite,1,ieo,ibl) = Gtmp( 1,isite,1,ieo,ibl) - 1.0_KR
         Gtmp( 9,isite,1,ieo,ibl) = Gtmp( 9,isite,1,ieo,ibl) - 1.0_KR
         Gtmp(17,isite,1,ieo,ibl) = Gtmp(17,isite,1,ieo,ibl) - 1.0_KR
        enddo ! isite
        call mm(icount,Gtmp(:,:,1,ieo,ibl),Gtmp(:,:,1,ieo,ibl),G(:,:,1,ieo,ibl))
        G = omega*Gtmp + 0.5*omega*(omega-1.0_KR)*G
        do isite = 1,nvhalf
         G( 1,isite,1,ieo,ibl) = 1.0_KR + G( 1,isite,1,ieo,ibl)
         G( 9,isite,1,ieo,ibl) = 1.0_KR + G( 9,isite,1,ieo,ibl)
         G(17,isite,1,ieo,ibl) = 1.0_KR + G(17,isite,1,ieo,ibl)
        enddo ! isite
       enddo ! ibl
      enddo ! ieo
      call fixsu3(G)
     endif

! Invoke the transformation onto the configuration.
     call rotlinks(u,G(:,:,1,:,:),myid,nn,ldiv,nms,lv,ib,lbd,iblv,MRT)

    enddo maindo

 end subroutine gfixL

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 subroutine fuzz(ttdir,u,nfuzz,epsfuzz,myid,nn,ldiv,nms,lv,ib,lbd,iblv,MRT)
! Perform fuzzing for the link variables "u", and output as "u".
! NOTE: The outputted "u" is projected back onto SU(3), but the projection
!       algorithm is not unique.
! NOTE: epsfuzz is used in each "spatial" direction, despite possible
!       differences in lattice spacings and tadpole factors.

    integer(kind=KI), intent(in)                     :: ttdir, nfuzz, myid, MRT
    real(kind=KR),    intent(inout), dimension(:,:,:,:,:) :: u
    real(kind=KR),    intent(in)                          :: epsfuzz
    integer(kind=KI), intent(in),    dimension(:,:)       :: nn, iblv
    logical,          intent(in),    dimension(:)         :: ldiv
    integer(kind=KI), intent(in),    dimension(:)         :: nms
    integer(kind=KI), intent(in),    dimension(:,:,:)     :: lv
    integer(kind=KI), intent(in),    dimension(:,:,:,:)   :: ib
    logical,          intent(in),    dimension(:,:)       :: lbd

    integer(kind=KI)                           :: ifuzz, ieo, ibl, mu, nu, isite
    real(kind=KR), dimension(18,nvhalf)        :: c
    real(kind=KR), dimension(18,nvhalf,4,2,16) :: ubit

    if (nfuzz>0) then
     do ifuzz = 1,nfuzz
      ubit = 0.0_KR
      do ieo = 1,2
       do ibl = 1,16
        do mu = 1,4
         if (mu/=ttdir) then
          do nu = 1,4
           if (nu/=ttdir.and.nu/=mu) then
            call loop(mu,nu,ieo,ibl,c,u,myid,nn,ldiv,nms,lv,ib,lbd,iblv,MRT)
            ubit(:,:,mu,ieo,ibl) = ubit(:,:,mu,ieo,ibl) + c(:,:)
           endif
          enddo ! nu
         endif
        enddo ! mu
       enddo ! ibl
      enddo ! ieo
      do mu = 1,4
       if (mu/=ttdir) then
        do isite = 1,nvhalf
         u(:,isite,mu,:,:) = u(:,isite,mu,:,:) + epsfuzz*ubit(:,isite,mu,:,:)
        enddo ! isite
       endif
      enddo ! mu
      call fixsu3(u)
     enddo ! ifuzz
    endif

 end subroutine fuzz

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 subroutine wilsonloop(ttdir,leni,pathi,lenf,pathf,wloop,u,myid,nn,ldiv, &
                       nms,lv,ib,lbd,iblv,MRT)
! Compute a Wilson loop with fixed "spatial" cross section, and "temporal"
! extents between 1 and the total "temporal" extent of the lattice.
! Each size of loop is averaged over all sites on the entire configuration.
! INPUT:
!   ttdir is the "time" direction, i.e. the direction in which extended
!         loops will be calculated.
!   leni is the number of links in the initial path (i.e. the line whose
!        endpoints will be extended in the ttdir direction to create a
!        Wilson loop.
!   lenf is the number of links in the final path (i.e. the line whose
!        endpoints will be extended in the ttdir direction to create a
!        Wilson loop.
!   pathi(i) is the spacetime direction of the i'th link of the initial path.
!   pathf(i) is the spacetime direction of the i'th link of the final path.
!   MRT is MPIREALTYPE.
! OUTPUT:
!   wloop(t) is the numerical value of the Wilson loop of "time" length t,
!            averaged over all sites.

    integer(kind=KI), intent(in)                :: ttdir, leni, lenf, myid, MRT
    integer(kind=KI), intent(in),    dimension(:)         :: pathi, pathf
    real(kind=KR),    intent(out),   dimension(:)         :: wloop
    real(kind=KR),    intent(inout), dimension(:,:,:,:,:) :: u
    integer(kind=KI), intent(in),    dimension(:,:)       :: nn, iblv
    logical,          intent(in),    dimension(:)         :: ldiv
    integer(kind=KI), intent(in),    dimension(:)         :: nms
    integer(kind=KI), intent(in),    dimension(:,:,:)     :: lv
    integer(kind=KI), intent(in),    dimension(:,:,:,:)   :: ib
    logical,          intent(in),    dimension(:,:)       :: lbd

    integer(kind=KI) :: ttmax, ttlen, ilen, i, idir, ieo, jeo, ibl, ieo1, &
                        jeo1, ibl1, ieo2, jeo2, ibl2, ieo3, jeo3, ibl3
    real(kind=KR), dimension(18,ntotal) :: uprod1, uprod2, uprod3, uprod4

! Define the maximum "temporal" extent of a Wilson loop.
    if (ttdir==1) ttmax=nx
    if (ttdir==2) ttmax=ny
    if (ttdir==3) ttmax=nz
    if (ttdir==4) ttmax=nt

    wloop = 0.0_KR
    do ieo = 1,2
     jeo = 3 - ieo
     do ibl = 1,16

! Initialize the two halves of the Wilson loop.
      uprod1(:,:) = u(:,:,pathi(1),ieo,ibl)
      uprod2(:,:) = u(:,:,ttdir,ieo,ibl)

! Define "ieo" and "ibl" indices that can move as the Wilson loop grows.
      ieo1 = ieo
      jeo1 = jeo
      ibl1 = ibl
      ieo2 = ieo
      jeo2 = jeo
      ibl2 = ibl

! Extend uprod1 to include the complete "initial time" path plus one link in
! the "time" direction.
      do ilen = 1,leni
       if (ilen<leni) then
        call setmmf(pathi(ilen),lbd(ibl1,pathi(ilen))                       &
                   ,u(:,:,pathi(ilen+1),ieo1,iblv(ibl1,pathi(ilen)))        &
                   ,ldiv(pathi(ilen)),ib(:,pathi(ilen),1,jeo1),uprod1       &
                   ,u(:,:,pathi(ilen+1),jeo1,iblv(ibl1,pathi(ilen))),uprod4 &
                   ,lv(:,pathi(ilen),ieo1),1,nms,myid,nn,MRT)
       else
        call setmmf(pathi(ilen),lbd(ibl1,pathi(ilen))                 &
                   ,u(:,:,ttdir,ieo1,iblv(ibl1,pathi(ilen)))          &
                   ,ldiv(pathi(ilen)),ib(:,pathi(ilen),1,jeo1),uprod1 &
                   ,u(:,:,ttdir,jeo1,iblv(ibl1,pathi(ilen))),uprod4   &
                   ,lv(:,pathi(ilen),ieo1),1,nms,myid,nn,MRT)
       endif
       if (lbd(ibl1,pathi(ilen))) then
        uprod1 = uprod4
       else
        do i = 1,nvhalf
         uprod1(:,lv(i,pathi(ilen),ieo1)) = uprod4(:,i)
        enddo ! i
        if (ldiv(pathi(ilen))) call slidematrix(uprod1 &
                  ,ib(:,pathi(ilen),1,jeo1),pathi(ilen),nms,myid,nn,MRT)
        ieo1 = jeo1
        jeo1 = 3 - ieo1
       endif
       ibl1 = iblv(ibl1,pathi(ilen))
      enddo ! ilen

      do ttlen = 1,ttmax

! Extend uprod1 and uprod2 into the ttdir direction by one additional link.
       if (ttlen>1) then
        call setmmf(ttdir,lbd(ibl1,ttdir)                     &
                   ,u(:,:,ttdir,ieo1,iblv(ibl1,ttdir))        &
                   ,ldiv(ttdir),ib(:,ttdir,1,jeo1),uprod1     &
                   ,u(:,:,ttdir,jeo1,iblv(ibl1,ttdir)),uprod4 &
                   ,lv(:,ttdir,ieo1),1,nms,myid,nn,MRT)
        call setmmf(ttdir,lbd(ibl2,ttdir)                     &
                   ,u(:,:,ttdir,ieo2,iblv(ibl2,ttdir))        &
                   ,ldiv(ttdir),ib(:,ttdir,1,jeo2),uprod2     &
                   ,u(:,:,ttdir,jeo2,iblv(ibl2,ttdir)),uprod3 &
                   ,lv(:,ttdir,ieo2),1,nms,myid,nn,MRT)
        if (lbd(ibl1,ttdir)) then
         uprod1 = uprod4
        else
         do i = 1,nvhalf
          uprod1(:,lv(i,ttdir,ieo1)) = uprod4(:,i)
         enddo ! i
         if (ldiv(ttdir)) call slidematrix(uprod1 &
                     ,ib(:,ttdir,1,jeo1),ttdir,nms,myid,nn,MRT)
         ieo1 = jeo1
         jeo1 = 3 - ieo1
        endif
        ibl1 = iblv(ibl1,ttdir)
        if (lbd(ibl2,ttdir)) then
         uprod2 = uprod3
        else
         do i = 1,nvhalf
          uprod2(:,lv(i,ttdir,ieo2)) = uprod3(:,i)
         enddo ! i
         if (ldiv(ttdir)) call slidematrix(uprod2 &
                     ,ib(:,ttdir,1,jeo2),ttdir,nms,myid,nn,MRT)
         ieo2 = jeo2
         jeo2 = 3 - ieo2
        endif
        ibl2 = iblv(ibl2,ttdir)
       endif

! Extend uprod2 (now called uprod3) to include the complete "final time" path.
       ieo3 = ieo2
       jeo3 = jeo2
       ibl3 = ibl2
       uprod3 = uprod2
       do ilen = 0,lenf-1
        if (ilen==0) then
         call setmmf(ttdir,lbd(ibl3,ttdir)                        &
                    ,u(:,:,pathf(1),ieo3,iblv(ibl3,ttdir))        &
                    ,ldiv(ttdir),ib(:,ttdir,1,jeo3),uprod3        &
                    ,u(:,:,pathf(1),jeo3,iblv(ibl3,ttdir)),uprod4 &
                    ,lv(:,ttdir,ieo3),1,nms,myid,nn,MRT)
         idir = ttdir
        else
         call setmmf(pathf(ilen),lbd(ibl3,pathf(ilen))                       &
                    ,u(:,:,pathf(ilen+1),ieo3,iblv(ibl3,pathf(ilen)))        &
                    ,ldiv(pathf(ilen)),ib(:,pathf(ilen),1,jeo3),uprod3       &
                    ,u(:,:,pathf(ilen+1),jeo3,iblv(ibl3,pathf(ilen))),uprod4 &
                    ,lv(:,pathf(ilen),ieo3),1,nms,myid,nn,MRT)
         idir = pathf(ilen)
        endif
        if (lbd(ibl3,idir)) then
         uprod3 = uprod4
        else
         do i = 1,nvhalf
          uprod3(:,lv(i,idir,ieo3)) = uprod4(:,i)
         enddo ! i
         if (ldiv(idir)) call slidematrix(uprod3 &
                     ,ib(:,idir,1,jeo3),idir,nms,myid,nn,MRT)
         ieo3 = jeo3
         jeo3 = 3 - ieo3
        endif
        ibl3 = iblv(ibl3,idir)
       enddo ! ilen

! Perform a final shift so uprod1 and uprod3 arguments are aligned.
       if (.not.lbd(ibl3,pathf(lenf))) then
        do i = 1,nvhalf
         uprod4(:,lv(i,pathf(lenf),ieo3)) = uprod3(:,i)
        enddo ! i
        if (ldiv(pathf(lenf))) call slidematrix(uprod4 &
                  ,ib(:,pathf(lenf),1,jeo3),pathf(lenf),nms,myid,nn,MRT)
        uprod3 = uprod4
       endif

! Close the Wilson loop and compute ReTr.
       call setmmf(ttdir,lbd(ibl1,ttdir),uprod3                  &
                  ,ldiv(ttdir),ib(:,ttdir,1,jeo1),uprod1         &
                  ,uprod3,uprod4,lv(:,ttdir,ieo1),2,nms,myid,nn,MRT)
       do i = 1,nvhalf
        wloop(ttlen) = wloop(ttlen)+uprod4(1,i)+uprod4(9,i)+uprod4(17,i)
       enddo ! i

      enddo ! ttlen

     enddo ! ibl
    enddo ! ieo
    wloop = wloop/real(3*nxyzt,KR)

 end subroutine wilsonloop

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 subroutine rotlinks(u,rot,myid,nn,ldiv,nms,lv,ib,lbd,iblv,MRT)
! Apply a gauge transformation to the link variables.
! INPUT:
!   rot() contains the SU(3) gauge transformation matrices.
!         expected size: rot(18,ntotal,2,16)
! OUTPUT:
!   The links, u(), are updated.
 
    real(kind=KR),    intent(inout), dimension(:,:,:,:,:) :: u
    real(kind=KR),    intent(inout), dimension(:,:,:,:)   :: rot
    integer(kind=KI), intent(in)                          :: myid, MRT
    integer(kind=KI), intent(in),    dimension(:,:)       :: nn
    logical,          intent(in),    dimension(:)         :: ldiv
    integer(kind=KI), intent(in),    dimension(:)         :: nms
    integer(kind=KI), intent(in),    dimension(:,:,:)     :: lv
    integer(kind=KI), intent(in),    dimension(:,:,:,:)   :: ib
    logical,          intent(in),    dimension(:,:)       :: lbd
    integer(kind=KI), intent(in),    dimension(:,:)       :: iblv
 
    real(kind=KR),   dimension(18,nvhalf) :: tmp1, tmp2
    integer(kind=KI)                      :: ieo, ibl, mu, jeo, isite, icount
 
    icount = nvhalf
    do ieo = 1,2
     jeo = 3 - ieo
     do ibl = 1,16
      do mu = 1,4
       call mm(icount,rot(:,:,ieo,ibl),u(:,:,mu,ieo,ibl),tmp1)
       if (lbd(ibl,mu)) then
        call mmd(icount,tmp1,rot(:,:,ieo,iblv(ibl,mu)),tmp2)
       else
        if (ldiv(mu)) call setlink(rot(:,:,jeo,iblv(ibl,mu)),ib(:,mu,1,jeo), &
                                   mu,nms,myid,nn,MRT)
        call mmdg(icount,tmp1,rot(:,:,jeo,iblv(ibl,mu)),tmp2,lv(:,mu,ieo))
       endif
       do isite = 1,nvhalf
        u(:,isite,mu,ieo,ibl) = tmp2(:,isite)
       enddo ! isite
      enddo ! mu
     enddo ! ibl
    enddo ! ieo
 
 end subroutine rotlinks

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 
 end module gaugetools
 
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
