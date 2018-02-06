 subroutine gmresproject(rwdir,b,xshift,GMRES,resmax,itermin,itercount,u,GeeGooinv, &
                    iflag,kappa,coact,bc,vecbl,vecblinv,myid,nn,ldiv,nms, &
                    lvbc,ib,lbd,iblv,MRT,MRT2,isignal,mvp,gdr)
! GMRES-PROJECT(n,k) matrix inverter of
! Ronald B. Morgan, SIAM Journal on Scientific Computing, 24, 20 (2002).
!
! Gmresproject projects over the approximate eigenvectors found in the 
! deflation section of gmresdr (gmresdrshift). These rojected evectors
! are used in the basis to solve the following right-hand side with the
! usual extraction methods.  
!
! INPUT:
!   b() is the source vector.
!   x() is used as an initial estimate of the true solution vector.
!   GMRES(1)=n in GMRES-DR(n,k): maximum dimension of the subspace.
!   GMRES(2)=k in GMRES-DR(n,k): number of approx eigenvectors kept at restart.
!   resmax is the stopping criterion for the iteration.
!   itermin is the minimum number of iteration required by the user.
!   u() contains the gauge fields for this sublattice.
!   GeeGooinv contains the clover/twisted-mass matrix G on globally-even sites
!             and its inverse on globally-odd sites.
!   iflag=-1 for Wilson or -2 for clover.
!   kappa(1) is 1/(8+2*m_0) with m_0 from eq.(1.1) of JHEP08(2001)058.
!   kappa(2) is mu_q from eq.(1.1) of JHEP08(2001)058.
!   coact(2,4,4) contains the coefficients from the action.
!   bc(mu) = 1,-1,0 for periodic,antiperiodic,fixed boundary conditions
!            respectively.
!   vecbl() defines the global checkerboarding of the lattice.
!           vecbl(1,:) means sites on blocks ibl=1,4,6,7,10,11,13,16.
!           vecbl(2,:) means sites on blocks ibl=2,3,5,8,9,12,14,15.
!   myid = number (i.e. single-integer address) of current process
!          (value = 0, 1, 2, ...).
!   nn(j,1) = single-integer address, on the grid of processes, of the
!             neighbour to myid in the +j direction.
!   nn(j,2) = single-integer address, on the grid of processes, of the
!             neighbour to myid in the -j direction.
!   ldiv(mu) is true if there is more than one process in the mu direction,
!            otherwise it is false.
!   nms(mu) is number of even (or odd) boundary sites per block between
!           processes in the mu'th direction IF there is more than one
!           process in this direction.
!   lvbc(ivhalf,mu,ieo) is the nearest neighbour site in +mu direction.
!                       If there is only one process in the mu direction,
!                       then lvbc still points to the buffer whenever
!                       non-periodic boundary conditions are needed.
!   ib(ibmax,mu,1,ieo) contains the boundary sites at the -mu edge of this
!                      block of the sublattice.
!   ib(ibmax,mu,2,ieo) contains the boundary sites at the +mu edge of this
!                      block of the sublattice.
!   lbd(ibl,mu) is true if ibl is the first of the two blocks in the mu
!               direction, otherwise it is false.
!   iblv(ibl,mu) is the number of the neighbouring block (to the block
!                numbered ibl) in the mu direction.
!                Both ibl and iblv() run from 1 through 16.
!   MRT is MPIREALTYPE.
!   MRT2 is MPIDBLETYPE.
! OUTPUT:
!   x() is the computed solution vector.
!   itercount is the number of iterations used by gmresproject.

! This is PROJ
 
    use shift

    character(len=*), intent(in),    dimension(:)           :: rwdir
    real(kind=KR),    intent(in),    dimension(:,:,:,:,:)   :: b 
    ! real(kind=KR),    intent(inout), dimension(:,:,:,:,:) :: x
    real(kind=KR),    intent(inout), dimension(:,:,:,:,:,:) :: xshift
    integer(kind=KI), intent(in),    dimension(:)         :: GMRES, bc, nms
    integer(kind=KI), intent(in)                          :: itermin, iflag, &
                                                             myid, MRT, MRT2
    real(kind=KR),    intent(in)                          :: resmax
    integer(kind=KI), intent(out)                         :: itercount
    real(kind=KR),    intent(in),    dimension(:,:,:,:,:) :: u
    real(kind=KR),    intent(in),    dimension(:,:,:,:,:) :: GeeGooinv
    real(kind=KR),    intent(in),    dimension(:)         :: kappa
    real(kind=KR),    intent(in),    dimension(:,:,:)     :: coact
    integer(kind=KI), intent(in),    dimension(:,:)       :: vecbl, vecblinv
    integer(kind=KI), intent(in),    dimension(:,:)       :: nn, iblv
    logical,          intent(in),    dimension(:)         :: ldiv
    integer(kind=KI), intent(in),    dimension(:,:,:)     :: lvbc
    integer(kind=KI), intent(in),    dimension(:,:,:,:)   :: ib
    logical,          intent(in),    dimension(:,:)       :: lbd
    ! real(kind=KR2),   intent(in),    dimension(:,:,:)     :: hcnew
    ! real(kind=KR2),   intent(in),    dimension(:,:,:,:,:,:) :: vtemp
!   real(kind=KR2),   intent(inout), dimension(:,:,:,:,:) :: beg
    integer(kind=KI), intent(in)                          :: isignal
    integer(kind=KI), intent(inout)                       :: mvp
    real(kind=KR2),   intent(inout), dimension(:,:)       :: gdr
 
    integer(kind=KI) :: icycle, i, j, k, jp1, jj, p, ir, irm1, is, ivb, ii, &
                        idis, idag, icri, nDR, kDR, ilo, ihi, ischur, &
                        id, ieo, ibleo, ikappa, nkappa, ifreq,temprhs !,nrhs
 
    ! real(kind=KR),    dimension(6,ntotal,4,2,8,nshifts)     :: xshift
    logical,          dimension(nmaxGMRES)                  :: myselect
    integer(kind=KI), dimension(nmaxGMRES)                  :: ipiv,ierr
    real(kind=KR2)                                          :: const, tval, &
                                                               con2, rv, amags
    real(kind=KR2)                                          :: rn, rinit, rnnn 
    integer(kind=KI)                                        :: ldh, ldz, ldhcht, &
                                                               lwork, lzwork, info,rnt
    real(kind=KR2),   dimension(2)                          :: beta, alpha, &
                                                               tv, tv1 ,tv2
    real(kind=KR2),   dimension(2)                          :: beta1, beta2, beta3
    real(kind=KR2),   dimension(2,nshifts)                  :: betashift
    real(kind=KR2),   dimension(kcyclim,nshifts)            :: rnale
    real(kind=KR2),   dimension(nmaxGMRES)                  :: mag
    real(kind=KR2),   dimension(nmaxGMRES)                  :: sr, si
    real(kind=KR2),   dimension(2,nmaxGMRES)                :: ss, gs, gc, w
    real(kind=KR2),   dimension(2,nmaxGMRES)                :: tau, work
    complex(kind=KCC), dimension(nmaxGMRES+1)                :: ztau
    complex(kind=KCC), dimension(nmaxGMRES+1)                :: zwork
    real(kind=KR2),   dimension(nshifts)                    :: sigma
    real(kind=KR2),   dimension(2, nshifts)                 :: alph, cmult
    real(kind=KR2),   dimension(2,nmaxGMRES,nshifts)        :: d
    complex(kind=KCC), dimension(nmaxGMRES,nshifts)          :: zd
    real(kind=KR2),   dimension(2,nmaxGMRES+1)              :: st
    real(kind=KR2),   dimension(2,nmaxGMRES+1,nshifts)      :: srv
    complex(kind=KCC), dimension(nmaxGMRES+1)                :: zsrvrot
    real(kind=KR2),   dimension(2,nmaxGMRES+1)              :: srvis
    complex(kind=KCC), dimension(nmaxGMRES+1)                :: zsrvis
    real(kind=KR2),   dimension(2,nmaxGMRES+1)              :: c, c2
    real(kind=KR2),   dimension(2,nmaxGMRES)                :: cmas
    complex(kind=KCC), dimension(nmaxGMRES+1)                :: zcrot
    real(kind=KR2),   dimension(2,nmaxGMRES,1)              :: em
    real(kind=KR2),   dimension(2,nmaxGMRES+1,nmaxGMRES+1)  :: z
    real(kind=KR2),   dimension(2,nmaxGMRES,nmaxGMRES+1)    :: hcht
    real(kind=KR2),   dimension(2,nmaxGMRES+1,nmaxGMRES)    :: hc, hc2, hc3
    real(kind=KR2),   dimension(2,nmaxGMRES+1,nmaxGMRES)    :: hcs2
    real(kind=KR2),   dimension(2,nmaxGMRES+1,nmaxGMRES+1)  :: hcs
    complex(kind=KCC), dimension(nmaxGMRES+1,nmaxGMRES+1)    :: zhcs
    complex(kind=KCC), dimension(nmaxGMRES+1,nmaxGMRES)      :: zhc2
    complex(kind=KCC), dimension(nmaxGMRES+1,nmaxGMRES)      :: zhcnew
    real(kind=KR2),   dimension(2,nmaxGMRES+1,nmaxGMRES)    :: rr
    complex(kind=KCC), dimension(nmaxGMRES+1,nmaxGMRES)      :: zrr
    real(kind=KR2),   dimension(2,nmaxGMRES+1,kmaxGMRES)    :: ws
    real(kind=KR2),   dimension(2,kmaxGMRES+1,kmaxGMRES)    :: gca, gsa
!    real(kind=KR2),   dimension(6,ntotal,4,2,8,nshifts)     :: xt,xbt
!    real(kind=KR2),   dimension(6,ntotal,4,2,8)             :: xb
    real(kind=KR2),   dimension(6,nvhalf,4,2,8)             :: r, h
    real(kind=KR2),   dimension(6,nvhalf,4,2,8,nshifts)       :: rshift
    !real(kind=KR2),   dimension(6,ntotal,4,2,8,nmaxGMRES+1) :: v
    real(kind=KR2),   dimension(6,nmaxGMRES+1)              :: vt
    ! real(kind=KR2),   dimension(2,nmaxGMRES)                :: gc, gs
    real(kind=KR2),   dimension(2)                          :: gam

    complex(kind=KCC)                                        :: ztemp1, ztemp2, ztemp3, zalpha

    integer(kind=KI)                                        ::  isite, icolorir, idirac, iblock, &
                                                                site, icolorr, irow, ishift

! This is still PROJ 

 
! Shift sigmamu to base mu above (mtmqcd(1,2))
    p = 6!order of the polynomial
    sigma = 0.0_KR
    y = 0.0_KR2!put y=0
    try = 0.0_KR2!put try=0
! DEAN ~ HEY! I need to take out the sigma(1) part because I am not 
!        shifting at all and the residuals should be r = b-Ax
 
! init all x to 0

      xshift = 0.0_KR
 
! Define some parameters.
    nDR = GMRES(1)
    kDR = GMRES(2)


    ldh = nmaxGMRES+ 1
    ldz = nmaxGMRES+ 1
    lzwork = nmaxGMRES+1

!   if (myid==0) then
!     print *, "nDR = ", nDR
!     print *, "kDR = ", kDR
!   endif ! myid

    icycle = 1
    ifreq = 1
    idag = 0
    alph = 0.0_KR
    zcrot = 0.0_KR
    zsrvrot = 0.0_KR
    cmas = 0.0_KR
    srv = 0.0_KR
    srvis = 0.0_KR
    ss = 0.0_KR
    st = 0.0_KR
    hc = 0.0_KR
    hc2 = 0.0_KR
    hc3 = 0.0_KR
    hcs = 0.0_KR
    hcht = 0.0_KR
    zrr = 0.0_KR
    c2 = 0.0_KR
    c = 0.0_KR
    v = 0.0_KR
    xb = 0.0_KR
    d = 0.0_KR
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!poly parameters!!!!!!!!!!!!!!!!!!!!!
    lsmat = 0.0_KR2
    co = 0.0_KR2
    cls = 0.0_KR2
    ipiv2 = 0.0_KI
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Initialize cmult
     
!   do is = 1,nshifts
       is =1
       cmult(1, is) = 1.0_KR
       cmult(2, is) = 0.0_KR
!   enddo ! is
 
    ! do is=1,nshifts
    !   xshift(:,:,:,:,:,is) = x(:,:,:,:,:,:)
    ! enddo ! is
 
! Compute r=b-M*x and v=r/|r| and beta=|r|.

    call vecdot(b(:,:,:,:,:), b(:,:,:,:,:), beta,MRT2)
   !if (myid == 0) then
   ! print *, "norm o b in project =", sqrt(beta(1))
   !endif

! do iblock =1,8
!   do ieo = 1,2
!     do idirac=1,4
!       do isite=1,nvhalf 
!         do icolorir=1,5,2
!                icolorr = icolorir/2 +1
!
! To print single rhs source vector use ..
!
!                irow = icolorr + 3*(isite -1) + 24*(idirac -1) + 96*(ieo -1) + 192*(iblock - 1)
!                       print *, irow, b(icolorir,isite,idirac,ieo,iblock), b(icolorir+1,isite,idirac,ieo,iblock)
!
!             enddo ! icolorir
!          enddo ! isite
!       enddo ! idirac 
!    enddo ! ieo
!  enddo ! iblock

    call Hdbletm(h,u,GeeGooinv,xshift(:,:,:,:,:,1),idag,coact,kappa,iflag,bc,vecbl, &
                    vecblinv,myid,nn,ldiv,nms,lvbc,ib,lbd,iblv,MRT)

! NOTE ~ I don't think that I need this ....

    do ii = 1,nvhalf
        rshift(:,ii,:,:,:,1) = b(:,ii,:,:,:) - h(:,ii,:,:,:) +sigma(1)*xshift(:,ii,:,:,:,1)
    enddo ! ii

! For error correction put error in one direction and use gmresproject to 
! correct and solve the later right hand sides.

    if (isignal == 1) then
!     do is=1,nshifts
       rshift(:,:,:,:,:,1) = beg(:,:,:,:,:)
!     enddo ! is
    else
      beg(:,:,:,:,:) = rshift(:,:,:,:,:,1)
     !beg(:,:,:,:,:) = b(:,:,:,:,:)
    endif
!   call checkNonZero(beg,nvhalf)
! Create rinit for logic passing used in gmresprojet    

   beta = 0.0_KR

   call vecdot(rshift(:,:,:,:,:,1), rshift(:,:,:,:,:,1), beta, MRT2)

   rinit = sqrt(beta(1))
   rn = rinit

   call Hdbletm(h,u,GeeGooinv,xshift(:,:,:,:,:,1),idag,coact,kappa,iflag,bc,vecbl, &
                   vecblinv,myid,nn,ldiv,nms,lvbc,ib,lbd,iblv,MRT)
 
! should b be vtemp because I am solving the error solution???
    
    do i = 1,nvhalf
!       r(:,i,:,:,:) = b(:,i,:,:,:) - h(:,i,:,:,:)

! NOTE~ need to add sigma(1)*xshift(...,1) to end of this even though sigma1=0

        r(:,i,:,:,:) = beg(:,i,:,:,:) - h(:,i,:,:,:) + sigma(1)*xshift(:,i,:,:,:,1)
!      v(:,i,:,:,:,:,1) = r(:,i,:,:,:,:)

!      do is=1,nshifts
!         xt(:,i,:,:,:,is) = x(:,i,:,:,:,:)
!      enddo ! is
    enddo ! i
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!generate the poly!!!!!!!!!!!!!!!!!!!!!!!!!!
    do k = 1,nvhalf
     vprime(:,k,:,:,:,1) = b(:,k,:,:,:)
    enddo !k
    do i = 1,p
     call Hdbletm(vprime(:,:,:,:,:,i+1),u,GeeGooinv,vprime(:,:,:,:,:,i),idag, &
                 coact,kappa,iflag,bc,vecbl,vecblinv,myid,nn, &
                 ldiv,nms,lvbc,ib,lbd,iblv,MRT)
    enddo !i
    
    do i=2,p+1
     do j=2,p+1
      call vecdot(vprime(:,:,:,:,:,i),vprime(:,:,:,:,:,j),beta,MRT2)
      lsmat(:,i-1,j-1) = beta(:)  !lsmat(2,p,p) ,cls(2,p,1)
!      print *, "i,j, lsmat(:,i,j)=", i-1,j-1, lsmat(:,i-1,j-1)
     enddo!j
    enddo!i
        
   do i=2,p+1
     call vecdot(vprime(:,:,:,:,:,i),phi(:,:,:,:,:),beta,MRT2)
     cls(:,i-1,1) = beta(:)
!     print *, "i,cls(:,i)=", i-1, cls(:,i-1,1)
   enddo!i
    
    call linearsolver(p,1,lsmat,ipiv2,cls)
    co(:,:) = cls(:,:,1)    
!    co = 0.0_KR2    
!    co(1,1) = 4
   if(myid==0) then
    do i=1,p
     print *, "i,result from the project(:,i)=", i, co(:,i)
    enddo!i  
   endif!myid

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
! Copy the initial resiudal into the initial residual for each shift

! r and rshift(...,1) are the source vector at this point

!   do is=2,nshifts
       is =1
       rshift(:,:,:,:,:,is) = r(:,:,:,:,:)
 !!!!!!!!!!!!!!!!!!!!!!!r=P(A)*r!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      do ii = 1,nvhalf
          try(:,ii,:,:,:,1) = rshift(:,ii,:,:,:,1)!initiate
      enddo!ii

      do icri=1,5,2
       do k=1,nvhalf
        y(icri,k,:,:,:,1) = co(1,1)*try(icri,k,:,:,:,1) &
                           -co(2,1)*try(icri+1,k,:,:,:,1)
        y(icri+1,k,:,:,:,1) = co(1,1)*try(icri+1,k,:,:,:,1) &
                             +co(2,1)*try(icri,k,:,:,:,1)
       enddo!k
      enddo!icri

      do i=1,p-1 
       call Hdbletm(try(:,:,:,:,:,i+1),u,GeeGooinv,try(:,:,:,:,:,i),idag, &
                    coact,kappa,iflag,bc,vecbl,vecblinv,myid,nn,ldiv,nms, &
                    lvbc,ib,lbd,iblv,MRT )  !z1=M*z
       do icri=1,5,2
        do k=1,nvhalf
         y(icri  ,k,:,:,:,1) = y(icri ,k,:,:,:,1) &
                               +co(1,i+1)*try(icri,k,:,:,:,i+1) &
                               -co(2,i+1)*try(icri+1,k,:,:,:,i+1)
         y(icri+1,k,:,:,:,1) = y(icri+1,k,:,:,:,1) &
                               +co(1,i+1)*try(icri+1,k,:,:,:,i+1) &
                               +co(2,i+1)*try(icri,k,:,:,:,i+1)   !y=P(A)*r
        enddo!k
       enddo!icri
      enddo!i
      do k=1,nvhalf
       rshift(:,k,:,:,:,1) = y(:,k,:,:,:,1) ! Let rshift=P(A)*rshift
      enddo!k
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!r=r*P(A)!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!      rshift(:,:,:,:,:,is) = rshift(:,:,:,:,:,1)
!   enddo ! is
 
! Need to zero out the first shift after creating the first res.
! so that the solution from the projection section is not initiated
! with something other than zeros.

! ~NOTE this should not be done here, but not effetcing hopefully

    xshift(:,:,:,:,:,1) = 0.0_KR
   !if (myid == 0) then
   ! print *, "xshift2  in proj="
!    call checkNonZero(xshift(:,:,:,:,:,2),ntotal)
   !endif

!....Need logic here to determine the cycle in which it leaves...

    k = kDR

    cycledo: do

     if (rn/rinit <= resmax .or. icycle > kcyclim ) exit cycledo

     !if (myid==0) then
    !!   print *, "At start of loop PROJ rn/rinit = ", rn, rinit, rn/rinit
    ! endif

! DEAN~ By uncommenting next line - take out projection 
!     if (icycle -1 ==-1 )then

! The next if statment allows the projection step to occur.
     if (icycle-1 == ((icycle-1)/ifreq)*ifreq) then

      do i=1,k+1
        call vecdot(vtemp(:,:,:,:,:,i),rshift(:,:,:,:,:,1),beta,MRT2)
         c(1,i) = beta(1)
         c(2,i) = beta(2)
         c2(1,i) = c(1,i)
         c2(2,i) = c(2,i)
      enddo ! i

      do ii=1,k+1
        do jj=1,k 
          hc2(1,ii,jj) = hcnew(1,ii,jj)
          hc2(2,ii,jj) = hcnew(2,ii,jj)
        enddo ! jj
      enddo ! ii
 
      do jj =1,k
       hc2(1,jj,jj) = hc2(1,jj,jj) - sigma(1)
      enddo ! jj

      do jj = 1,k
         do i = jj+1,k+1
            amags = hc2(1,jj,jj)**2 + hc2(2,jj,jj)**2
            con2 = 1.0_KR2/amags
            tv(1) = sqrt(amags+hc2(1,i,jj)**2+hc2(2,i,jj)**2)
            tv(2) = 0.0_KR2
            gca(1,i,jj) = sqrt(amags)/tv(1)
            gca(2,i,jj) = 0.0_KR2
            gsa(1,i,jj) = gca(1,i,jj)*con2 &
                          *(hc2(1,i,jj)*hc2(1,jj,jj)+hc2(2,i,jj)*hc2(2,jj,jj))
            gsa(2,i,jj) = gca(1,i,jj)*con2 &
                        *(hc2(2,i,jj)*hc2(1,jj,jj)-hc2(1,i,jj)*hc2(2,jj,jj))
            do j = jj,k
               tv1(1) = gca(1,i,jj)*hc2(1,jj,j) + gsa(1,i,jj)*hc2(1,i,j) &
                                             + gsa(2,i,jj)*hc2(2,i,j)
               tv1(2) = gca(1,i,jj)*hc2(2,jj,j) + gsa(1,i,jj)*hc2(2,i,j) &
                                             - gsa(2,i,jj)*hc2(1,i,j)
               tv2(1) = gca(1,i,jj)*hc2(1,i,j) - gsa(1,i,jj)*hc2(1,jj,j) &
                                            + gsa(2,i,jj)*hc2(2,jj,j)
               tv2(2) = gca(1,i,jj)*hc2(2,i,j) - gsa(1,i,jj)*hc2(2,jj,j) &
                                            - gsa(2,i,jj)*hc2(1,jj,j)
               hc2(:,jj,j) = tv1(:)
               hc2(:,i,j) = tv2(:)
            enddo ! j
            tv1(1) = gca(1,i,jj)*c(1,jj) + gsa(1,i,jj)*c(1,i) + gsa(2,i,jj)*c(2,i)
            tv1(2) = gca(1,i,jj)*c(2,jj) + gsa(1,i,jj)*c(2,i) - gsa(2,i,jj)*c(1,i)
            tv2(1) = gca(1,i,jj)*c(1,i) - gsa(1,i,jj)*c(1,jj) + gsa(2,i,jj)*c(2,jj)
            tv2(2) = gca(1,i,jj)*c(2,i) - gsa(1,i,jj)*c(2,jj) - gsa(2,i,jj)*c(1,jj)
            c(:,jj) = tv1(:)
            c(:,i) = tv2(:)
         enddo ! i
      enddo ! jj

    ! Solve linear equation
! ~ NOTE check dimension of ss...does it need to be zeroed out?

      do i = 1,k
         ss(:,i) = c(:,i)
      enddo ! i

      con2 = 1.0_KR2/(hc2(1,k,k)**2+hc2(2,k,k)**2)
      const = con2*(ss(1,k)*hc2(1,k,k)+ss(2,k)*hc2(2,k,k))
      ss(2,k) = con2*(ss(2,k)*hc2(1,k,k)-ss(1,k)*hc2(2,k,k))
      ss(1,k) = const

      if (k/=1) then
         do i = 1,k-1
            ir = k - i + 1
            irm1 = ir - 1
            do jj = 1,irm1
               const = ss(1,jj) - ss(1,ir)*hc2(1,jj,ir) + ss(2,ir)*hc2(2,jj,ir)
               ss(2,jj) = ss(2,jj) - ss(1,ir)*hc2(2,jj,ir) - ss(2,ir)*hc2(1,jj,ir)
               ss(1,jj) = const
            enddo ! jj
            con2 = 1.0_KR2/(hc2(1,irm1,irm1)**2+hc2(2,irm1,irm1)**2)
            const = con2*(ss(1,irm1)*hc2(1,irm1,irm1)+ss(2,irm1)*hc2(2,irm1,irm1))
            ss(2,irm1) = con2*(ss(2,irm1)*hc2(1,irm1,irm1)-ss(1,irm1)*hc2(2,irm1,irm1))
            ss(1,irm1) = const
         enddo ! i
      endif ! (k/=1)

    ! ... Define new variable d to assist in shifting masses
    ! d is the "short" solution vector to small problem

      do jj=1,k
        d(1,jj,1) = ss(1,jj)
        d(2,jj,1) = ss(2,jj)
      enddo ! jj

! Put this in to take out of algorithum...temporarily. ONCE
!   IT WORKS NEED TO TAKE OUT THE is LOOP!!!

!     is=1
!     if(is==0) then
!
!     do is=2,nshifts
!
!        do ii=1,k+1
!           do jj=1,k
!              hc2(1,ii,jj) = hcnew(1,ii,jj)
!              hc2(2,ii,jj) = hcnew(2,ii,jj)
!           enddo ! jj
!        enddo ! ii
!
!        do ii=1,k+1
!          st(:,ii) = 0.0_KR
!        enddo ! ii
!
! NOTE ~ if first shift is not zero we need hc2 to be shifted...
!
!        do ii = 1,k+1
!           do jj = 1,k
!             st(1,ii) = st(1,ii) + hc2(1,ii,jj)*d(1,jj,1) - hc2(2,ii,jj)*d(2,jj,1)
!             st(2,ii) = st(2,ii) + hc2(1,ii,jj)*d(2,jj,1) + hc2(2,ii,jj)*d(1,jj,1)
!           enddo ! jj
!           srv(1,ii,1) = c2(1,ii) - st(1,ii)
!           srv(2,ii,1) = c2(2,ii) - st(2,ii)
!        enddo ! ii
!
!        do jj=1,k
!          hc2(1,jj,jj) = hc2(1,jj,jj) - sigma(is)
!        enddo ! jj
!
! NOTE ~ d here is a "work" vector give differenet temp name..
!        NOT the short solution vector!
!
!        do ii=1,k+1
!           d(1,ii,is) = cmult(1,is)*st(1,ii) - cmult(2,is)*st(2,ii)
!           d(2,ii,is) = cmult(1,is)*st(2,ii) + cmult(2,is)*st(1,ii)
!        enddo ! ii
!
!       !if (myid == 0) then
!       ! print *, "is, ourrhshereis =", is, d(:,1:2,is)
!       ! print *, "hc2 before rotation =", hc2(:,1:3,1:2)
!       ! print *, " k before =", k
!       !endif ! myid
!
!        do jj = 1,k
!           do i = jj+1,k
!              amags = hc2(1,jj,jj)**2 + hc2(2,jj,jj)**2
!              con2 = 1.0_KR2/amags
!              tv(1) = sqrt(amags+hc2(1,i,jj)**2+hc2(2,i,jj)**2)
!              tv(2) = 0.0_KR2
!              gca(1,i,jj) = sqrt(amags)/tv(1)
!              gca(2,i,jj) = 0.0_KR2
!              gsa(1,i,jj) = gca(1,i,jj)*con2 &
!                            *(hc2(1,i,jj)*hc2(1,jj,jj)+hc2(2,i,jj)*hc2(2,jj,jj))
!              gsa(2,i,jj) = gca(1,i,jj)*con2 &
!                            *(hc2(2,i,jj)*hc2(1,jj,jj)-hc2(1,i,jj)*hc2(2,jj,jj))
!              do j = jj,k
!                 tv1(1) = gca(1,i,jj)*hc2(1,jj,j) + gsa(1,i,jj)*hc2(1,i,j) &
!                                                  + gsa(2,i,jj)*hc2(2,i,j)
!                 tv1(2) = gca(1,i,jj)*hc2(2,jj,j) + gsa(1,i,jj)*hc2(2,i,j) &
!                                                  - gsa(2,i,jj)*hc2(1,i,j)
!                 tv2(1) = gca(1,i,jj)*hc2(1,i,j) - gsa(1,i,jj)*hc2(1,jj,j) &
!                                                 + gsa(2,i,jj)*hc2(2,jj,j)
!                 tv2(2) = gca(1,i,jj)*hc2(2,i,j) - gsa(1,i,jj)*hc2(2,jj,j) &
!                                                 - gsa(2,i,jj)*hc2(1,jj,j)
!                 hc2(:,jj,j) = tv1(:)
!                 hc2(:,i,j) = tv2(:)
!              enddo ! j
!              tv1(1) = gca(1,i,jj)*d(1,jj,is) + gsa(1,i,jj)*d(1,i,is) + gsa(2,i,jj)*d(2,i,is)
!              tv1(2) = gca(1,i,jj)*d(2,jj,is) + gsa(1,i,jj)*d(2,i,is) - gsa(2,i,jj)*d(1,i,is)
!              tv2(1) = gca(1,i,jj)*d(1,i,is) - gsa(1,i,jj)*d(1,jj,is) + gsa(2,i,jj)*d(2,jj,is)
!              tv2(2) = gca(1,i,jj)*d(2,i,is) - gsa(1,i,jj)*d(2,jj,is) - gsa(2,i,jj)*d(1,jj,is)
!!             d(:,jj,is) = tv1(:)
!              d(:,i,is) = tv2(:)
!           enddo ! i
!        enddo ! jj
!
! NOTE~ why are we zeroing out lower triangluar part of hc2
!       when these values should be allready zero.
!
!        do jj=1,k-1
!           do ii=jj+1,k
!              hc2(1,ii,jj) = 0.0_KR2
!              hc2(2,ii,jj) = 0.0_KR2
!           enddo ! ii
!        enddo ! jj
!
!        do ii=1,k
!           hc2(1,k+1,ii) = 0.0_KR2
!           hc2(2,k+1,ii) = 0.0_KR2
!        enddo ! ii
!
!      ! if (myid == 0) then
!      !  print *, "is, d in between =", is, d(:,1:2,is)
!      !  print *, "k in da' middle =", k
!      ! endif ! myid 
! NOTE~ back solve these linear equations
!
!     con2 = 1.0_KR2/(hc2(1,k,k)**2+hc2(2,k,k)**2)
!     const = con2*(d(1,k,is)*hc2(1,k,k)+d(2,k,is)*hc2(2,k,k))
!     d(2,k,is) = con2*(d(2,k,is)*hc2(1,k,k)-d(1,k,is)*hc2(2,k,k))
!     d(1,k,is) = const
!
!     if (k/=1) then
!        do i = 1,k-1
!           ir = k - i + 1
!           irm1 = ir - 1
!           do jj = 1,irm1
!              const = d(1,jj,is) - d(1,ir,is)*hc2(1,jj,ir) + d(2,ir,is)*hc2(2,jj,ir)
!              d(2,jj,is) = d(2,jj,is) - d(1,ir,is)*hc2(2,jj,ir) - d(2,ir,is)*hc2(1,jj,ir)
!              d(1,jj,is) = const
!           enddo ! jj
!           con2 = 1.0_KR2/(hc2(1,irm1,irm1)**2+hc2(2,irm1,irm1)**2)
!           const = con2*(d(1,irm1,is)*hc2(1,irm1,irm1)+d(2,irm1,is)*hc2(2,irm1,irm1))
!           d(2,irm1,is) = con2*(d(2,irm1,is)*hc2(1,irm1,irm1)-d(1,irm1,is)*hc2(2,irm1,irm1))
!           d(1,irm1,is) = const
!        enddo ! i
!     endif ! (k/=1)
!
!    !if (myid == 0) then
!    ! print *, "is, oursolnvector =" , is, d(:,1:2,is)
!    ! print *, "k after =", k
!    !endif ! myid
!
! NOTE ~ took out least squares solution upon DR Morgan request..
!
!        call real2complex_mat(hc2, nmaxGMRES+1, nmaxGMRES, zhc2)
!
!        call zgels('N', k+1, k, 1, zhc2, nmaxGMRES+1, zd(1,is), nmaxGMRES, zwork, lzwork, info)
!
!        call complex2real_mat(zd,nmaxGMRES,nshifts,d)
!        call complex2real_mat(zhc2, nmaxGMRES+1,nmaxGMRES, hc2)
!        
!     enddo ! is
!     
!     endif ! (is==0)
     
! Form the approximate new solution x = xt + xb.
! First zero xb then xb = V*d(1,2,:,is), and x =xt +xb
! ... Need to keep track of each solution for each shift.
!     Loop over shifts while creating the soln vector.
 
! DEAN ~ ONCE THIS ALGORITHUM WORKS, NEED TO TAKE OUT THE
!        is  LOOP AND HARD CODE IN is==1

 !          do is = 1,nshifts
               is = 1
               do i=1,k
                  sr(i) = d(1,i,is)
                  si(i) = d(2,i,is)
               enddo ! i
               do jj=1,nvhalf
                  xb(:,jj,:,:,:) = 0.0_KR2
               enddo ! jj
               do jj = 1,k
                  do icri = 1,5,2 ! 6=nri*nc
                     do i = 1,nvhalf
                        xb(icri  ,i,:,:,:) = xb(icri  ,i,:,:,:) &
                                           + sr(jj)*vtemp(icri  ,i,:,:,:,jj) &
                                           - si(jj)*vtemp(icri+1,i,:,:,:,jj)
                        xb(icri+1,i,:,:,:) = xb(icri+1,i,:,:,:) &
                                           + si(jj)*vtemp(icri  ,i,:,:,:,jj) &
                                           + sr(jj)*vtemp(icri+1,i,:,:,:,jj)
                     enddo ! i
                  enddo ! icri
               enddo ! jj
               do i = 1,nvhalf
                  xshift(:,i,:,:,:,is) = xshift(:,i,:,:,:,is) + xb(:,i,:,:,:)
               enddo ! i

! This call to Hdbletm is purely for informational purpose.

!               call Hdbletm(h,u,GeeGooinv,xshift(:,:,:,:,:,is),idag,coact,kappa,iflag, &
!                               bc,vecbl,vecblinv,myid,nn,ldiv,nms,lvbc,ib,lbd,iblv,MRT)
!               do ii = 1,nvhalf
!                  rshift(:,ii,:,:,:,is) = beg(:,ii,:,:,:) - h(:,ii,:,:,:) &
!                                                +sigma(is)*xshift(:,ii,:,:,:,is)
!               enddo ! ii
!           enddo ! is

!DD note ~ equivqlent to exiting the 505 loop.

   endif ! (icycle-1 == ((icycle-1)/ifreq)*ifreq) 

  !if (myid == 0) then
  ! print *, "cycles of gmres"
  !endif

! Normalize rshift(:,:,:,:,:,1) and put into the first coln. of V, the 
! orthonormal matrix whose colns span the Krylov subspace.

! NOTE ~ may not need all the res norms for shifts >=2

   betashift = 0.0_KR

!    do is =1,nshifts
     is=1

       call vecdot(rshift(:,:,:,:,:,1),rshift(:,:,:,:,:,1), beta,MRT2)
       betashift(1,is) = beta(1)
       betashift(2,is) = beta(2)


!    enddo ! is

   const = 1.0_KR/sqrt(betashift(1,1))


! With the normalized residual form the first coln of V.
   do jj=1,nvhalf
     v(:,jj,:,:,:,1) = const*rshift(:,jj,:,:,:,1)
   enddo ! jj  

   ztemp1 = DCMPLX(betashift(1,1),betashift(2,1))
   ztemp1 = sqrt(ztemp1)
   c(1,1) = REAL(ztemp1)
   c(2,1) = AIMAG(ztemp1)
   c2(1,1) = c(1,1)
   c2(2,1) = c(2,1)
   
! Perform GMRES between projections...

     do j = 1,nDR
       jp1 = j + 1
       itercount = itercount + 1
 
! LEFT OFF HERE!
 
!*****MORGAN'S STEP 2A AND STEPS 7,8A: Apply standard GMRES(nDR).
! Generate V_(m+1) and Hbar_m with the Arnoldi iteration.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!v_j=P(A)*v_j!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    try = 0.0_KR2
    y = 0.0_KR2
    do i=1,nvhalf
     try(:,i,:,:,:,1) = v(:,i,:,:,:,j)
    enddo!i
    do icri=1,5,2
     do k=1,nvhalf
      y(icri,k,:,:,:,1) = co(1,1)*try(icri,k,:,:,:,1) &
                         -co(2,1)*try(icri+1,k,:,:,:,1)
      y(icri+1,k,:,:,:,1) = co(1,1)*try(icri+1,k,:,:,:,1) &
                           +co(2,1)*try(icri,k,:,:,:,1)
     enddo!k
    enddo!icri
    do i=1,p-1 
     call Hdbletm(try(:,:,:,:,:,i+1),u,GeeGooinv,try(:,:,:,:,:,i),idag,coact, &
                  kappa,iflag,bc,vecbl,vecblinv,myid,nn,ldiv,nms,lvbc,ib, &
                  lbd,iblv,MRT )  !z1=M*z
     do icri=1,5,2
      do k=1,nvhalf
       y(icri  ,k,:,:,:,1) = y(icri ,k,:,:,:,1) &
                             +co(1,i+1)*try(icri,k,:,:,:,i+1) &
                             -co(2,i+1)*try(icri+1,k,:,:,:,i+1)
       y(icri+1,k,:,:,:,1) = y(icri+1,k,:,:,:,1) &
                             +co(1,i+1)*try(icri+1,k,:,:,:,i+1) &
                             +co(2,i+1)*try(icri,k,:,:,:,i+1)   !y=P(A)*r
      enddo!k
     enddo!icri
    enddo!i

    call Hdbletm(v(:,:,:,:,:,jp1),u,GeeGooinv,y(:,:,:,:,:,1),idag,coact, &
                 kappa,iflag,bc,vecbl,vecblinv,myid,nn,ldiv,nms,lvbc,ib,lbd, &
                 iblv,MRT)

     mvp = mvp+1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!V_(j+1)=P(A)*A*V_(j)!!!!!!!!!!!!!!!!
!


       do i = 1,j
 
          call vecdot(v(:,:,:,:,:,i),v(:,:,:,:,:,jp1),beta,MRT2)

          ! if (myid==0) then
          !   print *,"after 2.beta", beta(1), beta(2)
          ! endif
 
          hc(:,i,j) = beta(:)
          hc2(:,i,j) = hc(:,i,j)
          hc3(:,i,j) = hc(:,i,j)
 
          do icri = 1,5,2 ! 6=nri*nc
             do jj = 1,nvhalf
                v(icri  ,jj,:,:,:,jp1) = v(icri  ,jj,:,:,:,jp1) &
                                        - beta(1)*v(icri  ,jj,:,:,:,i) &
                                        + beta(2)*v(icri+1,jj,:,:,:,i)
                v(icri+1,jj,:,:,:,jp1) = v(icri+1,jj,:,:,:,jp1) &
                                        - beta(2)*v(icri  ,jj,:,:,:,i) &
                                        - beta(1)*v(icri+1,jj,:,:,:,i)
             enddo ! jj 
          enddo ! icri
       enddo ! i
 
 
       hc(1,j,j) = hc(1,j,j) - sigma(1)
       hc2(:,j,j) = hc(:,j,j)
       hc3(:,j,j) = hc(:,j,j)

       call vecdot(v(:,:,:,:,:,jp1),v(:,:,:,:,:,jp1),beta,MRT2)
 
       hc(1,jp1,j) = sqrt(beta(1))
       hc(2,jp1,j) = 0.0_KR2
       hc2(:,jp1,j) = hc(:,jp1,j)
       hc3(:,jp1,j) = hc(:,jp1,j)

       const = 1.0_KR2/sqrt(beta(1))
       do jj=1,nvhalf
          v(:,jj,:,:,:,jp1) = const*v(:,jj,:,:,:,jp1)
       enddo

       c(:,jp1) = 0.0_KR2
       c2(:,jp1) = c(:,jp1)
 
! DD note ~ I need to find a way of doing Givens rotations
! in the pseudocomplex routines.....

       if (j /= 1) then

          do i = 1,j-1
             tv1(1) = gc(1,i)*hc(1,i,j) - gc(2,i)*hc(2,i,j) &
                    + gs(1,i)*hc(1,i+1,j) + gs(2,i)*hc(2,i+1,j)
             tv1(2) = gc(1,i)*hc(2,i,j) + gc(2,i)*hc(1,i,j) &
                    + gs(1,i)*hc(2,i+1,j) - gs(2,i)*hc(1,i+1,j)
             tv2(1) = gc(1,i)*hc(1,i+1,j) - gc(2,i)*hc(2,i+1,j) &
                    - gs(1,i)*hc(1,i,j) + gs(2,i)*hc(2,i,j)
             tv2(2) = gc(1,i)*hc(2,i+1,j) + gc(2,i)*hc(1,i+1,j) &
                    - gs(1,i)*hc(2,i,j) - gs(2,i)*hc(1,i,j)
             hc(:,i,j) = tv1(:)
             hc(:,i+1,j) = tv2(:)
          enddo ! i

       endif ! (j /= 1)

       amags = hc(1,j,j)**2 + hc(2,j,j)**2
       tv(1) = sqrt( amags + hc(1,j+1,j)**2 + hc(2,j+1,j)**2 )
       tv(2) = 0.0_KR2
       gc(1,j) = sqrt(amags)/tv(1)
       gc(2,j) = 0.0_KR2
       con2 = gc(1,j)/amags
       gs(1,j) = con2*(hc(1,j+1,j)*hc(1,j,j)+hc(2,j+1,j)*hc(2,j,j))
       gs(2,j) = con2*(hc(2,j+1,j)*hc(1,j,j)-hc(1,j+1,j)*hc(2,j,j))
       hc(1,j,j) = gc(1,j)*hc(1,j,j) + gs(1,j)*hc(1,j+1,j) + gs(2,j)*hc(2,j+1,j)
       hc(2,j,j) = gc(1,j)*hc(2,j,j) + gs(1,j)*hc(2,j+1,j) - gs(2,j)*hc(1,j+1,j)
       hc(:,j+1,j) = 0.0_KR2
       tv1(1) = gc(1,j)*c(1,j) + gs(1,j)*c(1,j+1) + gs(2,j)*c(2,j+1)
       tv1(2) = gc(1,j)*c(2,j) + gs(1,j)*c(2,j+1) - gs(2,j)*c(1,j+1)
       tv2(1) = gc(1,j)*c(1,j+1) - gs(1,j)*c(1,j) + gs(2,j)*c(2,j)
       tv2(2) = gc(1,j)*c(2,j+1) - gs(1,j)*c(2,j) - gs(2,j)*c(1,j)
       c(:,j) = tv1(:)
       c(:,j+1) = tv2(:)

       do i = 1,j
          ss(:,i) = c(:,i)
       enddo ! i

      con2 = 1.0_KR2/(hc(1,j,j)**2+hc(2,j,j)**2)
      const = con2*(ss(1,j)*hc(1,j,j)+ss(2,j)*hc(2,j,j))
      ss(2,j) = con2*(ss(2,j)*hc(1,j,j)-ss(1,j)*hc(2,j,j))
      ss(1,j) = const

      if (j/=1) then
         do i = 1,j-1
            ir = j - i + 1
            irm1 = ir - 1
            do jj = 1,irm1
               const = ss(1,jj) - ss(1,ir)*hc(1,jj,ir) + ss(2,ir)*hc(2,jj,ir)
               ss(2,jj) = ss(2,jj) - ss(1,ir)*hc(2,jj,ir) - ss(2,ir)*hc(1,jj,ir)
               ss(1,jj) = const
            enddo ! jj
            con2 = 1.0_KR2/(hc(1,irm1,irm1)**2+hc(2,irm1,irm1)**2)
            const = con2*(ss(1,irm1)*hc(1,irm1,irm1)+ss(2,irm1)*hc(2,irm1,irm1))
            ss(2,irm1) = con2*(ss(2,irm1)*hc(1,irm1,irm1)-ss(1,irm1)*hc(2,irm1,irm1))
            ss(1,irm1) = const
         enddo ! i
      endif ! (j/=1)

      ! ... Define new variable d to assist in shifting masses
      ! d is the "short" solution vector to small problem

      do i=1,j
         d(1,i,1) = ss(1,i)
         d(2,i,1) = ss(2,i)
      enddo ! i

      do i=1,jp1
        st(:,i) = 0.0_KR2
      enddo ! i
  
      do ii = 1,jp1
         do jj = 1,j
            st(1,ii) = st(1,ii) + hc2(1,ii,jj)*d(1,jj,1) - hc2(2,ii,jj)*d(2,jj,1)
            st(2,ii) = st(2,ii) + hc2(1,ii,jj)*d(2,jj,1) + hc2(2,ii,jj)*d(1,jj,1)
         enddo ! jj
         srv(1,ii,1) = c2(1,ii) - st(1,ii)
         srv(2,ii,1) = c2(2,ii) - st(2,ii)
      enddo ! ii

      beta(1) = 0.0_KR

      do jj = 1,j+1
         beta(1) = beta(1) + srv(1,jj,1)**2 + srv(2,jj,1)**2
      enddo ! jj

      ! ... form gdr

      gdr(mvp,1) = sqrt(beta(1))

       ! ... To keep subspaces parallel for different shifts vector
       !     needs to be rotated. Use a QR factorization to do the
       !     ortognal rotation.
       !     BIG Shifting loop


!      if (is==0) then
!
!      do is = 2,nshifts
!
!         do ii=1,jp1
!            do jj=1,j
!               hcs(1,ii,jj) = hc2(1,ii,jj)
!               hcs(2,ii,jj) = hc2(2,ii,jj)
!               hcs2(1,ii,jj) = hc2(1,ii,jj)
!               hcs2(2,ii,jj) = hc2(2,ii,jj)
!            enddo ! jj
!         enddo ! ii
!
!         do jj=1,j
!            hcs(1,jj,jj) = hc2(1,jj,jj) + sigma(1) - sigma(is)
!            hcs2(1,jj,jj) = hc2(1,jj,jj) + sigma(1) - sigma(is) 
!         enddo ! jj
!
!         ! Copy the 12complex arrays into true complex arrays for use with
!         ! lapack routines
!
!         call real2complex_mat(hcs, nmaxGMRES+1, nmaxGMRES+1, zhcs)
!
!         call zgeqrf(j+1,j,zhcs,ldh,ztau,zwork,lzwork,info)
!
!         ! ... store R (upper triangular) in rr
!
!         do ii=1,jp1 
!            do jj=1,j
!               zrr(ii,jj) = zhcs(ii,jj)
!            enddo ! jj
!         enddo ! ii
!
!         call zungqr(j+1,j+1,j,zhcs,ldh,ztau,zwork,lzwork,info)
!
!         ! Copy the complex zhcs array back to hcs 
!
!         call complex2real_mat(zhcs, nmaxGMRES+1, nmaxGMRES+1, hcs)
!
!         ! ... hcs after this call is the qq part of the qr factorization
!
!         ! Now zero out crot (keeps shifts parrallel) and srvrot
!
!          do ii=1,jp1
!             zcrot(ii)=0.0_KR
!             zsrvrot(ii)=0.0_KR
!          enddo ! ii
!
!         do ii=1,jp1
!            do jj=1,jp1
!
!               ztemp1 = DCMPLX(cmult(1,is), cmult(2,is))
!               ztemp2 = DCMPLX(c2(1,jj), c2(2,jj))
!               ztemp3 = DCMPLX(srv(1,jj,1), srv(2,jj,1))
!
!               zcrot(ii) = zcrot(ii) + ztemp1 * CONJG(zhcs(jj,ii)) * ztemp2
!               zsrvrot(ii) = zsrvrot(ii) + CONJG(zhcs(jj,ii)) * ztemp3
!
!            enddo ! jj
!         enddo ! ii
!
!         ! ... construct alpha
!
!           if ((myid == 0) .and. (is==3)) then
!             print *, "jp1, zcrot(jp1), zsrvrot(jp1) = ", jp1, zcrot(jp1), zsrvrot(jp1)
!           endif
!            
!         zalpha = zcrot(jp1)/zsrvrot(jp1)
!
!         alph(1, is) = REAL(zalpha)
!         alph(2, is) = AIMAG(zalpha)
!
!           if (myid==0)  then
!             print *,"is, zalpha, alph = ", is, zalpha, alph(:,is)
!           endif
!
!         do jj=1,j
!            ztemp1 = zalpha * zsrvrot(jj)
!            cmas(1,jj) = REAL(zcrot(jj)) - REAL(ztemp1)    ! alpha*srvrot(1,jj)
!            cmas(2,jj) = AIMAG(zcrot(jj)) - AIMAG(ztemp1)   ! alpha*srvrot(2,jj)
!         enddo ! jj
!
!           if ((myid==0) .and. (is==3))  then
!             print *,"cmas = ", cmas
!           endif
!
!         ! ... solve linear eqns problem d(1:j,is)=rr(1:j,1:j)\cmas
!
!         do ii=1,j
!            d(1,ii,is)=cmas(1,ii)
!            d(2,ii,is)=cmas(2,ii)
!         enddo ! ii
!
!         call real2complex_mat(d,nmaxGMRES,nshifts,zd)
!
!         zd(j,is) = zd(j,is)/zrr(j,j)
!
!         if (j /= 1) then 
!            do i=1,j-1
!               ir = j-i +1
!               irm1 = ir -1
!               call zaxpy(irm1,-zd(ir,is),zrr(1,ir),1,zd(1,is),1)
!               zd(irm1,is) = zd(irm1,is)/zrr(irm1,irm1)
!            enddo ! i
!         endif ! j/=1
!
!         call complex2real_mat(zd, nmaxGMRES, nshifts, d)
!
!         do ii=1,jp1
!            st(1,ii) = 0.0_KR 
!            st(2,ii) = 0.0_KR
!            srvis(1,ii) = 0.0_KR
!            srvis(2,ii) = 0.0_KR
!         enddo ! ii
!
!         do ii=1,jp1
!            do jj=1,j
!               st(1,ii) = st(1,ii) + hcs2(1,ii,jj)*d(1,jj,is) &
!                                   - hcs2(2,ii,jj)*d(2,jj,is)
!               st(2,ii) = st(2,ii) + hcs2(1,ii,jj)*d(2,jj,is) &
!                                   + hcs2(2,ii,jj)*d(1,jj,is)
!            enddo ! jj
!            srvis(1,ii) = cmult(1,is)*c2(1,ii)-cmult(2,is)*c2(2,ii) &
!                        - st(1,ii)
!            srvis(2,ii) = cmult(1,is)*c2(2,ii)+cmult(2,is)*c2(1,ii) &
!                        - st(2,ii)
!         enddo ! ii
!
!         ! ... form the norm of srvis and put in gdr
!
!         beta(1) =0.0_KR
!
!         do jj=1,j+1
!            beta(1) = beta(1) + srvis(1,jj)*srvis(1,jj) &
!                              + srvis(2,jj)*srvis(2,jj)
!         enddo ! jj
!
!         gdr(mvp,is) = sqrt(beta(1))
!
!         if(myid==0.and.is==2) then 
!           print *, "mvp, gdr(mvp,2)=", mvp, gdr(mvp,2)
!         endif ! myid
!
!      enddo ! BIG is loop
!
!      endif ! (is==0)

       if (j>=nDR) then 

        !if(myid==0) then
        !  print *,"j=nDR", j,nDR
        !endif ! myid

! DEAN~ NEED TO HARD CODE IN is==1

!         do is = 1,nshifts
          is=1
             ! cmult(1,is) = alph(1,is) 
             ! cmult(2,is) = alph(2,is) 
             do i=1,j
                sr(i) = d(1,i,is)
                si(i) = d(2,i,is)
             enddo ! i

             do jj=1,nvhalf
               xb(:,jj,:,:,:) = 0.0_KR2
             enddo ! jj

             do jj = 1,j
                do icri = 1,5,2 ! 6=nri*nc
                   do i = 1,nvhalf
                      xb(icri  ,i,:,:,:) = xb(icri  ,i,:,:,:) &
                                         + sr(jj)*v(icri  ,i,:,:,:,jj) &
                                         - si(jj)*v(icri+1,i,:,:,:,jj)
                      xb(icri+1,i,:,:,:) = xb(icri+1,i,:,:,:) &
                                         + si(jj)*v(icri  ,i,:,:,:,jj) &
                                         + sr(jj)*v(icri+1,i,:,:,:,jj)
                   enddo ! i
                enddo ! icri
             enddo ! jj

             do i = 1,nvhalf
                ! HEY!
                xshift(:,i,:,:,:,is) = xshift(:,i,:,:,:,is) + xb(:,i,:,:,:)
             enddo ! i

            !if (myid == 0) then
            ! print *, "idag in proj=", idag
            !endif

             call Hdbletm(h,u,GeeGooinv,xshift(:,:,:,:,:,is),idag,coact,kappa,&
                         iflag,bc,vecbl,vecblinv,myid,nn,ldiv, &
                         nms,lvbc,ib,lbd,iblv,MRT)

             do i=1,nvhalf
                r(:,i,:,:,:,is) = beg(:,i,:,:,:) - h(:,i,:,:,:) &
                                 + sigma(is)*xshift(:,i,:,:,:,is)
             enddo ! i 
           

!       enddo ! is
       endif ! (j>=nDR)
  
       if (j >= nDR) then
          betashift = 0.0_KR2

! DEAN HARD CODE IN is=1

!         do is=1,nshifts
          is=1
             call vecdot(r(:,:,:,:,:,is), r(:,:,:,:,:,is), beta, MRT2)
          
! ~ NOTE ... this should be a sqrt of beta.....

             betashift(1,is) = sqrt(beta(1))
             betashift(2,is) = 0.0_KR2
             rnale(icycle,is) = sqrt(beta(1))
!            cmult(1,is) = alph(1,is)
!            cmult(2,is) = alph(2,is)
!         enddo ! is

! use betashift (1,1) to check convergence of residual and ultimately exit.

         rn = betashift(1,1)

           !if(myid==0) then
           !  print *, "Print before write to LOG"
           !endif ! myid

          ! if (j >= nDR) then 
             if (myid==0) then
                open(unit=8,file=trim(rwdir(myid+1))//"CFGSPROPS.LOG",action="write", &
                     form="formatted",status="old",position="append")
                  write(unit=8,fmt=*) "gmresdrproject-rn",itercount,rn
                  write(unit=8,fmt=*) "gmresdrproject-gdr",itercount,gdr(mvp,1)/rinit!,gdr(mvp,2), gdr(mvp,3), gdr(mvp,4)
                 !print *, betashift(1,1), betashift(1,2), betashift(1,3)
                close(unit=8,status="keep")
             endif
          ! endif ! (j >= nDR)

      !   if (myid ==0) then
      !      open(unit=8,file=trim(rwdir(myid+1))//"GDR.LOG",action="write", &
      !           form="formatted",status="old",position="append")
!     !      write(unit=8,fmt=*) itercount,betashift(1,1),betashift(1,2),betashift(1,3)
      !      write(unit=8,fmt="(i9,a1,es17.10,a1,es17.10,a1,es17.10,a1,es17.10)") itercount," ",&
      !            gdr(mvp,1)/rinit!, " ", gdr(mvp,2), " ", gdr(mvp,3), " ", gdr(mvp,4)
      !      close(unit=8,status="keep")
      !   endif ! myid

          rshift = 0.0_KR2

!         if (myid == 0) then
!            print *,"Checking rshift, v = ", v
!         endif

          do ii=1,nDR+1
             do icri = 1,5,2 ! 6=nri*nc
                do jj=1,nvhalf
                 rshift(icri  ,jj,:,:,:,1) = rshift(icri  ,jj,:,:,:,1) & 
                                            + v(icri,jj,:,:,:,ii)*srv(1,ii,1) &
                                            - v(icri+1,jj,:,:,:,ii)*srv(2,ii,1)
                 rshift(icri+1,jj,:,:,:,1) = rshift(icri+1,jj,:,:,:,1) 
                                             +v(icri,jj,:,:,:,ii)*srv(2,ii,1) &
                                             + v(icri+1,jj,:,:,:,ii)*srv(1,ii,1)
                enddo ! jj
             enddo ! icri
          enddo ! ii

!          rn=gdr(mvp,1)

        !if (myid==0) then
        !   print *, "At bottom of loop PROJ rn/rinit = ", rn, rinit, rn/rinit, "j,mvp=",j,mvp
        !endif
        !if (myid==0) then
        !   print *, "end of gmres "
        !endif

! Don't need this vecdot

         call vecdot(rshift(:,:,:,:,:,1), rshift(:,:,:,:,:,1), beta, MRT2)

! NOTE~ since I took sqrt of betashift above don't need to here,

         const = 1.0_KR/betashift(1,1)

         do jj=1,nvhalf
            v(:,jj,:,:,:,1) = const*rshift(:,jj,:,:,:,1)
         enddo ! jj

         icycle = icycle + 1

      endif ! (j>=nDR)
 
    enddo ! j BIG J LOOP

 enddo cycledo

 !  if (isignal == 1) then 
 !    if (myid == 0) then
 !        print *, "vk+1 gdr ="
 !       do ii=1,mvp 
 !        print *, gdr(ii,:)
 !       enddo ! ii
 !    endif ! myid
 !  else ! isignal
 !    if (myid == 0) then
 !        print *, "RHS gdr ="
 !       do ii=1,mvp 
 !        print *, gdr(ii,:)
 !       enddo ! ii
 !    endif ! myid
 !  endif ! isiganl

  ! if(myid==0) then
  !   print *, "Leaveing PROJ"
  ! endif ! myid


 end subroutine gmresproject

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

