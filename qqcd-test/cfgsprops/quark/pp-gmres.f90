 subroutine PP-gmresdr(rwdir,phi,x,GMRES,resmax,itermin,itercount,u,GeeGooinv, &
                    iflag,kappa,coact,bc,vecbl,vecblinv,myid,nn,ldiv,nms, &
                    lvbc,ib,lbd,iblv,MRT,MRT2)
! GMRES-DR(n,k) with polynomial preconditioning
    use shift

    character(len=*), intent(in),    dimension(:)         :: rwdir
    real(kind=KR),    intent(in),    dimension(:,:,:,:,:) :: phi
    real(kind=KR),    intent(inout), dimension(:,:,:,:,:) :: x
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

    integer(kind=KI) :: icycle, i, j, p, k, jp1, jj, ir, irm1, is, ivb, ii, &
                        idis, idag, icri, nDR, kDR, nrhs, ilo, ihi, ischur, &
                        id, ieo, ibleo!p is the order of the polynomial
    logical,          dimension(nmaxGMRES)                  :: myselect
    integer(kind=KI), dimension(nmaxGMRES)                  :: ipiv
!!!!!!!!!!!!!!!!!!!!!!!!!!solving the coeffcients for the polynomial!!!!!!!!!!!!
    integer(kind=KI), dimension(4)                          :: ipiv2!order of po                                                                     lynomial
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=KR2)                                          :: const, tval, &
                                                               amags, con2, rv, &
                                                               normnum
    real(kind=KR2),   dimension(2)                          :: beta, tv1, tv2, &
                                                               tv
    real(kind=KR2),   dimension(nmaxGMRES)                  :: mag
    real(kind=KR2),   dimension(2,nmaxGMRES)                :: ss, gs, gc, &
                                                               tau, w, work
    real(kind=KR2),   dimension(2,nmaxGMRES+1)              :: c, c2, srv
    real(kind=KR2),   dimension(2,nmaxGMRES,1)              :: em
    real(kind=KR2),   dimension(2,nmaxGMRES+1,nmaxGMRES+1)  :: z,ztmp, q
    real(kind=KR2),   dimension(2,nmaxGMRES,nmaxGMRES+1)    :: hcht,matss
    real(kind=KR2),   dimension(2,nmaxGMRES+1,nmaxGMRES)    :: hc, hc2, hc3, hprint, t, hhh
   
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=KR2),   dimension(2,4,4)                      ::lsmat
    real(kind=KR2),   dimension(2,4,1)                      ::cls
    real(kind=KR2),   dimension(2,4)                        ::co!coefficients
!!!!!!!!!!!parameters in determing the preconditioning polynomial!!!!!!!!!!!!!!!

    real(kind=KR2),   dimension(2,nmaxGMRES+1,kmaxGMRES)    :: ws, ev
    real(kind=KR2),   dimension(2,kmaxGMRES+1,kmaxGMRES)    :: gca, gsa
    real(kind=KR2),   dimension(6,nvhalf,4,2,8)             :: r, h , xt
!!!!!!!!!!!Some intermediate parameters used for the preconditioning!!!!!!!!!!!!
    real(kind=KR2),   dimension(6,nvhalf,4,2,8)             :: w, y , z
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   real(kind=KR2),   dimension(6,ntotal,4,2,8)             :: xb
    !real(kind=KR2),   dimension(6,ntotal,4,2,8,nmaxGMRES+1) :: v
    !real(kind=KR2),   dimension(6,nvhalf,4,2,8,nmaxGMRES+1) :: vtemp
    !real(kind=KR2),   dimension(6,ntotal,4,2,8)             :: xb
    !real(kind=KR2),   dimension(2,nmaxGMRES+1,nmaxGMRES)    :: hcnew
    real(kind=KR2),   dimension(6,nmaxGMRES+1)              :: vt
    real(kind=KR2),   dimension(2,kmaxGMRES,kmaxGMRES)      :: greal
   
    real(kind=KR), dimension(6,nvhalf,4,2,8)   :: htemp, bopart
    real(kind=KR), dimension(6,ntotal,4,2,8,1)  :: getemp
    integer(kind=KI) :: iblock, isite, idirac,icolorir, site, icolorr, irow
    integer(kind=KI) :: didmaindo,roww
    integer(kind=KI), dimension(nmaxGMRES)              :: sortev

    real(kind=KR), dimension(nxyzt,3,4)   :: realz2noise, imagz2noise


! We need to allocate the array v because on 32 bit Linux (IA-32) very large
! lattice sizes (nxyzt) cause the data segment to be too large and the program
! won't run.


! Define some parameters.
!parameters adding for P.P: p,ipiv2,cls,lsmat,co,w,y,z.
    nDR = GMRES(1)
    kDR = GMRES(2)
    p = 4 !the order of the polynomial
    didmaindo = 0
    icycle = 1
    idag = 0
    ss = 0.0_KR
    hc = 0.0_KR
    hc2 = 0.0_KR
    hc3 = 0.0_KR
    hcht = 0.0_KR

 ! htemp = 0.0_KR
 ! site = 0.0_KR
 !do iblock =1,8
 !  do ieo = 1,2
 !    do isite=1,nvhalf
 !      do idirac=1,4
 !        do icolorir=1,5,2
 !               site = ieo + 16*(isite - 1) + 2*(iblock - 1)
 !               icolorr = icolorir/2 +1
 !              !print *, "site,icolorr =", site,icolorr
 !              !print *, "nvhalf, ntotal, nps =", nvhalf, ntotal, nps
!
 !               getemp = 0.0_KR
 !               getemp(icolorir   ,isite,idirac,ieo,iblock,1) = 1.0_KR
 !               getemp(icolorir +1,isite,idirac,ieo,iblock,1) = 0.0_KR
!
!                call Hdbletm(htemp,u,GeeGooinv,getemp(:,:,:,:,:,1),idag,coact,kappa,iflag,bc,vecbl, &
!                                vecblinv,myid,nn,ldiv,nms,lvbc,ib,lbd,iblv,MRT)
!
!                call checkNonZero(htemp(:,:,:,:,:), nvhalf,iblock,ieo,isite,idirac,icolorir,site,icolorr)
!
! To print single rhs source vector use ..
!
!               !irow = icolorr + 3*(isite -1) + 24*(idirac -1) + 96*(ieo -1) + 192*(iblock - 1)
!               !       print *, irow, phi(icolorir,isite,idirac,ieo,iblock), phi(icolorir+1,isite,idirac,ieo,iblock)
!
!             enddo ! icolorir
!          enddo ! idirac
!       enddo ! isite
!    enddo ! ieo
!  enddo ! iblock

!*****MORGAN'S STEP 1: Start.
! Compute r=phi-M*x and v=r/|r| and beta=|r|.
    call Hdbletm(h,u,GeeGooinv,x,idag,coact,kappa,iflag,bc,vecbl,vecblinv, &
                 myid,nn,ldiv,nms,lvbc,ib,lbd,iblv,MRT)
    do i = 1,nvhalf
     r(:,i,:,:,:) = phi(:,i,:,:,:) - h(:,i,:,:,:)
!     v(:,i,:,:,:,1) = r(:,i,:,:,:) we do this after the prconditioning.
     xt(:,i,:,:,:) = x(:,i,:,:,:)
    enddo ! i
!Determine the polynomial.!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    do k = 1,nvhalf
     vprime(:,k,:,:,:,1) = phi(:,k,:,:,:)
    enddo !k
    do i = 1,p
     call Hdbletm(vprime(:,:,:,:,:,i+1),u,GeeGooinv,vprime(:,:,:,:,:,i),idag, &
                 coact,kappa,iflag,bc,vecbl,vecblinv,myid,nn,ldiv,nms,lvbs,ib, &                 lbd,iblv,MRT)
    enddo !i
    
    do i=2,p+1
     do j=2,p+1
      call vecdot(vprime(:,:,:,:,:,i),vprime(:,:,:,:,:,j),beta,MRT2)
      lsmat(:,i-1,j-1) = beta(:)  !lsmat(2,p,p) ,cls(2,p,1)
     enddo!j
    enddo!i

    do i=2,p+1
     call vecdot(vprime(:,:,:,:,:,i),phi(:,:,:,:,:),beta,MRT2)
     cls(:,i-1,1) = beta(:)
    enddo!i
    
    call linearsolver(p,1,lsmat,ipiv2,cls)
    co(:,:) = cls(:,:,1)
    
!!!!
!!!!!Times the polynomial to the residue
    do k=1,nvhalf
     w(:,k,:,:,:) = r(:,k,:,:,:)
    enddo!k
    
    do icri=1,5,2
     do k=1,nvhalf
      y(icri,k,:,:,:) = co(1,1)*w(icri,k,:,:,:)-co(2,1)*w(icri+1,k,:,:,:)
      y(icri+1,k,:,:,:) = co(1,1)*w(icri+1,k,:,:,:)+co(2,1)*w(icri,k,:,:,:)
     enddo!k
    enddo!icri

    do k=1,nvhalf
     z(:,k,:,:,:) = w(:,k,:,:,:)
    enddo!k

    do i=2:p
     call Hdbletm(z(:,:,:,:,:),u,GeeGooinv,z(:,:,:,:,:),idag,coact,kappa, &
                  iflag,bc,vecbl,vecblinv,myid,nn,ldiv,nms,lvbc,ib,ldb,iblv, &
                  MRT )  !z=M*z
     do icri=1,5,2
      do k=1,nvhalf
       y(icri,k,:,:,:) = y(icri,k,:,:,:)+co(1,i)*z(icri,k,:,:,:) &
                         -co(2,i)*z(icri+1,k,:,:,:)
       y(icri+1,k,:,:,:) = y(icri+1,k,:,:,:)+co(1,i)*z(icri+1,k,:,:,:) &
                           +co(2,i)*z(icri,k,:,:,:)   !y=P(A)*r
      enddo!k
     enddo!icri
    enddo!i
    
    do k=1,nvhalf
     v(:,k,:,:,:,1) = y(:,k,:,:,:) ! Use y=P(A)*rto generate V_(m+1)&Hbar_m
    enddo!k
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
    call vecdot(v(:,:,:,:,:,1),v(:,:,:,:,:,1),beta,MRT2)
    beta(1) = sqrt(beta(1))
    normnum = beta(1)

    const = 1.0_KR2/beta(1)
    v(:,:,:,:,:,1) = const*v(:,:,:,:,:,1)
! For use in Morgan's step 2a, define c = beta*e_1.
    c(1,1) = beta(1)
    c(2,1) = 0.0_KR2
    c2(:,1) = c(:,1)

!*****The main loop.
    itercount = 0
    j = 0
    maindo: do




     if ( icycle > kcyclim) exit maindo

     j = j + 1
     jp1 = j + 1
     itercount = itercount + 1

!*****MORGAN'S STEP 2A AND STEPS 7,8A: Apply standard GMRES(nDR).
! Generate V_(m+1) and Hbar_m with the Arnoldi iteration.
     call Hdbletm(v(:,:,:,:,:,jp1),u,GeeGooinv,v(:,:,:,:,:,j),idag,coact, &
                  kappa,iflag,bc,vecbl,vecblinv,myid,nn,ldiv,nms,lvbc,ib,lbd, &
                  iblv,MRT)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!V_(j+1)=P(A)*A*V_(j)!!!!!!!!!!!!!!!!
     do k=1,nvhalf
      w(:,k,:,:,:) = v(:,k,:,:,:,jp1)
     enddo!k

     do icri=1,5,2
      do k=1,nvhalf
       y(icri,k,:,:,:) = co(1,1)*w(icri,k,:,:,:)-co(2,1)*w(icri+1,k,:,:,:)
       y(icri+1,k,:,:,:) = co(1,1)*w(icri+1,k,:,:,:)+co(2,1)*w(icri,k,:,:,:)
      enddo!k
     enddo!icri

     do k=1,nvhalf
      z(:,k,:,:,:) = w(:,k,:,:,:)
     enddo!k

     do i=2,p
      call Hdbletm(z(:,:,:,:,:),u,GeeGooinv,z(:,:,:,:,:),idag,coact,kappa, &
                   iflag,bc,vecbl,vecblinv,myid,nn,ldiv,nms,lvbc,ib,lbd,iblv, &
                   MRT)
      do icri = 1,5,2
       do k=1,nvhalf
        y(icri,k,:,:,:) = y(icri,k,:,:,:)+co(1,i)*z(icri,k,:,:,:) &
                          -co(2,i)*z(icri+1,k,:,:,:)
        y(icri+1,k,:,:,:) = y(icri+1,k,:,:,:)+co(1,i)*z(icri+1,k,:,:,:) &
                            +co(2,i)*z(icri,k,:,:,:) !y=P(A)*A*V_(j)
       enddo!k
      enddo!icri
     enddo!i
     do k=1,nvhalf
      v(:,k,:,:,:,jp1)=y(:,k,:,:,:)
     enddo!k

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     do i = 1,j
      call vecdot(v(:,:,:,:,:,i),v(:,:,:,:,:,jp1),beta,MRT2)
      hc(:,i,j) = beta(:)

!     print *, "i,j, hc(:,i,j)=", i,j, hc(:,i,j)

      hc2(:,i,j) = hc(:,i,j)
      hc3(:,i,j) = hc(:,i,j)
      hcht(1,j,i) = hc(1,i,j)
      hcht(2,j,i) = -hc(2,i,j)
      do icri = 1,5,2 ! 6=nri*nc
       do k = 1,nvhalf
        v(icri  ,k,:,:,:,jp1) = v(icri  ,k,:,:,:,jp1) &
                              - beta(1)*v(icri  ,k,:,:,:,i) &
                              + beta(2)*v(icri+1,k,:,:,:,i)
        v(icri+1,k,:,:,:,jp1) = v(icri+1,k,:,:,:,jp1) &
                              - beta(2)*v(icri  ,k,:,:,:,i) &
                              - beta(1)*v(icri+1,k,:,:,:,i)
       enddo ! k
      enddo ! icri
     enddo ! i



     call vecdot(v(:,:,:,:,:,jp1),v(:,:,:,:,:,jp1),beta,MRT2)
     hc(1,jp1,j) = sqrt(beta(1))
     hc(2,jp1,j) = 0.0_KR2
     hc2(:,jp1,j) = hc(:,jp1,j)
     hc3(:,jp1,j) = hc(:,jp1,j)
     hcht(1,j,jp1) = hc(1,jp1,j)
     hcht(2,j,jp1) = -hc(2,jp1,j)
     const = 1.0_KR2/sqrt(beta(1))
     v(:,:,:,:,:,jp1) = const*v(:,:,:,:,:,jp1)
     c(:,jp1) = 0.0_KR2
     c2(:,jp1) = c(:,jp1)


! Solve min|c-Hbar*ss| for ss, where c=beta*e_1.
     if (icycle/=1) then
      do jj = 1,kDR
       do i = jj+1,kDR+1
        tv1(1) = gca(1,i,jj)*hc(1,jj,j) - gca(2,i,jj)*hc(2,jj,j) &
               + gsa(1,i,jj)*hc(1,i,j) + gsa(2,i,jj)*hc(2,i,j)
        tv1(2) = gca(1,i,jj)*hc(2,jj,j) + gca(2,i,jj)*hc(1,jj,j) &
               + gsa(1,i,jj)*hc(2,i,j) - gsa(2,i,jj)*hc(1,i,j)
        tv2(1) = gca(1,i,jj)*hc(1,i,j) - gca(2,i,jj)*hc(2,i,j) &
               - gsa(1,i,jj)*hc(1,jj,j) + gsa(2,i,jj)*hc(2,jj,j)
        tv2(2) = gca(1,i,jj)*hc(2,i,j) + gca(2,i,jj)*hc(1,i,j) &
               - gsa(1,i,jj)*hc(2,jj,j) - gsa(2,i,jj)*hc(1,jj,j)
        hc(:,jj,j) = tv1(:)
        hc(:,i,j) = tv2(:)
       enddo ! i
      enddo ! jj
      if (j>kDR+1) then
       do i = kDR+1,j-1
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
      endif
     elseif (j/=1) then
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
     endif
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
     con2 = 1.0_KR2/(hc(1,j,j)**2+hc(2,j,j)**2
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
     endif
! Form the approximate new solution x = xt + V*ss.
     xb = 0.0_KR2
     do jj = 1,j
      do icri = 1,5,2 ! 6=nri*nc
       do i = 1,nvhalf
        xb(icri  ,i,:,:,:) = xb(icri  ,i,:,:,:) + ss(1,jj)*v(icri  ,i,:,:,:,jj)&
                                                - ss(2,jj)*v(icri+1,i,:,:,:,jj)
        xb(icri+1,i,:,:,:) = xb(icri+1,i,:,:,:) + ss(2,jj)*v(icri  ,i,:,:,:,jj)&
                                                + ss(1,jj)*v(icri+1,i,:,:,:,jj)
       enddo ! i
      enddo ! icri
     enddo ! jj
     do i = 1,nvhalf
      x(:,i,:,:,:) = xt(:,i,:,:,:) + xb(:,i,:,:,:)
     enddo ! i

! Define a small residual vector, srv = c-Hbar*ss, which corresponds to the
! kDR+1 column of the new V that will be formed.
     do i = 1,nDR+1
      srv(:,i) = c2(:,i)
     enddo ! i
     do jj = 1,nDR
      do i = 1,nDR+1
       srv(1,i) = srv(1,i) - ss(1,jj)*hc3(1,i,jj) + ss(2,jj)*hc3(2,i,jj)
       srv(2,i) = srv(2,i) - ss(1,jj)*hc3(2,i,jj) - ss(2,jj)*hc3(1,i,jj)
      enddo ! i
     enddo ! jj

!*****Only deflate after V_(m+1) and Hbar_m have been fully formed.
     if (j>=nDR) then

!*****MORGAN'S STEP 2B AND STEP 8B: Let xt=x and r=phi-M*x.

      call Hdbletm(h,u,GeeGooinv,xb,idag,coact,kappa,iflag,bc,vecbl,vecblinv, &
                   myid,nn,ldiv,nms,lvbc,ib,lbd,iblv,MRT)

      do i = 1,nvhalf
       r(:,i,:,:,:) = r(:,i,:,:,:) - h(:,i,:,:,:)
      enddo ! i

      beta =0.0_KR
      call vecdot(r,r,beta,MRT2)
      beta(1) = sqrt(beta(1))
      if (myid==0) then
       open(unit=8,file=trim(rwdir(myid+1))//"CFGSPROPS.LOG",action="write", &
            form="formatted",status="old",position="append")
!        write(unit=8,fmt="(a12,i9,es17.10)") "gmresdr",itercount,beta(1)
        write(unit=8,fmt="(a12,i9,es17.10)") "gmresdr-norm",itercount,beta(1)/normnum
       close(unit=8,status="keep")
      endif

      if ((beta(1)/normnum)<resmax .and. itercount>=itermin) then
         if(didmaindo==0) then
         print *, "Everyone needs more resmax!"
         endif ! didmaindo
         exit maindo
      endif ! beta
      didmaindo=1
      do i = 1,nvhalf
       xt(:,i,:,:,:) = x(:,i,:,:,:)
      enddo ! i
!*****MORGAN'S STEP 2C AND STEP 9: Compute the kDR smallest eigenpairs of
!                                  H + beta^2*H^(-dagger)*e_(nDR)*e_(nDR)^T.
! These eigenvalues are the harmonic Ritz values, and they are approximate
! eigenvalues for the large matrix.  (Not always accurate approximations.)
      do i = 1,nDR
       em(:,i,1) = 0.0_KR
      enddo ! i
      em(1,nDR,1) = 1.0_KR2
      nrhs = 1





      call linearsolver(nDR,nrhs,hcht,ipiv,em)


      do i = 1,nDR
       hc2(1,i,nDR) = hc2(1,i,nDR) + em(1,i,1)*hc2(1,nDR+1,nDR)**2 &
                    - em(1,i,1)*hc2(2,nDR+1,nDR)**2 &
                    - 2.0_KR2*em(2,i,1)*hc2(1,nDR+1,nDR)*hc2(2,nDR+1,nDR)
       hc2(2,i,nDR) = hc2(2,i,nDR) + em(2,i,1)*hc2(1,nDR+1,nDR)**2 &
                    + 2.0_KR2*em(1,i,1)*hc2(2,nDR+1,nDR)*hc2(1,nDR+1,nDR) &
                    - em(2,i,1)*hc2(2,nDR+1,nDR)**2
      enddo ! i

      ilo = 1
      ihi = nDR
      call hessenberg(hc2,nDR,ilo,ihi,tau,work)




      do jj = 1,nDR
       do i = jj,nDR
        z(:,i,jj) = hc2(:,i,jj)
       enddo ! i
      enddo ! jj
      call qgenerator(nDR,ilo,ihi,z,tau,work)
      ischur = 1
      call evalues(ischur,nDR,ilo,ihi,hc2,w,z,work)

      !   if (myid==0) then
      !      print *, "------------"
      !      print *, "evaluesgmresdr:", w
      !      print *, "------------"
      !   endif 

!*****MORGAN'S STEP 3: Orthonormalization of the first kDR vectors.
! Instead of using eigenvectors, the Schur vectors will be used
! -- see the sentence following equation 3.1 of Morgan -- 
! so reorder the harmonic Ritz values and rearrange the Schur form.
      mag(1) = sqrt(w(1,1)**2+w(2,1)**2)
      do i = 2,nDR
       mag(i) = sqrt(w(1,i)**2+w(2,i)**2)
       is = 0
       ritzloop: do
        is = is + 1
        if (is>i-1) exit ritzloop
        if (mag(i)<mag(is)) then
         tval = mag(i)
         do ivb = i-1,is,-1
          mag(ivb+1) = mag(ivb)
         enddo ! ivb
         mag(is) = tval
         exit ritzloop
        endif
       enddo ritzloop
      enddo ! i
      do i = 1,nDR
       myselect(i) = .false.
       if (sqrt(w(1,i)**2+w(2,i)**2)<=mag(kDR)) myselect(i)=.true.
      enddo ! i

      call orgschur(myselect,nDR,hc2,z,w,idis)
!*****MORGAN'S STEP 4: Orthonormalization of the kDR+1 vector.
! Orthonormalize the vector srv against the first kDR columns of z to form the
! kDR+1 column of z.
      do i = 1,kDR
       z(:,nDR+1,i) = 0.0_KR2
      enddo ! i
      do i = 1,nDR+1
       z(:,i,kDR+1) = srv(:,i)
      enddo ! i
      do jj = 1,kDR
       tv = 0.0_KR2
       do i = 1,nDR+1
        tv(1) = tv(1) + z(1,i,jj)*z(1,i,kDR+1) + z(2,i,jj)*z(2,i,kDR+1)
        tv(2) = tv(2) + z(1,i,jj)*z(2,i,kDR+1) - z(2,i,jj)*z(1,i,kDR+1)
       enddo ! i
       do i = 1,nDR+1
        z(1,i,kDR+1) = z(1,i,kDR+1) - tv(1)*z(1,i,jj) + tv(2)*z(2,i,jj)
        z(2,i,kDR+1) = z(2,i,kDR+1) - tv(1)*z(2,i,jj) - tv(2)*z(1,i,jj)
       enddo ! i
      enddo ! jj
      jj = nDR + 1
      rv = twonorm(z(:,:,kDR+1),jj)
      rv = 1.0_KR2/rv
      do i = 1,nDR+1
       z(:,i,kDR+1) = rv*z(:,i,kDR+1)
      enddo ! i

!*****MORGAN'S STEP 5: Form portions of the new H and V using the old H and V.
      do jj = 1,kDR
       do ii = 1,nDR+1
        ws(:,ii,jj) = 0.0_KR2
       enddo ! ii
       do ii = 1,nDR
        do i = 1,nDR+1
         ws(1,i,jj) = ws(1,i,jj) + z(1,ii,jj)*hc3(1,i,ii) &
                                 - z(2,ii,jj)*hc3(2,i,ii)
         ws(2,i,jj) = ws(2,i,jj) + z(1,ii,jj)*hc3(2,i,ii) &
                                 + z(2,ii,jj)*hc3(1,i,ii)
        enddo ! i
       enddo ! ii
      enddo ! jj
      do jj = 1,kDR
       do ii = 1,kDR+1
        hcnew(:,ii,jj) = 0.0_KR2
        do i = 1,nDR+1
         hcnew(1,ii,jj) = hcnew(1,ii,jj) + z(1,i,ii)*ws(1,i,jj) &
                                         + z(2,i,ii)*ws(2,i,jj)
         hcnew(2,ii,jj) = hcnew(2,ii,jj) + z(1,i,ii)*ws(2,i,jj) &
                                         - z(2,i,ii)*ws(1,i,jj)
        enddo ! i
       enddo ! ii
      enddo ! jj

      do jj = 1,nDR
       do ii = 1,nDR
        hcht(:,ii,jj) = 0.0_KR2
        hc2(:,ii,jj) = 0.0_KR2
       enddo ! ii
      enddo ! jj
      do jj = 1,kDR
       do ii = 1,kDR+1
        hc(:,ii,jj) = hcnew(:,ii,jj)
        hc2(:,ii,jj) = hcnew(:,ii,jj)
        hc3(:,ii,jj) = hcnew(:,ii,jj)
       enddo ! ii
       do ii = 1,kDR+1
        hcht(1,jj,ii) = hcnew(1,ii,jj)
        hcht(2,jj,ii) = -hcnew(2,ii,jj)
       enddo ! ii
      enddo ! jj
      do ii = 1,kDR+1
       c(:,ii) = 0.0_KR2
       do i = 1,nDR+1
        c(1,ii) = c(1,ii) + z(1,i,ii)*srv(1,i) + z(2,i,ii)*srv(2,i)
        c(2,ii) = c(2,ii) + z(1,i,ii)*srv(2,i) - z(2,i,ii)*srv(1,i)
       enddo ! i
       c2(:,ii) = c(:,ii)
      enddo ! ii

      vtemp = 0.0_KR

      do ibleo = 1,8
       do ieo = 1,2
        do id = 1,4
         do i = 1,nvhalf
          do jj = 1,kDR+1
           vt(:,jj) = 0.0_KR2
           do k = 1,nDR+1
            do icri = 1,5,2 ! 6=nri*nc
             vt(icri  ,jj) = vt(icri  ,jj) &
                           + z(1,k,jj)*v(icri  ,i,id,ieo,ibleo,k) &
                           - z(2,k,jj)*v(icri+1,i,id,ieo,ibleo,k)
             vt(icri+1,jj) = vt(icri+1,jj) &
                           + z(2,k,jj)*v(icri  ,i,id,ieo,ibleo,k) &
                           + z(1,k,jj)*v(icri+1,i,id,ieo,ibleo,k)
            enddo ! icri
           enddo ! k
          enddo ! jj
          do jj = 1,kDR+1
           v(:,i,id,ieo,ibleo,jj) = vt(:,jj)
           vtemp(:,i,id,ieo,ibleo,jj) = vt(:,jj)
          enddo ! jj
         enddo ! i
        enddo ! id
       enddo ! ieo
      enddo ! ibleo
         !if(myid==0) then
         ! print *, "vtemp in this shizat-only 1:kDR", vtemp
         !endif ! myid
!*****MORGAN'S STEP 6: Reorthogonalization of k+1 vector.
      do jj = 1,kDR
       call vecdot(v(:,:,:,:,:,jj),v(:,:,:,:,:,kDR+1),beta,MRT2)
       do icri = 1,5,2 ! 6=nri*nc
        do i = 1,nvhalf
         v(icri  ,i,:,:,:,kDR+1) = v(icri  ,i,:,:,:,kDR+1) &
                                 - beta(1)*v(icri  ,i,:,:,:,jj) &
                                 + beta(2)*v(icri+1,i,:,:,:,jj)
         v(icri+1,i,:,:,:,kDR+1) = v(icri+1,i,:,:,:,kDR+1) &
                                 - beta(2)*v(icri  ,i,:,:,:,jj) &
                                 - beta(1)*v(icri+1,i,:,:,:,jj)
        enddo ! i
       enddo ! icri
      enddo ! jj
      call vecdot(v(:,:,:,:,:,kDR+1),v(:,:,:,:,:,kDR+1),beta,MRT2)
      const = 1.0_KR2/sqrt(beta(1))
      do i = 1,nvhalf
       v(:,i,:,:,:,kDR+1)     = const*v(:,i,:,:,:,kDR+1)
      enddo ! i

! Need to have the vtemp vector for the gmresproj routine....

      do jj = 1,nvhalf
       vtemp(:,jj,:,:,:,kDR+1) = v(:,jj,:,:,:,kDR+1)
      enddo ! jj

! Rotations for newly formed hc() matrix.
      do jj = 1,kDR
       do i = jj+1,kDR+1
        amags = hc(1,jj,jj)**2 + hc(2,jj,jj)**2
        con2 = 1.0_KR2/amags
        tv(1) = sqrt(amags+hc(1,i,jj)**2+hc(2,i,jj)**2)
        tv(2) = 0.0_KR2
        gca(1,i,jj) = sqrt(amags)/tv(1)
        gca(2,i,jj) = 0.0_KR2
        gsa(1,i,jj) = gca(1,i,jj)*con2 &
                      *(hc(1,i,jj)*hc(1,jj,jj)+hc(2,i,jj)*hc(2,jj,jj))
        gsa(2,i,jj) = gca(1,i,jj)*con2 &
                      *(hc(2,i,jj)*hc(1,jj,jj)-hc(1,i,jj)*hc(2,jj,jj))
        do j = jj,kDR
         tv1(1) = gca(1,i,jj)*hc(1,jj,j) + gsa(1,i,jj)*hc(1,i,j) &
                                         + gsa(2,i,jj)*hc(2,i,j)
         tv1(2) = gca(1,i,jj)*hc(2,jj,j) + gsa(1,i,jj)*hc(2,i,j) &
                                         - gsa(2,i,jj)*hc(1,i,j)
         tv2(1) = gca(1,i,jj)*hc(1,i,j) - gsa(1,i,jj)*hc(1,jj,j) &
                                        + gsa(2,i,jj)*hc(2,jj,j)
         tv2(2) = gca(1,i,jj)*hc(2,i,j) - gsa(1,i,jj)*hc(2,jj,j) &
                                        - gsa(2,i,jj)*hc(1,jj,j)
         hc(:,jj,j) = tv1(:)
         hc(:,i,j) = tv2(:)
        enddo ! j
        tv1(1) = gca(1,i,jj)*c(1,jj) + gsa(1,i,jj)*c(1,i) + gsa(2,i,jj)*c(2,i)
        tv1(2) = gca(1,i,jj)*c(2,jj) + gsa(1,i,jj)*c(2,i) - gsa(2,i,jj)*c(1,i)
        tv2(1) = gca(1,i,jj)*c(1,i) - gsa(1,i,jj)*c(1,jj) + gsa(2,i,jj)*c(2,jj)
        tv2(2) = gca(1,i,jj)*c(2,i) - gsa(1,i,jj)*c(2,jj) - gsa(2,i,jj)*c(1,jj)
        c(:,jj) = tv1(:)
        c(:,i) = tv2(:)
       enddo ! i
      enddo ! jj
      j = kDR
      icycle = icycle + 1
     endif
    enddo maindo

     !     if (myid==0) then
     !        print *, "evaluesgmresdr:-Dean wuz here"
     !        do ii =1,kDR
     !          print *, w(:,ii)
     !        enddo ! ii
     !     endif 
 end subroutine PP-gmresdr

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
