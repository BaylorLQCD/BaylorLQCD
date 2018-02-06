! DEAN -(7-31-05) Tomarrow I need to set up the reads of the tmU and 
!                 tmD files to use in the nucleon construction. I also
!                 need to change the looping structure since the masses
!                 are looped on the outside of the subroutine in 
!                 generator.



! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
! quark.f90, Randy Lewis, randy.lewis@uregina.ca
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
! This module computes fermion propagators on existing configurations.
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 module quark

!    use MPI
    use kinds
    use latdims
    use basics
    use lattice
    use gaugetools
    use diracops
    use inverters
    use gmresrhs

    implicit none
    private

! Use the following line if the MPI module is not available.
   include 'mpif.h'

! Define access to subroutines.
    public  :: pointsource, SSTsource, fermprop, pionrho, SSTpion, mesoncorr, &
               corr, nucleoncorr, FPwrite, FPread, gamma5vector, fermpropR,&
               fermprop_wilson
    private :: makeG, bsource, xvector, check4ones, printlog

 contains

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 subroutine pointsource(be,bo,src,icsrc,idsrc,myid,iblv,vecblinv)
! Define the source vector to be used for propagator construction.
! INPUT:
!   src is the location of the point source on the global lattice.

    real(kind=KR),    intent(out), dimension(:,:,:,:,:) :: be, bo
    integer(kind=KI), intent(in),  dimension(:)         :: src
    integer(kind=KI), intent(in)                        :: icsrc, idsrc, myid
    integer(kind=KI), intent(in),  dimension(:,:)       :: iblv, vecblinv

    integer(kind=KI), dimension(4) :: np, ip, srcbit
    integer(kind=KI) :: intx, inty, intz, intt, isite, ieo, ibl, mu, ibleo

! Some initializations.
    np(1) = npx
    np(2) = npy
    np(3) = npz
    np(4) = npt
    call atoc(myid,np,ip)

! Initialize the source vector.
    be = 0.0_KR
    bo = 0.0_KR

! srcbit is the position of the source on this sublattice.
    srcbit(1) = src(1) - ip(1)*nx/npx
    srcbit(2) = src(2) - ip(2)*ny/npy
    srcbit(3) = src(3) - ip(3)*nz/npz
    srcbit(4) = src(4) - ip(4)*nt/npt

! Entries are left at zero unless the point source is on this process.
    if (srcbit(1)>0.and.srcbit(1)<nx/npx+1.and. &
        srcbit(2)>0.and.srcbit(2)<ny/npy+1.and. &
        srcbit(3)>0.and.srcbit(3)<nz/npz+1.and. &
        srcbit(4)>0.and.srcbit(4)<nt/npt+1) then
!-determine the isite of the point source.
     intx = (srcbit(1)-1)/4
     inty = (srcbit(2)-1)/2
     intz = (srcbit(3)-1)/2
     intt = (srcbit(4)-1)/2
     isite = 1 + intx + inty*nx/(4*npx) + intz*nx*ny/(8*npx*npy) &
           + intt*nx*ny*nz/(16*npx*npy*npz)
!-determine the ieo and ibl of the point source.
     ieo = 1
     ibl = 1
     do mu = 1,4
      if (modulo(srcbit(mu),4)==3.or.modulo(srcbit(mu),4)==0) ieo = 3 - ieo
      if (modulo(srcbit(mu),2)==0) ibl = iblv(ibl,mu)
     enddo ! mu
!-determine the ibleo of the point source, and set that entry to unity.
     ibleo = vecblinv(2,ibl)
     if (vecblinv(1,ibl)==1) then
!     print "(a34,6i3)", "myid,be(ic,isite,id,ieo,ibleo)=",&
!                         myid,2*icsrc-1,isite,idsrc,ieo,ibleo
      be(2*icsrc-1,isite,idsrc,ieo,ibleo) = 1.0_KR
     else
!     print "(a34,6i3)", "myid,bo(ic,isite,id,ieo,ibleo)=",&
!                         myid,2*icsrc-1,isite,idsrc,ieo,ibleo
      bo(2*icsrc-1,isite,idsrc,ieo,ibleo) = 1.0_KR
     endif
    endif

 end subroutine pointsource

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 subroutine SSTsource(xe,xo,src,opSST,itgbl,pSST,be,bo,bc,myid,iblv,vecblinv)
! In order to insert the current between two propagators, this subroutine
! appends an operator to an existing propagator.
! INPUT:
!   xe() contains one column of the initial propagator, on even lattice sites
!        with respect to the entire lattice (not blocked).
!   xo() contains one column of the initial propagator, on odd lattice sites
!        with respect to the entire lattice (not blocked).
!   src is the location of the first propagator's source on the global lattice.
!   opSST=1 and opSST=-1 are equivalent.  They denote a pseudoscalar operator.
!           (opSST=1 is for a local pseudoscalar and opSST=-1 for a smeared one,
!           but the smearing is NOT done within this subroutine.)
!   tgbl is the timeslice of the inserted operator on the global lattice.
!   pSST() is the 3-momentum of the inserted operator.
! OUTPUT:
!   be() contains one column of the initial propagator with the operator
!        appended, on even lattice sites with respect to the entire lattice
!        (not blocked).
!   bo() contains one column of the initial propagator with the operator
!        appended, on odd lattice sites with respect to the entire lattice
!        (not blocked).

    real(kind=KR),    intent(in),  dimension(:,:,:,:,:) :: xe, xo
    integer(kind=KI), intent(in),  dimension(:)         :: src, pSST, bc
    integer(kind=KI), intent(in)                        :: opSST, itgbl, myid
    real(kind=KR),    intent(out), dimension(:,:,:,:,:) :: be, bo
    integer(kind=KI), intent(in),  dimension(:,:)       :: iblv, vecblinv

    integer(kind=KI), dimension(4) :: np, ip, ierr
    integer(kind=KI)               :: ixgbl, iygbl, izgbl, ixloc, iyloc, &
                                      izloc, itloc, isite, ieo, ibl, ibleo, &
                                      ir, ixbit, iybit, izbit, itbit, lessnoise
    real(kind=KR)                  :: cospx, cospy, cospz, sinpx, sinpy, &
                                      sinpz, cospr, sinpr
    real(kind=KR2),   parameter    :: twoPi=6.283185307179586_KR2

! Initialize the output.
    be = 0.0_KR
    bo = 0.0_KR

! NOTICE: If the boundary conditions are periodic in ALL spatial directions,
!         then the terms which are odd in sinx, siny, sinz, sin2x, sin2y or
!         sin2z are omitted, since the integral of an even periodic function,
!         times sinx or sin2x etc., vanishes.  This should reduce the
!         statistical noise.
! WE ARE ASSUMING THAT THESE SST PROPAGATORS WILL BE USED IN 3-POINT FUNCTIONS.
     if (bc(1)==1.and.bc(2)==1.and.bc(3)==1) then
      lessnoise = 1
     else
      lessnoise = 0
     endif

    select case(opSST)
     case(1,-1) ! local or smeared pseudoscalar operator
      np(1) = npx
      np(2) = npy
      np(3) = npz
      np(4) = npt
      call atoc(myid,np,ip)
      itbit = itgbl - ip(4)*nt/npt
!-only calculate if the relevant timeslice is on this sublattice.
      if (itbit>0 .and. itbit<nt/npt+1) then
       itloc = (itbit-1)/2
       do izbit = 1,nz/npz
        izloc = (izbit-1)/2
        izgbl = izbit + ip(3)*nz/npz
        do iybit = 1,ny/npy
         iyloc = (iybit-1)/2
         iygbl = iybit + ip(2)*ny/npy
         do ixbit = 1,nx/npx
          ixloc = (ixbit-1)/4
          ixgbl = ixbit + ip(1)*nx/npx
!-determine isite,ieo,ibl,ibleo for this site.
          isite = 1 + ixloc + iyloc*nx/(4*npx) + izloc*nx*ny/(8*npx*npy) &
                + itloc*nx*ny*nz/(16*npx*npy*npz)
          ieo = 1
          ibl = 1
          if (modulo(ixbit,4)==3.or.modulo(ixbit,4)==0) ieo = 3 - ieo
          if (modulo(iybit,4)==3.or.modulo(iybit,4)==0) ieo = 3 - ieo
          if (modulo(izbit,4)==3.or.modulo(izbit,4)==0) ieo = 3 - ieo
          if (modulo(itbit,4)==3.or.modulo(itbit,4)==0) ieo = 3 - ieo
          if (modulo(ixbit,2)==0) ibl = iblv(ibl,1)
          if (modulo(iybit,2)==0) ibl = iblv(ibl,2)
          if (modulo(izbit,2)==0) ibl = iblv(ibl,3)
          if (modulo(itbit,2)==0) ibl = iblv(ibl,4)
          ibleo = vecblinv(2,ibl)
!-at this site, be,bo = gamma5 * exp(i*p.r) * xe,xo.
          cospx = cos(twoPi*real(pSST(1)*(ixgbl-src(1)),KR)/real(nx,KR))
          sinpx = sin(twoPi*real(pSST(1)*(ixgbl-src(1)),KR)/real(nx,KR))
          cospy = cos(twoPi*real(pSST(2)*(iygbl-src(2)),KR)/real(ny,KR))
          sinpy = sin(twoPi*real(pSST(2)*(iygbl-src(2)),KR)/real(ny,KR))
          cospz = cos(twoPi*real(pSST(3)*(izgbl-src(3)),KR)/real(nz,KR))
          sinpz = sin(twoPi*real(pSST(3)*(izgbl-src(3)),KR)/real(nz,KR))
          if (lessnoise==0) then
           cospr = cospx*cospy*cospz - cospx*sinpy*sinpz &
                 - sinpx*cospy*cospz - sinpx*sinpy*cospz
           sinpr = sinpx*cospy*cospz - sinpx*sinpy*sinpz &
                 + cospx*sinpy*cospz + cospx*cospy*sinpz
          else
           cospr = cospx*cospy*cospz
           sinpr = 0.0_KR
          endif
          if (vecblinv(1,ibl)==1) then
           do ir = 1,5,2
            be(ir  ,isite,1,ieo,ibleo) = cospr*xe(ir+1,isite,3,ieo,ibleo) &
                                       + sinpr*xe(ir  ,isite,3,ieo,ibleo)
            be(ir+1,isite,1,ieo,ibleo) = sinpr*xe(ir+1,isite,3,ieo,ibleo) &
                                       - cospr*xe(ir  ,isite,3,ieo,ibleo)
            be(ir  ,isite,2,ieo,ibleo) = cospr*xe(ir+1,isite,4,ieo,ibleo) &
                                       + sinpr*xe(ir  ,isite,4,ieo,ibleo)
            be(ir+1,isite,2,ieo,ibleo) = sinpr*xe(ir+1,isite,4,ieo,ibleo) &
                                       - cospr*xe(ir  ,isite,4,ieo,ibleo)
            be(ir  ,isite,3,ieo,ibleo) =-cospr*xe(ir+1,isite,1,ieo,ibleo) &
                                       - sinpr*xe(ir  ,isite,1,ieo,ibleo)
            be(ir+1,isite,3,ieo,ibleo) = cospr*xe(ir  ,isite,1,ieo,ibleo) &
                                       - sinpr*xe(ir+1,isite,1,ieo,ibleo)
            be(ir  ,isite,4,ieo,ibleo) =-cospr*xe(ir+1,isite,2,ieo,ibleo) &
                                       - sinpr*xe(ir  ,isite,2,ieo,ibleo)
            be(ir+1,isite,4,ieo,ibleo) = cospr*xe(ir  ,isite,2,ieo,ibleo) &
                                       - sinpr*xe(ir+1,isite,2,ieo,ibleo)
           enddo ! ir
          else
           do ir = 1,5,2
            bo(ir  ,isite,1,ieo,ibleo) = cospr*xo(ir+1,isite,3,ieo,ibleo) &
                                       + sinpr*xo(ir  ,isite,3,ieo,ibleo)
            bo(ir+1,isite,1,ieo,ibleo) = sinpr*xo(ir+1,isite,3,ieo,ibleo) &
                                       - cospr*xo(ir  ,isite,3,ieo,ibleo)
            bo(ir  ,isite,2,ieo,ibleo) = cospr*xo(ir+1,isite,4,ieo,ibleo) &
                                       + sinpr*xo(ir  ,isite,4,ieo,ibleo)
            bo(ir+1,isite,2,ieo,ibleo) = sinpr*xo(ir+1,isite,4,ieo,ibleo) &
                                       - cospr*xo(ir  ,isite,4,ieo,ibleo)
            bo(ir  ,isite,3,ieo,ibleo) =-cospr*xo(ir+1,isite,1,ieo,ibleo) &
                                       - sinpr*xo(ir  ,isite,1,ieo,ibleo)
            bo(ir+1,isite,3,ieo,ibleo) = cospr*xo(ir  ,isite,1,ieo,ibleo) &
                                       - sinpr*xo(ir+1,isite,1,ieo,ibleo)
            bo(ir  ,isite,4,ieo,ibleo) =-cospr*xo(ir+1,isite,2,ieo,ibleo) &
                                       - sinpr*xo(ir  ,isite,2,ieo,ibleo)
            bo(ir+1,isite,4,ieo,ibleo) = cospr*xo(ir  ,isite,2,ieo,ibleo) &
                                       - sinpr*xo(ir+1,isite,2,ieo,ibleo)
           enddo ! ir
          endif
         enddo ! ixbit
        enddo ! iybit
       enddo ! izbit
      endif
     case default
      open(unit=8,file="QUARK.ERROR",action="write",status="replace", &
           form="formatted")
       write(unit=8,fmt=*) "subroutine SSTsource: opSST =", opSST
      close(unit=8,status="keep")
      print *, "Calling MPI_ABORT near 'subroutine SSTsource: opSST'"
      call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
      stop
    end select

 end subroutine SSTsource

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
! This code was originally designed to do serial masses (or kappas)
! for the Wilson case, however, it is presently done by treating
! Wilson as a special case of twisted with mu=0, and using fermprop.
!
 subroutine fermprop_wilson(rwdir,be,bo,kappa,inverter,coact,bc,resmax, &
                     itermin,u,mGMRES,kGMRES,xe,xo,myid,nn,ldiv, &
                     nms,lvbc,lvcl,ib,lbd,iblv,vecbl,vecblinv,&
                     hoption,MRT,MRT2)
! Compute fermion propagators for wilson fermions with no multimass.
! INPUT:
!   be() contains the entries of the input Dirac vector that are on even
!        lattice sites, with respect to the entire lattice (not blocked).
!   bo() contains the entries of the input Dirac vector that are on odd
!        lattice sites, with respect to the entire lattice (not blocked).
!   kappa is the hopping parameter
!   inverter chooses the algorithm
!   bc(mu) = 1,-1,0 for periodic,antiperiodic,fixed boundary conditions
!            respectively.
!   resmax is the stopping criterion for the iteration.
!   itermin is the minimum number of iterations required by the user.
!   mGMRES is the maximum subspace dimension
!   kGMRES is the number of eigenvectors to be deflated
!   hoption=1 to solve gamma5*M and 2 to solve M^dagger*M
! OUTPUT:
!   xe(), xo() is the solution vector.

    character(len=*), intent(in),    dimension(:)           :: rwdir
    real(kind=KR),    intent(in),    dimension(:,:,:,:,:)   :: be
    real(kind=KR),    intent(inout), dimension(:,:,:,:,:)   :: bo
    integer(kind=KI), intent(in)                            :: hoption,mGMRES,&
                                                               inverter, kGMRES
    real(kind=KR),    intent(in)                            :: kappa
    real(kind=KR),    intent(in)                            :: resmax
    real(kind=KR),    intent(in),    dimension(:,:,:)       :: coact
    integer(kind=KI), intent(in),    dimension(:)           :: bc, nms
    integer(kind=KI), intent(in)                            :: itermin, myid, &
                                                               MRT, MRT2
    !real(kind=KR),    intent(inout), dimension(:,:,:,:,:)   :: u
    real(kind=KR),    intent(in), dimension(:,:,:,:,:)   :: u
    real(kind=KR),    intent(inout), dimension(:,:,:,:,:)   :: xe, xo
    integer(kind=KI), intent(in),    dimension(:,:)         :: nn, iblv, &
                                                               vecbl, vecblinv
    logical,          intent(in),    dimension(:)           :: ldiv
    integer(kind=KI), intent(in),    dimension(:,:,:)       :: lvbc, lvcl
    integer(kind=KI), intent(in),    dimension(:,:,:,:)     :: ib
    logical,          intent(in),    dimension(:,:)         :: lbd





    !local variables

    integer(kind=KI)                           :: isite, itercount, gblclr, &
                                                  iflag, ntrue, ierr, Gext, &
                                                  idag, LUcount,i,j
    integer(kind=KI), dimension(2)             :: nGMRES
    integer(kind=KI)                           :: isignal, is
    real(kind=KR)                              :: fac1, fac2
    real(kind=KR), dimension(2)                :: kappatm
    real(kind=KR), dimension(6,ntotal,4,2,8)   :: vecsrc
    real(kind=KR), dimension(18,nvhalf,8,2,16) :: GeeGooinv
    integer(kind=KI) :: k,l,m,n,proc
    
    print *,"hoption=",hoption   

    iflag=-1 !Wilson

    kappatm(1)=kappa
    kappatm(2)=0.0_KR
   
    GeeGooinv=0.0_KR
    vecsrc = 0.0_KR
    nGMRES(1) = mGMRES
    nGMRES(2) = kGMRES
     
     !-Wilson: vecsrc = b_e/kappa^2+ H_eo*b_o/kappa.
      gblclr = 1
      idag = 0
      call Hsingle(xo,u,bo,idag,coact,bc,gblclr,vecbl, &
                   vecblinv,myid,nn,ldiv,nms,lvbc,ib,lbd,iblv,MRT)

      fac1 = 1.0_KR/kappatm(1)**2
      fac2 = 1.0_KR/kappatm(1)
      do isite = 1,nvhalf
       vecsrc(:,isite,:,:,:) = fac1*be(:,isite,:,:,:) + fac2*xo(:,isite,:,:,:)
      enddo ! isite

      !Start with a zero initial guess
      xe = 0.0_KR
      xo=  0.0_KR
    
      itercount = 0
      !-perform the matrix inversion using chosen algorithm.
     if (inverter==0) then
      call bicgstab(rwdir,vecsrc,xe,resmax,itermin,itercount,u, &
                   GeeGooinv,iflag,kappatm,coact,bc,vecbl,vecblinv,myid,nn, &
                   ldiv,nms,lvbc,ib,lbd,iblv,MRT,MRT2)
     elseif (inverter==1) then
      call cgne(rwdir,vecsrc,xe,resmax,itermin,itercount,u, &
                GeeGooinv,iflag,kappatm,coact,bc,vecbl,vecblinv,myid,nn,ldiv, &
                nms,lvbc,ib,lbd,iblv,MRT,MRT2)
     elseif (inverter==2) then
      call qmr(rwdir,vecsrc,xe,resmax,itermin,itercount,u, &
               GeeGooinv,iflag,kappatm,coact,bc,vecbl,vecblinv,myid,nn,ldiv, &
               nms,lvbc,ib,lbd,iblv,MRT,MRT2)
     elseif (inverter==3) then
      call qmrgam5(rwdir,vecsrc,xe,resmax,itermin,itercount,u, &
                   GeeGooinv,iflag,kappatm,coact,bc,vecbl,vecblinv,myid,nn, &
                   ldiv,nms,lvbc,ib,lbd,iblv,MRT,MRT2)
     elseif (inverter==4) then
      call gmres(rwdir,vecsrc,xe,nGMRES(1),resmax,itermin, &
                 itercount,LUcount,u,GeeGooinv,iflag,kappatm,coact,bc,vecbl, &
                 vecblinv,myid,nn,ldiv,nms,lvbc,ib,lbd,iblv,MRT,MRT2)
     elseif (inverter==5) then
      call gmresdr(rwdir,vecsrc,xe,nGMRES,resmax,itermin, &
                   itercount,u,GeeGooinv,iflag,kappatm,coact,bc,vecbl, &
                   vecblinv,myid,nn,ldiv,nms,lvbc,ib,lbd,iblv,MRT,MRT2)
     elseif (inverter==6) then

      !call psc(rwdir,vecsrc,xe,nGMRES,resmax,itermin, &
      !             itercount,u,GeeGooinv,iflag,kappatm,coact,bc,vecbl, &
      !             vecblinv,myid,nn,ldiv,nms,lvbc,ib,lbd,iblv,MRT,MRT2)

      call  landr_proj(rwdir,vecsrc,xe,nGMRES,resmax,itermin, &
                       itercount,u,GeeGooinv,iflag,kappatm,coact,bc,vecbl, &
                       vecblinv,myid,nn,ldiv,nms,lvbc,ib,lbd,iblv,hoption,MRT,MRT2)
!     elseif (inverter==8) then
!      call PP-gmresdr(rwdir,vecsrc,xe,nGMRES,resmax,itermin, &
!                   itercount,u,GeeGooinv,iflag,kappatm,coact,bc,vecbl, &
!                   vecblinv,myid,nn,ldiv,nms,lvbc,ib,lbd,iblv,MRT,MRT2)

      !call  cg(rwdir,vecsrc,xe,resmax,itermin, &
      !                 itercount,u,GeeGooinv,iflag,kappatm,coact,bc,vecbl, &
      !                 vecblinv,myid,nn,ldiv,nms,lvbc,ib,lbd,iblv,hoption,MRT,MRT2)

     endif ! inverter


! Record the number of iterations used to invert this matrix.
    if (myid==0) then
     open(unit=8,file=trim(rwdir(myid+1))//"CFGSPROPS.LOG",action="write",&
          form="formatted",status="old",position="append")
      if (inverter==0) then
       write(unit=8,fmt="(a32,i10)") &
             "bicgstab iterations: itercount =", itercount
      elseif (inverter==1) then
       write(unit=8,fmt="(a28,i10)") &
             "cgne iterations: itercount =", itercount
      elseif (inverter==2) then
       write(unit=8,fmt="(a27,i10)") &
             "qmr iterations: itercount =", itercount
      elseif (inverter==3) then
       write(unit=8,fmt="(a31,i10)") &
             "qmrgam5 iterations: itercount =", itercount
      elseif (inverter==4) then
       write(unit=8,fmt="(a32,i10)") &
             "gmres(n) iterations: itercount =", itercount
       write(unit=8,fmt="(a33,i10)") &
             "gmres(n) LU inversions: LUcount =", LUcount
      elseif (inverter==5) then
       write(unit=8,fmt="(a36,i10)") &
             "gmresdr(n,k) iterations: itercount =", itercount
      elseif (inverter==6) then
       write(unit=8,fmt="(a36,i10)") &
             "Dean - psc(n,k) iterations: itercount =", itercount
!      elseif (inverter==8) then
!       write(unit=8,fmt="(a36,i10)") &
!             "PP-gmresdr(n,k) iterations: itercount =", itercount
      endif
     close(unit=8,status="keep")
    endif


!         x_o = b_o + kappa*H_oe*x_e 
    
    gblclr = 2
    idag = 0
    call Hsingle(xo,u,xe,idag,coact,bc, &
                 gblclr,vecbl,vecblinv,myid,nn,ldiv,nms,lvbc,ib,lbd,iblv,MRT)
     do isite = 1,nvhalf
      xo(:,isite,:,:,:) = bo(:,isite,:,:,:) &
                               + kappa*xo(:,isite,:,:,:)
     enddo ! isite
 
 end subroutine fermprop_wilson
! - - - - - - - - - - - - - - - - - - - - - - - -- - - - - - - - - - - - - -
 subroutine fermprop(rwdir,be,bo,nkappa,kappa,cSW,inverter,coact,bc,resmax, &
                     itermin,omegaM3R,u,vectype,kGMRES,xe,xo,myid,nn,ldiv, &
                     nms,lvbc,lvcl,ib,lbd,iblv,vecbl,vecblinv,ntmqcd, &
                     hoption,MRT,MRT2)
! Compute fermion propagators.
! INPUT:
!   be() contains the entries of the input Dirac vector that are on even
!        lattice sites, with respect to the entire lattice (not blocked).
!   bo() contains the entries of the input Dirac vector that are on odd
!        lattice sites, with respect to the entire lattice (not blocked).
!   For BiCGstab or CGNE or QMR or GMRES...
!       nkappa is -1 (for Wilson) or -2 (for clover).
!       kappa(1) is 1/(8+2*m_0) with m_0 from eq.(1.1) of JHEP08(2001)058.
!       kappa(2) is mu_q from eq.(1.1) of JHEP08(2001)058.
!       cSW is the clover coefficient in the fermion action.
!       (cSW is irrelevant for nkappa=-1).
!       inverter=0 for bicgstab; inverter=1 for CGNE, inverter=2 for QMR;
!       inverter=3 for QMR(gamma5); inverter=4 for GMRES(n);
!       inverter=5 for GMRES-DR(n,k); inverter=6 for GMRES-Drshift(n,k).
!       inverter=7 for GMRES-DREIG (n,k) (gmresdr equiv to gmresdr from matlab; 
!                                         generate eigmodes stored in xshift.f90)
!       kGMRES is the "k" in GMRES-DR(n,k) and is irrelevant unless inverter=5
!       or inverter=6.
!       Recall that QMR(gamma5) is only valid for mu_q=0.
!   For minimal residual...
!       nkappa is the number of hopping parameters.
!       kappa() is the set of hopping parameters.
!       cSW is irrelevant.  The Wilson action is used.
!       inverter is irrelevant.
!   bc(mu) = 1,-1,0 for periodic,antiperiodic,fixed boundary conditions
!            respectively.
!   resmax is the stopping criterion for the iteration.
!   itermin is the minimum number of iterations required by the user.
!   omegaM3R is the overrelaxation parameter for the minimal residual inverter.
!   For CGNE or QMR...
!       vectype is irrelevant
!   For BiCGstab...
!       if vectype=1, then xe() is set to some initial vector.
!                     otherwise, the input vector xe() is used.
!   For minimal residual...
!       if vectype=1, then be() is assumed to be zero.
!       if vectype=2, then bo() is assumed to be zero.
!       if vectype=3, then neither be() nor bo() is assumed to be zero.
!   For GMRES(n) and GMRES-DR(n,k)...
!       vectype is "n".
! OUTPUT:
!   xe(), xo() are the solution vectors.

    character(len=*), intent(in),    dimension(:)           :: rwdir
    real(kind=KR),    intent(in),    dimension(:,:,:,:,:)   :: be
    real(kind=KR),    intent(inout), dimension(:,:,:,:,:)   :: bo
    integer(kind=KI), intent(in)                            :: nkappa, hoption,&
                                                               inverter, kGMRES
    real(kind=KR),    intent(in), dimension(:)              :: kappa
    real(kind=KR),    intent(in)                            :: cSW, resmax
    real(kind=KR),    intent(in),    dimension(:,:,:)       :: coact
    integer(kind=KI), intent(in)                            :: ntmqcd
    integer(kind=KI), intent(in),    dimension(:)           :: bc, nms
    integer(kind=KI), intent(in)                    :: itermin, myid, MRT, MRT2
    real(kind=KR2),   intent(in)                            :: omegaM3R
    real(kind=KR),    intent(inout), dimension(:,:,:,:,:)   :: u
    integer(kind=KI), intent(inout)                         :: vectype
    real(kind=KR),    intent(inout), dimension(:,:,:,:,:,:) :: xe, xo
    integer(kind=KI), intent(in),    dimension(:,:) :: nn, iblv, vecbl, vecblinv
    logical,          intent(in),    dimension(:)           :: ldiv
    integer(kind=KI), intent(in),    dimension(:,:,:)       :: lvbc, lvcl

! Note: lvcl() is the nearest neighbour site in +mu direction to be used by
! subroutine makeclover.  It should probably be identified with lv() and not
! lvbc(), but I've given it a new name -- lvcl() -- so as to remind the user
! that the routine calling fermprop has made a choice.
    integer(kind=KI), intent(in),    dimension(:,:,:,:)     :: ib
    logical,          intent(in),    dimension(:,:)         :: lbd

    integer(kind=KI)                           :: isite, itercount, gblclr, &
                                                  ikappa, ntrue, ierr, Gext, &
                                                  idag, LUcount,i,j
    integer(kind=KI), dimension(2)             :: nGMRES
    integer(kind=KI)                           :: isignal, is
    real(kind=KR)                              :: fac1, fac2,time1,time2
    real(kind=KR), dimension(6,ntotal,4,2,8)   :: vecsrc
    real(kind=KR), dimension(18,nvhalf,8,2,16) :: GeeGooinv

    real(kind=KR), dimension(6,ntotal,4,2,8,2)     :: bvec
    real(kind=KR), dimension(6,ntotal,4,2,8,2,nk1) :: xwhole
    real(kind=KR), dimension(6,ntotal,2,2,8,1)     :: xtemp
    integer(kind=KI) :: k,l,m,n,proc


   ! print *, "nkappa inside fermprop =", nkappa
    !Timing fermprop
    call cpu_time(time1) !start

    vecsrc = 0.0_KR
!     print *,"nkappa ori=",nkappa
!*OPTION 1: BICGSTAB/CGNE/QMR/GMRES(n) FOR TWISTED MASS (WILSON OR CLOVER).
    if (nkappa<0) then
!-identify the value of "n" for GMRES(n).  It is harmless for other inverters.
!     print *,"nkappa ori=",nkappa           
     nGMRES(1) = vectype
     nGMRES(2) = kGMRES
     if (nkappa==-2) then
     if (myid == 0) then
       print *, "Shouldn't be in here right?"
     endif ! myid
!-clover: vecsrc = b_e/kappa^2 + H_eo*G_oo^(-1)*b_o/kappa.
!         where G_oo = 1-kappa*C_oo+2*kappa*mu*i*gam_5, where i=sqrt(-1).
!               C_oo is the clover term on the globally-odd lattice sites.
      fac1 = 1.0_KR/kappa(1)
      do isite = 1,nvhalf
       vecsrc(:,isite,:,:,:) = fac1*bo(:,isite,:,:,:)
      enddo ! isite
      call makeclover(GeeGooinv,coact,u,myid,nn,ldiv,nms,lvcl,ib,lbd,iblv,MRT)
      call makeG(GeeGooinv,kappa,cSW,vecbl,Gext)
      if (nps==1) then
       fac1 = real(Gext,KR)
      else
       fac2 = real(Gext,KR)
       call MPI_REDUCE(fac2,fac1,1,MRT,MPI_SUM,0,MPI_COMM_WORLD,ierr)
      endif
      if (myid==0 .and. fac1>0.0_KR) then
       open(unit=8,file=trim(rwdir(myid+1))//"CFGSPROPS.LOG",action="write", &
            form="formatted",status="old",position="append")
        write(unit=8,fmt="(a24,es17.10)") "extra makeG iterations: ", fac1
       close(unit=8,status="keep")
      endif
      gblclr = 2
      call multG(GeeGooinv,vecsrc,gblclr,vecbl)
      gblclr = 1
      idag = 0
      call Hsingle(xo(:,:,:,:,:,1),u,vecsrc,idag,coact,bc,gblclr,vecbl, &
                   vecblinv,myid,nn,ldiv,nms,lvbc,ib,lbd,iblv,MRT)
      fac1 = 1.0_KR/kappa(1)**2
      do isite = 1,nvhalf
       vecsrc(:,isite,:,:,:) = fac1*be(:,isite,:,:,:) + xo(:,isite,:,:,:,1)
      enddo ! isite

     elseif (nkappa==-1) then
     
!-Wilson: vecsrc = b_e/kappa^2
!                  + H_eo*(1-2*kappa*mu*i*gam_5)*b_o/kappa/(1+(2*kappa*mu)^2).

      fac1 = 2.0_KR*kappa(1)*kappa(2)
      do isite = 1,nvhalf
       vecsrc(:,isite,1,:,:) = bo(:,isite,1,:,:) - fac1*bo(:,isite,3,:,:)
       vecsrc(:,isite,2,:,:) = bo(:,isite,2,:,:) - fac1*bo(:,isite,4,:,:)
       vecsrc(:,isite,3,:,:) = bo(:,isite,3,:,:) + fac1*bo(:,isite,1,:,:)
       vecsrc(:,isite,4,:,:) = bo(:,isite,4,:,:) + fac1*bo(:,isite,2,:,:)
      enddo ! isite
!
      gblclr = 1
      idag = 0
      call Hsingle(xo(:,:,:,:,:,1),u,vecsrc,idag,coact,bc,gblclr,vecbl, &
                   vecblinv,myid,nn,ldiv,nms,lvbc,ib,lbd,iblv,MRT)

      fac1 = 1.0_KR/kappa(1)**2
      fac2 = 1.0_KR/( kappa(1)*(1.0_KR+(2.0_KR*kappa(1)*kappa(2))**2) )
      do isite = 1,nvhalf
       vecsrc(:,isite,:,:,:) = fac1*be(:,isite,:,:,:) + fac2*xo(:,isite,:,:,:,1)
      enddo ! isite
     endif


!-for cgne, qmr, gmres(n), gmresdr(n,k), gmresdrEIG(n,k) or 1st bicgstab call, make initial
! solution vector.
!-for gmresdr-twisted multi-mass (inverter=6) make both even and odd inital solution
! vector. 

! COMMENT ~ the check for isignal is for inverter==6. For GMRESDR we want to
!           get the initial solution vector and for following RHS use the
!           the previous solution vector as the starting vector.

     if (inverter>3) then
        xe = 0.0_KR
     elseif (inverter/=0 .or. vectype==1) then
        xe(1,:,:,:,:,1) = 1.0_KR
        xe(2,:,:,:,:,1) = 0.0_KR
        xe(3,:,:,:,:,1) = 1.0_KR
        xe(4,:,:,:,:,1) = 0.0_KR
        xe(5,:,:,:,:,1) = 1.0_KR
        xe(6,:,:,:,:,1) = 0.0_KR
        vectype=0
     endif

    !itercount = 0
!-perform the matrix inversion using chosen algorithm.
!      print *,"nkappa sec=",nkappa
     if (inverter==0) then
      call bicgstab(rwdir,vecsrc,xe(:,:,:,:,:,1),resmax,itermin,itercount,u, &
                   GeeGooinv,nkappa,kappa,coact,bc,vecbl,vecblinv,myid,nn, &
                   ldiv,nms,lvbc,ib,lbd,iblv,MRT,MRT2)
     elseif (inverter==1) then
      call cgne(rwdir,vecsrc,xe(:,:,:,:,:,1),resmax,itermin,itercount,u, &
                GeeGooinv,nkappa,kappa,coact,bc,vecbl,vecblinv,myid,nn,ldiv, &
                nms,lvbc,ib,lbd,iblv,MRT,MRT2)
     elseif (inverter==2) then
      call qmr(rwdir,vecsrc,xe(:,:,:,:,:,1),resmax,itermin,itercount,u, &
               GeeGooinv,nkappa,kappa,coact,bc,vecbl,vecblinv,myid,nn,ldiv, &
               nms,lvbc,ib,lbd,iblv,MRT,MRT2)
     elseif (inverter==3) then
      call qmrgam5(rwdir,vecsrc,xe(:,:,:,:,:,1),resmax,itermin,itercount,u, &
                   GeeGooinv,nkappa,kappa,coact,bc,vecbl,vecblinv,myid,nn, &
                   ldiv,nms,lvbc,ib,lbd,iblv,MRT,MRT2)
     elseif (inverter==4) then
      call gmres(rwdir,vecsrc,xe(:,:,:,:,:,1),nGMRES(1),resmax,itermin, &
                 itercount,LUcount,u,GeeGooinv,nkappa,kappa,coact,bc,vecbl, &
                 vecblinv,myid,nn,ldiv,nms,lvbc,ib,lbd,iblv,MRT,MRT2)
     elseif (inverter==5) then
      call gmresdr(rwdir,vecsrc,xe(:,:,:,:,:,1),nGMRES,resmax,itermin, &
                   itercount,u,GeeGooinv,nkappa,kappa,coact,bc,vecbl, &
                   vecblinv,myid,nn,ldiv,nms,lvbc,ib,lbd,iblv,MRT,MRT2)
     elseif (inverter==6) then
!
! If one wants to do a deflated GMRES followed by GMRES-Proj, then 
! uncomment (or make active) the psc call below.
! 
      call psc(rwdir,vecsrc,xe(:,:,:,:,:,1),nGMRES,resmax,itermin, &
                   itercount,u,GeeGooinv,nkappa,kappa,coact,bc,vecbl, &
                   vecblinv,myid,nn,ldiv,nms,lvbc,ib,lbd,iblv,MRT,MRT2)

   !   call  landr_proj(rwdir,vecsrc,xe(:,:,:,:,:,1),nGMRES,resmax,itermin, &
   !                    itercount,u,GeeGooinv,nkappa,kappa,coact,bc,vecbl, &
   !                    vecblinv,myid,nn,ldiv,nms,lvbc,ib,lbd,iblv,hoption,MRT,MRT2)

      !call  cg(rwdir,vecsrc,xe(:,:,:,:,:,1),resmax,itermin, &
      !         itercount,u,GeeGooinv,nkappa,kappa,coact,bc,vecbl, &
      !         vecblinv,myid,nn,ldiv,nms,lvbc,ib,lbd,iblv,hoption,MRT,MRT2)
     elseif (inverter==7) then
      call gmresdrEIG(rwdir,vecsrc,xe(:,:,:,:,:,1),nGMRES,resmax,itermin, &
                   itercount,u,GeeGooinv,nkappa,idag,kappa,coact,bc,vecbl, &
                   vecblinv,myid,nn,ldiv,nms,lvbc,ib,lbd,iblv,MRT,MRT2)
     elseif (inverter==8) then
      call ppgmresdr(rwdir,vecsrc,xe(:,:,:,:,:,1),nGMRES,resmax,itermin, &
                   itercount,u,GeeGooinv,nkappa,kappa,coact,bc,vecbl, &
                   vecblinv,myid,nn,ldiv,nms,lvbc,ib,lbd,iblv,MRT,MRT2)
!      print *,"nkappa thi=",nkappa
     endif ! inverter
!
! This "else" option below is for the case nkappa .ne. 0, which implies we
! want to multi-mass, and the only option for this is M3R.
!
!*OPTION 2: M3R FOR MULTI-MASS (WILSON).
    else ! M3R
!-perform the matrix inversion using a multi-mass minimal residual algorithm,
! ...for  b_e=0, b_o/=0
     if (vectype==1) then
      xe = 0.0_KR
      gblclr = 1
      idag = 0
      call Hsingle(vecsrc,u,bo,idag,coact,bc,gblclr,vecbl,vecblinv,myid,nn, &
                   ldiv,nms,lvbc,ib,lbd,iblv,MRT)
      call m3r(rwdir,vecsrc,xo,resmax,itermin,itercount,omegaM3R,u,nkappa, &
               kappa,coact,bc,vecbl,vecblinv,myid,nn,ldiv,nms,lvbc,ib,lbd, &
               iblv,MRT,MRT2)
! ...for  b_e/=0, b_o=0
     elseif (vectype==2) then
      call m3r(rwdir,be,xe,resmax,itermin,itercount,omegaM3R,u,nkappa,kappa, &
               coact,bc,vecbl,vecblinv,myid,nn,ldiv,nms,lvbc,ib,lbd,iblv,MRT, &
               MRT2)
      xo = 0.0_KR
! ...for  b_e/=0, b_o/=0
     else
      call m3r(rwdir,be,xe,resmax,itermin,itercount,omegaM3R,u,nkappa,kappa, &
               coact,bc,vecbl,vecblinv,myid,nn,ldiv,nms,lvbc,ib,lbd,iblv,MRT, &
               MRT2)
      gblclr = 1
      idag = 0
      call Hsingle(vecsrc,u,bo,idag,coact,bc,gblclr,vecbl,vecblinv,myid,nn, &
                   ldiv,nms,lvbc,ib,lbd,iblv,MRT)
      call m3r(rwdir,vecsrc,xo,resmax,itermin,itercount,omegaM3R,u,nkappa, &
               kappa,coact,bc,vecbl,vecblinv,myid,nn,ldiv,nms,lvbc,ib,lbd, &
               iblv,MRT,MRT2)
     endif
!-build the full result for xe().
     do ikappa = 1,nkappa
      print *,"nkappa=",nkappa                              
      do isite = 1,nvhalf
       xe(:,isite,:,:,:,ikappa) = xe(:,isite,:,:,:,ikappa)/kappa(ikappa)**2 &
                                + xo(:,isite,:,:,:,ikappa)/kappa(ikappa)
      enddo ! isite
     enddo ! ikappa
    endif ! M3R

! Record the number of iterations used to invert this matrix.
    if (myid==0) then
     open(unit=8,file=trim(rwdir(myid+1))//"CFGSPROPS.LOG",action="write",&
          form="formatted",status="old",position="append")
      if (nkappa<0 .and. inverter==0) then
       write(unit=8,fmt="(a32,i10)") &
             "bicgstab iterations: itercount =", itercount
      elseif (nkappa<0 .and. inverter==1) then
       write(unit=8,fmt="(a28,i10)") &
             "cgne iterations: itercount =", itercount
      elseif (nkappa<0 .and. inverter==2) then
       write(unit=8,fmt="(a27,i10)") &
             "qmr iterations: itercount =", itercount
      elseif (nkappa<0 .and. inverter==3) then
       write(unit=8,fmt="(a31,i10)") &
             "qmrgam5 iterations: itercount =", itercount
      elseif (nkappa<0 .and. inverter==4) then
       write(unit=8,fmt="(a32,i10)") &
             "gmres(n) iterations: itercount =", itercount
       write(unit=8,fmt="(a33,i10)") &
             "gmres(n) LU inversions: LUcount =", LUcount
      elseif (nkappa<0 .and. inverter==5) then
       write(unit=8,fmt="(a36,i10)") &
             "gmresdr(n,k) iterations: itercount =", itercount
      elseif (nkappa<0 .and. inverter==6) then
       write(unit=8,fmt="(a36,i10)") &
             "Dean - psc(n,k) iterations: itercount =", itercount
      elseif (nkappa<0 .and. inverter==8) then
       write(unit=8,fmt="(a36,i10)") &
             "Andy - ppgmres iterations: itercount =", itercount      
      elseif (nkappa>0) then
       write(unit=8,fmt="(a40,i10)") &
             "minimal residual iterations: itercount =", itercount
      endif
     close(unit=8,status="keep")
    endif

! For multi-mass Wilson, construct
!         x_o = b_o + kappa*H_oe*x_e.
! For twisted mass Wilson, construct
!         x_o = ( 1 - 2*kappa*mu*i*gam_5 ) / ( 1 + (2*kappa*mu)^2  )
!               * ( b_o + kappa*H_oe*x_e ).
! For twisted mass clover, construct
!         x_o = ( 1 - kappa*C_oo + 2*kappa*mu*i*gam_5 )^(-1)
!               * ( b_o + kappa*H_oe*x_e ).

    if (nkappa<0) then
     ntrue = 1
    else
     ntrue = nkappa
    endif
    gblclr = 2
    idag = 0


! All cases
    do ikappa = 1,ntrue
     call Hsingle(xo(:,:,:,:,:,ikappa),u,xe(:,:,:,:,:,ikappa),idag,coact,bc, &
                  gblclr,vecbl,vecblinv,myid,nn,ldiv,nms,lvbc,ib,lbd,iblv,MRT)
     do isite = 1,nvhalf
      xo(:,isite,:,:,:,ikappa) = bo(:,isite,:,:,:) &
                               + kappa(ikappa)*xo(:,isite,:,:,:,ikappa)
     enddo ! isite
    enddo ! ikappa


    if (nkappa==-1) then
! tmqcd case
     fac1 = 1.0_KR/( 1.0_KR + (2.0_KR*kappa(1)*kappa(2))**2 )
     do isite = 1,nvhalf
      vecsrc(:,isite,:,:,:) = fac1*xo(:,isite,:,:,:,1)
     enddo ! isite
     fac1 = 2.0_KR*kappa(1)*kappa(2)
     do isite = 1,nvhalf
      xo(:,isite,1,:,:,1) = vecsrc(:,isite,1,:,:) - fac1*vecsrc(:,isite,3,:,:)
      xo(:,isite,2,:,:,1) = vecsrc(:,isite,2,:,:) - fac1*vecsrc(:,isite,4,:,:)
      xo(:,isite,3,:,:,1) = vecsrc(:,isite,3,:,:) + fac1*vecsrc(:,isite,1,:,:)
      xo(:,isite,4,:,:,1) = vecsrc(:,isite,4,:,:) + fac1*vecsrc(:,isite,2,:,:)
     enddo ! isite
    elseif (nkappa==-2) then
! Clover case
     gblclr = 2
     call multG(GeeGooinv,xo(:,:,:,:,:,1),gblclr,vecbl)
    endif


  if (ntmqcd /= 0.and.abs(kappa(2)) > 0.0001) then
    call gamma5vector(xo(:,:,:,:,:,1),ntmqcd,myid)
    call gamma5vector(xe(:,:,:,:,:,1),ntmqcd,myid)
  endif ! kappa

  call cpu_time(time2)


  if (myid==0) then
    open(unit=8,file=trim(rwdir(myid+1))//"CFGSPROPS.LOG",action="write", &
            form="formatted",status="old",position="append")
    write(unit=8,fmt="(a50,es17.10)") "Time taken by fermprop in seconds: ", time2-time1
    close(unit=8,status="keep")
  endif

 end subroutine fermprop

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 subroutine fermpropR(rwdir,be,bo,nkappa,kappa,cSW,inverter,coact,bc,resmax, &
                     itermin,omegaM3R,u,vectype,kGMRES,xe,xo,myid,nn,ldiv, &
                     nms,lvbc,lvcl,ib,lbd,iblv,vecbl,vecblinv,ntmqcd,&
                     hoption,MRT,MRT2,ir)
! Compute fermion propagators.
! INPUT:
!   be() contains the entries of the input Dirac vector that are on even
!        lattice sites, with respect to the entire lattice (not blocked).
!   bo() contains the entries of the input Dirac vector that are on odd
!        lattice sites, with respect to the entire lattice (not blocked).
!   For BiCGstab or CGNE or QMR or GMRES...
!       nkappa is -1 (for Wilson) or -2 (for clover).
!       kappa(1) is 1/(8+2*m_0) with m_0 from eq.(1.1) of JHEP08(2001)058.
!       kappa(2) is mu_q from eq.(1.1) of JHEP08(2001)058.
!       cSW is the clover coefficient in the fermion action.
!       (cSW is irrelevant for nkappa=-1).
!       inverter=0 for bicgstab; inverter=1 for CGNE, inverter=2 for QMR;
!       inverter=3 for QMR(gamma5); inverter=4 for GMRES(n);
!       inverter=5 for GMRES-DR(n,k); inverter=6 for GMRES-Drshift(n,k).
!       kGMRES is the "k" in GMRES-DR(n,k) and is irrelevant unless inverter=5
!       or inverter=6.
!       Recall that QMR(gamma5) is only valid for mu_q=0.
!   For minimal residual...
!       nkappa is the number of hopping parameters.
!       kappa() is the set of hopping parameters.
!       cSW is irrelevant.  The Wilson action is used.
!       inverter is irrelevant.
!   bc(mu) = 1,-1,0 for periodic,antiperiodic,fixed boundary conditions
!            respectively.
!   resmax is the stopping criterion for the iteration.
!   itermin is the minimum number of iterations required by the user.
!   omegaM3R is the overrelaxation parameter for the minimal residual inverter.
!   For CGNE or QMR...
!       vectype is irrelevant
!   For BiCGstab...
!       if vectype=1, then xe() is set to some initial vector.
!                     otherwise, the input vector xe() is used.
!   For minimal residual...
!       if vectype=1, then be() is assumed to be zero.
!       if vectype=2, then bo() is assumed to be zero.
!       if vectype=3, then neither be() nor bo() is assumed to be zero.
!   For GMRES(n) and GMRES-DR(n,k)...
!       vectype is "n".
!   ir value is the Wilson r coefficient.
! OUTPUT:
!   xe(), xo() are the solution vectors.

    character(len=*), intent(in),    dimension(:)           :: rwdir
    real(kind=KR),    intent(in),    dimension(:,:,:,:,:)   :: be
    real(kind=KR),    intent(inout), dimension(:,:,:,:,:)   :: bo
    integer(kind=KI), intent(in)                            :: nkappa, hoption,&
                                                               inverter, kGMRES, ir
    real(kind=KR),    intent(in), dimension(:)              :: kappa
    real(kind=KR),    intent(in)                            :: cSW, resmax
    real(kind=KR),    intent(in),    dimension(:,:,:)       :: coact
    integer(kind=KI), intent(in)                            :: ntmqcd
    integer(kind=KI), intent(in),    dimension(:)           :: bc, nms
    integer(kind=KI), intent(in)                    :: itermin, myid, MRT, MRT2
    real(kind=KR2),   intent(in)                            :: omegaM3R
    real(kind=KR),    intent(inout), dimension(:,:,:,:,:)   :: u
    integer(kind=KI), intent(inout)                         :: vectype
    real(kind=KR),    intent(inout), dimension(:,:,:,:,:,:) :: xe, xo
    integer(kind=KI), intent(in),    dimension(:,:) :: nn, iblv, vecbl, vecblinv
    logical,          intent(in),    dimension(:)           :: ldiv
    integer(kind=KI), intent(in),    dimension(:,:,:)       :: lvbc, lvcl

! Note: lvcl() is the nearest neighbour site in +mu direction to be used by
! subroutine makeclover.  It should probably be identified with lv() and not
! lvbc(), but I've given it a new name -- lvcl() -- so as to remind the user
! that the routine calling fermprop has made a choice.
    integer(kind=KI), intent(in),    dimension(:,:,:,:)     :: ib
    logical,          intent(in),    dimension(:,:)         :: lbd
    real(kind=KR),                   dimension(4,4,4)       :: coactr
    integer(kind=KI)::i

    coactr = coact
    if (ir.eq.-1) then
     do i=1,4
      coactr(3,1,i) = -coact(3,1,i)
     enddo
    endif
    call fermprop(rwdir,be,bo,nkappa,kappa,cSW,inverter,coactr,bc,resmax, &
                     itermin,omegaM3R,u,vectype,kGMRES,xe,xo,myid,nn,ldiv, &
                     nms,lvbc,lvcl,ib,lbd,iblv,vecbl,vecblinv,ntmqcd,hoption,MRT,MRT2)

    end subroutine fermpropR

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 subroutine makeG(GeeGooinv,kappa,cSW,vecbl,Gext)
! INPUT:
!   GeeGooinv contains the clover term C.
!   kappa(1) is the hopping parameter.
!   kappa(2) is the twisted mass parameter, mu_q.
!   cSW is the clover coefficient.
! OUTPUT:
!   Gext is the number of extra LU inversions required.
!   GeeGooinv is returned as G_ee on globally-even sites and G_oo^(-1) on
!             globally-odd sites, where G = 1-kappa*C+2*kappa*mu*i*gam_5
!             and i = sqrt(-1).
!   Note that G = /   V  i*W \ = (1/2) /  1 1 \  /V+W 0 \  / 1  i \
!                 \ -i*W  V  /         \ -i i /  \ 0 V-W/  \ 1 -i /
!
!        and G^(-1) = (1/2) /  1 1 \  /(V+W)^(-1)      0    \  / 1  i \
!                           \ -i i /  \     0     (V-W)^(-1)/  \ 1 -i /
!
!   so the 12x12 matrix G can be inverted using two 6x6 inversions: V+W and V-W
!   as discussed by Gockeler et al., NPB(PS)53, 312 (1997).
!   Therefore, GeeGooinv(:,:,a,:,:) is
!   V_(11) for a=1, V_(12) for a=2, V_(21) for a=3, V_(22) for a=4,
!   W_(11) for a=5, W_(12) for a=6, W_(21) for a=7, W_(22) for a=8.

    real(kind=KR),    intent(inout), dimension(:,:,:,:,:) :: GeeGooinv
    real(kind=KR),    intent(in),    dimension(:)         :: kappa
    real(kind=KR),    intent(in)                          :: cSW
    integer(kind=KI), intent(in),    dimension(:,:)       :: vecbl
    integer(kind=KI), intent(out)                         :: Gext

    integer(kind=KI)                :: icri, isite, ieo, ibl, iri, icd, jcd, &
                                       ncd, ibleo, ipauli, jpauli
    real(kind=KR)                   :: bit
    real(kind=KR), dimension(2,6,6) :: Gbit, V, W, VpWinv, VmWinv

! Initialize matrices that will only get defined via a subroutine call.
    VpWinv = 0.0_KR
    VmWinv = 0.0_KR

! Construct the matrix G.
    bit = -kappa(1)*cSW
    GeeGooinv = bit*GeeGooinv
    GeeGooinv( 1,:,1,:,:) = GeeGooinv( 1,:,1,:,:) + 1.0_KR
    GeeGooinv( 9,:,1,:,:) = GeeGooinv( 9,:,1,:,:) + 1.0_KR
    GeeGooinv(17,:,1,:,:) = GeeGooinv(17,:,1,:,:) + 1.0_KR
    GeeGooinv( 1,:,4,:,:) = GeeGooinv( 1,:,4,:,:) + 1.0_KR
    GeeGooinv( 9,:,4,:,:) = GeeGooinv( 9,:,4,:,:) + 1.0_KR
    GeeGooinv(17,:,4,:,:) = GeeGooinv(17,:,4,:,:) + 1.0_KR
    bit = 2.0_KR*kappa(1)*kappa(2)
    GeeGooinv( 2,:,5,:,:) = GeeGooinv( 2,:,5,:,:) - bit
    GeeGooinv(10,:,5,:,:) = GeeGooinv(10,:,5,:,:) - bit
    GeeGooinv(18,:,5,:,:) = GeeGooinv(18,:,5,:,:) - bit
    GeeGooinv( 2,:,8,:,:) = GeeGooinv( 2,:,8,:,:) - bit
    GeeGooinv(10,:,8,:,:) = GeeGooinv(10,:,8,:,:) - bit
    GeeGooinv(18,:,8,:,:) = GeeGooinv(18,:,8,:,:) - bit

! Invert G on globally-odd sites.  (Leave globally-even sites intact.)
    ncd = 6
    do ibleo = 1,8
     ibl = vecbl(2,ibleo)
     do ieo = 1,2
      do isite = 1,nvhalf
!-step 1: invert V+W and V-W.
       do icri = 1,18
        iri = 1 + modulo(icri-1,2)
        do ipauli = 1,4
         jpauli = 4 + ipauli
         icd = 1 + (icri-1)/6 + 3*((ipauli-1)/2)
         jcd = 1 + modulo((icri-1)/2,3) + 3*modulo(ipauli-1,2)
         V(iri,icd,jcd) = GeeGooinv(icri,isite,ipauli,ieo,ibl)
         W(iri,icd,jcd) = GeeGooinv(icri,isite,jpauli,ieo,ibl)
        enddo ! ipauli
       enddo ! icri
       Gbit = V + W
       Gext = 0
       call LUinversion(Gbit,ncd,VpWinv,Gext)
       Gbit = V - W
       call LUinversion(Gbit,ncd,VmWinv,Gext)
!-step 2: construct G^(-1).
       V = 0.5_KR*(VpWinv+VmWinv)
       W = 0.5_KR*(VpWinv-VmWinv)
       do icri = 1,18
        iri = 1 + modulo(icri-1,2)
        do ipauli = 1,4
         jpauli = 4 + ipauli
         icd = 1 + (icri-1)/6 + 3*((ipauli-1)/2)
         jcd = 1 + modulo((icri-1)/2,3) + 3*modulo(ipauli-1,2)
         GeeGooinv(icri,isite,ipauli,ieo,ibl) = V(iri,icd,jcd)
         GeeGooinv(icri,isite,jpauli,ieo,ibl) = W(iri,icd,jcd)
        enddo ! ipauli
       enddo ! icri
      enddo ! isite
     enddo ! ieo
    enddo ! ibleo

 end subroutine makeG

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 subroutine pionrho(pion,rho,decay,xe,xo,idsrc,ntmqcd,myid,MRT)
! Using local operators, compute the pion and rho correlators and pion decay.
! Only the third component is used for the rho operator.
! Only the fourth component is used for axial operator in pion decay.
! NOTE: For maximally twisted mass QCD (alpha=pi/2), this "pion decay"
!       amplitude will vanish.  The pion then decays via a "vector" operator.
! NOTE: it=1 is uniquely defined by          1<=isite<=  factor and 1<=ibl<= 8
!       it=2 is uniquely defined by          1<=isite<=  factor and 9<=ibl<=16
!       it=3 is uniquely defined by   factor+1<=isite<=2*factor and 1<=ibl<= 8
!       it=4 is uniquely defined by   factor+1<=isite<=2*factor and 9<=ibl<=16
!       it=5 is uniquely defined by 2*factor+1<=isite<=3*factor and 1<=ibl<= 8
!       it=6 is uniquely defined by 2*factor+1<=isite<=3*factor and 9<=ibl<=16
!       etc, where factor = 2*nvhalf*npt/nt.

    real(kind=KR),    intent(inout), dimension(:)         :: pion, rho, decay
    real(kind=KR),    intent(inout), dimension(:,:,:,:,:) :: xe, xo
!   real(kind=KR),    intent(in),    dimension(:,:,:,:,:) :: xe, xo
    integer(kind=KI), intent(in)                          :: idsrc, myid, MRT
    integer(kind=KI), intent(in)                          :: ntmqcd

    integer(kind=KI) :: isite, ieo, ibleo, nsite, sitemin, sitemax, it, id, &
                        icri, itbit, ierr
    real(kind=KR), dimension(6,ntotal,4,2,8) :: gam5xe, gam5xo
    integer(kind=KI), dimension(4)           :: np, ip
    real(kind=KR2),   dimension(nt)          :: pbit, rbit, dbit
    real(kind=KR),    dimension(nt)          :: pcorr, rcorr, dcorr, psum, rsum, dsum
    real(kind=KR),    dimension(4)           :: rsgn

! COMMENT ~ Due to the twisted calculation that I am doing I put the next
!           x*gamma_5 multiplications in so that we get the correct neutral 
!           pion which is calculated here for the twisted mass basis! BE CAREFUL,
!           Dr Lewis analysis does this afterward and these next 10 lines 
!           may cause a potential problem.

! NOTE    ~ Only do this for the pseudo-scalar pion. The vectors and axial 
! currents mix as well, but that is not taken care of here.

    gam5xe = 0.0_KR
    gam5xo = 0.0_KR
 
!   if (ntmqcd/=0)then
!     do isite = 1,nvhalf
!       gam5xe(:,isite,1,:,:) = -xe(:,isite,3,:,:)
!       gam5xe(:,isite,2,:,:) = -xe(:,isite,4,:,:)
!       gam5xe(:,isite,3,:,:) =  xe(:,isite,1,:,:)
!       gam5xe(:,isite,4,:,:) =  xe(:,isite,2,:,:)
!
!       gam5xo(:,isite,1,:,:) = -xo(:,isite,3,:,:)
!       gam5xo(:,isite,2,:,:) = -xo(:,isite,4,:,:)
!       gam5xo(:,isite,3,:,:) =  xo(:,isite,1,:,:)
!       gam5xo(:,isite,4,:,:) =  xo(:,isite,2,:,:)
!     enddo ! isite
!   else ! ntmqcd==0
        gam5xe =  xe
        gam5xo =  xo
!   endif ! ntmqcd

! IMPORTANT!!!!!
! Note that the neutral pion Tr[D_u D_u] is now Tr[gamma_5 D_u gamma_5 (D_d)^dagger] in
! the twisted basis. However, we have returned to the pysical (Wilson) basis in the
! beggining of the calculation by rotating the source and solution vectors by maximal
! twist rotations. 

! Also, note that in the coding since there is an overall factor of (i)^2 the
! final solution should be multiplied by a minus sign. So the overall sign is correct
! since the original trace also had a negative sign.

! NEED TO FIX THE TRACE OF THE NEUTEAL PION....NEED TO MAKE ONE OF THE
!  TERMS THE TRANSPOSE OF XE,XO !!!!!!

    pbit = 0.0_KR
    rbit = 0.0_KR
    dbit = 0.0_KR
    rsgn(1:4) = (/ 1.0_KR, -1.0_KR, -1.0_KR, 1.0_KR /) 
    np(1) = npx
    np(2) = npy
    np(3) = npz
    np(4) = npt
    call atoc(myid,np,ip)
    nsite = 2*nvhalf*npt/nt
    sitemin = 1
    sitemax = nsite
    do itbit = 2,nt/npt,2
     it = itbit + ip(4)*nt/npt
     do isite = sitemin,sitemax
      do ibleo = 1,4
       do ieo = 1,2
        do icri = 1,6
         if (it==4) then
            gam5xe(icri,isite,1,ieo,ibleo) = -xe(icri,isite,3,ieo,ibleo)
            gam5xe(icri,isite,2,ieo,ibleo) = -xe(icri,isite,4,ieo,ibleo)
            gam5xe(icri,isite,3,ieo,ibleo) =  xe(icri,isite,1,ieo,ibleo)
            gam5xe(icri,isite,4,ieo,ibleo) =  xe(icri,isite,2,ieo,ibleo)

            gam5xo(icri,isite,1,ieo,ibleo) = -xo(icri,isite,3,ieo,ibleo)
            gam5xo(icri,isite,2,ieo,ibleo) = -xo(icri,isite,4,ieo,ibleo)
            gam5xo(icri,isite,3,ieo,ibleo) =  xo(icri,isite,1,ieo,ibleo)
            gam5xo(icri,isite,4,ieo,ibleo) =  xo(icri,isite,2,ieo,ibleo)
         elseif (it==nt) then
            gam5xe(icri,isite,1,ieo,ibleo+4) = -xe(icri,isite,3,ieo,ibleo+4)
            gam5xe(icri,isite,2,ieo,ibleo+4) = -xe(icri,isite,4,ieo,ibleo+4)
            gam5xe(icri,isite,3,ieo,ibleo+4) =  xe(icri,isite,1,ieo,ibleo+4)
            gam5xe(icri,isite,4,ieo,ibleo+4) =  xe(icri,isite,2,ieo,ibleo+4)

            gam5xo(icri,isite,1,ieo,ibleo+4) = -xo(icri,isite,3,ieo,ibleo+4)
            gam5xo(icri,isite,2,ieo,ibleo+4) = -xo(icri,isite,4,ieo,ibleo+4)
            gam5xo(icri,isite,3,ieo,ibleo+4) =  xo(icri,isite,1,ieo,ibleo+4)
            gam5xo(icri,isite,4,ieo,ibleo+4) =  xo(icri,isite,2,ieo,ibleo+4)
         endif ! it
         do id = 1,4
          pbit(it-1) = pbit(it-1) + real(gam5xe(icri,isite,id,ieo,ibleo),KR2)**2 &
                                  + real(gam5xo(icri,isite,id,ieo,ibleo),KR2)**2
          pbit(it) = pbit(it) + real(gam5xe(icri,isite,id,ieo,ibleo+4),KR2)**2 &
                              + real(gam5xo(icri,isite,id,ieo,ibleo+4),KR2)**2
         enddo ! id
         rbit(it-1) = rbit(it-1) + real(xe(icri,isite,1,ieo,ibleo),KR2)**2 &
                                 + real(xo(icri,isite,1,ieo,ibleo),KR2)**2 &
                                 - real(xe(icri,isite,2,ieo,ibleo),KR2)**2 &
                                 - real(xo(icri,isite,2,ieo,ibleo),KR2)**2 &
                                 - real(xe(icri,isite,3,ieo,ibleo),KR2)**2 &
                                 - real(xo(icri,isite,3,ieo,ibleo),KR2)**2 &
                                 + real(xe(icri,isite,4,ieo,ibleo),KR2)**2 &
                                 + real(xo(icri,isite,4,ieo,ibleo),KR2)**2
         rbit(it) = rbit(it) + real(xe(icri,isite,1,ieo,ibleo+4),KR2)**2 &
                             + real(xo(icri,isite,1,ieo,ibleo+4),KR2)**2 &
                             - real(xe(icri,isite,2,ieo,ibleo+4),KR2)**2 &
                             - real(xo(icri,isite,2,ieo,ibleo+4),KR2)**2 &
                             - real(xe(icri,isite,3,ieo,ibleo+4),KR2)**2 &
                             - real(xo(icri,isite,3,ieo,ibleo+4),KR2)**2 &
                             + real(xe(icri,isite,4,ieo,ibleo+4),KR2)**2 &
                             + real(xo(icri,isite,4,ieo,ibleo+4),KR2)**2
         dbit(it-1) = dbit(it-1) + real(xe(icri,isite,1,ieo,ibleo),KR2)**2 &
                                 + real(xo(icri,isite,1,ieo,ibleo),KR2)**2 &
                                 + real(xe(icri,isite,2,ieo,ibleo),KR2)**2 &
                                 + real(xo(icri,isite,2,ieo,ibleo),KR2)**2 &
                                 - real(xe(icri,isite,3,ieo,ibleo),KR2)**2 &
                                 - real(xo(icri,isite,3,ieo,ibleo),KR2)**2 &
                                 - real(xe(icri,isite,4,ieo,ibleo),KR2)**2 &
                                 - real(xo(icri,isite,4,ieo,ibleo),KR2)**2
         dbit(it) = dbit(it) + real(xe(icri,isite,1,ieo,ibleo+4),KR2)**2 &
                             + real(xo(icri,isite,1,ieo,ibleo+4),KR2)**2 &
                             + real(xe(icri,isite,2,ieo,ibleo+4),KR2)**2 &
                             + real(xo(icri,isite,2,ieo,ibleo+4),KR2)**2 &
                             - real(xe(icri,isite,3,ieo,ibleo+4),KR2)**2 &
                             - real(xo(icri,isite,3,ieo,ibleo+4),KR2)**2 &
                             - real(xe(icri,isite,4,ieo,ibleo+4),KR2)**2 &
                             - real(xo(icri,isite,4,ieo,ibleo+4),KR2)**2
        enddo ! icri
       enddo ! ieo
      enddo ! ibleo
     enddo ! isite
     sitemin = sitemin + nsite
     sitemax = sitemax + nsite
    enddo ! itbit
    pcorr = real(pbit,KR)
    rcorr = rsgn(idsrc)*real(rbit,KR)
    dcorr = real(dbit,KR)
    if (nps==1) then
     pion = pion + pcorr
     rho = rho + rcorr
     decay = decay + dcorr
    else
     call MPI_REDUCE(pcorr(1),psum(1),nt,MRT,MPI_SUM,0,MPI_COMM_WORLD,ierr)
     call MPI_BCAST(psum(1),nt,MRT,0,MPI_COMM_WORLD,ierr)
     call MPI_REDUCE(rcorr(1),rsum(1),nt,MRT,MPI_SUM,0,MPI_COMM_WORLD,ierr)
     call MPI_BCAST(rsum(1),nt,MRT,0,MPI_COMM_WORLD,ierr)
     call MPI_REDUCE(dcorr(1),dsum(1),nt,MRT,MPI_SUM,0,MPI_COMM_WORLD,ierr)
     call MPI_BCAST(dsum(1),nt,MRT,0,MPI_COMM_WORLD,ierr)

! The pion, rho, and pion-decay are summed here for each final dirac source
! that is looped on in the subroutine generator. So the output of this program
! for example, up to some normalization, is the physical 
! neutral pion u*ubar + d*dbar.

     pion = pion + psum
     rho = rho + rsum
     decay = decay + dsum
    endif

 end subroutine pionrho

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 subroutine SSTpion(pion,xe,xo,src,icsrc,idsrc,myid,iblv,vecblinv,MRT)
! Using a local source, compute the pion correlator from existing sequential
! source propagators.
! NOTE: it=1 is uniquely defined by          1<=isite<=  factor and 1<=ibl<= 8
!       it=2 is uniquely defined by          1<=isite<=  factor and 9<=ibl<=16
!       it=3 is uniquely defined by   factor+1<=isite<=2*factor and 1<=ibl<= 8
!       it=4 is uniquely defined by   factor+1<=isite<=2*factor and 9<=ibl<=16
!       it=5 is uniquely defined by 2*factor+1<=isite<=3*factor and 1<=ibl<= 8
!       it=6 is uniquely defined by 2*factor+1<=isite<=3*factor and 9<=ibl<=16
!       etc, where factor = 2*nvhalf*npt/nt.

    real(kind=KR),    intent(inout)                      :: pion
    real(kind=KR),    intent(in),   dimension(:,:,:,:,:) :: xe, xo
    integer(kind=KI), intent(in),   dimension(:)         :: src
    integer(kind=KI), intent(in)                         :: icsrc, idsrc, &
                                                            myid, MRT
    integer(kind=KI), intent(in),   dimension(:,:)       :: iblv, vecblinv

    integer(kind=KI)               :: isite, ieo, ibl, ibleo, mu, ix, iy, iz, &
                                      it, icri, ierr, id
    integer(kind=KI), dimension(4) :: agam5, np, ip, srcbit
    real(kind=KR),    dimension(4) :: bgam5
    real(kind=KR)                  :: pcorr, psum

! Some initializations.
    agam5(1:4) = (/ 3, 4, 1, 2 /) 
    bgam5(1:4) = (/ 1.0_KR, 1.0_KR, -1.0_KR, -1.0_KR /) 
    np(1) = npx
    np(2) = npy
    np(3) = npz
    np(4) = npt
    call atoc(myid,np,ip)

! srcbit is the position of the source on this sublattice.
    srcbit(1) = src(1) - ip(1)*nx/npx
    srcbit(2) = src(2) - ip(2)*ny/npy
    srcbit(3) = src(3) - ip(3)*nz/npz
    srcbit(4) = src(4) - ip(4)*nt/npt

! There is no contribution unless the point source is on this process.
    if (srcbit(1)>0.and.srcbit(1)<nx/npx+1.and. &
        srcbit(2)>0.and.srcbit(2)<ny/npy+1.and. &
        srcbit(3)>0.and.srcbit(3)<nz/npz+1.and. &
        srcbit(4)>0.and.srcbit(4)<nt/npt+1) then
!-determine the isite of the point source.
     ix = (srcbit(1)-1)/4
     iy = (srcbit(2)-1)/2
     iz = (srcbit(3)-1)/2
     it = (srcbit(4)-1)/2
     isite = 1 + ix + iy*nx/(4*npx) + iz*nx*ny/(8*npx*npy) &
           + it*nx*ny*nz/(16*npx*npy*npz)
!-determine the ieo and ibl of the point source.
     ieo = 1
     ibl = 1
     do mu = 1,4
      if (modulo(srcbit(mu),4)==3.or.modulo(srcbit(mu),4)==0) ieo = 3 - ieo
      if (modulo(srcbit(mu),2)==0) ibl = iblv(ibl,mu)
     enddo ! mu
!-determine the ibleo of the point source.
     ibleo = vecblinv(2,ibl)
!-the pion correlator is built by putting gamma5 at the source.
     icri = 2*icsrc
     id = agam5(idsrc)
     if (vecblinv(1,ibl)==1) then
      pcorr = bgam5(idsrc)*xe(icri,isite,id,ieo,ibleo)
     else
      pcorr = bgam5(idsrc)*xo(icri,isite,id,ieo,ibleo)
     endif
    endif
    if (nps==1) then
     pion = pion + pcorr
    else
     call MPI_REDUCE(pcorr,psum,1,MRT,MPI_SUM,0,MPI_COMM_WORLD,ierr)
     call MPI_BCAST(psum,1,MRT,0,MPI_COMM_WORLD,ierr)
     pion = pion + psum
    endif

 end subroutine SSTpion

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 subroutine mesoncorr(nkappa,ntmqcd,nSST,src,nsmrsrc,nsmrsnk,asmrsnk,bc, &
                      rwdir,cfgfile,u,vecbl,vecblinv,myid,nn,ldiv,nms,lvbc, &
                      ib,lbd,iblv,MRT,be,bo,xe,xo,icfgsave)
! Compute two-point correlation functions for operators of type qbar*Gamma*q.
! In this subroutine, scalar means psibar 1 psi,
!                     pseudo means psibar gamma_5 psi,
!                     vector means psibar gamma_mu psi,
!                     axial  means psibar gamma_mu gamma_5 psi,
!                     tensor means psibar i*sigma_(mu,nu) psi,
! where i=sqrt(-1).
! Three-point functions are computed by this same code if sequential source
! propagators are read from disk.
! 
! INPUT:
!   nkappa is the number of hopping parameters for the unimproved Wilson action.
!   ntmqcd is the number of (kappa,mu) pairs for the Wilson-clover action.
!   nSST(): abs(nSST(1)*nSST(2))=number of sequential source propagator pairs.
!           If nSST(1)>0, nSST(2)>0, then tmUtmU SST prop's are used.
!           If nSST(1)>0, nSST(2)<0, then tmUtmU,tmUtmD SST prop's are used.
!           If nSST(1)<0, nSST(2)>0, then tmUtmU,tmDtmU SST prop's are used.
!           If nSST(1)<0, nSST(2)<0, then tmUtmU,tmUtmD,tmDtmU,tmDtmD are used.
!   src(:) is the location of the point source on the global lattice.
!   nsmrsrc is the user-defined smearing parameter for the source (only used
!           here to determine filenames).
!   nsmrsnk is the user-defined smearing parameter for the sink.
!   asmrsnk is the user-defined smearing coefficient for the sink.
!   bc(mu) is 1,-1,0 for periodic,antiperiodic,fixed boundary conditions
!          respectively.
!   rwdir(i) is the read/write directory for the i'th process.
!   cfgfile is the filename for propagators, omitting some trailing characters.
!   u() contains the gauge fields for this sublattice.
!   vecbl(1,:) are the blocks of the lattice where sites are globally "even".
!   vecbl(2,:) are the blocks of the lattice where sites are globally "odd".
!   vecblinv(1,1:16) = (/ 1, 2, 2, 1, 2, 1, 1, 2, 2, 1, 1, 2, 1, 2, 2, 1 /)
!   vecblinv(2,1:16) = (/ 1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7, 8, 8 /)
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
!   MRT is MPI_REAL or MPI_DOUBLE_PRECISION.  It must match KR.
!   be,bo,xe,xo are used as scratch space.  Their input values are irrelevant.
!   expected size: be(6,ntotal,4,2,8), bo(6,ntotal,4,2,8),
!                  xe(6,nvhalf,4,2,8), xo(6,nvhalf,4,2,8)
! OUTPUT:
!   be,bo,xe,xo are used as scratch space.  Their output values are irrelevant.

    integer(kind=KI), intent(in)                          :: nkappa, ntmqcd, &
                                                             nsmrsrc, &
                                                             nsmrsnk, myid, MRT
    integer(kind=KI), intent(in),    dimension(:)         :: nSST, bc, nms
    integer(kind=KI), intent(in),    dimension(:)         :: src
    real(kind=KR),    intent(in)                          :: asmrsnk
    character(len=*), intent(in),    dimension(:)         :: rwdir
    character(len=*), intent(in)                          :: cfgfile
    real(kind=KR),    intent(in),    dimension(:,:,:,:,:) :: u
    integer(kind=KI), intent(in),    dimension(:,:)       :: vecbl, vecblinv
    integer(kind=KI), intent(in),    dimension(:,:)       :: nn, iblv
    logical,          intent(in),    dimension(:)         :: ldiv
    integer(kind=KI), intent(in),    dimension(:,:,:)     :: lvbc
    integer(kind=KI), intent(in),    dimension(:,:,:,:)   :: ib
    logical,          intent(in),    dimension(:,:)       :: lbd
    real(kind=KR),    intent(inout), dimension(:,:,:,:,:) :: be, bo, xe, xo
    integer(kind=KI), intent(in)                          :: icfgsave

    integer(kind=KI), dimension(2)           :: ix
    character(len=1), dimension(2)           :: iLSsymbol, jLSsymbol
    character(len=14)                        :: trailer
    character(len=11)                        :: trailerq, trailerqbar
    character(len=8)                         :: smrq, smrqbar
    character(len=23)                        :: labelqqbar
    character(len=128)                       :: propfile, opfile
    character(len=3), dimension(16)          :: opname
    real(kind=KR2),   dimension(nt,-4:4,-4:4,-4:4,28) :: corrRebit, corrImbit
    real(kind=KR),    dimension(nt,-4:4,-4:4,-4:4,28) :: corrReKR, corrImKR, &
                                                         corrRe, corrIm
    real(kind=KR),    dimension(3,4,2,3,4,2) :: q, qbar
    integer(kind=KI), dimension(4)           :: np, ip
    real(kind=KR2),   parameter              :: twoPi=6.283185307179586_KR2
    real(kind=KR2),   dimension(3,-4:4)      :: mycos, mysin
    real(kind=KR2)                           :: bitRe, bitIm, bittemp
    integer(kind=KI), dimension(16)          :: rootm1
    integer(kind=KI), dimension(2)           :: iri2, iri3
    integer(kind=KI), dimension(16,4)        :: jd2
    real(kind=KR),    dimension(16,4)        :: sgni, sgnj
    real(kind=KR)                            :: kapin, muin, cSWin
    integer(kind=KI), dimension(30)          :: jop, iop
    logical                                  :: docorr
    integer(kind=KI) :: nLSsrc, nLSsnk, iLS, jLS, ncalcs, ikapq, ikapqbar, &
                        jc, jd, ieo, ieo1, ieo2, ixbit, iybit, izbit, itbit, &
                        ixbit2, iybit2, izbit2, itbit2, ixbit3, ibl, iblbit, &
                        iy, iz, it, ifile, imom, ixmom, iymom, izmom, iri, &
                        ic, id, lessnoise, icin, idin, icfgin, &
                        myidin, nmesg, jx, imax1, imax2, imax3, &
                        imax4, imax5, imax6, sgntmqcd, sgnsst1, sgnsst2, ii, &
                        jj, iimax, ierr, icount, ntmin, ntmax

    integer(kind=KI) :: lessnoiseoverride
    logical          :: fileexists, corrsp
    real :: time1, time2

! Define some parameters.
    if (nsmrsrc==0) then
     nLSsrc = 1
     jLSsymbol(1) = "L"
    elseif (nsmrsrc>0) then
     nLSsrc = 1
     jLSsymbol(1) = "S"
    else
     nLSsrc = 2
     jLSsymbol(1) = "L"
     jLSsymbol(2) = "S"
    endif
    if (nsmrsnk==0) then
     nLSsnk = 1
     iLSsymbol(1) = "L"
    elseif (nsmrsnk>0) then
     nLSsnk = 1
     iLSsymbol(1) = "S"
    else
     nLSsnk = 2
     iLSsymbol(1) = "L"
     iLSsymbol(2) = "S"
    endif
    if (ntmqcd<0) then
     sgntmqcd = 2
    else
     sgntmqcd = 1
    endif
    if (nSST(1)<0) then
     sgnSST1 = 2
    else
     sgnSST1 = 1
    endif
    if (nSST(2)<0) then
     sgnSST2 = 2
    else
     sgnSST2 = 1
    endif
    imax1 = nkappa
    imax2 = imax1 + abs(ntmqcd)
    imax3 = imax2 + (sgntmqcd-1)*abs(ntmqcd)
    imax4 = imax3 + abs(nSST(1)*nSST(2))
    imax5 = imax4 + (sgnSST2-1)*abs(nSST(1)*nSST(2))
    imax6 = imax5 + (sgnSST1-1)*abs(nSST(1)*nSST(2))
    ncalcs = imax6 + (sgnSST1-1)*(sgnSST2-1)*abs(nSST(1)*nSST(2))

! Define names for the Dirac structures (to label entries in the output file).
    opname( 1) = "sca"
    opname( 2) = "pse"
    opname( 3) = "ve1"
    opname( 4) = "ve2"
    opname( 5) = "ve3"
    opname( 6) = "ve4"
    opname( 7) = "ax1"
    opname( 8) = "ax2"
    opname( 9) = "ax3"
    opname(10) = "ax4"
    opname(11) = "t12"
    opname(12) = "t13"
    opname(13) = "t23"
    opname(14) = "t14"
    opname(15) = "t24"
    opname(16) = "t34"

! Define the 16 possible Dirac structures for the source, corresponding to
! gamma_4*(operator)^dagger*gamma_4*gamma_5.
! "rootm1" is -1 if the operator is purely imaginary, 1 is it is purely real.
! "jd2(:,i)" is the non-zero column in the i'th row of the Dirac matrix.
    rootm1( 1) = -1 ! scalar       needs -gamma_5
    rootm1( 2) =  1 ! pseudoscalar needs 1
    rootm1( 3) = -1 ! vector_1     needs gamma_1*gamma_5
    rootm1( 4) =  1 ! vector_2     needs gamma_2*gamma_5
    rootm1( 5) = -1 ! vector_3     needs gamma_3*gamma_5
    rootm1( 6) = -1 ! vector_4     needs -gamma_4*gamma_5
    rootm1( 7) =  1 ! axial_1      needs gamma_1
    rootm1( 8) = -1 ! axial_2      needs gamma_2
    rootm1( 9) =  1 ! axial_3      needs gamma_3
    rootm1(10) =  1 ! axial_4      needs -gamma_4
    rootm1(11) =  1 ! tensor_12    needs gamma_1*gamma_2*gamma_5
    rootm1(12) = -1 ! tensor_13    needs gamma_1*gamma_3*gamma_5
    rootm1(13) =  1 ! tensor_23    needs gamma_2*gamma_3*gamma_5
    rootm1(14) = -1 ! tensor_14    needs -gamma_1*gamma_4*gamma_5
    rootm1(15) =  1 ! tensor_24    needs -gamma_2*gamma_4*gamma_5
    rootm1(16) = -1 ! tensor_34    needs -gamma_3*gamma_4*gamma_5
    jd2( 1,1:4) = (/ 3, 4, 1, 2 /) ! scalar       needs -gamma_5
    jd2( 2,1:4) = (/ 1, 2, 3, 4 /) ! pseudoscalar needs 1
    jd2( 3,1:4) = (/ 2, 1, 4, 3 /) ! vector_1     needs gamma_1*gamma_5
    jd2( 4,1:4) = (/ 2, 1, 4, 3 /) ! vector_2     needs gamma_2*gamma_5
    jd2( 5,1:4) = (/ 1, 2, 3, 4 /) ! vector_3     needs gamma_3*gamma_5
    jd2( 6,1:4) = (/ 3, 4, 1, 2 /) ! vector_4     needs -gamma_4*gamma_5
    jd2( 7,1:4) = (/ 4, 3, 2, 1 /) ! axial_1      needs gamma_1
    jd2( 8,1:4) = (/ 4, 3, 2, 1 /) ! axial_2      needs gamma_2
    jd2( 9,1:4) = (/ 3, 4, 1, 2 /) ! axial_3      needs gamma_3
    jd2(10,1:4) = (/ 1, 2, 3, 4 /) ! axial_4      needs -gamma_4
    jd2(11,1:4) = (/ 3, 4, 1, 2 /) ! tensor_12    needs gamma_1*gamma_2*gamma_5
    jd2(12,1:4) = (/ 4, 3, 2, 1 /) ! tensor_13    needs gamma_1*gamma_3*gamma_5
    jd2(13,1:4) = (/ 4, 3, 2, 1 /) ! tensor_23    needs gamma_2*gamma_3*gamma_5
    jd2(14,1:4) = (/ 2, 1, 4, 3 /) ! tensor_14    needs -gamma_1*gamma_4*gamma_5
    jd2(15,1:4) = (/ 2, 1, 4, 3 /) ! tensor_24    needs -gamma_2*gamma_4*gamma_5
    jd2(16,1:4) = (/ 1, 2, 3, 4 /) ! tensor_34    needs -gamma_3*gamma_4*gamma_5
    sgnj( 1,1:4) = (/  1.0_KR,  1.0_KR, -1.0_KR, -1.0_KR /) ! scalar
    sgnj( 2,1:4) = (/  1.0_KR,  1.0_KR,  1.0_KR,  1.0_KR /) ! pseudoscalar
    sgnj( 3,1:4) = (/  1.0_KR,  1.0_KR, -1.0_KR, -1.0_KR /) ! vector_1
    sgnj( 4,1:4) = (/ -1.0_KR,  1.0_KR,  1.0_KR, -1.0_KR /) ! vector_2
    sgnj( 5,1:4) = (/  1.0_KR, -1.0_KR, -1.0_KR,  1.0_KR /) ! vector_3
    sgnj( 6,1:4) = (/  1.0_KR,  1.0_KR,  1.0_KR,  1.0_KR /) ! vector_4
    sgnj( 7,1:4) = (/  1.0_KR,  1.0_KR,  1.0_KR,  1.0_KR /) ! axial_1
    sgnj( 8,1:4) = (/ -1.0_KR,  1.0_KR, -1.0_KR,  1.0_KR /) ! axial_2
    sgnj( 9,1:4) = (/  1.0_KR, -1.0_KR,  1.0_KR, -1.0_KR /) ! axial_3
    sgnj(10,1:4) = (/ -1.0_KR, -1.0_KR,  1.0_KR,  1.0_KR /) ! axial_4
    sgnj(11,1:4) = (/ -1.0_KR,  1.0_KR,  1.0_KR, -1.0_KR /) ! tensor_12
    sgnj(12,1:4) = (/  1.0_KR, -1.0_KR, -1.0_KR,  1.0_KR /) ! tensor_13
    sgnj(13,1:4) = (/ -1.0_KR, -1.0_KR,  1.0_KR,  1.0_KR /) ! tensor_23
    sgnj(14,1:4) = (/  1.0_KR,  1.0_KR,  1.0_KR,  1.0_KR /) ! tensor_14
    sgnj(15,1:4) = (/ -1.0_KR,  1.0_KR, -1.0_KR,  1.0_KR /) ! tensor_24
    sgnj(16,1:4) = (/  1.0_KR, -1.0_KR,  1.0_KR, -1.0_KR /) ! tensor_34

! Define the 16 possible Dirac structures for the sink, corresponding to
! gamma_5*operator.  Note: rootm1(1:16) and id2(1:16,1:4) are not defined
! since they would be identical to jd2() defined above.
! scalar       needs gamma_5
! pseudoscalar needs 1
! vector_1     needs -gamma_1*gamma_5
! vector_2     needs -gamma_2*gamma_5
! vector_3     needs -gamma_3*gamma_5
! vector_4     needs -gamma_4*gamma_5
! axial_1      needs -gamma_1
! axial_2      needs -gamma_2
! axial_3      needs -gamma_3
! axial_4      needs -gamma_4
! tensor_12    needs -gamma_1*gamma_2*gamma_5
! tensor_13    needs -gamma_1*gamma_3*gamma_5
! tensor_23    needs -gamma_2*gamma_3*gamma_5
! tensor_14    needs -gamma_1*gamma_4*gamma_5
! tensor_24    needs -gamma_2*gamma_4*gamma_5
! tensor_34    needs -gamma_3*gamma_4*gamma_5
    sgni( 1,1:4) = (/ -1.0_KR, -1.0_KR,  1.0_KR,  1.0_KR /) ! scalar
    sgni( 2,1:4) = (/  1.0_KR,  1.0_KR,  1.0_KR,  1.0_KR /) ! pseudoscalar
    sgni( 3,1:4) = (/ -1.0_KR, -1.0_KR,  1.0_KR,  1.0_KR /) ! vector_1
    sgni( 4,1:4) = (/  1.0_KR, -1.0_KR, -1.0_KR,  1.0_KR /) ! vector_2
    sgni( 5,1:4) = (/ -1.0_KR,  1.0_KR,  1.0_KR, -1.0_KR /) ! vector_3
    sgni( 6,1:4) = (/  1.0_KR,  1.0_KR,  1.0_KR,  1.0_KR /) ! vector_4
    sgni( 7,1:4) = (/ -1.0_KR, -1.0_KR, -1.0_KR, -1.0_KR /) ! axial_1
    sgni( 8,1:4) = (/  1.0_KR, -1.0_KR,  1.0_KR, -1.0_KR /) ! axial_2
    sgni( 9,1:4) = (/ -1.0_KR,  1.0_KR, -1.0_KR,  1.0_KR /) ! axial_3
    sgni(10,1:4) = (/ -1.0_KR, -1.0_KR,  1.0_KR,  1.0_KR /) ! axial_4
    sgni(11,1:4) = (/  1.0_KR, -1.0_KR, -1.0_KR,  1.0_KR /) ! tensor_12
    sgni(12,1:4) = (/ -1.0_KR,  1.0_KR,  1.0_KR, -1.0_KR /) ! tensor_13
    sgni(13,1:4) = (/  1.0_KR,  1.0_KR, -1.0_KR, -1.0_KR /) ! tensor_23
    sgni(14,1:4) = (/  1.0_KR,  1.0_KR,  1.0_KR,  1.0_KR /) ! tensor_14
    sgni(15,1:4) = (/ -1.0_KR,  1.0_KR, -1.0_KR,  1.0_KR /) ! tensor_24
    sgni(16,1:4) = (/  1.0_KR, -1.0_KR,  1.0_KR, -1.0_KR /) ! tensor_34

! snk(i): An operator at lattice site (x,y,z) will actually be computed using
!         propagators going to site (x+snk(1),y+snk(2),z+snk(3)).  This only
!         affects the momentum factors and is useful for derivative operators.

! Define the correlators of interest that do not involve an SST propagator.
! S->S, P->P, V_mu->V_mu, A_mu->A_mu
    do icount = 1,10
     jop(icount) = icount
     iop(icount) = icount
    enddo ! icount

! WARNING :: if corrsp = .true., we will neglect the P-A_4 2pt function
!            to calcualte the scalar-pseudoscalar 2pt-function. This is done because
!            when doing TWdisconloops, it is necessary to investigate "contamination" of
!            the scalar/pseudoscalar when using max-twist defined by V-A. 

  corrsp = .true. 

! P->A_4
    jop(11) = 2
    iop(11) = 10

  if (corrsp) then 
! S->P
    jop(11) = 1
    iop(11) = 2
  endif !

! P->V_4
    jop(12) = 2
    iop(12) = 6
! V_mu->A_mu
    jop(13:16) = (/ 3, 4, 5, 6 /)
    iop(13:16) = (/ 7, 8, 9, 10 /)
! A_mu->V_mu
    jop(17:20) = (/ 7, 8, 9, 10 /)
    iop(17:20) = (/ 3, 4, 5, 6 /)
! V_4->P
    jop(21) = 6
    iop(21) = 2

! A_4->P
    jop(22) = 10
    iop(22) = 2

  if (corrsp) then
! S->P
    jop(22) = 2
    iop(22) = 1
  endif

! V_mu->T_(mu,4)
    jop(23:25) = (/ 3, 4, 5 /)
    iop(23:25) = (/ 14, 15, 16 /)
! A_mu->T_(mu,4)
    jop(26:28) = (/ 7, 8, 9 /)
    iop(26:28) = (/ 14, 15, 16 /)

! Define the correlators of interest that involve one SST propagator.
! P->V_4
    jop(29) = 2
    iop(29) = 6
! P->A_4
    jop(30) = 2
    iop(30) = 10

! Consider local and/or smeared sources, as requested through "nsmrsrc".
    do jLS = 1,nLSsrc

! Consider local and/or smeared sinks, as requested through "nsmrsnk".
     do iLS = 1,nLSsnk

! Consider each quark in turn.
      do ikapq = 1,ncalcs
       if (ikapq>imax6) then
        trailerq = ".xxk"//jLSsymbol(jLS)//"tmDtmD"
        write(unit=trailerq(4:4),fmt="(i1.1)") ikapq - imax6
       elseif (ikapq>imax5) then
        trailerq = ".xxk"//jLSsymbol(jLS)//"tmDtmU"
        write(unit=trailerq(4:4),fmt="(i1.1)") ikapq - imax5
       elseif (ikapq>imax4) then
        trailerq = ".xxk"//jLSsymbol(jLS)//"tmUtmD"
        write(unit=trailerq(4:4),fmt="(i1.1)") ikapq - imax4
       elseif (ikapq>imax3) then
        trailerq = ".xxk"//jLSsymbol(jLS)//"tmUtmU"
        write(unit=trailerq(4:4),fmt="(i1.1)") ikapq - imax3
       elseif (ikapq>imax2) then
        trailerq = ".xxk"//jLSsymbol(jLS)//"tmD"
        write(unit=trailerq(4:4),fmt="(i1.1)") ikapq - imax2
       elseif (ikapq>imax1) then
        trailerq = ".xxk"//jLSsymbol(jLS)//"tmU"
        write(unit=trailerq(4:4),fmt="(i1.1)") ikapq - imax1
       else
        trailerq = ".xxk"//jLSsymbol(jLS)
        write(unit=trailerq(4:4),fmt="(i1.1)") ikapq
       endif

! Smear the quark at the sink, if desired.
       if ( nsmrsnk>0 .or. (nsmrsnk<0.and.nLSsnk==2) ) then
        smrq = "tempq"
        do jc = 1,3
         write(unit=trailerq(2:2),fmt="(i1.1)") jc
         do jd = 1,4
          write(unit=trailerq(3:3),fmt="(i1.1)") jd
          propfile = trim(rwdir(myid+1))//trim(cfgfile)//trim(trailerq)
          call FPread(xe,xo,vecblinv,propfile,kapin,muin,cSWin,icin,idin, &
                      icfgin,myidin)
          call smear(be,bo,u,xe,xo,asmrsnk,nsmrsnk,bc,vecbl,vecblinv,myid,nn, &
                     ldiv,nms,lvbc,ib,lbd,iblv,MRT)
          propfile = trim(propfile)//trim(smrq)
          call FPwrite(be,bo,vecblinv,propfile,kapin,muin,cSWin,icin,idin, &
                       icfgin,myidin)
         enddo ! jd
        enddo ! jc
       else
        smrq = ""
       endif

! Consider each antiquark in turn, and decide whether to calculate for this
! q-qbar pair.
       do ikapqbar = 1,ncalcs
        docorr = .false.
        if (ikapqbar>imax6) then
         trailerqbar = ".xxk"//jLSsymbol(jLS)//"tmDtmD"
         write(unit=trailerqbar(4:4),fmt="(i1.1)") ikapqbar - imax6
         if ( ikapq>imax1 .and. ikapq<imax3 ) docorr=.true.
        elseif (ikapqbar>imax5) then
         trailerqbar = ".xxk"//jLSsymbol(jLS)//"tmDtmU"
         write(unit=trailerqbar(4:4),fmt="(i1.1)") ikapqbar - imax5
         if ( ikapq>imax1 .and. ikapq<imax3 ) docorr=.true.
        elseif (ikapqbar>imax4) then
         trailerqbar = ".xxk"//jLSsymbol(jLS)//"tmUtmD"
         write(unit=trailerqbar(4:4),fmt="(i1.1)") ikapqbar - imax4
         if ( ikapq>imax1 .and. ikapq<imax3 ) docorr=.true.
        elseif (ikapqbar>imax3) then
         trailerqbar = ".xxk"//jLSsymbol(jLS)//"tmUtmU"
         write(unit=trailerqbar(4:4),fmt="(i1.1)") ikapqbar - imax3
         if ( ikapq>imax1 .and. ikapq<imax3 ) docorr=.true.
        elseif (ikapqbar>imax2) then
         trailerqbar = ".xxk"//jLSsymbol(jLS)//"tmD"
         write(unit=trailerqbar(4:4),fmt="(i1.1)") ikapqbar - imax2
         if ( ikapq>imax1 ) docorr=.true.
        elseif (ikapqbar>imax1) then
         trailerqbar = ".xxk"//jLSsymbol(jLS)//"tmU"
         write(unit=trailerqbar(4:4),fmt="(i1.1)") ikapqbar - imax1
         if ( ikapq>imax1 ) docorr=.true.
        else
         trailerqbar = ".xxk"//jLSsymbol(jLS)
         write(unit=trailerqbar(4:4),fmt="(i1.1)") ikapqbar
         docorr=.true.
        endif
        if (docorr) then
         if (ikapq>imax3) then
          iimax = 2
         else
          iimax = 28
         endif

! Smear the antiquark at the sink, if desired.
         if ( nsmrsnk>0 .or. (nsmrsnk<0.and.nLSsnk==2) ) then
          smrqbar = "tempqbar"
          do jc = 1,3
           write(unit=trailerqbar(2:2),fmt="(i1.1)") jc
           do jd = 1,4
            write(unit=trailerqbar(3:3),fmt="(i1.1)") jd
            propfile = trim(rwdir(myid+1))//trim(cfgfile)//trim(trailerqbar)
            call FPread(xe,xo,vecblinv,propfile,kapin,muin,cSWin,icin,idin, &
                        icfgin,myidin)
            call smear(be,bo,u,xe,xo,asmrsnk,nsmrsnk,bc,vecbl,vecblinv,myid, &
                       nn,ldiv,nms,lvbc,ib,lbd,iblv,MRT)
            propfile = trim(propfile)//trim(smrqbar)
            call FPwrite(be,bo,vecblinv,propfile,kapin,muin,cSWin,icin,idin, &
                         icfgin,myidin)
           enddo ! jd
          enddo ! jc
         else
          smrqbar = ""
         endif

! Open the quark and antiquark files.
         do jc = 1,3
          write(unit=trailerq(2:2),fmt="(i1.1)") jc
          write(unit=trailerqbar(2:2),fmt="(i1.1)") jc
          do jd = 1,4
           write(unit=trailerq(3:3),fmt="(i1.1)") jd
           write(unit=trailerqbar(3:3),fmt="(i1.1)") jd ! nt check
           propfile = trim(rwdir(myid+1))//trim(cfgfile)//trim(trailerq)// &
                      trim(smrq)
           ifile = 8 + jd + (jc-1)*4
           open(unit=ifile,file=trim(propfile),action="read",status="old", &
                position="rewind",form="unformatted")
           if (trailerqbar/=trailerq .or. smrqbar/=smrq) then
            propfile = trim(rwdir(myid+1))//trim(cfgfile)//trim(trailerqbar)// &
                       trim(smrqbar)
            ifile = 20 + jd + (jc-1)*4
            open(unit=ifile,file=trim(propfile),action="read",status="old", &
                 position="rewind",form="unformatted")
           endif
          enddo ! jd
         enddo ! jc

! Some initializations.
         corrRebit = 0.0_KR2
         corrImbit = 0.0_KR2
         np(1) = npx
         np(2) = npy
         np(3) = npz
         np(4) = npt
         call atoc(myid,np,ip)

! NOTICE: If the boundary conditions are periodic in ALL spatial directions,
!         then the terms which are odd in sinx, siny, sinz, sin2x, sin2y or
!         sin2z are omitted, since the integral of an even periodic function,
!         times sinx or sin2x etc., vanishes.  This should reduce the
!         statistical noise.

         lessnoiseoverride=1

         if (bc(1)==1.and.bc(2)==1.and.bc(3)==1.and.lessnoiseoverride==0) then
          lessnoise = 1
         else
          lessnoise = 0
         endif

        !call cpu_time(time1)
! Begin main loop, constructing ix,iy,iz,it from isite,ieo,ibl.
         ieo1 = 2
         ieo2 = 1
         do itbit = 2,nt/npt,2
          itbit2 = itbit + ip(4)*nt/npt
          ieo1 = 3 - ieo1
          ieo2 = 3 - ieo2
          do izbit = 2,nz/npz,2
           izbit2 = izbit + ip(3)*nz/npz
           ieo1 = 3 - ieo1
           ieo2 = 3 - ieo2
           do iybit = 2,ny/npy,2
            iybit2 = iybit + ip(2)*ny/npy
            ieo1 = 3 - ieo1
            ieo2 = 3 - ieo2
            do ixbit = 4,nx/npx,4
             ixbit2 = ixbit + ip(1)*nx/npx
             do ibl = 1,16
              if (ibl>8) then
               it = itbit2
              else
               it = itbit2 - 1
              endif
              iblbit = 1 + modulo(ibl-1,8)
              if (iblbit>4) then
               iz = izbit2
              else
               iz = izbit2 - 1
              endif
              iblbit = 1 + modulo(iblbit-1,4)
              if (iblbit>2) then
               iy = iybit2
              else
               iy = iybit2 - 1
              endif
              if (modulo(ibl,2)==1) then
               ixbit3 = ixbit2 - 1
              else
               ixbit3 = ixbit2
              endif
              ix(ieo1) = ixbit3 - 2
              ix(ieo2) = ixbit3

! Read in the 12 propagators for this pair of lattice sites (ieo=1 and ieo=2).
              do ic = 1,3
               do iri = 1,2
                do jc = 1,3
                 do jd = 1,4
                  ifile = 8 + jd + (jc-1)*4
                  read(unit=ifile) &
                      (q(jc,jd,iri,ic,id,1),q(jc,jd,iri,ic,id,2),id=1,4)
                  if (trailerqbar/=trailerq .or. smrqbar/=smrq) then
                   ifile = ifile + 12
                   read(unit=ifile) &
                       (qbar(jc,jd,iri,ic,id,1),qbar(jc,jd,iri,ic,id,2),id=1,4)
                   else
                    qbar(jc,jd,iri,ic,:,:) = q(jc,jd,iri,ic,:,:)
                  endif
                 enddo ! jd
                enddo ! jc
               enddo ! iri
              enddo ! ic

! Construct the correlators.
              do ieo = 1,2
               do id = 1,4
                do jd = 1,4
                 if (iimax==2) then
                  jj = 28
                 else
                  jj = 0
                 endif
                 do ii = 1,iimax
                  jj = jj + 1
                  if (jop(jj)<=2 .and.  iop(jj)<=2) then ! (P->P, S->P correlations)
                     if (jop(jj) /=1 .or. iop(jj)/=1) then ! (omit S->S correlations)
                  do imom = -4,4 
                   mycos(1,imom) = cos(twoPi*real(imom*(ix(ieo)-src(1)),KR2) &
                                   /real(nx,KR2))
                   mysin(1,imom) = sin(twoPi*real(imom*(ix(ieo)-src(1)),KR2) &
                                   /real(nx,KR2))
                   mycos(2,imom) = cos(twoPi*real(imom*(iy-src(2)),KR2) &
                                   /real(ny,KR2))
                   mysin(2,imom) = sin(twoPi*real(imom*(iy-src(2)),KR2) &
                                   /real(ny,KR2))
                   mycos(3,imom) = cos(twoPi*real(imom*(iz-src(3)),KR2) &
                                   /real(nz,KR2))
                   mysin(3,imom) = sin(twoPi*real(imom*(iz-src(3)),KR2) &
                                   /real(nz,KR2))
                  enddo ! imom
                  if (rootm1(jop(jj))*rootm1(iop(jj))==1) then
                   iri2 = (/ 1, 2 /)
                   iri3 = (/ 2, 1 /)
                  else
                   iri2 = (/ 2, 1 /)
                   iri3 = (/ 1, 2 /)
                  endif
                  bitRe = 0.0_KR2
                  bitIm = 0.0_KR2
                  do ic = 1,3
                   do jc = 1,3
                    bitRe = bitRe &
                          + real(q(jc,jd,1,ic,jd2(iop(jj),id),ieo) &
                                *qbar(jc,jd2(jop(jj),jd),1,ic,id,ieo) &
                                +q(jc,jd,2,ic,jd2(iop(jj),id),ieo) &
                                *qbar(jc,jd2(jop(jj),jd),2,ic,id,ieo),KR2)
                    bitIm = bitIm &
                          + real(q(jc,jd,2,ic,jd2(iop(jj),id),ieo) &
                                *qbar(jc,jd2(jop(jj),jd),1,ic,id,ieo) &
                                -q(jc,jd,1,ic,jd2(iop(jj),id),ieo) &
                                 *qbar(jc,jd2(jop(jj),jd),2,ic,id,ieo),KR2)
                   enddo ! jc
                  enddo ! ic
                  if (rootm1(jop(jj))==1 .and. rootm1(iop(jj))==1) then
                   bitRe = bitRe*sgnj(jop(jj),jd)*sgni(iop(jj),id)
                   bitIm = bitIm*sgnj(jop(jj),jd)*sgni(iop(jj),id)
                  elseif (rootm1(jop(jj))==-1 .and. rootm1(iop(jj))==-1) then
                   bitRe = -bitRe*sgnj(jop(jj),jd)*sgni(iop(jj),id)
                   bitIm = -bitIm*sgnj(jop(jj),jd)*sgni(iop(jj),id)
                  else
                   bittemp = bitIm*sgnj(jop(jj),jd)*sgni(iop(jj),id)
                   bitIm = bitRe*sgnj(jop(jj),jd)*sgni(iop(jj),id)
                   bitRe = bittemp
                  endif
                  do izmom = -4,4
                   do iymom = -4,4
                    do ixmom = -4,4
                     corrRebit(it,ixmom,iymom,izmom,ii) &
                          = corrRebit(it,ixmom,iymom,izmom,ii) &
                          + bitRe*mycos(1,ixmom)*mycos(2,iymom)*mycos(3,izmom)
                     corrImbit(it,ixmom,iymom,izmom,ii) &
                          = corrImbit(it,ixmom,iymom,izmom,ii) &
                          + bitIm*mycos(1,ixmom)*mycos(2,iymom)*mycos(3,izmom)
                    enddo ! ixmom
                   enddo ! iymom
                  enddo ! izmom
                  if (lessnoise==0) then
                   do izmom = -4,4
                    do iymom = -4,4
                     do ixmom = -4,4
                      corrRebit(it,ixmom,iymom,izmom,ii) &
                        = corrRebit(it,ixmom,iymom,izmom,ii) &
                        - bitIm*(mycos(1,ixmom)*mycos(2,iymom)*mysin(3,izmom)&
                                +mycos(1,ixmom)*mysin(2,iymom)*mycos(3,izmom)&
                                +mysin(1,ixmom)*mycos(2,iymom)*mycos(3,izmom))&
                        - bitRe*(mysin(1,ixmom)*mysin(2,iymom)*mycos(3,izmom)&
                                +mysin(1,ixmom)*mycos(2,iymom)*mysin(3,izmom)&
                                +mycos(1,ixmom)*mysin(2,iymom)*mysin(3,izmom))&
                        + bitIm*mysin(1,ixmom)*mysin(2,iymom)*mysin(3,izmom)
                      corrImbit(it,ixmom,iymom,izmom,ii) &
                        = corrImbit(it,ixmom,iymom,izmom,ii) &
                        + bitRe*(mycos(1,ixmom)*mycos(2,iymom)*mysin(3,izmom)&
                                +mycos(1,ixmom)*mysin(2,iymom)*mycos(3,izmom)&
                                +mysin(1,ixmom)*mycos(2,iymom)*mycos(3,izmom))&
                        - bitIm*(mysin(1,ixmom)*mysin(2,iymom)*mycos(3,izmom)&
                                +mysin(1,ixmom)*mycos(2,iymom)*mysin(3,izmom)&
                                +mycos(1,ixmom)*mysin(2,iymom)*mysin(3,izmom))&
                        - bitRe*mysin(1,ixmom)*mysin(2,iymom)*mysin(3,izmom)
                     enddo ! ixmom
                    enddo ! iymom
                   enddo ! izmom
                  endif
                    endif ! omit scalar
                  endif ! check Kaon
                 enddo ! ii
                enddo ! jd
               enddo ! id
              enddo ! ieo

! End main loop.
             enddo ! ibl
            enddo ! ixbit
           enddo ! iybit
          enddo ! izbit
         enddo ! itbit

! Close the quark and antiquark files.
         do jc = 1,3
          do jd = 1,4
           ifile = 8 + jd + (jc-1)*4
           if ( (nsmrsnk>0 .or. (nsmrsnk<0.and.nLSsnk==2)) ) then
            close(unit=ifile,status="delete")
            if (trailerqbar/=trailerq .or. smrqbar/=smrq) then
             ifile = ifile + 12
             close(unit=ifile,status="delete")
            endif
           else
            close(unit=ifile,status="keep")
            if (trailerqbar/=trailerq .or. smrqbar/=smrq) then
             ifile = ifile + 12
             close(unit=ifile,status="keep")
            endif
           endif
          enddo ! jd
         enddo ! jc

! Combine the results from all processes.
         corrReKR = real(corrRebit,KR)
         corrImKR = real(corrImbit,KR)
         if (nps==1) then
          corrRe = corrReKR
          corrIm = corrImKR
         else
          nmesg = 729*nt*iimax
          call MPI_REDUCE(corrReKR(1,-4,-4,-4,1),corrRe(1,-4,-4,-4,1),nmesg, &
                          MRT,MPI_SUM,0,MPI_COMM_WORLD,ierr)
          call MPI_BCAST(corrRe(1,-4,-4,-4,1),nmesg,MRT,0,MPI_COMM_WORLD,ierr)
          call MPI_REDUCE(corrImKR(1,-4,-4,-4,1),corrIm(1,-4,-4,-4,1),nmesg, &
                          MRT,MPI_SUM,0,MPI_COMM_WORLD,ierr)
          call MPI_BCAST(corrIm(1,-4,-4,-4,1),nmesg,MRT,0,MPI_COMM_WORLD,ierr)
         endif

! Record the results.
         if (myid==0) then
          labelqqbar = "q=nxxxyyy qbar=nxxxyyy "
          if (ikapq>imax6) then 
           write(unit=labelqqbar(3:3),fmt="(i1)") ikapq-imax6
           write(unit=labelqqbar(4:9),fmt="(a6)") "tmDtmD"
          elseif (ikapq>imax5) then 
           write(unit=labelqqbar(3:3),fmt="(i1)") ikapq-imax5
           write(unit=labelqqbar(4:9),fmt="(a6)") "tmDtmU"
          elseif (ikapq>imax4) then 
           write(unit=labelqqbar(3:3),fmt="(i1)") ikapq-imax4
           write(unit=labelqqbar(4:9),fmt="(a6)") "tmUtmD"
          elseif (ikapq>imax3) then 
           write(unit=labelqqbar(3:3),fmt="(i1)") ikapq-imax3
           write(unit=labelqqbar(4:9),fmt="(a6)") "tmUtmU"
          elseif (ikapq>imax2) then 
           write(unit=labelqqbar(3:3),fmt="(i1)") ikapq-imax2
           write(unit=labelqqbar(4:9),fmt="(a6)") "tmD   "
          elseif (ikapq>imax1) then 
           write(unit=labelqqbar(3:3),fmt="(i1)") ikapq-imax1
           write(unit=labelqqbar(4:9),fmt="(a6)") "tmU   "
          else
           write(unit=labelqqbar(3:3),fmt="(i1)") ikapq
           write(unit=labelqqbar(4:9),fmt="(a6)") "Wilson"
          endif
          if (ikapqbar>imax6) then 
           write(unit=labelqqbar(16:16),fmt="(i1)") ikapqbar-imax6
           write(unit=labelqqbar(17:22),fmt="(a6)") "tmDtmD"
          elseif (ikapqbar>imax5) then 
           write(unit=labelqqbar(16:16),fmt="(i1)") ikapqbar-imax5
           write(unit=labelqqbar(17:22),fmt="(a6)") "tmDtmU"
          elseif (ikapqbar>imax4) then 
           write(unit=labelqqbar(16:16),fmt="(i1)") ikapqbar-imax4
           write(unit=labelqqbar(17:22),fmt="(a6)") "tmUtmD"
          elseif (ikapqbar>imax3) then 
           write(unit=labelqqbar(16:16),fmt="(i1)") ikapqbar-imax3
           write(unit=labelqqbar(17:22),fmt="(a6)") "tmUtmU"
          elseif (ikapqbar>imax2) then 
           write(unit=labelqqbar(16:16),fmt="(i1)") ikapqbar-imax2
           write(unit=labelqqbar(17:22),fmt="(a6)") "tmD   "
          elseif (ikapqbar>imax1) then 
           write(unit=labelqqbar(16:16),fmt="(i1)") ikapqbar-imax1
           write(unit=labelqqbar(17:22),fmt="(a6)") "tmU   "
          else
           write(unit=labelqqbar(16:16),fmt="(i1)") ikapqbar
           write(unit=labelqqbar(17:22),fmt="(a6)") "Wilson"
          endif

          trailer = "xxxx.MESON.LOG"
          write(unit=trailer(1:4), fmt = "(1i4.4)") icfgsave
          opfile = trim(rwdir(myid+1))//trim(trailer)
          inquire(file=opfile, exist=fileexists)
          if (.not. fileexists) then
             open(unit=10,file=opfile, &
                  action="write",form="formatted",status="new")
          else
             open(unit=10,file=opfile, &
                  action="write",form="formatted",status="old",position="append")
          endif ! file check

           do iz = -4,4
            do iy = -4,4
             do jx = -4,4
              if (iimax==2) then
               jj = 28
              else
               jj = 0
              endif
              do ii = 1,iimax
               jj = jj + 1
               if (jop(jj) <= 2 .and. iop(jj) <= 2) then ! (P->P, S->P correlations)
                  if (jop(jj) /=1 .or. iop(jj)/=1) then ! (omit S->S correlations)
                  if (nt<49) then
                     write(unit=10, &
                           fmt="(a38,a1,a3,a1,a3,a8,i2,a1,i2,a1,i2,a5,48es17.10)") &
                           labelqqbar, jLSsymbol(jLS), opname(jop(jj)), &
                           iLSsymbol(iLS), opname(iop(jj)), &
                           " Re [p=(",jx,",",iy,",",iz,")] = ", corrRe(:,jx,iy,iz,ii)
                     write(unit=10, &
                           fmt="(a38,a1,a3,a1,a3,a8,i2,a1,i2,a1,i2,a5,48es17.10)") &
                           labelqqbar, jLSsymbol(jLS), opname(jop(jj)), &
                           iLSsymbol(iLS), opname(iop(jj)), &
                           " Im [p=(",jx,",",iy,",",iz,")] = ", corrIm(:,jx,iy,iz,ii)
                  else
                     ntmax = 48
                     ntmin = 49
                     write(unit=10, &
                           fmt="(a38,a1,a3,a1,a3,a8,i2,a1,i2,a1,i2,a5,48es17.10)") &
                           labelqqbar, jLSsymbol(jLS), opname(jop(jj)), &
                           iLSsymbol(iLS), opname(iop(jj)), &
                           " ReA[p=(",jx,",",iy,",",iz,")] = ", &
                           corrRe(1:ntmax,jx,iy,iz,ii)
                     write(unit=10, &
                           fmt="(a38,a1,a3,a1,a3,a8,i2,a1,i2,a1,i2,a5,48es17.10)") &
                           labelqqbar, jLSsymbol(jLS), opname(jop(jj)), &
                           iLSsymbol(iLS), opname(iop(jj)), &
                           " ReB[p=(",jx,",",iy,",",iz,")] = ", &
                           corrRe(ntmin:,jx,iy,iz,ii)
                     write(unit=10, &
                           fmt="(a38,a1,a3,a1,a3,a8,i2,a1,i2,a1,i2,a5,48es17.10)") &
                           labelqqbar, jLSsymbol(jLS), opname(jop(jj)), &
                           iLSsymbol(iLS), opname(iop(jj)), &
                           " ImA[p=(",jx,",",iy,",",iz,")] = ", &
                           corrIm(1:ntmax,jx,iy,iz,ii)
                     write(unit=10, &
                           fmt="(a38,a1,a3,a1,a3,a8,i2,a1,i2,a1,i2,a5,48es17.10)") &
                           labelqqbar, jLSsymbol(jLS), opname(jop(jj)), &
                           iLSsymbol(iLS), opname(iop(jj)), &
                           " ImB[p=(",jx,",",iy,",",iz,")] = ", &
                           corrIm(ntmin:,jx,iy,iz,ii)
                endif ! nt check
                endif ! omit scalar
               endif ! check Kaon
              enddo ! ii
             enddo ! jx
            enddo ! iy
           enddo ! iz
          close(unit=10,status="keep")
         endif
        endif
        !call cpu_time(time2)
        !if (myid==0) print *, "Time to do mesoncorr=", time2-time1
        !call printlog("Stopping in MESONCORR!!", myid,rwdir)
        !call MPI_BARRIER(MPI_COMM_WORLD,ierr)
        !stop
       enddo ! ikapqbar
      enddo ! ikapq
     enddo ! iLS
    enddo ! jLS

 end subroutine mesoncorr

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 subroutine corr(nkappa,ntmqcd,nsmear,src,bc,myid,rwdir,cfgfile,MRT)
!
! This subroutine will eventually get replaced by subroutine baryoncorr,
! but since that has not yet been written this old code is retained.
!
! WARNING! WARNING! WARNING! WARNING! WARNING! WARNING! WARNING! WARNING!
! This subroutine is only correct for cSW=0.0_KR
! WARNING! WARNING! WARNING! WARNING! WARNING! WARNING! WARNING! WARNING!
!
! Compute the pion, rho and nucleon correlators (Re and Im parts and all
! Dirac components) using local operators.
!
! This subroutine uses N = epsilon^(abc) q^a [ (q^b)^T C*gamma5 q^c ]
!
!                  / 0 -1  0 0 \
! where C*gamma5 = | 1  0  0 0 |
!                  | 0  0  0 1 |
!                  \ 0  0 -1 0 /
!
! The Dirac index for the creation operator is set equal to the Dirac index
! for the annihilation operator.  (All four components are computed.)

    integer(kind=KI), intent(in)               :: nkappa, ntmqcd, nsmear, &
                                                   myid, MRT
    integer(kind=KI), intent(in), dimension(:) :: src, bc
    character(len=*), intent(in), dimension(:) :: rwdir
    character(len=*), intent(in)               :: cfgfile

    character(len=128)                       :: propfile
    character(len=8)                         :: trailer
    integer(kind=KI), dimension(4)           :: np, ip, cgam5, dsgn
    integer(kind=KI), dimension(6,3)         :: col
    real(kind=KR2),   parameter              :: twoPi=6.283185307179586_KR2
    real(kind=KR2),   dimension(6)           :: csgn
    real(kind=KR2)                           :: bit3, bit4
    real(kind=KR2),   dimension(4,4)         :: bitRe, bitIm, bit1Re, bit1Im, &
                                                bit2Re, bit2Im
    real(kind=KR),    dimension(3,4,2,3,4,2) :: q
    integer(kind=KI), dimension(2)           :: ix
    real(kind=KR2),   dimension(3,0:4)       :: mycos, mysin
    real(kind=KR2),   dimension(nt,0:4,0:4,0:4)       :: pbit, rbit
    real(kind=KR),    dimension(nt,0:4,0:4,0:4)       :: pcorr, rcorr, pion, rho
    real(kind=KR2),   dimension(2,nt,4,4,0:2,0:2,0:2) :: nbit
    real(kind=KR),    dimension(2,nt,4,4,0:2,0:2,0:2) :: ncorr, nucl
    character(len=1), dimension(2)                    :: LSsymbol
    integer(kind=KI) :: jc, jd, ifile, itbit, jx, iy, iz, it, ibl, ic, id, &
                        iri, ieo, ieo1, ieo2, icol, jcol, ca, cb, cc, cd, ce, &
                        cf, da, db, dc, dd, nmesg, ierr, itbit2, iblbit, &
                        ixbit, iybit, izbit, ixbit2, iybit2, izbit2, ixbit3, &
                        lessnoise, ikappa, nLS, iLS, ncalcs, imom, ixmom, &
                        iymom, izmom

    real(kind=KR),   dimension(nt,2)                  :: znuc
    integer(kind=KI)  :: ittt

! Compute the nucleon correlator from a local operator.
! NOTE: it=1 is uniquely defined by          1<=isite<=  factor and 1<=ibl<= 8
!       it=2 is uniquely defined by          1<=isite<=  factor and 9<=ibl<=16
!       it=3 is uniquely defined by   factor+1<=isite<=2*factor and 1<=ibl<= 8
!       it=4 is uniquely defined by   factor+1<=isite<=2*factor and 9<=ibl<=16
!       it=5 is uniquely defined by 2*factor+1<=isite<=3*factor and 1<=ibl<= 8
!       it=6 is uniquely defined by 2*factor+1<=isite<=3*factor and 9<=ibl<=16
!       etc, where factor = 2*nvhalf*npt/nt.

! Some parameters.
     if (nsmear==0) then
      nLS = 1
      LSsymbol(1) = "L"
     elseif (nsmear>0) then
      nLS = 1
      LSsymbol(1) = "S"
     else
      nLS = 2
      LSsymbol(1) = "L"
      LSsymbol(2) = "S"
     endif
     if (ntmqcd>0) then
      ncalcs = nkappa + ntmqcd
     else
      ncalcs = nkappa + 2*abs(ntmqcd)
     endif

! Consider multiple quark masses SEPARATELY.  (Hopefully I can mix them soon.)
    do ikappa = 1,ncalcs
    do iLS = 1,nLS
       if (ikappa<=nkappa) then
        trailer = ".xxk"//LSsymbol(iLS)
        write(unit=trailer(4:4),fmt="(i1.1)") ikappa
       elseif (ikappa<=nkappa+abs(ntmqcd)) then
        trailer = ".xxk"//LSsymbol(iLS)//"tmU"
        write(unit=trailer(4:4),fmt="(i1.1)") ikappa - nkappa
       else
        trailer = ".xxk"//LSsymbol(iLS)//"tmD"
        write(unit=trailer(4:4),fmt="(i1.1)") ikappa - nkappa - abs(ntmqcd)
       endif

! Open the necessary files.
     do jc = 1,3
      do jd = 1,4
       write(unit=trailer(2:2),fmt="(i1.1)") jc
       write(unit=trailer(3:3),fmt="(i1.1)") jd
       propfile = trim(rwdir(myid+1))//trim(cfgfile)//trim(trailer)
       ifile = 100 + jd + (jc-1)*4
       open(unit=ifile,file=propfile,action="read",status="old", &
            position="rewind",form="unformatted")
      enddo ! jd
     enddo ! jc
 
! Some initializations.
     pbit = 0.0_KR2
     rbit = 0.0_KR2
     nbit = 0.0_KR2
     col(1,1:3) = (/ 1, 2, 3 /)
     col(2,1:3) = (/ 2, 3, 1 /)
     col(3,1:3) = (/ 3, 1, 2 /)
     col(4,1:3) = (/ 1, 3, 2 /)
     col(5,1:3) = (/ 2, 1, 3 /)
     col(6,1:3) = (/ 3, 2, 1 /)
     cgam5(1:4) = (/ 2, 1, 4, 3 /)
     csgn(1:6) = (/ 1.0_KR2, 1.0_KR2, 1.0_KR2, -1.0_KR2, -1.0_KR2, -1.0_KR2 /)
     dsgn(1:4) = (/ 1.0_KR2, -1.0_KR2, -1.0_KR2, 1.0_KR2 /)
     np(1) = npx
     np(2) = npy
     np(3) = npz
     np(4) = npt
     call atoc(myid,np,ip)
 
! NOTICE: If the boundary conditions are periodic in ALL spatial directions,
!         then the terms which are odd in sinx, siny, sinz, sin2x, sin2y or
!         sin2z are omitted, since the integral of an even periodic function,
!         times sinx or sin2x etc., vanishes.  This should reduce the
!         statistical noise.
     if (bc(1)==1.and.bc(2)==1.and.bc(3)==1) then
      lessnoise = 1
     else
      lessnoise = 0
     endif
 
! Begin main loop, constructing ix,iy,iz,it from isite,ieo,ibl.
     ieo1 = 2
     ieo2 = 1
     do itbit = 2,nt/npt,2
      itbit2 = itbit + ip(4)*nt/npt
      ieo1 = 3 - ieo1
      ieo2 = 3 - ieo2
      do izbit = 2,nz/npz,2
       izbit2 = izbit + ip(3)*nz/npz
       ieo1 = 3 - ieo1
       ieo2 = 3 - ieo2
       do iybit = 2,ny/npy,2
        iybit2 = iybit + ip(2)*ny/npy
        ieo1 = 3 - ieo1
        ieo2 = 3 - ieo2
        do ixbit = 4,nx/npx,4
         ixbit2 = ixbit + ip(1)*nx/npx
         do ibl = 1,16
          if (ibl>8) then
           it = itbit2
          else
           it = itbit2 - 1
          endif
          iblbit = 1 + modulo(ibl-1,8)
          if (iblbit>4) then
           iz = izbit2
          else
           iz = izbit2 - 1
          endif
          iblbit = 1 + modulo(iblbit-1,4)
          if (iblbit>2) then
           iy = iybit2
          else
           iy = iybit2 - 1
          endif
          if (modulo(ibl,2)==1) then
           ixbit3 = ixbit2 - 1
          else
           ixbit3 = ixbit2
          endif
          ix(ieo1) = ixbit3 - 2
          ix(ieo2) = ixbit3
   
! Factors for the smallest momenta in the y and z directions.
          do imom = 0,4
           mycos(2,imom) = cos(twoPi*real(imom*(iy-src(2)),KR2)/real(ny,KR2))
           mysin(2,imom) = sin(twoPi*real(imom*(iy-src(2)),KR2)/real(ny,KR2))
           mycos(3,imom) = cos(twoPi*real(imom*(iz-src(3)),KR2)/real(nz,KR2))
           mysin(3,imom) = sin(twoPi*real(imom*(iz-src(3)),KR2)/real(nz,KR2))
          enddo ! imom
   
! Read in the 12 propagators for this pair of lattice sites (ieo=1 and ieo=2).
          do ic = 1,3
           do iri = 1,2
            do jc = 1,3
             do jd = 1,4
              ifile = 100 + jd + (jc-1)*4
              read(unit=ifile)(q(jc,jd,iri,ic,id,1),q(jc,jd,iri,ic,id,2),id=1,4)
             enddo ! jd
            enddo ! jc
           enddo ! iri
          enddo ! ic
   
! Construct the pion and rho correlators.  (rho uses gamma_3 only.)
          do ieo = 1,2
           do imom = 0,4
            mycos(1,imom) = cos(twoPi*real(imom*(ix(ieo)-src(1)),KR2) &
                            /real(nx,KR2))
            mysin(1,imom) = sin(twoPi*real(imom*(ix(ieo)-src(1)),KR2) &
                            /real(nx,KR2))
           enddo ! imom
           do id = 1,4
            do jd = 1,4
             bit3 = 0.0_KR2
             do ic = 1,3
              do jc = 1,3
               bit3 = bit3 &
                    + real(q(jc,jd,1,ic,id,ieo)**2+q(jc,jd,2,ic,id,ieo)**2,KR2)
              enddo ! jc
             enddo ! ic
             bit4 = dsgn(id)*dsgn(jd)*bit3
             do izmom = 0,4
              do iymom = 0,4
               do ixmom = 0,4
                pbit(it,ixmom,iymom,izmom) = pbit(it,ixmom,iymom,izmom) &
                            + bit3*mycos(1,ixmom)*mycos(2,iymom)*mycos(3,izmom)
                rbit(it,ixmom,iymom,izmom) = rbit(it,ixmom,iymom,izmom) &
                            + bit4*mycos(1,ixmom)*mycos(2,iymom)*mycos(3,izmom)
               enddo ! ixmom
              enddo ! iymom
             enddo ! izmom
             if (lessnoise==0) then
              do izmom = 1,4
               do iymom = 1,4
                pbit(it,0,iymom,izmom) = pbit(it,0,iymom,izmom) &
                                       - bit3*mysin(2,iymom)*mysin(3,izmom)
                pbit(it,iymom,0,izmom) = pbit(it,iymom,0,izmom) &
                                       - bit3*mysin(1,iymom)*mysin(3,izmom)
                pbit(it,iymom,izmom,0) = pbit(it,iymom,izmom,0) &
                                       - bit3*mysin(1,iymom)*mysin(2,izmom)
                rbit(it,0,iymom,izmom) = rbit(it,0,iymom,izmom) &
                                       - bit4*mysin(2,iymom)*mysin(3,izmom)
                rbit(it,iymom,0,izmom) = rbit(it,iymom,0,izmom) &
                                       - bit4*mysin(1,iymom)*mysin(3,izmom)
                rbit(it,iymom,izmom,0) = rbit(it,iymom,izmom,0) &
                                       - bit4*mysin(1,iymom)*mysin(2,izmom)
                do ixmom = 1,4
                 pbit(it,ixmom,iymom,izmom) = pbit(it,ixmom,iymom,izmom) &
                           - bit3*(mysin(1,ixmom)*mysin(2,iymom)*mycos(3,izmom)&
                                  +mysin(1,ixmom)*mycos(2,iymom)*mysin(3,izmom)&
                                  +mycos(1,ixmom)*mysin(2,iymom)*mysin(3,izmom))
                 rbit(it,ixmom,iymom,izmom) = rbit(it,ixmom,iymom,izmom) &
                           - bit4*(mysin(1,ixmom)*mysin(2,iymom)*mycos(3,izmom)&
                                  +mysin(1,ixmom)*mycos(2,iymom)*mysin(3,izmom)&
                                  +mycos(1,ixmom)*mysin(2,iymom)*mysin(3,izmom))
                enddo ! ixmom
               enddo ! iymom
              enddo ! izmom
             endif
            enddo ! jd
           enddo ! id
          enddo ! ieo
   

! Construct the nucleon correlator.

! -the colours are ca,cb,cc,cd,ce,cf.
! -the internal Dirac indices are da,db,dc,dd.
! -the external Dirac indices are source=jd, sink=id.

! ASIDE: Typically, one only uses jd=id.  However, the off-diagonal entries
!        can be important when this nucleon two-point function is one piece
!        within a larger disconnected diagram, such as the magnetic form factor.

! NOTE ~DD: the dc and dd internal indices are doing the gamma5 multiplications
!           due to the verticies in the nucleon construction. (7-31-05)

          do ieo = 1,2
           mycos(1,1) = cos(twoPi*real(ix(ieo)-src(1),KR2)/real(nx,KR2))
           mysin(1,1) = sin(twoPi*real(ix(ieo)-src(1),KR2)/real(nx,KR2))
           mycos(1,2) = cos(twoPi*real(2*(ix(ieo)-src(1)),KR2)/real(nx,KR2))
           mysin(1,2) = sin(twoPi*real(2*(ix(ieo)-src(1)),KR2)/real(nx,KR2))
           bit1Re = 0.0_KR2
           bit1Im = 0.0_KR2
           bit2Re = 0.0_KR2
           bit2Im = 0.0_KR2
           do icol = 1,6
            ca = col(icol,1)
            cb = col(icol,2)
            cc = col(icol,3)
            do jcol = 1,6
             cd = col(jcol,1)
             ce = col(jcol,2)
             cf = col(jcol,3)
             do jd = 1,4
              do id = 1,4
               do da = 1,4
                dc = cgam5(da)
                do db = 1,4
                 dd = cgam5(db)
                 bit1Re(jd,id) = bit1Re(jd,id)                                 &
                 + csgn(icol)*csgn(jcol)*dsgn(da)*dsgn(db)*real(               &
                 q(cd,jd,1,ca,id,ieo)*q(ce,da,1,cb,db,ieo)*q(cf,dc,1,cc,dd,ieo)&
                -q(cd,jd,1,ca,id,ieo)*q(ce,da,2,cb,db,ieo)*q(cf,dc,2,cc,dd,ieo)&
                -q(cd,jd,2,ca,id,ieo)*q(ce,da,1,cb,db,ieo)*q(cf,dc,2,cc,dd,ieo)&
                -q(cd,jd,2,ca,id,ieo)*q(ce,da,2,cb,db,ieo)*q(cf,dc,1,cc,dd,ieo)&
                 ,KR2)
                 bit1Im(jd,id) = bit1Im(jd,id)                                 &
                 + csgn(icol)*csgn(jcol)*dsgn(da)*dsgn(db)*real(               &
                 q(cd,jd,2,ca,id,ieo)*q(ce,da,1,cb,db,ieo)*q(cf,dc,1,cc,dd,ieo)&
                +q(cd,jd,1,ca,id,ieo)*q(ce,da,2,cb,db,ieo)*q(cf,dc,1,cc,dd,ieo)&
                +q(cd,jd,1,ca,id,ieo)*q(ce,da,1,cb,db,ieo)*q(cf,dc,2,cc,dd,ieo)&
                -q(cd,jd,2,ca,id,ieo)*q(ce,da,2,cb,db,ieo)*q(cf,dc,2,cc,dd,ieo)&
                 ,KR2)
                 bit2Re(jd,id) = bit2Re(jd,id)                                 &
                 + csgn(icol)*csgn(jcol)*dsgn(da)*dsgn(db)*real(               &
                 q(cd,dc,1,ca,id,ieo)*q(ce,da,1,cb,db,ieo)*q(cf,jd,1,cc,dd,ieo)&
                -q(cd,dc,1,ca,id,ieo)*q(ce,da,2,cb,db,ieo)*q(cf,jd,2,cc,dd,ieo)&
                -q(cd,dc,2,ca,id,ieo)*q(ce,da,1,cb,db,ieo)*q(cf,jd,2,cc,dd,ieo)&
                -q(cd,dc,2,ca,id,ieo)*q(ce,da,2,cb,db,ieo)*q(cf,jd,1,cc,dd,ieo)&
                 ,KR2)
                 bit2Im(jd,id) = bit2Im(jd,id)                                 &
                 + csgn(icol)*csgn(jcol)*dsgn(da)*dsgn(db)*real(               &
                 q(cd,dc,2,ca,id,ieo)*q(ce,da,1,cb,db,ieo)*q(cf,jd,1,cc,dd,ieo)&
                +q(cd,dc,1,ca,id,ieo)*q(ce,da,2,cb,db,ieo)*q(cf,jd,1,cc,dd,ieo)&
                +q(cd,dc,1,ca,id,ieo)*q(ce,da,1,cb,db,ieo)*q(cf,jd,2,cc,dd,ieo)&
                -q(cd,dc,2,ca,id,ieo)*q(ce,da,2,cb,db,ieo)*q(cf,jd,2,cc,dd,ieo)&
                 ,KR2)
                enddo ! db
               enddo ! da
              enddo ! id
             enddo ! jd
            enddo ! jcol
           enddo ! icol
           bitRe(:,:) = bit1Re(:,:) + bit2Re(:,:)
           bitIm(:,:) = bit1Im(:,:) + bit2Im(:,:)
           nbit(1,it,:,:,0,0,0) = nbit(1,it,:,:,0,0,0) + bitRe(:,:)
           nbit(2,it,:,:,0,0,0) = nbit(2,it,:,:,0,0,0) + bitIm(:,:)
           nbit(1,it,:,:,1,0,0) = nbit(1,it,:,:,1,0,0) + bitRe(:,:)*mycos(1,1)
           nbit(2,it,:,:,1,0,0) = nbit(2,it,:,:,1,0,0) + bitIm(:,:)*mycos(1,1)
           nbit(1,it,:,:,0,1,0) = nbit(1,it,:,:,0,1,0) + bitRe(:,:)*mycos(2,1)
           nbit(2,it,:,:,0,1,0) = nbit(2,it,:,:,0,1,0) + bitIm(:,:)*mycos(2,1)
           nbit(1,it,:,:,0,0,1) = nbit(1,it,:,:,0,0,1) + bitRe(:,:)*mycos(3,1)
           nbit(2,it,:,:,0,0,1) = nbit(2,it,:,:,0,0,1) + bitIm(:,:)*mycos(3,1)
           nbit(1,it,:,:,1,1,0) = nbit(1,it,:,:,1,1,0) &
                                + bitRe(:,:)*mycos(1,1)*mycos(2,1)
           nbit(2,it,:,:,1,1,0) = nbit(2,it,:,:,1,1,0) &
                                + bitIm(:,:)*mycos(1,1)*mycos(2,1)
           nbit(1,it,:,:,1,0,1) = nbit(1,it,:,:,1,0,1) &
                                + bitRe(:,:)*mycos(1,1)*mycos(3,1)
           nbit(2,it,:,:,1,0,1) = nbit(2,it,:,:,1,0,1) &
                                + bitIm(:,:)*mycos(1,1)*mycos(3,1)
           nbit(1,it,:,:,0,1,1) = nbit(1,it,:,:,0,1,1) &
                                + bitRe(:,:)*mycos(2,1)*mycos(3,1)
           nbit(2,it,:,:,0,1,1) = nbit(2,it,:,:,0,1,1) &
                                + bitIm(:,:)*mycos(2,1)*mycos(3,1)
           nbit(1,it,:,:,1,1,1) = nbit(1,it,:,:,1,1,1) &
                                + bitRe(:,:)*mycos(1,1)*mycos(2,1)*mycos(3,1)
           nbit(2,it,:,:,1,1,1) = nbit(2,it,:,:,1,1,1) &
                                + bitIm(:,:)*mycos(1,1)*mycos(2,1)*mycos(3,1)
           nbit(1,it,:,:,2,0,0) = nbit(1,it,:,:,2,0,0) + bitRe(:,:)*mycos(1,2)
           nbit(2,it,:,:,2,0,0) = nbit(2,it,:,:,2,0,0) + bitIm(:,:)*mycos(1,2)
           nbit(1,it,:,:,0,2,0) = nbit(1,it,:,:,0,2,0) + bitRe(:,:)*mycos(2,2)
           nbit(2,it,:,:,0,2,0) = nbit(2,it,:,:,0,2,0) + bitIm(:,:)*mycos(2,2)
           nbit(1,it,:,:,0,0,2) = nbit(1,it,:,:,0,0,2) + bitRe(:,:)*mycos(3,2)
           nbit(2,it,:,:,0,0,2) = nbit(2,it,:,:,0,0,2) + bitIm(:,:)*mycos(3,2)
           if (lessnoise==0) then
            nbit(1,it,:,:,1,0,0) = nbit(1,it,:,:,1,0,0) - bitIm(:,:)*mysin(1,1)
            nbit(2,it,:,:,1,0,0) = nbit(2,it,:,:,1,0,0) + bitRe(:,:)*mysin(1,1)
            nbit(1,it,:,:,0,1,0) = nbit(1,it,:,:,0,1,0) - bitIm(:,:)*mysin(2,1)
            nbit(2,it,:,:,0,1,0) = nbit(2,it,:,:,0,1,0) + bitRe(:,:)*mysin(2,1)
            nbit(1,it,:,:,0,0,1) = nbit(1,it,:,:,0,0,1) - bitIm(:,:)*mysin(3,1)
            nbit(2,it,:,:,0,0,1) = nbit(2,it,:,:,0,0,1) + bitRe(:,:)*mysin(3,1)
            nbit(1,it,:,:,1,1,0) = nbit(1,it,:,:,1,1,0) &
                                 - bitRe(:,:)*mysin(1,1)*mysin(2,1) &
                                 - bitIm(:,:)*(mycos(1,1)*mysin(2,1) &
                                              +mysin(1,1)*mycos(2,1))
            nbit(2,it,:,:,1,1,0) = nbit(2,it,:,:,1,1,0) &
                                 - bitIm(:,:)*mysin(1,1)*mysin(2,1) &
                                 + bitRe(:,:)*(mycos(1,1)*mysin(2,1) &
                                              +mysin(1,1)*mycos(2,1))
            nbit(1,it,:,:,1,0,1) = nbit(1,it,:,:,1,0,1) &
                                 - bitRe(:,:)*mysin(1,1)*mysin(3,1) &
                                 - bitIm(:,:)*(mycos(1,1)*mysin(3,1) &
                                              +mysin(1,1)*mycos(3,1))
            nbit(2,it,:,:,1,0,1) = nbit(2,it,:,:,1,0,1) &
                                 - bitIm(:,:)*mysin(1,1)*mysin(3,1) &
                                 + bitRe(:,:)*(mycos(1,1)*mysin(3,1) &
                                              +mysin(1,1)*mycos(3,1))
            nbit(1,it,:,:,0,1,1) = nbit(1,it,:,:,0,1,1) &
                                 - bitRe(:,:)*mysin(2,1)*mysin(3,1) &
                                 - bitIm(:,:)*(mycos(2,1)*mysin(3,1) &
                                              +mysin(2,1)*mycos(3,1))
            nbit(2,it,:,:,0,1,1) = nbit(2,it,:,:,0,1,1) &
                                 - bitIm(:,:)*mysin(2,1)*mysin(3,1) &
                                 + bitRe(:,:)*(mycos(2,1)*mysin(3,1) &
                                              +mysin(2,1)*mycos(3,1))
            nbit(1,it,:,:,1,1,1) = nbit(1,it,:,:,1,1,1)                        &
                                - bitRe(:,:)*(mysin(1,1)*mysin(2,1)*mycos(3,1) &
                                             +mysin(1,1)*mycos(2,1)*mysin(3,1) &
                                             +mycos(1,1)*mysin(2,1)*mysin(3,1))&
                                - bitIm(:,:)*(mycos(1,1)*mycos(2,1)*mysin(3,1) &
                                             +mycos(1,1)*mysin(2,1)*mycos(3,1) &
                                             +mysin(1,1)*mycos(2,1)*mycos(3,1) &
                                             -mysin(1,1)*mysin(2,1)*mysin(3,1))
            nbit(2,it,:,:,1,1,1) = nbit(2,it,:,:,1,1,1)                        &
                                - bitIm(:,:)*(mysin(1,1)*mysin(2,1)*mycos(3,1) &
                                             +mysin(1,1)*mycos(2,1)*mysin(3,1) &
                                             +mycos(1,1)*mysin(2,1)*mysin(3,1))&
                                + bitRe(:,:)*(mycos(1,1)*mycos(2,1)*mysin(3,1) &
                                             +mycos(1,1)*mysin(2,1)*mycos(3,1) &
                                             +mysin(1,1)*mycos(2,1)*mycos(3,1) &
                                             -mysin(1,1)*mysin(2,1)*mysin(3,1))
            nbit(1,it,:,:,2,0,0) = nbit(1,it,:,:,2,0,0) - bitIm(:,:)*mysin(1,2)
            nbit(2,it,:,:,2,0,0) = nbit(2,it,:,:,2,0,0) + bitRe(:,:)*mysin(1,2)
            nbit(1,it,:,:,0,2,0) = nbit(1,it,:,:,0,2,0) - bitIm(:,:)*mysin(2,2)
            nbit(2,it,:,:,0,2,0) = nbit(2,it,:,:,0,2,0) + bitRe(:,:)*mysin(2,2)
            nbit(1,it,:,:,0,0,2) = nbit(1,it,:,:,0,0,2) - bitIm(:,:)*mysin(3,2)
            nbit(2,it,:,:,0,0,2) = nbit(2,it,:,:,0,0,2) + bitRe(:,:)*mysin(3,2)
           endif
          enddo ! ieo
 
! End main loop.
         enddo ! ibl
        enddo ! ixbit
       enddo ! iybit
      enddo ! izbit
     enddo ! itbit
 
! Close the necessary files.
     do jc = 1,3
      do jd = 1,4
       ifile = 100 + jd + (jc-1)*4
       close(unit=ifile,status="keep")
      enddo ! jd
     enddo ! jc
 
! Combine the results from all processes.
     pcorr = real(pbit,KR)
     rcorr = real(rbit,KR)
     ncorr = real(nbit,KR)
     if (nps==1) then
      pion = pcorr
      rho = rcorr
      nucl = ncorr
     else
      nmesg = 125*nt
      call MPI_REDUCE(pcorr(1,0,0,0),pion(1,0,0,0),nmesg,MRT,MPI_SUM,0, &
                      MPI_COMM_WORLD,ierr)
      call MPI_BCAST(pion(1,0,0,0),nmesg,MRT,0,MPI_COMM_WORLD,ierr)
      call MPI_REDUCE(rcorr(1,0,0,0),rho(1,0,0,0),nmesg,MRT,MPI_SUM,0, &
                      MPI_COMM_WORLD,ierr)
      call MPI_BCAST(rho(1,0,0,0),nmesg,MRT,0,MPI_COMM_WORLD,ierr)
      nmesg = 32*27*nt
      call MPI_REDUCE(ncorr(1,1,1,1,0,0,0),nucl(1,1,1,1,0,0,0),nmesg,MRT, &
                      MPI_SUM,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(nucl(1,1,1,1,0,0,0),nmesg,MRT,0,MPI_COMM_WORLD,ierr)
     endif
 
! Record the results.
     if (myid==0) then
      open(unit=8,file=trim(rwdir(myid+1))//"CFGSPROPS.LOG",action="write", &
           form="formatted",status="old",position="append")
       do iz = 0,4
        do iy = 0,4
         do jx = 0,4
         !write(unit=8,fmt="(i1,a14,i1,a1,i1,a1,i1,a5,500es17.10)") &
         !    ikappa,"pion      [p=(",jx,",",iy,",",iz,")] = ", pion(:,jx,iy,iz)
         !write(unit=8,fmt="(i1,a14,i1,a1,i1,a1,i1,a5,500es17.10)") &
         !    ikappa,"rho       [p=(",jx,",",iy,",",iz,")] = ", rho(:,jx,iy,iz)
          if (jx**2+iy**2+iz**2<5) then
           do jd = 1,4
            do id = 1,4
             write(unit=8,fmt="(i1,a7,2i1,a5,i1,a1,i1,a1,i1,a5,500es17.10)") &
                   ikappa,"nuclRe,",jd,id," [p=(",jx,",",iy,",",iz,")] = ",  &
                   nucl(1,:,jd,id,jx,iy,iz)
             write(unit=8,fmt="(i1,a7,2i1,a5,i1,a1,i1,a1,i1,a5,500es17.10)") &
                   ikappa,"nuclIm,",jd,id," [p=(",jx,",",iy,",",iz,")] = ",  &
                   nucl(2,:,jd,id,jx,iy,iz)
            enddo ! id
           enddo ! jd
          endif
         enddo ! jx
        enddo ! iy
       enddo ! iz
      close(unit=8,status="keep")
     endif

     if (myid==0) then
      open(unit=8,file=trim(rwdir(myid+1))//"DUBLIN.OUT",action="write", &
           form="formatted",status="old",position="append")
        znuc = 0.0_KR
        do jd=1,2
            do ittt = 1,nt
              znuc(ittt,1) = znuc(ittt,1) + nucl(1,ittt,jd,jd,0,0,0)
              znuc(ittt,2) = znuc(ittt,2) + nucl(2,ittt,jd,jd,0,0,0)
            enddo ! ittt
        enddo ! jd
        write(unit=8,fmt="(a12,1es17.10)") "quark mass=", ikappa
        do ittt=1,nt
          write(unit=8,fmt="(i5,2es17.10)")  ittt,znuc(ittt,1),znuc(ittt,2)
        enddo ! ittt
        close(unit=8,status="keep")
     endif

    enddo ! iLS
    enddo ! ikappa

 end subroutine corr

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 subroutine nucleoncorr(nkappa,ntmqcd,nsmear,src,bc,myid,rwdir,cfgfile,&
                        particle,MRT,icfgsave)
!
! This subroutine will eventually get replaced by subroutine baryoncorr,
! but since that has not yet been written this old code is retained.
!
! WARNING! WARNING! WARNING! WARNING! WARNING! WARNING! WARNING! WARNING!
! This subroutine is only correct for cSW=0.0_KR
! WARNING! WARNING! WARNING! WARNING! WARNING! WARNING! WARNING! WARNING!
!
! Compute the pion, rho and nucleon correlators (Re and Im parts and all
! Dirac components) using local operators.
!
! This subroutine uses N = epsilon^(abc) q^a [ (q^b)^T C*gamma5 q^c ]
!
!                  / 0 -1  0 0 \
! where C*gamma5 = | 1  0  0 0 |
!                  | 0  0  0 1 |
!                  \ 0  0 -1 0 /
!
! The Dirac index for the creation operator is set equal to the Dirac index
! for the annihilation operator.  (All four components are computed.)

    integer(kind=KI), intent(in)               :: nkappa, ntmqcd, nsmear, &
                                                  particle, myid, MRT
    integer(kind=KI), intent(in), dimension(:) :: src, bc
    integer(Kind=KI), intent(in)               :: icfgsave
    character(len=*), intent(in), dimension(:) :: rwdir
    character(len=*), intent(in)               :: cfgfile

    character(len=128)                       :: propfile, opfile
    character(len=12)                        :: particlestring
    character(len=8)                         :: trailer
    character(len=16)                         :: ntrailer
    integer(kind=KI), dimension(4)           :: np, ip, cgam5, dsgn
    integer(kind=KI), dimension(6,3)         :: col
    real(kind=KR2),   parameter              :: twoPi=6.283185307179586_KR2
    real(kind=KR2),   dimension(6)           :: csgn
    real(kind=KR2)                           :: bit3, bit4
    real(kind=KR2),   dimension(4,4)         :: bitRe, bitIm, bit1Re, bit1Im, &
                                                bit2Re, bit2Im
    real(kind=KR),    dimension(3,4,2,3,4,2) :: q, q1, q2 
    real(kind=KR),    dimension(3,4,2,3,4,2) :: qpion1, qpion2
    integer(kind=KI), dimension(2)           :: ix
    real(kind=KR2),   dimension(3,0:4)       :: mycos, mysin
    real(kind=KR2),   dimension(nt,0:4,0:4,0:4)       :: pbit, rbit
    real(kind=KR),    dimension(nt,0:4,0:4,0:4)       :: pcorr, rcorr, pion, rho
    real(kind=KR2),   dimension(2,nt,4,4,0:2,0:2,0:2) :: nbit
    real(kind=KR),    dimension(2,nt,4,4,0:2,0:2,0:2) :: ncorr, nucl
    character(len=1), dimension(2)                    :: LSsymbol
    integer(kind=KI) :: jc, jd, ifile, itbit, jx, iy, iz, it, ibl, ic, id, &
                        iri, ieo, ieo1, ieo2, icol, jcol, ca, cb, cc, cd, ce, &
                        cf, da, db, dc, dd, nmesg, ierr, itbit2, iblbit, &
                        ixbit, iybit, izbit, ixbit2, iybit2, izbit2, ixbit3, &
                        lessnoise, ikappa, nLS, iLS, ncalcs, imom, ixmom, &
                        iymom, izmom

    real(kind=KR),   dimension(2,nt)                  :: znuc
    real(kind=KR)                                     :: fac
    integer(kind=KI)  :: ittt, row, row2 
    logical           :: fileexists

! Compute the nucleon correlator from a local operator.
! NOTE: it=1 is uniquely defined by          1<=isite<=  factor and 1<=ibl<= 8
!       it=2 is uniquely defined by          1<=isite<=  factor and 9<=ibl<=16
!       it=3 is uniquely defined by   factor+1<=isite<=2*factor and 1<=ibl<= 8
!       it=4 is uniquely defined by   factor+1<=isite<=2*factor and 9<=ibl<=16
!       it=5 is uniquely defined by 2*factor+1<=isite<=3*factor and 1<=ibl<= 8
!       it=6 is uniquely defined by 2*factor+1<=isite<=3*factor and 9<=ibl<=16
!       etc, where factor = 2*nvhalf*npt/nt.

! Some parameters.
     if (nsmear==0) then
      nLS = 1
      LSsymbol(1) = "L"
     elseif (nsmear>0) then
      nLS = 1
      LSsymbol(1) = "S"
     else
      nLS = 2
      LSsymbol(1) = "L"
      LSsymbol(2) = "S"
     endif
!  need to change this to do the nucleon tmU and tmD propagators for
!  each mass independantly.

! ncalcs is the number of masses that are being looped over.

      ncalcs = abs(ntmqcd)

! Consider multiple quark masses SEPARATELY.  (Hopefully I can mix them soon.)
    do ikappa = 1,ncalcs
    do iLS = 1,nLS
!      if (ikappa<=nkappa) then
!       trailer = ".xxk"//LSsymbol(iLS)
!       write(unit=trailer(4:4),fmt="(i1.1)") ikappa
!      elseif (ikappa<=nkappa+abs(ntmqcd)) then
!       trailer = ".xxk"//LSsymbol(iLS)//"tmU"
!       write(unit=trailer(4:4),fmt="(i1.1)") ikappa - nkappa
!      else
!       trailer = ".xxk"//LSsymbol(iLS)//"tmD"
!       write(unit=trailer(4:4),fmt="(i1.1)") ikappa - nkappa - abs(ntmqcd)
!      endif

! Open the necessary tmU quark propagator files.
     do jc = 1,3
      do jd = 1,4
       trailer = ".xxk"//LSsymbol(iLS)//"tmU"
       write(unit=trailer(2:2),fmt="(i1.1)") jc
       write(unit=trailer(3:3),fmt="(i1.1)") jd
       write(unit=trailer(4:4),fmt="(i1.1)") ikappa
       propfile = trim(rwdir(myid+1))//trim(cfgfile)//trim(trailer)
       ifile = 100 + jd + (jc-1)*4
       open(unit=ifile,file=propfile,action="read",status="old", &
            position="rewind",form="unformatted")
      enddo ! jd
     enddo ! jc

! Open the necessary tmD quark propagator files.
     do jc = 1,3
      do jd = 1,4
       trailer = ".xxk"//LSsymbol(iLS)//"tmD"
       write(unit=trailer(2:2),fmt="(i1.1)") jc
       write(unit=trailer(3:3),fmt="(i1.1)") jd
       write(unit=trailer(4:4),fmt="(i1.1)") ikappa
       propfile = trim(rwdir(myid+1))//trim(cfgfile)//trim(trailer)
       ifile = 200 + jd + (jc-1)*4
       open(unit=ifile,file=propfile,action="read",status="old", &
            position="rewind",form="unformatted")
      enddo ! jd
     enddo ! jc

! Some initializations.
     pbit = 0.0_KR2
     rbit = 0.0_KR2
     nbit = 0.0_KR2
     col(1,1:3) = (/ 1, 2, 3 /)
     col(2,1:3) = (/ 2, 3, 1 /)
     col(3,1:3) = (/ 3, 1, 2 /)
     col(4,1:3) = (/ 1, 3, 2 /)
     col(5,1:3) = (/ 2, 1, 3 /)
     col(6,1:3) = (/ 3, 2, 1 /)


! NOTE ~ cgam5 = -1*(gam5)^-1 and (cgam5)^-1 = (cgam5)^T
     cgam5(1:4) = (/ 2, 1, 4, 3 /)

     csgn(1:6) = (/ 1.0_KR2, 1.0_KR2, 1.0_KR2, -1.0_KR2, -1.0_KR2, -1.0_KR2 /)
     dsgn(1:4) = (/ 1.0_KR2, -1.0_KR2, -1.0_KR2, 1.0_KR2 /)
     np(1) = npx
     np(2) = npy
     np(3) = npz
     np(4) = npt
     call atoc(myid,np,ip)

! NOTICE: If the boundary conditions are periodic in ALL spatial directions,
!         then the terms which are odd in sinx, siny, sinz, sin2x, sin2y or
!         sin2z are omitted, since the integral of an even periodic function,
!         times sinx or sin2x etc., vanishes.  This should reduce the
!         statistical noise.
     if (bc(1)==1.and.bc(2)==1.and.bc(3)==1) then
      lessnoise = 1
     else
      lessnoise = 0
     endif

! Begin main loop, constructing ix,iy,iz,it from isite,ieo,ibl.
     ieo1 = 2
     ieo2 = 1
     do itbit = 2,nt/npt,2
      itbit2 = itbit + ip(4)*nt/npt
      ieo1 = 3 - ieo1
      ieo2 = 3 - ieo2
      do izbit = 2,nz/npz,2
       izbit2 = izbit + ip(3)*nz/npz
       ieo1 = 3 - ieo1
       ieo2 = 3 - ieo2
       do iybit = 2,ny/npy,2
        iybit2 = iybit + ip(2)*ny/npy
        ieo1 = 3 - ieo1
        ieo2 = 3 - ieo2
        do ixbit = 4,nx/npx,4
         ixbit2 = ixbit + ip(1)*nx/npx
         do ibl = 1,16
          if (ibl>8) then
           it = itbit2
          else
           it = itbit2 - 1
          endif
          iblbit = 1 + modulo(ibl-1,8)
          if (iblbit>4) then
           iz = izbit2
          else
           iz = izbit2 - 1
          endif
          iblbit = 1 + modulo(iblbit-1,4)
          if (iblbit>2) then
           iy = iybit2
          else
           iy = iybit2 - 1
          endif
          if (modulo(ibl,2)==1) then
           ixbit3 = ixbit2 - 1
          else
           ixbit3 = ixbit2
          endif
          ix(ieo1) = ixbit3 - 2
          ix(ieo2) = ixbit3

! Factors for the smallest momenta in the y and z directions.
          do imom = 0,4
           mycos(2,imom) = cos(twoPi*real(imom*(iy-src(2)),KR2)/real(ny,KR2))
           mysin(2,imom) = sin(twoPi*real(imom*(iy-src(2)),KR2)/real(ny,KR2))
           mycos(3,imom) = cos(twoPi*real(imom*(iz-src(3)),KR2)/real(nz,KR2))
           mysin(3,imom) = sin(twoPi*real(imom*(iz-src(3)),KR2)/real(nz,KR2))
          enddo ! imom

! NOTE!!!! If you wish to calculate the proton (uud) then particle=1 and the 
!          subroutine will use tmU data twice for the nucleon construction  and only one tmD
!          otherwise for the neutron (udd) then the opposite is true.

! Read in the 12 tmU propagators for this pair of lattice sites (ieo=1 and ieo=2).
          do ic = 1,3
           do iri = 1,2
            do jc = 1,3
             do jd = 1,4
              ifile = 100 + jd + (jc-1)*4
               if (particle==1) then
                 read(unit=ifile)(q2(jc,jd,iri,ic,id,1),q2(jc,jd,iri,ic,id,2),id=1,4)
               else
                 read(unit=ifile)(q1(jc,jd,iri,ic,id,1),q1(jc,jd,iri,ic,id,2),id=1,4)
               endif ! particle
             enddo ! jd
            enddo ! jc
           enddo ! iri
          enddo ! ic

! Read in the 12 tmD propagators for this pair of lattice sites (ieo=1 and ieo=2).
          do ic = 1,3
           do iri = 1,2
            do jc = 1,3
             do jd = 1,4
              ifile = 200 + jd + (jc-1)*4
              if (particle == 1) then
                 read(unit=ifile)(q1(jc,jd,iri,ic,id,1),q1(jc,jd,iri,ic,id,2),id=1,4)
              else
                 read(unit=ifile)(q2(jc,jd,iri,ic,id,1),q2(jc,jd,iri,ic,id,2),id=1,4)
              endif ! particle
             enddo ! jd
            enddo ! jc
           enddo ! iri
          enddo ! ic

! Construct the pion and rho correlators.  (rho uses gamma_3 only.)

! This does the charged pion with tmU and tmD. Does not matter which one is 
! assigned to q1 and q2. The trace is real so we only need the real part of the
! multiplication.

! To make the charged pion we must multiply the "charged scalar" Tr[D_u D_d]
! by gamma_5 to give Tr[gamma_5 D_u gamma_5 (D_u)^dagger]


          if (ntmqcd/=0) then
            qpion1(:,1,:,:,:,:) = -q1(:,3,:,:,:,:) 
            qpion1(:,2,:,:,:,:) = -q1(:,4,:,:,:,:) 
            qpion1(:,3,:,:,:,:) =  q1(:,1,:,:,:,:) 
            qpion1(:,4,:,:,:,:) =  q1(:,2,:,:,:,:) 

            qpion2(:,1,:,:,:,:) = -q2(:,3,:,:,:,:) 
            qpion2(:,2,:,:,:,:) = -q2(:,4,:,:,:,:) 
            qpion2(:,3,:,:,:,:) =  q2(:,1,:,:,:,:) 
            qpion2(:,4,:,:,:,:) =  q2(:,2,:,:,:,:) 
          else
            qpion1 = q1
            qpion2 = q2
         endif ! ntmqcd

          do ieo = 1,2
           do imom = 0,4
            mycos(1,imom) = cos(twoPi*real(imom*(ix(ieo)-src(1)),KR2) &
                            /real(nx,KR2))
            mysin(1,imom) = sin(twoPi*real(imom*(ix(ieo)-src(1)),KR2) &
                            /real(nx,KR2))
           enddo ! imom
! NOTE ~ Since we want to calculate the charged pion with the 
!        gamma_5 rotations, we need the transpose of one of the 
!        propagators....(gamma_5*(D_u)^T)
           do id = 1,4
            do jd = 1,4
             bit3 = 0.0_KR2
             do ic = 1,3
              do jc = 1,3
               bit3 = bit3 &
                    + real(qpion1(jc,jd,1,ic,id,ieo)*qpion2(ic,id,1,jc,jd,ieo)+ &
                           qpion1(jc,jd,2,ic,id,ieo)*qpion2(ic,id,2,jc,jd,ieo),KR2)
              enddo ! jc
             enddo ! ic
             bit4 = dsgn(id)*dsgn(jd)*bit3
             do izmom = 0,4
              do iymom = 0,4
               do ixmom = 0,4
                pbit(it,ixmom,iymom,izmom) = pbit(it,ixmom,iymom,izmom) &
                            + bit3*mycos(1,ixmom)*mycos(2,iymom)*mycos(3,izmom)
                rbit(it,ixmom,iymom,izmom) = rbit(it,ixmom,iymom,izmom) &
                            + bit4*mycos(1,ixmom)*mycos(2,iymom)*mycos(3,izmom)
               enddo ! ixmom
              enddo ! iymom
             enddo ! izmom
             if (lessnoise==0) then
              do izmom = 1,4
               do iymom = 1,4
                pbit(it,0,iymom,izmom) = pbit(it,0,iymom,izmom) &
                                       - bit3*mysin(2,iymom)*mysin(3,izmom)
                pbit(it,iymom,0,izmom) = pbit(it,iymom,0,izmom) &
                                       - bit3*mysin(1,iymom)*mysin(3,izmom)
                pbit(it,iymom,izmom,0) = pbit(it,iymom,izmom,0) &
                                       - bit3*mysin(1,iymom)*mysin(2,izmom)
                rbit(it,0,iymom,izmom) = rbit(it,0,iymom,izmom) &
                                       - bit4*mysin(2,iymom)*mysin(3,izmom)
                rbit(it,iymom,0,izmom) = rbit(it,iymom,0,izmom) &
                                       - bit4*mysin(1,iymom)*mysin(3,izmom)
                rbit(it,iymom,izmom,0) = rbit(it,iymom,izmom,0) &
                                       - bit4*mysin(1,iymom)*mysin(2,izmom)
                do ixmom = 1,4
                 pbit(it,ixmom,iymom,izmom) = pbit(it,ixmom,iymom,izmom) &
                           - bit3*(mysin(1,ixmom)*mysin(2,iymom)*mycos(3,izmom)&
                                  +mysin(1,ixmom)*mycos(2,iymom)*mysin(3,izmom)&
                                  +mycos(1,ixmom)*mysin(2,iymom)*mysin(3,izmom))
                 rbit(it,ixmom,iymom,izmom) = rbit(it,ixmom,iymom,izmom) &
                           - bit4*(mysin(1,ixmom)*mysin(2,iymom)*mycos(3,izmom)&
                                  +mysin(1,ixmom)*mycos(2,iymom)*mysin(3,izmom)&
                                  +mycos(1,ixmom)*mysin(2,iymom)*mysin(3,izmom))
                enddo ! ixmom
               enddo ! iymom
              enddo ! izmom
             endif
            enddo ! jd
           enddo ! id
          enddo ! ieo


! Construct the nucleon correlator.

! -the colours are ca,cb,cc,cd,ce,cf.
! -the internal Dirac indices are da,db,dc,dd.
! -the external Dirac indices are source=jd, sink=id.

! ASIDE: Typically, one only uses jd=id.  However, the off-diagonal entries
!        can be important when this nucleon two-point function is one piece
!        within a larger disconnected diagram, such as the magnetic form factor.

! NOTE ~DD: the dc and dd internal indices are doing the gamma5 multiplications
!           due to the verticies in the nucleon construction. (7-31-05)

          do ieo = 1,2
           mycos(1,1) = cos(twoPi*real(ix(ieo)-src(1),KR2)/real(nx,KR2))
           mysin(1,1) = sin(twoPi*real(ix(ieo)-src(1),KR2)/real(nx,KR2))
           mycos(1,2) = cos(twoPi*real(2*(ix(ieo)-src(1)),KR2)/real(nx,KR2))
           mysin(1,2) = sin(twoPi*real(2*(ix(ieo)-src(1)),KR2)/real(nx,KR2))
           bit1Re = 0.0_KR2
           bit1Im = 0.0_KR2
           bit2Re = 0.0_KR2
           bit2Im = 0.0_KR2
           do icol = 1,6
            ca = col(icol,1)
            cb = col(icol,2)
            cc = col(icol,3)
            do jcol = 1,6
             cd = col(jcol,1)
             ce = col(jcol,2)
             cf = col(jcol,3)
             do jd = 1,4
              do id = 1,4
               do da = 1,4
                dc = cgam5(da)
                do db = 1,4
                 dd = cgam5(db)

! I need to make sure that the traces are on the correct objects.

                 bit1Re(jd,id) = bit1Re(jd,id)                                 &
                 + csgn(icol)*csgn(jcol)*dsgn(da)*dsgn(db)*real(               &
                 q2(cd,jd,1,ca,id,ieo)*q1(ce,da,1,cb,db,ieo)*q2(cf,dc,1,cc,dd,ieo)&
                -q2(cd,jd,1,ca,id,ieo)*q1(ce,da,2,cb,db,ieo)*q2(cf,dc,2,cc,dd,ieo)&
                -q2(cd,jd,2,ca,id,ieo)*q1(ce,da,1,cb,db,ieo)*q2(cf,dc,2,cc,dd,ieo)&
                -q2(cd,jd,2,ca,id,ieo)*q1(ce,da,2,cb,db,ieo)*q2(cf,dc,1,cc,dd,ieo)&
                 ,KR2)

                 bit1Im(jd,id) = bit1Im(jd,id)                                 &
                 + csgn(icol)*csgn(jcol)*dsgn(da)*dsgn(db)*real(               &
                 q2(cd,jd,2,ca,id,ieo)*q1(ce,da,1,cb,db,ieo)*q2(cf,dc,1,cc,dd,ieo)&
                +q2(cd,jd,1,ca,id,ieo)*q1(ce,da,2,cb,db,ieo)*q2(cf,dc,1,cc,dd,ieo)&
                +q2(cd,jd,1,ca,id,ieo)*q1(ce,da,1,cb,db,ieo)*q2(cf,dc,2,cc,dd,ieo)&
                -q2(cd,jd,2,ca,id,ieo)*q1(ce,da,2,cb,db,ieo)*q2(cf,dc,2,cc,dd,ieo)&
                 ,KR2)
 
! I think that bit2 is supposed to be q1*q2*q1...for the following reason,
! use tr[D_u Dbar_d D_u] = tr[D_u (cgam5 D_d cgam5^-1)^T D_u] 
!                        = tr[ (D_u*(cgam5^-1)^T (D_d)^T (cgam5)^T*D_u]


                 bit2Re(jd,id) = bit2Re(jd,id)                                 &
                 + csgn(icol)*csgn(jcol)*dsgn(da)*dsgn(db)*real(               &
                 q2(cd,dc,1,ca,id,ieo)*q1(ce,da,1,cb,db,ieo)*q2(cf,jd,1,cc,dd,ieo)&
                -q2(cd,dc,1,ca,id,ieo)*q1(ce,da,2,cb,db,ieo)*q2(cf,jd,2,cc,dd,ieo)&
                -q2(cd,dc,2,ca,id,ieo)*q1(ce,da,1,cb,db,ieo)*q2(cf,jd,2,cc,dd,ieo)&
                -q2(cd,dc,2,ca,id,ieo)*q1(ce,da,2,cb,db,ieo)*q2(cf,jd,1,cc,dd,ieo)&
                 ,KR2)

                 bit2Im(jd,id) = bit2Im(jd,id)                                 &
                 + csgn(icol)*csgn(jcol)*dsgn(da)*dsgn(db)*real(               &
                 q2(cd,dc,2,ca,id,ieo)*q1(ce,da,1,cb,db,ieo)*q2(cf,jd,1,cc,dd,ieo)&
                +q2(cd,dc,1,ca,id,ieo)*q1(ce,da,2,cb,db,ieo)*q2(cf,jd,1,cc,dd,ieo)&
                +q2(cd,dc,1,ca,id,ieo)*q1(ce,da,1,cb,db,ieo)*q2(cf,jd,2,cc,dd,ieo)&
                -q2(cd,dc,2,ca,id,ieo)*q1(ce,da,2,cb,db,ieo)*q2(cf,jd,2,cc,dd,ieo)&
                 ,KR2)

                enddo ! db
               enddo ! da
              enddo ! id
             enddo ! jd
            enddo ! jcol
           enddo ! icol


           bitRe(:,:) = bit1Re(:,:) + bit2Re(:,:)
           bitIm(:,:) = bit1Im(:,:) + bit2Im(:,:)
          
           nbit(1,it,:,:,0,0,0) = nbit(1,it,:,:,0,0,0) + bitRe(:,:)
           nbit(2,it,:,:,0,0,0) = nbit(2,it,:,:,0,0,0) + bitIm(:,:)
           nbit(1,it,:,:,1,0,0) = nbit(1,it,:,:,1,0,0) + bitRe(:,:)*mycos(1,1)
           nbit(2,it,:,:,1,0,0) = nbit(2,it,:,:,1,0,0) + bitIm(:,:)*mycos(1,1)
           nbit(1,it,:,:,0,1,0) = nbit(1,it,:,:,0,1,0) + bitRe(:,:)*mycos(2,1)
           nbit(2,it,:,:,0,1,0) = nbit(2,it,:,:,0,1,0) + bitIm(:,:)*mycos(2,1)
           nbit(1,it,:,:,0,0,1) = nbit(1,it,:,:,0,0,1) + bitRe(:,:)*mycos(3,1)
           nbit(2,it,:,:,0,0,1) = nbit(2,it,:,:,0,0,1) + bitIm(:,:)*mycos(3,1)
           nbit(1,it,:,:,1,1,0) = nbit(1,it,:,:,1,1,0) &
                                + bitRe(:,:)*mycos(1,1)*mycos(2,1)
           nbit(2,it,:,:,1,1,0) = nbit(2,it,:,:,1,1,0) &
                                + bitIm(:,:)*mycos(1,1)*mycos(2,1)
           nbit(1,it,:,:,1,0,1) = nbit(1,it,:,:,1,0,1) &
                                + bitRe(:,:)*mycos(1,1)*mycos(3,1)
           nbit(2,it,:,:,1,0,1) = nbit(2,it,:,:,1,0,1) &
                                + bitIm(:,:)*mycos(1,1)*mycos(3,1)
           nbit(1,it,:,:,0,1,1) = nbit(1,it,:,:,0,1,1) &
                                + bitRe(:,:)*mycos(2,1)*mycos(3,1)
           nbit(2,it,:,:,0,1,1) = nbit(2,it,:,:,0,1,1) &
                                + bitIm(:,:)*mycos(2,1)*mycos(3,1)
           nbit(1,it,:,:,1,1,1) = nbit(1,it,:,:,1,1,1) &
                                + bitRe(:,:)*mycos(1,1)*mycos(2,1)*mycos(3,1)
           nbit(2,it,:,:,1,1,1) = nbit(2,it,:,:,1,1,1) &
                                + bitIm(:,:)*mycos(1,1)*mycos(2,1)*mycos(3,1)
           nbit(1,it,:,:,2,0,0) = nbit(1,it,:,:,2,0,0) + bitRe(:,:)*mycos(1,2)
           nbit(2,it,:,:,2,0,0) = nbit(2,it,:,:,2,0,0) + bitIm(:,:)*mycos(1,2)
           nbit(1,it,:,:,0,2,0) = nbit(1,it,:,:,0,2,0) + bitRe(:,:)*mycos(2,2)
           nbit(2,it,:,:,0,2,0) = nbit(2,it,:,:,0,2,0) + bitIm(:,:)*mycos(2,2)
           nbit(1,it,:,:,0,0,2) = nbit(1,it,:,:,0,0,2) + bitRe(:,:)*mycos(3,2)
           nbit(2,it,:,:,0,0,2) = nbit(2,it,:,:,0,0,2) + bitIm(:,:)*mycos(3,2)
           if (lessnoise==0) then
            nbit(1,it,:,:,1,0,0) = nbit(1,it,:,:,1,0,0) - bitIm(:,:)*mysin(1,1)
            nbit(2,it,:,:,1,0,0) = nbit(2,it,:,:,1,0,0) + bitRe(:,:)*mysin(1,1)
            nbit(1,it,:,:,0,1,0) = nbit(1,it,:,:,0,1,0) - bitIm(:,:)*mysin(2,1)
            nbit(2,it,:,:,0,1,0) = nbit(2,it,:,:,0,1,0) + bitRe(:,:)*mysin(2,1)
            nbit(1,it,:,:,0,0,1) = nbit(1,it,:,:,0,0,1) - bitIm(:,:)*mysin(3,1)
            nbit(2,it,:,:,0,0,1) = nbit(2,it,:,:,0,0,1) + bitRe(:,:)*mysin(3,1)
            nbit(1,it,:,:,1,1,0) = nbit(1,it,:,:,1,1,0) &
                                 - bitRe(:,:)*mysin(1,1)*mysin(2,1) &
                                 - bitIm(:,:)*(mycos(1,1)*mysin(2,1) &
                                              +mysin(1,1)*mycos(2,1))
            nbit(2,it,:,:,1,1,0) = nbit(2,it,:,:,1,1,0) &
                                 - bitIm(:,:)*mysin(1,1)*mysin(2,1) &
                                 + bitRe(:,:)*(mycos(1,1)*mysin(2,1) &
                                              +mysin(1,1)*mycos(2,1))
            nbit(1,it,:,:,1,0,1) = nbit(1,it,:,:,1,0,1) &
                                 - bitRe(:,:)*mysin(1,1)*mysin(3,1) &
                                 - bitIm(:,:)*(mycos(1,1)*mysin(3,1) &
                                              +mysin(1,1)*mycos(3,1))
            nbit(2,it,:,:,1,0,1) = nbit(2,it,:,:,1,0,1) &
                                 - bitIm(:,:)*mysin(1,1)*mysin(3,1) &
                                 + bitRe(:,:)*(mycos(1,1)*mysin(3,1) &
                                              +mysin(1,1)*mycos(3,1))
            nbit(1,it,:,:,0,1,1) = nbit(1,it,:,:,0,1,1) &
                                 - bitRe(:,:)*mysin(2,1)*mysin(3,1) &
                                 - bitIm(:,:)*(mycos(2,1)*mysin(3,1) &
                                              +mysin(2,1)*mycos(3,1))
            nbit(2,it,:,:,0,1,1) = nbit(2,it,:,:,0,1,1) &
                                 - bitIm(:,:)*mysin(2,1)*mysin(3,1) &
                                 + bitRe(:,:)*(mycos(2,1)*mysin(3,1) &
                                              +mysin(2,1)*mycos(3,1))
            nbit(1,it,:,:,1,1,1) = nbit(1,it,:,:,1,1,1)                        &
                                - bitRe(:,:)*(mysin(1,1)*mysin(2,1)*mycos(3,1) &
                                             +mysin(1,1)*mycos(2,1)*mysin(3,1) &
                                             +mycos(1,1)*mysin(2,1)*mysin(3,1))&
                                - bitIm(:,:)*(mycos(1,1)*mycos(2,1)*mysin(3,1) &
                                             +mycos(1,1)*mysin(2,1)*mycos(3,1) &
                                             +mysin(1,1)*mycos(2,1)*mycos(3,1) &
                                             -mysin(1,1)*mysin(2,1)*mysin(3,1))
            nbit(2,it,:,:,1,1,1) = nbit(2,it,:,:,1,1,1)                        &
                                - bitIm(:,:)*(mysin(1,1)*mysin(2,1)*mycos(3,1) &
                                             +mysin(1,1)*mycos(2,1)*mysin(3,1) &
                                             +mycos(1,1)*mysin(2,1)*mysin(3,1))&
                                + bitRe(:,:)*(mycos(1,1)*mycos(2,1)*mysin(3,1) &
                                             +mycos(1,1)*mysin(2,1)*mycos(3,1) &
                                             +mysin(1,1)*mycos(2,1)*mycos(3,1) &
                                             -mysin(1,1)*mysin(2,1)*mysin(3,1))
            nbit(1,it,:,:,2,0,0) = nbit(1,it,:,:,2,0,0) - bitIm(:,:)*mysin(1,2)
            nbit(2,it,:,:,2,0,0) = nbit(2,it,:,:,2,0,0) + bitRe(:,:)*mysin(1,2)
            nbit(1,it,:,:,0,2,0) = nbit(1,it,:,:,0,2,0) - bitIm(:,:)*mysin(2,2)
            nbit(2,it,:,:,0,2,0) = nbit(2,it,:,:,0,2,0) + bitRe(:,:)*mysin(2,2)
            nbit(1,it,:,:,0,0,2) = nbit(1,it,:,:,0,0,2) - bitIm(:,:)*mysin(3,2)
            nbit(2,it,:,:,0,0,2) = nbit(2,it,:,:,0,0,2) + bitRe(:,:)*mysin(3,2)
           endif
          enddo ! ieo

! End main loop.
         enddo ! ibl
        enddo ! ixbit
       enddo ! iybit
      enddo ! izbit
     enddo ! itbit

! Close the necessary tmU files.
     do jc = 1,3
      do jd = 1,4
       ifile = 100 + jd + (jc-1)*4
       close(unit=ifile,status="keep")
      enddo ! jd
     enddo ! jc

! Close the necessary tmD files.
     do jc = 1,3
      do jd = 1,4
       ifile = 200 + jd + (jc-1)*4
       close(unit=ifile,status="keep")
      enddo ! jd
     enddo ! jc

! Combine the results from all processes.
     pcorr = real(pbit,KR)
     rcorr = real(rbit,KR)
     ncorr = real(nbit,KR)
     if (nps==1) then
      pion = pcorr
      rho = rcorr
      nucl = ncorr
     else
      nmesg = 125*nt
      call MPI_REDUCE(pcorr(1,0,0,0),pion(1,0,0,0),nmesg,MRT,MPI_SUM,0, &
                      MPI_COMM_WORLD,ierr)
      call MPI_BCAST(pion(1,0,0,0),nmesg,MRT,0,MPI_COMM_WORLD,ierr)
      call MPI_REDUCE(rcorr(1,0,0,0),rho(1,0,0,0),nmesg,MRT,MPI_SUM,0, &
                      MPI_COMM_WORLD,ierr)
      call MPI_BCAST(rho(1,0,0,0),nmesg,MRT,0,MPI_COMM_WORLD,ierr)
      nmesg = 32*27*nt
      call MPI_REDUCE(ncorr(1,1,1,1,0,0,0),nucl(1,1,1,1,0,0,0),nmesg,MRT, &
                      MPI_SUM,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(nucl(1,1,1,1,0,0,0),nmesg,MRT,0,MPI_COMM_WORLD,ierr)
     endif

! Record the results.
     if (myid==0) then
        ntrailer = "xxxx.NUCLEON.LOG"
        write(unit=ntrailer(1:4), fmt="(1i4.4)") icfgsave
        opfile = trim(rwdir(myid+1))//trim(ntrailer)
        inquire(file=opfile, exist=fileexists)
        if (.not. fileexists) then
           open(unit=11,file=opfile, &
                action="write",form="formatted",status="new")
        else
           open(unit=11,file=opfile, &
                action="write",form="formatted",status="old",position="append")
        endif ! file check

       do iz = 0,4
        do iy = 0,4
         do jx = 0,4
          write(unit=11,fmt="(i1,a22,i1,a1,i1,a1,i1,a5,500es17.10)") &
              ikappa,"charged pion      [p=(",jx,",",iy,",",iz,")] = ", pion(:,jx,iy,iz)
          write(unit=11,fmt="(i1,a14,i1,a1,i1,a1,i1,a5,500es17.10)") &
              ikappa,"rho       [p=(",jx,",",iy,",",iz,")] = ", rho(:,jx,iy,iz)

          if (jx**2+iy**2+iz**2<5) then
             if(particle==1) then
                particlestring = " uud "
             else
                particlestring = " ddu "
             endif ! particle
             do jd = 1,4
              do id = 1,4
             write(unit=11,fmt="(i1,a7,2i1,a6,a5,i1,a1,i1,a1,i1,a5,500es17.10)") &
                   ikappa,"nuclRe,",jd,id,trim(particlestring)," [p=(",jx,",",iy,",",iz,")] = ",  &
                   nucl(1,:,jd,id,jx,iy,iz)
             write(unit=11,fmt="(i1,a7,2i1,a6,a5,i1,a1,i1,a1,i1,a5,500es17.10)") &
                   ikappa,"nuclIm,",jd,id,trim(particlestring)," [p=(",jx,",",iy,",",iz,")] = ",  &
                   nucl(2,:,jd,id,jx,iy,iz)
            enddo ! id
           enddo ! jd
          endif
         enddo ! jx
        enddo ! iy
       enddo ! iz
      close(unit=11,status="keep")
     endif

! COMMENT ~ NUCLEON (only for debug!!!)
! We print out the jd=1,2 values corresponding to the foward parity 
! proton determined by the GAMMA_4 matrix (in the Wilson basis).
! Since we are doing the rotation from the twisted basis to the physical basis,
! we need to use (gamma_5*GAMMA_4*gamma_5) in the nucleon correlation function
! which makes use of jd=3,4....picture

! For the proton...

    if (particle==1) then 
!
!                                        / 1  1 \
!     (1+ gamma_5)*GAMMA_4*(1+gamma_5) = |      |
!                                        \-1 -1 /

! DEAN ~
! I printed out the values and did the calclation by hand and it agrees with the 
! GAMMA_4(new) multiplication done below!

!   if(myid==0) then
!     do id=1,4
!       do jd=1,4
!         print *, "nucl(1,4,id,jd,0,0,0)=",id,jd, nucl(1,4,id,jd,0,0,0) 
!       enddo ! jd
!     enddo ! id
!   endif ! myid

! This adds the first and third rows together as the above multiplication suggests. 
! Then I dump the CHANGED first column into the third since then only differ by a 
! minus sign. (This is not an error, but done intentionaly) 

! CAUTION:
! This is a test for the twMass on the nucleon. IF the propagator is NOT twisted
! you may twist the Gamma_4 matrix and see the same result. THIS PROGRAM IS TWISTING
! PROPAGATORS so do not use. 

 if (.false.) then
     do ittt =1,nt
       do row = 1,2
          row2 = row + 2

          nucl(1,ittt,row,:,0,0,0) = nucl(1,ittt,row,:,0,0,0) + nucl(1,ittt,row2,:,0,0,0)
          nucl(2,ittt,row,:,0,0,0) = nucl(2,ittt,row,:,0,0,0) + nucl(2,ittt,row2,:,0,0,0)
          nucl(1,ittt,row2,:,0,0,0) = -1.0_KR*nucl(1,ittt,row,:,0,0,0)
          nucl(2,ittt,row2,:,0,0,0) = -1.0_KR*nucl(2,ittt,row,:,0,0,0)

       enddo ! row
     enddo ! ittt
 
     fac = 0.5_KR
     nucl = fac*nucl
 endif ! .false.

     if (myid==0) then
      open(unit=12,file=trim(rwdir(myid+1))//"DUBLIN.OUT",action="write", &
           form="formatted",status="old",position="append")
        write(unit=12,fmt="(a12)") particlestring
        znuc = 0.0_KR
        do jd=1,2
          do ittt = 1,nt
              znuc(1,ittt) = znuc(1,ittt) + nucl(1,ittt,jd,jd,0,0,0)
              znuc(2,ittt) = znuc(2,ittt) + nucl(2,ittt,jd,jd,0,0,0)
          enddo ! ittt
        enddo ! jd

        do ittt=1,nt
          write(unit=12,fmt="(i5,2es17.10)")  ittt,znuc(1,ittt),znuc(2,ittt)
        enddo ! ittt
        close(unit=12,status="keep")
     endif

    else !  particle
! For the neutron
!
!                                        / 1 -1 \
!     (1- gamma_5)*GAMMA_4*(1-gamma_5) = |      |
!                                        \ 1 -1 /


 if (.false.) then
     do ittt =1,nt
       do row = 1,2
         row2 = row + 2

         nucl(1,ittt,row,:,0,0,0)    = nucl(1,ittt,row,:,0,0,0) - nucl(1,ittt,row2,:,0,0,0)
         nucl(2,ittt,row,:,0,0,0)    = nucl(2,ittt,row,:,0,0,0) - nucl(2,ittt,row2,:,0,0,0)
         nucl(1,ittt,row2,:,0,0,0) = nucl(1,ittt,row,:,0,0,0)
         nucl(2,ittt,row2,:,0,0,0) = nucl(2,ittt,row,:,0,0,0)

       enddo ! row
     enddo ! ittt

     fac = 0.5_KR
     nucl = fac*nucl
 endif ! .false.
 
     if (myid==0) then
      open(unit=12,file=trim(rwdir(myid+1))//"DUBLIN.OUT",action="write", &
           form="formatted",status="old",position="append")
        write(unit=12,fmt="(a12)") particlestring
        znuc = 0.0_KR
        do jd=1,2
            do ittt = 1,nt
              znuc(1,ittt) = znuc(1,ittt) + nucl(1,ittt,jd,jd,0,0,0)
              znuc(2,ittt) = znuc(2,ittt) + nucl(2,ittt,jd,jd,0,0,0)
            enddo ! ittt
        enddo ! jd

        do ittt=1,nt
          write(unit=12,fmt="(i5,2es17.10)")  ittt,znuc(1,ittt),znuc(2,ittt)
        enddo ! ittt
        close(unit=12,status="keep")
     endif

    endif ! particle

    enddo ! iLS
    enddo ! ikappa

 end subroutine nucleoncorr

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 subroutine FPwrite(xe,xo,vecblinv,filename,kappa,tmqcd,cSW,icsrc,idsrc,icfg, &
                    myid)
! Writes a fermion propagator to disk.
! isite is the outermost do loop, so the writing is organized by timeslice.
! (Actually by pairs of timeslices, since
!       it=1 is uniquely defined by        1<isite<  factor and 1<ibl< 8
!       it=2 is uniquely defined by        1<isite<  factor and 9<ibl<16
!       it=3 is uniquely defined by   factor<isite<2*factor and 1<ibl< 8
!       it=4 is uniquely defined by   factor<isite<2*factor and 9<ibl<16
!       it=5 is uniquely defined by 2*factor<isite<3*factor and 1<ibl< 8
!       it=6 is uniquely defined by 2*factor<isite<3*factor and 9<ibl<16
!       etc, where factor = 2*nvhalf*npt/nt.
! INPUT:
!   xe is the ibleo=1 half of the propagator.
!   xo is the ibleo=2 half of the propagator.
!   vecbl(1,1:8) = (/ 1, 4, 6, 7, 10, 11, 13, 16 /)
!   vecbl(2,1:8) = (/ 2, 3, 5, 8,  9, 12, 14, 15 /)
!   vecblinv(1,1:16) = (/ 1, 2, 2, 1, 2, 1, 1, 2, 2, 1, 1, 2, 1, 2, 2, 1 /)
!   vecblinv(2,1:16) = (/ 1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7, 8, 8 /)
!   filename is the file to be created.
!   kappa is the hopping parameter (or "m" for twisted mass QCD).
!   tmqcd is "mu" for twisted mass QCD (unused for standard QCD).
!   cSW is the clover coefficient in the fermion action.
!   icsrc is the colour of the point source.
!   idsrc is the Dirac index of the point source.
!   icfg is the integer which has been used to identify the gauge configuration.

    integer(kind=KI), intent(in)                          :: icsrc, idsrc, icfg, myid
    real(kind=KR),    intent(inout), dimension(:,:,:,:,:) :: xe, xo
    integer(kind=KI), intent(in),    dimension(:,:)       :: vecblinv
    character(len=*), intent(in)                          :: filename
    real(kind=KR),    intent(in)                          :: kappa, tmqcd, cSW

    integer(kind=KI) :: ibl, ivbl, ibleo, isite, icri, id

    open(unit=9,file=trim(filename),action="write",status="replace", &!change
         form="unformatted")
     do isite = 1,nvhalf
      do ibl = 1,16
       ivbl = vecblinv(1,ibl)
       ibleo = vecblinv(2,ibl)
       do icri = 1,6
        if (ivbl==1) then
         write(unit=9) (xe(icri,isite,id,1,ibleo), &
                        xe(icri,isite,id,2,ibleo),id=1,4)
        else
         write(unit=9) (xo(icri,isite,id,1,ibleo), &
                        xo(icri,isite,id,2,ibleo),id=1,4)
        endif
       enddo ! icri
      enddo ! ibl
     enddo ! isite
     write(unit=9) kappa, tmqcd, cSW, icsrc, idsrc, icfg, myid
    close(unit=9,status="keep")

 end subroutine FPwrite

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 subroutine FPread(xe,xo,vecblinv,filename,kappa,tmqcd,cSW,icsrc,idsrc,icfgin,&
                   myidin)
! Read a fermion propagator from disk.
! INPUT:
!   filename is the file to be read in.
! OUTPUT:
!   x is the fermion propagator.
!   kappa is the hopping parameter (or "m" for twisted mass QCD).
!   tmqcd is "mu" for twisted mass QCD (unused for standard QCD).
!   cSW is the clover coefficient in the fermion action.
!   icsrc is the colour of the point source.
!   idsrc is the Dirac index of the point source.
!   icfgin is the integer used to identify the gauge configuration.
!   myidin is "myid" for the process that created this partial configuration.

    real(kind=KR),    intent(inout), dimension(:,:,:,:,:) :: xe, xo
    integer(kind=KI), intent(in),    dimension(:,:)       :: vecblinv
    character(len=*), intent(in)                          :: filename
    real(kind=KR),    intent(out)                         :: kappa, tmqcd, cSW
    integer(kind=KI), intent(out)               :: icsrc, idsrc, icfgin, myidin

    integer(kind=KI) :: ibl, isite, icri, id, ivbl, ibleo

    open(unit=9,file=trim(filename),action="read",status="old", &
         position="rewind",form="unformatted")
     do isite = 1,nvhalf
      do ibl = 1,16
       ivbl = vecblinv(1,ibl)
       ibleo = vecblinv(2,ibl)
       do icri = 1,6
        if (ivbl==1) then
         read(unit=9) (xe(icri,isite,id,1,ibleo), &
                       xe(icri,isite,id,2,ibleo),id=1,4)
        else
         read(unit=9) (xo(icri,isite,id,1,ibleo), &
                       xo(icri,isite,id,2,ibleo),id=1,4)
        endif
       enddo ! icri
      enddo ! ibl
     enddo ! isite
     read(unit=9) kappa, tmqcd, cSW, icsrc, idsrc, icfgin, myidin
    close(unit=9,status="keep")

 end subroutine FPread

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 subroutine gamma5vector(a,ntmqcd,myid)

! Multiplies vector a by -i*gamma5.

! This subroutine is employed to move the "twisted operators" into
! the physical basis. For example we it can be shown in the 
! twisted formalism that the neutral pion and the charged pion need
! to be rotated by the exp[+/-igamma_5w/2] to have the same meaning
! that it had in the physical basis. The plus and minus sign
! are a result of the tau_3 doublet. 

! The form of the rotation is:

!        a -> (1+/-igamma_5)a

! where the +/- represnts the tmU and tmD.

! We are assuming MAIXIMAL TWIST with this (1 +/- i gamma_5) term
! such that the twist angle is exactly pi/4

! INPUT: vecsrc is the input e/o preconditioned RHS.
!        angle is the maximal twist angle
!        ntmqcd is the parameter that signals tmU or tmD
! OUTPUT: vecsrc -> exp[+/-igamma_5w/2]*vecsrc

 

    real(kind=KR),    intent(inout),   dimension(6,ntotal,4,2,8)   :: a
    integer(kind=KI), intent(in)                                   :: ntmqcd, myid

    real(kind=KR),                     dimension(6,ntotal,2,2,8)   :: temp
    real(kind=KR)                                                  :: fac
    integer(kind=KI)                                               :: isite


    temp = 0.0_KR
    fac = 1.0_KR/(sqrt(2.0_KR))


    if(ntmqcd>0) then
! This will do the tmU case and multiple a by (1+igamma_5)
      do isite = 1,ntotal
       temp(:,isite,1,:,:) = a(:,isite,1,:,:)
       temp(:,isite,2,:,:) = a(:,isite,2,:,:)
       a(:,isite,1,:,:) = a(:,isite,1,:,:) + a(:,isite,3,:,:)
       a(:,isite,2,:,:) = a(:,isite,2,:,:) + a(:,isite,4,:,:)
       a(:,isite,3,:,:) = a(:,isite,3,:,:) - temp(:,isite,1,:,:)
       a(:,isite,4,:,:) = a(:,isite,4,:,:) - temp(:,isite,2,:,:)
      enddo ! isite
    a = fac*a
    else if(ntmqcd<0) then
! This will do the tmD case and multiple a by (1-igamma_5)
      do isite = 1,ntotal
       temp(:,isite,1,:,:) = a(:,isite,1,:,:)
       temp(:,isite,2,:,:) = a(:,isite,2,:,:)
       a(:,isite,1,:,:) = a(:,isite,1,:,:) - a(:,isite,3,:,:)
       a(:,isite,2,:,:) = a(:,isite,2,:,:) - a(:,isite,4,:,:)
       a(:,isite,3,:,:) = a(:,isite,3,:,:) + temp(:,isite,1,:,:)
       a(:,isite,4,:,:) = a(:,isite,4,:,:) + temp(:,isite,2,:,:)
      enddo ! isite
    a = fac*a
    endif ! ntmqcd

 
 end subroutine gamma5vector
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 subroutine bsource(bvec,be,bo,kappa,myid)
! This subroutine combines the even and odd parts of the source vector
!    and outputs in bsource to form ....
!
!                                      /be\
!      bsource = -i*gamma5/(2*kappa*mu)|  | ,
!                                      \bo/
!
! bsource has the same first five indices as be/bo. The last index corresponds
!     to the even(1) and odd(2) source vectors.

    real(kind=KR),    intent(out),   dimension(:,:,:,:,:,:) :: bvec
    real(kind=KR),    intent(in),    dimension(:,:,:,:,:)   :: be
    real(kind=KR),    intent(in),    dimension(:,:,:,:,:)   :: bo
    real(kind=KR),    intent(in),    dimension(:)           :: kappa
    integer(kind=KI), intent(in)                            :: myid

    real(kind=KR),              dimension(6,nvhalf,2,2,8,2) :: temp
    real(kind=KR)                                           :: fac1
    integer(kind=KI)  :: ii

    bvec =0.0_KR
    do ii= 1,nvhalf
       bvec(:,ii,:,:,:,1) = be(:,ii,:,:,:)
       bvec(:,ii,:,:,:,2) = bo(:,ii,:,:,:)
    enddo ! ii

! Do a normilization by gamma5/2*i*kappa*mu or (-i*gamma5)/2.0*kappa*mu

    fac1 = 1.0_KR/(2.0_KR*kappa(1)*kappa(2))

    if(myid==0) then
     print *, "fac1=", fac1
    endif ! myid

    temp = 0.0_KR
    do ii =1,nvhalf
      temp(:,ii,1,:,:,:) = bvec(:,ii,1,:,:,:)
      temp(:,ii,2,:,:,:) = bvec(:,ii,2,:,:,:)

      bvec(:,ii,1,:,:,:) = -fac1*bvec(:,ii,3,:,:,:)
      bvec(:,ii,2,:,:,:) = -fac1*bvec(:,ii,4,:,:,:)
      bvec(:,ii,3,:,:,:) =  fac1*temp(:,ii,1,:,:,:)
      bvec(:,ii,4,:,:,:) =  fac1*temp(:,ii,2,:,:,:)
    enddo ! ii
 end subroutine bsource

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 subroutine xvector(xwhole,xe,xo)

! This subroutine combines the even and odd parts of the solutuion vector
!    and outputs in bsource to form ....
!
!                /xe\
!      xwhole =  |  | ,
!                \xo/
!
! xwhole has the same indices as xe/xo excpet the 6th index. That index corresponds
!     to the even(1) and odd(2) solution vectors.

    real(kind=KR),    intent(out), dimension(:,:,:,:,:,:,:)   :: xwhole
    real(kind=KR),    intent(in),  dimension(:,:,:,:,:,:)     :: xe, xo

    integer(kind=KI)     :: ii,is

!HEY, HEY ,HEY DEAN
! Dean - Put an explicit loop here on last index
! Only need to do first shift here....
    
     xwhole =0.0_KR

     do is = 1,nshifts
      do ii = 1,nvhalf
        xwhole(:,ii,:,:,:,1,is) = xe(:,ii,:,:,:,is)
        xwhole(:,ii,:,:,:,2,is) = xo(:,ii,:,:,:,is)
      enddo ! ii
     enddo ! is
 end subroutine xvector

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 subroutine check4ones(xe,myid)

    real(kind=KR), intent(in), dimension(6,ntotal,4,2,8)     :: xe
    integer(kind=KI), intent(in) ::  myid
    integer(kind=KI) :: i,j,k,l,m,ic, ierr

    if (myid == 0) then
      print *, "check4ones"
      do i =1,6
         do j =1,nvhalf
            do k=1,4
               do l=1,2
                  do m =1,8
                     if(abs(xe(i,j,k,l,m)) /= 0.0_KR) then
                        print "(a32,6i3,1es17.10)", "myid,xe(i,j,k,l,m)=", myid,i,j,k,l,m,xe(i,j,k,l,m)
                     endif
                  enddo ! m
               enddo ! l
           enddo ! k
         enddo ! j
      enddo ! i
     end if ! myid
       call MPI_BARRIER(MPI_COMM_WORLD, ierr)
       stop
   end subroutine check4ones

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 subroutine printlog(s,myid,rwdir)

 character(len=*), intent(in)                  :: s
 character(len=*), intent(in), dimension(:)    :: rwdir
 integer(kind=KI), intent(in)                  :: myid
      if (myid==0) then
        open(unit=8,file=trim(rwdir(myid+1))//"CFGSPROPS.LOG",     &
             action="write",form="formatted",status="old",position="append")
        write(unit=8,fmt=*)  s
        close(unit=8,status="keep")
       endif

  end subroutine printlog

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 end module quark
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

