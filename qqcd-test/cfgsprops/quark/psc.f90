 subroutine psc(rwdir,b,x,GMRES,resmax,itermin,itercount,u,GeeGooinv, &
                    iflag,kappa,coact,bc,vecbl,vecblinv,myid,nn,ldiv,nms, &
                    lvbc,ib,lbd,iblv,MRT,MRT2)
! GMRES-DR(n,k) matrix inverter of
! Ronald B. Morgan, SIAM Journal on Scientific Computing, 24, 20 (2002).
! Solves M*x=b for the vector x.
! INPUT:
!   b() is the source vector.
!   x() is used as an initial estimate of the true solution vector.
!   GMRES(1)=n in GMRES-DR(n,k): maximum dimension of the subspace.
!   GMRES(2)=k in GMRES-DR(n,k): number of approx eigenvectors kept at restart.
!   resmax is the stopping criterion for the iteration.
!   itermin is the minimum number of iteration required by the user.
!   u() contains the gauge fields for this sublattice.
!   GeeGooinv contains the clover/twisted-mass matrix G on globlly-even sites
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
!   itercount is the number of iterations used by gmresdrshift.
 
    use shift 

    character(len=*), intent(in),    dimension(:)             :: rwdir
    real(kind=KR),    intent(in),    dimension(:,:,:,:,:)     :: b 
    real(kind=KR),    intent(inout), dimension(:,:,:,:,:)     :: x
    integer(kind=KI), intent(in),    dimension(:)             :: GMRES, bc, nms
    integer(kind=KI), intent(in)                              :: itermin, iflag, &
                                                                 myid, MRT, MRT2
                                                             
    real(kind=KR),    intent(in)                              :: resmax
    integer(kind=KI), intent(out)                             :: itercount
    real(kind=KR),    intent(in),    dimension(:,:,:,:,:)     :: u
    real(kind=KR),    intent(in),    dimension(:,:,:,:,:)     :: GeeGooinv
    real(kind=KR),    intent(in),    dimension(:)             :: kappa
    real(kind=KR),    intent(in),    dimension(:,:,:)         :: coact
    integer(kind=KI), intent(in),    dimension(:,:)           :: vecbl, vecblinv
    integer(kind=KI), intent(in),    dimension(:,:)           :: nn, iblv
    logical,          intent(in),    dimension(:)             :: ldiv
    integer(kind=KI), intent(in),    dimension(:,:,:)         :: lvbc
    integer(kind=KI), intent(in),    dimension(:,:,:,:)       :: ib
    logical,          intent(in),    dimension(:,:)           :: lbd
  
    integer(kind=KI), dimension(2)                            :: GMRESproj
    real(kind=KR),    dimension(6,ntotal,4,2,8,nshifts)       :: tempxvshift
    real(kind=KR2)                                            :: rn, rnc, rnt, rrr
    real(kind=KR2),   dimension(2)                            :: gam
    real(kind=KR2),   dimension(2,nshifts)                    :: gamtemp
    integer(kind=KI)                                          :: mvp, nDR, kDR, is, &
                                                                 icri, jj, idag, i, k, &
                                                                 ishift, ierr, j,l,m,n
    real(kind=KR2),   dimension(6,nvhalf,4,2,8,nshifts)       :: rshift
    real(kind=KR),    dimension(6,ntotal,4,2,8,nshifts)       :: xshift

! NOTE~ the dimension of gdr must be big enough to hold all the norms
!       for the specified size of matrix

    real(kind=KR),    dimension(gdrsize,nshifts)              :: gdr 
!   real(kind=KR),    dimension(15000,nshifts)              :: gdr 
    real(kind=KR)                                             :: projresmax
    real(kind=KR),    dimension(2)                            :: beta
    real(kind=KR2),   dimension(nshifts)                      :: sigma

! This is still PSC

! Shift sigmamu to base mu above  (mtmqcd(1,2))
      print *,'from psc inside quark'

! Inititilizations  
 
! NOTE ~ BIG NOTE - might need to pass in idag into psc if that part of gamma5mult is needed

    idag = 0
    mvp = 0
    gdr = 0.0_KR
      v = 0.0_KR
    

    nDR = GMRES(1)
    kDR = GMRES(2)

    if (isignal == 1) then

      xshift(:,:,:,:,:,1) = x(:,:,:,:,:)
      vtemp = 0.0_KR2
      v = 0.0_KR


      call gmresdr(rwdir,b,xshift(:,:,:,:,:,1),GMRES,resmax,itermin,itercount,u,GeeGooinv, &
                    iflag,kappa,coact,bc,vecbl,vecblinv,myid,nn,ldiv,nms, &
                    lvbc,ib,lbd,iblv,MRT,MRT2)

    !BS change for polynomial substraction  call gmresdr(rwdir,b,xshift(:,:,:,:,:,1),GMRES,resmax,itermin,itercount,u,GeeGooinv, &
                !BS    iflag,kappa,coact,bc,vecbl,vecblinv,myid,nn,ldiv,nms, &
                !BS    lvbc,ib,lbd,iblv,MRT,MRT2)

      beg(:,1:nvhalf,:,:,:) = vtemp(:,1:nvhalf,:,:,:,kDR+1)
     
    else ! (isignal == 1)

    itercount = 0
    projresmax = 1.0e-7

    call vecdot(b(:,:,:,:,:), b(:,:,:,:,:), beta, MRT2)


! For comparison with gmresdr(m,k) we need to compare with  
!     gmres(m-k) - proj(k)

    GMRESproj(1) = GMRES(1) - GMRES(2)
    GMRESproj(2) = GMRES(2)

print *,'GMRESproj(1)',GMRES(1)-GMRES(2)
print *,'GMRESproj(2)',GMRES(2)
    call gmresproject(rwdir,b,xshift,GMRESproj,projresmax,itermin,itercount,u,GeeGooinv, &
                     iflag,kappa,coact,bc,vecbl,vecblinv,myid,nn,ldiv,nms, &
                     lvbc,ib,lbd,iblv,MRT,MRT2,isignal,mvp,gdr)

!BS for polynomial substraction    call gmresproject(rwdir,b,xshift,GMRESproj,projresmax,itermin,itercount,u,GeeGooinv, &
                   !BS  iflag,kappa,coact,bc,vecbl,vecblinv,myid,nn,ldiv,nms, &
                    !BS lvbc,ib,lbd,iblv,MRT,MRT2,isignal,mvp,gdr)

    endif ! (isignal /= 1)

    x(:,:,:,:,:) = xshift(:,:,:,:,:,1)

end subroutine psc

