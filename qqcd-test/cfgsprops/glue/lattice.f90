! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
! lattice.f90, Randy Lewis, randy.lewis@uregina.ca
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
! Routines for the initialization of a lattice of quenched QCD gauge fields.
!
! parts of this code are based on:
! *QCDMPI (version 5, 1996) by Shinji Hioki, Parallel Computing 22 (1997) 1335.
! and
! *QCDimMPI (version 0.9, 1998) by Astushi Nakamura and Shinji Hioki,
!           Nucl Phys B (Proc Suppl) 73 (1999) 895.
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 module lattice

    use kinds
    use latdims
    use basics
    implicit none
    private

! Define access to subroutines.
    public  :: initlat, atoc, ctoa, fixaction
    private :: setnodeid, siteindex, setblk, gblchker

 contains

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 subroutine initlat(myid,bc,nn,ldiv,nms,lv,lvbc,ib,lbd,iblv,vecbl,vecblinv)

    integer(kind=KI), intent(in)                      :: myid
    integer(kind=KI), intent(in),  dimension(:)       :: bc
    integer(kind=KI), intent(out), dimension(:,:)  :: nn, iblv, vecbl, vecblinv
    logical,          intent(out), dimension(:)       :: ldiv
    integer(kind=KI), intent(out), dimension(:)       :: nms
    integer(kind=KI), intent(out), dimension(:,:,:)   :: lv, lvbc
    integer(kind=KI), intent(out), dimension(:,:,:,:) :: ib
    logical,          intent(out), dimension(:,:)     :: lbd

    integer(kind=KI), dimension(4) :: ip

!*Identify the sublattice to be assigned to myid and neighbouring processes.
    call setnodeid(myid,ip,nn)
! ip(i) = coordinates of the current process on the grid of processes

!*Indexing.
! Set up the indexing of lattice sites, nearest neighbours in the positive
! spacetime directions, and process boundaries.
    call siteindex(ip,bc,ldiv,nms,lv,lvbc,ib)
! Set up an index for the various blocks of a sublattice.
    call setblk(lbd,iblv)
! Set up global checkboarding (used for fermion inversion; not for Monte Carlo).
    call gblchker(vecbl,vecblinv)

 end subroutine initlat

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 subroutine setnodeid(myid,ip,nn)
! Identify the sublattice to be assigned to myid.
! OUTPUT:
!   np(i) = number of processes in the i'th spacetime direction
!   ip(i) = coordinates of the current process on the grid of processes
!         expected size: ip(4)
!   nn(j,1) = single-integer address, on the grid of processes, of the
!             neighbour to myid in the +j direction.
!   nn(j,2) = single-integer address, on the grid of processes, of the
!             neighbour to myid in the -j direction.
!           expected size: nn(4,2)
 
    integer(kind=KI), intent(in)                  :: myid
    integer(kind=KI), intent(out), dimension(:)   :: ip
    integer(kind=KI), intent(out), dimension(:,:) :: nn

    integer(kind=KI)               :: ipsave, mu
    integer(kind=KI), dimension(4) :: np
 
    np(1) = npx
    np(2) = npy
    np(3) = npz
    np(4) = npt
    call atoc(myid,np,ip)
    do mu = 1,4
     ipsave = ip(mu)
     ip(mu) = modulo(ipsave+1,np(mu))
     call ctoa(ip,np,nn(mu,1))
     ip(mu) = modulo(ipsave-1+np(mu),np(mu))
     call ctoa(ip,np,nn(mu,2))
     ip(mu) = ipsave
    enddo ! mu
 
 end subroutine setnodeid
 
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 subroutine atoc(ia,nc,ic)
! Convert single-integer address "ia" of a lattice site to spacetime
! coordinates "ic" in a spacetime of size "nc"
! OR
! Convert single-integer address "ia" of a process to 4-D coordinates
! "ic" on a grid of processes of size "nc"
! NOTE: coordinates and addresses are numbered 0, 1, 2, ...
!       expected sizes: ic(4), nc(4)
 
    integer(kind=KI), intent(in)                :: ia
    integer(kind=KI), intent(in),  dimension(:) :: nc
    integer(kind=KI), intent(out), dimension(:) :: ic
 
    ic(1) = modulo(ia,nc(1))
    ic(2) = modulo(ia/nc(1),nc(2))
    ic(3) = modulo(ia/(nc(1)*nc(2)),nc(3))
    ic(4) = modulo(ia/(nc(1)*nc(2)*nc(3)),nc(4))
 
 end subroutine atoc
 
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 
 subroutine ctoa(ic,nc,ia)
! Convert spacetime coordinates "ic" in a spacetime of size "nc" to
! single-integer address "ia"
! OR
! Convert 4-D coordinates "ic" on a grid of processes of size "nc" to
! single-integer address "ia"
! NOTE: coordinates and addresses are numbered 0, 1, 2, ...
!       expected sizes: ic(4), nc(4)
 
    integer(kind=KI), intent(in), dimension(:) :: ic, nc
    integer(kind=KI), intent(out)              :: ia
 
    ia = ic(1) + ic(2)*nc(1) + ic(3)*nc(2)*nc(1) + ic(4)*nc(3)*nc(2)*nc(1)
 
 end subroutine ctoa
 
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 subroutine siteindex(ip,bc,ldiv,nms,lv,lvbc,ib)
! Create an index for the lattice sites within one block of a sublattice.
! INPUT:
!   ip(mu) = number between 0 and np(mu); position of this process in
!            mu'th direction.
!   bc(mu) = 1,-1,0 for periodic,antiperiodic,fixed boundary conditions
!            respectively.
! OUTPUT:
!   ns(1)=n1, ns(2)=n2, ns(3)=n3, ns(4)=n4
!   np(1)=npx, np(2)=npy, np(3)=npz, np(4)=npt
!   isg(iv) is the site index for this block of the sublattice.
!           It counts 0, 1, 2,... separately for the even and the odd sites.
!           expected size: isg(nv)
!   ldiv(mu) is true if there is more than one process in the mu'th direction,
!            otherwise it is false.
!   lv(ivhalf,mu,ieo) is the site of nearest neighbour in +mu direction.
!                     expected size: lv(nvhalf,ndim,neo)
!   lvbc(ivhalf,mu,ieo) is the site of the nearest neighbour in +mu direction.
!                       If there is only one process in the mu direction,
!                       then lvbc still points to the buffer whenever
!                       non-periodic boundary conditions are needed.
!   ib(ibmax,mu,1,ieo) contains the boundary sites at the -mu edge of this
!                      block of the sublattice.
!   ib(ibmax,mu,2,ieo) contains the boundary sites at the +mu edge of this
!                      block of the sublattice.
!                      expected size: ib(nbmax,4,2,2)
!
! NOTE: coordinates and addresses are numbered 0, 1, 2, ...
! NOTE: A site (x,y,z,t) is "even" (ieo=1) if x+y+z+t is even,
!       otherwise it is odd (ieo=2).
 
    integer(kind=KI), intent(in),  dimension(:)       :: ip, bc
    logical,          intent(out), dimension(:)       :: ldiv
    integer(kind=KI), intent(out), dimension(:)       :: nms
    integer(kind=KI), intent(out), dimension(:,:,:)   :: lv, lvbc
    integer(kind=KI), intent(out), dimension(:,:,:,:) :: ib

    integer(kind=KI)                      :: iv, mu, ig, ieo, issave, k
    integer(kind=KI), dimension(2)        :: iseo
    integer(kind=KI), dimension(4)        :: is, ns
    integer(kind=KI), dimension(2,4)      :: imin, imax
    integer(kind=KI), dimension(2*nvhalf) :: isg
 
!*PART I: Some definitions.
    ldiv(1) = npx/=1
    ldiv(2) = npy/=1
    ldiv(3) = npz/=1
    ldiv(4) = npt/=1
    ns(1) = nx/2/npx
    ns(2) = ny/2/npy
    ns(3) = nz/2/npz
    ns(4) = nt/2/npt
 
!*PART II: make local site index, isg().
    iseo = 0
    do iv = 1,2*nvhalf
     call atoc(iv-1,ns,is)
! ig sums x+y+z+t on one block of the sublattice
     ig = 0
     do mu = 1,4
      ig = ig + ip(mu)*ns(mu) + is(mu)
     enddo ! mu
! ieo=1 if ig is even, otherwise ieo=2.
     ieo = modulo(ig,2) + 1
! iseo(1) counts the number of even sites on one block of the sublattice
! iseo(2) counts the number of odd sites on one block of the sublattice
     iseo(ieo) = iseo(ieo) + 1
     isg(iv) = iseo(ieo)
    enddo ! iv
 
!*PART III: find nearest neighbour in positive direction lv(,,),
!           and boundary index ib(,,,).
    iseo = 0
    imin = 0
    imax = 0
    do iv = 1,2*nvhalf
     call atoc(iv-1,ns,is)
! ig sums x+y+z+t on one block of the sublattice
     ig = 0
     do mu = 1,4
      ig = ig + ip(mu)*ns(mu) + is(mu)
     enddo ! mu
! ieo = 1 if ig is even, otherwise ieo = 2.
     ieo = modulo(ig,2) + 1
! iseo(1) counts the number of even sites on one block of the sublattice
! iseo(2) counts the number of odd sites on one block of the sublattice
     iseo(ieo) = iseo(ieo) + 1
     do mu = 1,4
!-Begin calculation for multiple processes.
      if (ldiv(mu)) then
       if (is(mu)<ns(mu)-1) then
        issave = is(mu)
        is(mu) = is(mu) + 1
        call ctoa(is,ns,k)
! Assign site of nearest neighbour in +mu direction to lv(,,).
! (Note: '+1' in k+1 simply converts from 0,1,2,... convention to isg index.)
        lv(iseo(ieo),mu,ieo) = isg(k+1)
        lvbc(iseo(ieo),mu,ieo) = isg(k+1)
        is(mu) = issave
! At the sublattice boundary, identify nearest neighbour as the buffer portion
! of the link variable u(); i.e. choose lv()>nvhalf.
       elseif (is(mu)==ns(mu)-1) then
        imax(ieo,mu) = imax(ieo,mu) + 1
        lv(iseo(ieo),mu,ieo) = nvhalf + imax(ieo,mu)
        lvbc(iseo(ieo),mu,ieo) = nvhalf + imax(ieo,mu)
        ib(imax(ieo,mu),mu,2,ieo) = iseo(ieo)
       endif
       if (is(mu)==0) then
        imin(ieo,mu) = imin(ieo,mu) + 1
        ib(imin(ieo,mu),mu,1,ieo) = iseo(ieo)
       endif
!-Begin calculation for a single process with non-periodic fermion boundaries.
      elseif (bc(mu)/=1) then
       if (is(mu)<ns(mu)-1) then
        issave = is(mu)
        is(mu) = is(mu) + 1
        call ctoa(is,ns,k)
        lv(iseo(ieo),mu,ieo) = isg(k+1)
        lvbc(iseo(ieo),mu,ieo) = isg(k+1)
        is(mu) = issave
       elseif (is(mu)==ns(mu)-1) then
        issave = is(mu)
        is(mu) = 0
        call ctoa(is,ns,k)
        lv(iseo(ieo),mu,ieo) = isg(k+1)
        is(mu) = issave
        imax(ieo,mu) = imax(ieo,mu) + 1 
        lvbc(iseo(ieo),mu,ieo) = nvhalf + imax(ieo,mu)
        ib(imax(ieo,mu),mu,2,ieo) = iseo(ieo)
       endif
       if (is(mu)==0) then
        imin(ieo,mu) = imin(ieo,mu) + 1
        ib(imin(ieo,mu),mu,1,ieo) = iseo(ieo)
       endif
!-Begin calculation for a single process with periodic fermion boundaries.
      else
       issave = is(mu)
! Identify nearest neighbour in +mu direction.
       is(mu) = modulo(is(mu)+1,ns(mu))
       call ctoa(is,ns,k)
! Assign site of nearest neighbour in +mu direction to lv(,,).
! (Note: k+1 simply converts from 0,1,2,... convention to isg index.)
       lv(iseo(ieo),mu,ieo) = isg(k+1)
       lvbc(iseo(ieo),mu,ieo) = isg(k+1)
       is(mu) = issave
      endif
     enddo ! mu
    enddo ! iv
 
!*PART IV: Determine nms().
    do mu = 1,4
     nms(mu) = nvhalf/ns(mu)
    enddo ! mu
 
 end subroutine siteindex
 
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 
 subroutine setblk(lbd,iblv)
! Build an index for the various blocks of a sublattice.
! OUTPUT:
!   iblv(ibl,mu) is the number of the neighbouring block (to the block
!                numbered ibl) in the mu direction.
!                iblv() runs from 1 through 16.
!   lbd(ibl,mu) is true if ibl is the first of the two blocks in the mu
!               direction, otherwise it is false.
 
    logical,          intent(out), dimension(:,:) :: lbd
    integer(kind=KI), intent(out), dimension(:,:) :: iblv

    integer(kind=KI) :: ib0, ibx, iby, ibz, ibt
 
    ib0 = 0
    do ibt = 1,2
     do ibz = 1,2
      do iby = 1,2
       do ibx = 1,2
        ib0 = ib0 + 1
        lbd(ib0,1) = ibx==1
        lbd(ib0,2) = iby==1
        lbd(ib0,3) = ibz==1
        lbd(ib0,4) = ibt==1
        iblv(ib0,1) = ib0 +    3-2*ibx
        iblv(ib0,2) = ib0 + 2*(3-2*iby)
        iblv(ib0,3) = ib0 + 4*(3-2*ibz)
        iblv(ib0,4) = ib0 + 8*(3-2*ibt)
       enddo ! ibx
      enddo ! iby
     enddo ! ibz
    enddo ! ibt
 
 end subroutine setblk

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 
 subroutine gblchker(vecbl,vecblinv)
! Define the checkerboarding of the entire lattice (distinct from ieo,ibl).
! In particular, the blocks vecbl(1,:), including ieo=1 and ieo=2, form
! a checkerboard pattern.  So do the complement: vecbl(2,:).
! OUTPUT:
!   vecbl() defines the global checkerboarding of the lattice.
!           vecbl(1,:) means sites on blocks ibl=1,4,6,7,10,11,13,16.
!           vecbl(2,:) means sites on blocks ibl=2,3,5,8,9,12,14,15.
!           vecbl(1,:) is called "globally even", because ix+iy+iz+it is even.
!           vecbl(2,:) is called "globally odd", because ix+iy+iz+it is odd.
!   vecblinv() is the inverse of vecbl().
 
    integer(kind=KI), intent(out), dimension(:,:) :: vecbl, vecblinv

    vecbl(1,1:8) = (/ 1, 4, 6, 7, 10, 11, 13, 16 /)
    vecbl(2,1:8) = (/ 2, 3, 5, 8,  9, 12, 14, 15 /)
    vecblinv(1,1:16) = (/ 1, 2, 2, 1, 2, 1, 1, 2, 2, 1, 1, 2, 1, 2, 2, 1 /)
    vecblinv(2,1:16) = (/ 1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7, 8, 8 /)
 
 end subroutine gblchker

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 subroutine fixaction(gaction,beta,alat,tad,coact)
! Define the gauge and fermion actions to be employed.
! gaction(i) = 1 for nearest neighbour gauge action in i'th spacetime direction.
! gaction(i) = 2 for next-nearest neighbour gauge action in the i'th direction.
! In next-nearest neighbour directions, the leading classical discretization
! errors are eliminated from the gauge action.
! The action is fully anisotropic in any number of the 4 spacetime directions.
! Tadpole improvement is built in.
! The tadpole coefficients follow the notation of Aoki et al, hep-lat/0107009.
!
! gauge action = sum_(n,i,j>i) coact(1,i,j)*[1-(1/3)ReTr(Uplaq_(i,j))]
!              + sum_(n,i,j<>i) coact(2,i,j)*[1-(1/3)ReTr(Urect_(i,j))]
! and
! fermion action = sum_(m,n) psibar(m) [A(m,n)-kappa*B(m,n)] psi(n)
! with A(m,n) = delta_(m,n)
!               *[1+kappa*cSW*sum_(i,j>i)coact(4,i,j)*i*sigma_(i,j)*F_(i,j)(m)]
!      B(m,n) = sum_i [coact(3,1,i)*(U_i(m)*delta_(m+i,n)
!                                   +U_i^dagger(n)*delta_(m-i,n))
!                     +coact(3,2,i)*gamma_i*(-U_i(m)*delta_(m+i,n)
!                                           +U_i^dagger(n)*delta_(m-i,n))]
! where n is a lattice site
!       i and j are spacetime directions
!       Uplaq_(i,j) is an elementary plaquette in the ij plane
!       Urect_(i,j) is a rectangle with the long side in the i'th direction
!                   and short side in the j'th direction
!       F_(i,j)(m) is the field strength tensor WITHOUT tadpole improvement,
!                  since the tadpole improvement is within coact().

    integer(kind=KI), intent(in),  dimension(:)     :: gaction
    real(kind=KR),    intent(in)                    :: beta
    real(kind=KR),    intent(in),  dimension(:)     :: alat, tad
    real(kind=KR),    intent(out), dimension(:,:,:) :: coact

    real(kind=KR)                  :: bit
    integer(kind=KI)               :: i, j
    integer(kind=KI), dimension(6) :: a, b

! Initialize.
    coact = 0.0_KR

! For gauge action coefficients, consider the 6 "ab" planes.
    a(1:6) = (/ 1, 1, 1, 2, 2, 3 /)
    b(1:6) = (/ 2, 3, 4, 3, 4, 4 /)
    do i = 1,6
     j = 7 - i
     bit = beta*alat(a(j))*alat(b(j))/alat(a(i))/alat(b(i)) &
           /tad(a(i))**2/tad(b(i))**2
     if (gaction(a(i))==2.and.gaction(b(i))==2) then
      coact(1,a(i),b(i)) = 5.0_KR/3.0_KR*bit
      coact(2,a(i),b(i)) = -bit/12.0_KR/tad(a(i))**2
      coact(2,b(i),a(i)) = -bit/12.0_KR/tad(b(i))**2
     elseif (gaction(a(i))==2.and.gaction(b(i))==1) then
      coact(1,a(i),b(i)) = 4.0_KR/3.0_KR*bit
      coact(2,a(i),b(i)) = -bit/12.0_KR/tad(a(i))**2
      coact(2,b(i),a(i)) = 0.0_KR
     elseif (gaction(a(i))==1.and.gaction(b(i))==2) then
      coact(1,a(i),b(i)) = 4.0_KR/3.0_KR*bit
      coact(2,a(i),b(i)) = 0.0_KR
      coact(2,b(i),a(i)) = -bit/12.0_KR/tad(b(i))**2
     else
      coact(1,a(i),b(i)) = bit
      coact(2,a(i),b(i)) = 0.0_KR
      coact(2,b(i),a(i)) = 0.0_KR
     endif
     coact(1,b(i),a(i)) = coact(1,a(i),b(i))
    enddo ! i

! Fermion action coefficients.
    do i = 1,4
     coact(3,1,i) = 1.0_KR/tad(i)/alat(i)**2
     coact(3,2,i) = 1.0_KR/tad(i)/alat(i)
    enddo ! i
    do i = 1,3
     do j = i+1,4
      coact(4,i,j) = 1.0_KR/(alat(i)*alat(j)*(tad(i)*tad(j))**2)
      coact(4,j,i) = coact(4,i,j)
     enddo ! j
    enddo ! i

 end subroutine fixaction

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 
 end module lattice

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
