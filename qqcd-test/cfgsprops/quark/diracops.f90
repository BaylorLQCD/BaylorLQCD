! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
! diracops.f90, Randy Lewis, randy.lewis@uregina.ca
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
! This module contains subroutines which construct the Dirac operator, M,
! and multiply it onto any vector.
!
! CONVENTION FOR THE (WILSON) DIRAC ACTION:
!
! S = 1/(2 kappa) sum_{x,y} [ psibar(x) M(x,y) psi(y) ]
! where
! M(x,y) = delta_{x,y} - kappa sum_{mu=1,4} [ (1-gam_mu) U_mu(x) delta_{x,y-mu}
!                               + (1+gam_mu) U_mu^dagger(x-mu) delta_{x,y+mu} ]
! The factor 1/(2 kappa) is absorbed by psi(x) -> sqrt(2kappa) psi(x).
!
! CONVENTION FOR CHECKERBOARDING:
!
! The eigenvalue problem is M x = b.
! To make the kappa dependence explicit, one can define M = 1 - kappa * H,
! with
!        H_{x,x+mu} = (1-gam_mu) U_mu(x)
!        H_{x,x-mu} = (1+gam_mu) U_mu^dagger(x-mu),
! which emphasizes the checkerboard pattern of this operator.
! Notice,
!        M/kappa = | 1/kappa  -H_eo   | with x = | x_e |, b = | b_e |
!                  | -H_oe    1/kappa |          | x_o |      | b_o |.
! so that
!        x_o = b_o + kappa*H_oe*x_e;
!        (1/kappa^2 - H_eo*H_oe) * x_e = b_e/kappa^2 + H_eo*b_o/kappa.
! The final solution is thus { x_e, x_o }
! where x_o is obtained directly from x_e, as shown above.
!
! For multimass computations, it should be noted that the equations decouple:
!        (1/kappa^2 - H_eo * H_oe) * y_e = b_e
!        (1/kappa^2 - H_eo * H_oe) * z_e = H_eo * b_o
! The final solution is { x_e = y_e/kappa^2 + z_e/kappa, x_o }
! where x_o is obtained directly from x_e, as shown above.
!
! The checkerboarding of the entire lattice is not simply the index "ieo"
! (4th argument of link fields) due to the presence of "ibl" (5th argument
! of link fields).  The required global checkerboard colour will be
! denoted by "gblclr" such that
!        gblclr=1 means sites on blocks ibl=1,4,6,7,10,11,13,16 (ieo=1 and 2).
!        gblclr=2 means sites on blocks ibl=2,3,5,8, 9,12,14,15 (ieo=1 and 2).
! Thus, a complete Dirac vector follows the index structure of a link variable:
!        g(a,b,c,d,e), where
! a=1,6 runs over Re and Im parts of the colour vector as follows:
!            / v(1)+i*v(2) \
!        v = | v(3)+i*v(4) |, where i=sqrt(-1)
!            \ v(5)+i*v(6) /
! b=1,ntotal runs over even or odd sites on one block of a sublattice, plus
!            room for the largest boundary shared by two processes.
! c=1,4 runs over the 4 Dirac components of the spinor.
! d=1,2 runs over the even(1) and odd(2) lattice sites (ieo, NOT gblclr).
! e=1,16 runs over the blocks of a sublattice.
! while the notation for a checkerboarded Dirac vector is
!        g(a,b,c,d,ibleo) where a,b,c,d are as above.
!
! CONVENTION FOR DIRAC GAMMA MATRICES:
!
!      / 0 0 0 1 \       / 0  0 0 -i \       / 0  0 1  0 \       / 1 0  0  0 \
! gam1=| 0 0 1 0 |, gam2=| 0  0 i  0 |, gam3=| 0  0 0 -1 |, gam4=| 0 1  0  0 |.
!      | 0 1 0 0 |       | 0 -i 0  0 |       | 1  0 0  0 |       | 0 0 -1  0 |
!      \ 1 0 0 0 /       \ i  0 0  0 /       \ 0 -1 0  0 /       \ 0 0  0 -1 /
!
!                                        / 0 0 -i  0 \
! Note that gam5 = gam1*gam2*gam3*gam4 = | 0 0  0 -i |.
!                                        | i 0  0  0 |
!                                        \ 0 i  0  0 /
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 module diracops

 !   use MPI
    use kinds
    use latdims
    use basics
    use lattice
    implicit none
    private

! Use the following line if the MPI module is not available.
   include 'mpif.h'

! Define access to subroutines.
    public  :: Hdouble, Hdouble1, Hdbletm, Hsingle, multG, &
               smear, setmvf, mv, mdv, &
               setvector, gamma5mult, twMshift, zsmear, gammamult, gammamultshort
    private :: mvg, mdvg, slidevector,mulbac, mulfor

 contains

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 ! Multiply a Randy Lewis vector by any gamma matrix. The Gamma matrix 
 ! convention is that which is defined above. 
 ! INPUT: 
 !    evenIN/oddIN : input vectors (not changed in code)
 !    mu : 1=gamma1, 2=gamma2, 3=gamma3, 4=gamma4, 5=gamma5
 subroutine gammamult(evenIN, oddIN, evenOUT, oddOUT, mu,ieo,ibleo)
 real(kind=KR2), dimension(:,:,:,:,:), intent(in)  :: evenIN, oddIN
 real(kind=KR2), dimension(:,:,:,:,:), intent(inout) :: evenOUT, oddOUT
 integer(kind=KI), intent(in)                      :: mu,ieo,ibleo

 integer(kind=KI)                                  :: j

 if (mu == 5) then  ! gamma5
   do j=1,5,2
     evenOUT(j  ,:nvhalf,3,ieo,ibleo) = -evenIN(j+1,:nvhalf,1,ieo,ibleo)
     evenOUT(j+1,:nvhalf,3,ieo,ibleo) =  evenIN(j  ,:nvhalf,1,ieo,ibleo)
     oddOUT(j  ,:nvhalf,3,ieo,ibleo) = -oddIN(j+1,:nvhalf,1,ieo,ibleo)
     oddOUT(j+1,:nvhalf,3,ieo,ibleo) =  oddIN(j  ,:nvhalf,1,ieo,ibleo)

     evenOUT(j  ,:nvhalf,4,ieo,ibleo) = -evenIN(j+1,:nvhalf,2,ieo,ibleo)
     evenOUT(j+1,:nvhalf,4,ieo,ibleo) =  evenIN(j  ,:nvhalf,2,ieo,ibleo)
     oddOUT(j  ,:nvhalf,4,ieo,ibleo) = -oddIN(j+1,:nvhalf,2,ieo,ibleo)
     oddOUT(j+1,:nvhalf,4,ieo,ibleo) =  oddIN(j  ,:nvhalf,2,ieo,ibleo)

     evenOUT(j  ,:nvhalf,1,ieo,ibleo) =  evenIN(j+1,:nvhalf,3,ieo,ibleo)
     evenOUT(j+1,:nvhalf,1,ieo,ibleo) = -evenIN(j  ,:nvhalf,3,ieo,ibleo)
     oddOUT(j  ,:nvhalf,1,ieo,ibleo) =  oddIN(j+1,:nvhalf,3,ieo,ibleo)
     oddOUT(j+1,:nvhalf,1,ieo,ibleo) = -oddIN(j  ,:nvhalf,3,ieo,ibleo)

     evenOUT(j  ,:nvhalf,2,ieo,ibleo) =  evenIN(j+1,:nvhalf,4,ieo,ibleo)
     evenOUT(j+1,:nvhalf,2,ieo,ibleo) = -evenIN(j  ,:nvhalf,4,ieo,ibleo)
     oddOUT(j  ,:nvhalf,2,ieo,ibleo) =  oddIN(j+1,:nvhalf,4,ieo,ibleo)
     oddOUT(j+1,:nvhalf,2,ieo,ibleo) = -oddIN(j  ,:nvhalf,4,ieo,ibleo)
   enddo
 elseif (mu == 1) then  ! gamma1
   evenOUT(:6,:nvhalf,4,ieo,ibleo) = evenIN(:6,:nvhalf,1,ieo,ibleo)
   oddOUT(:6,:nvhalf,4,ieo,ibleo) = oddIN(:6,:nvhalf,1,ieo,ibleo)

   evenOUT(:6,:nvhalf,3,ieo,ibleo) = evenIN(:6,:nvhalf,2,ieo,ibleo)
   oddOUT(:6,:nvhalf,3,ieo,ibleo) = oddIN(:6,:nvhalf,2,ieo,ibleo)

   evenOUT(:6,:nvhalf,2,ieo,ibleo) = evenIN(:6,:nvhalf,3,ieo,ibleo)
   oddOUT(:6,:nvhalf,2,ieo,ibleo) = oddIN(:6,:nvhalf,3,ieo,ibleo)

   evenOUT(:6,:nvhalf,1,ieo,ibleo) = evenIN(:6,:nvhalf,4,ieo,ibleo)
   oddOUT(:6,:nvhalf,1,ieo,ibleo) = oddIN(:6,:nvhalf,4,ieo,ibleo)
 elseif (mu == 2) then  ! gamma2
   do j=1,5,2
     evenOUT(j  ,:nvhalf,4,ieo,ibleo) = -evenIN(j+1,:nvhalf,1,ieo,ibleo)
     evenOUT(j+1,:nvhalf,4,ieo,ibleo) =  evenIN(j  ,:nvhalf,1,ieo,ibleo)
     oddOUT(j  ,:nvhalf,4,ieo,ibleo) = -oddIN(j+1,:nvhalf,1,ieo,ibleo)
     oddOUT(j+1,:nvhalf,4,ieo,ibleo) =  oddIN(j  ,:nvhalf,1,ieo,ibleo)

     evenOUT(j  ,:nvhalf,3,ieo,ibleo) =  evenIN(j+1,:nvhalf,2,ieo,ibleo)
     evenOUT(j+1,:nvhalf,3,ieo,ibleo) = -evenIN(j  ,:nvhalf,2,ieo,ibleo)
     oddOUT(j  ,:nvhalf,3,ieo,ibleo) =  oddIN(j+1,:nvhalf,2,ieo,ibleo)
     oddOUT(j+1,:nvhalf,3,ieo,ibleo) = -oddIN(j  ,:nvhalf,2,ieo,ibleo)

     evenOUT(j  ,:nvhalf,2,ieo,ibleo) = -evenIN(j+1,:nvhalf,3,ieo,ibleo)
     evenOUT(j+1,:nvhalf,2,ieo,ibleo) =  evenIN(j  ,:nvhalf,3,ieo,ibleo)
     oddOUT(j  ,:nvhalf,2,ieo,ibleo) = -oddIN(j+1,:nvhalf,3,ieo,ibleo)
     oddOUT(j+1,:nvhalf,2,ieo,ibleo) =  oddIN(j  ,:nvhalf,3,ieo,ibleo)

     evenOUT(j  ,:nvhalf,1,ieo,ibleo) =  evenIN(j+1,:nvhalf,4,ieo,ibleo)
     evenOUT(j+1,:nvhalf,1,ieo,ibleo) = -evenIN(j  ,:nvhalf,4,ieo,ibleo)
     oddOUT(j  ,:nvhalf,1,ieo,ibleo) =  oddIN(j+1,:nvhalf,4,ieo,ibleo)
     oddOUT(j+1,:nvhalf,1,ieo,ibleo) = -oddIN(j  ,:nvhalf,4,ieo,ibleo)
   enddo
 elseif (mu == 3) then  ! gamma3
   evenOUT(:6,:nvhalf,3,ieo,ibleo) = evenIN(:6,:nvhalf,1,ieo,ibleo)
   oddOUT(:6,:nvhalf,3,ieo,ibleo) = oddIN(:6,:nvhalf,1,ieo,ibleo)

   evenOUT(:6,:nvhalf,4,ieo,ibleo) = -evenIN(:6,:nvhalf,2,ieo,ibleo)
   oddOUT(:6,:nvhalf,4,ieo,ibleo) = -oddIN(:6,:nvhalf,2,ieo,ibleo)

   evenOUT(:6,:nvhalf,1,ieo,ibleo) = evenIN(:6,:nvhalf,3,ieo,ibleo)
   oddOUT(:6,:nvhalf,1,ieo,ibleo) = oddIN(:6,:nvhalf,3,ieo,ibleo)

   evenOUT(:6,:nvhalf,2,ieo,ibleo) = -evenIN(:6,:nvhalf,4,ieo,ibleo)
   oddOUT(:6,:nvhalf,2,ieo,ibleo) = -oddIN(:6,:nvhalf,4,ieo,ibleo)
 elseif (mu == 4) then  ! gamma4
   evenOUT(:6,:nvhalf,1,ieo,ibleo) = evenIN(:6,:nvhalf,1,ieo,ibleo)
   oddOUT(:6,:nvhalf,1,ieo,ibleo) = oddIN(:6,:nvhalf,1,ieo,ibleo)

   evenOUT(:6,:nvhalf,2,ieo,ibleo) = evenIN(:6,:nvhalf,2,ieo,ibleo)
   oddOUT(:6,:nvhalf,2,ieo,ibleo) = oddIN(:6,:nvhalf,2,ieo,ibleo)

   evenOUT(:6,:nvhalf,3,ieo,ibleo) = -evenIN(:6,:nvhalf,3,ieo,ibleo)
   oddOUT(:6,:nvhalf,3,ieo,ibleo) = -oddIN(:6,:nvhalf,3,ieo,ibleo)

   evenOUT(:6,:nvhalf,4,ieo,ibleo) = -evenIN(:6,:nvhalf,4,ieo,ibleo)
   oddOUT(:6,:nvhalf,4,ieo,ibleo) = -oddIN(:6,:nvhalf,4,ieo,ibleo)
 else
   write(*,*) 'ERROR: mu value is not valid'
 endif
 end subroutine gammamult



 ! Multiply a Randy Lewis vector by any gamma matrix. The Gamma matrix 
 ! convention is that which is defined above. 
 ! INPUT: 
 !    evenIN/oddIN : input vectors (not changed in code)
 !    mu : 1=gamma1, 2=gamma2, 3=gamma3, 4=gamma4, 5=gamma5
 subroutine gammamultshort(vIN,vOUT, mu)
 real(kind=KR2), dimension(:,:,:), intent(in)    :: vIN
 real(kind=KR2), dimension(:,:,:), intent(inout) :: vOUT 
 integer(kind=KI), intent(in)                    :: mu

 integer(kind=KI)                                :: j

 if (mu == 5) then  ! gamma5
   do j=1,5,2
     vOUT(j  ,:nvhalf,3) = -vIN(j+1,:nvhalf,1)
     vOUT(j+1,:nvhalf,3) =  vIN(j  ,:nvhalf,1)

     vOUT(j  ,:nvhalf,4) = -vIN(j+1,:nvhalf,2)
     vOUT(j+1,:nvhalf,4) =  vIN(j  ,:nvhalf,2)

     vOUT(j  ,:nvhalf,1) =  vIN(j+1,:nvhalf,3)
     vOUT(j+1,:nvhalf,1) = -vIN(j  ,:nvhalf,3)

     vOUT(j  ,:nvhalf,2) =  vIN(j+1,:nvhalf,4)
     vOUT(j+1,:nvhalf,2) = -vIN(j  ,:nvhalf,4)
   enddo
 elseif (mu == 1) then  ! gamma1
   vOUT(:6,:nvhalf,4) = vIN(:6,:nvhalf,1)
   vOUT(:6,:nvhalf,3) = vIN(:6,:nvhalf,2)
   vOUT(:6,:nvhalf,2) = vIN(:6,:nvhalf,3)
   vOUT(:6,:nvhalf,1) = vIN(:6,:nvhalf,4)
 elseif (mu == 2) then  ! gamma2
   do j=1,5,2
     vOUT(j  ,:nvhalf,4) = -vIN(j+1,:nvhalf,1)
     vOUT(j+1,:nvhalf,4) =  vIN(j  ,:nvhalf,1)

     vOUT(j  ,:nvhalf,3) =  vIN(j+1,:nvhalf,2)
     vOUT(j+1,:nvhalf,3) = -vIN(j  ,:nvhalf,2)

     vOUT(j  ,:nvhalf,2) = -vIN(j+1,:nvhalf,3)
     vOUT(j+1,:nvhalf,2) =  vIN(j  ,:nvhalf,3)

     vOUT(j  ,:nvhalf,1) =  vIN(j+1,:nvhalf,4)
     vOUT(j+1,:nvhalf,1) = -vIN(j  ,:nvhalf,4)
   enddo
 elseif (mu == 3) then  ! gamma3
   vOUT(:6,:nvhalf,3) =  vIN(:6,:nvhalf,1)
   vOUT(:6,:nvhalf,4) = -vIN(:6,:nvhalf,2)
   vOUT(:6,:nvhalf,1) =  vIN(:6,:nvhalf,3)
   vOUT(:6,:nvhalf,2) = -vIN(:6,:nvhalf,4)
 elseif (mu == 4) then  ! gamma4
   vOUT(:6,:nvhalf,1) =  vIN(:6,:nvhalf,1)
   vOUT(:6,:nvhalf,2) =  vIN(:6,:nvhalf,2)
   vOUT(:6,:nvhalf,3) = -vIN(:6,:nvhalf,3)
   vOUT(:6,:nvhalf,4) = -vIN(:6,:nvhalf,4)
 else
   write(*,*) 'ERROR: mu value is not valid'
 endif
 end subroutine gammamultshort



!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

 subroutine Hdouble(he,u,ge,idag,coact,kappa,bc,vecbl,vecblinv,myid,nn,ldiv, &
                    nms,lvbc,ib,lbd,iblv,MRT)
! Multiply the operator (1/kappa^2-H_eo*H_oe) onto a given vector, ge().
! To get the Lagrangian term, one would premultiply this addition by gebar(x).
! For use in CGNE, idag=1 allows use of the operator's dagger instead.
! INPUT:
!   u() contains the gauge fields for this sublattice.
!   ge() contains the gblclr=1 half of the Dirac spinor for this sublattice.
!   idag=0 is for computing the operator; idag=1 is for computing H^dagger.
!   coact(2,4,4) contains the coefficients from the action.
!   kappa is the hopping parameter.
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
! OUTPUT:
!   he() is set to (1/kappa^2-H_eo*H_oe) * ge().
!   expected size: he(a,b,c),
!                  where a=1..6 is the real/imaginary and colour index.
!                        b=1,2,3,...,nvhalf combines colour and position.
!                        c=1,2,3,4 is the Dirac spinor index.

    real(kind=KR),    intent(out),   dimension(:,:,:,:,:) :: he
    real(kind=KR),    intent(in),    dimension(:,:,:,:,:) :: u
    real(kind=KR),    intent(inout), dimension(:,:,:,:,:) :: ge
    integer(kind=KI), intent(in)                          :: idag, myid, MRT
    real(kind=KR),    intent(in),    dimension(:,:,:)     :: coact
    real(kind=KR),    intent(in)                          :: kappa
    integer(kind=KI), intent(in),    dimension(:)         :: bc, nms
    integer(kind=KI), intent(in),    dimension(:,:)       :: vecbl, vecblinv
    integer(kind=KI), intent(in),    dimension(:,:)       :: nn, iblv
    logical,          intent(in),    dimension(:)         :: ldiv
    integer(kind=KI), intent(in),    dimension(:,:,:)     :: lvbc
    integer(kind=KI), intent(in),    dimension(:,:,:,:)   :: ib
    logical,          intent(in),    dimension(:,:)       :: lbd

    real(kind=KR), dimension(6,ntotal,4,2,8) :: htmp
    integer(kind=KI)                         :: gblclr, isite
    real(kind=KR)                            :: fac1

! First build H_oe * g_e.
    gblclr = 2
    call Hsingle(htmp,u,ge,idag,coact,bc,gblclr,vecbl,vecblinv,myid,nn,ldiv, &
                 nms,lvbc,ib,lbd,iblv,MRT)

! Next build H_eo * H_oe * g_e.
    gblclr = 1
    call Hsingle(he,u,htmp,idag,coact,bc,gblclr,vecbl,vecblinv,myid,nn,ldiv, &
                 nms,lvbc,ib,lbd,iblv,MRT)

! Finally, construct h_e = (1/kappa^2-H_eo*H_oe)*g_e.
    fac1 = 1.0_KR/kappa**2
    do isite = 1,nvhalf
     he(:,isite,:,:,:) = fac1*ge(:,isite,:,:,:) - he(:,isite,:,:,:)
    enddo ! isite

 end subroutine Hdouble

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
subroutine Hdouble1(ho,u,go,idag,coact,kappa,bc,vecbl,vecblinv,myid,nn,ldiv, &
                    nms,lvbc,ib,lbd,iblv,MRT)
! Multiply the operator (1/kappa^2-H_oe*H_eo) onto a given vector, go().
! To get the Lagrangian term, one would premultiply this addition by gebar(x).
! For use in CGNE, idag=1 allows use of the operator's dagger instead.
! INPUT:
!   u() contains the gauge fields for this sublattice.
!   ge() contains the gblclr=1 half of the Dirac spinor for this sublattice.
!   idag=0 is for computing the operator; idag=1 is for computing H^dagger.
!   coact(2,4,4) contains the coefficients from the action.
!   kappa is the hopping parameter.
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
! OUTPUT:
!   he() is set to (1/kappa^2-H_eo*H_oe) * ge().
!   expected size: he(a,b,c),
!                  where a=1..6 is the real/imaginary and colour index.
!                        b=1,2,3,...,nvhalf combines colour and position.
!                        c=1,2,3,4 is the Dirac spinor index.

    real(kind=KR),    intent(out),   dimension(:,:,:,:,:) :: ho
    real(kind=KR),    intent(in),    dimension(:,:,:,:,:) :: u
    real(kind=KR),    intent(inout), dimension(:,:,:,:,:) :: go
    integer(kind=KI), intent(in)                          :: idag, myid, MRT
    real(kind=KR),    intent(in),    dimension(:,:,:)     :: coact
    real(kind=KR),    intent(in)                          :: kappa
    integer(kind=KI), intent(in),    dimension(:)         :: bc, nms
    integer(kind=KI), intent(in),    dimension(:,:)       :: vecbl, vecblinv
    integer(kind=KI), intent(in),    dimension(:,:)       :: nn, iblv
    logical,          intent(in),    dimension(:)         :: ldiv
    integer(kind=KI), intent(in),    dimension(:,:,:)     :: lvbc
    integer(kind=KI), intent(in),    dimension(:,:,:,:)   :: ib
    logical,          intent(in),    dimension(:,:)       :: lbd

    real(kind=KR), dimension(6,ntotal,4,2,8) :: htmp
    integer(kind=KI)                         :: gblclr, isite
    real(kind=KR)                            :: fac1

! First build H_eo * g_o.
    gblclr = 1
    call Hsingle(htmp,u,go,idag,coact,bc,gblclr,vecbl,vecblinv,myid,nn,ldiv, &
                 nms,lvbc,ib,lbd,iblv,MRT)

! Next build H_oe * H_eo * g_o.
    gblclr = 2
    call Hsingle(ho,u,htmp,idag,coact,bc,gblclr,vecbl,vecblinv,myid,nn,ldiv, &
                 nms,lvbc,ib,lbd,iblv,MRT)

! Finally, construct h_o = (1/kappa^2-H_oe*H_eo)*g_o.
    fac1 = 1.0_KR/kappa**2
    do isite = 1,nvhalf
     ho(:,isite,:,:,:) = fac1*go(:,isite,:,:,:) - ho(:,isite,:,:,:)
    enddo ! isite

 end subroutine Hdouble1

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 subroutine Hdbletm(he,u,GeeGooinv,ge,idag,coact,kappa,iflag,bc,vecbl, &
                    vecblinv,myid,nn,ldiv,nms,lvbc,ib,lbd,iblv,MRT)
! Generalization of subroutine Hdouble to the case of twisted mass QCD.
! Multiply the operator G_ee - H_eo*G_oo^(-1)*H_oe onto a given vector, ge().
! where G = 1 - kappa*C + 2*kappa*mu*i*gam5.  [NOTE: C is the clover term.]
!   kappa(1) is 1/(8+2*m_0) with m_0 from eq.(1.1) of JHEP08(2001)058.
!   kappa(2) is mu_q from eq.(1.1) of JHEP08(2001)058.
! For use in CGNE, idag=1 allows use of the operator's dagger instead.
! INPUT: iflag=-1 for Wilson or -2 for clover. 


    real(kind=KR),    intent(out),   dimension(:,:,:,:,:) :: he
    real(kind=KR),    intent(in),    dimension(:,:,:,:,:) :: u
    real(kind=KR),    intent(in),    dimension(:,:,:,:,:) :: GeeGooinv
    real(kind=KR),    intent(inout), dimension(:,:,:,:,:) :: ge
    integer(kind=KI), intent(in)                      :: idag, iflag, myid, MRT
    real(kind=KR),    intent(in),    dimension(:,:,:)     :: coact
    real(kind=KR),    intent(in),    dimension(:)         :: kappa
    integer(kind=KI), intent(in),    dimension(:)         :: bc, nms
    integer(kind=KI), intent(in),    dimension(:,:)       :: vecbl, vecblinv
    integer(kind=KI), intent(in),    dimension(:,:)       :: nn, iblv
    logical,          intent(in),    dimension(:)         :: ldiv
    integer(kind=KI), intent(in),    dimension(:,:,:)     :: lvbc
    integer(kind=KI), intent(in),    dimension(:,:,:,:)   :: ib
    logical,          intent(in),    dimension(:,:)       :: lbd

    real(kind=KR), dimension(6,ntotal,4,2,8) :: htmp
    integer(kind=KI)                         :: gblclr, isite,i,j
    real(kind=KR)                            :: fac1, fac2, fac3

     !for kappa(2)=0 and no clover term use Hdouble to avoid roundoff errors
     !coming from un-necessary multiplications
     if( (kappa(2)<0.00000001) .and. (iflag==-1)) then

      !if(myid==0) then
      ! open(unit=8, file="/home/darnelld/qqcd/scratch/CFGSPROPS.LOG",action="write",&
      !      status="old",position="append")
      ! write(unit=8,fmt="(a30,2d20.10)")"using Hdouble, kappa=",kappa(1),kappa(2)
      ! close(unit=8,status="keep")
      !endif !myid==0

      call  Hdouble(he,u,ge,idag,coact,kappa(1),bc,vecbl,vecblinv,myid,nn,ldiv, &
                    nms,lvbc,ib,lbd,iblv,MRT)
     else



!*WILSON.
    if (iflag==-1) then
! First build H_oe * g_e.
     gblclr = 2
     call Hsingle(he,u,ge,idag,coact,bc,gblclr,vecbl,vecblinv,myid,nn,ldiv, &
                  nms,lvbc,ib,lbd,iblv,MRT)
! Next build (1-2*kappa*mu*i*gam5) * H_oe * g_e.
     fac1 = 2.0_KR*kappa(1)*kappa(2)
     do isite = 1,nvhalf
      htmp(:,isite,1,:,:) = he(:,isite,1,:,:) - fac1*he(:,isite,3,:,:)
      htmp(:,isite,2,:,:) = he(:,isite,2,:,:) - fac1*he(:,isite,4,:,:)
      htmp(:,isite,3,:,:) = he(:,isite,3,:,:) + fac1*he(:,isite,1,:,:)
      htmp(:,isite,4,:,:) = he(:,isite,4,:,:) + fac1*he(:,isite,2,:,:)
     enddo ! isite
! Next build H_eo * (1-2*kappa*mu*i*gam5) * H_oe * g_e.
     gblclr = 1
     call Hsingle(he,u,htmp,idag,coact,bc,gblclr,vecbl,vecblinv,myid,nn,ldiv, &
                  nms,lvbc,ib,lbd,iblv,MRT)
! Finally, construct h_e = [(1+2*kappa*mu*i*gam5)/kappa^2
!               - H_eo*(1-2*kappa*mu*i*gam5)*H_oe/(1+4*kappa^2*mu^2)] * g_e.
     fac1 = 1.0_KR/kappa(1)**2
     fac2 = 2.0_KR*kappa(2)/kappa(1)
     fac3 = 1.0_KR/( 1.0_KR + (2.0_KR*kappa(1)*kappa(2))**2 )
     do isite = 1,nvhalf
      he(:,isite,1,:,:) = fac1*ge(:,isite,1,:,:) + fac2*ge(:,isite,3,:,:) &
                        - fac3*he(:,isite,1,:,:)
      he(:,isite,2,:,:) = fac1*ge(:,isite,2,:,:) + fac2*ge(:,isite,4,:,:) &
                        - fac3*he(:,isite,2,:,:)
      he(:,isite,3,:,:) = fac1*ge(:,isite,3,:,:) - fac2*ge(:,isite,1,:,:) &
                        - fac3*he(:,isite,3,:,:)
      he(:,isite,4,:,:) = fac1*ge(:,isite,4,:,:) - fac2*ge(:,isite,2,:,:) &
                        - fac3*he(:,isite,4,:,:)
     enddo ! isite
!    if (myid==0) then
!    do j =1,35
!    do i =1,4
!    print *, "he", he(1,j,i,1,1), "j=", j, "i=", i
!    enddo ! i
!    enddo ! j
!    endif
!*CLOVER.
    elseif (iflag==-2) then
! First build H_oe * g_e.
     gblclr = 2
     call Hsingle(htmp,u,ge,idag,coact,bc,gblclr,vecbl,vecblinv,myid,nn,ldiv, &
                  nms,lvbc,ib,lbd,iblv,MRT)
! Next build G_oo^(-1) * H_oe * g_e.
     gblclr = 2
     call multG(GeeGooinv,htmp,gblclr,vecbl)
! Next build H_eo * G_oo^(-1) * H_oe * g_e.
     gblclr = 1
     call Hsingle(he,u,htmp,idag,coact,bc,gblclr,vecbl,vecblinv,myid,nn,ldiv, &
                  nms,lvbc,ib,lbd,iblv,MRT)
! Finally, construct h_e = [G_ee/kappa^2 - H_eo*G_oo^(-1)*H_oe] * g_e.
     do isite = 1,nvhalf
      htmp(:,isite,:,:,:) = ge(:,isite,:,:,:)
     enddo ! isite
     gblclr = 1
     call multG(GeeGooinv,htmp,gblclr,vecbl)
     fac1 = 1.0_KR/kappa(1)**2
     do isite = 1,nvhalf
      he(:,isite,:,:,:) = fac1*htmp(:,isite,:,:,:) - he(:,isite,:,:,:)
     enddo ! isite
    endif


    endif !if(kappa(2)<0.000001.and.....)

 end subroutine Hdbletm

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 subroutine Hsingle(gout,u,gin,idag,coact,bc,gblclr,vecbl,vecblinv,myid,nn, &
                    ldiv,nms,lvbc,ib,lbd,iblv,MRT)
! Perform the multiplication of H_eo or H_oe onto a given vector.
! INPUT:
!   u() contains the gauge fields for this sublattice.
!   gin() contains the Dirac spinor for this sublattice.
!   idag=0 is for computing H; idag=1 is for computing H^dagger.
!   coact(2,4,4) contains the coefficients from the action.
!   bc(mu) = 1,-1,0 for periodic,antiperiodic,fixed boundary conditions
!            respectively.
!   vecbl() defines the global checkerboarding of the lattice.
!           vecbl(1,:) means sites on blocks ibl=1,4,6,7,10,11,13,16.
!           vecbl(2,:) means sites on blocks ibl=2,3,5,8,9,12,14,15.
! OUTPUT:
!   gout() is set to H_eo * gin() if gblclr=1, or to H_oe * gin() if gblclr=2.
!   expected size: gout(a,b,c),
!                  where a=1..6 is the real/imaginary and colour index.
!                        b=1,2,3,...,nvhalf combines colour and position.
!                        c=1,2,3,4 is the Dirac spinor index.
    real(kind=KR),    intent(out),   dimension(:,:,:,:,:) :: gout
    real(kind=KR),    intent(in),    dimension(:,:,:,:,:) :: u
    real(kind=KR),    intent(inout), dimension(:,:,:,:,:) :: gin
    integer(kind=KI), intent(in),    dimension(:)         :: bc, nms
    integer(kind=KI), intent(in)                     :: idag, gblclr, myid, MRT
    real(kind=KR),    intent(in),    dimension(:,:,:)     :: coact
    integer(kind=KI), intent(in),    dimension(:,:)       :: vecbl, vecblinv
    integer(kind=KI), intent(in),    dimension(:,:)       :: nn, iblv
    logical,          intent(in),    dimension(:)         :: ldiv
    integer(kind=KI), intent(in),    dimension(:,:,:)     :: lvbc
    integer(kind=KI), intent(in),    dimension(:,:,:,:)   :: ib
    logical,          intent(in),    dimension(:,:)       :: lbd

    integer(kind=KI) :: mu

! Some initializations.
    gout = 0.0_KR

! First build H_eo * g_o (gblclr=1) or H_oe * g_e (gblclr=2).
    do mu = 1,4
     call mulbac(gout,u,gin,mu,idag,coact,bc,gblclr,vecbl,vecblinv,myid,nn, &
                 ldiv,nms,lvbc,ib,lbd,iblv,MRT)
     call mulfor(gout,u,gin,mu,idag,coact,bc,gblclr,vecbl,vecblinv,myid,nn, &
                 ldiv,nms,lvbc,ib,lbd,iblv,MRT)
    enddo ! mu

 end subroutine Hsingle

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 subroutine multG(GeeGooinv,spinor,gblclr,vecbl)
! Multiply the matrix G_ee or G_oo^(-1) by a spinor:
!                 G*S = /  V  iW \ / S_u \ = / V*S_u + iW*S_d \
!                       \ -iW  V / \ S_d /   \ V*S_d - iW*S_u /
! INPUT:
!   GeeGooinv() contains the matrix G_ee on globally-even sites and
!               the matrix G_oo^(-1) on globally-odd sites.
!   spinor() contains the Dirac spinor for globally-even or globally-odd sites.
!   gblclr is the "global checkerboard colour" of spinor().
!          NOTE: This is NOT ieo, because of the lattice's blocked structure.
!                gblclr=1 means sites on blocks 1,4,6,7,10,11,13,16.
!                gblclr=2 means sites on blocks 2,3,5,8,9,12,14,15.
!   vecbl() defines the checkerboarding of the entire lattice.
!           In particular, the blocks vecbl(1,:), including ieo=1 and ieo=2,
!           form a checkerboard pattern.  So do the complement: vecbl(2,:).
! OUTPUT:
!   spinor = G_ee*spinor for gblclr=1, or G_oo^(-1)*spinor for gblclr=2.

    real(kind=KR),    intent(inout), dimension(:,:,:,:,:) :: spinor
    real(kind=KR),    intent(in),    dimension(:,:,:,:,:) :: GeeGooinv
    integer(kind=KI), intent(in)                          :: gblclr
    integer(kind=KI), intent(in),    dimension(:,:)       :: vecbl

    integer(kind=KI)                     :: icount, ieo, ibl, ibleo, isite, ir
    real(kind=KR), dimension(6,nvhalf)   :: prod
    real(kind=KR), dimension(6,nvhalf,4) :: localspinor

    icount = nvhalf

! Multiply at each site having Dirac global checkerboard colour = gblclr.
    do ibleo = 1,8
     ibl = vecbl(gblclr,ibleo)
     do ieo = 1,2
      localspinor = 0.0_KR

! Compute V*S_u and add to upper components of result.
      call mv(icount,GeeGooinv(:,:,1,ieo,ibl),spinor(:,:,1,ieo,ibleo),prod)
      localspinor(:,:,1) = localspinor(:,:,1) + prod(:,:)
      call mv(icount,GeeGooinv(:,:,2,ieo,ibl),spinor(:,:,2,ieo,ibleo),prod)
      localspinor(:,:,1) = localspinor(:,:,1) + prod(:,:)
      call mv(icount,GeeGooinv(:,:,3,ieo,ibl),spinor(:,:,1,ieo,ibleo),prod)
      localspinor(:,:,2) = localspinor(:,:,2) + prod(:,:)
      call mv(icount,GeeGooinv(:,:,4,ieo,ibl),spinor(:,:,2,ieo,ibleo),prod)
      localspinor(:,:,2) = localspinor(:,:,2) + prod(:,:)

! Compute V*S_d and add to lower components of result.
      call mv(icount,GeeGooinv(:,:,1,ieo,ibl),spinor(:,:,3,ieo,ibleo),prod)
      localspinor(:,:,3) = localspinor(:,:,3) + prod(:,:)
      call mv(icount,GeeGooinv(:,:,2,ieo,ibl),spinor(:,:,4,ieo,ibleo),prod)
      localspinor(:,:,3) = localspinor(:,:,3) + prod(:,:)
      call mv(icount,GeeGooinv(:,:,3,ieo,ibl),spinor(:,:,3,ieo,ibleo),prod)
      localspinor(:,:,4) = localspinor(:,:,4) + prod(:,:)
      call mv(icount,GeeGooinv(:,:,4,ieo,ibl),spinor(:,:,4,ieo,ibleo),prod)
      localspinor(:,:,4) = localspinor(:,:,4) + prod(:,:)

! Compute i*W*S_u and subtract from lower components of result.
      call mv(icount,GeeGooinv(:,:,5,ieo,ibl),spinor(:,:,1,ieo,ibleo),prod)
      do ir = 1,5,2
       localspinor(ir  ,:,3) = localspinor(ir  ,:,3) + prod(ir+1,:)
       localspinor(ir+1,:,3) = localspinor(ir+1,:,3) - prod(ir  ,:)
      enddo ! ir
      call mv(icount,GeeGooinv(:,:,6,ieo,ibl),spinor(:,:,2,ieo,ibleo),prod)
      do ir = 1,5,2
       localspinor(ir  ,:,3) = localspinor(ir  ,:,3) + prod(ir+1,:)
       localspinor(ir+1,:,3) = localspinor(ir+1,:,3) - prod(ir  ,:)
      enddo ! ir
      call mv(icount,GeeGooinv(:,:,7,ieo,ibl),spinor(:,:,1,ieo,ibleo),prod)
      do ir = 1,5,2
       localspinor(ir  ,:,4) = localspinor(ir  ,:,4) + prod(ir+1,:)
       localspinor(ir+1,:,4) = localspinor(ir+1,:,4) - prod(ir  ,:)
      enddo ! ir
      call mv(icount,GeeGooinv(:,:,8,ieo,ibl),spinor(:,:,2,ieo,ibleo),prod)
      do ir = 1,5,2
       localspinor(ir  ,:,4) = localspinor(ir  ,:,4) + prod(ir+1,:)
       localspinor(ir+1,:,4) = localspinor(ir+1,:,4) - prod(ir  ,:)
      enddo ! ir
      
! Compute i*W*S_d and add to upper components of result.
      call mv(icount,GeeGooinv(:,:,5,ieo,ibl),spinor(:,:,3,ieo,ibleo),prod)
      do ir = 1,5,2
       localspinor(ir  ,:,1) = localspinor(ir  ,:,1) - prod(ir+1,:)
       localspinor(ir+1,:,1) = localspinor(ir+1,:,1) + prod(ir  ,:)
      enddo ! ir
      call mv(icount,GeeGooinv(:,:,6,ieo,ibl),spinor(:,:,4,ieo,ibleo),prod)
      do ir = 1,5,2
       localspinor(ir  ,:,1) = localspinor(ir  ,:,1) - prod(ir+1,:)
       localspinor(ir+1,:,1) = localspinor(ir+1,:,1) + prod(ir  ,:)
      enddo ! ir
      call mv(icount,GeeGooinv(:,:,7,ieo,ibl),spinor(:,:,3,ieo,ibleo),prod)
      do ir = 1,5,2
       localspinor(ir  ,:,2) = localspinor(ir  ,:,2) - prod(ir+1,:)
       localspinor(ir+1,:,2) = localspinor(ir+1,:,2) + prod(ir  ,:)
      enddo ! ir
      call mv(icount,GeeGooinv(:,:,8,ieo,ibl),spinor(:,:,4,ieo,ibleo),prod)
      do ir = 1,5,2
       localspinor(ir  ,:,2) = localspinor(ir  ,:,2) - prod(ir+1,:)
       localspinor(ir+1,:,2) = localspinor(ir+1,:,2) + prod(ir  ,:)
      enddo ! ir

      do isite = 1,nvhalf
       spinor(:,isite,:,ieo,ibleo) = localspinor(:,isite,:)
      enddo ! isite
     enddo ! ieo
    enddo ! ibleo

 end subroutine multG

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 subroutine smear(he,ho,u,ge,go,asmear,nsmear,bc,vecbl,vecblinv,myid,nn,ldiv, &
                  nms,lvbc,ib,lbd,iblv,MRT)
! Multiply the vector g(x) by a Gaussian smearing function as defined by
! equation (23) of Trottier, Phys. Rev. D55, 6844 (1997).
! INPUT:
!   u() contains the gauge fields for this sublattice.
!   ge() contains the Dirac spinor on globally even sites.
!   go() contains the Dirac spinor on globally odd sites.
!   asmear is the smearing "radius".
!   nsmear is the number of smearing factors applied.
! OUTPUT:
!   h(x) = (1 + asmear*[H_1(x,y)+H_2(x,y)+H_3(x,y)])^nsmear * g(y)
!        where H_mu = U_mu(x)*delta_(x,y-mu) + U_mu^dagger(x-mu)*delta_(x,y+mu)
!                   - 2*delta_(x,y)

    real(kind=KR),    intent(out), dimension(:,:,:,:,:) :: he, ho
    real(kind=KR),    intent(in),  dimension(:,:,:,:,:) :: u
    real(kind=KR),    intent(in),  dimension(:,:,:,:,:) :: ge, go
    real(kind=KR),    intent(in)                        :: asmear
    integer(kind=KI), intent(in)                        :: nsmear, myid, MRT
    integer(kind=KI), intent(in),  dimension(:)         :: bc, nms
    integer(kind=KI), intent(in),  dimension(:,:)       :: vecbl, vecblinv
    integer(kind=KI), intent(in),  dimension(:,:)       :: nn, iblv
    logical,          intent(in),  dimension(:)         :: ldiv
    integer(kind=KI), intent(in),  dimension(:,:,:)     :: lvbc
    integer(kind=KI), intent(in),  dimension(:,:,:,:)   :: ib
    logical,          intent(in),  dimension(:,:)       :: lbd

    integer(kind=KI)                          :: ismear, mu, ieo, jeo1, jeo2, &
                                                 ibl, jbl, ibleo, jbleo, id, &
                                                 i, j
    real(kind=KR), dimension(6,nvhalf,4,2,16) :: hebit, hobit
    real(kind=KR), dimension(6,nvhalf,4)      :: ha, hb
    real(kind=KR), dimension(6,ntotal,4)      :: hc
    real(kind=KR)                             :: bsmear

! A useful definition.
    bsmear = 1.0_KR - 6.0_KR*asmear

! Initialize the output vectors.
    do i = 1,nvhalf
     he(:,i,:,:,:) = ge(:,i,:,:,:)
     ho(:,i,:,:,:) = go(:,i,:,:,:)
    enddo ! i

! Sum over number of smearings.
    do ismear = 1,nsmear

! Compute the increment for each globally even site (gblclr=1).
     hebit = 0.0_KR
     do ieo = 1,2
      jeo1 = 3 - ieo
      do ibleo = 1,8
       ibl = vecbl(1,ibleo)
       do mu = 1,3
        jbl = iblv(ibl,mu)
        jbleo = vecblinv(2,jbl)
        if (lbd(jbl,mu)) then
         jeo2 = ieo
        else
         jeo2 = 3 - ieo
        endif
!-construct ha() = U_mu(x) * delta_{x,y-mu} * ho(y).
        do id = 1,4
         call setmvf(mu,bc,lbd(ibl,mu),ho(:,:,id,ieo,jbleo),ldiv(mu), &
                     ib(:,mu,1,jeo1),u(:,:,mu,ieo,ibl),ho(:,:,id,jeo1,jbleo), &
                     ha(:,:,id),lvbc(:,mu,ieo),1,nms,myid,nn,MRT)
!-construct hb() = U_mu^dagger(x-mu) * delta_{x,y+mu} * ho(y).
         call mdv(nvhalf,u(:,:,mu,jeo2,jbl),ho(:,:,id,jeo2,jbleo),hb(:,:,id))
!..next, move entries from hb(), defined at site y, to hc(), defined at site x.
!  If jbl is the second of two blocks in this direction, then move this
!  entry to its true owner (neighbouring process in the +mu direction)
!  using slidevector.  In any case, add the entry to the total.
         if (lbd(jbl,mu)) then
          do j = 1,nvhalf
           hc(:,j,id) = hb(:,j,id)
          enddo ! j
         else
          do j = 1,nvhalf
           hc(:,lvbc(j,mu,jeo2),id) = hb(:,j,id)
          enddo ! j
          if ( ldiv(mu) .or. ((.not.ldiv(mu)).and.bc(mu)/=1) ) &
           call slidevector(hc(:,:,id),ib(:,mu,1,ieo),mu,bc,nms, &
                            myid,nn,ldiv(mu),MRT)
         endif
        enddo ! id
        do i = 1,nvhalf
         hebit(:,i,:,ieo,ibleo) = hebit(:,i,:,ieo,ibleo) + ha(:,i,:) + hc(:,i,:)
        enddo ! i
       enddo ! mu
      enddo ! ibleo
     enddo ! ieo

! Compute the increment for each globally odd site (gblclr=2).
     hobit = 0.0_KR
     do ieo = 1,2
      jeo1 = 3 - ieo
      do ibleo = 1,8
       ibl = vecbl(2,ibleo)
       do mu = 1,3
        jbl = iblv(ibl,mu)
        jbleo = vecblinv(2,jbl)
        if (lbd(jbl,mu)) then
         jeo2 = ieo
        else
         jeo2 = 3 - ieo
        endif
!-construct ha() = U_mu(x) * delta_{x,y-mu} * he(y).
        do id = 1,4
         call setmvf(mu,bc,lbd(ibl,mu),he(:,:,id,ieo,jbleo),ldiv(mu), &
                     ib(:,mu,1,jeo1),u(:,:,mu,ieo,ibl),he(:,:,id,jeo1,jbleo), &
                     ha(:,:,id),lvbc(:,mu,ieo),1,nms,myid,nn,MRT)
!-construct hb() = U_mu^dagger(x-mu) * delta_{x,y+mu} * he(y).
         call mdv(nvhalf,u(:,:,mu,jeo2,jbl),he(:,:,id,jeo2,jbleo),hb(:,:,id))
!..next, move entries from hb(), defined at site y, to hc(), defined at site x.
!  If jbl is the second of two blocks in this direction, then move this
!  entry to its true owner (neighbouring process in the +mu direction)
!  using slidevector.  In any case, add the entry to the total.
         if (lbd(jbl,mu)) then
          do j = 1,nvhalf
           hc(:,j,id) = hb(:,j,id)
          enddo ! j
         else
          do j = 1,nvhalf
           hc(:,lvbc(j,mu,jeo2),id) = hb(:,j,id)
          enddo ! j
          if ( ldiv(mu) .or. ((.not.ldiv(mu)).and.bc(mu)/=1) ) &
           call slidevector(hc(:,:,id),ib(:,mu,1,ieo),mu,bc,nms, &
                            myid,nn,ldiv(mu),MRT)
         endif
        enddo ! id
        do i = 1,nvhalf
         hobit(:,i,:,ieo,ibleo) = hobit(:,i,:,ieo,ibleo) + ha(:,i,:) + hc(:,i,:)
        enddo ! i
       enddo ! mu
      enddo ! ibleo
     enddo ! ieo

! Update the quark fields for this smearing iteration.
     do i = 1,nvhalf
      he(:,i,:,:,:) = bsmear*he(:,i,:,:,:) + asmear*hebit(:,i,:,:,:)
      ho(:,i,:,:,:) = bsmear*ho(:,i,:,:,:) + asmear*hobit(:,i,:,:,:)
     enddo ! i

    enddo ! ismear

 end subroutine smear

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 subroutine zsmear(be,bo,src,icsrc,idsrc,myid,iblv,vecblinv)
! Define the source vector to be used for propagator construction.
! INPUT:
!   src is the location of the point source on the global lattice.

    real(kind=KR),    intent(out), dimension(:,:,:,:,:) :: be, bo
    integer(kind=KI), intent(in),  dimension(:)         :: src
    integer(kind=KI), intent(in)                        :: icsrc, idsrc, myid
    integer(kind=KI), intent(in),  dimension(:,:)       :: iblv, vecblinv

    integer(kind=KI), dimension(4) :: np, ip, srcbit
    integer(kind=KI) :: intx, inty, intz, intt, isite, ieo, ibl, mu, ibleo
    integer(kind=KI) :: ix,iy,iz
    integer(kind=KI) :: ierr, ididit

! Some initializations.
    np(1) = npx
    np(2) = npy
    np(3) = npz
    np(4) = npt
    call atoc(myid,np,ip)

    ididit = 0

!    print *, "myid,ip=",myid,ip

!    call MPI_BARRIER(MPI_COMM_WORLD, ierr)
!    stop

! Initialize the source vector.
    be = 0.0_KR
    bo = 0.0_KR

    do ix=1,nx
      do iy=1,ny
        do iz=1,nz

! srcbit is the position of the source on this sublattice.
    srcbit(1) = ix - ip(1)*nx/npx
    srcbit(2) = iy - ip(2)*ny/npy
    srcbit(3) = iz - ip(3)*nz/npz
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
      be(2*icsrc-1,isite,idsrc,ieo,ibleo) = 1.0_KR
     else
      bo(2*icsrc-1,isite,idsrc,ieo,ibleo) = 1.0_KR
     endif
     ididit = 1
    endif

    enddo ! iz
   enddo !iy
  enddo ! ix 

! if (ididit == 1) then
!    print *, "myididit", myid
! end if ! ididit

! call MPI_BARRIER(MPI_COMM_WORLD, ierr)
! stop

 end subroutine zsmear 

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


 subroutine mulbac(h,u,g,mu,idag,coact,bc,gblclr,vecbl,vecblinv,myid,nn,ldiv, &
                   nms,lvbc,ib,lbd,iblv,MRT)
! Add to h(x) the backward propagation term from the action:
! h(x) = h(x) + (1 - gamma_mu) U_mu(x) delta_{x,y-mu} g(y)
! To get the Lagrangian term, one would premultiply this addition by gbar(x).
! ***OR*** compute with the opposite sign on gamma_mu (for use in M^dagger):
! h(x) = h(x) + (1 + gamma_mu) U_mu(x) delta_{x,y-mu} g(y)
! INPUT:
!   h() contains the previously computed portion of the vector h=M*g.
!   u() contains the gauge fields for this sublattice.
!   g() contains the Dirac spinor for this sublattice.
!   mu is the direction of the gauge field (and Dirac gamma).
!   idag=0 is the Lagrangian term; idag=1 has the opposite sign for gamma_mu.
!   coact(2,4,4) contains the coefficients from the action.
!   bc(mu) = 1,-1,0 for periodic,antiperiodic,fixed boundary conditions
!            respectively.
!   gblclr is the "global checkerboard colour" of the gbar(x) spinor,
!          and thus of h(x).  g(y) spinor has the opposite checkerboard colour.
!          NOTE: This is NOT ieo, because of the lattice's blocked structure.
!                gblclr=1 means sites on blocks 1,4,6,7,10,11,13,16.
!                gblclr=2 means sites on blocks 2,3,5,8,9,12,14,15.
!   vecbl() defines the checkerboarding of the entire lattice.
!           In particular, the blocks vecbl(1,:), including ieo=1 and ieo=2,
!           form a checkerboard pattern.  So do the complement: vecbl(2,:).
! OUTPUT:
!   h(x) = h(x) + (1 - gamma_mu) U_mu(x) delta_{x,y-mu} g(y)
!   expected size: h(a,b,c),
!                  where a=1..6 is the real/imaginary and colour index.
!                        b=1,2,3,...,nvhalf combines colour and position.
!                        c=1,2,3,4 is the Dirac spinor index.

    real(kind=KR),    intent(inout), dimension(:,:,:,:,:) :: h, g
    real(kind=KR),    intent(in),    dimension(:,:,:,:,:) :: u
    integer(kind=KI), intent(in)                 :: mu, idag, gblclr, myid, MRT
    real(kind=KR),    intent(in),    dimension(:,:,:)     :: coact
    integer(kind=KI), intent(in),    dimension(:)         :: bc, nms
    integer(kind=KI), intent(in),    dimension(:,:)       :: vecbl, vecblinv
    integer(kind=KI), intent(in),    dimension(:,:)       :: nn, iblv
    logical,          intent(in),    dimension(:)         :: ldiv
    integer(kind=KI), intent(in),    dimension(:,:,:)     :: lvbc
    integer(kind=KI), intent(in),    dimension(:,:,:,:)   :: ib
    logical,          intent(in),    dimension(:,:)       :: lbd

    integer(kind=KI) :: ieo, jeo, ibl, jbl, ibleo, jbleo, id, i, ir, icount
    real(kind=KR), dimension(6,nvhalf,4) :: wx

    icount = nvhalf

!*Consider each site having Dirac global checkerboard colour = gblclr.
    do ieo = 1,2
     jeo = 3 - ieo
     do ibleo = 1,8
      ibl = vecbl(gblclr,ibleo)
      jbl = iblv(ibl,mu)
      jbleo = vecblinv(2,jbl)

!-first, construct wx() = U_mu(x) * delta_{x,y-mu} * g(y).
      do id = 1,4
       call setmvf(mu,bc,lbd(ibl,mu),g(:,:,id,ieo,jbleo),         &
                   ldiv(mu),ib(:,mu,1,jeo),u(:,:,mu,ieo,ibl),     &
                   g(:,:,id,jeo,jbleo),wx(:,:,id),lvbc(:,mu,ieo), &
                   1,nms,myid,nn,MRT)
      enddo ! id

!-next, construct h() = h() + [coact(3,1,mu) - coact(3,2,mu)*gamma_mu] * wx().
      if (idag==0) then
       select case(mu)
        case(1)
         do i = 1,nvhalf
          do ir = 1,6 ! 6=nri*nc
           h(ir,i,1,ieo,ibleo) = h(ir,i,1,ieo,ibleo) &
                          + coact(3,1,mu)*wx(ir,i,1) - coact(3,2,mu)*wx(ir,i,4)
           h(ir,i,2,ieo,ibleo) = h(ir,i,2,ieo,ibleo) &
                          + coact(3,1,mu)*wx(ir,i,2) - coact(3,2,mu)*wx(ir,i,3)
           h(ir,i,3,ieo,ibleo) = h(ir,i,3,ieo,ibleo) &
                          + coact(3,1,mu)*wx(ir,i,3) - coact(3,2,mu)*wx(ir,i,2)
           h(ir,i,4,ieo,ibleo) = h(ir,i,4,ieo,ibleo) &
                          + coact(3,1,mu)*wx(ir,i,4) - coact(3,2,mu)*wx(ir,i,1)
          enddo ! ir
         enddo ! i
        case(2)
         do i = 1,nvhalf
          do ir = 1,5,2 ! 6=nri*nc
           h(ir  ,i,1,ieo,ibleo) = h(ir  ,i,1,ieo,ibleo) &
                      + coact(3,1,mu)*wx(ir  ,i,1) - coact(3,2,mu)*wx(ir+1,i,4)
           h(ir+1,i,1,ieo,ibleo) = h(ir+1,i,1,ieo,ibleo) &
                      + coact(3,1,mu)*wx(ir+1,i,1) + coact(3,2,mu)*wx(ir  ,i,4)
           h(ir  ,i,2,ieo,ibleo) = h(ir  ,i,2,ieo,ibleo) &
                      + coact(3,1,mu)*wx(ir  ,i,2) + coact(3,2,mu)*wx(ir+1,i,3)
           h(ir+1,i,2,ieo,ibleo) = h(ir+1,i,2,ieo,ibleo) &
                      + coact(3,1,mu)*wx(ir+1,i,2) - coact(3,2,mu)*wx(ir  ,i,3)
           h(ir  ,i,3,ieo,ibleo) = h(ir  ,i,3,ieo,ibleo) &
                      + coact(3,1,mu)*wx(ir  ,i,3) - coact(3,2,mu)*wx(ir+1,i,2)
           h(ir+1,i,3,ieo,ibleo) = h(ir+1,i,3,ieo,ibleo) &
                      + coact(3,1,mu)*wx(ir+1,i,3) + coact(3,2,mu)*wx(ir  ,i,2)
           h(ir  ,i,4,ieo,ibleo) = h(ir  ,i,4,ieo,ibleo) &
                      + coact(3,1,mu)*wx(ir  ,i,4) + coact(3,2,mu)*wx(ir+1,i,1)
           h(ir+1,i,4,ieo,ibleo) = h(ir+1,i,4,ieo,ibleo) &
                      + coact(3,1,mu)*wx(ir+1,i,4) - coact(3,2,mu)*wx(ir  ,i,1)
          enddo ! ir
         enddo ! i
        case(3)
         do i = 1,nvhalf
          do ir = 1,6 ! 6=nri*nc
           h(ir,i,1,ieo,ibleo) = h(ir,i,1,ieo,ibleo) &
                          + coact(3,1,mu)*wx(ir,i,1) - coact(3,2,mu)*wx(ir,i,3)
           h(ir,i,2,ieo,ibleo) = h(ir,i,2,ieo,ibleo) &
                          + coact(3,1,mu)*wx(ir,i,2) + coact(3,2,mu)*wx(ir,i,4)
           h(ir,i,3,ieo,ibleo) = h(ir,i,3,ieo,ibleo) &
                          + coact(3,1,mu)*wx(ir,i,3) - coact(3,2,mu)*wx(ir,i,1)
           h(ir,i,4,ieo,ibleo) = h(ir,i,4,ieo,ibleo) &
                          + coact(3,1,mu)*wx(ir,i,4) + coact(3,2,mu)*wx(ir,i,2)
          enddo ! ir
         enddo ! i
        case(4)
         do i = 1,nvhalf
          do ir = 1,6 ! 6=nri*nc
           h(ir,i,1,ieo,ibleo) = h(ir,i,1,ieo,ibleo) &
                               + (coact(3,1,mu)-coact(3,2,mu))*wx(ir,i,1)
           h(ir,i,2,ieo,ibleo) = h(ir,i,2,ieo,ibleo) &
                               + (coact(3,1,mu)-coact(3,2,mu))*wx(ir,i,2)
           h(ir,i,3,ieo,ibleo) = h(ir,i,3,ieo,ibleo) &
                               + (coact(3,1,mu)+coact(3,2,mu))*wx(ir,i,3)
           h(ir,i,4,ieo,ibleo) = h(ir,i,4,ieo,ibleo) &
                               + (coact(3,1,mu)+coact(3,2,mu))*wx(ir,i,4)
          enddo ! ir
         enddo ! i
        case default
         open(unit=8,file="DIRACOPS.ERROR",action="write",status="replace", &
              form="formatted")
          write(unit=8,fmt=*) "subroutine mulbac: mu =", mu
         close(unit=8,status="keep")
         stop
       end select
!-or construct h() = h() + [coact(3,1,mu) + coact(3,2,mu)*gamma_mu] * wx().
      elseif (idag==1) then
       select case(mu)
        case(1)
         do i = 1,nvhalf
          do ir = 1,6 ! 6=nri*nc
           h(ir,i,1,ieo,ibleo) = h(ir,i,1,ieo,ibleo) &
                          + coact(3,1,mu)*wx(ir,i,1) + coact(3,2,mu)*wx(ir,i,4)
           h(ir,i,2,ieo,ibleo) = h(ir,i,2,ieo,ibleo) &
                          + coact(3,1,mu)*wx(ir,i,2) + coact(3,2,mu)*wx(ir,i,3)
           h(ir,i,3,ieo,ibleo) = h(ir,i,3,ieo,ibleo) &
                          + coact(3,1,mu)*wx(ir,i,3) + coact(3,2,mu)*wx(ir,i,2)
           h(ir,i,4,ieo,ibleo) = h(ir,i,4,ieo,ibleo) &
                          + coact(3,1,mu)*wx(ir,i,4) + coact(3,2,mu)*wx(ir,i,1)
          enddo ! ir
         enddo ! i
        case(2)
         do i = 1,nvhalf
          do ir = 1,5,2 ! 6=nri*nc
           h(ir  ,i,1,ieo,ibleo) = h(ir  ,i,1,ieo,ibleo) &
                      + coact(3,1,mu)*wx(ir  ,i,1) + coact(3,2,mu)*wx(ir+1,i,4)
           h(ir+1,i,1,ieo,ibleo) = h(ir+1,i,1,ieo,ibleo) &
                      + coact(3,1,mu)*wx(ir+1,i,1) - coact(3,2,mu)*wx(ir  ,i,4)
           h(ir  ,i,2,ieo,ibleo) = h(ir  ,i,2,ieo,ibleo) &
                      + coact(3,1,mu)*wx(ir  ,i,2) - coact(3,2,mu)*wx(ir+1,i,3)
           h(ir+1,i,2,ieo,ibleo) = h(ir+1,i,2,ieo,ibleo) &
                      + coact(3,1,mu)*wx(ir+1,i,2) + coact(3,2,mu)*wx(ir  ,i,3)
           h(ir  ,i,3,ieo,ibleo) = h(ir  ,i,3,ieo,ibleo) &
                      + coact(3,1,mu)*wx(ir  ,i,3) + coact(3,2,mu)*wx(ir+1,i,2)
           h(ir+1,i,3,ieo,ibleo) = h(ir+1,i,3,ieo,ibleo) &
                      + coact(3,1,mu)*wx(ir+1,i,3) - coact(3,2,mu)*wx(ir  ,i,2)
           h(ir  ,i,4,ieo,ibleo) = h(ir  ,i,4,ieo,ibleo) &
                      + coact(3,1,mu)*wx(ir  ,i,4) - coact(3,2,mu)*wx(ir+1,i,1)
           h(ir+1,i,4,ieo,ibleo) = h(ir+1,i,4,ieo,ibleo) &
                      + coact(3,1,mu)*wx(ir+1,i,4) + coact(3,2,mu)*wx(ir  ,i,1)
          enddo ! ir
         enddo ! i
        case(3)
         do i = 1,nvhalf
          do ir = 1,6 ! 6=nri*nc
           h(ir,i,1,ieo,ibleo) = h(ir,i,1,ieo,ibleo) &
                          + coact(3,1,mu)*wx(ir,i,1) + coact(3,2,mu)*wx(ir,i,3)
           h(ir,i,2,ieo,ibleo) = h(ir,i,2,ieo,ibleo) &
                          + coact(3,1,mu)*wx(ir,i,2) - coact(3,2,mu)*wx(ir,i,4)
           h(ir,i,3,ieo,ibleo) = h(ir,i,3,ieo,ibleo) &
                          + coact(3,1,mu)*wx(ir,i,3) + coact(3,2,mu)*wx(ir,i,1)
           h(ir,i,4,ieo,ibleo) = h(ir,i,4,ieo,ibleo) &
                          + coact(3,1,mu)*wx(ir,i,4) - coact(3,2,mu)*wx(ir,i,2)
          enddo ! ir
         enddo ! i
        case(4)
         do i = 1,nvhalf
          do ir = 1,6 ! 6=nri*nc
           h(ir,i,1,ieo,ibleo) = h(ir,i,1,ieo,ibleo) &
                               + (coact(3,1,mu)+coact(3,2,mu))*wx(ir,i,1)
           h(ir,i,2,ieo,ibleo) = h(ir,i,2,ieo,ibleo) &
                               + (coact(3,1,mu)+coact(3,2,mu))*wx(ir,i,2)
           h(ir,i,3,ieo,ibleo) = h(ir,i,3,ieo,ibleo) &
                               + (coact(3,1,mu)-coact(3,2,mu))*wx(ir,i,3)
           h(ir,i,4,ieo,ibleo) = h(ir,i,4,ieo,ibleo) &
                               + (coact(3,1,mu)-coact(3,2,mu))*wx(ir,i,4)
          enddo ! ir
         enddo ! i
        case default
         open(unit=8,file="DIRACOPS.ERROR",action="write",status="replace", &
              form="formatted")
          write(unit=8,fmt=*) "subroutine mulbac: mu =", mu
         close(unit=8,status="keep")
         stop
       end select
      endif

     enddo ! ibleo
    enddo ! ieo

 end subroutine mulbac

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 subroutine mulfor(h,u,g,mu,idag,coact,bc,gblclr,vecbl,vecblinv,myid,nn,ldiv, &
                   nms,lvbc,ib,lbd,iblv,MRT)
! Add to h(x) the forward propagation term from the action:
! h(x) = h(x) + (1 + gamma_mu) U_mu^dagger(x-mu) delta_{x,y+mu} g(y)
! To get the Lagrangian term, one would premultiply this addition by gbar(x).
! ***OR*** compute with the opposite sign on gamma_mu (for use in M^dagger):
! h(x) = h(x) + (1 - gamma_mu) U_mu^dagger(x-mu) delta_{x,y+mu} g(y)
! INPUT:
!   h() contains the previously computed portion of the vector h=M*g.
!   u() contains the gauge fields for this sublattice.
!   g() contains the Dirac spinor for this sublattice.
!   mu is the direction of the gauge field (and Dirac gamma).
!   idag=0 is the Lagrangian term; idag=1 has the opposite sign for gamma_mu.
!   coact(2,4,4) contains the coefficients from the action.
!   bc(mu) = 1,-1,0 for periodic,antiperiodic,fixed boundary conditions
!            respectively.
!   gblclr is the "global checkerboard colour" of the gbar(x) spinor,
!          and thus of h(x).  g(y) spinor has the opposite checkerboard colour.
!          NOTE: This is NOT ieo, because of the lattice's blocked structure.
!                gblclr=1 means sites on blocks 1,4,6,7,10,11,13,16.
!                gblclr=2 means sites on blocks 2,3,5,8,9,12,14,15.
!   vecbl() defines the checkerboarding of the entire lattice.
!           In particular, the blocks vecbl(1,:), including ieo=1 and ieo=2,
!           form a checkerboard pattern.  So do the complement: vecbl(2,:).
! OUTPUT:
!   h(x) = h(x) + (1 + gamma_mu) U_mu^dagger(x-mu) delta_{x,y+mu} g(y)
!   expected size: h(a,b,c),
!                  where a=1..6 is the real/imaginary and colour index.
!                        b=1,2,3,...,nvhalf combines colour and position.
!                        c=1,2,3,4 is the Dirac spinor index.

    real(kind=KR),    intent(inout), dimension(:,:,:,:,:) :: h
    real(kind=KR),    intent(in),    dimension(:,:,:,:,:) :: u, g
    integer(kind=KI), intent(in)                 :: mu, idag, gblclr, myid, MRT
    real(kind=KR),    intent(in),    dimension(:,:,:)     :: coact
    integer(kind=KI), intent(in),    dimension(:)         :: bc, nms
    integer(kind=KI), intent(in),    dimension(:,:)       :: vecbl, vecblinv
    integer(kind=KI), intent(in),    dimension(:,:)       :: nn, iblv
    logical,          intent(in),    dimension(:)         :: ldiv
    integer(kind=KI), intent(in),    dimension(:,:,:)     :: lvbc
    integer(kind=KI), intent(in),    dimension(:,:,:,:)   :: ib
    logical,          intent(in),    dimension(:,:)       :: lbd

    integer(kind=KI) :: ieo, jeo, ibl, jbl, ibleo, jbleo, id, i, j, ir, icount
    real(kind=KR), dimension(6,ntotal,4) :: wx
    real(kind=KR), dimension(6,nvhalf,4) :: wy

    icount = nvhalf

!*Consider each site having Dirac global checkerboard colour = gblclr.
    do ieo = 1,2
     do ibleo = 1,8
      ibl = vecbl(gblclr,ibleo)
      jbl = iblv(ibl,mu)
      jbleo = vecblinv(2,jbl)
      if (lbd(jbl,mu)) then
       jeo = ieo
      else
       jeo = 3 - ieo
      endif

!-consider each Dirac entry in turn.
      do id = 1,4
!..first, construct wy() = U_mu^dagger(x-mu) * delta_{x,y+mu} * g(y).
       call mdv(icount,u(:,:,mu,jeo,jbl),g(:,:,id,jeo,jbleo),wy(:,:,id))
!..next, move entries from wy(), defined at site y, to wx(), defined at site x.
!  If jbl is the second of two blocks in this direction, then move this
!  entry to its true owner (neighbouring process in the +mu direction)
!  using slidevector.  In any case, add the entry to the total.                 
       if (lbd(jbl,mu)) then
        do j = 1,nvhalf
         wx(:,j,id) = wy(:,j,id)
        enddo ! j
       else
        do j = 1,nvhalf
         wx(:,lvbc(j,mu,jeo),id) = wy(:,j,id)
        enddo ! j
        if ( ldiv(mu) .or. ((.not.ldiv(mu)).and.bc(mu)/=1) ) &
         call slidevector(wx(:,:,id),ib(:,mu,1,ieo),mu,bc,nms, &
                          myid,nn,ldiv(mu),MRT)
       endif
      enddo ! id

!-finally, construct h() = h() + (coact(3,1,mu)+coact(3,2,mu)*gamma_mu) * wx().
      if (idag==0) then
       select case(mu)
        case(1)
         do i = 1,nvhalf
          do ir = 1,6 ! 6=nri*nc
           h(ir,i,1,ieo,ibleo) = h(ir,i,1,ieo,ibleo) &
                          + coact(3,1,mu)*wx(ir,i,1) + coact(3,2,mu)*wx(ir,i,4)
           h(ir,i,2,ieo,ibleo) = h(ir,i,2,ieo,ibleo) &
                          + coact(3,1,mu)*wx(ir,i,2) + coact(3,2,mu)*wx(ir,i,3)
           h(ir,i,3,ieo,ibleo) = h(ir,i,3,ieo,ibleo) &
                          + coact(3,1,mu)*wx(ir,i,3) + coact(3,2,mu)*wx(ir,i,2)
           h(ir,i,4,ieo,ibleo) = h(ir,i,4,ieo,ibleo) &
                          + coact(3,1,mu)*wx(ir,i,4) + coact(3,2,mu)*wx(ir,i,1)
          enddo ! ir
         enddo ! i
        case(2)
         do i = 1,nvhalf
          do ir = 1,5,2 ! 6=nri*nc
           h(ir  ,i,1,ieo,ibleo) = h(ir  ,i,1,ieo,ibleo) &
                      + coact(3,1,mu)*wx(ir  ,i,1) + coact(3,2,mu)*wx(ir+1,i,4)
           h(ir+1,i,1,ieo,ibleo) = h(ir+1,i,1,ieo,ibleo) &
                      + coact(3,1,mu)*wx(ir+1,i,1) - coact(3,2,mu)*wx(ir  ,i,4)
           h(ir  ,i,2,ieo,ibleo) = h(ir  ,i,2,ieo,ibleo) &
                      + coact(3,1,mu)*wx(ir  ,i,2) - coact(3,2,mu)*wx(ir+1,i,3)
           h(ir+1,i,2,ieo,ibleo) = h(ir+1,i,2,ieo,ibleo) &
                      + coact(3,1,mu)*wx(ir+1,i,2) + coact(3,2,mu)*wx(ir  ,i,3)
           h(ir  ,i,3,ieo,ibleo) = h(ir  ,i,3,ieo,ibleo) &
                      + coact(3,1,mu)*wx(ir  ,i,3) + coact(3,2,mu)*wx(ir+1,i,2)
           h(ir+1,i,3,ieo,ibleo) = h(ir+1,i,3,ieo,ibleo) &
                      + coact(3,1,mu)*wx(ir+1,i,3) - coact(3,2,mu)*wx(ir  ,i,2)
           h(ir  ,i,4,ieo,ibleo) = h(ir  ,i,4,ieo,ibleo) &
                      + coact(3,1,mu)*wx(ir  ,i,4) - coact(3,2,mu)*wx(ir+1,i,1)
           h(ir+1,i,4,ieo,ibleo) = h(ir+1,i,4,ieo,ibleo) &
                      + coact(3,1,mu)*wx(ir+1,i,4) + coact(3,2,mu)*wx(ir  ,i,1)
          enddo ! ir
         enddo ! i
        case(3)
         do i = 1,nvhalf
          do ir = 1,6 ! 6=nri*nc
           h(ir,i,1,ieo,ibleo) = h(ir,i,1,ieo,ibleo) &
                          + coact(3,1,mu)*wx(ir,i,1) + coact(3,2,mu)*wx(ir,i,3)
           h(ir,i,2,ieo,ibleo) = h(ir,i,2,ieo,ibleo) &
                          + coact(3,1,mu)*wx(ir,i,2) - coact(3,2,mu)*wx(ir,i,4)
           h(ir,i,3,ieo,ibleo) = h(ir,i,3,ieo,ibleo) &
                          + coact(3,1,mu)*wx(ir,i,3) + coact(3,2,mu)*wx(ir,i,1)
           h(ir,i,4,ieo,ibleo) = h(ir,i,4,ieo,ibleo) &
                          + coact(3,1,mu)*wx(ir,i,4) - coact(3,2,mu)*wx(ir,i,2)
          enddo ! ir
         enddo ! i
        case(4)
         do i = 1,nvhalf
          do ir = 1,6 ! 6=nri*nc
           h(ir,i,1,ieo,ibleo) = h(ir,i,1,ieo,ibleo) &
                               + (coact(3,1,mu)+coact(3,2,mu))*wx(ir,i,1)
           h(ir,i,2,ieo,ibleo) = h(ir,i,2,ieo,ibleo) &
                               + (coact(3,1,mu)+coact(3,2,mu))*wx(ir,i,2)
           h(ir,i,3,ieo,ibleo) = h(ir,i,3,ieo,ibleo) &
                               + (coact(3,1,mu)-coact(3,2,mu))*wx(ir,i,3)
           h(ir,i,4,ieo,ibleo) = h(ir,i,4,ieo,ibleo) &
                               + (coact(3,1,mu)-coact(3,2,mu))*wx(ir,i,4)
          enddo ! ir
         enddo ! i
        case default
         open(unit=8,file="DIRACOPS.ERROR",action="write",status="replace", &
              form="formatted")
          write(unit=8,fmt=*) "subroutine mulfor: mu =", mu
         close(unit=8,status="keep")
         stop
       end select
!-or construct h() = h() + (coact(3,1,mu) - coact(3,2,mu)*gamma_mu) * wx().
      elseif (idag==1) then
       select case(mu)
        case(1)
         do i = 1,nvhalf
          do ir = 1,6 ! 6=nri*nc
           h(ir,i,1,ieo,ibleo) = h(ir,i,1,ieo,ibleo) &
                          + coact(3,1,mu)*wx(ir,i,1) - coact(3,2,mu)*wx(ir,i,4)
           h(ir,i,2,ieo,ibleo) = h(ir,i,2,ieo,ibleo) &
                          + coact(3,1,mu)*wx(ir,i,2) - coact(3,2,mu)*wx(ir,i,3)
           h(ir,i,3,ieo,ibleo) = h(ir,i,3,ieo,ibleo) &
                          + coact(3,1,mu)*wx(ir,i,3) - coact(3,2,mu)*wx(ir,i,2)
           h(ir,i,4,ieo,ibleo) = h(ir,i,4,ieo,ibleo) &
                          + coact(3,1,mu)*wx(ir,i,4) - coact(3,2,mu)*wx(ir,i,1)
          enddo ! ir
         enddo ! i
        case(2)
         do i = 1,nvhalf
          do ir = 1,5,2 ! 6=nri*nc
           h(ir  ,i,1,ieo,ibleo) = h(ir  ,i,1,ieo,ibleo) &
                      + coact(3,1,mu)*wx(ir  ,i,1) - coact(3,2,mu)*wx(ir+1,i,4)
           h(ir+1,i,1,ieo,ibleo) = h(ir+1,i,1,ieo,ibleo) &
                      + coact(3,1,mu)*wx(ir+1,i,1) + coact(3,2,mu)*wx(ir  ,i,4)
           h(ir  ,i,2,ieo,ibleo) = h(ir  ,i,2,ieo,ibleo) &
                      + coact(3,1,mu)*wx(ir  ,i,2) + coact(3,2,mu)*wx(ir+1,i,3)
           h(ir+1,i,2,ieo,ibleo) = h(ir+1,i,2,ieo,ibleo) &
                      + coact(3,1,mu)*wx(ir+1,i,2) - coact(3,2,mu)*wx(ir  ,i,3)
           h(ir  ,i,3,ieo,ibleo) = h(ir  ,i,3,ieo,ibleo) &
                      + coact(3,1,mu)*wx(ir  ,i,3) - coact(3,2,mu)*wx(ir+1,i,2)
           h(ir+1,i,3,ieo,ibleo) = h(ir+1,i,3,ieo,ibleo) &
                      + coact(3,1,mu)*wx(ir+1,i,3) + coact(3,2,mu)*wx(ir  ,i,2)
           h(ir  ,i,4,ieo,ibleo) = h(ir  ,i,4,ieo,ibleo) &
                      + coact(3,1,mu)*wx(ir  ,i,4) + coact(3,2,mu)*wx(ir+1,i,1)
           h(ir+1,i,4,ieo,ibleo) = h(ir+1,i,4,ieo,ibleo) &
                      + coact(3,1,mu)*wx(ir+1,i,4) - coact(3,2,mu)*wx(ir  ,i,1)
          enddo ! ir
         enddo ! i
        case(3)
         do i = 1,nvhalf
          do ir = 1,6 ! 6=nri*nc
           h(ir,i,1,ieo,ibleo) = h(ir,i,1,ieo,ibleo) &
                          + coact(3,1,mu)*wx(ir,i,1) - coact(3,2,mu)*wx(ir,i,3)
           h(ir,i,2,ieo,ibleo) = h(ir,i,2,ieo,ibleo) &
                          + coact(3,1,mu)*wx(ir,i,2) + coact(3,2,mu)*wx(ir,i,4)
           h(ir,i,3,ieo,ibleo) = h(ir,i,3,ieo,ibleo) &
                          + coact(3,1,mu)*wx(ir,i,3) - coact(3,2,mu)*wx(ir,i,1)
           h(ir,i,4,ieo,ibleo) = h(ir,i,4,ieo,ibleo) &
                          + coact(3,1,mu)*wx(ir,i,4) + coact(3,2,mu)*wx(ir,i,2)
          enddo ! ir
         enddo ! i
        case(4)
         do i = 1,nvhalf
          do ir = 1,6 ! 6=nri*nc
           h(ir,i,1,ieo,ibleo) = h(ir,i,1,ieo,ibleo) &
                               + (coact(3,1,mu)-coact(3,2,mu))*wx(ir,i,1)
           h(ir,i,2,ieo,ibleo) = h(ir,i,2,ieo,ibleo) &
                               + (coact(3,1,mu)-coact(3,2,mu))*wx(ir,i,2)
           h(ir,i,3,ieo,ibleo) = h(ir,i,3,ieo,ibleo) &
                               + (coact(3,1,mu)+coact(3,2,mu))*wx(ir,i,3)
           h(ir,i,4,ieo,ibleo) = h(ir,i,4,ieo,ibleo) &
                               + (coact(3,1,mu)+coact(3,2,mu))*wx(ir,i,4)
          enddo ! ir
         enddo ! i
        case default
         open(unit=8,file="DIRACOPS.ERROR",action="write",status="replace", &
              form="formatted")
          write(unit=8,fmt=*) "subroutine mulfor: mu =", mu
         close(unit=8,status="keep")
         stop
       end select
      endif

     enddo ! ibleo
    enddo ! ieo

 end subroutine mulfor

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 subroutine setmvf(mu,bc,l1,b1,l,ibbit,a,b,c,lvbcbit,iflg,nms,myid,nn,MRT)
! Get Dirac vectors from a neighbouring process, and multiply each of them
! by a link to form part of the Dirac operator.

    integer(kind=KI), intent(in)                    :: mu, iflg, myid, MRT
    integer(kind=KI), intent(in),    dimension(:)   :: bc, nms
    logical,          intent(in)                    :: l1, l
    real(kind=KR),    intent(in),    dimension(:,:) :: b1
    integer(kind=KI), intent(in),    dimension(:)   :: ibbit
    real(kind=KR),    intent(in),    dimension(:,:) :: a
    real(kind=KR),    intent(inout), dimension(:,:) :: b
    real(kind=KR),    intent(out),   dimension(:,:) :: c
    integer(kind=KI), intent(in),    dimension(:)   :: lvbcbit
    integer(kind=KI), intent(in),    dimension(:,:) :: nn

    integer(kind=KI) :: icount
    icount = nvhalf

!*If this is the first of the two blocks in this direction, then no vectors
! are needed from a neighbouring process.
    if (l1) then
     if (iflg==1) then
      call mv(icount,a,b1,c)
     elseif (iflg==2) then
      call mdv(icount,a,b1,c)
     endif
!*If this is the second of the two blocks in this direction, then get the
! vectors from a neighbouring process if...
! ...there are multiple processes in the mu direction [l=.true.],
! OR
! ...there is a single process in the mu direction but with non-periodic
!    boundary conditions [(.not.l).and.bc(mu)/=1].
    else
     if ( l .or. ((.not.l).and.bc(mu)/=1) ) &
      call setvector(b,ibbit,mu,bc,nms,myid,nn,l,MRT)
     if (iflg==1) then
      call mvg(icount,a,b,c,lvbcbit)
     elseif (iflg==2) then
      call mdvg(icount,a,b,c,lvbcbit)
     endif
    endif

 end subroutine setmvf

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 subroutine mv(n,y,v,w)
! Matrix*vector multiplier: calculate w(i) = y(i)*v(i)
! INPUT:
!   y() is a set of 3x3 matrices
!       (e.g. all links in mu direction on even sites),
!       expected size: y(18,n)
!   v() is a vector of the form     / v(1)+i*v(2) \
!                               v = | v(3)+i*v(4) |, where i=sqrt(-1).
!                                   \ v(5)+i*v(6) /
!
!       expected size: v(6,n)
!   n is the number of y matrices = number of v vectors.
! OUTPUT:
!   w(i) = y(i)*v(i) is the output.
!        expected size: w(6,n)

    integer(kind=KI), intent(in)                  :: n
    real(kind=KR),    intent(in),  dimension(:,:) :: y, v
    real(kind=KR),    intent(out), dimension(:,:) :: w

    integer(kind=KI) :: i

    do i = 1,n
     w(1,i) = y(1 ,i) * v(1,i) - y(2 ,i) * v(2,i) &
            + y(3 ,i) * v(3,i) - y(4 ,i) * v(4,i) &
            + y(5 ,i) * v(5,i) - y(6 ,i) * v(6,i)
     w(2,i) = y(1 ,i) * v(2,i) + y(2 ,i) * v(1,i) &
            + y(3 ,i) * v(4,i) + y(4 ,i) * v(3,i) &
            + y(5 ,i) * v(6,i) + y(6 ,i) * v(5,i)
     w(3,i) = y(7 ,i) * v(1,i) - y(8 ,i) * v(2,i) &
            + y(9 ,i) * v(3,i) - y(10,i) * v(4,i) &
            + y(11,i) * v(5,i) - y(12,i) * v(6,i)
     w(4,i) = y(7 ,i) * v(2,i) + y(8 ,i) * v(1,i) &
            + y(9 ,i) * v(4,i) + y(10,i) * v(3,i) &
            + y(11,i) * v(6,i) + y(12,i) * v(5,i)
     w(5,i) = y(13,i) * v(1,i) - y(14,i) * v(2,i) &
            + y(15,i) * v(3,i) - y(16,i) * v(4,i) &
            + y(17,i) * v(5,i) - y(18,i) * v(6,i)
     w(6,i) = y(13,i) * v(2,i) + y(14,i) * v(1,i) &
            + y(15,i) * v(4,i) + y(16,i) * v(3,i) &
            + y(17,i) * v(6,i) + y(18,i) * v(5,i)
    enddo ! i

 end subroutine mv

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 subroutine mdv(n,y,v,w)
! Matrix*vector multiplier: calculate w(i) = y^dagger(i)*v(i)
! INPUT:
!   y() is a set of 3x3 matrices
!       (e.g. all links in mu direction on even sites),
!       expected size: y(18,n)
!   v() is a vector of the form     / v(1)+i*v(2) \
!                               v = | v(3)+i*v(4) |, where i=sqrt(-1)
!                                   \ v(5)+i*v(6) /
!
!       expected size: v(6,n)
!   n is the number of y matrices = number of v vectors.
! OUTPUT:
!   w(i) = y^dagger(i)*v(i) is the output.
!        expected size: w(6,n)

    integer(kind=KI), intent(in)                  :: n
    real(kind=KR),    intent(in),  dimension(:,:) :: y, v
    real(kind=KR),    intent(out), dimension(:,:) :: w

    integer(kind=KI) :: i

    do i = 1,n
     w(1,i) = y(1 ,i) * v(1,i) + y(2 ,i) * v(2,i) &
            + y(7 ,i) * v(3,i) + y(8 ,i) * v(4,i) &
            + y(13,i) * v(5,i) + y(14,i) * v(6,i)
     w(2,i) = y(1 ,i) * v(2,i) - y(2 ,i) * v(1,i) &
            + y(7 ,i) * v(4,i) - y(8 ,i) * v(3,i) &
            + y(13,i) * v(6,i) - y(14,i) * v(5,i)
     w(3,i) = y(3 ,i) * v(1,i) + y(4 ,i) * v(2,i) &
            + y(9 ,i) * v(3,i) + y(10,i) * v(4,i) &
            + y(15,i) * v(5,i) + y(16,i) * v(6,i)
     w(4,i) = y(3 ,i) * v(2,i) - y(4 ,i) * v(1,i) &
            + y(9 ,i) * v(4,i) - y(10,i) * v(3,i) &
            + y(15,i) * v(6,i) - y(16,i) * v(5,i)
     w(5,i) = y(5 ,i) * v(1,i) + y(6 ,i) * v(2,i) &
            + y(11,i) * v(3,i) + y(12,i) * v(4,i) &
            + y(17,i) * v(5,i) + y(18,i) * v(6,i)
     w(6,i) = y(5 ,i) * v(2,i) - y(6 ,i) * v(1,i) &
            + y(11,i) * v(4,i) - y(12,i) * v(3,i) &
            + y(17,i) * v(6,i) - y(18,i) * v(5,i)
    enddo ! i

 end subroutine mdv

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 subroutine mvg(n,y,v,w,j)
! Matrix*vector multiplier: calculate w(i) = y(i)*v(j(i))
! INPUT:
!   y() is a set of 3x3 matrices
!       (e.g. all links in mu direction on even sites),
!       expected size: y(18,n)
!   v() is a vector of the form     / v(1)+i*v(2) \
!                               v = | v(3)+i*v(4) |, where i=sqrt(-1).
!                                   \ v(5)+i*v(6) /
!
!       expected size: v(6,m)
!   m is the number of v vectors (not needed in this subroutine),
!   j(i) is typically the neighbouring site (to i) in the appropriate
!        direction so w(i) will represent part of a staple,
!        expected size: j(n)
!   n is the number of y matrices = number of entries in j().
! OUTPUT:
!   w(i) = y(i)*v(j(i)) is the output.
!        expected size: w(6,n)

    integer(kind=KI), intent(in)                  :: n
    integer(kind=KI), intent(in),  dimension(:)   :: j
    real(kind=KR),    intent(in),  dimension(:,:) :: y, v
    real(kind=KR),    intent(out), dimension(:,:) :: w

    integer(kind=KI) :: i, k

    do i = 1,n
     k = j(i)
     w(1,i) = y(1 ,i) * v(1,k) - y(2 ,i) * v(2,k) &
            + y(3 ,i) * v(3,k) - y(4 ,i) * v(4,k) &
            + y(5 ,i) * v(5,k) - y(6 ,i) * v(6,k)
     w(2,i) = y(1 ,i) * v(2,k) + y(2 ,i) * v(1,k) &
            + y(3 ,i) * v(4,k) + y(4 ,i) * v(3,k) &
            + y(5 ,i) * v(6,k) + y(6 ,i) * v(5,k)
     w(3,i) = y(7 ,i) * v(1,k) - y(8 ,i) * v(2,k) &
            + y(9 ,i) * v(3,k) - y(10,i) * v(4,k) &
            + y(11,i) * v(5,k) - y(12,i) * v(6,k)
     w(4,i) = y(7 ,i) * v(2,k) + y(8 ,i) * v(1,k) &
            + y(9 ,i) * v(4,k) + y(10,i) * v(3,k) &
            + y(11,i) * v(6,k) + y(12,i) * v(5,k)
     w(5,i) = y(13,i) * v(1,k) - y(14,i) * v(2,k) &
            + y(15,i) * v(3,k) - y(16,i) * v(4,k) &
            + y(17,i) * v(5,k) - y(18,i) * v(6,k)
     w(6,i) = y(13,i) * v(2,k) + y(14,i) * v(1,k) &
            + y(15,i) * v(4,k) + y(16,i) * v(3,k) &
            + y(17,i) * v(6,k) + y(18,i) * v(5,k)
    enddo ! i

 end subroutine mvg

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 subroutine mdvg(n,y,v,w,j)
! Matrix*vector multiplier: calculate w(i) = y^dagger(i)*v(j(i))
! INPUT:
!   y() is a set of 3x3 matrices
!       (e.g. all links in mu direction on even sites),
!       expected size: y(18,n)
!   v() is a vector of the form     / v(1)+i*v(2) \
!                               v = | v(3)+i*v(4) |, where i=sqrt(-1).
!                                   \ v(5)+i*v(6) /
!
!       expected size: v(6,m)
!   m is the number of v vectors (not needed in this subroutine),
!   j(i) is typically the neighbouring site (to i) in the appropriate
!        direction so w(i) will represent part of a staple,
!        expected size: j(n)
!   n is the number of y matrices = number of entries in j().
! OUTPUT:
!   w(i) = y^dagger(i)*v(j(i)) is the output.
!        expected size: w(6,n)

    integer(kind=KI), intent(in)                  :: n
    integer(kind=KI), intent(in),  dimension(:)   :: j
    real(kind=KR),    intent(in),  dimension(:,:) :: y, v
    real(kind=KR),    intent(out), dimension(:,:) :: w

    integer(kind=KI) :: i, k

    do i = 1,n
     k = j(i)
     w(1,i) = y(1 ,i) * v(1,k) + y(2 ,i) * v(2,k) &
            + y(7 ,i) * v(3,k) + y(8 ,i) * v(4,k) &
            + y(13,i) * v(5,k) + y(14,i) * v(6,k)
     w(2,i) = y(1 ,i) * v(2,k) - y(2 ,i) * v(1,k) &
            + y(7 ,i) * v(4,k) - y(8 ,i) * v(3,k) &
            + y(13,i) * v(6,k) - y(14,i) * v(5,k)
     w(3,i) = y(3 ,i) * v(1,k) + y(4 ,i) * v(2,k) &
            + y(9 ,i) * v(3,k) + y(10,i) * v(4,k) &
            + y(15,i) * v(5,k) + y(16,i) * v(6,k)
     w(4,i) = y(3 ,i) * v(2,k) - y(4 ,i) * v(1,k) &
            + y(9 ,i) * v(4,k) - y(10,i) * v(3,k) &
            + y(15,i) * v(6,k) - y(16,i) * v(5,k)
     w(5,i) = y(5 ,i) * v(1,k) + y(6 ,i) * v(2,k) &
            + y(11,i) * v(3,k) + y(12,i) * v(4,k) &
            + y(17,i) * v(5,k) + y(18,i) * v(6,k)
     w(6,i) = y(5 ,i) * v(2,k) - y(6 ,i) * v(1,k) &
            + y(11,i) * v(4,k) - y(12,i) * v(3,k) &
            + y(17,i) * v(6,k) - y(18,i) * v(5,k)
    enddo ! i

 end subroutine mdvg

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 subroutine setvector(avector,ibbit,mu,bc,nms,myid,nn,ldivmu,MRT)
! Send/receive Dirac vectors to/from a neighbouring process which are required
! for the Dirac operator and/or enforce temporal boundary conditions.
! INPUT:
!   avector() is the array containing the vectors to be transferred.
!             expected size: avector(6,ntotal)
!   ibbit() is a numbering of boundary sites for this block of the sublattice
!           taken from EITHER the +mu edge OR the -mu edge of the
!           neighbouring sublattice.
!           expected size: ibbit(nbmax)
!   mu defined the lattice face of interest (i.e. mu=1 or mu=maximum).
!   bc(mu) = 1,-1,0 for periodic,antiperiodic,fixed boundary conditions
!            respectively.
!   MRT is MPIREALTYPE.
! OUTPUT:
!   avector() is updated.

    real(kind=KR),    intent(inout), dimension(:,:) :: avector
    integer(kind=KI), intent(in),    dimension(:)   :: ibbit, bc, nms
    integer(kind=KI), intent(in)                    :: mu, myid, MRT
    integer(kind=KI), intent(in),    dimension(:,:) :: nn
    logical,          intent(in)                    :: ldivmu

    integer(kind=KI)                             :: i, ierr
    integer(kind=KI), dimension(MPI_STATUS_SIZE) :: istatus
    integer(kind=KI), dimension(4)               :: msgcnt, np, ip

! Locate vectors needed by neighbouring process and put them into the buffer.
    do i = 1,nms(mu)
     avector(:,nvhalf+i) = avector(:,ibbit(i))
    enddo ! i

! Send buffer vectors to neighbour; receive buffer vectors from other neighbour.
    if (ldivmu) then
     msgcnt = 6*nms
    !if (abs(nn(mu,2)) > 7) then
    !   print *,"myid=",myid,"sending to nn(mu,2)",nn(mu,2),"mu=",mu
    !endif ! abs()

     call MPI_SENDRECV_REPLACE(avector(1,nvhalf+1),msgcnt(mu),MRT,nn(mu,2), &
                             myid,nn(mu,1),nn(mu,1),MPI_COMM_WORLD,istatus,ierr)
    endif

! Impose periodic, antiperiodic or fixed boundary conditions.
    np(1) = npx
    np(2) = npy
    np(3) = npz
    np(4) = npt
    call atoc(myid,np,ip)
    if (ip(mu)==np(mu)-1) then
     do i = 1,nms(mu)
      avector(:,nvhalf+i) = real(bc(mu),KR)*avector(:,nvhalf+i)
     enddo ! i
    endif

 end subroutine setvector

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 subroutine slidevector(w,ibbit,mu,bc,nms,myid,nn,ldivmu,MRT)
! Move the Dirac vectors in w from each process to a neighbouring process.
! INPUT:
!   w() is the set of Dirac vectors to be transferred.
!   ibbit() is a numbering of boundary sites for this block of the sublattice
!           in EITHER the +mu direction OR the -mu direction.
!   mu is the direction of the Dirac vectors to be transferred.
!   bc(mu) = 1,-1,0 for periodic,antiperiodic,fixed boundary conditions
!            respectively.
!   MRT is MPIREALTYPE.
! OUTPUT:
!   w() is updated.

    real(kind=KR),    intent(inout), dimension(:,:) :: w
    integer(kind=KI), intent(in),    dimension(:)   :: ibbit, bc, nms
    integer(kind=KI), intent(in)                    :: mu, myid, MRT
    integer(kind=KI), intent(in),    dimension(:,:) :: nn
    logical,          intent(in)                    :: ldivmu

    integer(kind=KI)                             :: i, ierr
    integer(kind=KI), dimension(MPI_STATUS_SIZE) :: istatus
    integer(kind=KI), dimension(4)               :: msgcnt, np, ip

! Send buffer entries to neighbour; receive buffer entries from other neighbour.
    if (ldivmu) then
     msgcnt = 6*nms
     call MPI_SENDRECV_REPLACE(w(1,nvhalf+1),msgcnt(mu),MRT,nn(mu,1),myid, &
                               nn(mu,2),nn(mu,2),MPI_COMM_WORLD,istatus,ierr)
    endif

! Impose periodic, antiperiodic or fixed boundary conditions.
    np(1) = npx
    np(2) = npy
    np(3) = npz
    np(4) = npt
    call atoc(myid,np,ip)
    if (ip(mu)==0) then
     do i = 1,nms(mu)
      w(:,nvhalf+i) = real(bc(mu),KR)*w(:,nvhalf+i)
     enddo ! i
    endif

! Move vectors from the buffer to their true position on the lattice.
    do i = 1,nms(mu)
     w(:,ibbit(i)) = w(:,nvhalf+i)
    enddo ! i

 end subroutine slidevector

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 subroutine gamma5mult(he,u,GeeGooinv,ge,idag,coact,kappa,iflag,bc,vecbl, &
                    vecblinv,myid,nn,ldiv,nms,lvbc,ib,lbd,iblv,MRT)

! NOTE: this subroutine has been repalced with the subroutine twMshift becasue 
!       our colaboration found that it is not possible to shift on the even/odd
!       twisted mass matrix. One must use the full matrix.

! Generalization of subroutine Hdbletm to the case of twisted mass QCD associated
! with multiple shifting. Multiply the operator -i*gamma5/2*kappa*mu*(G_ee - H_eo*G_oo^(-1)*H_oe)  
! onto a given vector, ge().where G = 1 + 2*kappa*mu*i*gam5. 
!   kappa(1) is 1/(8+2*m_0) with m_0 from eq.(1.1) of JHEP08(2001)058.
!   kappa(2) is mu_q from eq.(1.1) of JHEP08(2001)058.
! For use in CGNE, idag=1 allows use of the operator's dagger instead.
! INPUT: iflag=-1 (Only valid for Wilson loop, NO CLOVER TERM!)
! Notice that this program multiplies Hdbletm by -i*gamma5.
 
 
    real(kind=KR),    intent(out),   dimension(:,:,:,:,:) :: he
    real(kind=KR),    intent(in),    dimension(:,:,:,:,:) :: u
    real(kind=KR),    intent(in),    dimension(:,:,:,:,:) :: GeeGooinv
    real(kind=KR),    intent(inout), dimension(:,:,:,:,:) :: ge
    integer(kind=KI), intent(in)                      :: idag, iflag, myid, MRT
    real(kind=KR),    intent(in),    dimension(:,:,:)     :: coact
    real(kind=KR),    intent(in),    dimension(:)         :: kappa
    real(kind=KR)                                         :: basekappa
    real(kind=KR)                                         :: basemu
    integer(kind=KI), intent(in),    dimension(:)         :: bc, nms
    integer(kind=KI), intent(in),    dimension(:,:)       :: vecbl, vecblinv
    integer(kind=KI), intent(in),    dimension(:,:)       :: nn, iblv
    logical,          intent(in),    dimension(:)         :: ldiv
    integer(kind=KI), intent(in),    dimension(:,:,:)     :: lvbc
    integer(kind=KI), intent(in),    dimension(:,:,:,:)   :: ib
    logical,          intent(in),    dimension(:,:)       :: lbd
 
    real(kind=KR), dimension(6,ntotal,4,2,8) :: htmp,he1tmp,he2tmp
    integer(kind=KI)                         :: gblclr, isite,i,j,ierr
    real(kind=KR)                            :: fac1, fac2, fac3, fac4
 
    basekappa = kappa(1)
    basemu = kappa(2)
! For Wilson mass only
    if (iflag == -1) then
! First build H_oe * g_e.
     gblclr = 2
     call Hsingle(he,u,ge,idag,coact,bc,gblclr,vecbl,vecblinv,myid,nn,ldiv, &
                  nms,lvbc,ib,lbd,iblv,MRT)
! Next build (1-2*kappa*mu*i*gam5) * H_oe * g_e.
     fac1 = 2.0_KR*basekappa*basemu
     do isite = 1,nvhalf
      htmp(:,isite,1,:,:) = he(:,isite,1,:,:) - fac1*he(:,isite,3,:,:)
      htmp(:,isite,2,:,:) = he(:,isite,2,:,:) - fac1*he(:,isite,4,:,:)
      htmp(:,isite,3,:,:) = he(:,isite,3,:,:) + fac1*he(:,isite,1,:,:)
      htmp(:,isite,4,:,:) = he(:,isite,4,:,:) + fac1*he(:,isite,2,:,:)
     enddo ! isite
! Next build H_eo * (1-2*kappa*mu*i*gam5) * H_oe * g_e.
     gblclr = 1
     call Hsingle(he,u,htmp,idag,coact,bc,gblclr,vecbl,vecblinv,myid,nn,ldiv, &
                  nms,lvbc,ib,lbd,iblv,MRT)
! Finally, construct -i*gamma5/2*kappa*mu*h_e = 
!                                 -i*gamma5/2*kappa*mu*[(1+2*kappa*mu*i*gam5)/kappa^2
!               - H_eo*(1-2*kappa*mu*i*gam5)*H_oe/(1+4*kappa^2*mu^2)] * g_e.
! NOTE: the matrix he was multiplied by -i!!!
!       and scaled by 1/2*kappa(1)*kappa(2)
! giving -i*gamma5*h_e/2*kappa(1)*kappa(2)
! Shifting will occur on kappa(2) -> newmu = kappa(2)(1-sig)
     fac1 = 1.0_KR/(basekappa**2)
     fac2 = 2.0_KR*basemu/basekappa
     fac3 = 1.0_KR/( 1.0_KR + (2.0_KR*basekappa*basemu)**2 )
     fac4 = 1.0_KR/(2.0_KR*basekappa*basemu)
!    if (myid==0)then
!      print *, "fac1=", fac1
!      print *, "fac2=", fac2
!      print *, "fac3=", fac3
!      print *, "fac4=", fac4
!    endif
! Create temp variables for he(:,:,1,:,:) and he(:,:,2,:,:)
! so that first iteration of gamma5 multiplication does not
! change these values.
     do isite = 1,nvhalf
      he1tmp(:,isite,1,:,:) = he(:,isite,1,:,:) 
      he2tmp(:,isite,2,:,:) = he(:,isite,2,:,:)
     enddo ! isite
     do isite = 1,nvhalf
!     he(:,isite,1,:,:) = -fac1*ge(:,isite,3,:,:) + fac2*ge(:,isite,1,:,:) &
!                       + fac3*he(:,isite,3,:,:)
!     he(:,isite,2,:,:) = -fac1*ge(:,isite,4,:,:) + fac2*ge(:,isite,2,:,:) &
!                       + fac3*he(:,isite,4,:,:)
!     he(:,isite,3,:,:) = fac1*ge(:,isite,1,:,:) + fac2*ge(:,isite,3,:,:) &
!                       - fac3*he1tmp(:,isite,1,:,:)
!     he(:,isite,4,:,:) = fac1*ge(:,isite,2,:,:) + fac2*ge(:,isite,4,:,:) &
!                       - fac3*he2tmp(:,isite,2,:,:)
!
      he(:,isite,1,:,:) = -fac4*fac1*ge(:,isite,3,:,:) + fac1*ge(:,isite,1,:,:) &
                        + fac4*fac3*he(:,isite,3,:,:)
      he(:,isite,2,:,:) = -fac4*fac1*ge(:,isite,4,:,:) + fac1*ge(:,isite,2,:,:) &
                        + fac4*fac3*he(:,isite,4,:,:)
      he(:,isite,3,:,:) = fac4*fac1*ge(:,isite,1,:,:) + fac1*ge(:,isite,3,:,:) &
                        - fac4*fac3*he1tmp(:,isite,1,:,:)
      he(:,isite,4,:,:) = fac4*fac1*ge(:,isite,2,:,:) + fac1*ge(:,isite,4,:,:) &
                        - fac4*fac3*he2tmp(:,isite,2,:,:)
 
     enddo ! isite
    elseif (iflag==-2)then 
     open(unit=8,file="CFGSPROPS.ERROR",action="write",status="replace" &
          ,form="formatted")
      write(unit=8,fmt=*) "Remember to turn off Clover term with multiple shifts!"
     close(unit=8,status="keep")
     stop
    endif

  end subroutine gamma5mult

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 subroutine twMshift(hwhole,u,GeeGooinv,gwhole,idag,coact,kappa,iflag,bc,vecbl, &
                     vecblinv,myid,nn,ldiv,nms,lvbc,ib,lbd,iblv,MRT)
! This subroutine is used to combine twisted and multi-massing. An attempt
!  to combine these two techniques using the even/odd technique
!  proves difficult (subroutine gamma5mult) so me must multi-mass
!  on the entire matrix. The matrix that we wish to build and shift on 
!  is...~DD
!
!                                          / 1         kappa*h_eo \
!     M^(shift) = 1 - i*gamma5/2*kappa*mu* |                      |
!                                          \ kappa*h_oe     1     /
!
!
!  where the kappa and mu values are gathered from...
!   kappa(1) is 1/(8+2*m_0) with m_0 from eq.(1.1) of JHEP08(2001)058.
!   kappa(2) is mu_q from eq.(1.1) of JHEP08(2001)058.
!
! We will now be able to shift on the parameter mu(new) = mu(old)*(1-sigma)

! INPUT: iflag=-1 (Only valid for Wilson loop, NO CLOVER TERM!)
!        ge,go are the even/odd vectors passed in/out of this routine.
! OUTPUT: the matrix mentioned above multiplied be the input vectors
!         ge,go. The output is in h and the modified gin vectors.

 
 
    real(kind=KR),    intent(out),   dimension(:,:,:,:,:,:) :: hwhole
    real(kind=KR),    intent(inout), dimension(:,:,:,:,:,:) :: gwhole
    real(kind=KR),    intent(in),    dimension(:,:,:,:,:)    :: u
    real(kind=KR),    intent(in),    dimension(:,:,:,:,:)    :: GeeGooinv
    integer(kind=KI), intent(in)                      :: idag, iflag, myid, MRT
    real(kind=KR),    intent(in),    dimension(:,:,:)        :: coact
    real(kind=KR),    intent(in),    dimension(:)            :: kappa
    real(kind=KR)                                            :: basekappa
    real(kind=KR)                                            :: basemu
    integer(kind=KI), intent(in),    dimension(:)            :: bc, nms
    integer(kind=KI), intent(in),    dimension(:,:)          :: vecbl, vecblinv
    integer(kind=KI), intent(in),    dimension(:,:)          :: nn, iblv
    logical,          intent(in),    dimension(:)            :: ldiv
    integer(kind=KI), intent(in),    dimension(:,:,:)        :: lvbc
    integer(kind=KI), intent(in),    dimension(:,:,:,:)      :: ib
    logical,          intent(in),    dimension(:,:)          :: lbd
 
    real(kind=KR), dimension(6,ntotal,4,2,8)   :: he,ho
    real(kind=KR), dimension(6,ntotal,4,2,8)   :: htmp 
    integer(kind=KI)                           :: gblclr, isite,i,j,ierr
    real(kind=KR)                              :: fac1, fac2, fac3, fac4
 

!*SOME BASIC PARAMETERS THAT ARE DEFINED IN module basics:(copied)
! nps = total number of processes.
! nvhalf is the number of even (or odd) lattices sites per process per block.
! nbmax is the number of spacetime sites in the largest boundary
!       shared by two processes.  (The boundary is 3-D.)
! ntotal runs over even or odd sites on one block of a sublattice, plus
!            room for the largest boundary shared by two processes.
! nwhole is the number of even AND odd lattice sites per process per block.~DD

    basekappa = kappa(1)
    basemu = kappa(2)
! For Wilson mass only
    if (iflag == -1) then

! First make the even solution.
! Build h_eo * g_o.
     gblclr = 1
     call Hsingle(he,u,gwhole(:,:,:,:,:,2),idag,coact,bc,gblclr,vecbl,vecblinv,myid,nn,ldiv, &
                  nms,lvbc,ib,lbd,iblv,MRT)

! Next build (-i*gamma5*kappa/2*kappa*mu)* h_eo * g_o.

     fac1 = 1.0_KR/(2.0_KR*basemu)

     do isite = 1,nvhalf
      htmp(:,isite,1,:,:) = - fac1*he(:,isite,3,:,:)
      htmp(:,isite,2,:,:) = - fac1*he(:,isite,4,:,:)
      htmp(:,isite,3,:,:) =   fac1*he(:,isite,1,:,:)
      htmp(:,isite,4,:,:) =   fac1*he(:,isite,2,:,:)
     enddo ! isite

! Now construct the new even solutions
!    ~h_e = (1 - i*gamma5/2*kappa*mu)ge - (i*gamma5/2*mu)*h_eo*g_o

     fac2 = 1.0_KR/(2.0_KR*basekappa*basemu)

     do isite =1,nvhalf
       hwhole(:,isite,1,:,:,1) = gwhole(:,isite,1,:,:,1) - fac2*gwhole(:,isite,3,:,:,1) &
                         + htmp(:,isite,1,:,:)
       hwhole(:,isite,2,:,:,1) = gwhole(:,isite,2,:,:,1) - fac2*gwhole(:,isite,4,:,:,1) &
                         + htmp(:,isite,2,:,:)
       hwhole(:,isite,3,:,:,1) = gwhole(:,isite,3,:,:,1) + fac2*gwhole(:,isite,1,:,:,1) &
                         + htmp(:,isite,3,:,:)
       hwhole(:,isite,4,:,:,1) = gwhole(:,isite,4,:,:,1) + fac2*gwhole(:,isite,2,:,:,1) &
                         + htmp(:,isite,4,:,:)
     enddo ! isite

! Now make the odd solution.
! Build h_oe * g_e.
     gblclr = 2
     call Hsingle(ho,u,gwhole(:,:,:,:,:,1),idag,coact,bc,gblclr,vecbl,vecblinv,myid,nn,ldiv, &
                  nms,lvbc,ib,lbd,iblv,MRT)
! Next build (-i*gamma5*kappa/2*kappa*mu)* h_oe * g_e.
     do isite = 1,nvhalf
      htmp(:,isite,1,:,:) = - fac1*ho(:,isite,3,:,:)
      htmp(:,isite,2,:,:) = - fac1*ho(:,isite,4,:,:)
      htmp(:,isite,3,:,:) =   fac1*ho(:,isite,1,:,:)
      htmp(:,isite,4,:,:) =   fac1*ho(:,isite,2,:,:)
     enddo ! isite

! Now construct the new odd solutions ... 
!    ~h_o = (1 - i*gamma5/2*kappa*mu)go - (i*gamma5/2*mu)*h_oe*g_e

     do isite =1,nvhalf
       hwhole(:,isite,1,:,:,2) = gwhole(:,isite,1,:,:,2) - fac2*gwhole(:,isite,3,:,:,2) &
                         + htmp(:,isite,1,:,:)
       hwhole(:,isite,2,:,:,2) = gwhole(:,isite,2,:,:,2) - fac2*gwhole(:,isite,4,:,:,2) &
                         + htmp(:,isite,2,:,:)
       hwhole(:,isite,3,:,:,2) = gwhole(:,isite,3,:,:,2) + fac2*gwhole(:,isite,1,:,:,2) &
                         + htmp(:,isite,3,:,:)
       hwhole(:,isite,4,:,:,2) = gwhole(:,isite,4,:,:,2) + fac2*gwhole(:,isite,2,:,:,2) &
                         + htmp(:,isite,4,:,:)
     enddo ! isite
 
    elseif (iflag==-2)then 
     open(unit=8,file="CFGSPROPS.ERROR",action="write",status="replace" &
          ,form="formatted")
      write(unit=8,fmt=*) "Remember to turn off Clover term with multiple shifts!"
      print *, "Need to remove Clover term for this inverter!"
     close(unit=8,status="keep")
     stop
    endif

  end subroutine twMshift

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 end module diracops

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
