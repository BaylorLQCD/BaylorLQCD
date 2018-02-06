! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
! latdims.f90, Randy Lewis, randy.lewis@uregina.ca
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
!   nx = number of lattice sites in the x direction on the full lattice
!        (similarly for ny, nz, nt)
!   npx = number of processes in the x direction
!         (similarly for npy, npz, npt)
! REQUIRED: nx/4/npx must be an integer
!           where 4 = (2 blocks) x (2 checkerboard even/odd types)
!           (similarly for y,z,t directions)
!
!   nk1 = number of hopping parameters for multi-mass quark propagators,
!         BUT IF NO MULTI-MASS HOPPING PARAMETERS ARE DESIRED, CHOOSE nk1=1.
!
! WARNING ~ if inverter ==6 then nk1 is the number of shifts desired for
!           the mu value
!
!   twoRmax should be set to the maximum value of 2*Rmax(iwilo).
!           IF nwilo=0, THEN CHOOSE twoRmax=1.
!
!   npathsmax must be at least at large as npaths.
!             IF nwilo=0 OR npaths=0, THEN CHOOSE npathsmax=1.
!
!   nmaxGMRES is the maximum allowed value of the "n" in GMRES(n).
!             IF inverter/=4 THEN CHOOSE nmaxGMRES=1.
!
!   kmaxGMRES is the maximum allowed value of the "k" in GMRES-DR(n,k).
!             IF inverter/=4 THEN CHOOSE kmaxGMRES=1.
!
!   kmax is the number of hopping parameters
!
!   nshifts is the number of twisted masses used in the mulit-mass
!           calculation.
!
!   kcyclim is the maximum allowed cycles of gmresDRSHIFT with 
!           mulitiple right-hand sides(RHS) and multiple twisted-mass.
!           ONLY valid for inverter==6.
!
!   nirloop=1 means that the Wilson coefficient r=+1 
!   nirloop=2 means that the Wilson coefficient r=-1
!   nirloop=3 means that we average the disconnected part over r=+1 and -1.
!
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


     module latdims
     implicit none

!Quenched Wilson

     !beta=6.0 and beta=5.85
      !BS integer, public, parameter  :: nx = 8, ny = 8, nz = 8, nt = 16, &
                                      !npx = 2, npy = 2, npz = 2, npt = 1 ! BS
     ! integer, public, parameter  :: nx = 8, ny = 8, nz = 8, nt = 8, & 
        !                             npx = 2, npy = 2, npz = 2, npt =1 !BS changed from 8 to 4 for values of nx,ny,nz,nt TW 10/13/17        
                                 ! npx = 1, npy = 1, npz = 1, npt = 1
      integer, public, parameter  :: nx = 24, ny = 24, nz = 24, nt = 32, &
                                    npx = 3, npy = 3, npz = 2, npt = 8
                                   ! npx = 3, npy = 3, npz = 6, npt = 8
      

     integer, public, parameter  :: nk1=1,nirloop=1
     integer, public, parameter  :: twoRmax=1, npathsmax=1, nmaxGMRES=510, &
                                    kmaxGMRES=400

     integer, public, parameter  :: kmax=9,nshifts=1,kcyclim=5000!BS
   
     
! The below is old Dean testing code.  
!
! NOTE~ ONLY if inverter==6 then choose sigmamu  to desired values. nshifts is
!       the total number of values allowed. The number of sigmamu values assigned
!       must equal nshifts. Also, the number of propagtors, nk1,  corresponding to these
!       shifts must be the same value as nshifts.
!




     real(kind=8), public,  dimension(nshifts)  :: sigmamu
     real(kind=8), public,  dimension(nshifts)  :: kappamu

     public ::  twistedvalues

     contains
     subroutine twistedvalues
      use kinds

      sigmamu(1) =  -5.0_KR

!     sigmamu(1) =  0.0_KR
!     sigmamu(2) = -1.0_KR
!     sigmamu(3) = -2.0_KR
!     sigmamu(4) = -5.0_KR

      kappamu(1) = .15701

!     kappamu(1) = .15728
!     kappamu(2) = .15721
!     kappamu(3) = .15708
!     kappamu(4) = .15679

     end subroutine twistedvalues

     end module latdims

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
