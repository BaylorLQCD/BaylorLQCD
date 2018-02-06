!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
! cfgspropsmain.f90, Randy Lewis, randy.lewis@uregina.ca
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
! Main program for quenched lattice QCD parallel codes.
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  program cfgspropsmain

    use kinds
    use latdims
    use cfgsprops
    implicit none

    integer(kind=KI) :: newcalc, iseed, ncold, nsweep, n0sweep, n1sweep, nhit, &
                        tadupdate, tadpole, tadlandau, nwilo, npaths, nfuzz,   &
                        ifilesave, gitest, nsmear, nkappa, ntmqcd, &
                        itermin, docorr, inverter, numnoises, opSST, tSST, &
                        invSST, nkappacorr, ntmqcdcorr, ishift,  is, wilsonmm

   
    
    character(len=128), dimension(1024) :: rwdir
    character(len=128)                  :: cfgfile
    integer(kind=KI),  dimension(4)     :: gaction, landimp, bc
    real(kind=KR)                       :: beta, epsfuzz, alpha, omega, &
                                           thetamax, resmax, asmear, &
                                           cSWtmqcd, cSWloop, cSWSST
    real(kind=KR),     dimension(2)     :: kappaloop
    real(kind=KR),     dimension(4)     :: alat, tad
    real(kind=KR),     dimension(3)     :: Rmax
    integer(kind=KI),  dimension(3,4)   :: wlat
    integer(kind=KI),  dimension(0:9,4) :: src
    real(kind=KR),     dimension(9)     :: kappa
    real(kind=KR),     dimension(9,2)   :: mtmqcd, mSST
    real(kind=KR),     dimension(20,2)  :: mtmqcdloop
    real(kind=KR2)                      :: omegaM3R
    integer(kind=KI),  dimension(2)     :: nGMRES, nSST, nSSTcorr
    integer(kind=KI),  dimension(3)     :: pSST
    integer(kind=KI)                    :: nucleonparticle
    integer(kind=KI)                    :: ntmqcdloop,hoption

! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

! FILE MANAGEMENT:

! Define the read/write directory, "rwdir(i+1)", for each process, "i".
!-two processes per node on openSPACE
!  rwdir = "/home/darnelld/qqcd/scratch/52/"
!  do nhit = 0,25
!   write(unit=rwdir(nhit+1)(29:30),fmt="(i2.2)") 52+nhit/2
!  enddo ! nhit
!-one process per node on openSPACE
!   rwdir = "/scratch52/randy/temp/"
!   do nhit = 0,13
!    write(unit=rwdir(nhit+1)(9:10),fmt="(i2.2)") 52+nhit
!    write(unit=rwdir(nhit+14)(9:10),fmt="(i2.2)") 52+nhit
!   enddo ! nhit
!-single process mode
    rwdir = "/data/qcd/qqcd-hutcheson/scratch/"
 
! Choose newcalc=0 to use existing configurations, existing single propagators
!                  and existing double (SST) propagators.
! Choose newcalc=1 to use existing configurations and existing single
!                  propagators but create new double (SST) propagators.
! Choose newcalc=2 to use existing configurations but create new single
!                  propagators and new double (SST) propagators.
! Choose newcalc=3 to create new configurations and new single propagators and
!                  new double (SST) propagators.
    newcalc = 3

! Set ifilesave=4 to retain all configurations, all single propagators and all
!                 double (SST) propagators on disk.
! Set ifilesave=3 to retain all configurations and all single propagators on
!                 disk, but erase all double (SST) propagators.
! Set ifilesave=2 to retain all configurations and all double (SST) propagators
!                 on disk, but erase all single propagators.
! Set ifilesave=1 to retain all configurations on disk, but erase all
!                 single propagators and all double (SST) propagators.
! Set ifilesave=0 to erase all configurations and all single propagators and
!                 all double (SST) propagators.
! Note: The final configuration and its single and double propagators are
!       never erased.
! Note: The code will not erase existing files.

    ifilesave = 1

! Filenames must take the form AxxxxxyyyBCDE...
! where xxxxx is used by the program to number configurations
! and yyy is used by the program to number processes.
! Enter a filename template (if this is not a cold start, enter the exact
! value of xxxxx for the first filename to be read in):
      cfgfile = "cfgfiles/G00001000qwilson-16-32-b6"

! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

! RANDOM NUMBER SEED:

! Choose iseed to be a non-negative integer.
! (NOTE: For multiple processors running in single process mode, iseed
!  should be unique for each value of nhit.  For all other cases, iseed
!  can be a simple constant.)
    iseed = 13

! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

! GAUGE FIELD UPDATING:

! Thermalization from a cold start? [yes=1, no=0]
! (The value of ncold is irrelevant unless newcalc=3)
    ncold = 1

! IF newcalc=3 THEN
!   "nsweep" sweeps are performed, with the following subset written to disk:
!   n0sweep, n0sweep+n1sweep, n0sweep+2*n1sweep, n0sweep+3*n1sweep, ...
!   (n0sweep > 0 is required.)

!!!!!! COMMENT !!!!!!!

! NOTE ~ when starting from a specific configuration and newcalc=3, ncold=0,
!        nsweep=(# of desired configs from specified config) (7-31-05)

!!!!! END COMMENT !!!!!

! IF newcalc<3 THEN
!   "nsweep" configurations are read from the disk.  They must be number
!   sequentially, beginning with "cfgfile" as defined above.
! (The values of n0sweep and n1sweep are irrelevant unless newcalc=3)

! IF newcalc/=3 and you wish to use exsisting configurations then
! n0sweep and n1sweep are irrelevant and nsweep is the desired number
! of configurations. (7-31-05)

    nsweep  = 30000
    n0sweep = 20000
    n1sweep = 2000 

    !nsweep =  100
    !n0sweep = 10000
    !n1sweep = 2000

! Choose the number of pseudo-heatbath hits/site/sweep.
! (The value of nhit is irrelevant unless newcalc=3)
    nhit = 1

! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 3

! GAUGE ACTION:

! Nearest [1] or next-nearest [2] neighbour action?
! Enter x,y,z,t values, eg: 2,2,2,1 or 2 2 2 1
! (The values in action are irrelevant unless newcalc=3 and ncold=1,
! otherwise the true value is read from each configuration file.)
    gaction(1:4) = (/ 1, 1, 1, 1 /)

! Enter the ratios of lattice spacings: a_x/a_t , a_y/a_t , a_z/a_t , a_t/a_t
! (The values in alat are irrelevant unless newcalc=3 and ncold=1,
! otherwise the true value is read from each configuration file.)
    alat(1:4) = (/ 1.0_KR, 1.0_KR, 1.0_KR, 1.0_KR /)

! Choose the plaquette coupling.
! (The value of beta is irrelevant unless newcalc=3 and ncold=1,
! otherwise the true value is read from each configuration file.)
    beta = 6.0_KR
! Note: Although the codes are written from a tadpole improvement perspective,
!       they can be reinterpreted in other ways.  For two examples...

!********************************Iwasaki action*********************************
!*Based on QCD-TARO, hep-lat/9911033v2 page 6,                                 *
!*[see also R. Gupta, hep-lat/9807028 equations (12.14) and (12.26)],          *
!*it seems that the Iwasaki gauge action is attainable by using                *
!*         gaction(1:4) = (/ 2, 2, 2, 2 /),                                    *
!*         tad(1:4) = (/ 0.7424_KR, 0.7424_KR, 0.7424_KR, 0.7424_KR /), and    *
!*         tadupdate = 0                                                       *
!*The value of tad(1:4) comes from solving 1/(20u^2) = 0.09073.                *
!*(This action is mentioned by Luscher on page 203 of Lattice2001 proceedings.)*
!*NOTE: The value of beta is the coefficient of (5/3)*(1x1 plaquette).         *
!*******************************************************************************

!***********************************DBW2 action*********************************
!*Based on QCD-TARO, hep-lat/9911033v2 page 7,                                 *
!*[see also R. Gupta, hep-lat/9807028 equations (12.14) and (12.31)],          *
!*it seems that the DBW2 gauge action is attainable by using                   *
!*         gaction(1:4) = (/ 2, 2, 2, 2 /),                                    *
!*         tad(1:4) = (/ 0.6600_KR, 0.6600_KR, 0.6600_KR, 0.6600_KR /), and    *
!*         tadupdate = 0                                                       *
!*The value of tad(1:4) comes from solving 1/(20u^2) = 0.1148.                 *
!*NOTE: The value of beta is the coefficient of (5/3)*(1x1 plaquette).         *
!*      What is called "beta=0.87" in hep-lat/0211023 is obtainable by using   *
!*      beta = 1.394845*"beta" = 1.394845*0.87 = 1.213515                      *
!*******************************************************************************

! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

! AVERAGE PLAQUETTE AND TADPOLE FACTORS:

! Enter initial values for the tadpole factors: u_x, u_y, u_z, u_t
! (The values in tad are irrelevant unless newcalc=3 and ncold=1,
! otherwise the true value is read from each configuration file.)
!   tad(1:4) = (/ 0.6600_KR, 0.6600_KR, 0.6600_KR, 0.6600_KR /)
    tad(1:4) = (/ 1.0_KR, 1.0_KR, 1.0_KR, 1.0_KR /)

! If "tadupdate" is nonzero, then the following algorithm is used:
! 1. Hold the tadpole factors fixed for "tadupdate" sweeps.
! 2. Set the tadpole factors to the median of
!    (a) the tadpole factors from step 1
!    (b) the average of the computed tadpole factors from the previous 
!        "tadupdate" sweeps.
! 3. Return to step 1.
! (The value of tadupdate is irrelevant unless newcalc=3)
    tadupdate = 0

! Compute average plaquette and tadpole factors for the saved configurations?
! [yes=1, no=0] 
! Note: if newcalc=3 and tadupdate is greater than zero, then tadpole=1 is 
!       used regardless of its setting here.
    tadpole = 1

! Choose tadlandau = 1 for the mean link in Landau gauge.
! Choose tadlandau = 0 for the 4th root of an elementary plaquette.
! (irrelevant IF tadpole = 0)
    tadlandau = 0

! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

! LANDAU GAUGE FIXING (irrelevant if tadupdate = 0):

! Set landimp(mu)=1 for improved Landau gauge fixing in the mu direction.
! Otherwise choose landimp(mu)=0.
    landimp = 0

! alpha is a tunable parameter for Landau gauge fixing.
! [Davies et al, Phys Rev D37, 1581 (1988) suggest alpha=0.1]
    alpha = 0.1_KR

! omega is a tunable overrelaxation parameter for Landau gauge fixing.
! [Mandula and Ogilvie, Phys Lett B248, 156 (1990) suggest omega=1.7]
! If omega<1.0, then overrelaxation is not used here.
    omega = 1.7_KR

! thetamax is the limit for an acceptable accuracy in Landau gauge fixing.
    thetamax = 0.00001_KR

! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

! WILSON LOOPS (To turn off all Wilson loop calculations, choose nwilo=0):

! For saved configurations, Wilson loops will be computed for "abs(nwilo)"
! orientations.
! NOTE: To divide all links in Wilson loops by the corresponding tadpole
!       factors, choose nwilo>0.
!       To omit all tadpole factors from Wilson loops, choose nwilo<0.
    nwilo = 0

! wlat(iwilo,1) = the "x" direction for Wilson loop calculations.
!                 (used for 2, 3 and 4-dimensional Wilson loops)
! wlat(iwilo,2) = the "y" direction for Wilson loop calculations.
!                 (used for 3 and 4-dimensional Wilson loops)
! wlat(iwilo,3) = the "z" direction for Wilson loop calculations.
!                 (used for 4-dimensional Wilson loops)
! wlat(iwilo,4) = the "t" direction for Wilson loop calculations.
!                 (used for 2, 3 and 4-dimensional Wilson loops)
! Computed Wilson loops satisfy "nx" >= "ny" >= "nz".
    wlat(1,1) = 1
    wlat(1,2) = 2
    wlat(1,3) = 3
    wlat(1,4) = 4

! For Wilson loops, choose the maximum "spatial" radius (in lattice units)
! for each orientation "iwilo".
    Rmax(1) = 10.0_KR

! Each Wilson loop is actually an average over all "spatial" paths,
! at both in initial and final "times", of all equal-length paths.
! The number of different loops in the average is [(Rx+Ry+Rz)!/Rx!/Ry!/Rz!]^2.
! Only 3 and 4-dimensional off-axis Wilson loops that require up to [npaths]^2
! paths will be computed.
! (2-dimensional Wilson loops have a unique path, and are always computed
! unless Rmax(iwilo)=0.)
    npaths = 6

! For Wilson loops, choose the number of fuzzing iterations and the fuzzing 
! parameter.  These choices are irrelevant if nwilo=0.
    nfuzz = 0
    epsfuzz = 0.3_KR

! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

! FERMION PROPAGATOR CONSTRUCTION (general information):

! src(0,:) is the location of the point source on the global lattice.
    src(0,1) = 1
    src(0,2) = 1
    src(0,3) = 1
    src(0,4) = 4

! bc(mu) = 1,-1,0 for periodic,antiperiodic,fixed boundary conditions
!            respectively.

!   PERIODIC BC (time)

    bc(1) = 1
    bc(2) = 1
    bc(3) = 1
    bc(4) = 1

!   NON-PERIODIC BC (time)

    !bc(1) = 1
    !bc(2) = 1
    !bc(3) = 1
    !bc(4) = 0

! resmax is the stopping criterion for the iteration.
    resmax = 1.0e-30_KR

! itermin is the minimum number of iterations for fermion matrix inversion.
    itermin = 1


! Choose values for nsmear and asmear.
! Note: source smearing for propagators is optional.
!     If nsmear=0 then no source smearing is performed.
!     If nsmear>0 then all propagators are source smeared using (nsmear,asmear).
!     If nsmear<0 then all propagators are computed twice:
!                 one with local source and once using (abs(nsmear),asmear).
! Note: propagators that do not originate at the source (i.e. appended
!       propagators in SST) are never ``source''-smeared in these codes.
    nsmear = 0
    asmear = 0.1_KR

!*******************************fermions on DBW2********************************
!*As discussed above, DBW2 requires                                            *
!*         gaction(1:4) = (/ 2, 2, 2, 2 /),                                    *
!*         tad(1:4) = (/ 0.6600_KR, 0.6600_KR, 0.6600_KR, 0.6600_KR /), and    *
!*         tadupdate = 0                                                       *
!*To use Wilson fermions, choose kappa(i) = 0.140_KR*tad(1),                   *
!*                            or mtmqcd(i,1:2) = (/ 0.140_KR*tad(1), 0.0_KR /).*
!*                         where 0.140 denotes the true kappa value.           *
!*To include the tree-level clover term, choose cSWtmqcd = 1.0_KR*tad(1)**3.   *
!*******************************************************************************

! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

! WILSON FERMION PROPAGATOR CONSTRUCTION WITH MULTI-MASS (NO TWISTED MASS):
! NEW OPTION: CHOOSE EITHER TO MULTI-MASS OR NOT

! To use multi-mass use wilsonmm=1 and to solve serially use wilsonmm=0
! If not using multimassing, we'll treat the kappa values as twisted-mass values
! with zero twisted mass mu. Also in this case their will be rotation necessary
! in when calculating the fermion propagator. The fermion propgator will always
! be given in the physical basis. 

   wilsonmm = 0

! When transforming the problem into a hermetian system and using lan_dr and cg_proj,
! we have the option of solving gamma5*M or M^dagger*M. Note that this transformation
! is applied to the even-odd preconditoined system. So, when calling lan_dr or cg_proj
! to solve M*x=b, we have the option to solve either gamma5*M*x=gamma5*b (indefinite system) 
! or M^dagger*M*x=M^dagger*b (definite system). If hoption=1 means solve gamma5*M and
! hoption=2 means solve M^dagger*M.

   hoption= 2 

! NOTE: The clover coefficient (cSW) is forced to zero in the code, since
! multi-mass and checkerboarding cannot be simultaneously implemented otherwise.

! Choose the number of hopping parameters for propagator construction.
! The maximum is 9, due to file naming limitations, but this could easily be
! changed in the codes.
! (To omit propagator construction, choose nkappa=0.)
   nkappa = 1

! Choose the hopping parameters for propagator construction.
! REQUIRED: kappa(1) > kappa(2) > kappa(3) > kappa(4) > ...
! (irrelevant if nkappa=0.)
    kappa = 0.0_KR

    !beta=5.85
    !kappa(1) = 0.162_KR*tad(1)
    !kappa(1) = 0.166168_KR*tad(1)
    !kappa(2) = 0.157394_KR*tad(1)
    !kappa(1)=0.103664_KR*tad(1)
    !kappa(2)= 0.153_KR*tad(1)



    !beta=6.0
    !kappa(1) = 0.159441_KR*tad(1)
    !kappa(2) = 0.154460_KR*tad(1)


    !beta=6.20
    !kappa(1) = 0.154783_KR*tad(1)
    !kappa(2) = 0.151647_KR*tad(1)

    kappa(1) = 0.165_KR*tad(1)
    kappa(2) = 0.1628_KR*tad(1)
    kappa(3) = 0.161_KR*tad(1)
    
    !kappa(1) = 0.153_KR*tad(1)
    !kappa(2) = 0.155_KR*tad(1)
    !kappa(3) = 0.1558_KR*tad(1)
    kappa(4) = 0.165_KR*tad(1)
    kappa(5) = 0.1675_KR*tad(1)
    kappa(6) = 0.170_KR*tad(1)
    kappa(7) = 0.180_KR*tad(1)
    kappa(8) = 0.190_KR*tad(1)
    kappa(9) = 0.200_KR*tad(1)

! Choose overrelaxation parameter for the minimal residual algorithm.
! Ying,Dong,Liu,Nucl.Phys.B(Proc.Suppl.)53(1997)993 uses omegaM3R=1.1
! but I have found this not to converge for some configurations when
! the quark mass is light.  omegaM3R=1.0 means no overrelaxation.
    omegaM3R = 1.0_KR

! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

! FERMION PROPAGATOR CONSTRUCTION WITH CLOVER AND TWISTED MASS (NO MULTI-MASS):
! (NOW ALLOWS FOR TWISTED-DISCONNECTED LOOPS)

! tmQCD propagators come in pairs.
! If both propagators in each pair need to be constructed,
!    then choose ntmqcd to be (-1)*(the number of pairs).
! If only one propagator in each pair needs to be constructed,
!    then choose ntmqcd to be the the number of pairs.
! File naming limitations require -9 <= ntmqcd <= 9, but this could
! easily be changed in the codes.
! (To omit propagator construction, choose ntmqcd=0.)


    ntmqcd=0
    

! If a proton (uud) using the tmU and tmD quark propagaotors is desired
! choose nucleonparticle=1. If a neutron (udd) is desired using the same 
! basis choose nucleonparticle=2. If both the proton and neutron desired
! choose nucleonparticle=3
! (To use the Wilson basis and construct the nucleon set nucleonparticle=0)

   nucleonparticle = 3 

! Choose the kappa and mu_q values for each pair, where the Dirac operator is
! D_W + m_0 + i*mu_q*gamma_5*tau_3 [see eq.(3) of JHEP08(2001)058]
! and kappa=1/(8+2*m_0).
! In particular, mtmqcd(i,1:2) = (/ kappa(i), mu_q(i) /).
    mtmqcd = 0.0_KR
    
! MASS AVERAGING

  ! beta=6.0 

    mtmqcd(1,1:2)  =  (/ 0.1530_KR*tad(1), 0.00001_KR /)
    mtmqcd(2,1:2)  =  (/ 0.1550_KR*tad(1), 0.00001_KR /)
    mtmqcd(3,1:2)  =  (/ 0.1558_KR*tad(1), 0.00001_KR /)
    mtmqcd(4,1:2)  =  (/ 0.1650_KR*tad(1), 0.00001_KR /)
    mtmqcd(5,1:2)  =  (/ 0.1675_KR*tad(1), 0.00001_KR /)
    mtmqcd(6,1:2)  =  (/ 0.1700_KR*tad(1), 0.00001_KR /)
    mtmqcd(7,1:2)  =  (/ 0.1800_KR*tad(1), 0.00001_KR /)
    mtmqcd(8,1:2)  =  (/ 0.1900_KR*tad(1), 0.00001_KR /)
    mtmqcd(9,1:2)  =  (/ 0.2000_KR*tad(1), 0.00001_KR /)

 
! THESE ARE THE MASSES FOR THE "BIG" RUN
     
    !mtmqcd(1,1:2) = (/ 0.15679_KR*tad(1), 0.030_KR /)
    !mtmqcd(2,1:2) = (/ 0.15708_KR*tad(1), 0.015_KR /)
    !mtmqcd(3,1:2) = (/ 0.15721_KR*tad(1), 0.010_KR /)
    !mtmqcd(4,1:2) = (/ 0.15728_KR*tad(1), 0.005_KR /)
    


! Choose the clover coefficient.
! To use the tree-level coefficient,
!                   choose cSWtmqcd = 1.000_KR*tad(1)*tad(2)*tad(3)*tad(4).
! (To omit the clover term, choose -0.0001 < cSWtmqcd < 0.0001.)
!   cSWtmqcd = 1.000_KR*(tad(1)/0.9058_KR)**3
!   cSWtmqcd = 1.769_KR

!  Dean put this in to omit clover term

    cSWtmqcd = 0.000_KR

! Choose src(i,:) to be the location of the source for the i'th propagator.
! When inserting momenta into an operator, the source is always understood
! to be src(0,:).  However, a derivative operator at src(0,:) requires
! propagators computed from neighbouring sites.  In such a case, src(i,:)
! should be set to the neighbouring site of interest, otherwise simply choose
! src(i,:) = src(0,:).
! RECALL: src(0,:) was already defined above.  Here, define src(1:9,:).
    src(1,:) = src(0,:)
    src(2,:) = src(0,:)
    src(3,:) = src(0,:)
    src(4,:) = src(0,:)
    src(5,:) = src(0,:)
    src(6,:) = src(0,:)
    src(7,:) = src(0,:)
    src(8,:) = src(0,:)
    src(9,:) = src(0,:)

! Choose inverter=0 for bicgstab; inverter=1 for CGNE; inverter=2 for QMR;
!        inverter=3 for QMR(gamma5); inverter=4 for GMRES(n);
!        inverter=5 for GMRES-DR(n,k);
!        inverter=6 for GMRES-Drshift(n,k) and multiple RHS (this also does twisted mass)
! Note: QMR(gamma5) is not permitted if mu_q(i) > 0.
! Note: maximum "n" in GMRES(n) and GMRES-DR(n,k) is defined in module latdims.
    inverter = 6

! Choose nGMRES(1) to be the desired value of "n" in GMRES(n) or GMRES-DR(n,k).
! It cannot be larger than nmaxGMRES as defined in module latdims.
! (irrelevant if inverter<4.)
    nGMRES(1) =200

! Choose nGMRES(2) to be the desired value of the "k" in GMRES-DR(n,k).
! It cannot be larger than nGMRES as defined above.
! (irrelevant if inverter/=5.)
    nGMRES(2) =100

! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

! FERMION PROPAGATOR CONSTRUCTION WITH CONNECTED CURRENT INSERTIONS:

! Choose the operator to be inserted.
! opSST= 1 for local pseudoscalar operator;
! opSST=-1 for smeared pseudoscalar operator (uses nsmear,asmear defined above);
! (To omit connected current insertions entirely, choose opSST=0.)
    opSST = 0

! Choose the timeslice of the inserted operator on the global lattice.
! (irrelevant if opSST=0.)
    tSST = 44

! Choose the 3-momentum for the inserted operator, in units of the smallest
! lattice momentum.
! (irrelevant if opSST=0.)
    pSST(1:3) = (/ 0, 0, 0 /)

! The first "nSST(1)" tmQCD propagator pairs, from the set of "abs(ntmqcd)"
! pairs computed above (or already existing on disk), will each be given
! "nSST(2)" separate appended propagator pairs.
! If both propagators in each pair need to be used/constructed,
!    then choose nSST(1) and/or nSST(2) to be (-1)*(the number of pairs),
! otherwise only one from each pair will be considered.
! File naming limitations require -9 <= nSST(1)*nSST(2) <= 9, but this could
! easily be changed in the codes.
! (irrelevant if opSST=0.)
    nSST(1:2) = (/ 1, 1 /)

! Choose the kappa and mu_q values for each of the appended propagators.
! The format is mSST(i,1:2) = (/ kappa(i), mu_q(i) /), where i=1,nSST(2).
! (irrelevant if opSST=0.)
    mSST = 0.0_KR
    mSST(1,1:2) = (/ 0.157_KR*tad(1), 0.0_KR /)

! Choose the clover coefficient for an appended tmQCD propagator.
! (irrelevant if opSST=0.)
    cSWSST = cSWtmqcd

! Choose the inverter to be used for the appended propagators.
! invSST=0 is bicgstab; 1 is CGNE; 2 is QMR; 3 is QMR(gamma5); 4 is GMRES(n);
!        5 is GMRES-DR(n,k).
! Note: QMR(gamma5) is not permitted if mu_q(i) > 0.
! Note: "n" and "k" in GMRES(n),GMRES-DR(n,k) are nGMRES(1:2) as defined by
!       the user above.
! Note: the maximum "n" in GMRES(n),GMRES-DR(n,k) is defined in module latdims.
! (irrelevant if opSST=0.)
    invSST = 0

! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

! FERMION PROPAGATOR CONSTRUCTION FOR DISCONNECTED LOOPS:



! Set numnoises to the number of Z2 noise sources to be used.
! (To omit disconnected loop construction, choose numnoises <=0.)
    numnoises = 0

! To choose if tmU and tmD disconnected loops are to be calculated.
! If ntmqcdloop = 1 only the tmU loop is calcualted. If ntmqcdloop = -1
! then both the tmU and tmD loops are calcualted, each with a serpate
! set of numnioses.(Value is irrelavant unless inverter==6
    ntmqcdloop = 1  ! For tmU
    !ntmqcdloop = 0 ! Wilson



     !PCAC values
     !mtmqcdloop(1,1:2)  =(/0.157390_KR*tad(1),0.0_KR/) 
     !mtmqcdloop(2,1:2)  =(/0.157390_KR*tad(1),0.01_KR/)
     !mtmqcdloop(3,1:2)  =(/0.157390_KR*tad(1),0.02_KR/)
     !mtmqcdloop(4,1:2)  =(/0.157390_KR*tad(1),0.03_KR/)
     !mtmqcdloop(5,1:2)  =(/0.157390_KR*tad(1),0.04_KR/)


     !beta=5.85
     !mtmqcdloop(1,1:2)  = (/ 0.166168_KR*tad(1), 0.000_KR/)   
     !mtmqcdloop(2,1:2)  = (/ 0.157394_KR*tad(1), 0.000_KR /)   


     !beta=6.0
     !mtmqcdloop(1,1:2)  = (/ 0.159441_KR*tad(1), 0.000_KR/)   
     !mtmqcdloop(2,1:2)  = (/ 0.154460_KR*tad(1), 0.000_KR /)   


     !beta=6.20
     mtmqcdloop(1,1:2)  = (/ 0.152_KR*tad(1), 0.000_KR/)   
     mtmqcdloop(2,1:2)  = (/ 0.151647_KR*tad(1), 0.000_KR /)   







     !mtmqcdloop(1,1:2)  = (/ 0.15679_KR*tad(1), 0.005_KR /)
     !mtmqcdloop(3,1:2)  = (/ 0.15679_KR*tad(1), 0.010_KR /)
     !mtmqcdloop(4,1:2)  = (/ 0.15679_KR*tad(1), 0.015_KR /)

  
     !mtmqcdloop(1,1:2)  = (/ 0.15679_KR*tad(1), 0.025_KR /)   
     !mtmqcdloop(2,1:2)  = (/ 0.15679_KR*tad(1), 0.026_KR /)
     !mtmqcdloop(3,1:2)  = (/ 0.15679_KR*tad(1), 0.027_KR /)
     !mtmqcdloop(4,1:2)  = (/ 0.15679_KR*tad(1), 0.028_KR /)
     !mtmqcdloop(5,1:2)  = (/ 0.15679_KR*tad(1), 0.029_KR /)  
     !mtmqcdloop(6,1:2)  = (/ 0.15679_KR*tad(1), 0.030_KR /)
     !mtmqcdloop(7,1:2)  = (/ 0.15679_KR*tad(1), 0.031_KR /)
     !mtmqcdloop(8,1:2)  = (/ 0.15679_KR*tad(1), 0.032_KR /)
     !mtmqcdloop(9,1:2)  = (/ 0.15679_KR*tad(1), 0.033_KR /) 
     !mtmqcdloop(10,1:2) = (/ 0.15679_KR*tad(1), 0.034_KR /) 
     !mtmqcdloop(11,1:2) = (/ 0.15679_KR*tad(1), 0.035_KR /) 
! Choose the hopping parameter for propagator construction.
! (irrelevant if numnoises <=0.)

! IRRELVANT IF DOING TWISTED DISCONNECTED LOOPS!
! Hard code in the strange quark mass for the subroutine discon.

    kappaloop(1) = .152 ! Wislon
    kappaloop(2) = .154 ! Wislon
  

! If inverter==6 then hopping parameter must be paired with
!   the appropriate mu-value as required in the paper
!   (hep-lat/0503007)

! Choose the clover coefficient.
! (To omit the clover term, choose -0.0001 < cSWloop < 0.0001.)
    cSWloop = 0.0_KR

! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

! CORRELATION FUNCTION CONSTRUCTION:

! Set docorr=0 to construct no correlators.
! Set docorr=1 to construct meson correlators without sink smearing.
! Set docorr=2 to construct meson correlators with sink smearing.
! Set docorr=3 to construct meson correlators with and without sink smearing.
! Set docorr=4 to construct baryon correlators without sink smearing.
! Set docorr=5 to construct baryon correlators with sink smearing.
! Set docorr=6 to construct baryon correlators with and without sink smearing.
! Set docorr=7 to construct meson and baryon correlators without sink smearing.
! Set docorr=8 to construct meson and baryon correlators with sink smearing.
! Set docorr=9 to construct meson and baryon correlators with and without
!              sink smearing.
! Note: baryon sink smearing is not yet implemented in the code, so all 
!       baryon correlators are presently local, even for docorr=4,5,8,9.
!       Hopefully this can be generalized soon.
    docorr = 7

! Choose nkappacorr to be the number of Wilson quark propagators to be used
! in the computation of correlation functions.
! Note: If 0<=nkappacorr<=nkappa is NOT satisfied, then the extra propagators
!       must already exist on disk, since they will not be generated now.

   
    nkappacorr = nkappa 

! tmQCD propagators come in pairs.
! If both propagators in each pair need to be used in correlation functions,
!    then choose ntmqcdcorr to be (-1)*(the number of pairs).
! If only one propagator in each pair needs to be used in correlation functions,
!    then choose ntmqcdcorr to be the the number of pairs.
! Note: If -abs(ntmqcd) <= ntmqcdcorr <= abs(ntmqcd), if ntmqcd<0.
!                     0 <= ntmqcdcorr <=     ntmqcd , if ntmqcd>=0.
!       are NOT satisfied, then the extra propagators must already exist on
!       disk, since they will not be generated now.
    ntmqcdcorr = ntmqcd

! The first "nSST(1)" tmQCD propagator pairs have "nSST(2)" separate appended
! propagator pairs.  If both propagators in each pair need to be used in
! correlation functions, then choose nSSTcorr(1) and/or nSSTcorr(2) to be
! (-1)*(the number of pairs), otherwise only one from each pair will be
! considered.
! File naming limitations require -9 <= nSSTcorr(1)*nSSTcorr(2) <= 9, but this
! could easily be changed in the codes.
    nSSTcorr(1:2) = (/ 1, 1 /)
    nSSTcorr = 0

! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

! CODE TESTING:

! Choose gitest=1 to perform a local gauge transformation on all saved
!                 configurations.  (The transformed configuration is never
!                 written to disk, but the original is, if ifilesave>0 above.)
! Choose gitest=2 to explicitly perform the disconnected loop computations
!                 both before and after a local gauge transformation.
!                 (Requires numnoises>0 above.)
! Choose gitest=0 to omit code testing.
    gitest = 0

! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

    call generator(rwdir,newcalc,ifilesave,cfgfile,iseed,ncold,nsweep, &
                   n0sweep,n1sweep,nhit,gaction,alat,beta,tad,tadupdate, &
                   tadpole,tadlandau,landimp,alpha,omega,thetamax,nwilo, &
                   wlat,Rmax,npaths,nfuzz,epsfuzz,src,bc,resmax,itermin, &
                   nsmear,asmear,nkappa,kappa,omegaM3R,ntmqcd,mtmqcd, &
                   cSWtmqcd,inverter,nGMRES,opSST,tSST,pSST,nSST,mSST,cSWSST, &
                   invSST,numnoises,kappaloop,cSWloop,docorr,nkappacorr, &
                   ntmqcdcorr,ntmqcdloop,mtmqcdloop,nSSTcorr, & 
                   nucleonparticle,gitest,hoption,wilsonmm)

 end program cfgspropsmain

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
