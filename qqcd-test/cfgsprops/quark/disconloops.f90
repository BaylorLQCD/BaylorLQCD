! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
! disconloops.f90, Randy Lewis, randy.lewis@uregina.ca
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
! WARNING! WARNING! WARNING! WARNING! WARNING! WARNING! WARNING! WARNING! 
!
! This module is only correct for cSW=0.0_KR
!
! WARNING! WARNING! WARNING! WARNING! WARNING! WARNING! WARNING! WARNING! 
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
! This module computes disconnected quark loops, each containing a chosen
! operator, from existing gauge field configurations.
!
! It was inspired by nonparallel F77 code provided by Walter Wilcox.
!
! The operators "Jtotal(iri,it,isub,imom,iop)" are written to CFGSPROPS.LOG
! where it = timestep
!       isub=1 is no subtraction,
!       isub=2 is subtraction up to and including O(kappa^3)
!       isub=3 is subtraction up to and including O(kappa^4)
!       isub=4 is subtraction up to and including O(kappa^5)
!       imom=1 is zero-momentum results:
!                 (iri,iop)=(1,1) is Re(J_1)
!                 (iri,iop)=(2,1) is Im(J_1)
!                 (iri,iop)=(1,2) is Re(J_2)
!                 (iri,iop)=(2,2) is Im(J_2)
!                 (iri,iop)=(1,3) is Re(J_3)
!                 (iri,iop)=(2,3) is Im(J_3)
!                 (iri,iop)=(1,4) is Re(J_4)
!                 (iri,iop)=(2,4) is Im(J_4)
!                 (iri,iop)=(1,5) is Re(psibar psi)
!                 (iri,iop)=(2,5) is Im(psibar psi)
!       imom=2,3,4,5 are the smallest nonzero-momentum results that exist
!                    on the specified lattice.
!                    (iri,iop)=(1,1) is factor*Im(J_1)
!                    (iri,iop)=(2,1) is factor*Im(J_1)
!                    (iri,iop)=(1,2) is factor*Im(J_2)
!                    (iri,iop)=(2,2) is factor*Im(J_2)
!                    (iri,iop)=(1,3) is factor*Im(J_3)
!                    (iri,iop)=(2,3) is factor*Im(J_3)
!                    (iri,iop)=(1,4) is factor*Im(J_4)
!                    (iri,iop)=(2,4) is factor*Im(J_4)
!                    (iri,iop)=(1,5) is factor*Im(psibar psi)
!                    (iri,iop)=(2,5) is factor*Im(psibar psi)
!                    where each "factor" is different! (see subroutine momfacs)
! Note that "J_mu" are the conserved (point-split) currents,
! and "psibar psi" is the local scalar.
!
! INFORMATION: The disconnected electric form factor requires the correlation
!              of the operator (iri,iop)=(2,4) with the imaginary part of the
!              unpolarized nucleon two point function.
!              The magnetic form factor can only be extracted at nonzero
!              momentum, and it requires three different polarized two point
!              functions to be separately calculated and correlated with
!              various combinations of operators (iri,iop)=(1,1), (2,1), (1,2),
!              (2,2), (1,3), (2,3).
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

!---------Changes made to the module----------------------------------------
!
!July, 5th, 2007 changes made by Abdou
!
!In order to include the 4 components of the point-split axial current in the 
!loop calculation we need to change nop from 6 to 10 in twistdiscon where nop
!is coded as a parameter. The convention that will be used is that nop=7 for
!the time component and nop=8-10 for x,y,z components respectively.
!---------------------------------------------------------------------------






     module disconloops

!   use MPI
     use kinds
     use latdims
     use basics
     use lagfib
     use lattice
     use gaugetools
     use diracops
     use quark
     use gmresrhs
     use shift
     use vevbleft
     use analysis
     use inverters 
     use working
     use pseudolapack
     use printops
     implicit none
     private

! Use the following line if the MPI module is not available.
     include 'mpif.h'

! Define access to subroutines.
     public  :: discon, twistdiscon, z2source, eigdiscon, z4source, testFUNC
     private :: vev, average, averagePRIME, loopops, mulfor2,mulfor3,mulbac2,mulbac3, vv, vgv, momfacs, ppaverage2, ppaverage, ppaverage1,&
                spacesum, checknonzero, fakegauge,  eigaverage, eigloopops, eigloopops_original, eigspacesum, eigloopops_abdou,  &
                eigmodesofM, eigmodesofMprime, xigenerator, vevoutput, generalaverage,    &
                printz2noise, nsaverage, nsaverage_abdou

     contains

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

! general average seems to be doing the various subtraction combinations 
! at zero momentum and printing them out. It uses nsaverage, eigaverage, 
! and average. --AA

 subroutine generalaverage(Jave,numEV,numEVP,eigstart,eigstep,u,z2e,z2o,xe,xo,kappa, &
                           evecReven,evecRodd,evecLeven,evecLodd,scalarMultiplier,   &
                           evecRevenP,evecRoddP,scalarMultiplierP,                   &
                           xiinvscalarMultiplier,xiinvscalarMultiplierPRIME,xiinvPOLY, &
                           coact,bc,vecbl,vecblinv,myid,nn,ldiv,nms,lvbc,ib,lbd,     &
                           iblv,rwdir,MRT,MRT2)
 real(kind=KR2),    intent(out),   dimension(:,:,:,:,:,:) :: Jave
 integer(kind=KI), intent(in)                            :: numEV,numEVP
 integer(kind=KI), intent(in)                            :: eigstart, eigstep
 real(kind=KR),    intent(inout), dimension(:,:,:,:,:)   :: u
 real(kind=KR),    intent(in),    dimension(:)           :: kappa
 real(kind=KR),    intent(inout), dimension(:,:,:,:,:)   :: z2e, z2o
 real(kind=KR),    intent(inout), dimension(:,:,:,:,:)   :: xe, xo
 real(kind=KR2),   intent(in),    dimension(:,:,:,:,:,:) :: evecReven, evecRodd
 real(kind=KR2),   intent(in),    dimension(:,:,:,:,:,:) :: evecLeven, evecLodd
 real(kind=KR2),   intent(in),    dimension(:,:,:,:,:,:) :: evecRevenP, evecRoddP
 real(kind=KR2),   intent(in),    dimension(:,:)         :: scalarMultiplier,scalarMultiplierP
 real(kind=KR2),   intent(in),    dimension(:,:,:)       :: xiinvscalarMultiplier,xiinvscalarMultiplierPRIME,xiinvPOLY 
 real(kind=KR),    intent(in),    dimension(:,:,:)       :: coact
 integer(kind=KI), intent(in),    dimension(:)           :: bc, nms
 integer(kind=KI), intent(in),    dimension(:,:)         :: vecbl, vecblinv
 integer(kind=KI), intent(in),    dimension(:,:)         :: nn, iblv
 logical,          intent(in),    dimension(:)           :: ldiv
 integer(kind=KI), intent(in),    dimension(:,:,:)       :: lvbc
 integer(kind=KI), intent(in),    dimension(:,:,:,:)     :: ib
 logical,          intent(in),    dimension(:,:)         :: lbd
 integer(kind=KI), intent(in)                            :: myid, MRT2, MRT
 character(len=*), intent(in),    dimension(:)           :: rwdir

! integer(kind=KI), parameter                 :: nsub=6, nmom=5, nop=5, nsty=6
 integer(kind=KI), parameter                 :: nsub=6, nmom=5, nop=9, nsty=10
!BS integer(kind=KI), parameter                 :: nsub=6, nmom=5, nop=9, nsty=9
 real(kind=KR),    dimension(6,nvhalf,4,2,8) :: ZERO_E=0.0_KR2, ZERO_O=0.0_KR2
 !real(kind=KR),    dimension(2,ntotal,4,2,8) :: ZERO_E=0.0_KR2,!ZERO_O=0.0_KR2!BS

 real(kind=KR2),   dimension(2,nt,nmom,nop)                :: JaveNS                             
 real(kind=KR2),    dimension(2,nt,kmaxGMRES,nmom,nop)      :: JaveES
 !real(kind=KR),    dimension(2,nt,nmom,nop)                :: tempHFES!BS
 !real(kind=KR),    dimension(2,nt,nmom,nop)                :: tempHFES2!BS
 !real(kind=KR),    dimension(2,nt,nmom,nop)                :: tempES2!BS
 !real(kind=KR),    dimension(2,nt,nmom,nop)                :: tempES!BS
 real(kind=KR2),    dimension(2,nt,kmaxGMRES,nmom,nop)      :: JaveHFES
 real(kind=KR2),   dimension(2,nt,nsub,nmom,nop)           :: JavePS
 real(kind=KR2),   dimension(2,nt,nmom,nop)                :: tempESPS!BS
 real(kind=KR2),   dimension(2,nt,nsub,nmom,nop)           :: JavePP
 real(kind=KR2),   dimension(2,nt,nsub,nmom,nop)           :: JavePP7
 real(kind=KR2),   dimension(2,nt,nsub,nmom,nop)           :: JavePP4

 real(kind=KR2),   dimension(2,nt,nsub,kmaxGMRES,nmom,nop) :: JaveESPS
 real(kind=KR2),   dimension(2,nt,nsub,kmaxGMRES,nmom,nop) ::JaveHFESPS,JaveHFESPOLY!BS
 !real(kind=KR2),   dimension(2,nt,nmom,nop)                :: tempHFESPS!BS
 !real(kind=KR2),   dimension(2,nt,nmom,nop)                :: tempHFESPS2!BS

 integer(kind=KI)                                   :: numEV_MAX, ieo, ibleo
 integer(kind=KI)                                   :: isub, imom, eig, sty
 real(kind=KR2),  dimension(nvhalf,2,16,nmom-1,nop) :: momfac
 integer(kind=KI)                                   :: dobndry = 0
 real(kind=KR2),  dimension(6,ntotal,4,2,8)         :: gz2e,gz2o !added by Abdou


!calculate gamma5*z

!do ieo=1,2
!   do ibleo =1,8
!     call gammamult(z2e, z2o, gz2e, gz2o, 5,ieo,ibleo)
!   enddo
!enddo

 ! Initalize data
 momfac     = 0.0_KR2
 JaveNS     = 0.0_KR2
 JaveES     = 0.0_KR2
 JaveHFES   = 0.0_KR2
 JavePS     = 0.0_KR2
 JaveESPS   = 0.0_KR2
 JaveHFESPS = 0.0_KR2
 JaveHFESPOLY = 0.0_KR2
 JavePP     = 0.0_KR2
 JavePP4     = 0.0_KR2
 JavePP7     = 0.0_KR2
 !tempPS    =   0.0_KR2
 ! non subtraction (NOTE: using the x solutions,not x')
 call nsaverage(JaveNS,u,z2e,z2o,xe,xo,kappa(1),coact,bc,vecbl,       &
                vecblinv,myid,nn,ldiv,nms,lvbc,ib,lbd,iblv,MRT2,1_KI)

 ! eig subtraction
 call eigaverage(JaveES,numEV,eigstart,eigstep,u,z2e,z2o,xe,xo,kappa(1), &
                 evecReven,evecRodd,evecLeven,evecLodd,scalarMultiplier, &
                 coact,bc,vecbl,vecblinv,myid,nn,ldiv,nms,lvbc,ib,lbd,   &
                 iblv,MRT2,1_KI)

!BSif (myid ==0) then
!BSprint *, "JaveES=",JaveES
!BSendif
 ! hermitian forced eig subtraction
 call eigaverage(JaveHFES,numEVP,eigstart,eigstep,u,z2e,z2o,xe,xo,kappa(1),   &
                 evecRevenP,evecRoddP,evecRevenP,evecRoddP,scalarMultiplierP, &
                 coact,bc,vecbl,vecblinv,myid,nn,ldiv,nms,lvbc,ib,lbd,        &
                 iblv,MRT2,2_KI)

!if (myid ==0) then
!BSprint *, "JaveHFES=",JaveHFES
!endif

 ! pert subtraction
 ! inorder to account for the "different" solution vector x' => gamma5Mx' = z
 ! instead of x => /Mx=z a zero vector is being passed in to calulate the subtraction 
 ! correction term of the perturbative case. The subtraction is done explicity in the
 ! next step. Note that the minus is in the correction calculation.

!Here I am testing the Mx=z case,not g5M case. 
call average(JavePS,nsub,nmom,nop,momfac,u,z2e,z2o,ZERO_E,ZERO_O,kappa(1), &
              dobndry,coact,bc,vecbl,vecblinv,myid,nn,ldiv,nms,lvbc,ib,lbd, &
              iblv,rwdir,MRT2)!To use ppaverage,add MRT to the list

 call ppaverage2(JavePP,nsub,nmom,nop,momfac,u,z2e,z2o,ZERO_E,ZERO_O,kappa(1), &
              dobndry,coact,bc,vecbl,vecblinv,myid,nn,ldiv,nms,lvbc,ib,lbd, &
              iblv,rwdir,MRT,MRT2)!To use ppaverage,add MRT to the list
 call ppaverage1(JavePP7,nsub,nmom,nop,momfac,u,z2e,z2o,ZERO_E,ZERO_O,kappa(1),&
              dobndry,coact,bc,vecbl,vecblinv,myid,nn,ldiv,nms,lvbc,ib,lbd, &
              iblv,rwdir,MRT,MRT2)!To use ppaverage,add MRT to the list
 call ppaverage(JavePP4,nsub,nmom,nop,momfac,u,z2e,z2o,ZERO_E,ZERO_O,kappa(1), &
              dobndry,coact,bc,vecbl,vecblinv,myid,nn,ldiv,nms,lvbc,ib,lbd, &
              iblv,rwdir,MRT,MRT2)!To use ppaverage,add MRT to the list



 JavePS = -JavePS 
 JavePP = -JavePP
 JavePP4 = -JavePP4
 JavePP7 = -JavePP7
 do isub=1,nsub
   ! eig subtraction for ES+PS
   call eigaverage(JaveESPS(:,:,isub,:,:,:),numEV,eigstart,eigstep,u,z2e,z2o,xe,xo,kappa(1), &
                   evecReven,evecRodd,evecLeven,evecLodd,xiinvscalarMultiplier(:,isub,:),    &
                   coact,bc,vecbl,vecblinv,myid,nn,ldiv,nms,lvbc,ib,lbd,iblv,MRT2,1_KI)

   ! hermitian forced eig subtraction for HFES+PS
   call eigaverage(JaveHFESPS(:,:,isub,:,:,:),numEVP,eigstart,eigstep,u,z2e,z2o,xe,xo,kappa(1),    &
                   evecRevenP,evecRoddP,evecRevenP,evecRoddP,xiinvscalarMultiplierPRIME(:,isub,:), &
                   coact,bc,vecbl,vecblinv,myid,nn,ldiv,nms,lvbc,ib,lbd,iblv,MRT2,2_KI)


   !BS 4/25/2016  polynomial deflation HFES+POLY
   call eigaverage(JaveHFESPOLY(:,:,isub,:,:,:),numEVP,eigstart,eigstep,u,z2e,z2o,xe,xo,kappa(1),&
                   evecRevenP,evecRoddP,evecRevenP,evecRoddP,xiinvPOLY(:,isub,:), &
                   coact,bc,vecbl,vecblinv,myid,nn,ldiv,nms,lvbc,ib,lbd,iblv,MRT2,2_KI)










 enddo

 !     JaveESPS = -JaveESPS!BS
  !    JaveHFESPS = - JaveHFESPS!BS
   !   JaveHFESPS = 0.0_KR2!BS

 ! COPY INTO FINAL HOLDER
 imom = 1 ! doing 0 momentum calculations
 numEV_MAX = max(numEV,numEVP)
 do isub=1,nsub
   do eig=1,numEV_MAX
     do sty=1,10!BS changed from 9 to 12
       ! sty = 1  NS
       if (sty==1 .and. eig==1 .and. isub==1) then
         Jave(:,:,isub,eig,sty,:) = JaveNS(:,:,imom,:)
       ! sty = 2  PS
       else if (sty==2 .and. eig==1) then
         Jave(:,:,isub,eig,sty,:) = JaveNS(:,:,imom,:) - JavePS(:,:,isub,imom,:)
       ! sty = 3  ES
       else if (sty==3 .and. isub==1) then
         !tempES(:,:,imom,:) = tempES(:,:,imom,:) + JaveES(:,:,eig,imom,:)   
         Jave(:,:,isub,eig,sty,:) = JaveNS(:,:,imom,:) - JaveES(:,:,eig,imom,:)
     
       ! sty = 4  HFES
       else if (sty==4 .and. isub==1) then
         Jave(:,:,isub,eig,sty,:) = JaveNS(:,:,imom,:) - JaveHFES(:,:,eig,imom,:)!BS
      
        ! sty = 5  ES+PS
       else if (sty==5) then
         
         Jave(:,:,isub,eig,sty,:) = JaveNS(:,:,imom,:) - JaveES(:,:,eig,imom,:)    &
                                  - (JavePS(:,:,isub,imom,:) - JaveESPS(:,:,isub,eig,imom,:)) 
       ! sty = 6  HFES+PS
       else if (sty==6) then
         Jave(:,:,isub,eig,sty,:) = JaveNS(:,:,imom,:)  - JaveHFES(:,:,eig,imom,:) &
                                  - (JavePS(:,:,isub,imom,:) - JaveHFESPS(:,:,isub,eig,imom,:))
       !sty = 7,8,9 polysubtraction
       else if (sty==7 .and. eig==1) then
         Jave(:,:,isub,eig,sty,:) = JaveNS(:,:,imom,:) - JavePP(:,:,isub,imom,:)
       else if (sty==8 .and. eig==1) then
         Jave(:,:,isub,eig,sty,:) = JaveNS(:,:,imom,:) -JavePP4(:,:,isub,imom,:)
       else if (sty==9 .and. eig==1) then
         Jave(:,:,isub,eig,sty,:) = JaveNS(:,:,imom,:) -JavePP7(:,:,isub,imom,:)
       else if (sty==10) then
         Jave(:,:,isub,eig,sty,:) = JaveNS(:,:,imom,:)  - JaveHFES(:,:,eig,imom,:) &
                                  - (JavePP7(:,:,isub,imom,:) -JaveHFESPOLY(:,:,isub,eig,imom,:))
      !BS else if (sty==11) then
         !BSJave(:,:,isub,eig,sty,:) =  JaveNS(:,:,imom,:)  - JaveHFES(:,:,eig,imom,:) & 
            !BS                      - (JavePP7(:,:,isub,imom,:) - JaveHFESPS(:,:,isub,eig,imom,:)) 
       !BS else if (sty==12) then
         !BS Jave(:,:,isub,eig,sty,:) = JaveHFES(:,:,eig,imom,:) !BS sty==10 added

       endif
     enddo
   enddo
 enddo

 end subroutine generalaverage



! vevoutput seems to be just a function to printout the output for loop calculations -AA

 subroutine vevoutput(Jtotal,numEV,numEVPRIME,eigstart,eigstep,rwdir,myid)
 real(kind=KR2), intent(in),   dimension(:,:,:,:,:,:) :: Jtotal
 character(len=*), intent(in), dimension(:)           :: rwdir
 integer(kind=KI), intent(in) :: numEV,numEVPRIME,eigstart,eigstep
 integer(kind=KI), intent(in) :: myid

 integer(kind=KI) :: eig,sty,iop,it,isub,numEV_MAX
! integer(kind=KI), parameter :: ksub=6, nop=5, nsub=6;
 integer(kind=KI), parameter :: ksub=6, nop=9, nsub=6;


 if (myid==0) then
   numEV_MAX = max(numEV,numEVPRIME)

   open(unit=21,file=trim(rwdir(myid+1))//"ns.dat",status="old",     &
        action="write",form="formatted",position="append")
   open(unit=22,file=trim(rwdir(myid+1))//"ps.dat",status="old",     &
        action="write",form="formatted",position="append")
   open(unit=23,file=trim(rwdir(myid+1))//"es.dat",status="old",     &
        action="write",form="formatted",position="append")
   open(unit=24,file=trim(rwdir(myid+1))//"hfes.dat",status="old",   &
        action="write",form="formatted",position="append")
   open(unit=25,file=trim(rwdir(myid+1))//"esps.dat",status="old",   &
        action="write",form="formatted",position="append")
   open(unit=26,file=trim(rwdir(myid+1))//"hfesps.dat",status="old", &
        action="write",form="formatted",position="append")
   open(unit=27,file=trim(rwdir(myid+1))//"pp.dat",status="old",     &
        action="write",form="formatted",position="append")
   open(unit=28,file=trim(rwdir(myid+1))//"pp4.dat",status="old",     &
        action="write",form="formatted",position="append")
   open(unit=29,file=trim(rwdir(myid+1))//"pp7.dat",status="old",     &
        action="write",form="formatted",position="append")
   open(unit=30,file=trim(rwdir(myid+1))//"hfespoly.dat",status="old",     &
        action="write",form="formatted",position="append")!BS
!BS   open(unit=31,file=trim(rwdir(myid+1))//"NSPSHFESPS.dat",status="old",     &
   !     action="write",form="formatted",position="append")!BS
  !BS open(unit=32,file=trim(rwdir(myid+1))//"NSPSHFES.dat",status="old",     &
     !   action="write",form="formatted",position="append")!BS


   !do sty=1,9!BS
   do sty=1,10
     do iop=1,nop
       do isub=1,nsub
         do eig=1,numEV_MAX
           do it=1,nt
             ! sty = 1  NS
             if (sty==1 .and. eig==1 .and. isub==1) then
               write(unit=21,fmt='(i3,i3,i3,i5,i3,f25.15,a3,f25.15)') sty,iop,isub,eig,it, &
                     Jtotal(1,it,isub,eig,sty,iop),'   ',Jtotal(2,it,isub,eig,sty,iop)
             ! sty = 2  PS
             else if (sty==2 .and. eig==1) then
               write(unit=22,fmt='(i3,i3,i3,i5,i3,f25.15,a3,f25.15)') sty,iop,isub,eig,it, &
                     Jtotal(1,it,isub,eig,sty,iop),'   ',Jtotal(2,it,isub,eig,sty,iop)
             ! sty = 3  ES
             else if (sty==3 .and. isub==1 .and. ((eig == eigstart) .or. (eig > eigstart .and. 0 == mod(eig-eigstart,eigstep)) .or. (eig == numEV))) then
               write(unit=23,fmt='(i3,i3,i3,i5,i3,f25.15,a3,f25.15)') sty,iop,isub,eig,it, &
                     Jtotal(1,it,isub,eig,sty,iop),'   ',Jtotal(2,it,isub,eig,sty,iop)
             ! sty = 4  HFES
             else if (sty==4 .and. isub==1 .and. ((eig == eigstart) .or. (eig > eigstart .and. 0 == mod(eig-eigstart,eigstep)) .or. (eig == numEVPRIME))) then
               write(unit=24,fmt='(i3,i3,i3,i5,i3,f25.15,a3,f25.15)') sty,iop,isub,eig,it, &
                     Jtotal(1,it,isub,eig,sty,iop),'   ',Jtotal(2,it,isub,eig,sty,iop)
             ! sty = 5  ES+PS
             else if (sty==5 .and. ((eig == eigstart) .or. (eig > eigstart .and. 0 == mod(eig-eigstart,eigstep)) .or. (eig == numEV))) then
               write(unit=25,fmt='(i3,i3,i3,i5,i3,f25.15,a3,f25.15)') sty,iop,isub,eig,it, &
                     Jtotal(1,it,isub,eig,sty,iop),'   ',Jtotal(2,it,isub,eig,sty,iop)
             ! sty = 6  HFES+PS
             else if (sty==6 .and. ((eig == eigstart) .or. (eig > eigstart .and. 0 == mod(eig-eigstart,eigstep)) .or. (eig == numEVPRIME))) then
               write(unit=26,fmt='(i3,i3,i3,i5,i3,f25.15,a3,f25.15)') sty,iop,isub,eig,it, &
                     Jtotal(1,it,isub,eig,sty,iop),'   ',Jtotal(2,it,isub,eig,sty,iop)
              ! sty = 7,8,9 PP, PP4, PP7
             else if (sty==7 .and. eig==1) then
               write(unit=27,fmt='(i3,i3,i3,i5,i3,f25.15,a3,f25.15)') sty,iop,isub,eig,it, &
                     Jtotal(1,it,isub,eig,sty,iop),'   ',Jtotal(2,it,isub,eig,sty,iop)
             else if (sty==8 .and. eig==1) then
               write(unit=28,fmt='(i3,i3,i3,i5,i3,f25.15,a3,f25.15)') sty,iop,isub,eig,it, &
                     Jtotal(1,it,isub,eig,sty,iop),'   ',Jtotal(2,it,isub,eig,sty,iop)
             else if (sty==9 .and. eig==1) then
               write(unit=29,fmt='(i3,i3,i3,i5,i3,f25.15,a3,f25.15)') sty,iop,isub,eig,it, &
                     Jtotal(1,it,isub,eig,sty,iop),'   ',Jtotal(2,it,isub,eig,sty,iop)
             else if (sty==10 .and. ((eig == eigstart) .or. (eig > eigstart .and. 0 == mod(eig-eigstart,eigstep)) .or. (eig == numEVPRIME))) then
               write(unit=30,fmt='(i3,i3,i3,i5,i3,f25.15,a3,f25.15)') sty,iop,isub,eig,it, &
                     Jtotal(1,it,isub,eig,sty,iop),'   ',Jtotal(2,it,isub,eig,sty,iop)!BS
             !else if (sty==11 .and. ((eig == eigstart) .or. (eig > eigstart.and. 0 == mod(eig-eigstart,eigstep)) .or. (eig == numEVPRIME))) then
              ! write(unit=31,fmt='(i3,i3,i3,i5,i3,f25.15,a3,f25.15)')sty,iop,isub,eig,it, &
               !      Jtotal(1,it,isub,eig,sty,iop),'',Jtotal(2,it,isub,eig,sty,iop)!BS
            ! else if (sty==12 .and. ((eig == eigstart) .or. (eig > eigstart.and. 0 == mod(eig-eigstart,eigstep)) .or. (eig == numEVPRIME))) then
             !  write(unit=32,fmt='(i3,i3,i3,i5,i3,f25.15,a3,f25.15)')sty,iop,isub,eig,it, &
              !       Jtotal(1,it,isub,eig,sty,iop),'',Jtotal(2,it,isub,eig,sty,iop)!BS





   !BS          else 
!BSprint *,"didnot print eig eigstart eigstep numEV ",eig,eigstart,eigstep,numEV

            endif
           enddo
         enddo
       enddo
     enddo
   enddo

   close(unit=21,status="keep")
   close(unit=22,status="keep")
   close(unit=23,status="keep")
   close(unit=24,status="keep")
   close(unit=25,status="keep")
   close(unit=26,status="keep")
   close(unit=27,status="keep")
   close(unit=28,status="keep")
   close(unit=29,status="keep")
   close(unit=30,status="keep")!BS
  ! close(unit=31,status="keep")!BS
  ! close(unit=32,status="keep")!BS
 endif
 end subroutine vevoutput




 ! OUTPUT: 
 !   xiinv(2,6,numEV)
 ! gammaid = 1 standard problem
 ! gammaid = anythingelse prime problem
 subroutine xigenerator(xiinv,numEV,eigstart,eigstep,u,kappa, &
                        evecReven,evecRodd,evecLeven,evecLodd, &
                        coact,bc,vecbl,vecblinv,myid,nn,ldiv,nms,lvbc,ib,lbd, &
                        iblv,MRT2,gammaid)
 real(kind=KR2),   intent(out), dimension(:,:,:)       :: xiinv
 integer(kind=KI), intent(in)                          :: numEV,gammaid
 integer(kind=KI), intent(in)                          :: eigstart, eigstep
 real(kind=KR),    intent(inout), dimension(:,:,:,:,:) :: u
 real(kind=KR),    intent(in)                          :: kappa
 real(kind=KR2),   intent(in), dimension(:,:,:,:,:,:)  :: evecReven, evecRodd
 real(kind=KR2),   intent(in), dimension(:,:,:,:,:,:)  :: evecLeven, evecLodd
 real(kind=KR),    intent(in), dimension(:,:,:)        :: coact
 integer(kind=KI), intent(in), dimension(:)            :: bc, nms
 integer(kind=KI), intent(in), dimension(:,:)          :: vecbl, vecblinv
 integer(kind=KI), intent(in), dimension(:,:)          :: nn, iblv
 logical,          intent(in), dimension(:)            :: ldiv
 integer(kind=KI), intent(in), dimension(:,:,:)        :: lvbc
 integer(kind=KI), intent(in), dimension(:,:,:,:)      :: ib
 logical,          intent(in), dimension(:,:)          :: lbd
 integer(kind=KI), intent(in)                          :: myid, MRT2

 real(kind=KR2), dimension(2)              :: tmp1, tmp2, tmp3, tmp4
 integer(kind=KI)                          :: eig,mu,gblclr,ieo,ibleo,icri
 integer(kind=KI)                          :: eigindx,sty,eigstyprime
 real(kind=KR2), dimension(6,ntotal,4,2,8) :: xe,xo
 real(kind=KR2), dimension(6,ntotal,4,2,8) :: eRe,eRo
 real(kind=KR2), dimension(6,ntotal,4,2,8) :: be,bo
 real(kind=KR2), dimension(6,ntotal,4,2,8) :: tmpE,tmpO
 real(kind=KR2), dimension(6,ntotal,4,2,8) :: sube,subo,tempsube,tempsubo
 real(kind=KR2), dimension(6,ntotal,4,2,8) :: sub1e,sub1o,temptempsube
 real(kind=KR2), dimension(6,ntotal,4,2,8) :: sub2e,sub2o,temptempsubo
 integer(kind=KI)   :: isite, iri, ibl, itbit, gblcrl
 integer(kind=KI)   :: isub, ksub, id,nsub, idag
 real(kind=KR2)     :: xk
 integer(kind=KI) :: exists

 idag = 0
 nsub = 6



! BS 4/25/2016 Finding coefficients for polynomials

























!BS start making changes for gamma5 multiplication

 do eig = 1, numEV

  !BS   if (gammaid == 1_KI) gamma5 not done here so its same for both gammaid
      be = evecReven(:,:,:,:,:,eig)
      bo = evecRodd(:,:,:,:,:,eig)
! !   else ! gamma5 mult 
    





if(.false.)then
  if (gammaid==2_KI) then
   if (eig==1) then
     if (myid==0) then
         inquire(file=trim('/data/barals/qqcd/scratch/evecR.dat'), exist=exists)
         if (.not. exists) then
            print *, "File does not exist. Creating it."
            open(unit=36,file=trim('/data/barals/qqcd/scratch/evecR.dat'),status="new",     &
                action="write",form="formatted")
            close(unit=36,status="keep")
         endif




       open(unit=36,file=trim('/data/barals/qqcd/scratch/evecR.dat'),status="old",     &
         action="write",form="formatted",position="append")

    
          do ibleo = 1,8
            do ieo = 1,2
              do id = 1,4
                do isite = 1,nvhalf
                  do icri = 1,6

          write(unit=36,fmt="(i5,a5,i5,a5,i5,a5,i5,a5,i5,a5,F13.10,a5,F13.10)")  ibleo,"",ieo,"",id,"",isite,"",icri,"",be(icri,isite,id,ieo,ibleo),"",bo(icri,isite,id,ieo,ibleo)
       
                  enddo ! icri
                enddo ! isite
              enddo ! id
            enddo ! ieo
          enddo ! ibleo
    
    
    

       close(unit=36,status="keep")
     endif!myid
  endif!eig
endif!gammaid

endif!true or false












 !    do ieo=1,2
  !      do ibleo=1,8
   !       call gammamult(evecReven(:,:,:,:,:,eig), evecRodd(:,:,:,:,:,eig), be, bo, 5,ieo,ibleo)
    !    enddo
    !  enddo
!    endif
! BS end comment out.......

    ! Compute the operators for each level of subtraction.
    do isub = 1,nsub
      select case(isub)
        case(1) ! no subtraction
          ksub = 0
          sube = 0.0_KR
          subo = 0.0_KR
        case(2) ! subtraction of O(kappa), O(kappa^2) and O(kappa^3)
          ksub = 1
          gblclr = 1
          call Hsingle(sub1e,u,bo,idag,coact,bc,gblclr,vecbl,vecblinv,myid, &
                       nn,ldiv,nms,lvbc,ib,lbd,iblv,MRT2)
          gblclr = 2
          call Hsingle(sub1o,u,be,idag,coact,bc,gblclr,vecbl,vecblinv,myid, &
                       nn,ldiv,nms,lvbc,ib,lbd,iblv,MRT2)
          do ibleo = 1,8
            do ieo = 1,2
              do id = 1,4
                do isite = 1,nvhalf
                  do icri = 1,6
                    sube(icri,isite,id,ieo,ibleo) = be(icri,isite,id,ieo,ibleo) &
                                         + kappa*sub1e(icri,isite,id,ieo,ibleo)
                    subo(icri,isite,id,ieo,ibleo) = bo(icri,isite,id,ieo,ibleo) &
                                         + kappa*sub1o(icri,isite,id,ieo,ibleo)
                  enddo ! icri
                enddo ! isite
              enddo ! id
            enddo ! ieo
          enddo ! ibleo
          gblclr = 1
          call Hsingle(sub2e,u,sub1o,idag,coact,bc,gblclr,vecbl,vecblinv,myid, &
                       nn,ldiv,nms,lvbc,ib,lbd,iblv,MRT2)
          gblclr = 2
          call Hsingle(sub2o,u,sub1e,idag,coact,bc,gblclr,vecbl,vecblinv,myid, &
                       nn,ldiv,nms,lvbc,ib,lbd,iblv,MRT2)
          xk = kappa**2
          do ibleo = 1,8
            do ieo = 1,2
              do id = 1,4
                do isite = 1,nvhalf
                  do icri = 1,6
                    sube(icri,isite,id,ieo,ibleo) = sube(icri,isite,id,ieo,ibleo) &
                                              + xk*sub2e(icri,isite,id,ieo,ibleo)
                    subo(icri,isite,id,ieo,ibleo) = subo(icri,isite,id,ieo,ibleo) &
                                              + xk*sub2o(icri,isite,id,ieo,ibleo)
                  enddo ! icri
                enddo ! isite
              enddo ! id
            enddo ! ieo
          enddo ! ibleo
          gblclr = 1
          call Hsingle(sub1e,u,sub2o,idag,coact,bc,gblclr,vecbl,vecblinv,myid, &
                       nn,ldiv,nms,lvbc,ib,lbd,iblv,MRT2)
          gblclr = 2
          call Hsingle(sub1o,u,sub2e,idag,coact,bc,gblclr,vecbl,vecblinv,myid, &
                       nn,ldiv,nms,lvbc,ib,lbd,iblv,MRT2)
          xk = kappa**3
          do ibleo = 1,8
            do ieo = 1,2
              do id = 1,4
                do isite = 1,nvhalf
                  do icri = 1,6
                    sube(icri,isite,id,ieo,ibleo) = sube(icri,isite,id,ieo,ibleo) &
                                              + xk*sub1e(icri,isite,id,ieo,ibleo)
                    subo(icri,isite,id,ieo,ibleo) = subo(icri,isite,id,ieo,ibleo) &
                                              + xk*sub1o(icri,isite,id,ieo,ibleo)
                  enddo ! icri
                enddo ! isite
              enddo ! id
            enddo ! ieo
          enddo ! ibleo
        case(3) ! subtraction of O(kappa^4)
          gblclr = 1
          call Hsingle(sub2e,u,sub1o,idag,coact,bc,gblclr,vecbl,vecblinv,myid, &
                       nn,ldiv,nms,lvbc,ib,lbd,iblv,MRT2)
          gblclr = 2
          call Hsingle(sub2o,u,sub1e,idag,coact,bc,gblclr,vecbl,vecblinv,myid, &
                       nn,ldiv,nms,lvbc,ib,lbd,iblv,MRT2)
          xk = kappa**4
          do ibleo = 1,8
            do ieo = 1,2
              do id = 1,4
                do isite = 1,nvhalf
                  do icri = 1,6
                    sube(icri,isite,id,ieo,ibleo) = sube(icri,isite,id,ieo,ibleo) &
                                              + xk*sub2e(icri,isite,id,ieo,ibleo)
                    subo(icri,isite,id,ieo,ibleo) = subo(icri,isite,id,ieo,ibleo) &
                                              + xk*sub2o(icri,isite,id,ieo,ibleo)
                  enddo ! icri
                enddo ! isite
              enddo ! id
            enddo ! ieo
          enddo ! ibleo
        case(4) ! subtraction of O(kappa^5)
          gblclr = 1
          call Hsingle(sub1e,u,sub2o,idag,coact,bc,gblclr,vecbl,vecblinv,myid, &
                       nn,ldiv,nms,lvbc,ib,lbd,iblv,MRT2)
          gblclr = 2
          call Hsingle(sub1o,u,sub2e,idag,coact,bc,gblclr,vecbl,vecblinv,myid, &
                       nn,ldiv,nms,lvbc,ib,lbd,iblv,MRT2)
          xk = kappa**5
          do ibleo = 1,8
            do ieo = 1,2
              do id = 1,4
                do isite = 1,nvhalf
                  do icri = 1,6
                    sube(icri,isite,id,ieo,ibleo) = sube(icri,isite,id,ieo,ibleo) &
                                              + xk*sub1e(icri,isite,id,ieo,ibleo)
                    subo(icri,isite,id,ieo,ibleo) = subo(icri,isite,id,ieo,ibleo) &
                                              + xk*sub1o(icri,isite,id,ieo,ibleo)
                  enddo ! icri
                enddo ! isite
              enddo ! id
            enddo ! ieo
          enddo ! ibleo
    !*****************BS  start commenting out******************88
    !    case(5) ! subtraction of O(kappa^6)
     !     gblclr = 1
      !    call Hsingle(sub1e,u,sub2o,idag,coact,bc,gblclr,vecbl,vecblinv,myid, &
       !                nn,ldiv,nms,lvbc,ib,lbd,iblv,MRT2)
        !  gblclr = 2
         ! call Hsingle(sub1o,u,sub2e,idag,coact,bc,gblclr,vecbl,vecblinv,myid, &
          !             nn,ldiv,nms,lvbc,ib,lbd,iblv,MRT2)
        !  xk = kappa**6
         ! do ibleo = 1,8
          !  do ieo = 1,2
           !   do id = 1,4
            !    do isite = 1,nvhalf
                 ! do icri = 1,6
                !    sube(icri,isite,id,ieo,ibleo) = sube(icri,isite,id,ieo,ibleo) &
               !                               + xk*sub1e(icri,isite,id,ieo,ibleo)
              !      subo(icri,isite,id,ieo,ibleo) = subo(icri,isite,id,ieo,ibleo) &
             !                                 + xk*sub1o(icri,isite,id,ieo,ibleo)
            !      enddo ! icri
           !     enddo ! isite
          !    enddo ! id
         !   enddo ! ieo
        !  enddo ! ibleo
      !*********BS end commenting out*****************
 
    ! *******BS 12/23/2015 start writiing*************8

        case(5) ! subtraction of O(kappa^6)
          gblclr = 1
          call Hsingle(sub2e,u,sub1o,idag,coact,bc,gblclr,vecbl,vecblinv,myid, &
                       nn,ldiv,nms,lvbc,ib,lbd,iblv,MRT2)
          gblclr = 2
          call Hsingle(sub2o,u,sub1e,idag,coact,bc,gblclr,vecbl,vecblinv,myid, &
                       nn,ldiv,nms,lvbc,ib,lbd,iblv,MRT2)
          xk = kappa**6
          do ibleo = 1,8
            do ieo = 1,2
              do id = 1,4
                do isite = 1,nvhalf
                  do icri = 1,6
                    sube(icri,isite,id,ieo,ibleo) =sube(icri,isite,id,ieo,ibleo) &
                                              +xk*sub2e(icri,isite,id,ieo,ibleo)
                    subo(icri,isite,id,ieo,ibleo) =subo(icri,isite,id,ieo,ibleo) &
                                              +xk*sub2o(icri,isite,id,ieo,ibleo)
                  enddo ! icri
                enddo ! isite
              enddo ! id
            enddo ! ieo
          enddo ! ibleo

       !**********BS end writing**************

        case(6) ! subtraction of O(kappa^7)
          gblclr = 1
          call Hsingle(sub1e,u,sub2o,idag,coact,bc,gblclr,vecbl,vecblinv,myid, &
                       nn,ldiv,nms,lvbc,ib,lbd,iblv,MRT2)
          gblclr = 2
          call Hsingle(sub1o,u,sub2e,idag,coact,bc,gblclr,vecbl,vecblinv,myid, &
                       nn,ldiv,nms,lvbc,ib,lbd,iblv,MRT2)
          xk = kappa**7!BS 12/23/2015 changed 6 to 7
          do ibleo = 1,8
            do ieo = 1,2
              do id = 1,4
                do isite = 1,nvhalf
                  do icri = 1,6
                    sube(icri,isite,id,ieo,ibleo) = sube(icri,isite,id,ieo,ibleo) &
                                              + xk*sub1e(icri,isite,id,ieo,ibleo)
                    subo(icri,isite,id,ieo,ibleo) = subo(icri,isite,id,ieo,ibleo) &
                                              + xk*sub1o(icri,isite,id,ieo,ibleo)
                  enddo ! icri
                enddo ! isite
              enddo ! id
            enddo ! ieo
          enddo ! ibleo
        case default





          open(unit=8,file="DISCONLOOPS.ERROR",action="write",status="replace", &
               form="formatted")
          write(unit=8,fmt=*) "subroutine average: isub =", isub
          close(unit=8,status="keep")
          stop
      end select

 
        temptempsube=sube !BS 4/25/2016
        temptempsubo=subo !BS 4/25/2016


     !BS start changes...flip gamma5 ...see vic's dissertation equation 4.83
     !BS  if (gammaid == 1_KI) then 
       if (gammaid .NE. 1_KI) then 
           do ieo=1,2
             do ibleo=1,8
               call gammamult(sube, subo, tempsube,tempsubo, 5,ieo,ibleo)
             enddo
           enddo
        sube = tempsube
        subo = tempsubo
      endif
    ! BS end changes here 

      call vecdot(evecLeven(:,:,:,:,:,eig),sube,tmp1,MRT2)
      call vecdot(evecLodd(:,:,:,:,:,eig),subo,tmp2,MRT2)
      xiinv(:2,isub, eig) = tmp1(:2) + tmp2(:2)

   sube=temptempsube
   subo=temptempsubo
    enddo
 enddo

   
 end subroutine xigenerator 








subroutine xipolygenerator(xiinvPOLY,numEV,eigstart,eigstep,u,kappa, &
                        evecReven,evecRodd,evecLeven,evecLodd,z2e,z2o, &
                        coact,bc,vecbl,vecblinv,myid,nn,ldiv,nms,lvbc,ib,lbd,&
                        iblv,MRT2,gammaid)
 real(kind=KR2),   intent(out), dimension(:,:,:)       :: xiinvPOLY
 integer(kind=KI), intent(in)                          :: numEV,gammaid
 integer(kind=KI), intent(in)                          :: eigstart,eigstep
 real(kind=KR),    intent(inout), dimension(:,:,:,:,:) :: u
 real(kind=KR),    intent(in)                          :: kappa
 real(kind=KR2),   intent(in), dimension(:,:,:,:,:,:)  :: evecReven,evecRodd
 real(kind=KR2),   intent(in), dimension(:,:,:,:,:,:)  :: evecLeven,evecLodd
  real(kind=KR),    intent(in), dimension(:,:,:,:,:) :: z2e, z2o

 real(kind=KR),    intent(in), dimension(:,:,:)        :: coact
 integer(kind=KI), intent(in), dimension(:)            :: bc, nms
 integer(kind=KI), intent(in), dimension(:,:)          :: vecbl,vecblinv
 integer(kind=KI), intent(in), dimension(:,:)          :: nn, iblv
 logical,          intent(in), dimension(:)            :: ldiv
 integer(kind=KI), intent(in), dimension(:,:,:)        :: lvbc
 integer(kind=KI), intent(in), dimension(:,:,:,:)      :: ib
 logical,          intent(in), dimension(:,:)          :: lbd
 integer(kind=KI), intent(in)                          :: myid, MRT2

 real(kind=KR2), dimension(2)              :: tmp1, tmp2, tmp3, tmp4
 integer(kind=KI)                          ::eig,mu,gblclr,ieo,ibleo,icri
 integer(kind=KI)                          :: eigindx,sty,eigstyprime
 real(kind=KR2), dimension(6,ntotal,4,2,8) :: xe,xo
 real(kind=KR2), dimension(6,ntotal,4,2,8) :: eRe,eRo
 real(kind=KR2), dimension(6,ntotal,4,2,8) :: be,bo
 real(kind=KR2), dimension(6,ntotal,4,2,8) :: tmpE,tmpO,tue,tuo
 real(kind=KR2), dimension(6,ntotal,4,2,8) ::sube,subo,tempsube,tempsubo
 real(kind=KR2), dimension(6,ntotal,4,2,8) :: sub1e,sub1o,temptempsube
 real(kind=KR2), dimension(6,ntotal,4,2,8) :: sub2e,sub2o,temptempsubo
 integer(kind=KI)   :: isite, iri, ibl, itbit, gblcrl
 integer(kind=KI)   :: isub, ksub, id,nsub, idag
 real(kind=KR2)     :: xk
 integer(kind=KI) :: exists
 integer(kind=KI) ::   p, i, j
 real(kind=KR2),   dimension(2)                          :: beta, beta1



 real(kind=KR2),   dimension(2,8,8)        ::  lsmat,lsmat1
 real(kind=KR2),   dimension(2,8,1)        ::  cls1
 real(kind=KR2),   dimension(2,8)          ::  co2


 idag = 0
 nsub = 6



!copy from ppaverage1 begins now


    p = 8



  if (gammaid .NE. 1_KI) then
         vprime = 0.0_KR2
         vprime1 = 0.0_KR2
        do isite = 1,nvhalf
           vprime(:,isite,:,:,:,1) = z2e(:,isite,:,:,:)
        enddo !isite

        do isite = 1,nvhalf
           vprime1(:,isite,:,:,:,1) = z2o(:,isite,:,:,:)
        enddo !isite


        do i = 1,p
           gblclr = 1
           call Hsingle(tue,u,vprime1(:,:,:,:,:,i),idag,coact,bc,gblclr,vecbl, &
                       vecblinv,myid,nn,ldiv,nms,lvbc,ib,lbd,iblv,MRT2)
           gblclr = 2
           call Hsingle(tuo,u,vprime(:,:,:,:,:,i),idag,coact,bc,gblclr,vecbl, &
                       vecblinv,myid,nn,ldiv,nms,lvbc,ib,lbd,iblv,MRT2)
!          
                do isite = 1,nvhalf
                   vprime(:,isite,:,:,:,i+1) = vprime(:,isite,:,:,:,i) &
                                               - kappa*tue(:,isite,:,:,:)
                   vprime1(:,isite,:,:,:,i+1) = vprime1(:,isite,:,:,:,i) &
                                                - kappa*tuo(:,isite,:,:,:)
                enddo ! isite
        enddo!i



        do i=2,p+1
          do j=2,p+1
             call vecdot(vprime(:,:,:,:,:,i),vprime(:,:,:,:,:,j),beta,MRT2)
             call vecdot(vprime1(:,:,:,:,:,i),vprime1(:,:,:,:,:,j),beta1,MRT2)
                  lsmat1(1,i-1,j-1) = beta(1) + beta1(1)  !lsmat(2,p,p) ,cls(2,p,1)
                  lsmat1(2,i-1,j-1) = beta(2) + beta1(2)
!      if(myid==0) then
!        print *, "i,j, lsmat(:,i,j)=", i-1,j-1, lsmat(:,i-1,j-1)
!      endif
        enddo!j
           enddo!i

       do i=2,p+1
          call vecdot(vprime(:,:,:,:,:,i),vprime(:,:,:,:,:,1),beta,MRT2)
          call vecdot(vprime1(:,:,:,:,:,i),vprime1(:,:,:,:,:,1),beta1,MRT2)
               cls1(1,i-1,1) = beta(1) + beta1(1)
               cls1(2,i-1,1) = beta(2) + beta1(2)
!      if(myid==0) then
!        print *, "i,cls(:,i)=", i-1, cls(:,i-1,1)
!      endif
       enddo!i

    call linearsolver(p,1,lsmat1,ipiv3,cls1)
    co2(:,:) = cls1(:,:,1)
!    co = 0.0_KR2
!    co(1,1) = 4
    
  else
       co2(1,:)=1.0_KR2
       co2(2,:)=0.0_KR2
  endif !gammaid
    
    ! if(myid==0) then
     !  do i=1,p
      !   print *, "fortran coefficients i co2(:,i)=", i, co2(:,i)
       !enddo!i
    ! endif!myid

do eig = 1, numEV



      be = 0.0_KR2
      bo = 0.0_KR2

      be = evecReven(:,:,:,:,:,eig)
      bo = evecRodd(:,:,:,:,:,eig)


 do isub = 1,nsub

 if (.true.) then


         vprime = 0.0_KR2
         vprime1 = 0.0_KR2
        do isite = 1,nvhalf
           vprime(:,isite,:,:,:,1) = be(:,isite,:,:,:)
        enddo !isite

        do isite = 1,nvhalf
           vprime1(:,isite,:,:,:,1) = bo(:,isite,:,:,:)
        enddo !isite


        do i = 1,p
           gblclr = 1
           call Hsingle(tue,u,vprime1(:,:,:,:,:,i),idag,coact,bc,gblclr,vecbl, &
                       vecblinv,myid,nn,ldiv,nms,lvbc,ib,lbd,iblv,MRT2)
           gblclr = 2
           call Hsingle(tuo,u,vprime(:,:,:,:,:,i),idag,coact,bc,gblclr,vecbl, &
                       vecblinv,myid,nn,ldiv,nms,lvbc,ib,lbd,iblv,MRT2)
!
                do isite = 1,nvhalf
                   vprime(:,isite,:,:,:,i+1) = vprime(:,isite,:,:,:,i) &
                                               - kappa*tue(:,isite,:,:,:)
                   vprime1(:,isite,:,:,:,i+1) = vprime1(:,isite,:,:,:,i) &
                                                - kappa*tuo(:,isite,:,:,:)
                enddo ! isite
        enddo!i












    ! Compute the operators for each level of subtraction.
   select case(isub)
     case(1) ! no subtraction
       ksub = 0
       sube = 0.0_KR
       subo = 0.0_KR
     case(2) ! subtraction of order1,order2 and order3
       ksub = 1
        do icri=1,5,2
         do isite=1,nvhalf
          sube(icri,isite,:,:,:) =  sube(icri,isite,:,:,:) &
                                    +co2(1,1)*vprime(icri,isite,:,:,:,1) &
                                    -co2(2,1)*vprime(icri+1,isite,:,:,:,1)
          sube(icri+1,isite,:,:,:) = sube(icri+1,isite,:,:,:) &
                                     +co2(1,1)*vprime(icri+1,isite,:,:,:,1) &
                                     +co2(2,1)*vprime(icri,isite,:,:,:,1)

         enddo!isite
        enddo!icri

      do icri=1,5,2
       do isite=1,nvhalf
        subo(icri,isite,:,:,:) =  subo(icri,isite,:,:,:) &
                                  +co2(1,1)*vprime1(icri,isite,:,:,:,1) &
                                  -co2(2,1)*vprime1(icri+1,isite,:,:,:,1)
       subo(icri+1,isite,:,:,:) = subo(icri+1,isite,:,:,:) &
                                  +co2(1,1)*vprime1(icri+1,isite,:,:,:,1) &
                                  +co2(2,1)*vprime1(icri,isite,:,:,:,1)
       enddo!isite
      enddo!icri
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!order0!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        do icri=1,5,2
         do isite=1,nvhalf
          sube(icri,isite,:,:,:) =  sube(icri,isite,:,:,:) &
                                    +co2(1,2)*vprime(icri,isite,:,:,:,2) &
                                    -co2(2,2)*vprime(icri+1,isite,:,:,:,2)
          sube(icri+1,isite,:,:,:) = sube(icri+1,isite,:,:,:) &
                                    +co2(1,2)*vprime(icri+1,isite,:,:,:,2) &
                                    +co2(2,2)*vprime(icri,isite,:,:,:,2)

         enddo!isite
        enddo!icri

      do icri=1,5,2
       do isite=1,nvhalf
        subo(icri,isite,:,:,:) =  subo(icri,isite,:,:,:) &
                                  +co2(1,2)*vprime1(icri,isite,:,:,:,2) &
                                  -co2(2,2)*vprime1(icri+1,isite,:,:,:,2)
        subo(icri+1,isite,:,:,:) = subo(icri+1,isite,:,:,:) &
                                   +co2(1,2)*vprime1(icri+1,isite,:,:,:,2)&
                                   +co2(2,2)*vprime1(icri,isite,:,:,:,2)



       enddo!isite
      enddo!icri

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1order!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        do icri=1,5,2
         do isite=1,nvhalf
          sube(icri,isite,:,:,:) =  sube(icri,isite,:,:,:) &
                                    +co2(1,3)*vprime(icri,isite,:,:,:,3) &
                                    -co2(2,3)*vprime(icri+1,isite,:,:,:,3)
          sube(icri+1,isite,:,:,:) = sube(icri+1,isite,:,:,:) &
                                     +co2(1,3)*vprime(icri+1,isite,:,:,:,3) &
                                     +co2(2,3)*vprime(icri,isite,:,:,:,3)

         enddo!isite
        enddo!icri

      do icri=1,5,2
       do isite=1,nvhalf
        subo(icri,isite,:,:,:) =  subo(icri,isite,:,:,:) &
                                  +co2(1,3)*vprime1(icri,isite,:,:,:,3) &
                                  -co2(2,3)*vprime1(icri+1,isite,:,:,:,3)
        subo(icri+1,isite,:,:,:) = subo(icri+1,isite,:,:,:) &
                                  +co2(1,3)*vprime1(icri+1,isite,:,:,:,3)&
                                  +co2(2,3)*vprime1(icri,isite,:,:,:,3)
       enddo!isite
      enddo!icri


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!2order!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        do icri=1,5,2
         do isite=1,nvhalf
          sube(icri,isite,:,:,:) =  sube(icri,isite,:,:,:)&
                                   +co2(1,4)*vprime(icri,isite,:,:,:,4) &
                                   -co2(2,4)*vprime(icri+1,isite,:,:,:,4)
          sube(icri+1,isite,:,:,:) = sube(icri+1,isite,:,:,:)&
                                    +co2(1,4)*vprime(icri+1,isite,:,:,:,4) &
                                    +co2(2,4)*vprime(icri,isite,:,:,:,4)

         enddo!isite
        enddo!icri

      do icri=1,5,2
       do isite=1,nvhalf
        subo(icri,isite,:,:,:) =  subo(icri,isite,:,:,:) &
                                  +co2(1,4)*vprime1(icri,isite,:,:,:,4)&
                                  -co2(2,4)*vprime1(icri+1,isite,:,:,:,4)
        subo(icri+1,isite,:,:,:) = subo(icri+1,isite,:,:,:) &
                                   +co2(1,4)*vprime1(icri+1,isite,:,:,:,4)&
                                   +co2(2,4)*vprime1(icri,isite,:,:,:,4)
       enddo!isite
      enddo!icri

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!3order!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        case(3) ! subtraction of O(kappa^4)

        do icri=1,5,2
         do isite=1,nvhalf
          sube(icri,isite,:,:,:) =  sube(icri,isite,:,:,:)&
                                   +co2(1,5)*vprime(icri,isite,:,:,:,5) &
                                   -co2(2,5)*vprime(icri+1,isite,:,:,:,5)
          sube(icri+1,isite,:,:,:) = sube(icri+1,isite,:,:,:)&
                                    +co2(1,5)*vprime(icri+1,isite,:,:,:,5) &
                                    +co2(2,5)*vprime(icri,isite,:,:,:,5)

         enddo!isite
        enddo!icri

      do icri=1,5,2
       do isite=1,nvhalf
        subo(icri,isite,:,:,:) =  subo(icri,isite,:,:,:) &
                                  +co2(1,5)*vprime1(icri,isite,:,:,:,5) &
                                  -co2(2,5)*vprime1(icri+1,isite,:,:,:,5)
        subo(icri+1,isite,:,:,:) = subo(icri+1,isite,:,:,:) &
                                   +co2(1,5)*vprime1(icri+1,isite,:,:,:,5)&
                                   +co2(2,5)*vprime1(icri,isite,:,:,:,5)
       enddo!isite
      enddo!icri

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!4order!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        case(4) ! subtraction of O(kappa^5)
         do icri=1,5,2
         do isite=1,nvhalf
          sube(icri,isite,:,:,:) =  sube(icri,isite,:,:,:)&
                                   +co2(1,6)*vprime(icri,isite,:,:,:,6) &
                                   -co2(2,6)*vprime(icri+1,isite,:,:,:,6)
          sube(icri+1,isite,:,:,:) = sube(icri+1,isite,:,:,:)&
                                    +co2(1,6)*vprime(icri+1,isite,:,:,:,6) &
                                    +co2(2,6)*vprime(icri,isite,:,:,:,6)

         enddo!isite
        enddo!icri

      do icri=1,5,2
       do isite=1,nvhalf
        subo(icri,isite,:,:,:) =  subo(icri,isite,:,:,:) &
                                  +co2(1,6)*vprime1(icri,isite,:,:,:,6) &
                                  -co2(2,6)*vprime1(icri+1,isite,:,:,:,6)
        subo(icri+1,isite,:,:,:) = subo(icri+1,isite,:,:,:) &
                                   +co2(1,6)*vprime1(icri+1,isite,:,:,:,6)&
                                   +co2(2,6)*vprime1(icri,isite,:,:,:,6)
       enddo!isite
      enddo!icri
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!5order!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        case(5) ! subtraction of O(kappa^6)
        do icri=1,5,2
         do isite=1,nvhalf
          sube(icri,isite,:,:,:) =  sube(icri,isite,:,:,:)&
                                   +co2(1,7)*vprime(icri,isite,:,:,:,7) &
                                   -co2(2,7)*vprime(icri+1,isite,:,:,:,7)
          sube(icri+1,isite,:,:,:) = sube(icri+1,isite,:,:,:)&
                                    +co2(1,7)*vprime(icri+1,isite,:,:,:,7) &
                                    +co2(2,7)*vprime(icri,isite,:,:,:,7)

         enddo!isite
        enddo!icri

      do icri=1,5,2
       do isite=1,nvhalf
        subo(icri,isite,:,:,:) =  subo(icri,isite,:,:,:) &
                                  +co2(1,7)*vprime1(icri,isite,:,:,:,7) &
                                  -co2(2,7)*vprime1(icri+1,isite,:,:,:,7)
        subo(icri+1,isite,:,:,:) = subo(icri+1,isite,:,:,:) &
                                   +co2(1,7)*vprime1(icri+1,isite,:,:,:,7)&
                                   +co2(2,7)*vprime1(icri,isite,:,:,:,7)
       enddo!isite
      enddo!icri
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!6order!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        case(6) ! subtraction of O(kappa^7)
        do icri=1,5,2
         do isite=1,nvhalf
          sube(icri,isite,:,:,:) =  sube(icri,isite,:,:,:)&
                                   +co2(1,8)*vprime(icri,isite,:,:,:,8) &
                                   -co2(2,8)*vprime(icri+1,isite,:,:,:,8)
          sube(icri+1,isite,:,:,:) = sube(icri+1,isite,:,:,:)&
                                    +co2(1,8)*vprime(icri+1,isite,:,:,:,8) &
                                    +co2(2,8)*vprime(icri,isite,:,:,:,8)

         enddo!isite
        enddo!icri

      do icri=1,5,2
       do isite=1,nvhalf
        subo(icri,isite,:,:,:) =  subo(icri,isite,:,:,:) &
                                  +co2(1,8)*vprime1(icri,isite,:,:,:,8) &
                                  -co2(2,8)*vprime1(icri+1,isite,:,:,:,8)
        subo(icri+1,isite,:,:,:) = subo(icri+1,isite,:,:,:) &
                                   +co2(1,8)*vprime1(icri+1,isite,:,:,:,8)&
                                   +co2(2,8)*vprime1(icri,isite,:,:,:,8)
       enddo!isite
      enddo!icri
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!7order!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        case default
          open(unit=8,file="DISCONLOOPS.ERROR",action="write",status="replace",&
               form="formatted")
          write(unit=8,fmt=*) "subroutine average: isub =", isub
          close(unit=8,status="keep")
          stop
      end select



endif !true or false




if (.false.) then



  ! Compute the operators for each level of subtraction.
      select case(isub)
        case(1) ! no subtraction
          ksub = 0
          sube = 0.0_KR
          subo = 0.0_KR
        case(2) ! subtraction of O(kappa), O(kappa^2) and O(kappa^3)
          ksub = 1
          gblclr = 1
          call Hsingle(sub1e,u,bo,idag,coact,bc,gblclr,vecbl,vecblinv,myid, &
                       nn,ldiv,nms,lvbc,ib,lbd,iblv,MRT2)
          gblclr = 2
          call Hsingle(sub1o,u,be,idag,coact,bc,gblclr,vecbl,vecblinv,myid, &
                       nn,ldiv,nms,lvbc,ib,lbd,iblv,MRT2)
          do ibleo = 1,8
            do ieo = 1,2
              do id = 1,4
                do isite = 1,nvhalf
                  do icri = 1,5,2
                    sube(icri,isite,id,ieo,ibleo) =co2(1,1)*be(icri,isite,id,ieo,ibleo)&
                                                  - co2(2,1)*be(icri+1,isite,id,ieo,ibleo)&
                                         +co2(1,2)*kappa*sub1e(icri,isite,id,ieo,ibleo)&
                                         -co2(2,2)*kappa*sub1e(icri+1,isite,id,ieo,ibleo)
                     sube(icri+1,isite,id,ieo,ibleo) =co2(1,1)*be(icri+1,isite,id,ieo,ibleo)&
                                                     +co2(2,1)*be(icri,isite,id,ieo,ibleo)&
                                         +co2(1,2)*kappa*sub1e(icri+1,isite,id,ieo,ibleo)&
                                         +co2(2,2)*kappa*sub1e(icri,isite,id,ieo,ibleo)


                    subo(icri,isite,id,ieo,ibleo) =co2(1,1)*bo(icri,isite,id,ieo,ibleo)&
                                                  -co2(2,1)*bo(icri+1,isite,id,ieo,ibleo)&
                                         +co2(1,2)*kappa*sub1o(icri,isite,id,ieo,ibleo)&
                                         -co2(2,2)*kappa*sub1o(icri+1,isite,id,ieo,ibleo)
                     subo(icri+1,isite,id,ieo,ibleo)= co2(1,1)*bo(icri+1,isite,id,ieo,ibleo)&
                                                     +co2(2,1)*bo(icri,isite,id,ieo,ibleo)&
                                         +co2(1,2)*kappa*sub1o(icri+1,isite,id,ieo,ibleo)&
                                         +co2(2,2)*kappa*sub1o(icri,isite,id,ieo,ibleo)

                  enddo ! icri
                enddo ! isite
              enddo ! id
            enddo ! ieo
          enddo ! ibleo

          gblclr = 1
          call Hsingle(sub2e,u,sub1o,idag,coact,bc,gblclr,vecbl,vecblinv,myid, &
                       nn,ldiv,nms,lvbc,ib,lbd,iblv,MRT2)
          gblclr = 2
          call Hsingle(sub2o,u,sub1e,idag,coact,bc,gblclr,vecbl,vecblinv,myid, &
                       nn,ldiv,nms,lvbc,ib,lbd,iblv,MRT2)
          xk = kappa**2
          do ibleo = 1,8
            do ieo = 1,2
              do id = 1,4
                do isite = 1,nvhalf
                  do icri = 1,5,2

                    sube(icri,isite,id,ieo,ibleo) =sube(icri,isite,id,ieo,ibleo)&
                                              +co2(1,3)*xk*sub2e(icri,isite,id,ieo,ibleo)&
                                              -co2(2,3)*xk*sub2e(icri+1,isite,id,ieo,ibleo)

                    sube(icri+1,isite,id,ieo,ibleo)=sube(icri+1,isite,id,ieo,ibleo)&
                                              +co2(1,3)*xk*sub2e(icri+1,isite,id,ieo,ibleo)&
                                              +co2(2,3)*xk*sub2e(icri,isite,id,ieo,ibleo)

                    subo(icri,isite,id,ieo,ibleo) =subo(icri,isite,id,ieo,ibleo)&
                                              +co2(1,3)*xk*sub2o(icri,isite,id,ieo,ibleo)&
                                              -co2(2,3)*xk*sub2e(icri+1,isite,id,ieo,ibleo)

                    subo(icri+1,isite,id,ieo,ibleo)=subo(icri+1,isite,id,ieo,ibleo)&
                                              +co2(1,3)*xk*sub2o(icri+1,isite,id,ieo,ibleo)&
                                              +co2(2,3)*xk*sub2e(icri,isite,id,ieo,ibleo)

                  enddo ! icri
                enddo ! isite
              enddo ! id
            enddo ! ieo
          enddo ! ibleo

          gblclr = 1
          call Hsingle(sub1e,u,sub2o,idag,coact,bc,gblclr,vecbl,vecblinv,myid, &
                       nn,ldiv,nms,lvbc,ib,lbd,iblv,MRT2)
          gblclr = 2
          call Hsingle(sub1o,u,sub2e,idag,coact,bc,gblclr,vecbl,vecblinv,myid, &
                       nn,ldiv,nms,lvbc,ib,lbd,iblv,MRT2)
          xk = kappa**3
          do ibleo = 1,8
            do ieo = 1,2
              do id = 1,4
                do isite = 1,nvhalf
                   do icri = 1,5,2

                    sube(icri,isite,id,ieo,ibleo)=sube(icri,isite,id,ieo,ibleo)&
                                              +co2(1,4)*xk*sub1e(icri,isite,id,ieo,ibleo)&
                                              -co2(2,4)*xk*sub1e(icri+1,isite,id,ieo,ibleo)

                    sube(icri+1,isite,id,ieo,ibleo)=sube(icri+1,isite,id,ieo,ibleo)&
                                              +co2(1,4)*xk*sub1e(icri+1,isite,id,ieo,ibleo)&
                                              +co2(2,4)*xk*sub1e(icri,isite,id,ieo,ibleo)

                    subo(icri,isite,id,ieo,ibleo)=subo(icri,isite,id,ieo,ibleo)&
                                              +co2(1,4)*xk*sub1o(icri,isite,id,ieo,ibleo)&
                                              -co2(2,4)*xk*sub1o(icri+1,isite,id,ieo,ibleo)

                    subo(icri+1,isite,id,ieo,ibleo)=subo(icri+1,isite,id,ieo,ibleo)&
                                              +co2(1,4)*xk*sub1o(icri+1,isite,id,ieo,ibleo)&
                                              +co2(2,4)*xk*sub1o(icri,isite,id,ieo,ibleo)



                  enddo ! icri
                enddo ! isite
              enddo ! id
            enddo ! ieo
          enddo ! ibleo




        case(3) ! subtraction of O(kappa^4)
          gblclr = 1
          call Hsingle(sub2e,u,sub1o,idag,coact,bc,gblclr,vecbl,vecblinv,myid, &
                       nn,ldiv,nms,lvbc,ib,lbd,iblv,MRT2)
          gblclr = 2
          call Hsingle(sub2o,u,sub1e,idag,coact,bc,gblclr,vecbl,vecblinv,myid, &
                       nn,ldiv,nms,lvbc,ib,lbd,iblv,MRT2)
          xk = kappa**4
          do ibleo = 1,8
            do ieo = 1,2
              do id = 1,4
                do isite = 1,nvhalf
                  do icri = 1,5,2

                    sube(icri,isite,id,ieo,ibleo)=sube(icri,isite,id,ieo,ibleo)&
                                              +co2(1,5)*xk*sub2e(icri,isite,id,ieo,ibleo)&
                                              -co2(2,5)*xk*sub2e(icri+1,isite,id,ieo,ibleo)

                    sube(icri+1,isite,id,ieo,ibleo)=sube(icri+1,isite,id,ieo,ibleo)&
                                              +co2(1,5)*xk*sub2e(icri+1,isite,id,ieo,ibleo)&
                                              +co2(2,5)*xk*sub2e(icri,isite,id,ieo,ibleo)

                    subo(icri,isite,id,ieo,ibleo)=subo(icri,isite,id,ieo,ibleo)&
                                              +co2(1,5)*xk*sub2o(icri,isite,id,ieo,ibleo)&
                                              -co2(2,5)*xk*sub2o(icri+1,isite,id,ieo,ibleo)


                    subo(icri+1,isite,id,ieo,ibleo)=subo(icri+1,isite,id,ieo,ibleo)&
                                              +co2(1,5)*xk*sub2o(icri+1,isite,id,ieo,ibleo)&
                                              +co2(2,5)*xk*sub2o(icri,isite,id,ieo,ibleo)

                  enddo ! icri
                enddo ! isite
              enddo ! id
            enddo ! ieo
          enddo ! ibleo






        case(4) ! subtraction of O(kappa^5)
          gblclr = 1
          call Hsingle(sub1e,u,sub2o,idag,coact,bc,gblclr,vecbl,vecblinv,myid, &
                       nn,ldiv,nms,lvbc,ib,lbd,iblv,MRT2)
          gblclr = 2
          call Hsingle(sub1o,u,sub2e,idag,coact,bc,gblclr,vecbl,vecblinv,myid, &
                       nn,ldiv,nms,lvbc,ib,lbd,iblv,MRT2)
          xk = kappa**5
          do ibleo = 1,8
            do ieo = 1,2
              do id = 1,4
                do isite = 1,nvhalf
                  do icri = 1,5,2

                    sube(icri,isite,id,ieo,ibleo)=sube(icri,isite,id,ieo,ibleo)&
                                              +co2(1,6)*xk*sub1e(icri,isite,id,ieo,ibleo)&
                                              -co2(2,6)*xk*sub1e(icri+1,isite,id,ieo,ibleo)

                    sube(icri+1,isite,id,ieo,ibleo)=sube(icri+1,isite,id,ieo,ibleo)&
                                              +co2(1,6)*xk*sub1e(icri+1,isite,id,ieo,ibleo)&
                                              +co2(2,6)*xk*sub1e(icri,isite,id,ieo,ibleo)

                    subo(icri,isite,id,ieo,ibleo)=subo(icri,isite,id,ieo,ibleo)&
                                              +co2(1,6)*xk*sub1o(icri,isite,id,ieo,ibleo)&
                                              -co2(2,6)*xk*sub1o(icri+1,isite,id,ieo,ibleo)

                    subo(icri+1,isite,id,ieo,ibleo)=subo(icri+1,isite,id,ieo,ibleo)&
                                              +co2(1,6)*xk*sub1o(icri+1,isite,id,ieo,ibleo)&
                                              +co2(2,6)*xk*sub1o(icri,isite,id,ieo,ibleo)



                  enddo ! icri
                enddo ! isite
              enddo ! id
            enddo ! ieo
          enddo ! ibleo





       case(5) ! subtraction of O(kappa^6)
          gblclr = 1
          call Hsingle(sub2e,u,sub1o,idag,coact,bc,gblclr,vecbl,vecblinv,myid, &
                       nn,ldiv,nms,lvbc,ib,lbd,iblv,MRT2)
          gblclr = 2
          call Hsingle(sub2o,u,sub1e,idag,coact,bc,gblclr,vecbl,vecblinv,myid, &
                       nn,ldiv,nms,lvbc,ib,lbd,iblv,MRT2)
          xk = kappa**6
          do ibleo = 1,8
            do ieo = 1,2
              do id = 1,4
                do isite = 1,nvhalf
                   do icri = 1,5,2

                    sube(icri,isite,id,ieo,ibleo)=sube(icri,isite,id,ieo,ibleo)&
                                              +co2(1,7)*xk*sub2e(icri,isite,id,ieo,ibleo)&
                                              -co2(2,7)*xk*sub2e(icri+1,isite,id,ieo,ibleo)

                    sube(icri+1,isite,id,ieo,ibleo)=sube(icri+1,isite,id,ieo,ibleo)&
                                              +co2(1,7)*xk*sub2e(icri+1,isite,id,ieo,ibleo)&
                                              +co2(2,7)*xk*sub2e(icri,isite,id,ieo,ibleo)

                    subo(icri,isite,id,ieo,ibleo)=subo(icri,isite,id,ieo,ibleo)&
                                              +co2(1,7)*xk*sub2o(icri,isite,id,ieo,ibleo)&
                                              -co2(2,7)*xk*sub2o(icri+1,isite,id,ieo,ibleo)

                    subo(icri+1,isite,id,ieo,ibleo)=subo(icri+1,isite,id,ieo,ibleo)&
                                              +co2(1,7)*xk*sub2o(icri+1,isite,id,ieo,ibleo)&
                                              +co2(2,7)*xk*sub2o(icri,isite,id,ieo,ibleo)




                  enddo ! icri
                enddo ! isite
              enddo ! id
            enddo ! ieo
          enddo ! ibleo




        case(6) ! subtraction of O(kappa^7)
          gblclr = 1
          call Hsingle(sub1e,u,sub2o,idag,coact,bc,gblclr,vecbl,vecblinv,myid, &
                       nn,ldiv,nms,lvbc,ib,lbd,iblv,MRT2)
          gblclr = 2
          call Hsingle(sub1o,u,sub2e,idag,coact,bc,gblclr,vecbl,vecblinv,myid, &
                       nn,ldiv,nms,lvbc,ib,lbd,iblv,MRT2)
          xk = kappa**7!BS 12/23/2015 changed 6 to 7
          do ibleo = 1,8
            do ieo = 1,2
              do id = 1,4
                do isite = 1,nvhalf
                  do icri = 1,5,2

                    sube(icri,isite,id,ieo,ibleo)=sube(icri,isite,id,ieo,ibleo)&
                                              +co2(1,8)*xk*sub1e(icri,isite,id,ieo,ibleo)&
                                              -co2(2,8)*xk*sub1e(icri+1,isite,id,ieo,ibleo)

                    sube(icri+1,isite,id,ieo,ibleo)=sube(icri+1,isite,id,ieo,ibleo)&
                                              +co2(1,8)*xk*sub1e(icri+1,isite,id,ieo,ibleo)&
                                              +co2(2,8)*xk*sub1e(icri,isite,id,ieo,ibleo)

                    subo(icri,isite,id,ieo,ibleo)=subo(icri,isite,id,ieo,ibleo)&
                                              +co2(1,8)*xk*sub1o(icri,isite,id,ieo,ibleo)&
                                              -co2(2,8)*xk*sub1o(icri+1,isite,id,ieo,ibleo)

                    subo(icri+1,isite,id,ieo,ibleo)=subo(icri+1,isite,id,ieo,ibleo)&
                                              +co2(1,8)*xk*sub1o(icri+1,isite,id,ieo,ibleo)&
                                              +co2(2,8)*xk*sub1o(icri,isite,id,ieo,ibleo)

                  enddo ! icri
                enddo ! isite
              enddo ! id
            enddo ! ieo
          enddo ! ibleo




      case default





          open(unit=8,file="DISCONLOOPS.ERROR",action="write",status="replace",&
               form="formatted")
          write(unit=8,fmt=*) "subroutine average: isub =", isub
          close(unit=8,status="keep")
          stop
      end select
endif !true or false


        temptempsube=sube !BS 4/25/2016
        temptempsubo=subo !BS 4/25/2016






    !BS start changes...flip gamma5 ...see vic's dissertation equation 4.83
     !BS  if (gammaid == 1_KI) then
     !  if (gammaid .NE. 1_KI) then
           do ieo=1,2
             do ibleo=1,8
               call gammamult(sube, subo, tempsube,tempsubo, 5,ieo,ibleo)
             enddo
           enddo
        sube = tempsube
        subo = tempsubo
      !endif
    ! BS end changes here

      call vecdot(evecLeven(:,:,:,:,:,eig),sube,tmp1,MRT2)
      call vecdot(evecLodd(:,:,:,:,:,eig),subo,tmp2,MRT2)
      xiinvPOLY(:2,isub, eig) = tmp1(:2) + tmp2(:2)


   sube=temptempsube
   subo=temptempsubo
    enddo
 enddo


 end subroutine xipolygenerator








 subroutine eigmodesofM(evectorRevn,evectorRodd,evectorLevn,evectorLodd,eigValExpand,numEV, &
                  u,kappa,coact,bc,vecbl,vecblinv,myid,nn,ldiv,nms,lvbc,             &
                  lvcl,ib,lbd,iblv,rwdir,MRT,MRT2,nGMRES)
 real(kind=KR2),   intent(out), dimension(:,:,:,:,:,:) :: evectorRevn,evectorRodd
 real(kind=KR2),   intent(out), dimension(:,:,:,:,:,:) :: evectorLevn,evectorLodd
 real(kind=KR2),   intent(out), dimension(:,:)         :: eigValExpand
 integer(kind=KI), intent(out)                         :: numEV
 ! Identify the number of modes to try to calculate
 integer(kind=KI), intent(in),  dimension(:)           :: nGMRES

 real(kind=KR),    intent(inout), dimension(:,:,:,:,:) :: u
 real(kind=KR),    intent(in),    dimension(:)         :: kappa
 character(len=*), intent(in),    dimension(:)         :: rwdir
 integer(kind=KI), intent(in),    dimension(:)         :: bc, nms
 integer(kind=KI), intent(in)                          :: myid, MRT,MRT2
 real(kind=KR),    intent(in),    dimension(:,:,:)     :: coact
 integer(kind=KI), intent(in),    dimension(:,:)       :: vecbl, vecblinv
 integer(kind=KI), intent(in),    dimension(:,:)       :: nn, iblv
 logical,          intent(in),    dimension(:)         :: ldiv
 integer(kind=KI), intent(in),    dimension(:,:,:)     :: lvbc, lvcl
 integer(kind=KI), intent(in),    dimension(:,:,:,:)   :: ib
 logical,          intent(in),    dimension(:,:)       :: lbd

 real(kind=KR2),  dimension(6,ntotal,4,2,8,kmaxGMRES) :: evectorRevnTMP, evectorRoddTMP
 real(kind=KR2),  dimension(6,ntotal,4,2,8,kmaxGMRES) :: evectorLevnTMP, evectorLoddTMP

 real(kind=KR2),  dimension(2,kmaxGMRES)      :: eigValExpandTMP
 real(kind=KR2),  dimension(kmaxGMRES)        :: evalueNORM


 real(kind=KR2),  dimension(6,ntotal,4,2,8)   :: tmpE, tmpO
 real(kind=KR2),  dimension(6,ntotal,4,2,8)   :: bepart, bopart
 real(kind=KR2),  dimension(6,ntotal,4,2,8)   :: b, xe,xo
 integer(kind=KI)                         :: idag
 real(kind=KR2),  dimension(2) :: tmp1,tmp2,tmp3,tmp4
 integer(kind=KI)  :: eig,eigpos,eigneg,eigrel,eigpair,eigprime
 integer(kind=KI)  :: iaa,j
 real(kind=KR2), dimension(2,kmaxGMRES) :: evalrel,evalpos,evalneg
 integer(kind=KI), dimension(kmaxGMRES) :: indrel,indpos,indneg,ind,ind2

 ! GMRES Variables
 real(kind=KR), dimension(18,nvhalf,8,2,16) :: GeeGooinv
 real(kind=KR2)     :: resmax
 integer(kind=KI)   :: itermin,iflag,itercount

 real(kind=KR2),   parameter        :: eigtol = 1e-8_KR2
 real(kind=KR2),   parameter        :: eigtoltemp = 1e-8_KR2
  integer(kind=KI) :: exists 
 ! initialize/declare vairable parameters
 idag = 0; ! guarentees Mx not Mdagx
 itermin = 10
 iflag = -1  ! Wilson
 GeeGooinv = 0.0_KR2

 ! GMRES parameters
 resmax = 1e-155_KR2

 !!!! -- Generate Eig-info from Mreduced -- !!!
 b(:6,:ntotal,:4,:2,:8) = 0.0_KR2   ! Arbitrary (non-zero) rhs
 b(1,1,1,1,1) = 1.0_KR2             ! Arbitrary (non-zero) rhs
 xe(:6,:ntotal,:4,:2,:8) = 0.0_KR2  ! Initial Guess
 call gmresdrEIGBS(rwdir,b,xe,nGMRES,resmax,itermin,itercount,u,GeeGooinv, &
                 iflag,idag,kappa,coact,bc,vecbl,vecblinv,myid,nn,ldiv,nms, &
                 lvbc,ib,lbd,iblv,MRT,MRT2)

!print *,'evaleig',evalue

 ! Separate modes into 3 catagories: {a+bi, a-bi, c (purly real)}
 eigpos = 0  ! number of a+bi modes
 eigneg = 0  ! number of a-bi modes
 eigrel = 0  ! number of c modes (purly real)
 do eig=1,nGMRES(2)
   if (abs(evalue(2,eig)) < eigtol) then ! c purly real
     eigrel = eigrel + 1
     indrel(eigrel) = eig
   elseif (evalue(2,eig) > 0) then       ! a+bi
     eigpos = eigpos + 1
     indpos(eigpos) = eig
   else                                  ! a-bi
     eigneg = eigneg + 1
     indneg(eigneg) = eig
   endif
 enddo

 ! Idenfity {a+bi,a-bi} pairs which are needed for L/R eigenmode calculations
 ! eigpair will hold the # of acceptable pair modes
 eigpair = 0 ! number of pairs found
 do eig=1,eigpos
   do eigprime=1,eigneg
     ! tmp1 = evalpos(eig) - conj(evalneg(eigprime))
     tmp1(1) = evalue(1,indpos(eig)) - evalue(1,indneg(eigprime))
     tmp1(2) = evalue(2,indpos(eig)) + evalue(2,indneg(eigprime))
     ! tmp1(1) = |tmp1| (abs norm)
     tmp1(1) = sqrt(tmp1(1)*tmp1(1) + tmp1(2)*tmp1(2))

     ! check if tol is made
     if (tmp1(1) <= eigtol) then
       eigpair = eigpair + 1
       ind(eigpair) = indpos(eig)
       eigpair = eigpair + 1
       ind(eigpair) = indneg(eigprime)
       exit ! break out of eigprime loop
     endif
   enddo
 enddo

 ! --- Form Full Right evectors & evaules ---
 do eig=1,eigpair
   evectorRevnTMP(:6,:ntotal,:4,:2,:8,eig) = evector(:6,:ntotal,:4,:2,:8,ind(eig))

   ! begin forming odd part of right evector
   call Hsingle(bopart,u,evectorRevnTMP(:,:,:,:,:,eig),idag,coact,bc,2_KI, &
                vecbl,vecblinv,myid,nn,ldiv,nms,lvbc,ib,lbd,iblv,MRT)
   bopart(:6,:nvhalf,:4,:2,:8)  = kappa(1)*bopart(:6,:nvhalf,:4,:2,:8) 
   tmp1(:2) = -kappa(1)*kappa(1)*evalue(:2,ind(eig))
   tmp1(1) = 1.0_KR2 + tmp1(1)
   call cmpxsqrt(tmp1,tmp2)

   ! eigenvalue is here
   eigValExpandTMP(:2,eig) = -tmp2(:2)
   eigValExpandTMP(1,eig) = 1.0_KR2 + eigValExpandTMP(1,eig)

   ! tmp1 = 1 / sqrt(1 - k*k*lamhat)
   call oneover(tmp2,tmp1)

   ! odd eigenvector
   do j = 1,5,2
     evectorRoddTMP(j  ,:nvhalf,:4,:2,:8,eig) = tmp1(1)*bopart(j  ,:nvhalf,:4,:2,:8) &
                                               -tmp1(2)*bopart(j+1,:nvhalf,:4,:2,:8)
     evectorRoddTMP(j+1,:nvhalf,:4,:2,:8,eig) = tmp1(2)*bopart(j  ,:nvhalf,:4,:2,:8) &
                                               +tmp1(1)*bopart(j+1,:nvhalf,:4,:2,:8)
   enddo ! j
 enddo

 !!! -- Form Full LEFT Evectors for Pairs -- !!!
 ! compute left eigenvectors; evL = gamma5*evR
 do eig=1,eigpair,2
   do j=1,5,2
     evectorLevnTMP(j  ,:ntotal,1,:2,:8,eig) =  evectorRevnTMP(j+1,:ntotal,3,:2,:8,eig+1)
     evectorLevnTMP(j+1,:ntotal,1,:2,:8,eig) = -evectorRevnTMP(j  ,:ntotal,3,:2,:8,eig+1)
     evectorLoddTMP(j  ,:ntotal,1,:2,:8,eig) =  evectorRoddTMP(j+1,:ntotal,3,:2,:8,eig+1)
     evectorLoddTMP(j+1,:ntotal,1,:2,:8,eig) = -evectorRoddTMP(j  ,:ntotal,3,:2,:8,eig+1)

     evectorLevnTMP(j  ,:ntotal,1,:2,:8,eig+1) =  evectorRevnTMP(j+1,:ntotal,3,:2,:8,eig)
     evectorLevnTMP(j+1,:ntotal,1,:2,:8,eig+1) = -evectorRevnTMP(j  ,:ntotal,3,:2,:8,eig)
     evectorLoddTMP(j  ,:ntotal,1,:2,:8,eig+1) =  evectorRoddTMP(j+1,:ntotal,3,:2,:8,eig)
     evectorLoddTMP(j+1,:ntotal,1,:2,:8,eig+1) = -evectorRoddTMP(j  ,:ntotal,3,:2,:8,eig)


     evectorLevnTMP(j  ,:ntotal,2,:2,:8,eig) =  evectorRevnTMP(j+1,:ntotal,4,:2,:8,eig+1)
     evectorLevnTMP(j+1,:ntotal,2,:2,:8,eig) = -evectorRevnTMP(j  ,:ntotal,4,:2,:8,eig+1)
     evectorLoddTMP(j  ,:ntotal,2,:2,:8,eig) =  evectorRoddTMP(j+1,:ntotal,4,:2,:8,eig+1)
     evectorLoddTMP(j+1,:ntotal,2,:2,:8,eig) = -evectorRoddTMP(j  ,:ntotal,4,:2,:8,eig+1)

     evectorLevnTMP(j  ,:ntotal,2,:2,:8,eig+1) =  evectorRevnTMP(j+1,:ntotal,4,:2,:8,eig)
     evectorLevnTMP(j+1,:ntotal,2,:2,:8,eig+1) = -evectorRevnTMP(j  ,:ntotal,4,:2,:8,eig)
     evectorLoddTMP(j  ,:ntotal,2,:2,:8,eig+1) =  evectorRoddTMP(j+1,:ntotal,4,:2,:8,eig)
     evectorLoddTMP(j+1,:ntotal,2,:2,:8,eig+1) = -evectorRoddTMP(j  ,:ntotal,4,:2,:8,eig)


     evectorLevnTMP(j  ,:ntotal,3,:2,:8,eig) = -evectorRevnTMP(j+1,:ntotal,1,:2,:8,eig+1)
     evectorLevnTMP(j+1,:ntotal,3,:2,:8,eig) =  evectorRevnTMP(j  ,:ntotal,1,:2,:8,eig+1)
     evectorLoddTMP(j  ,:ntotal,3,:2,:8,eig) = -evectorRoddTMP(j+1,:ntotal,1,:2,:8,eig+1)
     evectorLoddTMP(j+1,:ntotal,3,:2,:8,eig) =  evectorRoddTMP(j  ,:ntotal,1,:2,:8,eig+1)

     evectorLevnTMP(j  ,:ntotal,3,:2,:8,eig+1) = -evectorRevnTMP(j+1,:ntotal,1,:2,:8,eig)
     evectorLevnTMP(j+1,:ntotal,3,:2,:8,eig+1) =  evectorRevnTMP(j  ,:ntotal,1,:2,:8,eig)
     evectorLoddTMP(j  ,:ntotal,3,:2,:8,eig+1) = -evectorRoddTMP(j+1,:ntotal,1,:2,:8,eig)
     evectorLoddTMP(j+1,:ntotal,3,:2,:8,eig+1) =  evectorRoddTMP(j  ,:ntotal,1,:2,:8,eig)


     evectorLevnTMP(j  ,:ntotal,4,:2,:8,eig) = -evectorRevnTMP(j+1,:ntotal,2,:2,:8,eig+1)
     evectorLevnTMP(j+1,:ntotal,4,:2,:8,eig) =  evectorRevnTMP(j  ,:ntotal,2,:2,:8,eig+1)
     evectorLoddTMP(j  ,:ntotal,4,:2,:8,eig) = -evectorRoddTMP(j+1,:ntotal,2,:2,:8,eig+1)
     evectorLoddTMP(j+1,:ntotal,4,:2,:8,eig) =  evectorRoddTMP(j  ,:ntotal,2,:2,:8,eig+1)

     evectorLevnTMP(j  ,:ntotal,4,:2,:8,eig+1) = -evectorRevnTMP(j+1,:ntotal,2,:2,:8,eig)
     evectorLevnTMP(j+1,:ntotal,4,:2,:8,eig+1) =  evectorRevnTMP(j  ,:ntotal,2,:2,:8,eig)
     evectorLoddTMP(j  ,:ntotal,4,:2,:8,eig+1) = -evectorRoddTMP(j+1,:ntotal,2,:2,:8,eig)
     evectorLoddTMP(j+1,:ntotal,4,:2,:8,eig+1) =  evectorRoddTMP(j  ,:ntotal,2,:2,:8,eig)
   enddo
 enddo

 ! Normalize R/L eigenvectors
 do eig=1,eigpair
   call vecdot(evectorRevnTMP(:,:,:,:,:,eig),evectorLevnTMP(:,:,:,:,:,eig),tmp1,MRT)
   call vecdot(evectorRoddTMP(:,:,:,:,:,eig),evectorLoddTMP(:,:,:,:,:,eig),tmp2,MRT)
   tmp1 = tmp1 + tmp2

   call oneover(tmp1,tmp2)

   xe = evectorRevnTMP(:,:,:,:,:,eig)
   xo = evectorRoddTMP(:,:,:,:,:,eig)
   do j=1,5,2
     evectorRevnTMP(j  ,:,:,:,:,eig) = tmp2(1) * xe(j  ,:,:,:,:) &
                                     + tmp2(2) * xe(j+1,:,:,:,:)
     evectorRevnTMP(j+1,:,:,:,:,eig) = tmp2(1) * xe(j+1,:,:,:,:) &
                                     - tmp2(2) * xe(j  ,:,:,:,:)
     evectorRoddTMP(j  ,:,:,:,:,eig) = tmp2(1) * xo(j  ,:,:,:,:) &
                                     + tmp2(2) * xo(j+1,:,:,:,:)
     evectorRoddTMP(j+1,:,:,:,:,eig) = tmp2(1) * xo(j+1,:,:,:,:) &
                                     - tmp2(2) * xo(j  ,:,:,:,:)
   enddo
 enddo

 ! At this point the R/L mode are fully created, but only for the eigpairs, not
 ! for the purly real modes. They are also not sorted from smallest to greatest
 ! in magnitude. 

 !!! -- Check Residual of Purely Real Modes -- !!!
 ! Since their is no pair to compare the evalue to I just compute the residual
 ! for evector, if it is allowable, then the left is computed and compared!

 ! **********This needes to be added later********
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 ! measure accuracy of eigenvector/value combo, here there is no need it identify "pairs"
 ! since each left/right odd/even part has already been formed. 
 numEV = 0 ! number of accurate modes
 do eig=1,eigpair
   ! x=M*evector
   call Hsingle(xo,u,evectorRevnTMP(:,:,:,:,:,eig),idag,coact,bc,2_KI, &
                vecbl,vecblinv,myid,nn,ldiv,nms,lvbc,ib,lbd,iblv,MRT)

   call Hsingle(xe,u,evectorRoddTMP(:,:,:,:,:,eig),idag,coact,bc,1_KI, &
                vecbl,vecblinv,myid,nn,ldiv,nms,lvbc,ib,lbd,iblv,MRT)

   xe(:6,:nvhalf,:4,:2,:8) = evectorRevnTMP(:6,:nvhalf,:4,:2,:8,eig) - kappa(1)*xe(:6,:nvhalf,:4,:2,:8)
   xo(:6,:nvhalf,:4,:2,:8) = evectorRoddTMP(:6,:nvhalf,:4,:2,:8,eig) - kappa(1)*xo(:6,:nvhalf,:4,:2,:8)

   
   ! b=eigval * evector
   do j=1,5,2
     bepart(j  ,:nvhalf,:4,:2,:8) = eigValExpandTMP(1,eig) * evectorRevnTMP(j  ,:nvhalf,:4,:2,:8,eig) &
                                  - eigValExpandTMP(2,eig) * evectorRevnTMP(j+1,:nvhalf,:4,:2,:8,eig)
     bepart(j+1,:nvhalf,:4,:2,:8) = eigValExpandTMP(1,eig) * evectorRevnTMP(j+1,:nvhalf,:4,:2,:8,eig) &
                                  + eigValExpandTMP(2,eig) * evectorRevnTMP(j  ,:nvhalf,:4,:2,:8,eig)
     bopart(j  ,:nvhalf,:4,:2,:8) = eigValExpandTMP(1,eig) * evectorRoddTMP(j  ,:nvhalf,:4,:2,:8,eig) &
                                  - eigValExpandTMP(2,eig) * evectorRoddTMP(j+1,:nvhalf,:4,:2,:8,eig)
     bopart(j+1,:nvhalf,:4,:2,:8) = eigValExpandTMP(1,eig) * evectorRoddTMP(j+1,:nvhalf,:4,:2,:8,eig) &
                                  + eigValExpandTMP(2,eig) * evectorRoddTMP(j  ,:nvhalf,:4,:2,:8,eig)
   enddo

   ! x = x - b
   xe(:6,:nvhalf,:4,:2,:8) = xe(:6,:nvhalf,:4,:2,:8) - bepart(:6,:nvhalf,:4,:2,:8)
   xo(:6,:nvhalf,:4,:2,:8) = xo(:6,:nvhalf,:4,:2,:8) - bopart(:6,:nvhalf,:4,:2,:8)

   call vecdot(xe,xe,tmp1,MRT2)
   call vecdot(xo,xo,tmp2,MRT2)
   tmp1 = tmp1 + tmp2
   call cmpxsqrt(tmp1,tmp2)

!BS*********printing of tmp2.dat starts*************
if (.false.) then 

 if (myid==0) then

       inquire(file=trim(rwdir(myid+1))//"tmp2.dat", exist=exists)
           if (.not. exists) then
               open(unit=43,file=trim(rwdir(myid+1))//"tmp2.dat",status="new",&
                   action="write",form="formatted")
               close(unit=43,status="keep")
           endif



           open(unit=43,file=trim(rwdir(myid+1))//"tmp2.dat",status="old",&
                  action="write",form="formatted",position="append")
                              write(unit=43,fmt="(i5,a7,F17.13,a7,F17.13)") eig," ",tmp2(1),"  ",tmp2(2)

            close(unit=43,status="keep")
 endif


endif!true or false

!BS *******printing of tmp2.dat ends here








   if (tmp2(1) < eigtoltemp) then
     numEV = numEV + 1
     ind(numEV) = eig
   endif
 enddo

 !!! -- Sort Modes (smallest in magnitude to largest) -- !!!
 do eig=1,numEV
   evalueNORM(eig) = sqrt(eigValExpandTMP(1,ind(eig))*eigValExpandTMP(1,ind(eig)) &
                        + eigValExpandTMP(2,ind(eig))*eigValExpandTMP(2,ind(eig)))
 enddo
 call sort(evalueNORM,numEV,ind2)

 ! Put eigmodes in final holding variable
 do eig=1,numEV
   eigValExpand(:,eig) = eigValExpandTMP(:,ind2(ind(eig)))

   evectorRevn(:,:,:,:,:,eig) = evectorRevnTMP(:,:,:,:,:,ind2(ind(eig)))
   evectorRodd(:,:,:,:,:,eig) = evectorRoddTMP(:,:,:,:,:,ind2(ind(eig)))

   evectorLevn(:,:,:,:,:,eig) = evectorLevnTMP(:,:,:,:,:,ind2(ind(eig)))
   evectorLodd(:,:,:,:,:,eig) = evectorLoddTMP(:,:,:,:,:,ind2(ind(eig)))
 enddo
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 ! OUTPUT EIG INFO TO FILE
 if (myid==0) then
   open(unit=12,file=trim(rwdir(myid+1))//"EIGEN_VALS.LOG",status="old",&
        action="write",form="formatted",position="append")
   write(unit=12,fmt='(a60)')"-------GMRES-DR (MATLAB ver) (full computed)---------"
   write(unit=12,fmt='(a22,i4)') "               numEV: ",numEV
   do eig=1,numEV
     write(unit=12,fmt='(i5,a5,f25.15,a5,f25.15)')  eig,"     ",eigValExpand(1,eig), & 
                                                        "     ",eigValExpand(2,eig)
   enddo
   close(unit=12,status="keep")
 endif

 !{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{
 ! DEBUGGING PURPOSES
 ! *set to .true. to test L/R orthoganlity for all full eigenvectors (small and large)
 ! *only non-zero values are outputed to screen which would represent normality or
 !      non orthoganality
 ! *output is sent to standard out
 if (.true.) then
   if (myid==0) write(*,*) 'L/R Orthoganlity test for modes of M'
   do j=1,numEV
     do eig=1,numEV
       call vecdot(evectorRevn(:,:,:,:,:,eig),evectorLevn(:,:,:,:,:,j),tmp1,MRT)
       call vecdot(evectorRodd(:,:,:,:,:,eig),evectorLodd(:,:,:,:,:,j),tmp2,MRT)
       tmp3(:2) = tmp1(:2) + tmp2(:2)
       if (myid==0) then
         if (abs(tmp3(1)) > 1e-6_KR2 .or. abs(tmp3(2)) > 1e-6_KR2) write(*,*) eig,j,tmp3
       endif
     enddo
   enddo
 endif
 !}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}
 end subroutine eigmodesofM



 subroutine printz2noise(z2e,z2o,nvhalf,numprocs,rwdir,MRT,myid)
 real(kind=KR2),   intent(in), dimension(:,:,:,:,:) :: z2e,z2o
 character(len=*), intent(in), dimension(:)         :: rwdir
 integer(kind=KI), intent(in) :: nvhalf
 integer(kind=KI), intent(in) :: numprocs, MRT, myid

 real(kind=KR2), dimension(nxyzt,3,4) :: Rz2,Iz2
 integer(kind=KI) :: i,j,k
 logical :: file_exists

 call changevector(Rz2,Iz2,z2e,z2o(:,:nvhalf,:,:,:),numprocs,MRT,myid)

 if (myid == 0) then
   inquire(FILE=trim(rwdir(1))//"z2noise.dat", EXIST=file_exists)

   if (file_exists) then
     open(unit=40,file=trim(rwdir(myid+1))//"z2noise.dat",status="old", &
          action="write",form="formatted",position="append")
   else
     open(unit=40,file=trim(rwdir(myid+1))//"z2noise.dat",status="new", &
          action="write",form="formatted",position="append")
   endif

   do i=1,nxyzt
     do j=1,3
       do k=1,4
         write(unit=40,fmt='(f25.15,a3,f25.15)') Rz2(i,j,k),'   ',Iz2(i,j,k)
       enddo
     enddo
   enddo
   write(unit=40,fmt='(a3)') 'X'
   close(unit=40,status="keep")
 endif
 end subroutine





 subroutine testFUNC(u,kappa,nrhs,coact,bc,vecbl,vecblinv,myid,ntmqcd,nn,ldiv, &
                     nms,lvbc,lvcl,ib,lbd,iblv,rwdir,MRT,MRT2,numprocs)
 real(kind=KR),    intent(inout), dimension(:,:,:,:,:) :: u
 real(kind=KR),    intent(in),    dimension(:)         :: kappa
 integer(kind=KI), intent(in)                          :: nrhs
 character(len=*), intent(in),    dimension(:)         :: rwdir
 integer(kind=KI), intent(in),    dimension(:)         :: bc, nms
 integer(kind=KI), intent(in)                          :: myid, MRT,MRT2,numprocs
 real(kind=KR),    intent(in),    dimension(:,:,:)     :: coact
 integer(kind=KI), intent(in),    dimension(:,:)       :: vecbl, vecblinv
 integer(kind=KI), intent(in),    dimension(:,:)       :: nn, iblv
 logical,          intent(in),    dimension(:)         :: ldiv
 integer(kind=KI), intent(in),    dimension(:,:,:)     :: lvbc, lvcl
 integer(kind=KI), intent(in),    dimension(:,:,:,:)   :: ib
 logical,          intent(in),    dimension(:,:)       :: lbd
 integer(kind=KI), intent(in)                          :: ntmqcd

 integer(kind=KI)                             :: ieo,ibleo,idag,ierr
 integer(kind=KI)                             :: eig,irhs,gblclr
 real(kind=KR2),  dimension(2)                :: tmp1, tmp2

 real(kind=KR2),  dimension(6,ntotal,4,2,8)   :: tmpE,tmpO
 real(kind=KR2),  dimension(6,ntotal,4,2,8)   :: gxe,gxo,gxetmp,gxotmp!BS last2
 real(kind=KR2),  dimension(6,ntotal,4,2,8)   :: z2e, z2o,gz2e,gz2o,be,bo
 real(kind=KR2),  dimension(6,ntotal,4,2,8,1) :: xeTMP,xoTMP

 ! GMRES Variables
 integer(kind=KI), dimension(2)                :: nGMRES,nGMRESeig
 integer(kind=KI)                              :: itermin, iflag, inverter, isignal, mvp
 real(kind=KR),    dimension(2,kmaxGMRES)      :: gdr
 real(kind=KR)                                 :: resmax
 real(kind=KR),    dimension(18,nvhalf,8,2,16) :: GeeGooinv

! integer(kind=KI), parameter :: nop=5, nsty=6, eigstart=1, eigstep=5, ksub=6
! integer(kind=KI), parameter :: nop=9, nsty=9, eigstart=1, eigstep=5, ksub=6
 integer(kind=KI), parameter :: nop=9, nsty=10, eigstart=1, eigstep=1, ksub=6!BS
! integer(kind=KI), parameter :: nop=9, nsty=10, eigstart=1, eigstep=5, ksub=6

 real(kind=KR),   dimension(6,ntotal,4,2,8,kmaxGMRES) :: evectorRevnPRIME, evectorRoddPRIME
 real(kind=KR),   dimension(6,ntotal,4,2,8,kmaxGMRES) :: evectorRevn, evectorRodd
 real(kind=KR),   dimension(6,ntotal,4,2,8,kmaxGMRES) :: evectorLevn, evectorLodd
 real(kind=KR),   dimension(kmaxGMRES)                :: eigValExpandPRIME
 real(kind=KR),   dimension(2,kmaxGMRES)              :: eigValExpand
 integer(kind=KI)                                     :: numEVPRIME,numEV
 real(kind=KR),   dimension(2,kmaxGMRES)              :: inv_lam, inv_xi
 real(kind=KR),   dimension(2,kmaxGMRES)              :: inv_lamPRIME, inv_xiPRIME
 real(kind=KR2),  dimension(2,ksub,kmaxGMRES)         ::xiinv,xiinvPRIME,xiinvPOLY
 
 real(kind=KR2),  dimension(2,nt,ksub,kmaxGMRES,nsty,nop) :: Jvev,Jave,Jtotal
 !BS trying to get rid of error generated from already existing .dat files
  integer(kind=KI) :: i,isub,id,isite,icri
  integer(kind=KI) :: exists
character(len=128) :: noisefile
integer(kind=KI)   :: noiseunitnum

  real(kind=KR),    dimension(6,ntotal,4,2,8) ::   vecsrc, vecsrc1,temper !BS
  real(kind=KR)                               ::  fac1, fac2



!print *,"myid=",myid

                                                 
if (myid==0) then
!BS inquire if the file exists and create new if not for ns.dat
        inquire(file=trim(rwdir(myid+1))//"ns.dat", exist=exists)
  if (.not. exists) then

open(unit=21,file=trim(rwdir(myid+1))//"ns.dat",status="new",     &
        action="write",form="formatted")
   close(unit=21,status="keep")
end if


!BS inquire if the file exists and create new if not for ps.dat

 inquire(file=trim(rwdir(myid+1))//"ps.dat", exist=exists)
  if (.not. exists) then

open(unit=22,file=trim(rwdir(myid+1))//"ps.dat",status="new",     &
        action="write",form="formatted")
   close(unit=22,status="keep")
end if

!BS inquire if the file exists and create new if not for es.dat
inquire(file=trim(rwdir(myid+1))//"es.dat", exist=exists)

if (.not.exists) then

open(unit=23,file=trim(rwdir(myid+1))//"es.dat",status="new",     &
        action="write",form="formatted")
   close(unit=23,status="keep")
end if

!BS inquire if the file exists and create new if not for hfes.dat
inquire(file=trim(rwdir(myid+1))//"hfes.dat", exist=exists)
  if (.not. exists) then

open(unit=24,file=trim(rwdir(myid+1))//"hfes.dat",status="new",     &
        action="write",form="formatted")
   close(unit=24,status="keep")
end if
 
!BS inquire if the file exists and create new if not for esps.dat
inquire(file=trim(rwdir(myid+1))//"esps.dat", exist=exists)
  if (.not. exists) then

open(unit=25,file=trim(rwdir(myid+1))//"esps.dat",status="new",     &
        action="write",form="formatted")
   close(unit=25,status="keep")
end if
 

!BS inquire if the file exists and create new if not for hfesps.dat
inquire(file=trim(rwdir(myid+1))//"hfesps.dat", exist=exists)
  if (.not. exists) then

open(unit=26,file=trim(rwdir(myid+1))//"hfesps.dat",status="new",     &
        action="write",form="formatted")
   close(unit=26,status="keep")
end if
 

!inquire if the file exists and create new if not for pp.dat
inquire(file=trim(rwdir(myid+1))//"pp.dat", exist=exists)
  if (.not. exists) then

open(unit=27,file=trim(rwdir(myid+1))//"pp.dat",status="new",     &
        action="write",form="formatted")
   close(unit=27,status="keep")
end if
 
!inquire if the file exists and create new if not for pp4.dat
inquire(file=trim(rwdir(myid+1))//"pp4.dat", exist=exists)
  if (.not. exists) then

open(unit=28,file=trim(rwdir(myid+1))//"pp4.dat",status="new",     &
        action="write",form="formatted")
   close(unit=28,status="keep")
end if

!inquire if the file exists and create new if not for pp7.dat
inquire(file=trim(rwdir(myid+1))//"pp7.dat", exist=exists)
  if (.not. exists) then

open(unit=29,file=trim(rwdir(myid+1))//"pp7.dat",status="new",     &
        action="write",form="formatted")
   close(unit=29,status="keep")
end if
!BS added this hfespoly
 inquire(file=trim(rwdir(myid+1))//"hfespoly.dat", exist=exists)
  if (.not. exists) then
         print *, "File does not exist. Creating it."

open(unit=30,file=trim(rwdir(myid+1))//"hfespoly.dat",status="new",     &
        action="write",form="formatted")
   close(unit=30,status="keep")
end if
!BS added this 
! inquire(file=trim(rwdir(myid+1))//"NSPSHFESPS.dat", exist=exists)
 ! if (.not. exists) then
  !       print *, "File does not exist. Creating it."

!open(unit=31,file=trim(rwdir(myid+1))//"NSPSHFESPS.dat",status="new",     &
 !       action="write",form="formatted")
  ! close(unit=31,status="keep")
!end if
!BS added this 
 !inquire(file=trim(rwdir(myid+1))//"NSPSHFES.dat", exist=exists)
  !if (.not. exists) then
   !      print *, "File does not exist. Creating it."

!open(unit=32,file=trim(rwdir(myid+1))//"NSPSHFES.dat",status="new",     &
 !       action="write",form="formatted")
  ! close(unit=32,status="keep")
!end if

!BS Now change the status to old. Actually there was new inplace of old. All
!these inquire were put there  to stop the crashing

  open(unit=21,file=trim(rwdir(myid+1))//"ns.dat",status="old",     &

        action="write",form="formatted")
   close(unit=21,status="keep")
   open(unit=22,file=trim(rwdir(myid+1))//"ps.dat",status="old",     &
        action="write",form="formatted")
   close(unit=22,status="keep")
   open(unit=23,file=trim(rwdir(myid+1))//"es.dat",status="old",     &
        action="write",form="formatted")
   close(unit=23,status="keep")
   open(unit=24,file=trim(rwdir(myid+1))//"hfes.dat",status="old",   &
        action="write",form="formatted")
   close(unit=24,status="keep")
   open(unit=25,file=trim(rwdir(myid+1))//"esps.dat",status="old",   &
        action="write",form="formatted")
   close(unit=25,status="keep")
   open(unit=26,file=trim(rwdir(myid+1))//"hfesps.dat",status="old", &
        action="write",form="formatted")
   close(unit=26,status="keep")
   open(unit=27,file=trim(rwdir(myid+1))//"pp.dat",status="old", &
        action="write",form="formatted")
   close(unit=27,status="keep")
   open(unit=28,file=trim(rwdir(myid+1))//"pp4.dat",status="old", &
        action="write",form="formatted")
   close(unit=28,status="keep")
   open(unit=29,file=trim(rwdir(myid+1))//"pp7.dat",status="old", &
        action="write",form="formatted")
   close(unit=29,status="keep")
endif

 ! Idenfity number of Eigmodes to try to calculate
 nGMRESeig(1) =200  !400
 nGMRESeig(2) =160  !400!BS
 !  nGMRESeig(1) = 70 !400
  !nGMRESeig(2) = 50  !350
  ! nGMRESeig(2) = 160  !350!BS
 !BS nGMRESeig(2) = 50  !350
!testF
 !!!!       --- Gather eigenmode information ---    
 ! eigenmode info for M' = gamma5 M                     
 if (myid==0) then
 print *,'calling eigmodesofMprime'
 endif                                                   
 call eigmodesofMprime(evectorRevnPRIME,evectorRoddPRIME,eigValExpandPRIME,numEVPRIME, &
                       u,kappa,coact,bc,vecbl,vecblinv,myid,nn,ldiv,nms,lvbc,          &
                       lvcl,ib,lbd,iblv,rwdir,MRT,MRT2,nGMRESeig)





!!!!!!!!********BS Printing of eigenvector starts here .. To print just change
!TRUE to FALSE
!!!!!!!!!!!!!!!!*********************!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


if (.false.) then

 do eig = 1, numEVPRIME
   if (eig==1) then


      be = evectorRevnPRIME(:,:,:,:,:,eig)
      bo = evectorRoddPRIME(:,:,:,:,:,eig)
     if (myid==0) then
         inquire(file=trim('/data/barals/qqcd/scratch/evecR1.dat'), exist=exists)
         if (.not. exists) then
            print *, "File does not exist. Creating evecR1.dat."
            open(unit=38,file=trim('/data/barals/qqcd/scratch/evecR1.dat'),status="new",&
                action="write",form="formatted")
            close(unit=38,status="keep")
         endif



       open(unit=38,file=trim('/data/barals/qqcd/scratch/evecR1.dat'),status="old",&
         action="write",form="formatted",position="append")


          do ibleo = 1,8
            do ieo = 1,2
              do id = 1,4
                do isite = 1,nvhalf
                  do icri = 1,6

          write(unit=38,fmt="(i5,a5,i5,a5,i5,a5,i5,a5,i5,a5,F13.10,a5,F13.10)") ibleo,"",ieo,"",id,"",isite,"",icri,"",be(icri,isite,id,ieo,ibleo),"",bo(icri,isite,id,ieo,ibleo)

                  enddo ! icri
                enddo ! isite
              enddo ! id
            enddo ! ieo
          enddo ! ibleo




       close(unit=38,status="keep")
     endif!myid
  endif!eig
 enddo!eig
endif!TRUE FALSE


!!!!!!!***********BS printing of eigenvectorends here **********!!!!!!!!!!!!!!
!!!!!!!*********************************************************!!!!!!!!!!!!!!

  if(.false.) then
    do eig=1,numEVPRIME
     if(eig<10) then
        write(noisefile,"(a4,i1)")"eig-",eig
     elseif(eig<100) then
        write(noisefile,"(a4,i2)")"eig-",eig
     else
        write(noisefile,"(a4,i3)")"eig-",eig
     endif
         if (myid==0) then
             inquire(file=trim(rwdir(myid+1))//trim(noisefile)//".LOG",exist=exists)
             if (.not. exists) then
                open(unit=80,file=trim(rwdir(myid+1))//trim(noisefile)//".LOG",status="new",&
                action="write",form="formatted")
                close(unit=80,status="keep")
             endif
         endif



      be = evectorRevnPRIME(:,:,:,:,:,eig)
      bo = evectorRoddPRIME(:,:,:,:,:,eig)

      call printferm(be,bo,rwdir,noisefile,myid,iblv,vecblinv,MRT2)!BS printing 2/4/2016

    enddo
  endif










 ! eigenmode info for M
 call eigmodesofM(evectorRevn,evectorRodd,evectorLevn, evectorLodd,eigValExpand,numEV, &
                  u,kappa,coact,bc,vecbl,vecblinv,myid,nn,ldiv,nms,lvbc,               &
                  lvcl,ib,lbd,iblv,rwdir,MRT,MRT2,nGMRESeig)
 ! make 1/xi 

 call xigenerator(xiinv,numEV,eigstart,eigstep,u,kappa(1),                             &
                  evectorRevn,evectorRodd,evectorLevn,evectorLodd,                     &
                  coact,bc,vecbl,vecblinv,myid,nn,ldiv,nms,lvbc,ib,lbd,iblv,MRT2,1_KI)

 ! make 1/xiprime 
 call xigenerator(xiinvPRIME,numEVPRIME,eigstart,eigstep,u,kappa(1),                   &
                  evectorRevnPRIME,evectorRoddPRIME,evectorRevnPRIME,evectorRoddPRIME, &
                  coact,bc,vecbl,vecblinv,myid,nn,ldiv,nms,lvbc,ib,lbd,iblv,MRT2,2_KI)
! call xipolygenerator(xiinvPRIME,numEVPRIME,eigstart,eigstep,u,kappa(1),&
 !                 evectorRevnPRIME,evectorRoddPRIME,evectorRevnPRIME,evectorRoddPRIME,&
  !                z2e,z2o,coact,bc,vecbl,vecblinv,myid,nn,ldiv,nms,lvbc,ib,lbd,iblv,MRT2,1_KI)

!BS ********printing of xiinvprime.dat starts **********

if (.false.) then 
 if (myid==0) then

       inquire(file=trim(rwdir(myid+1))//"xiinvprime.dat", exist=exists)
           if (.not. exists) then
               open(unit=35,file=trim(rwdir(myid+1))//"xiinvprime.dat",status="new",     &
                   action="write",form="formatted")
               close(unit=35,status="keep")
           endif



           open(unit=35,file=trim(rwdir(myid+1))//"xiinvprime.dat",status="old",     &
                  action="write",form="formatted",position="append")
                        do eig =1,numEVPRIME
                           do isub=1,6
                              write(unit=35,fmt="(i5,a5,i5,a5,F13.10,a5,F13.10)") eig,"   ",isub,"   ",xiinvPRIME(1,isub,eig),"  ",xiinvPRIME(2,isub,eig)
                           enddo
                        enddo     
 
            close(unit=35,status="keep")
 endif

endif!true or false
!BS***********printing of xiinvprime ends************

 ! --  output -- !
 if (myid==0) then
   write(*,*) 'NUMEV:         ',numEV
   write(*,*) 'NUMEVPRIME:    ',numEVPRIME
   write(*,*) 'eigstart:      ',eigstart
   write(*,*) 'eigstep:       ',eigstep
 !  write(*,*)
 !  write(*,*)
 !  write(*,*)
 !  write(*,*) 'xiinv:       ',xiinv!BS
 !  write(*,*)
 !  write(*,*)
 !  write(*,*)
 !  write(*,*) 'xiinvprime:       ',xiinVPRIME!BS
   write(*,*)
   write(*,*)
   write(*,*) '--- Eigenvalues ---'
   do eig=1,numEV
     write(*,*) eig,eigValExpand(:,eig)
   enddo
   write(*,*)
   write(*,*)
   write(*,*) '--- Eigenvalues Gamma5 ---'
   do eig=1,numEV
     write(*,*) eig,eigValExpandPRIME(eig)
   enddo

 endif
 !!!!!!!!!!!!!!!!!


 ! GMRES Parameters
!BS nGMRES(1) = 60 ! must be smaller than nGMRESeig(1)
 !nGMRES(1) = 20!BS ! must be smaller than nGMRESeig(1)
 nGMRES(1) = 80!BS ! must be smaller than nGMRESeig(1)
 !nGMRES(2) = 10 ! must be smaller than nGMRESeig(2)
 nGMRES(2) = 40 ! must be smaller than nGMRESeig(2)
 !nGMRES(2) = 40 ! must be smaller than nGMRESeig(2)
 itermin = 10
 iflag = -1   ! Wilson
 GeeGooinv = 0.0_KR2
 resmax = 1e-7_KR2
! inverter = 6 ! use gmresdr/gmresproj
 inverter = 6!use bicgstab,abdou's suggestion,it seems psc doesn't converge
 isignal = 2  ! used to identify gmresproj, 1st rhs is called in eigmodeofM
 mvp = 0
 gdr = 0.0_KR2

 ! sets rhs for gmresdr in fermprop
 call setrhsnoise(2)

 ! Calculate scalar muliplier term needed in computation
 do eig=1,numEV
   call oneover(eigValExpand(:,eig),inv_lam(:,eig))
 enddo
 do eig=1,numEVPRIME
   inv_lamPRIME(1,eig) = 1.0_KR2/eigValExpandPRIME(eig)
   inv_lamPRIME(2,eig) = 0.0_KR2
 enddo

irhs=1
call z4source(z2e,z2o,nvhalf,myid)
 call xipolygenerator(xiinvPOLY,numEVPRIME,eigstart,eigstep,u,kappa(1),&
                  evectorRevnPRIME,evectorRoddPRIME,evectorRevnPRIME,evectorRoddPRIME,&
                  z2e,z2o,coact,bc,vecbl,vecblinv,myid,nn,ldiv,nms,lvbc,ib,lbd,iblv,MRT2,2_KI)

!BS ************printing of xiinvpoly starts********8888

if (.false.) then 


 if (myid==0) then

       inquire(file=trim(rwdir(myid+1))//"xiinvpoly.dat", exist=exists)
           if (.not. exists) then
               open(unit=41,file=trim(rwdir(myid+1))//"xiinvpoly.dat",status="new",&
                   action="write",form="formatted")
               close(unit=41,status="keep")
           endif



           open(unit=41,file=trim(rwdir(myid+1))//"xiinvpoly.dat",status="old",&
                  action="write",form="formatted",position="append")
                        do eig =1,numEVPRIME
                           do isub=1,6
                              write(unit=41,fmt="(i5,a5,i5,a5,F17.10,a5,F17.10)")eig,"   ",isub,"   ",xiinvPOLY(1,isub,eig),"  ",xiinvPOLY(2,isub,eig)
                           enddo
                        enddo

            close(unit=41,status="keep")
 endif




endif !true or false

!BS **********printing of xiinvPOLY ends ***************BS




 !!!!       --- Calculate Tr(\tilde{M}) ---      !!!!
 ! Since this term does not contribute to the error analysis it will not be implemented 
 ! at the moment. For a full calculation the appropriate function calls need to be implemented
 ! here. 
 Jvev = 0.0_KR2 ! vev is commented out since only error bars are being looked at ******
psignal = 0!the switch used to determine whether calculate the coeffecients for            !other right hands

 !drhs=1,o irhs=1,numnoises !200BS
 do irhs=1,nrhs !200
   Jtotal = 0.0_KR2
   Jave = 0.0_KR2
     ! print *,"checking if it passes line of creating noisefile111"

   ! generate noise vector
   call z4source(z2e,z2o,nvhalf,myid)






  if(.false.) then
     if(irhs<10) then
        write(noisefile,"(a6,i1)")"noise-",irhs
     elseif(irhs<100) then
        write(noisefile,"(a6,i2)")"noise-",irhs
     else
        write(noisefile,"(a6,i3)")"noise-",irhs
     endif
         if (myid==0) then
             inquire(file=trim(rwdir(myid+1))//trim(noisefile)//".LOG", exist=exists)
             if (.not. exists) then
                open(unit=80,file=trim(rwdir(myid+1))//trim(noisefile)//".LOG",status="new",     &
        action="write",form="formatted")
                close(unit=80,status="keep")
             endif
         endif
  endif 


 !call printferm(z2e,z2o,rwdir,noisefile,myid,iblv,vecblinv,MRT2)!BS printing 2/4/2016
   ! solve for solution vector (gxe,gxo): g5Mx'=z => Mx'=g5z
   ! form modified noise vector (gze,gzo) = gamma5*z 
!****Andy is here.I want to test the Mx=z case.***********!
!   do ieo=1,2
!     do ibleo=1,8
!       call gammamult(z2e, z2o, gz2e, gz2o, 5,ieo,ibleo)
!     enddo
!   enddo
!******************************end**************************!







   call fermprop(rwdir,z2e,z2o,iflag,kappa,0.0_KR2,inverter,coact,bc,resmax,       &
                 itermin,0.0_KR2,u,nGMRES(1),nGMRES(2),xeTMP,xoTMP,myid,nn,ldiv,nms, &
                 lvbc,lvcl,ib,lbd,iblv,vecbl,vecblinv,0,0,MRT,MRT2)
   gxe = xeTMP(:,:,:,:,:,1)!Mx=z,not g5Mx'=z
   gxo = xoTMP(:,:,:,:,:,1)


!to do the g5M case,change z2 to gz2!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! ---- DEBUG CODE ---- !
   ! Set to .true. to check residual of each solution explicitly
   ! NOTE: This is being done for the gamma5Mx'=z problem
   ! NOTE: This is written for the Wilson action only
   idag=0!abdou's change
   if (.true.) then


     call Hsingle(tmpO,u,gxe,idag,coact,bc,2_KI, &
                  vecbl,vecblinv,myid,nn,ldiv,nms,lvbc,ib,lbd,iblv,MRT)

      call Hsingle(tmpE,u,gxo,idag,coact,bc,1_KI, &
                   vecbl,vecblinv,myid,nn,ldiv,nms,lvbc,ib,lbd,iblv,MRT)

     tmpE(:6,:nvhalf,:4,:2,:8) = gxe(:6,:nvhalf,:4,:2,:8) - kappa(1)*tmpE(:6,:nvhalf,:4,:2,:8)
     tmpO(:6,:nvhalf,:4,:2,:8) = gxo(:6,:nvhalf,:4,:2,:8) - kappa(1)*tmpO(:6,:nvhalf,:4,:2,:8)

     tmpE = tmpE - z2e!Andy is here.To check g5M,change z2 to gz2.
     tmpO = tmpO - z2o

     call vecdot(tmpE,tmpE,tmp1,MRT2)
     call vecdot(tmpO,tmpO,tmp2,MRT2)
     tmp1 = tmp1 + tmp2
     call cmpxsqrt(tmp1, tmp2)

    !BS if (myid==0) write(*,*) 'residual(inv=0): ',tmp2(1)
     call vecdot(z2e,z2e,tmp1,MRT2)!Andy is here.
     call vecdot(z2o,z2o,tmp2,MRT2)!Andy is here.
     tmp1 = tmp1 + tmp2
     call cmpxsqrt(tmp1,tmp2)

     !BS if (myid==0) write(*,*) 'initial residual: ',tmp2(1)
   endif

   
   ! Eig Subtraction Corrections (all possible eigen combinations/operators)
   call generalaverage(Jave,numEV,numEVPRIME,eigstart,eigstep,u,z2e,z2o,xeTMP(:,:,:,:,:,1),xoTMP(:,:,:,:,:,1), &
                       kappa,evectorRevn,evectorRodd,evectorLevn, evectorLodd,inv_lam,                         &
                       evectorRevnPRIME, evectorRoddPRIME, inv_lamPRIME,xiinv, xiinvPRIME,xiinvPOLY,           &
                       coact,bc,vecbl,vecblinv,myid,nn,ldiv,nms,lvbc,ib,lbd,iblv,rwdir,MRT,MRT2)

   Jtotal = Jave + Jvev
   Jtotal = Jtotal/real(nx*ny*nz,KR2) ! Normalize by volume of latticea
   ! ####Is it correct to divide by the spatial volume every time a contribution from a noise vector gets computed?### -AA


   call vevoutput(Jtotal,numEV,numEVPRIME,eigstart,eigstep,rwdir,myid)

   psignal = psignal + 1!if want to caculate coeff for each rhs,comment out this
 enddo

 !call MPI_FINALIZE(ierr)
 end subroutine testFUNC



 subroutine eigmodesofMprime(evectorRevnPRIME,evectorRoddPRIME,eigValExpandPRIME,numEVPRIME, &
                             u,kappa,coact,bc,vecbl,vecblinv,myid,nn,ldiv,nms,lvbc,          &
                             lvcl,ib,lbd,iblv,rwdir,MRT,MRT2,nGMRES)
   use shift
   real(kind=KR),    intent(out), dimension(:,:,:,:,:,:) :: evectorRevnPRIME, evectorRoddPRIME
   real(kind=KR),    intent(out), dimension(:)           :: eigValExpandPRIME
   integer(kind=KI), intent(out)                         :: numEVPRIME
   real(kind=KR),    intent(inout), dimension(:,:,:,:,:) :: u
   real(kind=KR),    intent(in),    dimension(:)         :: kappa
   !! --- Identify number of eigenmodes to be calculated --- !!
   integer(kind=KI), intent(in),    dimension(:)         :: nGMRES
   

   character(len=*), intent(in),    dimension(:)         :: rwdir
   integer(kind=KI), intent(in),    dimension(:)         :: bc, nms
   integer(kind=KI), intent(in)                          :: myid, MRT,MRT2
   real(kind=KR),    intent(in),    dimension(:,:,:)     :: coact
   integer(kind=KI), intent(in),    dimension(:,:)       :: vecbl, vecblinv
   integer(kind=KI), intent(in),    dimension(:,:)       :: nn, iblv
   logical,          intent(in),    dimension(:)         :: ldiv
   integer(kind=KI), intent(in),    dimension(:,:,:)     :: lvbc, lvcl
   integer(kind=KI), intent(in),    dimension(:,:,:,:)   :: ib
   logical,          intent(in),    dimension(:,:)       :: lbd


   real(kind=KR), dimension(kmaxGMRES)          :: evalueNORM
   integer(kind=KI)                             :: irhs,j,ii,nv,gblclr
   integer(kind=KI)                             :: idag, ierr
   real(kind=KR2),  dimension(2)                :: dottmp1,dottmp2
   real(kind=KR2),  dimension(6,ntotal,4,2,8)   :: xe, xo
   real(kind=KR2),  dimension(6,ntotal,4,2,8)   :: b,tmpE,bepart,tmpEE !BS
   real(kind=KR2),  dimension(6,ntotal,4,2,8)   :: z2e, z2o, tmpO,tmpOO !BS 
   real(kind=KR2),  dimension(6,ntotal,4,2,8) :: xeTMP,xoTMP
   real(kind=KR2),  dimension(6,ntotal,4,2,8,2*kmaxGMRES) :: evectorODD,evectorODDL
   real(kind=KR2),  dimension(6,ntotal,4,2,8,2*kmaxGMRES) :: evectorL,evectorR
   real(kind=KR2),  dimension(2,kmaxGMRES)      :: evalueTMP
   real(kind=KR2),  dimension(2,2*kmaxGMRES)    :: eigValExpand,eigValExpandLEFT
   real(kind=KR2),  dimension(2,kmaxGMRES*2)    :: scalarMultiplier
   integer,         dimension(kmaxGMRES)        :: ind,ind2
   real(kind=KR2),  dimension(2)                :: tmp1, tmp2, tmp3, tmp4

   ! GMRES Variables
   integer(kind=KI)                           :: itermin, iflag, inverter
   real(kind=KR)                              :: resmax,constt
   integer(kind=KI)                           :: itercount, LUcount
   real(kind=KR), dimension(18,nvhalf,8,2,16) :: GeeGooinv

   real(kind=KR), dimension(6,ntotal,4,2,8)   :: bopart
   integer(kind=KI)                           :: iaa,ibb,icc,eigtest,eig,eigKEEP,eigVecKEEP,id
   integer(kind=KI)                           :: ieo,iop,sty,it,styindx,ibleo
   integer(kind=KI)                           :: numEV,evStrt,evEnd
   real(kind=KR)                              :: tmp,tmpReal,tmpImag
   integer(kind=KI)                           :: dim1, dim2,dim3,dim4
   integer(kind=KI)                           :: fullruns,styprint



if (myid==0) then
 print *, "ntotal aaj k din=", ntotal
endif
   ! initialize/declare vairable parameters
   idag = 0; ! guarentees Mx not Mdagx
   itermin = 10
   iflag = -1  ! Wilson
   GeeGooinv = 0.0_KR2
    if (myid==0) then
 print *,'just entered  eigmodesofMprime'
 endif

   ! a throw away right hand side
   call z4source(z2e,z2o,nvhalf,myid)

   ! Solve Mx=b
   xe(:6,:ntotal,:4,:2,:8) = 0.0_KR2  ! Initial Guess
   xo(:6,:ntotal,:4,:2,:8) = 0.0_KR2  ! Initial Guess
   call gmresdr5EIG(rwdir,z2e,z2o,xe,xo,nGMRES,resmax,itermin,itercount,u,GeeGooinv, &
                 iflag,idag,kappa,coact,bc,vecbl,vecblinv,myid,nn,ldiv,nms, &
                 lvbc,ib,lbd,iblv,MRT,MRT2)

   if (myid ==0) then
            print *,"just done with gmresdr5EIG"!BS
         !   print *,"evalue(1,3)=",evalue(1,ii)!BS
        endif



  ! **Debug** Use to verify residual
   if (.false.) then
     !!!!! check residual r=Mx-b
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!*****BS added for flipping gamma5 multiplication********
 ! gamma5*M
     do ieo=1,2
       do ibleo=1,8
         call gammamult(xe, xo, tmpE, tmpO, 5,ieo,ibleo)
       enddo
     enddo

 ! First build H_oe * v_e.
     gblclr = 2
     call Hsingle(tmpOO,u,tmpE,idag,coact,bc,gblclr,vecbl,vecblinv,myid,nn,ldiv, &
                  nms,lvbc,ib,lbd,iblv,MRT)
     ! Next build H_eo * v_o.
     gblclr = 1
     call Hsingle(tmpEE,u,tmpO,idag,coact,bc,gblclr,vecbl,vecblinv,myid,nn,ldiv, &
                  nms,lvbc,ib,lbd,iblv,MRT)
     tmpE = tmpE - kappa(1) * tmpEE
     tmpO = tmpO - kappa(1) * tmpOO

!****BS add ends here**********
!**********BS start commenting out wrong gamma5 multiply**********
     ! First build H_oe * v_e.
    ! gblclr = 2
    ! call Hsingle(tmpO,u,xe,idag,coact,bc,gblclr,vecbl,vecblinv,myid,nn,ldiv, &
     !             nms,lvbc,ib,lbd,iblv,MRT)
     ! Next build H_eo * v_o.
    ! gblclr = 1
    ! call Hsingle(tmpE,u,xo,idag,coact,bc,gblclr,vecbl,vecblinv,myid,nn,ldiv, &
     
!             nms,lvbc,ib,lbd,iblv,MRT)
 !    xe = xe - kappa(1) * tmpE 
  !   xo = xo - kappa(1) * tmpO


  ! gamma5*M 
   !  do ieo=1,2
    !   do ibleo=1,8
     !    call gammamult(xe, xo, tmpE, tmpO, 5,ieo,ibleo)
      ! enddo
    ! enddo
!*********BS comment out ends here*********

     tmpE = tmpE - z2e
     tmpO = tmpO - z2o

     call vecdot(tmpE,tmpE,tmp1,MRT2)
     call vecdot(tmpO,tmpO,tmp2,MRT2)
     tmp1 = tmp1 + tmp2
     call cmpxsqrt(tmp1, tmp2)
!BS if (myid==0) write(*,*) 'false is called '
     if (myid==0) write(*,*) 'wrong printing false residual: ',tmp2(1)
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   endif  !if false


   !!! check for acceptable eigenvectors
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   numEVprime = 0  ! count number of acceptable modes
   do ii=1,nGMRES(2)
     ! normalize eigenvectors
     call vecdot(evectore(:,:,:,:,:,ii),evectore(:,:,:,:,:,ii),tmp1,MRT2)
     call vecdot(evectoro(:,:,:,:,:,ii),evectoro(:,:,:,:,:,ii),tmp2,MRT2)
     tmp1 = tmp1 + tmp2
     tmp1(1) = sqrt(tmp1(1))
     evectore(:,:,:,:,:,ii) = (1.0_KR2/tmp1(1)) * evectore(:,:,:,:,:,ii)
     evectoro(:,:,:,:,:,ii) = (1.0_KR2/tmp1(1)) * evectoro(:,:,:,:,:,ii)

      ! if ((myid ==0).AND.(ii==2)) then
       !   print *," evectore(:,:,:,:,:,2)=", evectore(:,:,:,:,:,2)!BS
      ! endif



     ! holder for ease of use
     xe = evectore(:,:,:,:,:,ii)
     xo = evectoro(:,:,:,:,:,ii)
!*****BS start flipping the gamma5 calculation 
 do ieo=1,2
       do ibleo=1,8
         call gammamult(xe,xo, tmpE, tmpO,5,ieo,ibleo)
       enddo
     enddo
  ! Perform Mx (recall this is explicitly WILSON)
     ! First build H_oe * v_e.
     gblclr = 2
     call Hsingle(tmpOO,u,tmpE,idag,coact,bc,gblclr,vecbl,vecblinv,myid,nn,ldiv, &
                  nms,lvbc,ib,lbd,iblv,MRT)
     ! Next build H_eo * v_o.
     gblclr = 1
     call Hsingle(tmpEE,u,tmpO,idag,coact,bc,gblclr,vecbl,vecblinv,myid,nn,ldiv, &
                  nms,lvbc,ib,lbd,iblv,MRT)
     xeTMP = tmpE - kappa(1) * tmpEE
     xoTMP = tmpO - kappa(1) * tmpOO
!if ((myid ==0).AND.(ii==1)) then
!print *,"tmpEE1=",tmpEE!BS
!endif
!if ((myid ==0).AND.(ii==2)) then
!print *,"tmpEE2=",tmpEE!BS
!endif

  ! tmp = gamma5*M*x - eval*x
     tmpE = xeTMP - (evalue(1,ii)*xe)!BS should it be xe or gamma5xe
     tmpO = xoTMP - (evalue(1,ii)*xo)
       ! if (myid ==0) then
        !    print *,"evalue(1,1)=",evalue(1,1)!BS
         !   print *,"evalue(1,3)=",evalue(1,ii)!BS
       ! endif
!******BS end flipping gamma5 multiplication

! **********BS comment out starts here**********
     
    ! Perform Mx (recall this is explicitly WILSON)
     ! First build H_oe * v_e.
    ! gblclr = 2
    ! call Hsingle(tmpO,u,xe,idag,coact,bc,gblclr,vecbl,vecblinv,myid,nn,ldiv, &
     !             nms,lvbc,ib,lbd,iblv,MRT)
     ! Next build H_eo * v_o.
    ! gblclr = 1
    ! call Hsingle(tmpE,u,xo,idag,coact,bc,gblclr,vecbl,vecblinv,myid,nn,ldiv, &
    !              nms,lvbc,ib,lbd,iblv,MRT)
   !  tmpE = xe - kappa(1) * tmpE 
   !  tmpO = xo - kappa(1) * tmpO

     ! multiply by gamma5: gamma5*M*x
   !  do ieo=1,2  
    !   do ibleo=1,8
     !    call gammamult(tmpE, tmpO, xeTMP, xoTMP,5,ieo,ibleo)
     !  enddo
    ! enddo
     ! tmp = gamma5*M*x - eval*x
   !  tmpE = xeTMP - (evalue(1,ii)*xe)!BS should it be xe or gamma5xe
   !  tmpO = xoTMP - (evalue(1,ii)*xo)

! *********BS comment out ends here

     ! calculate norm of tmp2=|gamma5*M*x-eval*x|
     call vecdot(tmpE,tmpE,tmp1,MRT2)
     call vecdot(tmpO,tmpO,tmp2,MRT2)
     tmp1 = tmp1 + tmp2
     call cmpxsqrt(tmp1, tmp2)
       !if (myid ==0) then
       ! print *,"tmp2(1)=",tmp2(1)
       !endif
      !if (myid ==0) then
      !  print *,"numEVprime=",numEVprime
     !endif


     ! check for accuracy !BS dropped the criterion from 8 to 5
     if (tmp2(1) < 1e-5) then
       numEVprime = numEVprime + 1
       ind(numEVprime) = ii  ! list of accurate modes
     endif
   enddo
     !if (myid ==0) then 
      !  print *,"numEVprime=",numEVprime
     !endif 
   ! generate a list of values by absolute magnitude
   do eig=1,numEVprime
     evalueNORM(eig) = abs(evalue(1,ind(eig)))
   enddo

   ! sort modes by smallest-largest evalue
   call sort(evalueNORM,numEVprime,ind2)

   ! place data in final holders
   do eig=1,numEVprime
     eigvalExpandPRIME(eig) = evalue(1,ind(ind2(eig)))  ! note: only real part

     evectorRevnPRIME(:,:,:,:,:,eig) = evectore(:,:,:,:,:,ind(ind2(eig)))
     evectorRoddPRIME(:,:,:,:,:,eig) = evectoro(:,:,:,:,:,ind(ind2(eig)))
   enddo
    
 !BS if (myid==0) write(*,*) 'eigvalExpandPRIME1:',eigvalExpandPRIME(:,:,:,:,:,eig)
   !BS  if (myid==0) write(*,*) 'evectorRevnPRIME2:',evectorRevnPRIME(:,:,:,:,:,2)

   !{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{
   ! DEBUGGING PURPOSES
   ! *set to .true. to test L/R orthoganlity for all full eigenvectors 
   ! *only non-zero values are outputed to screen which would represent normality or
   !      non orthoganality
   ! *output is sent to standard out
   if (.true.) then
     ! check LR orthogonality
     if (myid==0) write(*,*) 'L/R Orthoganality test for Mprime'
     do j=1,numEVPRIME
       do eig=1,numEVPRIME
         call vecdot(evectorRevnPRIME(:,:,:,:,:,eig),evectorRevnPRIME(:,:,:,:,:,j),tmp1,MRT)
         call vecdot(evectorRoddPRIME(:,:,:,:,:,eig),evectorRoddPRIME(:,:,:,:,:,j),tmp2,MRT)
         tmp3(:2) = tmp1(:2) + tmp2(:2)
         if (myid==0) then
           if (abs(tmp3(1)) > 1e-6_KR2 .or. abs(tmp3(2)) > 1e-6_KR2) write(*,*) eig,j,tmp3
         endif
       enddo
     enddo
   endif
   !}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}


   !!! output acceptable eigenvalues
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   if (myid==0) then
     open(unit=12,file=trim(rwdir(myid+1))//"EIGEN_VALS.LOG",status="old",&
          action="write",form="formatted",position="append")
     write(unit=12,fmt='(a60)')"-------GMRES-DR5EIG (MATLAB ver) (full computed of Mprime)---------"
     write(unit=12,fmt='(a28,i4)') "               numEVprime: ",numEVprime
     do eig=1,numEVPRIME
       write(unit=12,fmt='(i5,a5,f25.15)')  eig,"     ",eigvalExpandPRIME(eig)
     enddo
     close(unit=12,status="keep")
   endif
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 end subroutine eigmodesofMprime







! Victor's Contribution
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 ! traceCorrection(:2,:styles,:operators,:eigs,:rhs)
 !   styles: 1) only small eigenmodes
 !           2) only large eigenmodes
 !           3) both small/large eigenmodes
 !           0) no subtraction
 !  operator: 1-4) vector <psibar gamma_mu psi>
 !            5) pseudoscalar <psibar gamma5 psi>
 !            6-9) point-split <psibar+a_mu (1+gamma_mu)U_mu^dag psi> - <psi (1-gamma_mu)U_mu psibar+a_mu>
 !            10) scalar <psibar psi>
 !  eigs: 1-numEV; all possible eigenmodes
 !        eigstart modes used then every eigstep after that is used
 !  rhs: 1-numnoises; all possible rhs done 
 subroutine eigdiscon(u,kappa,nrhs,coact,bc,vecbl,vecblinv,myid,ntmqcd,nn,ldiv, &
                     nms,lvbc,lvcl,ib,lbd,iblv,rwdir,MRT,MRT2,numprocs)
   use shift
   real(kind=KR),    intent(inout), dimension(:,:,:,:,:) :: u
   real(kind=KR),    intent(in),    dimension(:)         :: kappa
   integer(kind=KI), intent(in)                          :: nrhs

   character(len=*), intent(in),    dimension(:)         :: rwdir
   integer(kind=KI), intent(in),    dimension(:)         :: bc, nms
   integer(kind=KI), intent(in)                          :: myid, MRT,MRT2,numprocs
   real(kind=KR),    intent(in),    dimension(:,:,:)     :: coact
   integer(kind=KI), intent(in),    dimension(:,:)       :: vecbl, vecblinv
   integer(kind=KI), intent(in),    dimension(:,:)       :: nn, iblv
   logical,          intent(in),    dimension(:)         :: ldiv
   integer(kind=KI), intent(in),    dimension(:,:,:)     :: lvbc, lvcl
   integer(kind=KI), intent(in),    dimension(:,:,:,:)   :: ib
   logical,          intent(in),    dimension(:,:)       :: lbd

   integer(kind=KI), intent(in) :: ntmqcd

   integer(kind=KI)                             :: irhs,j,ii,nv
   integer(kind=KI)                             :: idag
   real(kind=KR2),  dimension(2)                :: dottmp1,dottmp2
   real(kind=KR2),  dimension(6,ntotal,4,2,8)   :: xe, xo
   real(kind=KR2),  dimension(6,ntotal,4,2,8)   :: b,tmpE,bepart
   real(kind=KR2),  dimension(6,ntotal,4,2,8)   :: z2e, z2o, tmpO
   real(kind=KR2),  dimension(6,ntotal,4,2,8,1) :: xeTMP,xoTMP
   real(kind=KR2),  dimension(6,ntotal,4,2,8,2*kmaxGMRES) :: evectorODD,evectorODDL
   real(kind=KR2),  dimension(6,ntotal,4,2,8,2*kmaxGMRES) :: evectorL,evectorR
   real(kind=KR2),  dimension(6,ntotal,4,2,8,kmaxGMRES) :: evectorRprime, evectorRODDprime
   real(kind=KR2),  dimension(2,kmaxGMRES)      :: evalueTMP
   real(kind=KR2),  dimension(2,2*kmaxGMRES)    :: eigValExpand,eigValExpandLEFT
   real(kind=KR2),  dimension(kmaxGMRES)        :: eigValExpandprime
   real(kind=KR2),  dimension(2,kmaxGMRES*2)    :: scalarMultiplier
   integer,         dimension(kmaxGMRES)        :: ind,ind2
   real(kind=KR2),  dimension(kmaxGMRES)        :: evalueNORM
   real(kind=KR2),  dimension(2)                :: tmp1, tmp2, tmp3, tmp4

!   integer(kind=KI), parameter :: nop=5, ksub=6
   integer(kind=KI), parameter :: nop=9, ksub=6
   real(kind=KR2),   dimension(2,nt,nop,4,kmaxGMRES,ksub)   :: Jvev,Jave,Jtotal


   ! GMRES Variables
   integer(kind=KI), dimension(2)             :: nGMRES
   integer(kind=KI)                           :: itermin, iflag, inverter
   real(kind=KR)                              :: resmax,constt
   integer(kind=KI)                           :: itercount, LUcount
   real(kind=KR), dimension(18,nvhalf,8,2,16) :: GeeGooinv

   real(kind=KR), dimension(6,ntotal,4,2,8)   :: bopart
   integer(kind=KI)                           :: iaa,ibb,icc,eigtest,eig,eigKEEP,eigVecKEEP,id
   integer(kind=KI)                           :: ieo,iop,sty,it,styindx
   integer(kind=KI)                           :: numEV,numEVprime,evStrt,evEnd
   real(kind=KR)                              :: tmp,tmpReal,tmpImag
   integer(kind=KI)                           :: dim1, dim2,dim3,dim4
   integer(kind=KI)                           :: fullruns,styprint

   integer(kind=KI), parameter        :: eigstart = 1_KI
   integer(kind=KI), parameter        :: eigstep = 10_KI
   integer(kind=KI), parameter        :: numruns = 1_KI
   integer(kind=KI), parameter        :: eigsty = 1_KI

   real(kind=KR2), dimension(2,ksub,kmaxGMRES) :: xiinv

   ! initialize/declare vairable parameters
   idag = 0; ! guarentees Mx not Mdagx
   itermin = 10
   iflag = -1  ! Wilson
   GeeGooinv = 0.0_KR2

   if (myid == 0) then
     write(*,*) 'kappa:',kappa(1)
   endif

   ! generate eigen modes of the full problem
 !  call eigmodesofM(evectorR,evectorODD,evectorL,evectorODDL,eigValExpand,numEV, &
 !                   u,kappa,coact,bc,vecbl,vecblinv,myid,nn,ldiv,nms,lvbc,       &
 !                   lvcl,ib,lbd,iblv,rwdir,MRT,MRT2)

   ! generate eigen modes of the full problem M'=gamma5*M (Wilson only)
!   call eigmodesofMprime(evectorRprime,evectorRODDprime,eigValExpandprime,numEVprime, &
!                         u,kappa,coact,bc,vecbl,vecblinv,myid,nn,ldiv,nms,lvbc,       &
!                         lvcl,ib,lbd,iblv,rwdir,MRT,MRT2)

   ! calculate "scalarMultiplier" which is 1/lam
   do eig=1,numEV+numEV
     call oneover(eigValExpand(:2,eig),scalarMultiplier(:2,eig))
   enddo

   ! calculate 1/xi, the representation of the low modes of Mpert
!   call xigenerator(xiinv,numEV,eigstart,eigstep,u,kappa(1), &
!                    evectorR,evectorODD,evectorL,evectorODDL, &
!                    coact,bc,vecbl,vecblinv,myid,nn,ldiv,nms,lvbc,ib,lbd, &
!                    iblv,MRT2)

   ! Prime files for output
   if (myid==0) then
     open(unit=20,file=trim(rwdir(myid+1))//"ns.dat",status="new", &
          action="write",form="formatted")
     close(unit=20,status="keep")
     open(unit=21,file=trim(rwdir(myid+1))//"es.dat",status="new", &
          action="write",form="formatted")
     close(unit=21,status="keep")
     open(unit=22,file=trim(rwdir(myid+1))//"ps.dat",status="new", &
          action="write",form="formatted")
     close(unit=22,status="keep")
     open(unit=23,file=trim(rwdir(myid+1))//"espsc.dat",status="new", &
          action="write",form="formatted")
     close(unit=23,status="keep")
   endif

!!!!!!! START ALGORITHM  !!!!!!!!!!!!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
   ! *** compute tr(OMtilde) ***
   Jvev = 0.0_KR2 ! vev is commented out since only error bars are being looked at ******
 !  call generalvev(Jvev,numEV,eigstart,eigstep,eigsty,evectorR,evectorODD, evectorL,evectorODDL,&
 !                  scalarMultiplier,xiinv,kappa,ntmqcd,vecbl,vecblinv,&
 !                  u,rwdir,coact,nms,bc,nn,iblv,ldiv,lvbc,ib,lbd,myid,MRT,MRT2,numprocs)

   ! INITIATE RUNS
   ! set parameters for GMRESDR
   !nGMRES(1) = 100
   !nGMRES(2) = 60
   nGMRES(1) = 50
   nGMRES(2) = 30

   resmax = 1E-14_KR2
   inverter = 6 ! use non-matlab-gmresdr and proj from fermprop

   ! Full runs: 1 full run will have nrhs comptued
   do fullruns=1,numruns
     ! initialize data holder per full run
     Jtotal = 0.0_KR2

     ! all right hand sides
     do irhs=1,nrhs
       if (myid==0) then
        open(unit=8,file=trim(rwdir(myid+1))//"CFGSPROPS.LOG",action="write", &
           form="formatted",status="old",position="append")
        write(unit=8,fmt="(a5,i3,a9,i3)") "run: ",fullruns,"   irhs: ",irhs
        close(unit=8,status="keep")
       endif

       ! sets rhs for gmresdr in fermprop
       call setrhsnoise(irhs + fullruns - 1)

       ! generate z2 noise
       !call z2source(z2e,z2o,nvhalf)
       call z4source(z2e,z2o,nvhalf,myid)
       call printz2noise(z2e,z2o,nvhalf,numprocs,rwdir,MRT2,myid)

       ! solve Mx=z2
       call fermprop(rwdir,z2e,z2o,iflag,kappa,0.0_KR2,inverter,coact,bc,resmax,   &
                     itermin,0.0_KR2,u,nGMRES(1),nGMRES(2),xeTMP,xoTMP,myid,nn,ldiv,nms, &
                     lvbc,lvcl,ib,lbd,iblv,vecbl,vecblinv,0,0,MRT,MRT2)

       ! Eig Subtraction Corrections (all possible eigen combinations/operators)
   !    call generalaverage(Jave,numEV,eigstart,eigstep,u,z2e,z2o,xeTMP(:,:,:,:,:,1),xoTMP(:,:,:,:,:,1), &
   !                    kappa,evectorR,evectorODD,evectorL,&
   !                    evectorODDL,scalarMultiplier,xiinv, &
   !                    coact,bc,vecbl,vecblinv,myid,nn,ldiv,nms,lvbc,ib,lbd, &
   !                    iblv,rwdir,MRT2)

       Jtotal = Jtotal + Jave
     enddo
     Jtotal = Jtotal/real(nrhs,KR2) + Jvev
     Jtotal = Jtotal/real(nx*ny*nz,KR2) ! Normalize by volume of latticea

     ! Output data for review
   !  call vevoutput(Jtotal,numEV,eigstart,eigstep,rwdir,myid)
   enddo ! full runs

   end subroutine eigdiscon




 ! Calculate the NS term in the Tr(O M^-1) calculation
 ! Tr(O M^-1) = z.Ox           (x =M^-1 z)         ESHFES = 1
 ! Tr(O M^-1) = z.gamma5 Ox'   (x'=M^-1 gamma5 z)  ESHFES = else
 ! Jtime(rc,nt,mom,iop)
 subroutine nsaverage(Jtime,u,z2e,z2o,xe,xo,kappa,coact,bc,vecbl, &
                      vecblinv,myid,nn,ldiv,nms,lvbc,ib,lbd,      &
                      iblv,MRT2,ESHFES)
 real(kind=KR),    intent(out),   dimension(:,:,:,:)   :: Jtime
 real(kind=KR),    intent(inout), dimension(:,:,:,:,:) :: u
 real(kind=KR),    intent(in)                          :: kappa
 real(kind=KR),    intent(inout), dimension(:,:,:,:,:) :: z2e, z2o
 real(kind=KR),    intent(inout), dimension(:,:,:,:,:) :: xe, xo
 real(kind=KR),    intent(in),    dimension(:,:,:)     :: coact
 integer(kind=KI), intent(in),    dimension(:)         :: bc, nms
 integer(kind=KI), intent(in),    dimension(:,:)       :: vecbl, vecblinv
 integer(kind=KI), intent(in),    dimension(:,:)       :: nn, iblv
 logical,          intent(in),    dimension(:)         :: ldiv
 integer(kind=KI), intent(in),    dimension(:,:,:)     :: lvbc
 integer(kind=KI), intent(in),    dimension(:,:,:,:)   :: ib
 logical,          intent(in),    dimension(:,:)       :: lbd
 integer(kind=KI), intent(in)                          :: myid, MRT2, ESHFES 

! integer(kind=KI), parameter                         :: nmom=5, nop=5
 integer(kind=KI), parameter                         :: nmom=5, nop=9

 integer(kind=KI)                                    :: ieo,ibleo

 real(kind=KR2),   dimension(2,nvhalf,nop)           :: Je, Jo

 real(kind=KR2),   dimension(nvhalf,2,16,nmom-1,nop) :: momfac

 ! initialize data holders
 Jtime  = 0.0_KR2
 momfac = 0.0_KR2

 ! ---------------- Non-Subtraction ----------------------------!
 if (ESHFES == 1) then ! non PRIME method (M no gamma5)
   ! apply operators, via time slices
   do ibleo = 1,8
     do ieo = 1,2
       call loopops(0_KI,xe,xo,u,z2e,z2o,Je,Jo,ieo,ibleo,bc,vecbl,   &
                    vecblinv,myid,nn,ldiv,nms,lvbc,ib,lbd,iblv,MRT2)
       call spacesum(Jtime(:,:,:,:),Je,Jo,ieo,ibleo,momfac,nmom,nop,    &
                     myid,vecbl,vecblinv,MRT2)
     enddo ! ieo
   enddo ! ibleo
 else ! PRIME method (M' = gamma5 M)
   ! apply operators, via time slices
   do ibleo = 1,8
     do ieo = 1,2

      !BS 5/8/015 call eigloopops(0_KI,xe,xo,u,z2e,z2o,Je,Jo,ieo,ibleo,bc,vecbl,   &
       !                vecblinv,myid,nn,ldiv,nms,lvbc,ib,lbd,iblv,MRT2)
       call eigloopops_original(0_KI,xe,xo,u,z2e,z2o,Je,Jo,ieo,ibleo,bc,vecbl,   &
                       vecblinv,myid,nn,ldiv,nms,lvbc,ib,lbd,iblv,MRT2)


       call spacesum(Jtime(:,:,:,:),Je,Jo,ieo,ibleo,momfac,nmom,nop,    &
                     myid,vecbl,vecblinv,MRT2)
     enddo ! ieo
   enddo ! ibleo
 endif                  

 ! kappa normalization for J_mu current
 Jtime(:,:,:,1) = kappa*Jtime(:,:,:,1)
 Jtime(:,:,:,2) = kappa*Jtime(:,:,:,2)
 Jtime(:,:,:,3) = kappa*Jtime(:,:,:,3)
 Jtime(:,:,:,4) = kappa*Jtime(:,:,:,4)
 end subroutine nsaverage
   

!------------------------------------------------------------------------------
subroutine nsaverage_abdou(Jtime,u,z2e,z2o,xe,xo,kappa,coact,bc,vecbl, &
                      vecblinv,myid,nn,ldiv,nms,lvbc,ib,lbd,      &
                      iblv,MRT2,ESHFES)
 real(kind=KR),    intent(out),   dimension(:,:,:,:)   :: Jtime
 real(kind=KR),    intent(inout), dimension(:,:,:,:,:) :: u
 real(kind=KR),    intent(in)                          :: kappa
 real(kind=KR),    intent(inout), dimension(:,:,:,:,:) :: z2e, z2o
 real(kind=KR),    intent(inout), dimension(:,:,:,:,:) :: xe, xo
 real(kind=KR),    intent(in),    dimension(:,:,:)     :: coact
 integer(kind=KI), intent(in),    dimension(:)         :: bc, nms
 integer(kind=KI), intent(in),    dimension(:,:)       :: vecbl, vecblinv
 integer(kind=KI), intent(in),    dimension(:,:)       :: nn, iblv
 logical,          intent(in),    dimension(:)         :: ldiv
 integer(kind=KI), intent(in),    dimension(:,:,:)     :: lvbc
 integer(kind=KI), intent(in),    dimension(:,:,:,:)   :: ib
 logical,          intent(in),    dimension(:,:)       :: lbd
 integer(kind=KI), intent(in)                          :: myid, MRT2, ESHFES 

! integer(kind=KI), parameter                         :: nmom=5, nop=5
 integer(kind=KI), parameter                         :: nmom=5, nop=9

 integer(kind=KI)                                    :: ieo,ibleo

 real(kind=KR2),   dimension(2,nvhalf,nop)           :: Je, Jo

 real(kind=KR2),   dimension(nvhalf,2,16,nmom-1,nop) :: momfac

 ! initialize data holders
 Jtime  = 0.0_KR2
 momfac = 0.0_KR2

 ! ---------------- Non-Subtraction ----------------------------!
 if (ESHFES == 1) then ! non PRIME method (M no gamma5)
   ! apply operators, via time slices
   do ibleo = 1,8
     do ieo = 1,2
       call loopops(0_KI,xe,xo,u,z2e,z2o,Je,Jo,ieo,ibleo,bc,vecbl,   &
                    vecblinv,myid,nn,ldiv,nms,lvbc,ib,lbd,iblv,MRT2)
       call spacesum(Jtime(:,:,:,:),Je,Jo,ieo,ibleo,momfac,nmom,nop,    &
                     myid,vecbl,vecblinv,MRT2)
     enddo ! ieo
   enddo ! ibleo
 else ! PRIME method (M' = gamma5 M)
   ! apply operators, via time slices
   do ibleo = 1,8
     do ieo = 1,2
       call eigloopops_abdou(0_KI,xe,xo,u,z2e,z2o,Je,Jo,ieo,ibleo,bc,vecbl,   &
                       vecblinv,myid,nn,ldiv,nms,lvbc,ib,lbd,iblv,MRT2)
       call spacesum(Jtime(:,:,:,:),Je,Jo,ieo,ibleo,momfac,nmom,nop,    &
                     myid,vecbl,vecblinv,MRT2)
     enddo ! ieo
   enddo ! ibleo
 endif                  

 ! kappa normalization for J_mu current
 Jtime(:,:,:,1) = kappa*Jtime(:,:,:,1)
 Jtime(:,:,:,2) = kappa*Jtime(:,:,:,2)
 Jtime(:,:,:,3) = kappa*Jtime(:,:,:,3)
 Jtime(:,:,:,4) = kappa*Jtime(:,:,:,4)
 end subroutine nsaverage_abdou

 !------------------------------------------------------------------------------ 






 ! Purpose of this is simply to calculate averages of operators using the
 ! eigensubtraction technique. Calculation is CORRECTION ONLY and does not have
 ! the ns part in the subtraction
 ! INPUT:
 !   nsub is the number of results, having different orders of subtraction,
 !        to be considered (including the case of order=0).
 !   z2e() contains the source vector on even (gblclr=1) lattice sites.
 !   z2o() contains the source vector on odd (gblclr=2) lattice sites.
 !        expected size: be(6,ntotal,4,2,8), bo(6,ntotal,4,2,8)
 !                       where the first entry is real/imaginary and colour
 !                       and the 3rd entry is the Dirac index
 !                       and the other three entries give the lattice site.
 !   xe() contains the sink vector on even (gblclr=1) lattice sites.
 !   xo() contains the sink vector on odd (gblclr=2) lattice sites.
 !        expected size: xe(6,ntotal,4,2,8), xo(6,ntotal,4,2,8)
 !                       where the first entry is real/imaginary and colour
 !                       and the 3rd entry is the Dirac index
 !                       and the other three entries give the lattice site.
 !   ESHFES is a flag to determin whether to do Eigenspectrum Subtraction or
 !          Hermitian Forced Eigenspectrum Subtraction
 !          ESHFES = 1 : ES Method    ESHFES = 2 : HFES Method
 ! OUTPUT:
 !   Jtime(iri,it,eig,mom,iop) are the averaged operators.
 !        it=1..nt is the time slice on the full (multi-process) lattice.
 subroutine eigaverage(Jtime,numEV,eigstart,eigstep,u,z2e,z2o,xe,xo,kappa,     &
                       evecReven,evecRodd,evecLeven,evecLodd,scalarMultiplier, &
                       coact,bc,vecbl,vecblinv,myid,nn,ldiv,nms,lvbc,ib,lbd,   &
                       iblv,MRT2,ESHFES)
 real(kind=KR2),    intent(out),   dimension(:,:,:,:,:)    :: Jtime
 integer(kind=KI), intent(in)                             :: numEV
 integer(kind=KI), intent(in)                             :: eigstart, eigstep
 real(kind=KR),    intent(inout), dimension(:,:,:,:,:)    :: u
 real(kind=KR),    intent(in)                             :: kappa
 real(kind=KR),    intent(inout), dimension(:,:,:,:,:)    :: z2e, z2o
 real(kind=KR),    intent(inout), dimension(:,:,:,:,:)    :: xe, xo
 real(kind=KR2),   intent(in),    dimension(:,:,:,:,:,:)  :: evecReven, evecRodd
 real(kind=KR2),   intent(in),    dimension(:,:,:,:,:,:)  :: evecLeven, evecLodd
 real(kind=KR2),   intent(in),    dimension(:,:)          :: scalarMultiplier
 real(kind=KR),    intent(in),    dimension(:,:,:)        :: coact
 integer(kind=KI), intent(in),    dimension(:)            :: bc, nms
 integer(kind=KI), intent(in),    dimension(:,:)          :: vecbl, vecblinv
 integer(kind=KI), intent(in),    dimension(:,:)          :: nn, iblv
 logical,          intent(in),    dimension(:)            :: ldiv
 integer(kind=KI), intent(in),    dimension(:,:,:)        :: lvbc
 integer(kind=KI), intent(in),    dimension(:,:,:,:)      :: ib
 logical,          intent(in),    dimension(:,:)          :: lbd
 integer(kind=KI), intent(in)                             :: myid, MRT2, ESHFES 

 real(kind=KR2), dimension(2)              :: tmp1, tmp2
 integer(kind=KI)                          :: eig,ieo,ibleo,icri
 real(kind=KR2), dimension(6,ntotal,4,2,8) :: xprimeE,xprimeO
 real(kind=KR2), dimension(6,ntotal,4,2,8) :: eigcor_q_even,eigcor_q_odd

! integer(kind=KI), parameter                         ::  nmom=5, nop=5
 integer(kind=KI), parameter                         ::  nmom=5, nop=9
 real(kind=KR2),   dimension(2,nvhalf,nop)           :: Je, Jo
 real(kind=KR2),   dimension(nvhalf,2,16,nmom-1,nop) :: momfac


 ! initialize data holders
 eigcor_q_even = 0.0_KR2
 eigcor_q_odd  = 0.0_KR2
 Jtime  = 0.0_KR2
 momfac = 0.0_KR2

 ! ---------------- Eigen Subtraction -------------------------!
 do eig=1,numEV ! only do small modes
   ! tmp2 = L^dag (dot) z2
   call vecdot(evecLeven(:,:,:,:,:,eig),z2e,tmp1,MRT2)
   call vecdot( evecLodd(:,:,:,:,:,eig),z2o,tmp2,MRT2)
   tmp2(:2) = tmp1(:2) + tmp2(:2)

   ! tmp1 = (1/lam) * [L^dag (dot) z2] = scalarMultiplier(:2,eig) * tmp2
   call cmpxmult(tmp2,scalarMultiplier(:2,eig),tmp1)

   ! eigcor_q = eigcor_q + (R * 1/lam * [L^dat (dot) z2]) = eigcor_q + (eigcor_q * tmp1)
   do icri = 1,5,2
     eigcor_q_even(icri  ,:nvhalf,:4,:2,:8) = eigcor_q_even(icri  ,:nvhalf,:4,:2,:8)          &
                                            + tmp1(1)*evecReven(icri  ,:nvhalf,:4,:2,:8,eig)  &
                                            - tmp1(2)*evecReven(icri+1,:nvhalf,:4,:2,:8,eig)
     eigcor_q_even(icri+1,:nvhalf,:4,:2,:8) = eigcor_q_even(icri+1,:nvhalf,:4,:2,:8)          &
                                            + tmp1(1)*evecReven(icri+1,:nvhalf,:4,:2,:8,eig)  &
                                            + tmp1(2)*evecReven(icri  ,:nvhalf,:4,:2,:8,eig)  
  
     eigcor_q_odd(icri  ,:nvhalf,:4,:2,:8) = eigcor_q_odd(icri  ,:nvhalf,:4,:2,:8)            &
                                           + tmp1(1)*evecRodd(icri  ,:nvhalf,:4,:2,:8,eig)    &
                                           - tmp1(2)*evecRodd(icri+1,:nvhalf,:4,:2,:8,eig)  
     eigcor_q_odd(icri+1,:nvhalf,:4,:2,:8) = eigcor_q_odd(icri+1,:nvhalf,:4,:2,:8)            &
                                           + tmp1(1)*evecRodd(icri+1,:nvhalf,:4,:2,:8,eig)    &
                                           + tmp1(2)*evecRodd(icri  ,:nvhalf,:4,:2,:8,eig)
   enddo ! icri

   ! only compute result for specific eigmodes
   if ( (eig == eigstart) .or. (eig > eigstart .and. 0 == mod(eig-eigstart,eigstep)) .or. (eig == numEV) ) then
     xprimeE(:6,:nvhalf,:4,:2,:8) = eigcor_q_even(:6,:nvhalf,:4,:2,:8)
     xprimeO(:6,:nvhalf,:4,:2,:8) =  eigcor_q_odd(:6,:nvhalf,:4,:2,:8)
     if (ESHFES  == 1) then ! ES
       ! apply operators, via time slices
       do ibleo = 1,8
         do ieo = 1,2
           call loopops(0_KI,xprimeE,xprimeO,u,z2e,z2o,Je,Jo,ieo,ibleo,bc,vecbl, &
                        vecblinv,myid,nn,ldiv,nms,lvbc,ib,lbd,iblv,MRT2)
           call spacesum(Jtime(:,:,eig,:,:),Je,Jo,ieo,ibleo,momfac,nmom,nop,     &
                         myid,vecbl,vecblinv,MRT2)
         enddo ! ieo
       enddo ! ibleo
     else ! HFES
       ! apply operators, via time slices
       do ibleo = 1,8
         do ieo = 1,2
           call eigloopops_original(0_KI,xprimeE,xprimeO,u,z2e,z2o,Je,Jo,ieo,ibleo,bc,vecbl, &
                           vecblinv,myid,nn,ldiv,nms,lvbc,ib,lbd,iblv,MRT2)
           call spacesum(Jtime(:,:,eig,:,:),Je,Jo,ieo,ibleo,momfac,nmom,nop,        &
                         myid,vecbl,vecblinv,MRT2)
         enddo ! ieo
       enddo ! ibleo
     endif
   endif ! specific eig modes
 enddo ! eig

 ! kappa normalization for J_mu current
 do eig=1,numEV
   Jtime(:,:,eig,:,1) = kappa*Jtime(:,:,eig,:,1)
   Jtime(:,:,eig,:,2) = kappa*Jtime(:,:,eig,:,2)
   Jtime(:,:,eig,:,3) = kappa*Jtime(:,:,eig,:,3)
   Jtime(:,:,eig,:,4) = kappa*Jtime(:,:,eig,:,4)
 enddo
 end subroutine eigaverage












      subroutine eigspacesum(opsout,opsine,opsino,ieo,ibleo,myid, &
                     vecbl,vecblinv,MRT2)
      ! Compute the observables per timeslice.
      ! NOTE: it=1 is uniquely defined by          1<=isite<=  factor and 1<=ibl<= 8
      !       it=2 is uniquely defined by          1<=isite<=  factor and 9<=ibl<=16
      !       it=3 is uniquely defined by   factor+1<=isite<=2*factor and 1<=ibl<= 8
      !       it=4 is uniquely defined by   factor+1<=isite<=2*factor and 9<=ibl<=16
      !       it=5 is uniquely defined by 2*factor+1<=isite<=3*factor and 1<=ibl<= 8
      !       it=6 is uniquely defined by 2*factor+1<=isite<=3*factor and 9<=ibl<=16
      !       etc, where factor = 2*nvhalf*npt/nt.
      real(kind=KR),    intent(inout), dimension(:,:,:)   :: opsout
      real(kind=KR),    intent(in),    dimension(:,:,:)     :: opsine, opsino
      integer(kind=KI), intent(in)                          :: ieo, ibleo,myid, MRT2
      integer(kind=KI), intent(in),    dimension(:,:)       :: vecbl,vecblinv
 
      integer(kind=KI), parameter :: nop=10 ! number of operators

      integer(kind=KI) :: isite, ibl, jbl, jbleo, nsite, sitemin
      integer(kind=KI) :: sitemax, it, itbit, ierr, iop, n, iri, itconstant
      integer(kind=KI), dimension(4)             :: np, ip

      real(kind=KR2),   dimension(2,nt,nop) :: opsbit
      real(kind=KR2),   dimension(2,nt,nop) :: opscorr, opssum

      ! Sum the data according to timestep.
      ibl = vecbl(1,ibleo)
      jbl = vecbl(2,ibleo)
      jbleo = vecblinv(2,jbl)
      opsbit = 0.0_KR
      np(1) = npx
      np(2) = npy
      np(3) = npz
      np(4) = npt
      call atoc(myid,np,ip)
      itconstant = ip(4)*nt/npt - 1 + (ibleo-1)/4
      nsite = 2*nvhalf*npt/nt
      sitemin = 1
      sitemax = nsite
      do itbit = 2,nt/npt,2
        it = itbit + itconstant
        do isite = sitemin,sitemax
          do iop = 1,nop
            do iri = 1,2
              opsbit(iri,it,iop) = opsbit(iri,it,iop) &
                             + real(opsine(iri,isite,iop),KR2) &
                             + real(opsino(iri,isite,iop),KR2)
            enddo ! iri
          enddo ! iop
        enddo ! isite
        sitemin = sitemin + nsite
        sitemax = sitemax + nsite
      enddo ! itbit

      ! Collect the final results into opsout.
      opscorr = real(opsbit,KR)
      if (nps==1) then
        do iop = 1,nop
          do it = 1,nt
            do iri = 1,2
              opsout(iri,it,iop) = opsout(iri,it,iop) &
                                 + opscorr(iri,it,iop)
            enddo ! iri
          enddo ! it
        enddo ! iop
      else
        n = 2*nt*nop
        call MPI_REDUCE(opscorr(1,1,1),opssum(1,1,1),n,MRT2,MPI_SUM,0, &
                     MPI_COMM_WORLD,ierr)
        call MPI_BCAST(opssum(1,1,1),n,MRT2,0,MPI_COMM_WORLD,ierr)
        do iop = 1,nop
          do it = 1,nt
            do iri = 1,2
              opsout(iri,it,iop) = opsout(iri,it,iop) &
                                 + opssum(iri,it,iop)
            enddo ! iri
          enddo ! it
        enddo ! iop
      endif
 
      end subroutine eigspacesum

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!






     ! The inputted be and bo are only used if numnoises=0 (intended for debugging).
     ! Expected sizes: be(6,ntotal,4,2,8), bo(6,ntotal,4,2,8)
     subroutine discon(numnoises,be,bo,rwdir,iflag,kappa,cSW,coact,bc,resmax, &
                       itermin,vectype,u,myid,nn,ldiv,nms,lvbc,lvcl,ib,lbd,iblv, &
                       vecbl,vecblinv,hoption,MRT,MRT2)
     integer(kind=KI), intent(in) :: hoption,numnoises, iflag, itermin, myid,MRT, MRT2
     real(kind=KR),    intent(inout), dimension(:,:,:,:,:) :: be, bo, u
     character(len=*), intent(in), dimension(:)            :: rwdir
     real(kind=KR),    intent(in)                          :: kappa, cSW, resmax
     real(kind=KR),    intent(in), dimension(:,:,:)        :: coact
     integer(kind=KI), intent(in), dimension(:)            :: bc, nms
     integer(kind=KI), intent(inout)                       :: vectype
     integer(kind=KI), intent(in), dimension(:,:)   :: nn, iblv, vecbl, vecblinv
     logical,          intent(in), dimension(:)            :: ldiv
     integer(kind=KI), intent(in), dimension(:,:,:)        :: lvbc, lvcl
     integer(kind=KI), intent(in), dimension(:,:,:,:)      :: ib
     logical,          intent(in), dimension(:,:)          :: lbd

     !integer(kind=KI), parameter                         :: nsub=4, nmom=5, nop=5
     integer(kind=KI), parameter                         :: nsub=6, nmom=5, nop=5   !I think nop has to be 9 after we add the local vector current -AA
     integer(kind=KI)                                    :: is
     integer(kind=KI), dimension(3)                      :: io
     integer(kind=KI), dimension(4)                      :: np, ip
     real(kind=KR),    dimension(2,nt,nmom,nop)          :: Jvev
     real(kind=KR),    dimension(2,nt,nsub,nmom,nop)     :: Jave, Jtotal
     real(kind=KR),    dimension(2,nt,nsub,nmom,nop,nshifts) :: Jtemp
     real(kind=KR),    dimension(nvhalf,2,16,nmom-1,nop) :: momfac
     integer(kind=KI) :: inverter, iplug, inoise, iri, isub, imom, iop, dobndry
     real(kind=KR),    dimension(6,ntotal,4,2,8,1)       :: xe
     real(kind=KR),    dimension(6,nvhalf,4,2,8,1)       :: xo !why not ntotal -AA
     real(kind=KR),    dimension(2)                      :: kapvec
     real(kind=KR2)                                      :: omegaM3R
     integer(kind=KI) :: ierr,it

     integer(kind=KI),  dimension(0:9,4) :: ptsrc
     integer(kind=KI) :: icolor,isite,idirac,ieo,iblock

     integer(kind=KI) :: iruns
     integer(kind=KI), parameter :: numruns = 20_KI

     integer(kind=KI) :: i
     integer(kind=KI) :: exists


     if (myid==0) write(*,*) 'in discon (not twist) - me'
     ! Only a single kappa is used in the propagator.
     ! (Be sure NOT to use the multi-mass algorithm, since it requires TWO
     ! inversions for a random source.)
     kapvec(1) = kappa
     kapvec(2) = 0.0_KR
     omegaM3R  = 0.0_KR2
     inverter  = 5

     vectype = 80 !this should be input -AA
     iplug   = 40

     ! initalize
     Jave = 0.0_KR2
     Jtotal = 0.0_KR2
     Jvev = 0.0_KR2
 
     ! Define the (x,y,z) coordinates of the origin.
     io(1) = 1
     io(2) = 1
     io(3) = 1

     ! Identify the location of my process.
     np(1) = npx
     np(2) = npy
     np(3) = npz
     np(4) = npt
     call atoc(myid,np,ip)
     if (ip(4)==np(4)-1) then
       dobndry = 1
     else
       dobndry = 0
     endif

     ! Calculate the average plaquette for certain operators.
     call vev(Jvev,u,rwdir,io,kappa,dobndry,myid,nn,ldiv,nms,lvbc,ib,lbd, &
              iblv,MRT)


     if (myid == 0) then
 inquire(file=trim(rwdir(myid+1))//"pertvev.dat", exist=exists)
  if (.not. exists) then
         print *, "File does not exist. Creating it."

open(unit=21,file=trim(rwdir(myid+1))//"pertvev.dat",status="new",     &
        action="write",form="formatted")
   close(unit=21,status="keep")
end if

end if

       open(unit=21,file=trim(rwdir(myid+1))//"pertvev.dat",status="old",&
            action="write",form="formatted")
     
     do iruns=1,numruns
       Jtotal = 0.0_KR2
       Jave = 0.0_KR2

       do isub = 2,nsub
         Jtotal(:,:,isub,:,:) = Jvev(:,:,:,:)
       enddo ! isub

       ! Initialize.
       call momfacs(io,momfac,myid,iblv)

       ! Loop over noises. (numnoises=0 is intended for debugging purposes only.)
       if (numnoises==0) then
         call fermprop(rwdir,be,bo,iflag,kapvec,cSW,inverter,coact,bc,resmax, &
                       itermin,omegaM3R,u,vectype,iplug,xe,xo,myid,nn,ldiv,nms, &
                       lvbc,lvcl,ib,lbd,iblv,vecbl,vecblinv,0,hoption,MRT,MRT2)
         call average(Jave,nsub,nmom,nop,momfac,u,be,bo,xe(:,:,:,:,:,1), &
                      xo(:,:,:,:,:,1),kappa,dobndry,coact,bc,vecbl,vecblinv,myid, &
                      nn,ldiv,nms,lvbc,ib,lbd,iblv,rwdir,MRT)
         Jtotal = Jtotal + Jave
       else
         do inoise = 1,numnoises
           call z2source(be,bo,nvhalf)

           call printlog("Finished after z2", myid,rwdir)
           call setrhsnoise(iruns + inoise - 1)

           ! inverter=2 is hard coded for gmresdr-proj
           call fermprop(rwdir,be,bo,iflag,kapvec,cSW,2_KI,coact,bc,resmax, &
                         itermin,omegaM3R,u,vectype,iplug,xe,xo,myid,nn,ldiv,nms, &
                         lvbc,lvcl,ib,lbd,iblv,vecbl,vecblinv,0,hoption,MRT,MRT2)

           ! Javg IS cumulative.
           call average(Jave,nsub,nmom,nop,momfac,u,be,bo,xe(:,:,:,:,:,1), &
                        xo(:,:,:,:,:,1),kappa,dobndry,coact,bc,vecbl,vecblinv, &
                        myid,nn,ldiv,nms,lvbc,ib,lbd,iblv,rwdir,MRT)

           Jtotal = Jtotal + Jave

         enddo ! inoise
         Jtotal = Jtotal/real(numnoises,KR)

         ! TEMPORARY: rescale for comparison to Walter's normalization!!!!!!!!
         Jtotal = Jtotal/real(nx*ny*nz,KR)

         call printlog("finished x-current", myid,rwdir)

         ! Save the results.
         call printlog("about to print",myid,rwdir)
         if (myid==0) then
           do iop = 1,nop
             do imom = 1,1 !nmom
               do isub = 1,nsub
                 do it=1,nt
                   !BS write(unit=21,fmt="(a29,4i3,a3,es17.10,a4,es17.10)") &
                    !  "it,isub,imom,iop, discon:",it,isub,imom,iop,"   ",Jtotal(1,it,isub,imom,iop),"    ",Jtotal(2,it,isub,imom,iop)
                   write(unit=21,fmt="(a29,4i3,a3,es17.10,a4,es17.10)") &
                      " ",it,isub,imom,iop,"   ",Jtotal(1,it,isub,imom,iop),"    ",Jtotal(2,it,isub,imom,iop)
                 enddo ! it
               enddo ! isub
             enddo ! imom
           enddo ! iop
         endif
       endif ! (numnoises = 0)

     enddo !numruns
     if (myid==0) close(unit=21,status="keep")
     end subroutine discon

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      subroutine twistdiscon(numnoises,bedebug,bodebug,rwdir,iflag,kappa,mu,cSW,coact,bc,resmax, &
                   itermin,vectype,u,myid,nn,ldiv,nms,lvbc,lvcl,ib,lbd,iblv, &
                   vecbl,vecblinv,hoption,MRT,MRT2,numprocs,delta,shiftnum,ntmqcd,kGMRES,&
                   sweepnum,nsweep,icfgsave,ir)

! The inputted be and bo are only used if numnoises=0 (intended for debugging).
! Expected sizes: be(6,ntotal,4,2,8), bo(6,ntotal,4,2,8)


     integer(kind=KI), intent(in) :: hoption,numnoises, iflag, itermin,numprocs, myid, MRT, MRT2, ir
!    real(kind=KR),    intent(inout), dimension(:,:,:,:,:)             :: be, bo, u
     real(kind=KR),    intent(inout), dimension(6,ntotal,4,2,8)        :: bedebug, bodebug
     real(kind=KR),    intent(inout), dimension(18,ntotal,4,2,16)       :: u
     character(len=*), intent(in),    dimension(:)                     :: rwdir
     real(kind=KR),    intent(in)                                      :: kappa,mu,cSW,resmax
     real(kind=KR),    intent(in)                                      :: delta  
     real(kind=KR),    intent(in),    dimension(:,:,:)                 :: coact
     integer(kind=KI), intent(in),    dimension(:)                     :: bc, nms
     integer(kind=KI), intent(inout)                                   :: vectype
     integer(kind=KI), intent(in),    dimension(:,:)                   :: nn, iblv, vecbl, vecblinv
     logical,          intent(in),    dimension(:)                     :: ldiv
     integer(kind=KI), intent(in),    dimension(:,:,:)                 :: lvbc, lvcl
     integer(kind=KI), intent(in),    dimension(:,:,:,:)               :: ib
     integer(kind=KI), intent(in)                                      :: shiftnum,sweepnum,nsweep
     integer(kind=KI), intent(in)                                      :: ntmqcd
     integer(kind=KI), intent(in)                                      :: icfgsave
     integer(kind=KI), intent(in)                                      :: kGMRES
     logical,          intent(in),    dimension(:,:)                   :: lbd

     integer(kind=KI), parameter                                       :: nsub=6, nmom=5, nop=6
     !integer(kind=KI), parameter                                       :: nsub=6, nmom=5, nop=10
     integer(kind=KI), dimension(3)                                    :: io
     integer(kind=KI), dimension(4)                                    :: np, ip
     real(kind=KR),    dimension(2,nt,6,nmom,nop)                      :: Jtotal, Jvev, Jave,&
                                                                         Jconu,Jcond,Jnoise
     real(kind=KR),    dimension(2,6,nmom,nop)              :: Jp1, Jp2
     !real(kind=KR),    dimension(2,nt,6,nmom,nop,numnoises)            :: Jtemp
     real(kind=KR),    dimension(:,:,:,:,:,:), allocatable            :: Jtemp
     real(kind=KR),    dimension(2,nt,6,nmom,nop)                      :: Jtvev
     real(kind=KR),    dimension(nt,6,nmom,nop)                        :: sigma
     real(kind=KR),    dimension(nvhalf,2,16,nmom-1,nop)               :: momfac
     integer(kind=KI) :: inverter, iplug, inoise, iri, isub, imom, iop, dobndry, iii, usub, irloop,&
                        irloop1, irloop2
!    real(kind=KR),    dimension(6,ntotal,4,2,8,nshifts)               :: xe
!    real(kind=KR),    dimension(6,nvhalf,4,2,8,nshifts)               :: xo
     real(kind=KR),allocatable,dimension(:,:,:,:,:,:)               :: xe
     real(kind=KR),allocatable,dimension(:,:,:,:,:,:)               :: xo
     real(kind=KR),    dimension(2)                                    :: kapvec
!    real(kind=KR),    dimension(6,ntotal,4,2,8)                       :: be, bo, ze, zo
     real(kind=KR),allocatable,dimension(:,:,:,:,:)                    :: be, bo, ze, zo
     real(kind=KR2)                                                    :: omegaM3R
     real(kind=KR)                                                     :: noisecount, sweepcount
     real(kind=KR)                                                     :: fac
     character(len=128)                                                :: opfile,tempfile
     character(len=2)                                                  :: trailer
     character(len=3)                                                  ::looptrailer,irtrailer
     character(len=15)                                                 :: ltrailer
     real(kind=KR)                                                     :: time1,time2,time3,time4

     integer(kind=KI),  dimension(0:9,4) :: ptsrc
     integer(kind=KI) :: it, isweep , icolor, idirac, isite, ieo, iblock, icsrc, idsrc
     integer(kind=KI),  dimension(0:9,4) :: src
     logical          :: fileexists

! DEAN~ note that the integer tempvec is to debug the raw propagator
!       and needs to be removed.
     integer(kind=KI)                                                  :: tempvec, ierr,i,j,k,l,m,n,ittt

if (myid==0) write(*,*) 'in twistdiscon - me'

     !Timing for twistdiscon
     call cpu_time(time1)

     time3=time1
     if (myid==0) then
      open(unit=8,file=trim(rwdir(myid+1))//"CFGSPROPS.LOG",action="write", &
         form="formatted",status="old",position="append")
      write(unit=8,fmt="(a50,es17.10,a10)") "Started twistdiscon at:",time3,"secs"
      close(unit=8,status="keep")
     endif

     allocate(xe(6,ntotal,4,2,8,nshifts))
     allocate(xo(6,nvhalf,4,2,8,nshifts))
     allocate(be(6,ntotal,4,2,8))
     allocate(bo(6,ntotal,4,2,8))
     allocate(ze(6,ntotal,4,2,8))
     allocate(zo(6,ntotal,4,2,8))


     if(ir.eq.1) irloop=1
     if(ir.eq.-1) irloop=2

     allocate(Jtemp(2,nt,6,nmom,nop,numnoises))

! Only a single kappa is used in the propagator.
! (Be sure NOT to use the multi-mass algorithm, since it requires TWO
! inversions for a random source.)

! NOTE: We need the fermion prop to have both kappa and mu value since we have implemented
! twisted mass into the disconnected part of the calcualtion.
     kapvec(1) = kappa
     kapvec(2) = mu ! When using this in the Wilson case, mu is hardcoded to zero.

     if (.false.) then
       call printlog("WARNING :: Using DEBUG kappa, mu!!!", myid,rwdir)
      !kapvec(1) = .00001
       kapvec(1) = .154
       kapvec(2) = 0.00_KR
     endif ! .false.

     omegaM3R = 0.0_KR2
     inverter = 6
!   iplug = 0
      iplug = kGMRES 

      sigma = 0.0_KR

! Define the (x,y,z) coordinates of the origin.
    io(1) = 1
    io(2) = 1
    io(3) = 1


! Identify the location of my process.
     np(1) = npx
     np(2) = npy
     np(3) = npz
     np(4) = npt
     call atoc(myid,np,ip)
     if (ip(4)==np(4)-1) then
      dobndry = 0
     else
      dobndry = 1
     endif

! Calculate the average plaquette for certain operators.

! Dean - NOte this kappa needs to be changed to mtmqcd(1)-Did it!

      Jtotal = 0.0_KR
      Jvev   = 0.0_KR
      Jave   = 0.0_KR
      Jtemp  = 0.0_KR
      Jtvev  = 0.0_KR
      
! Don't nee the above calucaltion for dobndry for our twvev and twAverage subroutines.
! dobndry = 1 (fixbc=.true.)
      
      !dobndry = 1           !Original code (Abdou)
      if(bc(4)==0) then   ! Allows for periodic or antiperiodic boundry conditions. 
       dobndry=1          ! fixbc is true only if bc(4)=0. -Abdou
      else
       dobndry=0
      endif
      
! twvev is not cumulative in the sense of adding "subtraction levels". 
! For highest order subtraction we must combine it with the first-
! nontrivial subtraction level.

      call cpu_time(time3)
  
       if (myid==0) then
        open(unit=8,file=trim(rwdir(myid+1))//"CFGSPROPS.LOG",action="write", &
            form="formatted",status="old",position="append")
        write(unit=8,fmt="(a50,es17.10,a10)") "Started twvev at:",time3,"secs"
        close(unit=8,status="keep")
       endif
     call twvev(Jvev,delta,io,kapvec(1),u,dobndry,numprocs,MRT2,&
               rwdir,myid,nsub,nmom,nop,shiftnum,ntmqcd,ir)

      call cpu_time(time4)

      if (myid==0) then
      open(unit=8,file=trim(rwdir(myid+1))//"CFGSPROPS.LOG",action="write", &
          form="formatted",status="old",position="append")
      write(unit=8,fmt="(a50,es17.10,a10)") "Time taken by twvev:",time4-time3,"secs"
      close(unit=8,status="keep")
     endif



! This is to match Randy's normilization
      do iop = 1,nop
         if (iop /= 2 .and. iop /= 6) then
            Jvev(:,:,:,:,iop) = kapvec(1)*Jvev(:,:,:,:,iop)
         endif
      enddo ! iop

! Combine highest order subtraction with first order.
    if (nsub == 6) then
      do iop = 1,nop
         Jvev(:,:,6,:,iop) = Jvev(:,:,4,:,iop) + Jvev(:,:,6,:,iop)
      enddo ! iop
    endif ! nsub


      if (myid==0) then
         inquire(file=trim(rwdir(myid+1))//"vev.file", exist=fileexists)
         if (.not. fileexists) then
            open(unit=118,file=trim(rwdir(myid+1))//"vev.file",action="write", form="formatted",status="new")
         else
            open(unit=118,file=trim(rwdir(myid+1))//"vev.file",action="write", form="formatted",status="old",position="append")
         endif ! file check

        do iop = 1,nop
           do imom = 1,nmom
              do it = 1,nt
                 write(unit=118, fmt="('isub, iop, imom, it, Real value, Imag value :: '4i3,2es17.10)")&
                                           4, iop, imom, it, Jvev(1,it,4,imom,iop)/real(nx*ny*nz,KR), Jvev(2,it,4,imom,iop)/real(nx*ny*nz,KR)
              enddo !iop
           enddo ! isub
        enddo ! it       
        do iop = 1,nop
           do imom = 1,nmom
              do it = 1,nt
                 write(unit=118, fmt="('isub, iop, imom, it, Real value, Imag value :: ',4i3,2es17.10)")&
                                           6, iop, imom, it, Jvev(1,it,6,imom,iop)/real(nx*ny*nz,KR), Jvev(2,it,6,imom,iop)/real(nx*ny*nz,KR)
              enddo !iop
           enddo ! isub
        enddo ! it       

        close (unit=118,status="keep")
      endif ! myid


! Copy Jvev into  Jtotal so that it can be added to the average oper. values.
! NEED TO NORMALIZE THE VEV BY NUMBER OF NOISES!!!
      Jtotal = numnoises*Jvev
      !if(myid==0) print *,"I am at A"
! Unlike Jtotal, Jtvev will be added to Jtemp which has noise information. I 
! have choosen to add in the vev for EACH noise with the array Jtvev.
      Jtvev(:,:,:,:,:) = Jvev(:,:,:,:,:)
      !if(myid==0) print *, "I am at B"
 if (.false.) then
  if (myid==0) then
    do it=1,nt
      !print "(a48,2i2,2es17.10)", "OUR: isub,it, t-Jvev=",1,it, Jvev(1,it,4,1,1),Jvev(2,it,4,1,1)
    enddo ! it
    do it=1,nt
      ! print "(a48,2i2,2es17.10)", "OUR: isub,it, x-Jvev=",4,it, Jvev(1,it,4,1,3),Jvev(2,it,4,1,3)
    enddo ! it
    do it=1,nt
      !print "(a48,2i2,2es17.10)", "OUR: isub,it, y-Jvev=",4,it, Jvev(1,it,4,1,4),Jvev(2,it,4,1,4)
    enddo ! it
    do it=1,nt
      !print "(a48,2i2,2es17.10)", "OUR: isub,it, z-Jvev=",4,it, Jvev(1,it,4,1,5),Jvev(2,it,4,1,5)
    enddo ! it
    do it=1,nt
      print "(a48,2i2,2es17.10)", "OUR: isub,it, psi-Jvev=",4,it, Jvev(1,it,4,1,2)/real(nx*ny*nz,KR),&
                                                                  Jvev(2,it,4,1,2)/real(nx*ny*nz,KR)
    enddo ! it
  endif ! myid
 endif ! .false.
     !call printlog("STOPPING!!!",myid,rwdir)
     !call MPI_BARRIER(MPI_COMM_WORLD,ierr)
     !stop      

! Initialize.
! I do this inside twAverage.......
!   call momfacs(io,momfac,myid,iblv)

! Loop over noises. (numnoises=0 is intended for debugging purposes only.)

    if (numnoises==0) then
     inverter=0
     call fermprop(rwdir,bedebug,bodebug,iflag,kapvec,cSW,inverter,coact,bc,resmax, &
                   itermin,omegaM3R,u,vectype,iplug,xe,xo,myid,nn,ldiv,nms, &
                   lvbc,lvcl,ib,lbd,iblv,vecbl,vecblinv,0,hoption,MRT,MRT2)

     call twAverage(Jave,be,bo,xe,xo,delta,&
                    io,kapvec(1),u,dobndry,numprocs,MRT,&
                    myid,nsub,nmom,nop,1,shiftnum,ntmqcd,&
                    rwdir,ir)


     Jtotal = Jtotal + Jave

! Save the results.
      if (myid==0) then
         opfile = trim(rwdir(myid+1))//"Discon.Total.LOG"
         inquire(file=opfile, exist=fileexists) 
         if (.not. fileexists) then
            open(unit=18,file=trim(rwdir(myid+1))//"Discon.Total.LOG",action="write", &
                 form="formatted",status="new")
         else
            open(unit=18,file=trim(rwdir(myid+1))//"Discon.Total.LOG",action="write", &
                 form="formatted",status="old",position="append")
         endif ! file check
         do iop = 1,nop
          do imom = 1,nmom
           do isub = 1,nsub
            do iri = 1,2
               write(unit=18,fmt="(a6,4i3,500es17.10)") &
                     "discon: iri,isub,imom,iop,shiftnum,Jtotal=" &
                      ,iri,isub,imom,iop,shiftnum,Jtotal(iri,:,isub,imom,iop)
            enddo ! iri
           enddo ! isub
          enddo ! imom
         enddo ! iop
         close(unit=18,status="keep")
       endif ! myid

    else

        if(sweepnum==1.and.ntmqcd.ge.0) then
          Jconu= 0.0_KR
        elseif (sweepnum==1.and.ntmqcd.lt.0) then
          Jcond=0.0_KR
        endif ! sweepnum

! NOISE LOOP!!!!!
        do inoise = 1,numnoises
if (myid==0) then
           open(unit=8,file=trim(rwdir(myid+1))//"CFGSPROPS.LOG",action="write", &
                form="formatted",status="old",position="append")
           write(unit=8,fmt="(a50,i5)") "enter inoise ",inoise
           close(unit=8,status="keep")
endif


           !print *, "myid=",myid," Inside noise loop doing noise number ",inoise
           if (inverter==6) then
              call setrhsnoise(inoise)
           endif ! inverter==6

           !print *,"myid=",myid," finished with setrhsnoise"


! Note to self: put pointsource here.
!         src(0,1)=1
!         src(0,2)=1
!         src(0,3)=1
!         src(0,4)=4
          !print *,"myid=",myid," before z2sourc"
if (myid==0) then
           open(unit=8,file=trim(rwdir(myid+1))//"CFGSPROPS.LOG",action="write", &
                form="formatted",status="old",position="append")
           write(unit=8,fmt="(a50)") "enter z2"
           close(unit=8,status="keep")
endif
            call z2source(be,bo,nvhalf)
if (myid==0) then
           open(unit=8,file=trim(rwdir(myid+1))//"CFGSPROPS.LOG",action="write", &
                form="formatted",status="old",position="append")
           write(unit=8,fmt="(a50)") "exiting z2"
           close(unit=8,status="keep")
endif
          !print *,"myid=",myid, " finished with z2source"
          !call MPI_BARRIER(MPI_COMM_WORLD, ierr)
                
          if(myid.eq.0.and..false.) then 
           open(unit=8,file=trim(rwdir(myid+1))//"CFGSPROPS.LOG",action="write", &
            form="formatted",status="old",position="append")
           write(unit=8,fmt=*)"ir=",ir,"inoise=",inoise,"kappa=",kappa,"mu=",mu
           write(unit=8,fmt=*)"icri=1,id=1,ieo=1,ibl=1"
           write(unit=8,fmt=*)"be="
           do i=1,10
            write(unit=8,fmt=*)i,be(1,i,1,1,1)
           enddo !i

           write(unit=8,fmt=*)"bo="
           do i=1,10
            write(unit=8,fmt=*)i,bo(1,i,1,1,1)
           enddo !i
          
          endif


!         icsrc=1
!         idsrc=1
!         call pointsource(be,bo,src(0,:),icsrc,idsrc,myid,iblv,vecblinv)

! Rotate the input vector so that the propagators come from Fermprop
! in the physical basis.

          !call gamma5vector(bo,ntmqcd,myid)
          !call gamma5vector(be,ntmqcd,myid)

! Note: irloop again on r Wilson value.

if (myid==0) then
           open(unit=8,file=trim(rwdir(myid+1))//"CFGSPROPS.LOG",action="write", &
                form="formatted",status="old",position="append")
           write(unit=8,fmt="(a50)") "enter ferm prop"
           close(unit=8,status="keep")
endif
          call fermpropR(rwdir,be,bo,iflag,kapvec,cSW,inverter,coact,bc,resmax, &
                        itermin,omegaM3R,u,vectype,iplug,xe,xo,myid,nn,ldiv,nms, &
                        lvbc,lvcl,ib,lbd,iblv,vecbl,vecblinv,0,hoption,MRT,MRT2,ir)
if (myid==0) then
           open(unit=8,file=trim(rwdir(myid+1))//"CFGSPROPS.LOG",action="write", &
                form="formatted",status="old",position="append")
           write(unit=8,fmt="(a50)") "exiting fermprop"
           close(unit=8,status="keep")
endif

! ABDOU ~ subroutine changevector in twAverage probably needs to be changed.
         
          
          call cpu_time(time3)

           if (myid==0) then
            open(unit=8,file=trim(rwdir(myid+1))//"CFGSPROPS.LOG",action="write", &
                 form="formatted",status="old",position="append")
            write(unit=8,fmt="(a50,es17.10,a10)") "Started twAverage at:",time3,"secs"
            close(unit=8,status="keep")
           endif

          call twAverage(Jave,be,bo,xe,xo,delta,&
                         io,kapvec(1),u,dobndry,numprocs,MRT,&
                         myid,nsub,nmom,nop,inoise,shiftnum,ntmqcd,&
                         rwdir,ir)

          !call twAverage_axial(Jave,be,bo,xe,xo,delta,&
          !               io,kapvec(1),u,dobndry,numprocs,MRT,&
          !               myid,nsub,nmom,nop,inoise,shiftnum,ntmqcd,&
          !               rwdir,ir)
          call cpu_time(time4)
          if (myid==0) then
           open(unit=8,file=trim(rwdir(myid+1))//"CFGSPROPS.LOG",action="write", &
                form="formatted",status="old",position="append")
           write(unit=8,fmt="(a50,es17.10,a10)") "Finishd twAvergae in:",time4-time3,"secs"
           close(unit=8,status="keep")
          endif

! ABDOU ~ need to have mass indicies on all "J-arrays" that have multi-mass
!         solutions. So, Jtotal, Jave, Jtvev, etc...
!         These are the arrays that ultimately output the data into our Signal.dat files



! To match Randy's normilization multiply the currents by kapvec(1)

          do iop=1,nop
             if (iop/=2 .and. iop/=6) then
                 Jave(:,:,:,:,iop) = kapvec(1)*Jave(:,:,:,:,iop)
             endif ! iop
          enddo ! iop
 
! For the first noise Jtotal=Jvev. Each following noise, Jtotal
! contains the Jave+Jvev sum from the previous noise. 

          Jtotal = Jtotal + Jave

! The way it was:
! Jtotal = Jtotal + Jave + Jvev, with Jtotal originally zeroed out.

! Print statements for debugging.
    if (.false.) then
          if (myid==0) then
           do it=1,nt
              print  "(a48,2i2,2es17.10)", "OUR: isub,it, t-Jave=", 1,it, Jave(1,it,1,1,1),&
                                                                          Jave(2,it,1,1,1) 
           enddo ! it
          endif
    
          if (myid==0) then
           do it=1,nt
              print  "(a48,2i2,2es17.10)", "OUR: isub,it, x-Jave=", 1,it, Jave(1,it,1,1,3),&
                                                                          Jave(2,it,1,1,3) 
           enddo ! it
          endif

          if (myid==0) then
           do it =1,nt
             print  "(a48,2i2,2es17.10)", "OUR: isub,it, y-Jave=", 1,it, Jave(1,it,1,1,4),&
                                                                         Jave(2,it,1,1,4) 
           enddo ! it
          endif

          if (myid==0) then
           do it =1,nt
             print  "(a48,2i2,2es17.10)", "OUR: isub,it, z-Jave=", 1,it, Jave(1,it,1,1,5),&
                                                                         Jave(2,it,1,1,5) 
           enddo ! it
          endif

          if (myid==0) then
           do it =1,nt
             print  "(a48,2i2,2es17.10)", "OUR: isub,it, psi-Jave=", 1,it, Jave(1,it,1,1,2),&
                                                                           Jave(2,it,1,1,2) 
           enddo ! it
          endif

          if (myid==0) then
           do it =1,nt
             print  "(a48,2i2,2es17.10)", "OUR: isub,it, psi-Jave=", 1,it, Jave(1,it,1,1,6),&
                                                                           Jave(2,it,1,1,6) 
           enddo ! it
          endif
     endif ! .false.

       if (myid==0 .and. .false.) then
          Jp1 = 0.0_KR
          Jp2 = 0.0_KR
          do ittt=1,nt
             Jp1(:,:,:,:) = Jp1(:,:,:,:) +Jave(:,ittt,:,:,:)/nt
             Jp2(:,:,:,:) = Jp2(:,:,:,:) +Jvev(:,ittt,:,:,:)/nt
          enddo ! ittt
          if (nsub==0) then
             call printarray1("rho=",myid,rwdir,nsub+1,Jp1(1,:,1,1))
!            call printarray1("Jtotal1=",myid,rwdir,nsub+1,Jp2(1,:,1,1))
             call printarray1("scalar=",myid,rwdir,nsub+1,Jp1(1,:,1,2)*nt)
!            call printarray1("Jtotal2=",myid,rwdir,nsub+1,Jp2(1,:,1,2)*nt)
!            call printarray1("Jtotal1 level 4",myid,rwdir,nt,Jtotal(1,:,4,1,1))
!            call printarray1("Jtotal1 level 6",myid,rwdir,nt,Jtotal(1,:,6,1,1))
          else
             call printarray1("Jnew1=",myid,rwdir,nsub,Jp1(1,:,1,1))
             call printarray1("Jtotal1=",myid,rwdir,nsub,Jp2(1,:,1,1))
             call printarray1("Jnew2=",myid,rwdir,nsub,Jp1(1,:,1,2)*nt)
             call printarray1("Jtotal2=",myid,rwdir,nsub,Jp2(1,:,1,2)*nt)
!            call printarray1("Jtotal1 level 4",myid,rwdir,nt,Jtotal(1,:,4,1,1))
!            call printarray1("Jtotal1 level 6",myid,rwdir,nt,Jtotal(1,:,6,1,1))
          endif ! nsub==0
        endif ! myid


! Add in the vev for each noise when doing Jtemp. (Differnet than Jtotal which has
! no noise information and N*Jvev is added on the first noise.)
          do iop=1,nop
             Jtemp(:,:,:,:,iop,inoise) = Jave(:,:,:,:,iop) + Jtvev(:,:,:,:,iop)
          enddo ! iop

       if (.false.) then
         if (myid==0) then
           fac = 1.0_KR/real(nx*ny*nz,KR)
           do it=1,nt
             write(unit=6, fmt="(a20,3i3,2es24.10)")  "zero Order", 1,it,1,fac*Jtemp(1,it,1,1,5,inoise),fac*Jtemp(2,it,1,1,5,inoise)
           enddo ! it
           do it=1,nt
             write(unit=6, fmt="(a20,3i3,2es24.10)")  "first Order", 1,it,1,fac*Jtemp(1,it,4,1,5,inoise),fac*Jtemp(2,it,4,1,5,inoise)
           enddo ! it
           do it=1,nt
             write(unit=6, fmt="(a20,3i3,2es24.10)")  "Highest Order", 1,it,1,fac*Jtemp(1,it,6,1,5,inoise),fac*Jtemp(2,it,6,1,5,inoise)
           enddo ! it
         endif ! myid
       endif ! .false.
if (myid==0) then
           open(unit=8,file=trim(rwdir(myid+1))//"CFGSPROPS.LOG",action="write", &
                form="formatted",status="old",position="append")
           write(unit=8,fmt="(a50,i10)") "exiting inoise ",inoise
           close(unit=8,status="keep")
endif


      enddo ! inoise


! Normalize by the noises.
        Jtotal= Jtotal/real(numnoises,KR)

! Normalize by the lattice space points.
        Jtotal = Jtotal/real(nx*ny*nz,KR)

        if(myid==0) then
          call printarray1("scalar=",myid,rwdir,nsub,Jtotal(1,1,:,1,2))
        endif ! myid
 
        if (ntmqcd.ge.0) then
           Jconu  = Jconu + Jtotal
        else if (ntmqcd.lt.0) then
           Jcond  = Jcond + Jtotal
        endif ! ntmqcd

        if(myid==0) then
          if(ntmqcd.gt.0) then
            call printarray1("Jconu1=",myid,rwdir,nsub,Jconu(1,1,:,1,2))
           else if (ntmqcd.lt.0) then
            call printarray1("Jcond1=",myid,rwdir,nsub,Jcond(1,1,:,1,2))
          endif ! ntmqcd
        endif ! myid
       
! Save the results.
        call printlog("about to print",myid,rwdir)
        if (myid==0) then
           opfile = trim(rwdir(myid+1))//"Discon.Total.LOG"
           inquire(file=opfile, exist=fileexists)
           if (.not. fileexists) then
              open(unit=18,file=trim(rwdir(myid+1))//"Discon.Total.LOG",action="write", &
                   form="formatted",status="new")
           else
              open(unit=18,file=trim(rwdir(myid+1))//"Discon.Total.LOG",action="write", &
                   form="formatted",status="old",position="append")
           endif ! file check
           do iop = 1,nop
            do imom = 1,nmom
             do isub = 1,nsub
              do it=1,nt
                 if (ntmqcd >= 0) then
                    if (Jconu(1,it,isub,imom,iop) /= 0.0_KR .or. Jconu(2,it,isub,imom,iop) /= 0.0_KR) then
                        write(unit=18,fmt="(a30,4i3,2es17.10)") &
                              "it,isub,imom,iop, tmUloop:",it,isub,imom,iop,Jconu(1,it,isub,imom,iop),Jconu(2,it,isub,imom,iop)
                    endif ! check for zero
                 else if (ntmqcd < 0) then
                    if (Jcond(1,it,isub,imom,iop) /= 0.0_KR .or. Jcond(2,it,isub,imom,iop) /= 0.0_KR) then
                        write(unit=18,fmt="(a30,4i3,2es17.10)") &
                              "it,isub,imom,iop, tmDloop:",it,isub,imom,iop,Jcond(1,it,isub,imom,iop),Jcond(2,it,isub,imom,iop)
                    endif ! check for zero again
                 endif ! ntmqcd
              enddo ! it
             enddo ! isub
            enddo ! imom
           enddo ! iop
           close(unit=18,status="keep")
         endif ! myid

! Numnoise and Nsweep both equal to one is not allowed.

! ABDOU ~ these output files need a mass index (if you are multi-massing)

       trailer= ".x"
       looptrailer=".xx"
       irtrailer=".xx"
       write(unit=trailer(2:2),fmt="(i1.1)") shiftnum
       write(unit=looptrailer(2:3),fmt="(i2.2)") shiftnum
       write(unit=irtrailer(2:3),fmt="(i2.2)")irloop

       if (myid==0) then 
          fac = 1.0_KR/real(nx*ny*nz,KR)
! Ave of noises
          Jtemp= fac* Jtemp
             ltrailer = "xxxx.Signal.dat"
             write(unit=ltrailer(1:4), fmt="(1i4.4)") icfgsave
             if (ntmqcd > 0) then
                opfile = trim(rwdir(myid+1))//trim("datafiles/")//trim(ltrailer)//&
                         trim(".loopmass")//trim(looptrailer)//trim(".irloop")//&
                         trim(irtrailer)//trim(".tmU")
             else if (ntmqcd < 0) then
                opfile = trim(rwdir(myid+1))//trim("datafiles/")//trim(ltrailer)//&
                         trim(".loopmass.")//trim(looptrailer)//trim(".irloop")//&
                         trim(irtrailer)//trim(".tmD")
             else if (ntmqcd==0) then
                if (shiftnum==1) opfile = trim(rwdir(myid+1))//trim("datafiles/")//trim(ltrailer)//trim(".152")
                if (shiftnum==2) opfile = trim(rwdir(myid+1))//trim("datafiles/")//trim(ltrailer)//trim(".154")
             endif ! ntmqcd
             tempfile = opfile
             looptrailer = ".xx"


             do iop=1,nop
                write(unit=looptrailer(2:3),fmt="(i2.2)") iop
                opfile = trim(tempfile)//trim(looptrailer)
                !print *, opfile
                inquire(file=opfile, exist=fileexists)
                if (.not.  fileexists ) then
                   open(unit=8,file=opfile,action="write",form="formatted",status="new")
                else
                   open(unit=8,file=opfile,action="write",form="formatted",status="old",position="append")
                endif ! open
                !write(unit=8, fmt="(3i3,2es24.10)")  1, 2, 3, 4.0, 5.0
                !close(unit=8,status="keep")
                !call MPI_BARRIER(MPI_COMM_WORLD, ierr)
                !stop
                if (nsub==0) then 
                   do inoise=1,numnoises
                      do imom=1,nmom
                         do it=1,nt
                            write(unit=8, fmt="(3i3,2es24.10)")  1,it,imom,Jtemp(1,it,1,imom,iop,inoise),Jtemp(2,it,1,imom,iop,inoise)
                         enddo ! it
                      enddo ! imom
                   enddo ! inoise
                endif ! nsub==0

                if (nsub==4) then 
                   do inoise=1,numnoises
                      do imom=1,nmom
                         do it=1,nt
                            write(unit=8, fmt="(3i3,2es24.10)")  1,it,imom,Jtemp(1,it,1,imom,iop,inoise),Jtemp(2,it,1,imom,iop,inoise)
                         enddo ! it
                      enddo ! imom
                   enddo ! inoise
                   do inoise=1,numnoises
                      do imom=1,nmom
                         do it=1,nt
                            write(unit=8, fmt="(3i3,2es24.10)")  4,it,imom,Jtemp(1,it,4,imom,iop,inoise),Jtemp(2,it,4,imom,iop,inoise)
                         enddo ! it
                      enddo ! imom
                   enddo ! inoise
                endif ! nsub==0

                if (nsub==6) then
! Look here, Carl.
                  !call printlog("Just before operator print",myid,rwdir)
                  !print *, numnoises, nmom, nt
                   do inoise=1,numnoises
                      do imom=1,nmom
                         do it=1,nt
                            !print *, 1,it,imom,Jtemp(1,it,1,imom,iop,inoise),Jtemp(2,it,1,imom,iop,inoise)
                            !write(unit=8, fmt="(3i3,2es24.10)")  1, 2, 3, 4.0, 5.0
                            write(unit=8, fmt="(3i3,2es24.10)")  1,it,imom,Jtemp(1,it,1,imom,iop,inoise),Jtemp(2,it,1,imom,iop,inoise)
                         enddo ! it
                      enddo ! imom
                   enddo ! inoise
                   do inoise=1,numnoises
                      do imom=1,nmom
                         do it=1,nt
                            write(unit=8, fmt="(3i3,2es24.10)")  4,it,imom,Jtemp(1,it,4,imom,iop,inoise),Jtemp(2,it,4,imom,iop,inoise)
                         enddo ! it
                      enddo ! imom
                   enddo ! inoise
                   do inoise=1,numnoises
                      do imom=1,nmom
                         do it=1,nt
                            write(unit=8, fmt="(3i3,2es24.10)")  6,it,imom,Jtemp(1,it,6,imom,iop,inoise),Jtemp(2,it,6,imom,iop,inoise)
                         enddo ! it
                      enddo ! imom
                   enddo ! inoise
                endif ! nsub==0

                close(unit=8,status="keep")
             enddo ! iop
         !endif ! sweepnum
       endif ! myid
     !call MPI_BARRIER(MPI_COMM_WORLD, ierr)                
     !stop

     endif ! (numnoises = 0)
     call printlog("Done with Signal.dat",myid,rwdir)

  deallocate(Jtemp)
  deallocate(xe)
  deallocate(xo)
  deallocate(be)
  deallocate(bo)
  deallocate(ze)
  deallocate(zo)


  call cpu_time(time2)

  if (myid==0) then
   open(unit=8,file=trim(rwdir(myid+1))//"CFGSPROPS.LOG",action="write", &
       form="formatted",status="old",position="append")
   write(unit=8,fmt="(a50,es17.10,a10)") "Time taken by twistdiscon:",time2-time1,"secs"
   close(unit=8,status="keep")
  endif

!    if (myid==0) print *, "Time to do noises in twistdiscon=", time2 - time1

 end subroutine twistdiscon
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -






 subroutine z4source(be,bo,nvhalf,myid)
! Define a source vector containing real Z4 noise.

    real(kind=KR),    intent(out), dimension(:,:,:,:,:) :: be, bo
    integer(kind=KI), intent(in)                        :: nvhalf,myid

    real(kind=KR), dimension(192*nvhalf) :: rn
    integer(kind=KI)                   :: icount, ic, icri, isite, id, ieo, ibl,ierr


    ! z8 nosie
    if (.false.) then
    be = 0.0_KR
    icount = 192*nvhalf
    call rnds(icount,rn)
    icount = 0
    do ic = 1,3
     icri = 2*ic - 1
     do isite = 1,nvhalf
      do id = 1,4
       do ieo = 1,2
        do ibl = 1,8
         icount = icount + 1

         if (rn(icount) > .875_KR) then
           be(icri  ,isite,id,ieo,ibl) = 1.0_KR
           be(icri+1,isite,id,ieo,ibl) = 0.0_KR
         elseif(rn(icount) > .75_KR) then
           be(icri  ,isite,id,ieo,ibl) = 0.707106781187_KR
           be(icri+1,isite,id,ieo,ibl) = 0.707106781187_KR
         elseif(rn(icount) > .625_KR) then
           be(icri  ,isite,id,ieo,ibl) = 0.0_KR
           be(icri+1,isite,id,ieo,ibl) = 1.0_KR
         elseif(rn(icount) > .5_KR) then
           be(icri  ,isite,id,ieo,ibl) = -0.707106781187_KR
           be(icri+1,isite,id,ieo,ibl) =  0.707106781187_KR
         elseif(rn(icount) > .375_KR) then
           be(icri  ,isite,id,ieo,ibl) = -1.0_KR
           be(icri+1,isite,id,ieo,ibl) =  0.0_KR
         elseif(rn(icount) > .25_KR) then
           be(icri  ,isite,id,ieo,ibl) = -0.707106781187_KR
           be(icri+1,isite,id,ieo,ibl) = -0.707106781187_KR
         elseif(rn(icount) > .125_KR) then
           be(icri  ,isite,id,ieo,ibl) =  0.0_KR
           be(icri+1,isite,id,ieo,ibl) = -1.0_KR
         else
           be(icri  ,isite,id,ieo,ibl) =  0.707106781187_KR
           be(icri+1,isite,id,ieo,ibl) = -0.707106781187_KR
         endif

        enddo ! ibl
       enddo ! ieo
      enddo ! id
     enddo ! isite
    enddo ! ic

    bo = 0.0_KR
    icount = 192*nvhalf
    call rnds(icount,rn)
    icount = 0
    do ic = 1,3
     icri = 2*ic - 1
     do isite = 1,nvhalf
      do id = 1,4
       do ieo = 1,2
        do ibl = 1,8
         icount = icount + 1

         if (rn(icount) > .875_KR) then
           bo(icri  ,isite,id,ieo,ibl) = 1.0_KR
           bo(icri+1,isite,id,ieo,ibl) = 0.0_KR
         elseif(rn(icount) > .75_KR) then
           bo(icri  ,isite,id,ieo,ibl) = 0.707106781187_KR
           bo(icri+1,isite,id,ieo,ibl) = 0.707106781187_KR
         elseif(rn(icount) > .625_KR) then
           bo(icri  ,isite,id,ieo,ibl) = 0.0_KR
           bo(icri+1,isite,id,ieo,ibl) = 1.0_KR
         elseif(rn(icount) > .5_KR) then
           bo(icri  ,isite,id,ieo,ibl) = -0.707106781187_KR
           bo(icri+1,isite,id,ieo,ibl) =  0.707106781187_KR
         elseif(rn(icount) > .375_KR) then
           bo(icri  ,isite,id,ieo,ibl) = -1.0_KR
           bo(icri+1,isite,id,ieo,ibl) =  0.0_KR
         elseif(rn(icount) > .25_KR) then
           bo(icri  ,isite,id,ieo,ibl) = -0.707106781187_KR
           bo(icri+1,isite,id,ieo,ibl) = -0.707106781187_KR
         elseif(rn(icount) > .125_KR) then
           bo(icri  ,isite,id,ieo,ibl) =  0.0_KR
           bo(icri+1,isite,id,ieo,ibl) = -1.0_KR
         else
           bo(icri  ,isite,id,ieo,ibl) =  0.707106781187_KR
           bo(icri+1,isite,id,ieo,ibl) = -0.707106781187_KR
         endif

        enddo ! ibl
       enddo ! ieo
      enddo ! id
     enddo ! isite
    enddo ! ic
    endif



    ! z3 nosie
    if (.false.) then
    be = 0.0_KR
    icount = 192*nvhalf
    call rnds(icount,rn)
    icount = 0
    do ic = 1,3
     icri = 2*ic - 1
     do isite = 1,nvhalf
      do id = 1,4
       do ieo = 1,2
        do ibl = 1,8
         icount = icount + 1
         if (rn(icount) > 0.66666666_KR) then
          be(icri,isite,id,ieo,ibl) = 1.0_KR
         elseif (rn(icount) > 0.33333333_KR) then
          be(icri  ,isite,id,ieo,ibl) = -0.5_KR
          be(icri+1,isite,id,ieo,ibl) =  0.866025403784
         else
          be(icri  ,isite,id,ieo,ibl) = -0.5_KR
          be(icri+1,isite,id,ieo,ibl) = -0.866025403784
         endif
        enddo ! ibl
       enddo ! ieo
      enddo ! id
     enddo ! isite
    enddo ! ic

    bo = 0.0_KR
    icount = 192*nvhalf
    call rnds(icount,rn)
    icount = 0
    do ic = 1,3
     icri = 2*ic - 1
     do isite = 1,nvhalf
      do id = 1,4
       do ieo = 1,2
        do ibl = 1,8
         icount = icount + 1
         if (rn(icount) > 0.66666666_KR) then
          bo(icri,isite,id,ieo,ibl) = 1.0_KR
         elseif (rn(icount) > 0.33333333_KR) then
          bo(icri  ,isite,id,ieo,ibl) = -0.5_KR
          bo(icri+1,isite,id,ieo,ibl) =  0.866025403784_KR
         else
          bo(icri  ,isite,id,ieo,ibl) = -0.5_KR
          bo(icri+1,isite,id,ieo,ibl) = -0.866025403784_KR
         endif
        enddo ! ibl
       enddo ! ieo
      enddo ! id
     enddo ! isite
    enddo ! ic
    endif

    ! z4 noise
    if (.true.) then
    be = 0.0_KR
    icount = 192*nvhalf
    call rnds(icount,rn)
    icount = 0
    do ic = 1,3
     icri = 2*ic - 1
     do isite = 1,nvhalf
      do id = 1,4
       do ieo = 1,2
        do ibl = 1,8
         icount = icount + 1
         if (rn(icount) > 0.75_KR) then
          be(icri,isite,id,ieo,ibl) = 1.0_KR
         elseif (0.5_KR < rn(icount) .and. rn(icount) < 0.75_KR) then
          be(icri,isite,id,ieo,ibl) = -1.0_KR
         elseif (.25_KR < rn(icount) .and. rn(icount) < 0.5_KR) then
          be(icri+1,isite,id,ieo,ibl) = 1.0_KR
         else
          be(icri+1,isite,id,ieo,ibl) = -1.0_KR
         endif
        enddo ! ibl
       enddo ! ieo
      enddo ! id
     enddo ! isite
    enddo ! ic

    bo = 0.0_KR
    icount = 192*nvhalf
    call rnds(icount,rn)
    icount = 0
    do ic = 1,3
     icri = 2*ic - 1
     do isite = 1,nvhalf
      do id = 1,4
       do ieo = 1,2
        do ibl = 1,8
         icount = icount + 1
         if (rn(icount) > 0.75_KR) then
          bo(icri,isite,id,ieo,ibl) = 1.0_KR
         elseif (0.5_KR < rn(icount) .and. rn(icount) < 0.75_KR) then
          bo(icri,isite,id,ieo,ibl) = -1.0_KR
         elseif (.25_KR < rn(icount) .and. rn(icount) < 0.5_KR) then
          bo(icri+1,isite,id,ieo,ibl) = 1.0_KR
         else
          bo(icri+1,isite,id,ieo,ibl) = -1.0_KR
         endif
        enddo ! ibl
       enddo ! ieo
      enddo ! id
     enddo ! isite
    enddo ! ic
    endif
 end subroutine z4source











 subroutine z2source(be,bo,nvhalf)
! Define a source vector containing real Z2 noise.

    real(kind=KR),    intent(out), dimension(:,:,:,:,:) :: be, bo
    integer(kind=KI), intent(in)                        :: nvhalf

    real(kind=KR), dimension(192*nvhalf) :: rn
    integer(kind=KI)                   :: icount, ic, icri, isite, id, ieo, ibl,ierr

    be = 0.0_KR
    icount = 192*nvhalf
    !print *,"Inside z2 source, before rnds"
    call rnds(icount,rn)
    !print *,"Inside z2 source , after rnds"


    !call MPI_BARRIER(MPI_COMM_WORLD, ierr)  
    !print *,"ierr=",ierr

    icount = 0
    do ic = 1,3
     icri = 2*ic - 1
     do isite = 1,nvhalf
      do id = 1,4
       do ieo = 1,2
        do ibl = 1,8
         icount = icount + 1
         if (rn(icount)>0.5_KR) then
          be(icri,isite,id,ieo,ibl) = 1.0_KR
         else
          be(icri,isite,id,ieo,ibl) = -1.0_KR
         endif
        enddo ! ibl
       enddo ! ieo
      enddo ! id
     enddo ! isite
    enddo ! ic

    bo = 0.0_KR
    icount = 192*nvhalf
    call rnds(icount,rn)
    icount = 0
    do ic = 1,3
     icri = 2*ic - 1
     do isite = 1,nvhalf
      do id = 1,4
       do ieo = 1,2
        do ibl = 1,8
         icount = icount + 1
         if (rn(icount)>0.5_KR) then
          bo(icri,isite,id,ieo,ibl) = 1.0_KR
         else
          bo(icri,isite,id,ieo,ibl) = -1.0_KR
         endif
        enddo ! ibl
       enddo ! ieo
      enddo ! id
     enddo ! isite
    enddo ! ic


    !print *,"Before exiting z2 source"

    !call MPI_BARRIER(MPI_COMM_WORLD, ierr)  
    !print *,"ierr=",ierr
              
 end subroutine z2source

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 subroutine vev(Jvev,u,rwdir,io,kappa,dobndry,myid,nn,ldiv,nms,lvbc,ib,lbd, &
                iblv,MRT)
 ! Calculates the average plaquette value of certain operators for a given
 ! configuration.

 real(kind=KR),    intent(out),   dimension(:,:,:,:)   :: Jvev
 integer(kind=KI), intent(in)                          :: dobndry, myid, MRT
 real(kind=KR),    intent(inout), dimension(:,:,:,:,:) :: u
 character(len=*), intent(in),    dimension(:)         :: rwdir
 integer(kind=KI), intent(in),    dimension(:)         :: io
 real(kind=KR),    intent(in)                          :: kappa
 integer(kind=KI), intent(in),    dimension(:)         :: nms
 integer(kind=KI), intent(in),    dimension(:,:)       :: nn, iblv
 logical,          intent(in),    dimension(:)         :: ldiv
 integer(kind=KI), intent(in),    dimension(:,:,:)     :: lvbc
 integer(kind=KI), intent(in),    dimension(:,:,:,:)   :: ib
 logical,          intent(in),    dimension(:,:)       :: lbd

 integer(kind=KI) :: i, mu, nu, isite, ieo, ibl, imom, ieo1, &
                     ieo2, ixbit, iybit, izbit, itbit, ixbit2, iybit2, &
                     izbit2, itbit2, iblbit, iy, iz, it, ixbit3, icount, &
                     nmesg, ierr, ipp, itmin, jeo
 integer(kind=KI), dimension(4)           :: np, ip, muvec
 integer(kind=KI), dimension(2)           :: ix
 integer(kind=KI), parameter              :: nmom=5, nop=5
 integer(kind=KI), parameter              :: itstep=2*nvhalf*npt/nt
 real(kind=KR)                            :: fac1, fac2, fac3, fac4, psi, &
                                             plaq, psitemp, plaqtemp, factor
 real(kind=KR),  dimension(18,itstep,2,8) :: ubndry
 real(kind=KR),  dimension(nvhalf)        :: cosx, cosy, cosz, &
                                             cos2x, cos2y, cos2z, &
                                             sinx, siny, sinz, &
                                             sin2x, sin2y, sin2z
 real(kind=KR),  dimension(3,nvhalf)      :: pcosx, pcosy, pcosz, &
                                             pcos2x, pcos2y, pcos2z, &
                                             psinx, psiny, psinz, &
                                             psin2x, psin2y, psin2z
 real(kind=KR),  dimension(18,ntotal)     :: staple
 real(kind=KR),  dimension(18,nvhalf)     :: staptemp, elemloop
 real(kind=KR),  dimension(nt)            :: tsum, ssum, j1t, j2t, j3t
 real(kind=KR),  dimension(4,nt)          :: mtsum, mssum, j4sum, j12s, &
                                             j13s, j21s, j23s, j31s, j32s, &
                                             j12t, j13t, j21t, j23t, j31t, &
                                             j32t
 real(kind=KR),  dimension(2,nt,nmom,nop) :: Jvevtemp
 real(kind=KR2), parameter                :: twoPi=6.283185307179586

 ! This local gaugelink is used for debugging
 real(kind=KR),  dimension(18,ntotal,4,2,16) :: uout

 ! Some initializations.
 plaq = 0.0_KR
 tsum = 0.0_KR
 ssum = 0.0_KR
 j1t = 0.0_KR
 j2t = 0.0_KR
 j3t = 0.0_KR
 mtsum = 0.0_KR
 mssum = 0.0_KR
 j4sum = 0.0_KR
 j12s = 0.0_KR
 j13s = 0.0_KR
 j21s = 0.0_KR
 j23s = 0.0_KR
 j31s = 0.0_KR
 j32s = 0.0_KR
 j12t = 0.0_KR
 j13t = 0.0_KR
 j21t = 0.0_KR
 j23t = 0.0_KR
 j31t = 0.0_KR
 j32t = 0.0_KR

 ! This routine is used for debugging.
 if (.false.) then ! debug
   call printlog("Debugging vev",myid,rwdir)
   call fakegauge(uout,myid,rwdir,MRT)
   u = uout
 endif ! debug

 ! Identify the location of my process.
 np(1) = npx
 np(2) = npy
 np(3) = npz
 np(4) = npt
 call atoc(myid,np,ip)
 itmin = 1 + ip(4)*nt/npt

 ! Set temporal links at the maximal timestep (on the global lattice) to zero.
 if (.false. .and. dobndry==1) then
   itbit = nvhalf-itstep
   do ibl = 9,16
     icount = ibl - 8
     do isite = itbit+1,nvhalf
       ubndry(:,isite-itbit,:,icount) = u(:,isite,4,:,ibl)
       u(:,isite,4,:,ibl) = 0.0_KR
     enddo ! isite
   enddo ! ibl
 endif

 do mu = 1,4
   ! Put the four lattice directions into an array.
   ! Note that muvec() simply runs over all 4 directions, beginning with mu.
   do i = 1,4
     muvec(i) = 1 + modulo(mu+i-2,4)
   enddo ! i

   do ieo = 1,2
     jeo = 3 - ieo
     do ibl = 1,16
       ! Now do the sum of the six plaquettes associated with a mu-direction link.
       do nu = 2,5-mu
         !-multiply first two links (from initial site go in +muvec(nu) dir, then +mu).
         call setmmf(muvec(nu),lbd(ibl,muvec(nu)),u(:,:,mu,ieo,iblv(ibl, &
                     muvec(nu))),ldiv(muvec(nu)),ib(:,muvec(nu),1,jeo), &
                     u(:,:,muvec(nu),ieo,ibl),u(:,:,mu,jeo,iblv(ibl, &
                     muvec(nu))),staptemp,lvbc(:,muvec(nu),ieo),1,nms,myid,nn, &
                     MRT)
         !-multiply in the third link (dagger of neighbour's +muvec(nu) directed link).
         call setmmf(mu,lbd(ibl,mu),u(:,:,muvec(nu),ieo,iblv(ibl,mu)), &
                     ldiv(mu),ib(:,mu,1,jeo),staptemp,u(:,:,muvec(nu),jeo, &
                     iblv(ibl,mu)),staple,lvbc(:,mu,ieo),2,nms,myid,nn,MRT)
         !-multiply in the final link (from initial site go in +mu direction).
         call mmd(nvhalf,u(:,:,mu,ieo,ibl),staple,elemloop)
         elemloop(1,:) = elemloop(1,:) + elemloop(9,:) + elemloop(17,:)
         elemloop(2,:) = elemloop(2,:) + elemloop(10,:) + elemloop(18,:)

         ! Define the phase factors based upon the location and plane of the plaquette.
         isite = 0
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
                 isite = isite + 1
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
                 cosx(isite) = cos(twoPi*real(ix(ieo)-io(1),KR2)/real(nx,KR2))
                 sinx(isite) = sin(twoPi*real(ix(ieo)-io(1),KR2)/real(nx,KR2))
                 cos2x(isite) = cos(twoPi*real(2*(ix(ieo)-io(1)),KR2)/real(nx,KR2))
                 sin2x(isite) = sin(twoPi*real(2*(ix(ieo)-io(1)),KR2)/real(nx,KR2))
                 cosy(isite) = cos(twoPi*real(iy-io(2),KR2)/real(ny,KR2))
                 siny(isite) = sin(twoPi*real(iy-io(2),KR2)/real(ny,KR2))
                 cos2y(isite) = cos(twoPi*real(2*(iy-io(2)),KR2)/real(ny,KR2))
                 sin2y(isite) = sin(twoPi*real(2*(iy-io(2)),KR2)/real(ny,KR2))
                 cosz(isite) = cos(twoPi*real(iz-io(3),KR2)/real(nz,KR2))
                 sinz(isite) = sin(twoPi*real(iz-io(3),KR2)/real(nz,KR2))
                 cos2z(isite) = cos(twoPi*real(2*(iz-io(3)),KR2)/real(nz,KR2))
                 sin2z(isite) = sin(twoPi*real(2*(iz-io(3)),KR2)/real(nz,KR2))

                 !-ipp identifies which vertex is the relevant one.
                 do ipp = 1,3
                   pcosx(ipp,isite) = cosx(isite)
                   psinx(ipp,isite) = sinx(isite)
                   pcos2x(ipp,isite) = cos2x(isite)
                   psin2x(ipp,isite) = sin2x(isite)
                   pcosy(ipp,isite) = cosy(isite)
                   psiny(ipp,isite) = siny(isite)
                   pcos2y(ipp,isite) = cos2y(isite)
                   psin2y(ipp,isite) = sin2y(isite)
                   pcosz(ipp,isite) = cosz(isite)
                   psinz(ipp,isite) = sinz(isite)
                   pcos2z(ipp,isite) = cos2z(isite)
                   psin2z(ipp,isite) = sin2z(isite)
                 enddo ! ipp
                 select case(mu)
                   case(1)
                     pcosx(1,isite) = cos(twoPi*real(1+ix(ieo)-io(1),KR2)/real(nx,KR2))
                     psinx(1,isite) = sin(twoPi*real(1+ix(ieo)-io(1),KR2)/real(nx,KR2))
                     pcos2x(1,isite) = cos(twoPi*real(2*(1+ix(ieo)-io(1)),KR2) &
                                      /real(nx,KR2))
                     psin2x(1,isite) = sin(twoPi*real(2*(1+ix(ieo)-io(1)),KR2) &
                                      /real(nx,KR2))
                     pcosx(3,isite) = pcosx(1,isite)
                     psinx(3,isite) = psinx(1,isite)
                     pcos2x(3,isite) = pcos2x(1,isite)
                     psin2x(3,isite) = psin2x(1,isite)
                   case(2)
                     pcosy(1,isite) = cos(twoPi*real(1+iy-io(2),KR2)/real(ny,KR2))
                     psiny(1,isite) = sin(twoPi*real(1+iy-io(2),KR2)/real(ny,KR2))
                     pcos2y(1,isite) = cos(twoPi*real(2*(1+iy-io(2)),KR2)/real(ny,KR2))
                     psin2y(1,isite) = sin(twoPi*real(2*(1+iy-io(2)),KR2)/real(ny,KR2))
                     pcosy(3,isite) = pcosy(1,isite)
                     psiny(3,isite) = psiny(1,isite)
                     pcos2y(3,isite) = pcos2y(1,isite)
                     psin2y(3,isite) = psin2y(1,isite)
                   case(3)
                     pcosz(1,isite) = cos(twoPi*real(1+iz-io(3),KR2)/real(nz,KR2))
                     psinz(1,isite) = sin(twoPi*real(1+iz-io(3),KR2)/real(nz,KR2))
                     pcos2z(1,isite) = cos(twoPi*real(2*(1+iz-io(3)),KR2)/real(nz,KR2))
                     psin2z(1,isite) = sin(twoPi*real(2*(1+iz-io(3)),KR2)/real(nz,KR2))
                     pcosz(3,isite) = pcosz(1,isite)
                     psinz(3,isite) = psinz(1,isite)
                     pcos2z(3,isite) = pcos2z(1,isite)
                     psin2z(3,isite) = psin2z(1,isite)
                   case(4) ! do nothing since these are never used.
                   case default
                     open(unit=8,file="DISCONLOOPS.ERROR",action="write", &
                          status="replace",form="formatted")
                     write(unit=8,fmt=*) "subroutine vev: mu =", mu
                     close(unit=8,status="keep")
                     stop
                 end select
                 select case(muvec(nu))
                   case(1)
                     pcosx(2,isite) = cos(twoPi*real(1+ix(ieo)-io(1),KR2)/real(nx,KR2))
                     psinx(2,isite) = sin(twoPi*real(1+ix(ieo)-io(1),KR2)/real(nx,KR2))
                     pcos2x(2,isite) = cos(twoPi*real(2*(1+ix(ieo)-io(1)),KR2) &
                                       /real(nx,KR2))
                     psin2x(2,isite) = sin(twoPi*real(2*(1+ix(ieo)-io(1)),KR2) &
                                       /real(nx,KR2))
                     pcosx(3,isite) = pcosx(2,isite)
                     psinx(3,isite) = psinx(2,isite)
                     pcos2x(3,isite) = pcos2x(2,isite)
                     psin2x(3,isite) = psin2x(2,isite)
                   case(2)
                     pcosy(2,isite) = cos(twoPi*real(1+iy-io(2),KR2)/real(ny,KR2))
                     psiny(2,isite) = sin(twoPi*real(1+iy-io(2),KR2)/real(ny,KR2))
                     pcos2y(2,isite) = cos(twoPi*real(2*(1+iy-io(2)),KR2)/real(ny,KR2))
                     psin2y(2,isite) = sin(twoPi*real(2*(1+iy-io(2)),KR2)/real(ny,KR2))
                     pcosy(3,isite) = pcosy(2,isite)
                     psiny(3,isite) = psiny(2,isite)
                     pcos2y(3,isite) = pcos2y(2,isite)
                     psin2y(3,isite) = psin2y(2,isite)
                   case(3)
                     pcosz(2,isite) = cos(twoPi*real(1+iz-io(3),KR2)/real(nz,KR2))
                     psinz(2,isite) = sin(twoPi*real(1+iz-io(3),KR2)/real(nz,KR2))
                     pcos2z(2,isite) = cos(twoPi*real(2*(1+iz-io(3)),KR2)/real(nz,KR2))
                     psin2z(2,isite) = sin(twoPi*real(2*(1+iz-io(3)),KR2)/real(nz,KR2))
                     pcosz(3,isite) = pcosz(2,isite)
                     psinz(3,isite) = psinz(2,isite)
                     pcos2z(3,isite) = pcos2z(2,isite)
                     psin2z(3,isite) = psin2z(2,isite)
                   case(4) ! do nothing since these are never used.
                   case default
                     open(unit=8,file="DISCONLOOPS.ERROR",action="write", &
                          status="replace",form="formatted")
                     write(unit=8,fmt=*) "subroutine vev: muvec(nu) =", muvec(nu)
                     close(unit=8,status="keep")
                     stop
                 end select

               enddo ! ixbit
             enddo ! iybit
           enddo ! izbit
         enddo ! itbit

         ! Zero momentum scalar.
         if (muvec(nu)==4.or.mu==4) then
           do i = 1,nvhalf
             it = itmin + 2*((i-1)/itstep) + (ibl-1)/8
             tsum(it) = tsum(it) + elemloop(1,i)
           enddo ! i
         else
           do i = 1,nvhalf
             it = itmin + 2*((i-1)/itstep) + (ibl-1)/8
             ssum(it) = ssum(it) + elemloop(1,i)
           enddo ! i
         endif

         ! Nonzero momentum scalar.
         ! Put in the phase factors based upon the location and plane of the plaquette.
         if (muvec(nu)==4) then
           do i = 1,nvhalf
             it = itmin + 2*((i-1)/itstep) + (ibl-1)/8
             fac1 = (cosx(i) + cosy(i) + cosz(i) &
                  + pcosx(1,i) + pcosy(1,i) + pcosz(1,i) &
                  + pcosx(2,i) + pcosy(2,i) + pcosz(2,i) &
                  + pcosx(3,i) + pcosy(3,i) + pcosz(3,i))/3.0_KR
             fac2 = (cosx(i)*cosy(i) + cosx(i)*cosz(i) + cosy(i)*cosz(i) &
                  + pcosx(1,i)*pcosy(1,i) + pcosx(1,i)*pcosz(1,i) &
                  + pcosy(1,i)*pcosz(1,i) + pcosx(2,i)*pcosy(2,i) &
                  + pcosx(2,i)*pcosz(2,i) + pcosy(2,i)*pcosz(2,i) &
                  + pcosx(3,i)*pcosy(3,i) + pcosx(3,i)*pcosz(3,i) &
                  + pcosy(3,i)*pcosz(3,i))/3.0_KR
             fac3 = cosx(i)*cosy(i)*cosz(i) &
                  + pcosx(1,i)*pcosy(1,i)*pcosz(1,i) &
                  + pcosx(2,i)*pcosy(2,i)*pcosz(2,i) &
                  + pcosx(3,i)*pcosy(3,i)*pcosz(3,i)
             fac4 = (cos2x(i) + cos2y(i) + cos2z(i) &
                  + pcos2x(1,i) + pcos2y(1,i) + pcos2z(1,i) &
                  + pcos2x(2,i) + pcos2y(2,i) + pcos2z(2,i) &
                  + pcos2x(3,i) + pcos2y(3,i) + pcos2z(3,i))/3.0_KR
             mtsum(1,it) = mtsum(1,it) + fac1*elemloop(1,i)
             mtsum(2,it) = mtsum(2,it) + fac2*elemloop(1,i)
             mtsum(3,it) = mtsum(3,it) + fac3*elemloop(1,i)
             mtsum(4,it) = mtsum(4,it) + fac4*elemloop(1,i)
           enddo ! i
         else
           do i = 1,nvhalf
             it = itmin + 2*((i-1)/itstep) + (ibl-1)/8
             fac1 = (cosx(i) + cosy(i) + cosz(i) &
                  + pcosx(1,i) + pcosy(1,i) + pcosz(1,i) &
                  + pcosx(2,i) + pcosy(2,i) + pcosz(2,i) &
                  + pcosx(3,i) + pcosy(3,i) + pcosz(3,i))/3.0_KR
             fac2 = (cosx(i)*cosy(i) + cosx(i)*cosz(i) + cosy(i)*cosz(i) &
                  + pcosx(1,i)*pcosy(1,i) + pcosx(1,i)*pcosz(1,i) &
                  + pcosy(1,i)*pcosz(1,i) + pcosx(2,i)*pcosy(2,i) &
                  + pcosx(2,i)*pcosz(2,i) + pcosy(2,i)*pcosz(2,i) &
                  + pcosx(3,i)*pcosy(3,i) + pcosx(3,i)*pcosz(3,i) &
                  + pcosy(3,i)*pcosz(3,i))/3.0_KR
             fac3 = cosx(i)*cosy(i)*cosz(i) &
                  + pcosx(1,i)*pcosy(1,i)*pcosz(1,i) &
                  + pcosx(2,i)*pcosy(2,i)*pcosz(2,i) &
                  + pcosx(3,i)*pcosy(3,i)*pcosz(3,i)
             fac4 = (cos2x(i) + cos2y(i) + cos2z(i) &
                  + pcos2x(1,i) + pcos2y(1,i) + pcos2z(1,i) &
                  + pcos2x(2,i) + pcos2y(2,i) + pcos2z(2,i) &
                  + pcos2x(3,i) + pcos2y(3,i) + pcos2z(3,i))/3.0_KR
             mssum(1,it) = mssum(1,it) + fac1*elemloop(1,i)
             mssum(2,it) = mssum(2,it) + fac2*elemloop(1,i)
             mssum(3,it) = mssum(3,it) + fac3*elemloop(1,i)
             mssum(4,it) = mssum(4,it) + fac4*elemloop(1,i)
           enddo ! i
         endif

         ! The electric operator is J_4.
         ! The zero momentum result exactly vanishes when summed over all space points
         ! at a given time.
         ! Nonzero momentum:
         if (muvec(nu)==4) then
           do i = 1,nvhalf
             it = itmin + 2*((i-1)/itstep) + (ibl-1)/8
             fac1 = (cosx(i)+cosy(i)+cosz(i)-pcosx(1,i)-pcosy(1,i)-pcosz(1,i)) &
                  /3.0_KR
             fac2 = (cosx(i)*cosy(i)+cosx(i)*cosz(i)+cosy(i)*cosz(i) &
                  -pcosx(1,i)*pcosy(1,i)-pcosx(1,i)*pcosz(1,i) &
                  -pcosy(1,i)*pcosz(1,i))/3.0_KR
             fac3 = cosx(i)*cosy(i)*cosz(i) - pcosx(1,i)*pcosy(1,i)*pcosz(1,i)
             fac4 = (cos2x(i)+cos2y(i)+cos2z(i) &
                  -pcos2x(1,i)-pcos2y(1,i)-pcos2z(1,i))/3.0_KR
             j4sum(1,it) = j4sum(1,it) - fac1*elemloop(2,i)
             j4sum(2,it) = j4sum(2,it) - fac2*elemloop(2,i)
             j4sum(3,it) = j4sum(3,it) - fac3*elemloop(2,i)
             j4sum(4,it) = j4sum(4,it) - fac4*elemloop(2,i)
           enddo ! i
         endif

         ! The magnetic operator is built from J_1, J_2 and J_3.
         ! (Only nonzero momentum imaginary parts contribute).

         ! J_1 spacetime part:
         if (mu==1.and.muvec(nu)==4) then
           do i = 1,nvhalf
             it = itmin + 2*((i-1)/itstep) + (ibl-1)/8
             j1t(it) = j1t(it) + elemloop(2,i)
             fac1 = sinz(i)
             fac2 = sinz(i)*(cosx(i)+cosy(i))/2.0_KR
             fac3 = sinz(i)*cosx(i)*cosy(i)
             fac4 = sin2z(i)
             j13t(1,it) = j13t(1,it) + fac1*elemloop(2,i)
             j13t(2,it) = j13t(2,it) + fac2*elemloop(2,i)
             j13t(3,it) = j13t(3,it) + fac3*elemloop(2,i)
             j13t(4,it) = j13t(4,it) + fac4*elemloop(2,i)
             fac1 = siny(i)
             fac2 = siny(i)*(cosx(i)+cosz(i))/2.0_KR
             fac3 = siny(i)*cosx(i)*cosz(i)
             fac4 = sin2y(i)
             j12t(1,it) = j12t(1,it) + fac1*elemloop(2,i)
             j12t(2,it) = j12t(2,it) + fac2*elemloop(2,i)
             j12t(3,it) = j12t(3,it) + fac3*elemloop(2,i)
             j12t(4,it) = j12t(4,it) + fac4*elemloop(2,i)
           enddo ! i
         endif

         ! J_2 spacetime part:
         if (mu==2.and.muvec(nu)==4) then
           do i = 1,nvhalf
             it = itmin + 2*((i-1)/itstep) + (ibl-1)/8
             j2t(it) = j2t(it) + elemloop(2,i)
             fac1 = sinx(i)
             fac2 = sinx(i)*(cosy(i)+cosz(i))/2.0_KR
             fac3 = sinx(i)*cosy(i)*cosz(i)
             fac4 = sin2x(i)
             j21t(1,it) = j21t(1,it) + fac1*elemloop(2,i)
             j21t(2,it) = j21t(2,it) + fac2*elemloop(2,i)
             j21t(3,it) = j21t(3,it) + fac3*elemloop(2,i)
             j21t(4,it) = j21t(4,it) + fac4*elemloop(2,i)
             fac1 = sinz(i)
             fac2 = sinz(i)*(cosx(i)+cosy(i))/2.0_KR
             fac3 = sinz(i)*cosx(i)*cosy(i)
             fac4 = sin2z(i)
             j23t(1,it) = j23t(1,it) + fac1*elemloop(2,i)
             j23t(2,it) = j23t(2,it) + fac2*elemloop(2,i)
             j23t(3,it) = j23t(3,it) + fac3*elemloop(2,i)
             j23t(4,it) = j23t(4,it) + fac4*elemloop(2,i)
           enddo ! i
         endif

         ! J_3 spacetime part:
         if (mu==3.and.muvec(nu)==4) then
           do i = 1,nvhalf
             it = itmin + 2*((i-1)/itstep) + (ibl-1)/8
             j3t(it) = j3t(it) + elemloop(2,i)
             fac1 = siny(i)
             fac2 = siny(i)*(cosx(i)+cosz(i))/2.0_KR
             fac3 = siny(i)*cosx(i)*cosz(i)
             fac4 = sin2y(i)
             j32t(1,it) = j32t(1,it) + fac1*elemloop(2,i)
             j32t(2,it) = j32t(2,it) + fac2*elemloop(2,i)
             j32t(3,it) = j32t(3,it) + fac3*elemloop(2,i)
             j32t(4,it) = j32t(4,it) + fac4*elemloop(2,i)
             fac1 = sinx(i)
             fac2 = sinx(i)*(cosy(i)+cosz(i))/2.0_KR
             fac3 = sinx(i)*cosy(i)*cosz(i)
             fac4 = sin2x(i)
             j31t(1,it) = j31t(1,it) + fac1*elemloop(2,i)
             j31t(2,it) = j31t(2,it) + fac2*elemloop(2,i)
             j31t(3,it) = j31t(3,it) + fac3*elemloop(2,i)
             j31t(4,it) = j31t(4,it) + fac4*elemloop(2,i)
           enddo ! i
         endif

         ! J_1 spacespace parts:
         if (mu==1.and.muvec(nu)/=4) then
           do i = 1,nvhalf
             it = itmin + 2*((i-1)/itstep) + (ibl-1)/8
             fac1 = sinz(i) - psinz(2,i)
             fac2 = sinz(i)*(cosx(i)+cosy(i))/2.0_KR &
                  - psinz(2,i)*(pcosx(2,i)+pcosy(2,i))/2.0_KR
             fac3 = sinz(i)*cosx(i)*cosy(i) - psinz(2,i)*pcosx(2,i)*pcosy(2,i)
             fac4 = sin2z(i) - psin2z(2,i)
             j13s(1,it) = j13s(1,it) + fac1*elemloop(2,i)
             j13s(2,it) = j13s(2,it) + fac2*elemloop(2,i)
             j13s(3,it) = j13s(3,it) + fac3*elemloop(2,i)
             j13s(4,it) = j13s(4,it) + fac4*elemloop(2,i)
             fac1 = siny(i) - psiny(2,i)
             fac2 = siny(i)*(cosx(i)+cosz(i))/2.0_KR &
                  - psiny(2,i)*(pcosx(2,i)+pcosz(2,i))/2.0_KR
             fac3 = siny(i)*cosx(i)*cosz(i) - psiny(2,i)*pcosx(2,i)*pcosz(2,i)
             fac4 = sin2y(i) - psin2y(2,i)
             j12s(1,it) = j12s(1,it) + fac1*elemloop(2,i)
             j12s(2,it) = j12s(2,it) + fac2*elemloop(2,i)
             j12s(3,it) = j12s(3,it) + fac3*elemloop(2,i)
             j12s(4,it) = j12s(4,it) + fac4*elemloop(2,i)
           enddo ! i
         endif

         ! J_2 spacespace part:
         if (mu==2.and.muvec(nu)==3) then
           do i = 1,nvhalf
             it = itmin + 2*((i-1)/itstep) + (ibl-1)/8
             fac1 = sinx(i) - psinx(2,i)
             fac2 = sinx(i)*(cosy(i)+cosz(i))/2.0_KR &
                  - psinx(2,i)*(pcosy(2,i)+pcosz(2,i))/2.0_KR
             fac3 = sinx(i)*cosy(i)*cosz(i) - psinx(2,i)*pcosy(2,i)*pcosz(2,i)
             fac4 = sin2x(i) - psin2x(2,i)
             j21s(1,it) = j21s(1,it) + fac1*elemloop(2,i)
             j21s(2,it) = j21s(2,it) + fac2*elemloop(2,i)
             j21s(3,it) = j21s(3,it) + fac3*elemloop(2,i)
             j21s(4,it) = j21s(4,it) + fac4*elemloop(2,i)
             fac1 = sinz(i) - psinz(2,i)
             fac2 = sinz(i)*(cosx(i)+cosy(i))/2.0_KR &
                  - psinz(2,i)*(pcosx(2,i)+pcosy(2,i))/2.0_KR
             fac3 = sinz(i)*cosx(i)*cosy(i) - psinz(2,i)*pcosx(2,i)*pcosy(2,i)
             fac4 = sin2z(i) - psin2z(2,i)
             j23s(1,it) = j23s(1,it) + fac1*elemloop(2,i)
             j23s(2,it) = j23s(2,it) + fac2*elemloop(2,i)
             j23s(3,it) = j23s(3,it) + fac3*elemloop(2,i)
             j23s(4,it) = j23s(4,it) + fac4*elemloop(2,i)
           enddo ! i
         elseif (mu==1.and.muvec(nu)==2) then
           do i = 1,nvhalf
             it = itmin + 2*((i-1)/itstep) + (ibl-1)/8
             fac1 = sinx(i) - psinx(1,i)
             fac2 = sinx(i)*(cosy(i)+cosz(i))/2.0_KR &
                  - psinx(1,i)*(pcosy(1,i)+pcosz(1,i))/2.0_KR
             fac3 = sinx(i)*cosy(i)*cosz(i) - psinx(1,i)*pcosy(1,i)*pcosz(1,i)
             fac4 = sin2x(i) - psin2x(1,i)
             j21s(1,it) = j21s(1,it) - fac1*elemloop(2,i)
             j21s(2,it) = j21s(2,it) - fac2*elemloop(2,i)
             j21s(3,it) = j21s(3,it) - fac3*elemloop(2,i)
             j21s(4,it) = j21s(4,it) - fac4*elemloop(2,i)
             fac1 = sinz(i) - psinz(1,i)
             fac2 = sinz(i)*(cosx(i)+cosy(i))/2.0_KR &
                  - psinz(1,i)*(pcosx(1,i)+pcosy(1,i))/2.0_KR
             fac3 = sinz(i)*cosx(i)*cosy(i) - psinz(1,i)*pcosx(1,i)*pcosy(1,i)
             fac4 = sin2z(i) - psin2z(1,i)
             j23s(1,it) = j23s(1,it) - fac1*elemloop(2,i)
             j23s(2,it) = j23s(2,it) - fac2*elemloop(2,i)
             j23s(3,it) = j23s(3,it) - fac3*elemloop(2,i)
             j23s(4,it) = j23s(4,it) - fac4*elemloop(2,i)
           enddo ! i
         endif

         ! J_3 spacespace part:
         if (muvec(nu)==3) then
           do i = 1,nvhalf
             it = itmin + 2*((i-1)/itstep) + (ibl-1)/8
             fac1 = siny(i) - psiny(1,i)
             fac2 = siny(i)*(cosx(i)+cosz(i))/2.0_KR &
                  - psiny(1,i)*(pcosx(1,i)+pcosz(1,i))/2.0_KR
             fac3 = siny(i)*cosx(i)*cosz(i) - psiny(1,i)*pcosx(1,i)*pcosz(1,i)
             fac4 = sin2y(i) - psin2y(1,i)
             j32s(1,it) = j32s(1,it) - fac1*elemloop(2,i)
             j32s(2,it) = j32s(2,it) - fac2*elemloop(2,i)
             j32s(3,it) = j32s(3,it) - fac3*elemloop(2,i)
             j32s(4,it) = j32s(4,it) - fac4*elemloop(2,i)
             fac1 = sinx(i) - psinx(1,i)
             fac2 = sinx(i)*(cosy(i)+cosz(i))/2.0_KR &
                  - psinx(1,i)*(pcosy(1,i)+pcosz(1,i))/2.0_KR
             fac3 = sinx(i)*cosy(i)*cosz(i) - psinx(1,i)*pcosy(1,i)*pcosz(1,i)
             fac4 = sin2x(i) - psin2x(1,i)
             j31s(1,it) = j31s(1,it) - fac1*elemloop(2,i)
             j31s(2,it) = j31s(2,it) - fac2*elemloop(2,i)
             j31s(3,it) = j31s(3,it) - fac3*elemloop(2,i)
             j31s(4,it) = j31s(4,it) - fac4*elemloop(2,i)
           enddo ! i
         endif
       enddo ! nu

     enddo ! ibl
   enddo ! ieo
 enddo ! mu

 ! The factor of 8 is for the Dirac trace.
 factor = 8.0_KR*kappa**4
 ! Walter's factor of 1.0_KR/real(nx*ny*nz,KR) was removed from Jvev by
 ! removing it from "factor".
 do it = 1,nt
   plaq = plaq + tsum(it) + ssum(it)
 enddo ! it
 plaq = plaq/real(18*nxyzt,KR)
 tsum = tsum*4.0_KR*factor
 ssum = ssum*8.0_KR*factor
 j1t  = j1t*2.0_KR*factor
 j2t  = j2t*2.0_KR*factor
 j3t  = j3t*2.0_KR*factor
 ! The numerical factor besides "factor" is for the
 ! number of times a given space-space plaquette is counted at
 ! a given location in the operator psibar-psi.
 mtsum = mtsum*factor
 mssum = mssum*2.0_KR*factor
 j4sum = j4sum*2.0_KR*factor
 j12s = j12s*2.0_KR*factor
 j13s = j13s*2.0_KR*factor
 j21s = j21s*2.0_KR*factor
 j23s = j23s*2.0_KR*factor
 j31s = j31s*2.0_KR*factor
 j32s = j32s*2.0_KR*factor
 j12t = j12t*2.0_KR*factor
 j13t = j13t*2.0_KR*factor
 j21t = j21t*2.0_KR*factor
 j23t = j23t*2.0_KR*factor
 j31t = j31t*2.0_KR*factor
 j32t = j32t*2.0_KR*factor

 ! For interest's sake, compute <psibar psi>.
 ! Since Walter's factor of 1.0_KR/real(nx*ny*nz,KR) was explicitly removed
 ! from ssum and tsum (via "factor"), the 1.0_KR/real(nx*ny*nz,KR) factor is
 ! now explicitly multiplied into psi.
 psi = ssum(1) + tsum(1)
 do it = 2,nt
   psi = psi + ssum(it) + tsum(it) + tsum(it-1)
   if (myid==0) then
     print *, "it, psi(it)=", it, psi
   endif
 enddo ! it
 psi = psi/real(nxyzt,KR)

 ! Return temporal links at the maximal timestep to their true nonzero values.
if (.false.) then
 if (dobndry==1) then
   itbit = nvhalf-itstep
   do ibl = 9,16
     icount = ibl - 8
     do isite = itbit+1,nvhalf
       u(:,isite,4,:,ibl) = ubndry(:,isite-itbit,:,icount)
     enddo ! isite
   enddo ! ibl
 endif
endif

 ! Initialize Jvev.
 Jvev = 0.0_KR

 ! Load the local scalar current into Jvev.
 Jvev(1,1,1,5) = ssum(1) + tsum(1)
 do it = 2,nt
   Jvev(1,it,1,5) = ssum(it) + tsum(it) + tsum(it-1)
 enddo ! it
 do imom = 1,4
   Jvev(1,1,imom+1,5) = mssum(imom,1) + mtsum(imom,1)
   do it = 2,nt
     Jvev(1,it,imom+1,5) = mssum(imom,it) + mtsum(imom,it) + mtsum(imom,it-1)
   enddo ! it
 enddo ! imom

 ! Load the temporal component of the conserved vector current into Jvev.
 do imom = 1,4
   do it = 1,nt
     Jvev(2,it,imom+1,4) = j4sum(imom,it)
   enddo ! it
 enddo ! imom

 ! Load the spatial components of the conserved vector current into Jvev.
 Jvev(2,1,1,1) = j1t(1)
 Jvev(2,1,1,2) = j2t(1)
 Jvev(2,1,1,3) = j3t(1)
 do it = 2,nt
   Jvev(2,it,1,1) = j1t(it) - j1t(it-1)
   Jvev(2,it,1,2) = j2t(it) - j2t(it-1)
   Jvev(2,it,1,3) = j3t(it) - j3t(it-1)
 enddo ! it
 do imom = 1,4
   Jvev(1,1,imom+1,1) = j13s(imom,1) + j13t(imom,1)
   Jvev(2,1,imom+1,1) = j21s(imom,1) + j21t(imom,1)
   Jvev(1,1,imom+1,2) = j32s(imom,1) + j32t(imom,1)
   Jvev(2,1,imom+1,2) = j12s(imom,1) + j12t(imom,1)
   Jvev(1,1,imom+1,3) = j23s(imom,1) + j23t(imom,1)
   Jvev(2,1,imom+1,3) = j31s(imom,1) + j31t(imom,1)
   do it = 2,nt
     Jvev(1,it,imom+1,1) = j13s(imom,it) - j13t(imom,it-1) + j13t(imom,it)
     Jvev(2,it,imom+1,1) = j21s(imom,it) - j21t(imom,it-1) + j21t(imom,it)
     Jvev(1,it,imom+1,2) = j32s(imom,it) - j32t(imom,it-1) + j32t(imom,it)
     Jvev(2,it,imom+1,2) = j12s(imom,it) - j12t(imom,it-1) + j12t(imom,it)
     Jvev(1,it,imom+1,3) = j23s(imom,it) - j23t(imom,it-1) + j23t(imom,it)
     Jvev(2,it,imom+1,3) = j31s(imom,it) - j31t(imom,it-1) + j31t(imom,it)
   enddo ! it
 enddo ! imom

 ! Tack on an extra minus sign to agree with output of zing.f for example.
 Jvev = -Jvev

 ! Sum the results from all processes.
 if (nps/=1) then
   plaqtemp = plaq
   call MPI_REDUCE(plaqtemp,plaq,1,MRT,MPI_SUM,0,MPI_COMM_WORLD,ierr)
   call MPI_BCAST(plaq,1,MRT,0,MPI_COMM_WORLD,ierr)
   psitemp = psi
   call MPI_REDUCE(psitemp,psi,1,MRT,MPI_SUM,0,MPI_COMM_WORLD,ierr)
   call MPI_BCAST(psi,1,MRT,0,MPI_COMM_WORLD,ierr)
   Jvevtemp = Jvev
   nmesg = 2*nt*nmom*nop
   call MPI_REDUCE(Jvevtemp(1,1,1,1),Jvev(1,1,1,1),nmesg,MRT,MPI_SUM,0, &
                   MPI_COMM_WORLD,ierr)
   call MPI_BCAST(Jvev(1,1,1,1),nmesg,MRT,0,MPI_COMM_WORLD,ierr)
 endif

 ! Write some secondary results to file.
 ! (They are otherwise lost upon completion of this subroutine.)
 if (myid==0) then
   open(unit=8,file=trim(rwdir(myid+1))//"CFGSPROPS.LOG",     &
        action="write",form="formatted",status="old",position="append")
   write(unit=8,fmt="(a55,es17.10)") &
         "disconloops vev: the average plaquette (mod boundary) =", plaq
   write(unit=8,fmt="(a55,es17.10)") &
         "disconloops vev: <psibar psi> (mod boundary)          =", psi
   close(unit=8,status="keep")
 endif

 end subroutine vev

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 subroutine average(Jtime,nsub,nmom,nop,momfac,u,be,bo,xe,xo,kappa,dobndry, &
                    coact,bc,vecbl,vecblinv,myid,nn,ldiv,nms,lvbc,ib,lbd, &
                    iblv,rwdir,MRT)
! Purpose of this is simply to calculate averages of operators.
! Due to limitations in subroutine vev, local operators are corrected
! only to 5th order in kappa, and nonlocal ones to 4th order.
! INPUT:
!   nsub is the number of results, having different orders of subtraction,
!        to be considered (including the case of order=0).
!   nmom is the number of momenta to be considered (including momentum=0).
!   nop is the number of operators to be considered.
!   momfac(isite,ieo,ibl,imom,iop) are the momentum factors
!          where the first 3 entries are the lattice site
!          and the 4th entry is the momentum
!          and the last entry denotes one of the nop operators.
!   be() contains the source vector on even (gblclr=1) lattice sites.
!   bo() contains the source vector on odd (gblclr=2) lattice sites.
!        expected size: be(6,ntotal,4,2,8), bo(6,ntotal,4,2,8)
!                       where the first entry is real/imaginary and colour
!                       and the 3rd entry is the Dirac index
!                       and the other three entries give the lattice site.
!   xe() contains the sink vector on even (gblclr=1) lattice sites.
!   xo() contains the sink vector on odd (gblclr=2) lattice sites.
!        expected size: xe(6,ntotal,4,2,8), xo(6,nvhalf,4,2,8)
!                       where the first entry is real/imaginary and colour
!                       and the 3rd entry is the Dirac index
!                       and the other three entries give the lattice site.
! OUTPUT:
!   Jtime(iri,it,isub,imom,iop) are the averaged operators.
!        it=1..nt is the time slice on the full (multi-process) lattice.
!        isub=1..nsub, where nsub is the number of subtractions
!               considered, including "no subtractions".
!        imom=1..nmom, where nmom is the number of momenta considered.
!        For imom=1, iri=1..2 is the real/imaginary index and iop=1..nop.
!        For imom>1, iri=1..2 and iop=1..nop together label Im operators.

    real(kind=KR),    intent(out),   dimension(:,:,:,:,:) :: Jtime
    character(len=*), intent(in),    dimension(:)         :: rwdir
    real(kind=KR),    intent(in),    dimension(:,:,:,:,:) :: momfac
    integer(kind=KI), intent(in)         :: nsub, nmom, nop, dobndry, myid, MRT
    real(kind=KR),    intent(in)                         :: kappa
    real(kind=KR),    intent(inout), dimension(:,:,:,:,:) :: u
    real(kind=KR),    intent(inout), dimension(:,:,:,:,:) :: be, bo
    real(kind=KR),    intent(inout), dimension(:,:,:,:,:) :: xe, xo
    real(kind=KR),    intent(in),    dimension(:,:,:)     :: coact
    integer(kind=KI), intent(in),    dimension(:)         :: bc, nms
    integer(kind=KI), intent(in),    dimension(:,:)       :: vecbl, vecblinv
    integer(kind=KI), intent(in),    dimension(:,:)       :: nn, iblv
    logical,          intent(in),    dimension(:)         :: ldiv
    integer(kind=KI), intent(in),    dimension(:,:,:)     :: lvbc
    integer(kind=KI), intent(in),    dimension(:,:,:,:)   :: ib
    logical,          intent(in),    dimension(:,:)       :: lbd

    integer(kind=KI) :: gblclr, icri, id, isub, ksub, iri, isite, ibl, itbit, &
                        ibleo, ieo, idag
    real(kind=KR)                               :: xk
    real(kind=KR),    dimension(6,ntotal,4,2,8) :: sub1e, sub1o, sub2e, sub2o,&
                                                   sub3e, sub3o
    real(kind=KR),    dimension(6,nvhalf,4,2,8) :: sube, subo
    integer(kind=KI), parameter                 :: itstep=2*nvhalf*npt/nt
    real(kind=KR),    dimension(18,itstep,2,8)  :: ubndry
    !real(kind=KR),    dimension(2,nvhalf,5)     :: Je, Jo !5 should be replaced by nop
    real(kind=KR),    dimension(2,nvhalf,nop)     :: Je, Jo !5 should be replaced by nop
    integer(kind=KI) :: ierr, it, iop
    integer(kind=KI) :: icolor, idirac, mu, ibl1, gblclr1, jbl, jbleo, jeo

    ! This local gaugelink is used for debugging
    real(kind=KR),  dimension(18,ntotal,4,2,16) :: uout

    ! Initialization.
    idag = 0
    Jtime = 0.0_KR

    ! This routine is used for debugging by setting logical to .true.
    !      in if statement below.
    if (.false.) then ! debug
      call printlog("Debugging average",myid,rwdir)
      call printlog("**CHANGE u's BACK**",myid,rwdir)
      call fakegauge(uout,myid,rwdir,MRT)
      u = uout
    endif ! debug

    ! Set temporal links at the maximal timestep (on the global lattice) to zero.
    if (.false. .and. dobndry==1) then
      itbit = nvhalf - itstep
      do ibl = 9,16
        iri = ibl - 8
        do isite = itbit+1,nvhalf
          ubndry(:,isite-itbit,:,iri) = u(:,isite,4,:,ibl)
          u(:,isite,4,:,ibl) = 0.0_KR
        enddo ! isite
      enddo ! ibl
    endif

    ! Compute the operators for each level of subtraction.
    do isub = 1,nsub
      select case(isub)
        case(1) ! no subtraction
          ksub = 0
          sube = 0.0_KR
          subo = 0.0_KR
        case(2) ! subtraction of O(kappa), O(kappa^2) and O(kappa^3)
          ksub = 1
          gblclr = 1
          call Hsingle(sub1e,u,bo,idag,coact,bc,gblclr,vecbl,vecblinv,myid, &
                       nn,ldiv,nms,lvbc,ib,lbd,iblv,MRT)
          gblclr = 2
          call Hsingle(sub1o,u,be,idag,coact,bc,gblclr,vecbl,vecblinv,myid, &
                       nn,ldiv,nms,lvbc,ib,lbd,iblv,MRT)
          do ibleo = 1,8
            do ieo = 1,2
              do id = 1,4
                do isite = 1,nvhalf
                  do icri = 1,6
                    sube(icri,isite,id,ieo,ibleo) = be(icri,isite,id,ieo,ibleo) &
                                       + kappa*sub1e(icri,isite,id,ieo,ibleo)
                    subo(icri,isite,id,ieo,ibleo) = bo(icri,isite,id,ieo,ibleo) &
                                       + kappa*sub1o(icri,isite,id,ieo,ibleo)
                  enddo ! icri
                enddo ! isite
              enddo ! id
            enddo ! ieo
          enddo ! ibleo
          gblclr = 1
          call Hsingle(sub2e,u,sub1o,idag,coact,bc,gblclr,vecbl,vecblinv,myid, &
                       nn,ldiv,nms,lvbc,ib,lbd,iblv,MRT)
          gblclr = 2
          call Hsingle(sub2o,u,sub1e,idag,coact,bc,gblclr,vecbl,vecblinv,myid, &
                       nn,ldiv,nms,lvbc,ib,lbd,iblv,MRT)
          xk = kappa**2
          do ibleo = 1,8
            do ieo = 1,2
              do id = 1,4
                do isite = 1,nvhalf
                  do icri = 1,6
                    sube(icri,isite,id,ieo,ibleo) = sube(icri,isite,id,ieo,ibleo) &
                                              + xk*sub2e(icri,isite,id,ieo,ibleo)
                    subo(icri,isite,id,ieo,ibleo) = subo(icri,isite,id,ieo,ibleo) &
                                              + xk*sub2o(icri,isite,id,ieo,ibleo)
                  enddo ! icri
                enddo ! isite
              enddo ! id
            enddo ! ieo
          enddo ! ibleo
          gblclr = 1
          call Hsingle(sub1e,u,sub2o,idag,coact,bc,gblclr,vecbl,vecblinv,myid, &
                       nn,ldiv,nms,lvbc,ib,lbd,iblv,MRT)
          gblclr = 2
          call Hsingle(sub1o,u,sub2e,idag,coact,bc,gblclr,vecbl,vecblinv,myid, &
                       nn,ldiv,nms,lvbc,ib,lbd,iblv,MRT)
          xk = kappa**3
          do ibleo = 1,8
            do ieo = 1,2
              do id = 1,4
                do isite = 1,nvhalf
                  do icri = 1,6
                    sube(icri,isite,id,ieo,ibleo) = sube(icri,isite,id,ieo,ibleo) &
                                              + xk*sub1e(icri,isite,id,ieo,ibleo)
                    subo(icri,isite,id,ieo,ibleo) = subo(icri,isite,id,ieo,ibleo) &
                                              + xk*sub1o(icri,isite,id,ieo,ibleo)
                  enddo ! icri
                enddo ! isite
              enddo ! id
            enddo ! ieo
          enddo ! ibleo
        case(3) ! subtraction of O(kappa^4)
          gblclr = 1
          call Hsingle(sub2e,u,sub1o,idag,coact,bc,gblclr,vecbl,vecblinv,myid, &
                       nn,ldiv,nms,lvbc,ib,lbd,iblv,MRT)
          gblclr = 2
          call Hsingle(sub2o,u,sub1e,idag,coact,bc,gblclr,vecbl,vecblinv,myid, &
                       nn,ldiv,nms,lvbc,ib,lbd,iblv,MRT)
          xk = kappa**4
          do ibleo = 1,8
            do ieo = 1,2
              do id = 1,4
                do isite = 1,nvhalf
                  do icri = 1,6
                    sube(icri,isite,id,ieo,ibleo) = sube(icri,isite,id,ieo,ibleo) &
                                              + xk*sub2e(icri,isite,id,ieo,ibleo)
                    subo(icri,isite,id,ieo,ibleo) = subo(icri,isite,id,ieo,ibleo) &
                                              + xk*sub2o(icri,isite,id,ieo,ibleo)
                  enddo ! icri
                enddo ! isite
              enddo ! id
            enddo ! ieo
          enddo ! ibleo
        case(4) ! subtraction of O(kappa^5)
          gblclr = 1
          call Hsingle(sub1e,u,sub2o,idag,coact,bc,gblclr,vecbl,vecblinv,myid, &
                       nn,ldiv,nms,lvbc,ib,lbd,iblv,MRT)
          gblclr = 2
          call Hsingle(sub1o,u,sub2e,idag,coact,bc,gblclr,vecbl,vecblinv,myid, &
                       nn,ldiv,nms,lvbc,ib,lbd,iblv,MRT)
          xk = kappa**5
          do ibleo = 1,8
            do ieo = 1,2
              do id = 1,4
                do isite = 1,nvhalf
                  do icri = 1,6
                    sube(icri,isite,id,ieo,ibleo) = sube(icri,isite,id,ieo,ibleo) &
                                              + xk*sub1e(icri,isite,id,ieo,ibleo)
                    subo(icri,isite,id,ieo,ibleo) = subo(icri,isite,id,ieo,ibleo) &
                                              + xk*sub1o(icri,isite,id,ieo,ibleo)
                  enddo ! icri
                enddo ! isite
              enddo ! id
            enddo ! ieo
          enddo ! ibleo
        case(5) ! subtraction of O(kappa^6)!change here
          gblclr = 1
          call Hsingle(sub2e,u,sub1o,idag,coact,bc,gblclr,vecbl,vecblinv,myid, &
                       nn,ldiv,nms,lvbc,ib,lbd,iblv,MRT)
          gblclr = 2
          call Hsingle(sub2o,u,sub1e,idag,coact,bc,gblclr,vecbl,vecblinv,myid, &
                       nn,ldiv,nms,lvbc,ib,lbd,iblv,MRT)
          xk = kappa**6
          do ibleo = 1,8
            do ieo = 1,2
              do id = 1,4
                do isite = 1,nvhalf
                  do icri = 1,6
                    sube(icri,isite,id,ieo,ibleo) = sube(icri,isite,id,ieo,ibleo) &
                                              + xk*sub2e(icri,isite,id,ieo,ibleo) !BS 12/23/2015 sub1e changed to sub2e
                    subo(icri,isite,id,ieo,ibleo) = subo(icri,isite,id,ieo,ibleo) &
                                              + xk*sub2o(icri,isite,id,ieo,ibleo)  !BS 12/23/2015 sub1o changed to sub2o
                  enddo ! icri
                enddo ! isite
              enddo ! id
            enddo ! ieo
          enddo ! ibleo
        case(6) ! subtraction of O(kappa^7)
          gblclr = 1
          call Hsingle(sub1e,u,sub2o,idag,coact,bc,gblclr,vecbl,vecblinv,myid, &
                       nn,ldiv,nms,lvbc,ib,lbd,iblv,MRT)
          gblclr = 2
          call Hsingle(sub1o,u,sub2e,idag,coact,bc,gblclr,vecbl,vecblinv,myid, &
                       nn,ldiv,nms,lvbc,ib,lbd,iblv,MRT)
          xk = kappa**7
          do ibleo = 1,8
            do ieo = 1,2
              do id = 1,4
                do isite = 1,nvhalf
                  do icri = 1,6
                    sube(icri,isite,id,ieo,ibleo) = sube(icri,isite,id,ieo,ibleo) &
                                              + xk*sub1e(icri,isite,id,ieo,ibleo)
                    subo(icri,isite,id,ieo,ibleo) = subo(icri,isite,id,ieo,ibleo) &
                                              + xk*sub1o(icri,isite,id,ieo,ibleo)
                  enddo ! icri
                enddo ! isite
              enddo ! id
            enddo ! ieo
          enddo ! ibleo

!          gblclr = 1
!          call Hsingle(sub2e,u,sub1o,idag,coact,bc,gblclr,vecbl,vecblinv,myid, &
!                       nn,ldiv,nms,lvbc,ib,lbd,iblv,MRT)
!          gblclr = 2
!          call Hsingle(sub2o,u,sub1e,idag,coact,bc,gblclr,vecbl,vecblinv,myid, &
!                       nn,ldiv,nms,lvbc,ib,lbd,iblv,MRT)
!          xk = kappa**8
!          do ibleo = 1,8
!            do ieo = 1,2
!              do id = 1,4
!                do isite = 1,nvhalf
!                  do icri = 1,6
!                    sube(icri,isite,id,ieo,ibleo) = sube(icri,isite,id,ieo,ibleo) &
!                                              + xk*sub2e(icri,isite,id,ieo,ibleo)
!                    subo(icri,isite,id,ieo,ibleo) = subo(icri,isite,id,ieo,ibleo) &
!                                              + xk*sub2o(icri,isite,id,ieo,ibleo)
!                  enddo ! icri
!                enddo ! isite
!              enddo ! id
!            enddo ! ieo
!          enddo ! ibleo
!          gblclr = 1
!          call Hsingle(sub1e,u,sub2o,idag,coact,bc,gblclr,vecbl,vecblinv,myid, &
!                       nn,ldiv,nms,lvbc,ib,lbd,iblv,MRT)
!          gblclr = 2
!          call Hsingle(sub1o,u,sub2e,idag,coact,bc,gblclr,vecbl,vecblinv,myid, &
!                       nn,ldiv,nms,lvbc,ib,lbd,iblv,MRT)
!          xk = kappa**9
!          do ibleo = 1,8
!            do ieo = 1,2
!              do id = 1,4
!                do isite = 1,nvhalf
!                  do icri = 1,6
!                    sube(icri,isite,id,ieo,ibleo) = sube(icri,isite,id,ieo,ibleo) &
!                                              + xk*sub1e(icri,isite,id,ieo,ibleo)
!                    subo(icri,isite,id,ieo,ibleo) = subo(icri,isite,id,ieo,ibleo) &
!                                              + xk*sub1o(icri,isite,id,ieo,ibleo)
!                  enddo ! icri
!                enddo ! isite
!              enddo ! id
!            enddo ! ieo
!          enddo ! ibleo
!          gblclr = 1
!          call Hsingle(sub2e,u,sub1o,idag,coact,bc,gblclr,vecbl,vecblinv,myid, &
!                       nn,ldiv,nms,lvbc,ib,lbd,iblv,MRT)
!          gblclr = 2
!          call Hsingle(sub2o,u,sub1e,idag,coact,bc,gblclr,vecbl,vecblinv,myid, &
!                       nn,ldiv,nms,lvbc,ib,lbd,iblv,MRT)
!          xk = kappa**10
!          do ibleo = 1,8
!            do ieo = 1,2
!              do id = 1,4
!                do isite = 1,nvhalf
!                  do icri = 1,6
!                    sube(icri,isite,id,ieo,ibleo) = sube(icri,isite,id,ieo,ibleo) &
!                                              + xk*sub2e(icri,isite,id,ieo,ibleo)
!                    subo(icri,isite,id,ieo,ibleo) = subo(icri,isite,id,ieo,ibleo) &
!                                              + xk*sub2o(icri,isite,id,ieo,ibleo)
!                  enddo ! icri
!                enddo ! isite
!              enddo ! id
!            enddo ! ieo
!          enddo ! ibleo


        case default
          open(unit=8,file="DISCONLOOPS.ERROR",action="write",status="replace", &
               form="formatted")
          write(unit=8,fmt=*) "subroutine average: isub =", isub
          close(unit=8,status="keep")
          stop
      end select
      if (isub==2) then
        do ibleo = 1,8
          do ieo = 1,2
            do id = 1,4
              do isite = 1,nvhalf
                do icri = 1,6
                  sub2e(icri,isite,id,ieo,ibleo) = xe(icri,isite,id,ieo,ibleo) &
                                               - sube(icri,isite,id,ieo,ibleo)
                  sub2o(icri,isite,id,ieo,ibleo) = xo(icri,isite,id,ieo,ibleo) &
                                               - subo(icri,isite,id,ieo,ibleo)
                enddo ! icri
              enddo ! isite
            enddo ! id
          enddo ! ieo
        enddo ! ibleo
        do ibleo = 1,8
          do ieo = 1,2
            call loopops(ksub,sub2e,sub2o,u,be,bo,Je,Jo,ieo,ibleo,bc,vecbl, &
                         vecblinv,myid,nn,ldiv,nms,lvbc,ib,lbd,iblv,MRT)
            call spacesum(Jtime(:,:,isub,:,:),Je,Jo,ieo,ibleo,momfac,nmom,nop, &
                          myid,vecbl,vecblinv,MRT)
          enddo ! ieo
        enddo ! ibleo
      else
        do ibleo = 1,8
          do ieo = 1,2
            do id = 1,4
              do isite = 1,nvhalf
                do icri = 1,6
                  sub3e(icri,isite,id,ieo,ibleo) = xe(icri,isite,id,ieo,ibleo) &
                                               - sube(icri,isite,id,ieo,ibleo)
                  sub3o(icri,isite,id,ieo,ibleo) = xo(icri,isite,id,ieo,ibleo) &
                                               - subo(icri,isite,id,ieo,ibleo)
                enddo ! icri
              enddo ! isite
            enddo ! id
          enddo ! ieo
        enddo ! ibleo
        do ibleo = 1,8
          do ieo = 1,2
            call loopops(ksub,sub3e,sub3o,u,be,bo,Je,Jo,ieo,ibleo,bc,vecbl, &
                         vecblinv,myid,nn,ldiv,nms,lvbc,ib,lbd,iblv,MRT)
            call spacesum(Jtime(:,:,isub,:,:),Je,Jo,ieo,ibleo,momfac,nmom,nop, &
                          myid,vecbl,vecblinv,MRT)
          enddo ! ieo
        enddo ! ibleo
      endif
    enddo ! isub

    ! Return temporal links at the maximal timestep to their true nonzero values.
    if (.false. .and. dobndry==1) then
      do ibl = 9,16
        iri = ibl - 8
        do isite = itbit+1,nvhalf
          u(:,isite,4,:,ibl) = ubndry(:,isite-itbit,:,iri)
        enddo ! isite
      enddo ! ibl
    endif
    
    ! A final normalization factor.
    Jtime(:,:,:,:,1) = kappa*Jtime(:,:,:,:,1)
    Jtime(:,:,:,:,2) = kappa*Jtime(:,:,:,:,2)
    Jtime(:,:,:,:,3) = kappa*Jtime(:,:,:,:,3)
    Jtime(:,:,:,:,4) = kappa*Jtime(:,:,:,:,4)
    end subroutine average
!------------------------------------------------------------------------------ 

!------------------------------------------------------------------------------ 
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 subroutine ppaverage2(Jtime,nsub,nmom,nop,momfac,u,be,bo,xe,xo,kappa,dobndry, &
                    coact,bc,vecbl,vecblinv,myid,nn,ldiv,nms,lvbc,ib,lbd, &
                    iblv,rwdir,MRT,MRT2)
! Purpose of this is simply to calculate averages of operators.
! Due to limitations in subroutine vev, local operators are corrected
! only to 5th order in kappa, and nonlocal ones to 4th order.
! INPUT:
!   nsub is the number of results, having different orders of subtraction,
!        to be considered (including the case of order=0).
!   nmom is the number of momenta to be considered (including momentum=0).
!   nop is the number of operators to be considered.
!   momfac(isite,ieo,ibl,imom,iop) are the momentum factors
!          where the first 3 entries are the lattice site
!          and the 4th entry is the momentum
!          and the last entry denotes one of the nop operators.
!   be() contains the source vector on even (gblclr=1) lattice sites.
!   bo() contains the source vector on odd (gblclr=2) lattice sites.
!        expected size: be(6,ntotal,4,2,8), bo(6,ntotal,4,2,8)
!                       where the first entry is real/imaginary and colour
!                       and the 3rd entry is the Dirac index
!                       and the other three entries give the lattice site.
!   xe() contains the sink vector on even (gblclr=1) lattice sites.
!   xo() contains the sink vector on odd (gblclr=2) lattice sites.
!        expected size: xe(6,ntotal,4,2,8), xo(6,nvhalf,4,2,8)
!                       where the first entry is real/imaginary and colour
!                       and the 3rd entry is the Dirac index
!                       and the other three entries give the lattice site.
! OUTPUT:
!   Jtime(iri,it,isub,imom,iop) are the averaged operators.
!        it=1..nt is the time slice on the full (multi-process) lattice.
!        isub=1..nsub, where nsub is the number of subtractions
!               considered, including "no subtractions".
!        imom=1..nmom, where nmom is the number of momenta considered.
!        For imom=1, iri=1..2 is the real/imaginary index and iop=1..nop.
!        For imom>1, iri=1..2 and iop=1..nop together label Im operators.

    real(kind=KR),    intent(out),   dimension(:,:,:,:,:) :: Jtime
    character(len=*), intent(in),    dimension(:)         :: rwdir
    real(kind=KR),    intent(in),    dimension(:,:,:,:,:) :: momfac
    integer(kind=KI), intent(in)    :: nsub, nmom, nop, dobndry, myid, MRT, MRT2
    real(kind=KR),    intent(in)                         :: kappa
    real(kind=KR),    intent(inout), dimension(:,:,:,:,:) :: u
    real(kind=KR),    intent(inout), dimension(:,:,:,:,:) :: be, bo
    real(kind=KR),    intent(inout), dimension(:,:,:,:,:) :: xe, xo
    real(kind=KR),    intent(in),    dimension(:,:,:)     :: coact
    integer(kind=KI), intent(in),    dimension(:)         :: bc, nms
    integer(kind=KI), intent(in),    dimension(:,:)       :: vecbl, vecblinv
    integer(kind=KI), intent(in),    dimension(:,:)       :: nn, iblv
    logical,          intent(in),    dimension(:)         :: ldiv
    integer(kind=KI), intent(in),    dimension(:,:,:)     :: lvbc
    integer(kind=KI), intent(in),    dimension(:,:,:,:)   :: ib
    logical,          intent(in),    dimension(:,:)       :: lbd

    integer(kind=KI) :: gblclr, icri, id, isub, ksub, iri, isite, ibl, itbit, &
                        ibleo, ieo, idag, p, i, j
    real(kind=KR)                               :: xk, fac1, fac2
    real(kind=KR),    dimension(6,ntotal,4,2,8) :: sub1e, sub1o, sub2e, sub2o,&
                                                   vecsrc, temp, tue, tuo, &
                                                   sub3e, sub3o
    real(kind=KR),    dimension(6,nvhalf,4,2,8) :: sube, subo
    integer(kind=KI), parameter                 :: itstep=2*nvhalf*npt/nt
    real(kind=KR),    dimension(18,itstep,2,8)  :: ubndry
    !real(kind=KR),    dimension(2,nvhalf,5)     :: Je, Jo !5 should be replaced by nop
    real(kind=KR),    dimension(2,nvhalf,nop)     :: Je, Jo !5 should be replaced by nop
    integer(kind=KI) :: ierr, it, iop
    real(kind=KR2),   dimension(2)                          :: beta, beta1
    integer(kind=KI) :: icolor, idirac, mu, ibl1, gblclr1, jbl, jbleo, jeo

    ! This local gaugelink is used for debugging
    real(kind=KR),  dimension(18,ntotal,4,2,16) :: uout

    ! Initialization.
    idag = 0
    Jtime = 0.0_KR
    p = 11
    vprime = 0.0_KR2
    vprime1 = 0.0_KR2
    ! This routine is used for debugging by setting logical to .true.
    !      in if statement below.
    if (.false.) then ! debug
      call printlog("Debugging average",myid,rwdir)
      call printlog("**CHANGE u's BACK**",myid,rwdir)
      call fakegauge(uout,myid,rwdir,MRT)
      u = uout
    endif ! debug

    ! Set temporal links at the maximal timestep (on the global lattice) to zero.
    if (.false. .and. dobndry==1) then
      itbit = nvhalf - itstep
      do ibl = 9,16
        iri = ibl - 8
        do isite = itbit+1,nvhalf
          ubndry(:,isite-itbit,:,iri) = u(:,isite,4,:,ibl)
          u(:,isite,4,:,ibl) = 0.0_KR
        enddo ! isite
      enddo ! ibl
    endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!define the right hand side for xe
!      gblclr = 1
!      idag = 0
!      call Hsingle(temp,u,bo,idag,coact,bc,gblclr,vecbl, &
!                   vecblinv,myid,nn,ldiv,nms,lvbc,ib,lbd,iblv,MRT)

!      fac1 = 1.0_KR/kappa**2
!      fac2 = 1.0_KR/kappa
!      do isite = 1,nvhalf
!       vecsrc(:,isite,:,:,:) = fac1*be(:,isite,:,:,:) + fac2*temp(:,isite,:,:,:)
!      enddo ! isite
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Determine the polynomial for M.!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    do isite = 1,nvhalf
     vprime(:,isite,:,:,:,1) = be(:,isite,:,:,:)
    enddo !isite

    do isite = 1,nvhalf
     vprime1(:,isite,:,:,:,1) = bo(:,isite,:,:,:)
    enddo !isite

!    do isite = 1,nvhalf
!     tue(:,isite,:,:,:) = be(:,isite,:,:,:)
!     tuo(:,isite,:,:,:) = bo(:,isite,:,:,:)
!    enddo !isite


    do i = 1,p
          gblclr = 1
          call Hsingle(tue,u,vprime1(:,:,:,:,:,i),idag,coact,bc,gblclr,vecbl, &
                       vecblinv,myid,nn,ldiv,nms,lvbc,ib,lbd,iblv,MRT)
          gblclr = 2
           call Hsingle(tuo,u,vprime(:,:,:,:,:,i),idag,coact,bc,gblclr,vecbl, &
                       vecblinv,myid,nn,ldiv,nms,lvbc,ib,lbd,iblv,MRT)
!          do ibleo = 1,8
!            do ieo = 1,2
!              do id = 1,4
                do isite = 1,nvhalf
!                  do icri = 1,6
                   vprime(:,isite,:,:,:,i+1) = vprime(:,isite,:,:,:,i) &
                                               - kappa*tue(:,isite,:,:,:)
                   vprime1(:,isite,:,:,:,i+1) = vprime1(:,isite,:,:,:,i) &
                                                - kappa*tuo(:,isite,:,:,:)
!                  enddo ! icri
                enddo ! isite
!              enddo ! id
!            enddo ! ieo
!          enddo ! ibleo
 !      do isite = 1,nvhalf
!        vprime(:,isite,:,:,:,i+1) = tue(:,isite,:,:,:)
!       enddo !isite

!       do isite = 1,nvhalf
!        vprime1(:,isite,:,:,:,i+1) = tuo(:,isite,:,:,:)
!       enddo !isite
  enddo !i

  if(psignal==0) then!determent whether to calculate coeff for different rhs  
    do i=2,p+1
     do j=2,p+1
      call vecdot(vprime(:,:,:,:,:,i),vprime(:,:,:,:,:,j),beta,MRT2)
      call vecdot(vprime1(:,:,:,:,:,i),vprime1(:,:,:,:,:,j),beta1,MRT2)
      lsmat(1,i-1,j-1) = beta(1) + beta1(1)  !lsmat(2,p,p) ,cls(2,p,1)
      lsmat(2,i-1,j-1) = beta(2) + beta1(2)
     !BS if(myid==0) then
        !print *, "i,j, lsmat(:,i,j)=", i-1,j-1, lsmat(:,i-1,j-1)
      !endif
     enddo!j
    enddo!i
        



   do i=2,p+1
     call vecdot(vprime(:,:,:,:,:,i),vprime(:,:,:,:,:,1),beta,MRT2)
     call vecdot(vprime1(:,:,:,:,:,i),vprime1(:,:,:,:,:,1),beta1,MRT2)
     cls(1,i-1,1) = beta(1) + beta1(1)
     cls(2,i-1,1) = beta(2) + beta1(2)
    !BS if(myid==0) then
      ! print *, "i,cls(:,i)=", i-1, cls(:,i-1,1)
     !endif
   enddo!i
    
    call linearsolver(p,1,lsmat,ipiv2,cls)
    co(:,:) = cls(:,:,1)    
 endif!psignal
!    co = 0.0_KR2    
!    co(1,1) = 4
 !BS  if(myid==0) then
   ! do i=1,p
    ! print *, "i,result(:,i)=", i, co(:,i)
    !enddo!i  
   !endif!myid
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! Compute the operators for each level of subtraction.
 do isub = 1,nsub
   select case(isub)
     case(1) ! no subtraction
       ksub = 0
       sube = 0.0_KR
       subo = 0.0_KR
     case(2) ! subtraction of order1,order2 and order3
       ksub = 1
        do icri=1,5,2
         do isite=1,nvhalf
          sube(icri,isite,:,:,:) =  sube(icri,isite,:,:,:) &
                                    +co(1,1)*vprime(icri,isite,:,:,:,1) &
                                    -co(2,1)*vprime(icri+1,isite,:,:,:,1)
          sube(icri+1,isite,:,:,:) = sube(icri+1,isite,:,:,:) & 
                                     +co(1,1)*vprime(icri+1,isite,:,:,:,1) &
                                     +co(2,1)*vprime(icri,isite,:,:,:,1)

         enddo!isite
        enddo!icri

      do icri=1,5,2
       do isite=1,nvhalf
        subo(icri,isite,:,:,:) =  subo(icri,isite,:,:,:) &
                                  +co(1,1)*vprime1(icri,isite,:,:,:,1) &
                                  -co(2,1)*vprime1(icri+1,isite,:,:,:,1)
       subo(icri+1,isite,:,:,:) = subo(icri+1,isite,:,:,:) &
                                  +co(1,1)*vprime1(icri+1,isite,:,:,:,1) &
                                  +co(2,1)*vprime1(icri,isite,:,:,:,1)
       enddo!isite
      enddo!icri
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!order0!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        do icri=1,5,2
         do isite=1,nvhalf
          sube(icri,isite,:,:,:) =  sube(icri,isite,:,:,:) &
                                    +co(1,2)*vprime(icri,isite,:,:,:,2) &
                                    -co(2,2)*vprime(icri+1,isite,:,:,:,2)
          sube(icri+1,isite,:,:,:) = sube(icri+1,isite,:,:,:) & 
                                    +co(1,2)*vprime(icri+1,isite,:,:,:,2) &
                                    +co(2,2)*vprime(icri,isite,:,:,:,2)

         enddo!isite
        enddo!icri

      do icri=1,5,2
       do isite=1,nvhalf
        subo(icri,isite,:,:,:) =  subo(icri,isite,:,:,:) &
                                  +co(1,2)*vprime1(icri,isite,:,:,:,2) &
                                  -co(2,2)*vprime1(icri+1,isite,:,:,:,2)
        subo(icri+1,isite,:,:,:) = subo(icri+1,isite,:,:,:) &
                                   +co(1,2)*vprime1(icri+1,isite,:,:,:,2)&
                                   +co(2,2)*vprime1(icri,isite,:,:,:,2)
       enddo!isite
      enddo!icri
       


           
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1order!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        do icri=1,5,2
         do isite=1,nvhalf
          sube(icri,isite,:,:,:) =  sube(icri,isite,:,:,:) &
                                    +co(1,3)*vprime(icri,isite,:,:,:,3) &
                                    -co(2,3)*vprime(icri+1,isite,:,:,:,3)
          sube(icri+1,isite,:,:,:) = sube(icri+1,isite,:,:,:) & 
                                     +co(1,3)*vprime(icri+1,isite,:,:,:,3) &
                                     +co(2,3)*vprime(icri,isite,:,:,:,3)

         enddo!isite
        enddo!icri

      do icri=1,5,2
       do isite=1,nvhalf
        subo(icri,isite,:,:,:) =  subo(icri,isite,:,:,:) &
                                  +co(1,3)*vprime1(icri,isite,:,:,:,3) &
                                  -co(2,3)*vprime1(icri+1,isite,:,:,:,3)
        subo(icri+1,isite,:,:,:) = subo(icri+1,isite,:,:,:) &
                                  +co(1,3)*vprime1(icri+1,isite,:,:,:,3)&
                                  +co(2,3)*vprime1(icri,isite,:,:,:,3)
       enddo!isite
      enddo!icri

          
                                                  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!2order!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        do icri=1,5,2
         do isite=1,nvhalf
          sube(icri,isite,:,:,:) =  sube(icri,isite,:,:,:)&
                                   +co(1,4)*vprime(icri,isite,:,:,:,4) &
                                   -co(2,4)*vprime(icri+1,isite,:,:,:,4)
          sube(icri+1,isite,:,:,:) = sube(icri+1,isite,:,:,:)& 
                                    +co(1,4)*vprime(icri+1,isite,:,:,:,4) &
                                    +co(2,4)*vprime(icri,isite,:,:,:,4)

         enddo!isite
        enddo!icri

      do icri=1,5,2
       do isite=1,nvhalf
        subo(icri,isite,:,:,:) =  subo(icri,isite,:,:,:) &
                                  +co(1,4)*vprime1(icri,isite,:,:,:,4)&
                                  -co(2,4)*vprime1(icri+1,isite,:,:,:,4)
        subo(icri+1,isite,:,:,:) = subo(icri+1,isite,:,:,:) &
                                   +co(1,4)*vprime1(icri+1,isite,:,:,:,4)&
                                   +co(2,4)*vprime1(icri,isite,:,:,:,4)
       enddo!isite
      enddo!icri

              
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!3order!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      
        case(3) ! subtraction of O(kappa^4)
          
        do icri=1,5,2
         do isite=1,nvhalf
          sube(icri,isite,:,:,:) =  sube(icri,isite,:,:,:)&
                                   +co(1,5)*vprime(icri,isite,:,:,:,5) &
                                   -co(2,5)*vprime(icri+1,isite,:,:,:,5)
          sube(icri+1,isite,:,:,:) = sube(icri+1,isite,:,:,:)& 
                                    +co(1,5)*vprime(icri+1,isite,:,:,:,5) &
                                    +co(2,5)*vprime(icri,isite,:,:,:,5)

         enddo!isite
        enddo!icri

      do icri=1,5,2
       do isite=1,nvhalf
        subo(icri,isite,:,:,:) =  subo(icri,isite,:,:,:) &
                                  +co(1,5)*vprime1(icri,isite,:,:,:,5) &
                                  -co(2,5)*vprime1(icri+1,isite,:,:,:,5)
        subo(icri+1,isite,:,:,:) = subo(icri+1,isite,:,:,:) &
                                   +co(1,5)*vprime1(icri+1,isite,:,:,:,5)&
                                   +co(2,5)*vprime1(icri,isite,:,:,:,5)
       enddo!isite
      enddo!icri

            
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!4order!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        case(4) ! subtraction of O(kappa^5)
         do icri=1,5,2
         do isite=1,nvhalf
          sube(icri,isite,:,:,:) =  sube(icri,isite,:,:,:)&
                                   +co(1,6)*vprime(icri,isite,:,:,:,6) &
                                   -co(2,6)*vprime(icri+1,isite,:,:,:,6)
          sube(icri+1,isite,:,:,:) = sube(icri+1,isite,:,:,:)& 
                                    +co(1,6)*vprime(icri+1,isite,:,:,:,6) &
                                    +co(2,6)*vprime(icri,isite,:,:,:,6)

         enddo!isite
        enddo!icri

      do icri=1,5,2
       do isite=1,nvhalf
        subo(icri,isite,:,:,:) =  subo(icri,isite,:,:,:) &
                                  +co(1,6)*vprime1(icri,isite,:,:,:,6) &
                                  -co(2,6)*vprime1(icri+1,isite,:,:,:,6)
        subo(icri+1,isite,:,:,:) = subo(icri+1,isite,:,:,:) &
                                   +co(1,6)*vprime1(icri+1,isite,:,:,:,6)&
                                   +co(2,6)*vprime1(icri,isite,:,:,:,6)
       enddo!isite
      enddo!icri
        
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!5order!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
        case(5) ! subtraction of O(kappa^6)
        do icri=1,5,2
         do isite=1,nvhalf
          sube(icri,isite,:,:,:) =  sube(icri,isite,:,:,:)&
                                   +co(1,7)*vprime(icri,isite,:,:,:,7) &
                                   -co(2,7)*vprime(icri+1,isite,:,:,:,7)
          sube(icri+1,isite,:,:,:) = sube(icri+1,isite,:,:,:)& 
                                    +co(1,7)*vprime(icri+1,isite,:,:,:,7) &
                                    +co(2,7)*vprime(icri,isite,:,:,:,7)

         enddo!isite
        enddo!icri

      do icri=1,5,2
       do isite=1,nvhalf
        subo(icri,isite,:,:,:) =  subo(icri,isite,:,:,:) &
                                  +co(1,7)*vprime1(icri,isite,:,:,:,7) &
                                  -co(2,7)*vprime1(icri+1,isite,:,:,:,7)
        subo(icri+1,isite,:,:,:) = subo(icri+1,isite,:,:,:) &
                                   +co(1,7)*vprime1(icri+1,isite,:,:,:,7)&
                                   +co(2,7)*vprime1(icri,isite,:,:,:,7)
       enddo!isite
      enddo!icri

    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!6order!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
        case(6) ! subtraction of O(kappa^7)
        do icri=1,5,2
         do isite=1,nvhalf
          sube(icri,isite,:,:,:) =  sube(icri,isite,:,:,:)&
                                   +co(1,8)*vprime(icri,isite,:,:,:,8) &
                                   -co(2,8)*vprime(icri+1,isite,:,:,:,8)
          sube(icri+1,isite,:,:,:) = sube(icri+1,isite,:,:,:)& 
                                    +co(1,8)*vprime(icri+1,isite,:,:,:,8) &
                                    +co(2,8)*vprime(icri,isite,:,:,:,8)

         enddo!isite
        enddo!icri

      do icri=1,5,2
       do isite=1,nvhalf
        subo(icri,isite,:,:,:) =  subo(icri,isite,:,:,:) &
                                  +co(1,8)*vprime1(icri,isite,:,:,:,8) &
                                  -co(2,8)*vprime1(icri+1,isite,:,:,:,8)
        subo(icri+1,isite,:,:,:) = subo(icri+1,isite,:,:,:) &
                                   +co(1,8)*vprime1(icri+1,isite,:,:,:,8)&
                                   +co(2,8)*vprime1(icri,isite,:,:,:,8)
       enddo!isite
      enddo!icri

         do icri=1,5,2
         do isite=1,nvhalf
          sube(icri,isite,:,:,:) =  sube(icri,isite,:,:,:)&
                                   +co(1,9)*vprime(icri,isite,:,:,:,9) &
                                   -co(2,9)*vprime(icri+1,isite,:,:,:,9)
          sube(icri+1,isite,:,:,:) = sube(icri+1,isite,:,:,:)& 
                                    +co(1,9)*vprime(icri+1,isite,:,:,:,9) &
                                    +co(2,9)*vprime(icri,isite,:,:,:,9)

         enddo!isite
        enddo!icri

      do icri=1,5,2
       do isite=1,nvhalf
        subo(icri,isite,:,:,:) =  subo(icri,isite,:,:,:) &
                                  +co(1,9)*vprime1(icri,isite,:,:,:,9) &
                                  -co(2,9)*vprime1(icri+1,isite,:,:,:,9)
        subo(icri+1,isite,:,:,:) = subo(icri+1,isite,:,:,:) &
                                   +co(1,9)*vprime1(icri+1,isite,:,:,:,9)&
                                   +co(2,9)*vprime1(icri,isite,:,:,:,9)
       enddo!isite
      enddo!icri
        
         do icri=1,5,2
         do isite=1,nvhalf
          sube(icri,isite,:,:,:) =  sube(icri,isite,:,:,:)&
                                   +co(1,10)*vprime(icri,isite,:,:,:,10) &
                                   -co(2,10)*vprime(icri+1,isite,:,:,:,10)
          sube(icri+1,isite,:,:,:) = sube(icri+1,isite,:,:,:)& 
                                    +co(1,10)*vprime(icri+1,isite,:,:,:,10) &
                                    +co(2,10)*vprime(icri,isite,:,:,:,10)

         enddo!isite
        enddo!icri

      do icri=1,5,2
       do isite=1,nvhalf
        subo(icri,isite,:,:,:) =  subo(icri,isite,:,:,:) &
                                  +co(1,10)*vprime1(icri,isite,:,:,:,10) &
                                  -co(2,10)*vprime1(icri+1,isite,:,:,:,10)
        subo(icri+1,isite,:,:,:) = subo(icri+1,isite,:,:,:) &
                                   +co(1,10)*vprime1(icri+1,isite,:,:,:,10)&
                                   +co(2,10)*vprime1(icri,isite,:,:,:,10)
       enddo!isite
      enddo!icri
                                       
         do icri=1,5,2
         do isite=1,nvhalf
          sube(icri,isite,:,:,:) =  sube(icri,isite,:,:,:)&
                                   +co(1,11)*vprime(icri,isite,:,:,:,11) &
                                   -co(2,11)*vprime(icri+1,isite,:,:,:,11)
          sube(icri+1,isite,:,:,:) = sube(icri+1,isite,:,:,:)& 
                                    +co(1,11)*vprime(icri+1,isite,:,:,:,11) &
                                    +co(2,11)*vprime(icri,isite,:,:,:,11)

         enddo!isite
        enddo!icri

      do icri=1,5,2
       do isite=1,nvhalf
        subo(icri,isite,:,:,:) =  subo(icri,isite,:,:,:) &
                                  +co(1,11)*vprime1(icri,isite,:,:,:,11) &
                                  -co(2,11)*vprime1(icri+1,isite,:,:,:,11)
        subo(icri+1,isite,:,:,:) = subo(icri+1,isite,:,:,:) &
                                   +co(1,11)*vprime1(icri+1,isite,:,:,:,11)&
                                   +co(2,11)*vprime1(icri,isite,:,:,:,11)
       enddo!isite
      enddo!icri
                                       

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!7order!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        case default
          open(unit=8,file="DISCONLOOPS.ERROR",action="write",status="replace", &
               form="formatted")
          write(unit=8,fmt=*) "subroutine average: isub =", isub
          close(unit=8,status="keep")
          stop
      end select
      if (isub==2) then
        do ibleo = 1,8
          do ieo = 1,2
            do id = 1,4
              do isite = 1,nvhalf
                do icri = 1,6
                  sub2e(icri,isite,id,ieo,ibleo) = xe(icri,isite,id,ieo,ibleo) &
                                               - sube(icri,isite,id,ieo,ibleo)
                  sub2o(icri,isite,id,ieo,ibleo) = xo(icri,isite,id,ieo,ibleo) &
                                               - subo(icri,isite,id,ieo,ibleo)
                enddo ! icri
              enddo ! isite
            enddo ! id
          enddo ! ieo
        enddo ! ibleo
        do ibleo = 1,8
          do ieo = 1,2
            call loopops(ksub,sub2e,sub2o,u,be,bo,Je,Jo,ieo,ibleo,bc,vecbl, &
                         vecblinv,myid,nn,ldiv,nms,lvbc,ib,lbd,iblv,MRT)
            call spacesum(Jtime(:,:,isub,:,:),Je,Jo,ieo,ibleo,momfac,nmom,nop, &
                          myid,vecbl,vecblinv,MRT)
          enddo ! ieo
        enddo ! ibleo
      else
        do ibleo = 1,8
          do ieo = 1,2
            do id = 1,4
              do isite = 1,nvhalf
                do icri = 1,6
                  sub1e(icri,isite,id,ieo,ibleo) = xe(icri,isite,id,ieo,ibleo) &
                                               - sube(icri,isite,id,ieo,ibleo)
                  sub1o(icri,isite,id,ieo,ibleo) = xo(icri,isite,id,ieo,ibleo) &
                                               - subo(icri,isite,id,ieo,ibleo)
                enddo ! icri
              enddo ! isite
            enddo ! id
          enddo ! ieo
        enddo ! ibleo
        do ibleo = 1,8
          do ieo = 1,2
            call loopops(ksub,sub1e,sub1o,u,be,bo,Je,Jo,ieo,ibleo,bc,vecbl, &
                         vecblinv,myid,nn,ldiv,nms,lvbc,ib,lbd,iblv,MRT)
            call spacesum(Jtime(:,:,isub,:,:),Je,Jo,ieo,ibleo,momfac,nmom,nop, &
                          myid,vecbl,vecblinv,MRT)
          enddo ! ieo
        enddo ! ibleo
      endif
    enddo ! isub

    ! Return temporal links at the maximal timestep to their true nonzero values.
    if (.false. .and. dobndry==1) then
      do ibl = 9,16
        iri = ibl - 8
        do isite = itbit+1,nvhalf
          u(:,isite,4,:,ibl) = ubndry(:,isite-itbit,:,iri)
        enddo ! isite
      enddo ! ibl
    endif
    
    ! A final normalization factor.
    Jtime(:,:,:,:,1) = kappa*Jtime(:,:,:,:,1)
    Jtime(:,:,:,:,2) = kappa*Jtime(:,:,:,:,2)
    Jtime(:,:,:,:,3) = kappa*Jtime(:,:,:,:,3)
    Jtime(:,:,:,:,4) = kappa*Jtime(:,:,:,:,4)
    end subroutine ppaverage2
!---------------------------------------------------------------------------
 subroutine ppaverage1(Jtime,nsub,nmom,nop,momfac,u,be,bo,xe,xo,kappa,dobndry, &
                    coact,bc,vecbl,vecblinv,myid,nn,ldiv,nms,lvbc,ib,lbd, &
                    iblv,rwdir,MRT,MRT2)
!This is the 7th order subtraction.
! Purpose of this is simply to calculate averages of operators.
! Due to limitations in subroutine vev, local operators are corrected
! only to 5th order in kappa, and nonlocal ones to 4th order.
! INPUT:
!   nsub is the number of results, having different orders of subtraction,
!        to be considered (including the case of order=0).
!   nmom is the number of momenta to be considered (including momentum=0).
!   nop is the number of operators to be considered.
!   momfac(isite,ieo,ibl,imom,iop) are the momentum factors
!          where the first 3 entries are the lattice site
!          and the 4th entry is the momentum
!          and the last entry denotes one of the nop operators.
!   be() contains the source vector on even (gblclr=1) lattice sites.
!   bo() contains the source vector on odd (gblclr=2) lattice sites.
!        expected size: be(6,ntotal,4,2,8), bo(6,ntotal,4,2,8)
!                       where the first entry is real/imaginary and colour
!                       and the 3rd entry is the Dirac index
!                       and the other three entries give the lattice site.
!   xe() contains the sink vector on even (gblclr=1) lattice sites.
!   xo() contains the sink vector on odd (gblclr=2) lattice sites.
!        expected size: xe(6,ntotal,4,2,8), xo(6,nvhalf,4,2,8)
!                       where the first entry is real/imaginary and colour
!                       and the 3rd entry is the Dirac index
!                       and the other three entries give the lattice site.
! OUTPUT:
!   Jtime(iri,it,isub,imom,iop) are the averaged operators.
!        it=1..nt is the time slice on the full (multi-process) lattice.
!        isub=1..nsub, where nsub is the number of subtractions
!               considered, including "no subtractions".
!        imom=1..nmom, where nmom is the number of momenta considered.
!        For imom=1, iri=1..2 is the real/imaginary index and iop=1..nop.
!        For imom>1, iri=1..2 and iop=1..nop together label Im operators.

    real(kind=KR),    intent(out),   dimension(:,:,:,:,:) :: Jtime
    character(len=*), intent(in),    dimension(:)         :: rwdir
    real(kind=KR),    intent(in),    dimension(:,:,:,:,:) :: momfac
    integer(kind=KI), intent(in)    :: nsub, nmom, nop, dobndry, myid, MRT, MRT2
    real(kind=KR),    intent(in)                         :: kappa
    real(kind=KR),    intent(inout), dimension(:,:,:,:,:) :: u
    real(kind=KR),    intent(inout), dimension(:,:,:,:,:) :: be, bo
    real(kind=KR),    intent(inout), dimension(:,:,:,:,:) :: xe, xo
    real(kind=KR),    intent(in),    dimension(:,:,:)     :: coact
    integer(kind=KI), intent(in),    dimension(:)         :: bc, nms
    integer(kind=KI), intent(in),    dimension(:,:)       :: vecbl, vecblinv
    integer(kind=KI), intent(in),    dimension(:,:)       :: nn, iblv
    logical,          intent(in),    dimension(:)         :: ldiv
    integer(kind=KI), intent(in),    dimension(:,:,:)     :: lvbc
    integer(kind=KI), intent(in),    dimension(:,:,:,:)   :: ib
    logical,          intent(in),    dimension(:,:)       :: lbd

    integer(kind=KI) :: gblclr, icri, id, isub, ksub, iri, isite, ibl, itbit, &
                        ibleo, ieo, idag, p, i, j
    real(kind=KR)                               :: xk, fac1, fac2
    real(kind=KR),    dimension(6,ntotal,4,2,8) :: sub1e, sub1o, sub2e, sub2o,&
                                                   vecsrc, temp, tue, tuo, &
                                                   sub3e, sub3o
    real(kind=KR),    dimension(6,nvhalf,4,2,8) :: sube, subo
    integer(kind=KI), parameter                 :: itstep=2*nvhalf*npt/nt
    real(kind=KR),    dimension(18,itstep,2,8)  :: ubndry
    !real(kind=KR),    dimension(2,nvhalf,5)     :: Je, Jo !5 should be replaced by nop
    real(kind=KR),    dimension(2,nvhalf,nop)     :: Je, Jo !5 should be replaced by nop
    integer(kind=KI) :: ierr, it, iop
    real(kind=KR2),   dimension(2)                          :: beta, beta1
    integer(kind=KI) :: icolor, idirac, mu, ibl1, gblclr1, jbl, jbleo, jeo

    ! This local gaugelink is used for debugging
    real(kind=KR),  dimension(18,ntotal,4,2,16) :: uout
    integer(kind=KI) :: exists!BS 3/2/2016
    character(len=128) :: startingvector !BS 3/2/2016

    ! Initialization.
    idag = 0
    Jtime = 0.0_KR
    p = 8
    vprime = 0.0_KR2
    vprime1 = 0.0_KR2
    ! This routine is used for debugging by setting logical to .true.
    !      in if statement below.
    if (.false.) then ! debug
      call printlog("Debugging average",myid,rwdir)
      call printlog("**CHANGE u's BACK**",myid,rwdir)
      call fakegauge(uout,myid,rwdir,MRT)
      u = uout
    endif ! debug

    ! Set temporal links at the maximal timestep (on the global lattice) to zero.
    if (.false. .and. dobndry==1) then
      itbit = nvhalf - itstep
      do ibl = 9,16
        iri = ibl - 8
        do isite = itbit+1,nvhalf
          ubndry(:,isite-itbit,:,iri) = u(:,isite,4,:,ibl)
          u(:,isite,4,:,ibl) = 0.0_KR
        enddo ! isite
      enddo ! ibl
    endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!define the right hand side for xe
!      gblclr = 1
!      idag = 0
!      call Hsingle(temp,u,bo,idag,coact,bc,gblclr,vecbl, &
!                   vecblinv,myid,nn,ldiv,nms,lvbc,ib,lbd,iblv,MRT)

!      fac1 = 1.0_KR/kappa**2
!      fac2 = 1.0_KR/kappa
!      do isite = 1,nvhalf
!       vecsrc(:,isite,:,:,:) = fac1*be(:,isite,:,:,:) + fac2*temp(:,isite,:,:,:)
!      enddo ! isite
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Determine the polynomial for M.!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!BS added the following lines to printout starting vector for matlab

if (.false.) then

write(startingvector,"(a14,i3)")"startingvector"

if (myid==0) then
        inquire(file=trim(rwdir(myid+1))//trim(startingvector)//".LOG",exist=exists)
  if (.not. exists) then

         print *, "File does not exist. Creating it inside noise file."
open(unit=80,file=trim(rwdir(myid+1))//trim(startingvector)//".LOG",status="new",&
        action="write",form="formatted")
   close(unit=80,status="keep")
end if

end if

 call printferm(be,bo,rwdir,startingvector,myid,iblv,vecblinv,MRT2)!BS 3/2/2016


endif!true or false



















    do isite = 1,nvhalf
     vprime(:,isite,:,:,:,1) = be(:,isite,:,:,:)
    enddo !isite

    do isite = 1,nvhalf
     vprime1(:,isite,:,:,:,1) = bo(:,isite,:,:,:)
    enddo !isite

!    do isite = 1,nvhalf
!     tue(:,isite,:,:,:) = be(:,isite,:,:,:)
!     tuo(:,isite,:,:,:) = bo(:,isite,:,:,:)
!    enddo !isite


    do i = 1,p
          gblclr = 1
          call Hsingle(tue,u,vprime1(:,:,:,:,:,i),idag,coact,bc,gblclr,vecbl, &
                       vecblinv,myid,nn,ldiv,nms,lvbc,ib,lbd,iblv,MRT)
          gblclr = 2
           call Hsingle(tuo,u,vprime(:,:,:,:,:,i),idag,coact,bc,gblclr,vecbl, &
                       vecblinv,myid,nn,ldiv,nms,lvbc,ib,lbd,iblv,MRT)
!          do ibleo = 1,8
!            do ieo = 1,2
!              do id = 1,4
                do isite = 1,nvhalf
!                  do icri = 1,6
                   vprime(:,isite,:,:,:,i+1) = vprime(:,isite,:,:,:,i) &
                                               - kappa*tue(:,isite,:,:,:)
                   vprime1(:,isite,:,:,:,i+1) = vprime1(:,isite,:,:,:,i) &
                                                - kappa*tuo(:,isite,:,:,:)
!                  enddo ! icri
                enddo ! isite
!              enddo ! id
!            enddo ! ieo
!          enddo ! ibleo
 !      do isite = 1,nvhalf
!        vprime(:,isite,:,:,:,i+1) = tue(:,isite,:,:,:)
!       enddo !isite

!       do isite = 1,nvhalf
!        vprime1(:,isite,:,:,:,i+1) = tuo(:,isite,:,:,:)
!       enddo !isite
  enddo !i

  if(psignal==0) then!determent whether to calculate coeff for different rhs  
    do i=2,p+1
     do j=2,p+1
      call vecdot(vprime(:,:,:,:,:,i),vprime(:,:,:,:,:,j),beta,MRT2)
      call vecdot(vprime1(:,:,:,:,:,i),vprime1(:,:,:,:,:,j),beta1,MRT2)
      lsmat1(1,i-1,j-1) = beta(1) + beta1(1)  !lsmat(2,p,p) ,cls(2,p,1)
      lsmat1(2,i-1,j-1) = beta(2) + beta1(2)
!      if(myid==0) then
!        print *, "i,j, lsmat(:,i,j)=", i-1,j-1, lsmat(:,i-1,j-1)
!      endif
     enddo!j
    enddo!i
        



   do i=2,p+1
     call vecdot(vprime(:,:,:,:,:,i),vprime(:,:,:,:,:,1),beta,MRT2)
     call vecdot(vprime1(:,:,:,:,:,i),vprime1(:,:,:,:,:,1),beta1,MRT2)
     cls1(1,i-1,1) = beta(1) + beta1(1)
     cls1(2,i-1,1) = beta(2) + beta1(2)
!     if(myid==0) then
!       print *, "i,cls(:,i)=", i-1, cls(:,i-1,1)
!     endif
   enddo!i
    
    call linearsolver(p,1,lsmat1,ipiv3,cls1)
    co2(:,:) = cls1(:,:,1)    
 endif!psignal
!    co = 0.0_KR2    
!    co(1,1) = 4
  !BS if(myid==0) then
    !do i=1,p
    ! print *, "i,result(:,i)=", i, co2(:,i)
    !enddo!i  
   !endif!myid
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! Compute the operators for each level of subtraction.
 do isub = 1,nsub
   select case(isub)
     case(1) ! no subtraction
       ksub = 0
       sube = 0.0_KR
       subo = 0.0_KR
     case(2) ! subtraction of order1,order2 and order3
       ksub = 1
        do icri=1,5,2
         do isite=1,nvhalf
          sube(icri,isite,:,:,:) =  sube(icri,isite,:,:,:) &
                                    +co2(1,1)*vprime(icri,isite,:,:,:,1) &
                                    -co2(2,1)*vprime(icri+1,isite,:,:,:,1)
          sube(icri+1,isite,:,:,:) = sube(icri+1,isite,:,:,:) & 
                                     +co2(1,1)*vprime(icri+1,isite,:,:,:,1) &
                                     +co2(2,1)*vprime(icri,isite,:,:,:,1)

         enddo!isite
        enddo!icri

      do icri=1,5,2
       do isite=1,nvhalf
        subo(icri,isite,:,:,:) =  subo(icri,isite,:,:,:) &
                                  +co2(1,1)*vprime1(icri,isite,:,:,:,1) &
                                  -co2(2,1)*vprime1(icri+1,isite,:,:,:,1)
       subo(icri+1,isite,:,:,:) = subo(icri+1,isite,:,:,:) &
                                  +co2(1,1)*vprime1(icri+1,isite,:,:,:,1) &
                                  +co2(2,1)*vprime1(icri,isite,:,:,:,1)
       enddo!isite
      enddo!icri
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!order0!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        do icri=1,5,2
         do isite=1,nvhalf
          sube(icri,isite,:,:,:) =  sube(icri,isite,:,:,:) &
                                    +co2(1,2)*vprime(icri,isite,:,:,:,2) &
                                    -co2(2,2)*vprime(icri+1,isite,:,:,:,2)
          sube(icri+1,isite,:,:,:) = sube(icri+1,isite,:,:,:) & 
                                    +co2(1,2)*vprime(icri+1,isite,:,:,:,2) &
                                    +co2(2,2)*vprime(icri,isite,:,:,:,2)

         enddo!isite
        enddo!icri

      do icri=1,5,2
       do isite=1,nvhalf
        subo(icri,isite,:,:,:) =  subo(icri,isite,:,:,:) &
                                  +co2(1,2)*vprime1(icri,isite,:,:,:,2) &
                                  -co2(2,2)*vprime1(icri+1,isite,:,:,:,2)
        subo(icri+1,isite,:,:,:) = subo(icri+1,isite,:,:,:) &
                                   +co2(1,2)*vprime1(icri+1,isite,:,:,:,2)&
                                   +co2(2,2)*vprime1(icri,isite,:,:,:,2)
       enddo!isite
      enddo!icri
       


           
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1order!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        do icri=1,5,2
         do isite=1,nvhalf
          sube(icri,isite,:,:,:) =  sube(icri,isite,:,:,:) &
                                    +co2(1,3)*vprime(icri,isite,:,:,:,3) &
                                    -co2(2,3)*vprime(icri+1,isite,:,:,:,3)
          sube(icri+1,isite,:,:,:) = sube(icri+1,isite,:,:,:) & 
                                     +co2(1,3)*vprime(icri+1,isite,:,:,:,3) &
                                     +co2(2,3)*vprime(icri,isite,:,:,:,3)

         enddo!isite
        enddo!icri

      do icri=1,5,2
       do isite=1,nvhalf
        subo(icri,isite,:,:,:) =  subo(icri,isite,:,:,:) &
                                  +co2(1,3)*vprime1(icri,isite,:,:,:,3) &
                                  -co2(2,3)*vprime1(icri+1,isite,:,:,:,3)
        subo(icri+1,isite,:,:,:) = subo(icri+1,isite,:,:,:) &
                                  +co2(1,3)*vprime1(icri+1,isite,:,:,:,3)&
                                  +co2(2,3)*vprime1(icri,isite,:,:,:,3)
       enddo!isite
      enddo!icri

          
                                                  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!2order!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        do icri=1,5,2
         do isite=1,nvhalf
          sube(icri,isite,:,:,:) =  sube(icri,isite,:,:,:)&
                                   +co2(1,4)*vprime(icri,isite,:,:,:,4) &
                                   -co2(2,4)*vprime(icri+1,isite,:,:,:,4)
          sube(icri+1,isite,:,:,:) = sube(icri+1,isite,:,:,:)& 
                                    +co2(1,4)*vprime(icri+1,isite,:,:,:,4) &
                                    +co2(2,4)*vprime(icri,isite,:,:,:,4)

         enddo!isite
        enddo!icri

      do icri=1,5,2
       do isite=1,nvhalf
        subo(icri,isite,:,:,:) =  subo(icri,isite,:,:,:) &
                                  +co2(1,4)*vprime1(icri,isite,:,:,:,4)&
                                  -co2(2,4)*vprime1(icri+1,isite,:,:,:,4)
        subo(icri+1,isite,:,:,:) = subo(icri+1,isite,:,:,:) &
                                   +co2(1,4)*vprime1(icri+1,isite,:,:,:,4)&
                                   +co2(2,4)*vprime1(icri,isite,:,:,:,4)
       enddo!isite
      enddo!icri

              
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!3order!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      
        case(3) ! subtraction of O(kappa^4)
          
        do icri=1,5,2
         do isite=1,nvhalf
          sube(icri,isite,:,:,:) =  sube(icri,isite,:,:,:)&
                                   +co2(1,5)*vprime(icri,isite,:,:,:,5) &
                                   -co2(2,5)*vprime(icri+1,isite,:,:,:,5)
          sube(icri+1,isite,:,:,:) = sube(icri+1,isite,:,:,:)& 
                                    +co2(1,5)*vprime(icri+1,isite,:,:,:,5) &
                                    +co2(2,5)*vprime(icri,isite,:,:,:,5)

         enddo!isite
        enddo!icri

      do icri=1,5,2
       do isite=1,nvhalf
        subo(icri,isite,:,:,:) =  subo(icri,isite,:,:,:) &
                                  +co2(1,5)*vprime1(icri,isite,:,:,:,5) &
                                  -co2(2,5)*vprime1(icri+1,isite,:,:,:,5)
        subo(icri+1,isite,:,:,:) = subo(icri+1,isite,:,:,:) &
                                   +co2(1,5)*vprime1(icri+1,isite,:,:,:,5)&
                                   +co2(2,5)*vprime1(icri,isite,:,:,:,5)
       enddo!isite
      enddo!icri

            
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!4order!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        case(4) ! subtraction of O(kappa^5)
         do icri=1,5,2
         do isite=1,nvhalf
          sube(icri,isite,:,:,:) =  sube(icri,isite,:,:,:)&
                                   +co2(1,6)*vprime(icri,isite,:,:,:,6) &
                                   -co2(2,6)*vprime(icri+1,isite,:,:,:,6)
          sube(icri+1,isite,:,:,:) = sube(icri+1,isite,:,:,:)& 
                                    +co2(1,6)*vprime(icri+1,isite,:,:,:,6) &
                                    +co2(2,6)*vprime(icri,isite,:,:,:,6)

         enddo!isite
        enddo!icri

      do icri=1,5,2
       do isite=1,nvhalf
        subo(icri,isite,:,:,:) =  subo(icri,isite,:,:,:) &
                                  +co2(1,6)*vprime1(icri,isite,:,:,:,6) &
                                  -co2(2,6)*vprime1(icri+1,isite,:,:,:,6)
        subo(icri+1,isite,:,:,:) = subo(icri+1,isite,:,:,:) &
                                   +co2(1,6)*vprime1(icri+1,isite,:,:,:,6)&
                                   +co2(2,6)*vprime1(icri,isite,:,:,:,6)
       enddo!isite
      enddo!icri
        
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!5order!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
        case(5) ! subtraction of O(kappa^6)
        do icri=1,5,2
         do isite=1,nvhalf
          sube(icri,isite,:,:,:) =  sube(icri,isite,:,:,:)&
                                   +co2(1,7)*vprime(icri,isite,:,:,:,7) &
                                   -co2(2,7)*vprime(icri+1,isite,:,:,:,7)
          sube(icri+1,isite,:,:,:) = sube(icri+1,isite,:,:,:)& 
                                    +co2(1,7)*vprime(icri+1,isite,:,:,:,7) &
                                    +co2(2,7)*vprime(icri,isite,:,:,:,7)

         enddo!isite
        enddo!icri

      do icri=1,5,2
       do isite=1,nvhalf
        subo(icri,isite,:,:,:) =  subo(icri,isite,:,:,:) &
                                  +co2(1,7)*vprime1(icri,isite,:,:,:,7) &
                                  -co2(2,7)*vprime1(icri+1,isite,:,:,:,7)
        subo(icri+1,isite,:,:,:) = subo(icri+1,isite,:,:,:) &
                                   +co2(1,7)*vprime1(icri+1,isite,:,:,:,7)&
                                   +co2(2,7)*vprime1(icri,isite,:,:,:,7)
       enddo!isite
      enddo!icri

    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!6order!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
        case(6) ! subtraction of O(kappa^7)
        do icri=1,5,2
         do isite=1,nvhalf
          sube(icri,isite,:,:,:) =  sube(icri,isite,:,:,:)&
                                   +co2(1,8)*vprime(icri,isite,:,:,:,8) &
                                   -co2(2,8)*vprime(icri+1,isite,:,:,:,8)
          sube(icri+1,isite,:,:,:) = sube(icri+1,isite,:,:,:)& 
                                    +co2(1,8)*vprime(icri+1,isite,:,:,:,8) &
                                    +co2(2,8)*vprime(icri,isite,:,:,:,8)

         enddo!isite
        enddo!icri

      do icri=1,5,2
       do isite=1,nvhalf
        subo(icri,isite,:,:,:) =  subo(icri,isite,:,:,:) &
                                  +co2(1,8)*vprime1(icri,isite,:,:,:,8) &
                                  -co2(2,8)*vprime1(icri+1,isite,:,:,:,8)
        subo(icri+1,isite,:,:,:) = subo(icri+1,isite,:,:,:) &
                                   +co2(1,8)*vprime1(icri+1,isite,:,:,:,8)&
                                   +co2(2,8)*vprime1(icri,isite,:,:,:,8)
       enddo!isite
      enddo!icri

!         do icri=1,5,2
!         do isite=1,nvhalf
!          sube(icri,isite,:,:,:) =  sube(icri,isite,:,:,:)&
!                                   +co(1,9)*vprime(icri,isite,:,:,:,9) &
!                                   -co(2,9)*vprime(icri+1,isite,:,:,:,9)
!          sube(icri+1,isite,:,:,:) = sube(icri+1,isite,:,:,:)& 
!                                    +co(1,9)*vprime(icri+1,isite,:,:,:,9) &
!                                    +co(2,9)*vprime(icri,isite,:,:,:,9)

!         enddo!isite
!        enddo!icri

!      do icri=1,5,2
!       do isite=1,nvhalf
!        subo(icri,isite,:,:,:) =  subo(icri,isite,:,:,:) &
!                                  +co(1,9)*vprime1(icri,isite,:,:,:,9) &
!                                  -co(2,9)*vprime1(icri+1,isite,:,:,:,9)
!        subo(icri+1,isite,:,:,:) = subo(icri+1,isite,:,:,:) &
!                                   +co(1,9)*vprime1(icri+1,isite,:,:,:,9)&
!                                   +co(2,9)*vprime1(icri,isite,:,:,:,9)
!       enddo!isite
!      enddo!icri
        
!         do icri=1,5,2
!         do isite=1,nvhalf
!          sube(icri,isite,:,:,:) =  sube(icri,isite,:,:,:)&
!                                   +co(1,10)*vprime(icri,isite,:,:,:,10) &
!                                   -co(2,10)*vprime(icri+1,isite,:,:,:,10)
!          sube(icri+1,isite,:,:,:) = sube(icri+1,isite,:,:,:)& 
!                                    +co(1,10)*vprime(icri+1,isite,:,:,:,10) &
!                                    +co(2,10)*vprime(icri,isite,:,:,:,10)

!         enddo!isite
!        enddo!icri

!      do icri=1,5,2
!       do isite=1,nvhalf
!        subo(icri,isite,:,:,:) =  subo(icri,isite,:,:,:) &
!                                  +co(1,10)*vprime1(icri,isite,:,:,:,10) &
!                                  -co(2,10)*vprime1(icri+1,isite,:,:,:,10)
!        subo(icri+1,isite,:,:,:) = subo(icri+1,isite,:,:,:) &
!                                   +co(1,10)*vprime1(icri+1,isite,:,:,:,10)&
!                                   +co(2,10)*vprime1(icri,isite,:,:,:,10)
!       enddo!isite
!      enddo!icri
                                       
!         do icri=1,5,2
!         do isite=1,nvhalf
!          sube(icri,isite,:,:,:) =  sube(icri,isite,:,:,:)&
!                                   +co(1,11)*vprime(icri,isite,:,:,:,11) &
!                                   -co(2,11)*vprime(icri+1,isite,:,:,:,11)
!          sube(icri+1,isite,:,:,:) = sube(icri+1,isite,:,:,:)& 
!                                    +co(1,11)*vprime(icri+1,isite,:,:,:,11) &
!                                    +co(2,11)*vprime(icri,isite,:,:,:,11)

!         enddo!isite
!        enddo!icri

!      do icri=1,5,2
!       do isite=1,nvhalf
!        subo(icri,isite,:,:,:) =  subo(icri,isite,:,:,:) &
!                                  +co(1,11)*vprime1(icri,isite,:,:,:,11) &
!                                  -co(2,11)*vprime1(icri+1,isite,:,:,:,11)
!        subo(icri+1,isite,:,:,:) = subo(icri+1,isite,:,:,:) &
!                                   +co(1,11)*vprime1(icri+1,isite,:,:,:,11)&
!                                   +co(2,11)*vprime1(icri,isite,:,:,:,11)
!       enddo!isite
!      enddo!icri
                                       

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!7order!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        case default
          open(unit=8,file="DISCONLOOPS.ERROR",action="write",status="replace", &
               form="formatted")
          write(unit=8,fmt=*) "subroutine average: isub =", isub
          close(unit=8,status="keep")
          stop
      end select
      if (isub==2) then
        do ibleo = 1,8
          do ieo = 1,2
            do id = 1,4
              do isite = 1,nvhalf
                do icri = 1,6
                  sub2e(icri,isite,id,ieo,ibleo) = xe(icri,isite,id,ieo,ibleo) &
                                               - sube(icri,isite,id,ieo,ibleo)
                  sub2o(icri,isite,id,ieo,ibleo) = xo(icri,isite,id,ieo,ibleo) &
                                               - subo(icri,isite,id,ieo,ibleo)
                enddo ! icri
              enddo ! isite
            enddo ! id
          enddo ! ieo
        enddo ! ibleo
        do ibleo = 1,8
          do ieo = 1,2
            call loopops(ksub,sub2e,sub2o,u,be,bo,Je,Jo,ieo,ibleo,bc,vecbl, &
                         vecblinv,myid,nn,ldiv,nms,lvbc,ib,lbd,iblv,MRT)
            call spacesum(Jtime(:,:,isub,:,:),Je,Jo,ieo,ibleo,momfac,nmom,nop, &
                          myid,vecbl,vecblinv,MRT)
          enddo ! ieo
        enddo ! ibleo
      else
        do ibleo = 1,8
          do ieo = 1,2
            do id = 1,4
              do isite = 1,nvhalf
                do icri = 1,6
                  sub1e(icri,isite,id,ieo,ibleo) = xe(icri,isite,id,ieo,ibleo) &
                                               - sube(icri,isite,id,ieo,ibleo)
                  sub1o(icri,isite,id,ieo,ibleo) = xo(icri,isite,id,ieo,ibleo) &
                                               - subo(icri,isite,id,ieo,ibleo)
                enddo ! icri
              enddo ! isite
            enddo ! id
          enddo ! ieo
        enddo ! ibleo
        do ibleo = 1,8
          do ieo = 1,2
            call loopops(ksub,sub1e,sub1o,u,be,bo,Je,Jo,ieo,ibleo,bc,vecbl, &
                         vecblinv,myid,nn,ldiv,nms,lvbc,ib,lbd,iblv,MRT)
            call spacesum(Jtime(:,:,isub,:,:),Je,Jo,ieo,ibleo,momfac,nmom,nop, &
                          myid,vecbl,vecblinv,MRT)
          enddo ! ieo
        enddo ! ibleo
      endif
    enddo ! isub

    ! Return temporal links at the maximal timestep to their true nonzero values.
    if (.false. .and. dobndry==1) then
      do ibl = 9,16
        iri = ibl - 8
        do isite = itbit+1,nvhalf
          u(:,isite,4,:,ibl) = ubndry(:,isite-itbit,:,iri)
        enddo ! isite
      enddo ! ibl
    endif
    
    ! A final normalization factor.
    Jtime(:,:,:,:,1) = kappa*Jtime(:,:,:,:,1)
    Jtime(:,:,:,:,2) = kappa*Jtime(:,:,:,:,2)
    Jtime(:,:,:,:,3) = kappa*Jtime(:,:,:,:,3)
    Jtime(:,:,:,:,4) = kappa*Jtime(:,:,:,:,4)
    end subroutine ppaverage1
!---------------------------------------------------------------------------
 subroutine ppaverage(Jtime,nsub,nmom,nop,momfac,u,be,bo,xe,xo,kappa,dobndry, &
                    coact,bc,vecbl,vecblinv,myid,nn,ldiv,nms,lvbc,ib,lbd, &
                    iblv,rwdir,MRT,MRT2)
!This is the 4th order subtraction.
! Purpose of this is simply to calculate averages of operators.
! Due to limitations in subroutine vev, local operators are corrected
! only to 5th order in kappa, and nonlocal ones to 4th order.
! INPUT:
!   nsub is the number of results, having different orders of subtraction,
!        to be considered (including the case of order=0).
!   nmom is the number of momenta to be considered (including momentum=0).
!   nop is the number of operators to be considered.
!   momfac(isite,ieo,ibl,imom,iop) are the momentum factors
!          where the first 3 entries are the lattice site
!          and the 4th entry is the momentum
!          and the last entry denotes one of the nop operators.
!   be() contains the source vector on even (gblclr=1) lattice sites.
!   bo() contains the source vector on odd (gblclr=2) lattice sites.
!        expected size: be(6,ntotal,4,2,8), bo(6,ntotal,4,2,8)
!                       where the first entry is real/imaginary and colour
!                       and the 3rd entry is the Dirac index
!                       and the other three entries give the lattice site.
!   xe() contains the sink vector on even (gblclr=1) lattice sites.
!   xo() contains the sink vector on odd (gblclr=2) lattice sites.
!        expected size: xe(6,ntotal,4,2,8), xo(6,nvhalf,4,2,8)
!                       where the first entry is real/imaginary and colour
!                       and the 3rd entry is the Dirac index
!                       and the other three entries give the lattice site.
! OUTPUT:
!   Jtime(iri,it,isub,imom,iop) are the averaged operators.
!        it=1..nt is the time slice on the full (multi-process) lattice.
!        isub=1..nsub, where nsub is the number of subtractions
!               considered, including "no subtractions".
!        imom=1..nmom, where nmom is the number of momenta considered.
!        For imom=1, iri=1..2 is the real/imaginary index and iop=1..nop.
!        For imom>1, iri=1..2 and iop=1..nop together label Im operators.

    real(kind=KR),    intent(out),   dimension(:,:,:,:,:) :: Jtime
    character(len=*), intent(in),    dimension(:)         :: rwdir
    real(kind=KR),    intent(in),    dimension(:,:,:,:,:) :: momfac
    integer(kind=KI), intent(in)    :: nsub, nmom, nop, dobndry, myid, MRT, MRT2
    real(kind=KR),    intent(in)                         :: kappa
    real(kind=KR),    intent(inout), dimension(:,:,:,:,:) :: u
    real(kind=KR),    intent(inout), dimension(:,:,:,:,:) :: be, bo
    real(kind=KR),    intent(inout), dimension(:,:,:,:,:) :: xe, xo
    real(kind=KR),    intent(in),    dimension(:,:,:)     :: coact
    integer(kind=KI), intent(in),    dimension(:)         :: bc, nms
    integer(kind=KI), intent(in),    dimension(:,:)       :: vecbl, vecblinv
    integer(kind=KI), intent(in),    dimension(:,:)       :: nn, iblv
    logical,          intent(in),    dimension(:)         :: ldiv
    integer(kind=KI), intent(in),    dimension(:,:,:)     :: lvbc
    integer(kind=KI), intent(in),    dimension(:,:,:,:)   :: ib
    logical,          intent(in),    dimension(:,:)       :: lbd

    integer(kind=KI) :: gblclr, icri, id, isub, ksub, iri, isite, ibl, itbit, &
                        ibleo, ieo, idag, p, i, j
    real(kind=KR)                               :: xk, fac1, fac2
    real(kind=KR),    dimension(6,ntotal,4,2,8) :: sub1e, sub1o, sub2e, sub2o,&
                                                   vecsrc, temp, tue, tuo, &
                                                   sub3e, sub3o
    real(kind=KR),    dimension(6,nvhalf,4,2,8) :: sube, subo
    integer(kind=KI), parameter                 :: itstep=2*nvhalf*npt/nt
    real(kind=KR),    dimension(18,itstep,2,8)  :: ubndry
    !real(kind=KR),    dimension(2,nvhalf,5)     :: Je, Jo !5 should be replaced by nop
    real(kind=KR),    dimension(2,nvhalf,nop)     :: Je, Jo !5 should be replaced by nop
    integer(kind=KI) :: ierr, it, iop
    real(kind=KR2),   dimension(2)                          :: beta, beta1
    integer(kind=KI) :: icolor, idirac, mu, ibl1, gblclr1, jbl, jbleo, jeo

    ! This local gaugelink is used for debugging
    real(kind=KR),  dimension(18,ntotal,4,2,16) :: uout

    ! Initialization.
    idag = 0
    Jtime = 0.0_KR
    p = 5
    vprime = 0.0_KR2
    vprime1 = 0.0_KR2
    ! This routine is used for debugging by setting logical to .true.
    !      in if statement below.
    if (.false.) then ! debug
      call printlog("Debugging average",myid,rwdir)
      call printlog("**CHANGE u's BACK**",myid,rwdir)
      call fakegauge(uout,myid,rwdir,MRT)
      u = uout
    endif ! debug

    ! Set temporal links at the maximal timestep (on the global lattice) to zero.
    if (.false. .and. dobndry==1) then
      itbit = nvhalf - itstep
      do ibl = 9,16
        iri = ibl - 8
        do isite = itbit+1,nvhalf
          ubndry(:,isite-itbit,:,iri) = u(:,isite,4,:,ibl)
          u(:,isite,4,:,ibl) = 0.0_KR
        enddo ! isite
      enddo ! ibl
    endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!define the right hand side for xe
!      gblclr = 1
!      idag = 0
!      call Hsingle(temp,u,bo,idag,coact,bc,gblclr,vecbl, &
!                   vecblinv,myid,nn,ldiv,nms,lvbc,ib,lbd,iblv,MRT)

!      fac1 = 1.0_KR/kappa**2
!      fac2 = 1.0_KR/kappa
!      do isite = 1,nvhalf
!       vecsrc(:,isite,:,:,:) = fac1*be(:,isite,:,:,:) + fac2*temp(:,isite,:,:,:)
!      enddo ! isite
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Determine the polynomial for M.!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    do isite = 1,nvhalf
     vprime(:,isite,:,:,:,1) = be(:,isite,:,:,:)
    enddo !isite

    do isite = 1,nvhalf
     vprime1(:,isite,:,:,:,1) = bo(:,isite,:,:,:)
    enddo !isite

!    do isite = 1,nvhalf
!     tue(:,isite,:,:,:) = be(:,isite,:,:,:)
!     tuo(:,isite,:,:,:) = bo(:,isite,:,:,:)
!    enddo !isite


    do i = 1,p
          gblclr = 1
          call Hsingle(tue,u,vprime1(:,:,:,:,:,i),idag,coact,bc,gblclr,vecbl, &
                       vecblinv,myid,nn,ldiv,nms,lvbc,ib,lbd,iblv,MRT)
          gblclr = 2
           call Hsingle(tuo,u,vprime(:,:,:,:,:,i),idag,coact,bc,gblclr,vecbl, &
                       vecblinv,myid,nn,ldiv,nms,lvbc,ib,lbd,iblv,MRT)
!          do ibleo = 1,8
!            do ieo = 1,2
!              do id = 1,4
                do isite = 1,nvhalf
!                  do icri = 1,6
                   vprime(:,isite,:,:,:,i+1) = vprime(:,isite,:,:,:,i) &
                                               - kappa*tue(:,isite,:,:,:)
                   vprime1(:,isite,:,:,:,i+1) = vprime1(:,isite,:,:,:,i) &
                                                - kappa*tuo(:,isite,:,:,:)
!                  enddo ! icri
                enddo ! isite
!              enddo ! id
!            enddo ! ieo
!          enddo ! ibleo
 !      do isite = 1,nvhalf
!        vprime(:,isite,:,:,:,i+1) = tue(:,isite,:,:,:)
!       enddo !isite

!       do isite = 1,nvhalf
!        vprime1(:,isite,:,:,:,i+1) = tuo(:,isite,:,:,:)
!       enddo !isite
  enddo !i

  if(psignal==0) then!determent whether to calculate coeff for different rhs  
    do i=2,p+1
     do j=2,p+1
      call vecdot(vprime(:,:,:,:,:,i),vprime(:,:,:,:,:,j),beta,MRT2)
      call vecdot(vprime1(:,:,:,:,:,i),vprime1(:,:,:,:,:,j),beta1,MRT2)
      lsmat2(1,i-1,j-1) = beta(1) + beta1(1)  !lsmat(2,p,p) ,cls(2,p,1)
      lsmat2(2,i-1,j-1) = beta(2) + beta1(2)
!      if(myid==0) then
!        print *, "i,j, lsmat(:,i,j)=", i-1,j-1, lsmat(:,i-1,j-1)
!      endif
     enddo!j
    enddo!i
        



   do i=2,p+1
     call vecdot(vprime(:,:,:,:,:,i),vprime(:,:,:,:,:,1),beta,MRT2)
     call vecdot(vprime1(:,:,:,:,:,i),vprime1(:,:,:,:,:,1),beta1,MRT2)
     cls2(1,i-1,1) = beta(1) + beta1(1)
     cls2(2,i-1,1) = beta(2) + beta1(2)
!     if(myid==0) then
!       print *, "i,cls(:,i)=", i-1, cls(:,i-1,1)
!     endif
   enddo!i
    
    call linearsolver(p,1,lsmat2,ipiv4,cls2)
    co4(:,:) = cls2(:,:,1)    
 endif!psignal
!    co = 0.0_KR2    
!    co(1,1) = 4
  !BS if(myid==0) then
   ! do i=1,p
   !  print *, "i,result(:,i)=", i, co4(:,i)
   ! enddo!i  
   !endif!myid
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! Compute the operators for each level of subtraction.
 do isub = 1,nsub
   select case(isub)
     case(1) ! no subtraction
       ksub = 0
       sube = 0.0_KR
       subo = 0.0_KR
     case(2) ! subtraction of order1,order2 and order3
       ksub = 1
        do icri=1,5,2
         do isite=1,nvhalf
          sube(icri,isite,:,:,:) =  sube(icri,isite,:,:,:) &
                                    +co4(1,1)*vprime(icri,isite,:,:,:,1) &
                                    -co4(2,1)*vprime(icri+1,isite,:,:,:,1)
          sube(icri+1,isite,:,:,:) = sube(icri+1,isite,:,:,:) & 
                                     +co4(1,1)*vprime(icri+1,isite,:,:,:,1) &
                                     +co4(2,1)*vprime(icri,isite,:,:,:,1)

         enddo!isite
        enddo!icri

      do icri=1,5,2
       do isite=1,nvhalf
        subo(icri,isite,:,:,:) =  subo(icri,isite,:,:,:) &
                                  +co4(1,1)*vprime1(icri,isite,:,:,:,1) &
                                  -co4(2,1)*vprime1(icri+1,isite,:,:,:,1)
       subo(icri+1,isite,:,:,:) = subo(icri+1,isite,:,:,:) &
                                  +co4(1,1)*vprime1(icri+1,isite,:,:,:,1) &
                                  +co4(2,1)*vprime1(icri,isite,:,:,:,1)
       enddo!isite
      enddo!icri
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!order0!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        do icri=1,5,2
         do isite=1,nvhalf
          sube(icri,isite,:,:,:) =  sube(icri,isite,:,:,:) &
                                    +co4(1,2)*vprime(icri,isite,:,:,:,2) &
                                    -co4(2,2)*vprime(icri+1,isite,:,:,:,2)
          sube(icri+1,isite,:,:,:) = sube(icri+1,isite,:,:,:) & 
                                    +co4(1,2)*vprime(icri+1,isite,:,:,:,2) &
                                    +co4(2,2)*vprime(icri,isite,:,:,:,2)

         enddo!isite
        enddo!icri

      do icri=1,5,2
       do isite=1,nvhalf
        subo(icri,isite,:,:,:) =  subo(icri,isite,:,:,:) &
                                  +co4(1,2)*vprime1(icri,isite,:,:,:,2) &
                                  -co4(2,2)*vprime1(icri+1,isite,:,:,:,2)
        subo(icri+1,isite,:,:,:) = subo(icri+1,isite,:,:,:) &
                                   +co4(1,2)*vprime1(icri+1,isite,:,:,:,2)&
                                   +co4(2,2)*vprime1(icri,isite,:,:,:,2)
       enddo!isite
      enddo!icri
       


           
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1order!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        do icri=1,5,2
         do isite=1,nvhalf
          sube(icri,isite,:,:,:) =  sube(icri,isite,:,:,:) &
                                    +co4(1,3)*vprime(icri,isite,:,:,:,3) &
                                    -co4(2,3)*vprime(icri+1,isite,:,:,:,3)
          sube(icri+1,isite,:,:,:) = sube(icri+1,isite,:,:,:) & 
                                     +co4(1,3)*vprime(icri+1,isite,:,:,:,3) &
                                     +co4(2,3)*vprime(icri,isite,:,:,:,3)

         enddo!isite
        enddo!icri

      do icri=1,5,2
       do isite=1,nvhalf
        subo(icri,isite,:,:,:) =  subo(icri,isite,:,:,:) &
                                  +co4(1,3)*vprime1(icri,isite,:,:,:,3) &
                                  -co4(2,3)*vprime1(icri+1,isite,:,:,:,3)
        subo(icri+1,isite,:,:,:) = subo(icri+1,isite,:,:,:) &
                                  +co4(1,3)*vprime1(icri+1,isite,:,:,:,3)&
                                  +co4(2,3)*vprime1(icri,isite,:,:,:,3)
       enddo!isite
      enddo!icri

          
                                                  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!2order!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        do icri=1,5,2
         do isite=1,nvhalf
          sube(icri,isite,:,:,:) =  sube(icri,isite,:,:,:)&
                                   +co4(1,4)*vprime(icri,isite,:,:,:,4) &
                                   -co4(2,4)*vprime(icri+1,isite,:,:,:,4)
          sube(icri+1,isite,:,:,:) = sube(icri+1,isite,:,:,:)& 
                                    +co4(1,4)*vprime(icri+1,isite,:,:,:,4) &
                                    +co4(2,4)*vprime(icri,isite,:,:,:,4)

         enddo!isite
        enddo!icri

      do icri=1,5,2
       do isite=1,nvhalf
        subo(icri,isite,:,:,:) =  subo(icri,isite,:,:,:) &
                                  +co4(1,4)*vprime1(icri,isite,:,:,:,4)&
                                  -co4(2,4)*vprime1(icri+1,isite,:,:,:,4)
        subo(icri+1,isite,:,:,:) = subo(icri+1,isite,:,:,:) &
                                   +co4(1,4)*vprime1(icri+1,isite,:,:,:,4)&
                                   +co4(2,4)*vprime1(icri,isite,:,:,:,4)
       enddo!isite
      enddo!icri

              
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!3order!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      
        case(3) ! subtraction of O(kappa^4)
          
        do icri=1,5,2
         do isite=1,nvhalf
          sube(icri,isite,:,:,:) =  sube(icri,isite,:,:,:)&
                                   +co4(1,5)*vprime(icri,isite,:,:,:,5) &
                                   -co4(2,5)*vprime(icri+1,isite,:,:,:,5)
          sube(icri+1,isite,:,:,:) = sube(icri+1,isite,:,:,:)& 
                                    +co4(1,5)*vprime(icri+1,isite,:,:,:,5) &
                                    +co4(2,5)*vprime(icri,isite,:,:,:,5)

         enddo!isite
        enddo!icri

      do icri=1,5,2
       do isite=1,nvhalf
        subo(icri,isite,:,:,:) =  subo(icri,isite,:,:,:) &
                                  +co4(1,5)*vprime1(icri,isite,:,:,:,5) &
                                  -co4(2,5)*vprime1(icri+1,isite,:,:,:,5)
        subo(icri+1,isite,:,:,:) = subo(icri+1,isite,:,:,:) &
                                   +co4(1,5)*vprime1(icri+1,isite,:,:,:,5)&
                                   +co4(2,5)*vprime1(icri,isite,:,:,:,5)
       enddo!isite
      enddo!icri

            
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!4order!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        case(4) ! subtraction of O(kappa^5)
         do icri=1,5,2
         do isite=1,nvhalf
          sube(icri,isite,:,:,:) =  sube(icri,isite,:,:,:)&
                                   +co4(1,5)*vprime(icri,isite,:,:,:,6) &
                                   -co4(2,5)*vprime(icri+1,isite,:,:,:,6)
          sube(icri+1,isite,:,:,:) = sube(icri+1,isite,:,:,:)& 
                                    +co4(1,5)*vprime(icri+1,isite,:,:,:,6) &
                                    +co4(2,5)*vprime(icri,isite,:,:,:,6)

         enddo!isite
        enddo!icri

      do icri=1,5,2
       do isite=1,nvhalf
        subo(icri,isite,:,:,:) =  subo(icri,isite,:,:,:) &
                                  +co4(1,5)*vprime1(icri,isite,:,:,:,6) &
                                  -co4(2,5)*vprime1(icri+1,isite,:,:,:,6)
        subo(icri+1,isite,:,:,:) = subo(icri+1,isite,:,:,:) &
                                   +co4(1,5)*vprime1(icri+1,isite,:,:,:,6)&
                                   +co4(2,5)*vprime1(icri,isite,:,:,:,6)
       enddo!isite
      enddo!icri
        
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!5order!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
        case(5) ! subtraction of O(kappa^6)
        do icri=1,5,2
         do isite=1,nvhalf
          sube(icri,isite,:,:,:) =  sube(icri,isite,:,:,:)&
                                   +co4(1,5)*vprime(icri,isite,:,:,:,7) &
                                   -co4(2,5)*vprime(icri+1,isite,:,:,:,7)
          sube(icri+1,isite,:,:,:) = sube(icri+1,isite,:,:,:)& 
                                    +co4(1,5)*vprime(icri+1,isite,:,:,:,7) &
                                    +co4(2,5)*vprime(icri,isite,:,:,:,7)

         enddo!isite
        enddo!icri

      do icri=1,5,2
       do isite=1,nvhalf
        subo(icri,isite,:,:,:) =  subo(icri,isite,:,:,:) &
                                  +co4(1,5)*vprime1(icri,isite,:,:,:,7) &
                                  -co4(2,5)*vprime1(icri+1,isite,:,:,:,7)
        subo(icri+1,isite,:,:,:) = subo(icri+1,isite,:,:,:) &
                                   +co4(1,5)*vprime1(icri+1,isite,:,:,:,7)&
                                   +co4(2,5)*vprime1(icri,isite,:,:,:,7)
       enddo!isite
      enddo!icri

    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!6order!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
        case(6) ! subtraction of O(kappa^7)
        do icri=1,5,2
         do isite=1,nvhalf
          sube(icri,isite,:,:,:) =  sube(icri,isite,:,:,:)&
                                   +co4(1,5)*vprime(icri,isite,:,:,:,8) &
                                   -co4(2,5)*vprime(icri+1,isite,:,:,:,8)
          sube(icri+1,isite,:,:,:) = sube(icri+1,isite,:,:,:)& 
                                    +co4(1,5)*vprime(icri+1,isite,:,:,:,8) &
                                    +co4(2,5)*vprime(icri,isite,:,:,:,8)

         enddo!isite
        enddo!icri

      do icri=1,5,2
       do isite=1,nvhalf
        subo(icri,isite,:,:,:) =  subo(icri,isite,:,:,:) &
                                  +co4(1,5)*vprime1(icri,isite,:,:,:,8) &
                                  -co4(2,5)*vprime1(icri+1,isite,:,:,:,8)
        subo(icri+1,isite,:,:,:) = subo(icri+1,isite,:,:,:) &
                                   +co4(1,5)*vprime1(icri+1,isite,:,:,:,8)&
                                   +co4(2,5)*vprime1(icri,isite,:,:,:,8)
       enddo!isite
      enddo!icri

!         do icri=1,5,2
!         do isite=1,nvhalf
!          sube(icri,isite,:,:,:) =  sube(icri,isite,:,:,:)&
!                                   +co(1,9)*vprime(icri,isite,:,:,:,9) &
!                                   -co(2,9)*vprime(icri+1,isite,:,:,:,9)
!          sube(icri+1,isite,:,:,:) = sube(icri+1,isite,:,:,:)& 
!                                    +co(1,9)*vprime(icri+1,isite,:,:,:,9) &
!                                    +co(2,9)*vprime(icri,isite,:,:,:,9)

!         enddo!isite
!        enddo!icri

!      do icri=1,5,2
!       do isite=1,nvhalf
!        subo(icri,isite,:,:,:) =  subo(icri,isite,:,:,:) &
!                                  +co(1,9)*vprime1(icri,isite,:,:,:,9) &
!                                  -co(2,9)*vprime1(icri+1,isite,:,:,:,9)
!        subo(icri+1,isite,:,:,:) = subo(icri+1,isite,:,:,:) &
!                                   +co(1,9)*vprime1(icri+1,isite,:,:,:,9)&
!                                   +co(2,9)*vprime1(icri,isite,:,:,:,9)
!       enddo!isite
!      enddo!icri
        
!         do icri=1,5,2
!         do isite=1,nvhalf
!          sube(icri,isite,:,:,:) =  sube(icri,isite,:,:,:)&
!                                   +co(1,10)*vprime(icri,isite,:,:,:,10) &
!                                   -co(2,10)*vprime(icri+1,isite,:,:,:,10)
!          sube(icri+1,isite,:,:,:) = sube(icri+1,isite,:,:,:)& 
!                                    +co(1,10)*vprime(icri+1,isite,:,:,:,10) &
!                                    +co(2,10)*vprime(icri,isite,:,:,:,10)

!         enddo!isite
!        enddo!icri

!      do icri=1,5,2
!       do isite=1,nvhalf
!        subo(icri,isite,:,:,:) =  subo(icri,isite,:,:,:) &
!                                  +co(1,10)*vprime1(icri,isite,:,:,:,10) &
!                                  -co(2,10)*vprime1(icri+1,isite,:,:,:,10)
!        subo(icri+1,isite,:,:,:) = subo(icri+1,isite,:,:,:) &
!                                   +co(1,10)*vprime1(icri+1,isite,:,:,:,10)&
!                                   +co(2,10)*vprime1(icri,isite,:,:,:,10)
!       enddo!isite
!      enddo!icri
                                       
!         do icri=1,5,2
!         do isite=1,nvhalf
!          sube(icri,isite,:,:,:) =  sube(icri,isite,:,:,:)&
!                                   +co(1,11)*vprime(icri,isite,:,:,:,11) &
!                                   -co(2,11)*vprime(icri+1,isite,:,:,:,11)
!          sube(icri+1,isite,:,:,:) = sube(icri+1,isite,:,:,:)& 
!                                    +co(1,11)*vprime(icri+1,isite,:,:,:,11) &
!                                    +co(2,11)*vprime(icri,isite,:,:,:,11)

!         enddo!isite
!        enddo!icri

!      do icri=1,5,2
!       do isite=1,nvhalf
!        subo(icri,isite,:,:,:) =  subo(icri,isite,:,:,:) &
!                                  +co(1,11)*vprime1(icri,isite,:,:,:,11) &
!                                  -co(2,11)*vprime1(icri+1,isite,:,:,:,11)
!        subo(icri+1,isite,:,:,:) = subo(icri+1,isite,:,:,:) &
!                                   +co(1,11)*vprime1(icri+1,isite,:,:,:,11)&
!                                   +co(2,11)*vprime1(icri,isite,:,:,:,11)
!       enddo!isite
!      enddo!icri
                                       

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!7order!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        case default
          open(unit=8,file="DISCONLOOPS.ERROR",action="write",status="replace", &
               form="formatted")
          write(unit=8,fmt=*) "subroutine average: isub =", isub
          close(unit=8,status="keep")
          stop
      end select
      if (isub==2) then
        do ibleo = 1,8
          do ieo = 1,2
            do id = 1,4
              do isite = 1,nvhalf
                do icri = 1,6
                  sub2e(icri,isite,id,ieo,ibleo) = xe(icri,isite,id,ieo,ibleo) &
                                               - sube(icri,isite,id,ieo,ibleo)
                  sub2o(icri,isite,id,ieo,ibleo) = xo(icri,isite,id,ieo,ibleo) &
                                               - subo(icri,isite,id,ieo,ibleo)
                enddo ! icri
              enddo ! isite
            enddo ! id
          enddo ! ieo
        enddo ! ibleo
        do ibleo = 1,8
          do ieo = 1,2
            call loopops(ksub,sub2e,sub2o,u,be,bo,Je,Jo,ieo,ibleo,bc,vecbl, &
                         vecblinv,myid,nn,ldiv,nms,lvbc,ib,lbd,iblv,MRT)
            call spacesum(Jtime(:,:,isub,:,:),Je,Jo,ieo,ibleo,momfac,nmom,nop, &
                          myid,vecbl,vecblinv,MRT)
          enddo ! ieo
        enddo ! ibleo
      else
        do ibleo = 1,8
          do ieo = 1,2
            do id = 1,4
              do isite = 1,nvhalf
                do icri = 1,6
                  sub1e(icri,isite,id,ieo,ibleo) = xe(icri,isite,id,ieo,ibleo) &
                                               - sube(icri,isite,id,ieo,ibleo)
                  sub1o(icri,isite,id,ieo,ibleo) = xo(icri,isite,id,ieo,ibleo) &
                                               - subo(icri,isite,id,ieo,ibleo)
                enddo ! icri
              enddo ! isite
            enddo ! id
          enddo ! ieo
        enddo ! ibleo
        do ibleo = 1,8
          do ieo = 1,2
            call loopops(ksub,sub1e,sub1o,u,be,bo,Je,Jo,ieo,ibleo,bc,vecbl, &
                         vecblinv,myid,nn,ldiv,nms,lvbc,ib,lbd,iblv,MRT)
            call spacesum(Jtime(:,:,isub,:,:),Je,Jo,ieo,ibleo,momfac,nmom,nop, &
                          myid,vecbl,vecblinv,MRT)
          enddo ! ieo
        enddo ! ibleo
      endif
    enddo ! isub

    ! Return temporal links at the maximal timestep to their true nonzero values.
    if (.false. .and. dobndry==1) then
      do ibl = 9,16
        iri = ibl - 8
        do isite = itbit+1,nvhalf
          u(:,isite,4,:,ibl) = ubndry(:,isite-itbit,:,iri)
        enddo ! isite
      enddo ! ibl
    endif
    
    ! A final normalization factor.
    Jtime(:,:,:,:,1) = kappa*Jtime(:,:,:,:,1)
    Jtime(:,:,:,:,2) = kappa*Jtime(:,:,:,:,2)
    Jtime(:,:,:,:,3) = kappa*Jtime(:,:,:,:,3)
    Jtime(:,:,:,:,4) = kappa*Jtime(:,:,:,:,4)
    end subroutine ppaverage

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 subroutine averagePRIME(Jtime,nsub,nmom,nop,momfac,u,be,bo,xe,xo,kappa,dobndry, &
                    coact,bc,vecbl,vecblinv,myid,nn,ldiv,nms,lvbc,ib,lbd, &
                    iblv,rwdir,MRT)

! This subroutine performs the perturbative subtraction similar to the subroutine
! average. However, it uses eigloopops instead of loopops because testFUNC subroutine
! gives only the solution of M*x'=gamma5*z. 
! -AA 11/26/2012 

! Purpose of this is simply to calculate averages of operators.
! Due to limitations in subroutine vev, local operators are corrected
! only to 5th order in kappa, and nonlocal ones to 4th order.
! INPUT:
!   nsub is the number of results, having different orders of subtraction,
!        to be considered (including the case of order=0).
!   nmom is the number of momenta to be considered (including momentum=0).
!   nop is the number of operators to be considered.
!   momfac(isite,ieo,ibl,imom,iop) are the momentum factors
!          where the first 3 entries are the lattice site
!          and the 4th entry is the momentum
!          and the last entry denotes one of the nop operators.
!   be() contains the source vector on even (gblclr=1) lattice sites.
!   bo() contains the source vector on odd (gblclr=2) lattice sites.
!        expected size: be(6,ntotal,4,2,8), bo(6,ntotal,4,2,8)
!                       where the first entry is real/imaginary and colour
!                       and the 3rd entry is the Dirac index
!                       and the other three entries give the lattice site.
!   xe() contains the sink vector on even (gblclr=1) lattice sites.
!   xo() contains the sink vector on odd (gblclr=2) lattice sites.
!        expected size: xe(6,ntotal,4,2,8), xo(6,nvhalf,4,2,8)
!                       where the first entry is real/imaginary and colour
!                       and the 3rd entry is the Dirac index
!                       and the other three entries give the lattice site.
! OUTPUT:
!   Jtime(iri,it,isub,imom,iop) are the averaged operators.
!        it=1..nt is the time slice on the full (multi-process) lattice.
!        isub=1..nsub, where nsub is the number of subtractions
!               considered, including "no subtractions".
!        imom=1..nmom, where nmom is the number of momenta considered.
!        For imom=1, iri=1..2 is the real/imaginary index and iop=1..nop.
!        For imom>1, iri=1..2 and iop=1..nop together label Im operators.

    real(kind=KR),    intent(out),   dimension(:,:,:,:,:) :: Jtime
    character(len=*), intent(in),    dimension(:)         :: rwdir
    real(kind=KR),    intent(in),    dimension(:,:,:,:,:) :: momfac
    integer(kind=KI), intent(in)         :: nsub, nmom, nop, dobndry, myid, MRT
    real(kind=KR),    intent(in)                          :: kappa
    real(kind=KR),    intent(inout), dimension(:,:,:,:,:) :: u
    real(kind=KR),    intent(inout), dimension(:,:,:,:,:) :: be, bo
    real(kind=KR),    intent(inout), dimension(:,:,:,:,:) :: xe, xo
    real(kind=KR),    intent(in),    dimension(:,:,:)     :: coact
    integer(kind=KI), intent(in),    dimension(:)         :: bc, nms
    integer(kind=KI), intent(in),    dimension(:,:)       :: vecbl, vecblinv
    integer(kind=KI), intent(in),    dimension(:,:)       :: nn, iblv
    logical,          intent(in),    dimension(:)         :: ldiv
    integer(kind=KI), intent(in),    dimension(:,:,:)     :: lvbc
    integer(kind=KI), intent(in),    dimension(:,:,:,:)   :: ib
    logical,          intent(in),    dimension(:,:)       :: lbd

    integer(kind=KI) :: gblclr, icri, id, isub, ksub, iri, isite, ibl, itbit, &
                        ibleo, ieo, idag
    real(kind=KR)                               :: xk
    real(kind=KR),    dimension(6,ntotal,4,2,8) :: sub1e, sub1o, sub2e, sub2o
    real(kind=KR),    dimension(6,nvhalf,4,2,8) :: sube, subo
    integer(kind=KI), parameter                 :: itstep=2*nvhalf*npt/nt
    real(kind=KR),    dimension(18,itstep,2,8)  :: ubndry
    !real(kind=KR),    dimension(2,nvhalf,5)     :: Je, Jo !5 should be replaced by nop
    real(kind=KR),    dimension(2,nvhalf,nop)     :: Je, Jo !5 should be replaced by nop
    integer(kind=KI) :: ierr, it, iop
    integer(kind=KI) :: icolor, idirac, mu, ibl1, gblclr1, jbl, jbleo, jeo

    ! This local gaugelink is used for debugging
    real(kind=KR),  dimension(18,ntotal,4,2,16) :: uout

    ! Initialization.
    idag = 0
    Jtime = 0.0_KR

    ! This routine is used for debugging by setting logical to .true.
    !      in if statement below.
    if (.false.) then ! debug
      call printlog("Debugging average",myid,rwdir)
      call printlog("**CHANGE u's BACK**",myid,rwdir)
      call fakegauge(uout,myid,rwdir,MRT)
      u = uout
    endif ! debug

    ! Set temporal links at the maximal timestep (on the global lattice) to zero.
    if (.false. .and. dobndry==1) then
      itbit = nvhalf - itstep
      do ibl = 9,16
        iri = ibl - 8
        do isite = itbit+1,nvhalf
          ubndry(:,isite-itbit,:,iri) = u(:,isite,4,:,ibl)
          u(:,isite,4,:,ibl) = 0.0_KR
        enddo ! isite
      enddo ! ibl
    endif

    ! Compute the operators for each level of subtraction.
    do isub = 1,nsub
      select case(isub)
        case(1) ! no subtraction
          ksub = 0
          sube = 0.0_KR
          subo = 0.0_KR
        case(2) ! subtraction of O(kappa), O(kappa^2) and O(kappa^3)
          ksub = 1
          gblclr = 1
          call Hsingle(sub1e,u,bo,idag,coact,bc,gblclr,vecbl,vecblinv,myid, &
                       nn,ldiv,nms,lvbc,ib,lbd,iblv,MRT)
          gblclr = 2
          call Hsingle(sub1o,u,be,idag,coact,bc,gblclr,vecbl,vecblinv,myid, &
                       nn,ldiv,nms,lvbc,ib,lbd,iblv,MRT)
          do ibleo = 1,8
            do ieo = 1,2
              do id = 1,4
                do isite = 1,nvhalf
                  do icri = 1,6
                    sube(icri,isite,id,ieo,ibleo) = be(icri,isite,id,ieo,ibleo) &
                                         + kappa*sub1e(icri,isite,id,ieo,ibleo)
                    subo(icri,isite,id,ieo,ibleo) = bo(icri,isite,id,ieo,ibleo) &
                                         + kappa*sub1o(icri,isite,id,ieo,ibleo)
                  enddo ! icri
                enddo ! isite
              enddo ! id
            enddo ! ieo
          enddo ! ibleo
          gblclr = 1
          call Hsingle(sub2e,u,sub1o,idag,coact,bc,gblclr,vecbl,vecblinv,myid, &
                       nn,ldiv,nms,lvbc,ib,lbd,iblv,MRT)
          gblclr = 2
          call Hsingle(sub2o,u,sub1e,idag,coact,bc,gblclr,vecbl,vecblinv,myid, &
                       nn,ldiv,nms,lvbc,ib,lbd,iblv,MRT)
          xk = kappa**2
          do ibleo = 1,8
            do ieo = 1,2
              do id = 1,4
                do isite = 1,nvhalf
                  do icri = 1,6
                    sube(icri,isite,id,ieo,ibleo) = sube(icri,isite,id,ieo,ibleo) &
                                              + xk*sub2e(icri,isite,id,ieo,ibleo)
                    subo(icri,isite,id,ieo,ibleo) = subo(icri,isite,id,ieo,ibleo) &
                                              + xk*sub2o(icri,isite,id,ieo,ibleo)
                  enddo ! icri
                enddo ! isite
              enddo ! id
            enddo ! ieo
          enddo ! ibleo
          gblclr = 1
          call Hsingle(sub1e,u,sub2o,idag,coact,bc,gblclr,vecbl,vecblinv,myid, &
                       nn,ldiv,nms,lvbc,ib,lbd,iblv,MRT)
          gblclr = 2
          call Hsingle(sub1o,u,sub2e,idag,coact,bc,gblclr,vecbl,vecblinv,myid, &
                       nn,ldiv,nms,lvbc,ib,lbd,iblv,MRT)
          xk = kappa**3
          do ibleo = 1,8
            do ieo = 1,2
              do id = 1,4
                do isite = 1,nvhalf
                  do icri = 1,6
                    sube(icri,isite,id,ieo,ibleo) = sube(icri,isite,id,ieo,ibleo) &
                                              + xk*sub1e(icri,isite,id,ieo,ibleo)
                    subo(icri,isite,id,ieo,ibleo) = subo(icri,isite,id,ieo,ibleo) &
                                              + xk*sub1o(icri,isite,id,ieo,ibleo)
                  enddo ! icri
                enddo ! isite
              enddo ! id
            enddo ! ieo
          enddo ! ibleo
        case(3) ! subtraction of O(kappa^4)
          gblclr = 1
          call Hsingle(sub2e,u,sub1o,idag,coact,bc,gblclr,vecbl,vecblinv,myid, &
                       nn,ldiv,nms,lvbc,ib,lbd,iblv,MRT)
          gblclr = 2
          call Hsingle(sub2o,u,sub1e,idag,coact,bc,gblclr,vecbl,vecblinv,myid, &
                       nn,ldiv,nms,lvbc,ib,lbd,iblv,MRT)
          xk = kappa**4
          do ibleo = 1,8
            do ieo = 1,2
              do id = 1,4
                do isite = 1,nvhalf
                  do icri = 1,6
                    sube(icri,isite,id,ieo,ibleo) = sube(icri,isite,id,ieo,ibleo) &
                                              + xk*sub2e(icri,isite,id,ieo,ibleo)
                    subo(icri,isite,id,ieo,ibleo) = subo(icri,isite,id,ieo,ibleo) &
                                              + xk*sub2o(icri,isite,id,ieo,ibleo)
                  enddo ! icri
                enddo ! isite
              enddo ! id
            enddo ! ieo
          enddo ! ibleo
        case(4) ! subtraction of O(kappa^5)
          gblclr = 1
          call Hsingle(sub1e,u,sub2o,idag,coact,bc,gblclr,vecbl,vecblinv,myid, &
                       nn,ldiv,nms,lvbc,ib,lbd,iblv,MRT)
          gblclr = 2
          call Hsingle(sub1o,u,sub2e,idag,coact,bc,gblclr,vecbl,vecblinv,myid, &
                       nn,ldiv,nms,lvbc,ib,lbd,iblv,MRT)
          xk = kappa**5
          do ibleo = 1,8
            do ieo = 1,2
              do id = 1,4
                do isite = 1,nvhalf
                  do icri = 1,6
                    sube(icri,isite,id,ieo,ibleo) = sube(icri,isite,id,ieo,ibleo) &
                                              + xk*sub1e(icri,isite,id,ieo,ibleo)
                    subo(icri,isite,id,ieo,ibleo) = subo(icri,isite,id,ieo,ibleo) &
                                              + xk*sub1o(icri,isite,id,ieo,ibleo)
                  enddo ! icri
                enddo ! isite
              enddo ! id
            enddo ! ieo
          enddo ! ibleo
        case(5) ! subtraction of O(kappa^6)
          gblclr = 1
          call Hsingle(sub1e,u,sub2o,idag,coact,bc,gblclr,vecbl,vecblinv,myid, &
                       nn,ldiv,nms,lvbc,ib,lbd,iblv,MRT)
          gblclr = 2
          call Hsingle(sub1o,u,sub2e,idag,coact,bc,gblclr,vecbl,vecblinv,myid, &
                       nn,ldiv,nms,lvbc,ib,lbd,iblv,MRT)
          xk = kappa**6
          do ibleo = 1,8
            do ieo = 1,2
              do id = 1,4
                do isite = 1,nvhalf
                  do icri = 1,6
                    sube(icri,isite,id,ieo,ibleo) = sube(icri,isite,id,ieo,ibleo) &
                                              + xk*sub1e(icri,isite,id,ieo,ibleo)
                    subo(icri,isite,id,ieo,ibleo) = subo(icri,isite,id,ieo,ibleo) &
                                              + xk*sub1o(icri,isite,id,ieo,ibleo)
                  enddo ! icri
                enddo ! isite
              enddo ! id
            enddo ! ieo
          enddo ! ibleo
        case(6) ! subtraction of O(kappa^7)
          gblclr = 1
          call Hsingle(sub1e,u,sub2o,idag,coact,bc,gblclr,vecbl,vecblinv,myid, &
                       nn,ldiv,nms,lvbc,ib,lbd,iblv,MRT)
          gblclr = 2
          call Hsingle(sub1o,u,sub2e,idag,coact,bc,gblclr,vecbl,vecblinv,myid, &
                       nn,ldiv,nms,lvbc,ib,lbd,iblv,MRT)
          xk = kappa**6
          do ibleo = 1,8
            do ieo = 1,2
              do id = 1,4
                do isite = 1,nvhalf
                  do icri = 1,6
                    sube(icri,isite,id,ieo,ibleo) = sube(icri,isite,id,ieo,ibleo) &
                                              + xk*sub1e(icri,isite,id,ieo,ibleo)
                    subo(icri,isite,id,ieo,ibleo) = subo(icri,isite,id,ieo,ibleo) &
                                              + xk*sub1o(icri,isite,id,ieo,ibleo)
                  enddo ! icri
                enddo ! isite
              enddo ! id
            enddo ! ieo
          enddo ! ibleo
        case default
          open(unit=8,file="DISCONLOOPS.ERROR",action="write",status="replace", &
               form="formatted")
          write(unit=8,fmt=*) "subroutine average: isub =", isub
          close(unit=8,status="keep")
          stop
      end select
      if (isub==2) then
        do ibleo = 1,8
          do ieo = 1,2
            do id = 1,4
              do isite = 1,nvhalf
                do icri = 1,6
                  sub2e(icri,isite,id,ieo,ibleo) = xe(icri,isite,id,ieo,ibleo) &
                                               - sube(icri,isite,id,ieo,ibleo)
                  sub2o(icri,isite,id,ieo,ibleo) = xo(icri,isite,id,ieo,ibleo) &
                                               - subo(icri,isite,id,ieo,ibleo)
                enddo ! icri
              enddo ! isite
            enddo ! id
          enddo ! ieo
        enddo ! ibleo
        do ibleo = 1,8
          do ieo = 1,2
            call eigloopops(ksub,sub2e,sub2o,u,be,bo,Je,Jo,ieo,ibleo,bc,vecbl, &
                         vecblinv,myid,nn,ldiv,nms,lvbc,ib,lbd,iblv,MRT)
            call spacesum(Jtime(:,:,isub,:,:),Je,Jo,ieo,ibleo,momfac,nmom,nop, &
                          myid,vecbl,vecblinv,MRT)
          enddo ! ieo
        enddo ! ibleo
      else
        do ibleo = 1,8
          do ieo = 1,2
            do id = 1,4
              do isite = 1,nvhalf
                do icri = 1,6
                  sub1e(icri,isite,id,ieo,ibleo) = xe(icri,isite,id,ieo,ibleo) &
                                               - sube(icri,isite,id,ieo,ibleo)
                  sub1o(icri,isite,id,ieo,ibleo) = xo(icri,isite,id,ieo,ibleo) &
                                               - subo(icri,isite,id,ieo,ibleo)
                enddo ! icri
              enddo ! isite
            enddo ! id
          enddo ! ieo
        enddo ! ibleo
        do ibleo = 1,8
          do ieo = 1,2
            call eigloopops(ksub,sub1e,sub1o,u,be,bo,Je,Jo,ieo,ibleo,bc,vecbl, &
                         vecblinv,myid,nn,ldiv,nms,lvbc,ib,lbd,iblv,MRT)
            call spacesum(Jtime(:,:,isub,:,:),Je,Jo,ieo,ibleo,momfac,nmom,nop, &
                          myid,vecbl,vecblinv,MRT)
          enddo ! ieo
        enddo ! ibleo
      endif
    enddo ! isub

    ! Return temporal links at the maximal timestep to their true nonzero values.
    if (.false. .and. dobndry==1) then
      do ibl = 9,16
        iri = ibl - 8
        do isite = itbit+1,nvhalf
          u(:,isite,4,:,ibl) = ubndry(:,isite-itbit,:,iri)
        enddo ! isite
      enddo ! ibl
    endif
    
    ! A final normalization factor.
    Jtime(:,:,:,:,1) = kappa*Jtime(:,:,:,:,1)
    Jtime(:,:,:,:,2) = kappa*Jtime(:,:,:,:,2)
    Jtime(:,:,:,:,3) = kappa*Jtime(:,:,:,:,3)
    Jtime(:,:,:,:,4) = kappa*Jtime(:,:,:,:,4)
    end subroutine averagePRIME

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 subroutine loopops(ksub,xe,xo,u,be,bo,Je,Jo,ieo,ibleo,bc,vecbl,vecblinv, &
                    myid,nn,ldiv,nms,lvbc,ib,lbd,iblv,MRT)
! Construct the operators to be included in disconnected loop contributions.

! INPUT:
!   ksub = 0 if the propagator has NOT been subtracted.
!   ksub = 1 if the propagator HAS been subtracted.
!   be() contains the source vector on even (gblclr=1) lattice sites.
!   bo() contains the source vector on odd (gblclr=2) lattice sites.
!        expected size: be(6,ntotal,4,2,8), bo(6,ntotal,4,2,8)
!                       where the first entry is real/imaginary and colour
!                       and the 3rd entry is the Dirac index
!                       and the other three entries give the lattice site.
!   xe() contains the sink vector on even (gblclr=1) lattice sites.
!   xo() contains the sink vector on odd (gblclr=2) lattice sites.
!        expected size: xe(6,ntotal,4,2,8), xo(6,ntotal,4,2,8)
!                       where the first entry is real/imaginary and colour
!                       and the 3rd entry is the Dirac index
!                       and the other three entries give the lattice site.
! OUTPUT:
!   Je() contains the operators on even (gblclr=1) lattice sites.
!   Jo() contains the operators on odd (gblclr=2) lattice sites.
!        expected size: Je(2,nvhalf,5), Jo(2,nvhalf,5)
!                       where the first entry is 1=real or 2=imaginary
!                       and the middle entry gives the lattice site
!                       and the last entry denotes one of the five operators.

    integer(kind=KI), intent(in)                 :: ksub, ieo, ibleo, myid, MRT
    real(kind=KR),    intent(inout), dimension(:,:,:,:,:) :: xe, xo
    real(kind=KR),    intent(in),    dimension(:,:,:,:,:) :: u
    real(kind=KR),    intent(inout), dimension(:,:,:,:,:) :: be, bo
    real(kind=KR),    intent(out),   dimension(:,:,:)     :: Je, Jo
    integer(kind=KI), intent(in),    dimension(:)         :: bc, nms
    integer(kind=KI), intent(in),    dimension(:,:)       :: vecbl, vecblinv
    integer(kind=KI), intent(in),    dimension(:,:)       :: nn, iblv
    logical,          intent(in),    dimension(:)         :: ldiv
    integer(kind=KI), intent(in),    dimension(:,:,:)     :: lvbc
    integer(kind=KI), intent(in),    dimension(:,:,:,:)   :: ib
    logical,          intent(in),    dimension(:,:)       :: lbd

    integer(kind=KI)                     :: gblclr, mu, ierr
!* tempe, tempo new temporary arrays (added by WW).
    real(kind=KR), dimension(6,nvhalf,4) :: temp,tempe,tempo
    integer(kind=KI) :: id, isite, ibl
!*Initializations.
    Je = 0.0_KR
    Jo = 0.0_KR

!*Point-split vector operators:

   do mu = 1,4

 ! DEAN !~ I put this in
   if (.false.) then
     gblclr = 1
     ibl = vecbl(gblclr,ibleo)
       do isite=1,nvhalf
          do id=1,4
          if (be(1,isite,id,ieo,ibleo)/=0.0_KR) then
             print "(a16,1i3,2es17.10)","myid: u,be=", myid, u(1,isite,mu,ieo,ibl),be(1,isite,id,ieo,ibleo)
          endif
          enddo ! id
       enddo! isite
   endif ! .false.
 ! end DEAN

! Je_mu(x)
! = Je_mu(x) + bebar(x) (1-gamma_mu) U_mu(x) xo(x+mu)
! = Je_mu(x) + [ xo^dagger(x+mu) (1-gamma_mu) U_mu^dagger(x) bebar^dagger(x) ]^*
! (The complex conjugation will be done at the end of this subroutine.)
     gblclr = 1
     call mulfor2(temp,u,be,mu,gblclr,ieo,ibleo,vecbl)
     call vgv(xo,temp,Je(:,:,mu),gblclr,ieo,ibleo,mu,bc,vecbl,vecblinv, &
              lbd,ldiv,ib,lvbc,nms,myid,nn,iblv,MRT)
 
! Jo_mu(x)
! = Jo_mu(x) + bobar(x) (1-gamma_mu) U_mu(x) xe(x+mu)
! = Jo_mu(x) + [ xe^dagger(x+mu) (1-gamma_mu) U_mu^dagger(x) bobar^dagger(x) ]^*
! (The complex conjugation will be done at the end of this subroutine.)
     gblclr = 2
     call mulfor2(temp,u,bo,mu,gblclr,ieo,ibleo,vecbl)
     call vgv(xe,temp,Jo(:,:,mu),gblclr,ieo,ibleo,mu,bc,vecbl,vecblinv, &
              lbd,ldiv,ib,lvbc,nms,myid,nn,iblv,MRT)

! Je_mu(x)
! = Je_mu(x) - bobar(x+mu) (1+gamma_mu) U_mu^dagger(x) xe(x)
! = Je_mu(x) - [ xe^dagger(x) (1+gamma_mu) U_mu(x) bobar^dagger(x+mu) ]^*
! (The complex conjugation will be done at the end of this subroutine.)
     gblclr = 1
     call mulbac2(temp,u,bo,mu,bc,gblclr,ieo,ibleo,vecbl,vecblinv,myid,nn, &
                  ldiv,nms,lvbc,ib,lbd,iblv,MRT)
     call vv(xe(:,:,:,ieo,ibleo),temp,Je(:,:,mu))
 
! Jo_mu(x)
! = Jo_mu(x) - bebar(x+mu) (1+gamma_mu) U_mu^dagger(x) xo(x)
! = Jo_mu(x) - [ xo^dagger(x) (1+gamma_mu) U_mu(x) bebar^dagger(x+mu) ]^*
! (The complex conjugation will be done at the end of this subroutine.)
     gblclr = 2
     call mulbac2(temp,u,be,mu,bc,gblclr,ieo,ibleo,vecbl,vecblinv,myid,nn, &
                  ldiv,nms,lvbc,ib,lbd,iblv,MRT)
     call vv(xo(:,:,:,ieo,ibleo),temp,Jo(:,:,mu))
 
    enddo ! mu

!*Local scalar operator [ psibar(x) psi(x) ] :
    call vv(xe(:,:,:,ieo,ibleo),be(:,:,:,ieo,ibleo),Je(:,:,5))
    call vv(xo(:,:,:,ieo,ibleo),bo(:,:,:,ieo,ibleo),Jo(:,:,5))
    if (ksub==1) then
     Je(1,:,5) = Je(1,:,5) + 12.0_KR ! 12 = 3 colours times 4 Dirac indices.
     Jo(1,:,5) = Jo(1,:,5) + 12.0_KR ! 12 = 3 colours times 4 Dirac indices.
    endif

!* Local vector operators (added by WW)


!##########I un-commented out this part to include the local vector current --Abdou 10/22/12#####
  ! Wilcox modification - Victor commented out
  do mu = 1,4
!* First multiply be and bo by gamma_mu and then call vv. We need to define Je and Jo
!* which have components which go from 6 to 9. How do we do this, Vic??

    ! I believe this is correct (7/9/12) -- Victor
    !call gammamult(be(:,:,:,ieo,ibleo),bo(:,:,:,ieo,ibleo),tempe,tempo,mu,ieo,ibleo) !changed gammamult into gammamultshort -AA
    call gammamultshort(be(:,:,:,ieo,ibleo),tempe,mu) !changed gammamult into gammamultshort -AA
    call gammamultshort(bo(:,:,:,ieo,ibleo),tempo,mu) !changed gammamult into gammamultshort -AA
    call vv(xe(:,:,:,ieo,ibleo),tempe,Je(:,:,mu+5))
    call vv(xo(:,:,:,ieo,ibleo),tempo,Jo(:,:,mu+5))


   enddo ! mu
!################End un-commenting out by Abdou#################################################

!*Now do the complex conjugation.
    Je(2,:,:) = -Je(2,:,:)
    Jo(2,:,:) = -Jo(2,:,:)

!   print *, "myid,Je(1,:,1)=", myid,Je(1,:,1)

!   call MPI_BARRIER(MPI_COMM_WORLD, ierr)
!   stop

 end subroutine loopops




! THIS FUNCTION WORKS EXACTLY LIKE loopops() EXCEPT THAT THE SOURCE VECTOR
! IS EXPECTED TO HOLD THE SOLUTION OF gamma5Mx'=z NOT Mx=z FOR USE WITH THE
! HERMITIAN FORCED EIGENSPECTRUM SUBTRACTION TECHNIQUE
! NOTE: source vector is "b" and noise vector is "x"
 subroutine eigloopops_original(ksub,xe,xo,u,be,bo,Je,Jo,ieo,ibleo,bc,vecbl,vecblinv, &
                    myid,nn,ldiv,nms,lvbc,ib,lbd,iblv,MRT)
! Construct the operators to be included in disconnected loop contributions.

! INPUT:
!   ksub = 0 if the propagator has NOT been subtracted.
!   ksub = 1 if the propagator HAS been subtracted.
!   be() contains the source vector on even (gblclr=1) lattice sites.
!   bo() contains the source vector on odd (gblclr=2) lattice sites.
!        expected size: be(6,ntotal,4,2,8), bo(6,ntotal,4,2,8)
!                       where the first entry is real/imaginary and colour
!                       and the 3rd entry is the Dirac index
!                       and the other three entries give the lattice site.
!   xe() contains the sink vector on even (gblclr=1) lattice sites.
!   xo() contains the sink vector on odd (gblclr=2) lattice sites.
!        expected size: xe(6,ntotal,4,2,8), xo(6,ntotal,4,2,8)
!                       where the first entry is real/imaginary and colour
!                       and the 3rd entry is the Dirac index
!                       and the other three entries give the lattice site.
! OUTPUT:
!   Je() contains the operators on even (gblclr=1) lattice sites.
!   Jo() contains the operators on odd (gblclr=2) lattice sites.
!        expected size: Je(2,nvhalf,5), Jo(2,nvhalf,5)
!                       where the first entry is 1=real or 2=imaginary
!                       and the middle entry gives the lattice site
!                       and the last entry denotes one of the five operators.

    integer(kind=KI), intent(in)                 :: ksub, ieo, ibleo, myid, MRT
    real(kind=KR),    intent(inout), dimension(:,:,:,:,:) :: xe, xo
    real(kind=KR),    intent(in),    dimension(:,:,:,:,:) :: u
    real(kind=KR),    intent(inout), dimension(:,:,:,:,:) :: be, bo
    real(kind=KR),    intent(out),   dimension(:,:,:)     :: Je, Jo
    integer(kind=KI), intent(in),    dimension(:)         :: bc, nms
    integer(kind=KI), intent(in),    dimension(:,:)       :: vecbl, vecblinv
    integer(kind=KI), intent(in),    dimension(:,:)       :: nn, iblv
    logical,          intent(in),    dimension(:)         :: ldiv
    integer(kind=KI), intent(in),    dimension(:,:,:)     :: lvbc
    integer(kind=KI), intent(in),    dimension(:,:,:,:)   :: ib
    logical,          intent(in),    dimension(:,:)       :: lbd

    integer(kind=KI)                     :: gblclr, mu, ierr
    real(kind=KR), dimension(6,nvhalf,4) :: temp
    real(kind=KR), dimension(6,nvhalf,4) :: gtemp
    integer(kind=KI) :: id, isite, ibl
!*Initializations.
    Je = 0.0_KR
    Jo = 0.0_KR

!*Point-split vector operators:

   do mu = 1,4

! Je_mu(x)
! = Je_mu(x) + bebar(x) (1-gamma_mu) U_mu(x) xo(x+mu)
! = Je_mu(x) + [ xo^dagger(x+mu) (1-gamma_mu) U_mu^dagger(x) bebar^dagger(x) ]^*
! (The complex conjugation will be done at the end of this subroutine.)
     gblclr = 1
     call mulfor2(temp,u,be,mu,gblclr,ieo,ibleo,vecbl)
     call gammamultshort(temp, gtemp, 5) ! multiply by gamma5
     call vgv(xo,gtemp,Je(:,:,mu),gblclr,ieo,ibleo,mu,bc,vecbl,vecblinv, &
              lbd,ldiv,ib,lvbc,nms,myid,nn,iblv,MRT)
 
! Jo_mu(x)
! = Jo_mu(x) + bobar(x) (1-gamma_mu) U_mu(x) xe(x+mu)
! = Jo_mu(x) + [ xe^dagger(x+mu) (1-gamma_mu) U_mu^dagger(x) bobar^dagger(x) ]^*
! (The complex conjugation will be done at the end of this subroutine.)
     gblclr = 2
     call mulfor2(temp,u,bo,mu,gblclr,ieo,ibleo,vecbl)
     call gammamultshort(temp, gtemp, 5) ! multiply by gamma5
     call vgv(xe,gtemp,Jo(:,:,mu),gblclr,ieo,ibleo,mu,bc,vecbl,vecblinv, &
              lbd,ldiv,ib,lvbc,nms,myid,nn,iblv,MRT)

! Je_mu(x)
! = Je_mu(x) - bobar(x+mu) (1+gamma_mu) U_mu^dagger(x) xe(x)
! = Je_mu(x) - [ xe^dagger(x) (1+gamma_mu) U_mu(x) bobar^dagger(x+mu) ]^*
! (The complex conjugation will be done at the end of this subroutine.)
     gblclr = 1
     call mulbac2(temp,u,bo,mu,bc,gblclr,ieo,ibleo,vecbl,vecblinv,myid,nn, &
                  ldiv,nms,lvbc,ib,lbd,iblv,MRT)
     call gammamultshort(temp, gtemp, 5) ! multiply by gamma5
     call vv(xe(:,:,:,ieo,ibleo),gtemp,Je(:,:,mu))
 
! Jo_mu(x)
! = Jo_mu(x) - bebar(x+mu) (1+gamma_mu) U_mu^dagger(x) xo(x)
! = Jo_mu(x) - [ xo^dagger(x) (1+gamma_mu) U_mu(x) bebar^dagger(x+mu) ]^*
! (The complex conjugation will be done at the end of this subroutine.)
     gblclr = 2
     call mulbac2(temp,u,be,mu,bc,gblclr,ieo,ibleo,vecbl,vecblinv,myid,nn, &
                  ldiv,nms,lvbc,ib,lbd,iblv,MRT)
     call gammamultshort(temp, gtemp, 5) ! multiply by gamma5
     call vv(xo(:,:,:,ieo,ibleo),gtemp,Jo(:,:,mu))
 
    enddo ! mu

!*Local scalar operator [ psibar(x) psi(x) ] :
    call gammamultshort(be(:,:,:,ieo,ibleo), gtemp, 5) ! multiply by gamma5
    call vv(xe(:,:,:,ieo,ibleo),gtemp,Je(:,:,5))
    call gammamultshort(bo(:,:,:,ieo,ibleo), gtemp, 5) ! multiply by gamma5
    call vv(xo(:,:,:,ieo,ibleo),gtemp,Jo(:,:,5))
    if (ksub==1) then
     Je(1,:,5) = Je(1,:,5) + 12.0_KR ! 12 = 3 colours times 4 Dirac indices.
     Jo(1,:,5) = Jo(1,:,5) + 12.0_KR ! 12 = 3 colours times 4 Dirac indices.
    endif

! Local vector operator [ psibar(x) gamma_mu psi(x) ]:
   do mu = 1,4
     call gammamultshort(be(:,:,:,ieo,ibleo), gtemp, mu) ! multiply by gamma_mu
     call gammamultshort(gtemp, temp, 5)                 ! multiply by gamma5
     call vv(xe(:,:,:,ieo,ibleo),temp,Je(:,:,mu+5))
     call gammamultshort(bo(:,:,:,ieo,ibleo), gtemp, mu) ! multiply by gamma_mu
     call gammamultshort(gtemp, temp, 5)                 ! multiply by gamma5
     call vv(xo(:,:,:,ieo,ibleo),temp,Jo(:,:,mu+5))
   enddo ! mu


!*Now do the complex conjugation.
    Je(2,:,:) = -Je(2,:,:)
    Jo(2,:,:) = -Jo(2,:,:)

 end subroutine eigloopops_original
!-------------------------------------------------------------------------------------------------------


! THIS FUNCTION WORKS EXACTLY LIKE loopops() EXCEPT THAT THE SOURCE VECTOR
! IS EXPECTED TO HOLD THE SOLUTION OF gamma5Mx'=z NOT Mx=z FOR USE WITH THE
! HERMITIAN FORCED EIGENSPECTRUM SUBTRACTION TECHNIQUE
! NOTE: source vector is "b" and noise vector is "x"
 subroutine eigloopops(ksub,xe,xo,u,be,bo,Je,Jo,ieo,ibleo,bc,vecbl,vecblinv, &
                    myid,nn,ldiv,nms,lvbc,ib,lbd,iblv,MRT)
! This is a modification of the original version written by Vic because I think 
! gamma_5 multiplication has to be done before the multiplication with the operator.
!
! what we want is b^dagger Gamma_5 *O * xprime  where M*xprime=gamma_5*b
! For coding reasons we calculate the hermitian conjugate of this quantity
! then take the complex conjugate. When we take the Hermitian conjugate,
! we get xprime^dagge* O^dagger*Gamma_5*b
!
! -AA 11/28/2012


! Construct the operators to be included in disconnected loop contributions.

! INPUT:
!   ksub = 0 if the propagator has NOT been subtracted.
!   ksub = 1 if the propagator HAS been subtracted.
!   be() contains the source vector on even (gblclr=1) lattice sites.
!   bo() contains the source vector on odd (gblclr=2) lattice sites.
!        expected size: be(6,ntotal,4,2,8), bo(6,ntotal,4,2,8)
!                       where the first entry is real/imaginary and colour
!                       and the 3rd entry is the Dirac index
!                       and the other three entries give the lattice site.
!   xe() contains the sink vector on even (gblclr=1) lattice sites.
!   xo() contains the sink vector on odd (gblclr=2) lattice sites.
!        expected size: xe(6,ntotal,4,2,8), xo(6,ntotal,4,2,8)
!                       where the first entry is real/imaginary and colour
!                       and the 3rd entry is the Dirac index
!                       and the other three entries give the lattice site.
! OUTPUT:
!   Je() contains the operators on even (gblclr=1) lattice sites.
!   Jo() contains the operators on odd (gblclr=2) lattice sites.
!        expected size: Je(2,nvhalf,5), Jo(2,nvhalf,5)
!                       where the first entry is 1=real or 2=imaginary
!                       and the middle entry gives the lattice site
!                       and the last entry denotes one of the five operators.

    integer(kind=KI), intent(in)                 :: ksub, ieo, ibleo, myid, MRT
    real(kind=KR),    intent(inout), dimension(:,:,:,:,:) :: xe, xo
    real(kind=KR),    intent(in),    dimension(:,:,:,:,:) :: u
    real(kind=KR),    intent(inout), dimension(:,:,:,:,:) :: be, bo
    real(kind=KR),    intent(out),   dimension(:,:,:)     :: Je, Jo
    integer(kind=KI), intent(in),    dimension(:)         :: bc, nms
    integer(kind=KI), intent(in),    dimension(:,:)       :: vecbl, vecblinv
    integer(kind=KI), intent(in),    dimension(:,:)       :: nn, iblv
    logical,          intent(in),    dimension(:)         :: ldiv
    integer(kind=KI), intent(in),    dimension(:,:,:)     :: lvbc
    integer(kind=KI), intent(in),    dimension(:,:,:,:)   :: ib
    logical,          intent(in),    dimension(:,:)       :: lbd

    integer(kind=KI)                     :: gblclr, mu, ierr
    real(kind=KR), dimension(6,nvhalf,4) :: temp
    real(kind=KR), dimension(6,nvhalf,4) :: gtemp
    integer(kind=KI) :: id, isite, ibl
!*Initializations.
    Je = 0.0_KR
    Jo = 0.0_KR

!*Point-split vector operators:

   do mu = 1,4

! Je_mu(x)
! = Je_mu(x) + bebar(x) (1-gamma_mu) U_mu(x) xo(x+mu)
! = Je_mu(x) + [ xo^dagger(x+mu) (1-gamma_mu) U_mu^dagger(x) bebar^dagger(x) ]^*
! (The complex conjugation will be done at the end of this subroutine.)
     gblclr = 1
     !call mulfor2(temp,u,be,mu,gblclr,ieo,ibleo,vecbl)
     !call gammamultshort(temp, gtemp, 5) ! multiply by gamma5
     !call gammamultshort(be(:,:,:,ieo,ibleo), temp, 5) ! multiply by gamma5
     !call mulfor2(gtemp,u,temp,mu,gblclr,ieo,ibleo,vecbl)
     call mulfor2(temp,u,be,mu,gblclr,ieo,ibleo,vecbl)!BS 12/3/2015
     !temp(:,:,:)=-temp(:,:,:)BS 12/3/2015
     call gammamultshort(temp, gtemp, 5) ! multiply by gamma5
     call vgv(xo,gtemp,Je(:,:,mu),gblclr,ieo,ibleo,mu,bc,vecbl,vecblinv, &
              lbd,ldiv,ib,lvbc,nms,myid,nn,iblv,MRT)
 
! Jo_mu(x)
! = Jo_mu(x) + bobar(x) (1-gamma_mu) U_mu(x) xe(x+mu)
! = Jo_mu(x) + [ xe^dagger(x+mu) (1-gamma_mu) U_mu^dagger(x) bobar^dagger(x) ]^*
! (The complex conjugation will be done at the end of this subroutine.)
     gblclr = 2
     !call mulfor2(temp,u,bo,mu,gblclr,ieo,ibleo,vecbl)
     !call gammamultshort(temp, gtemp, 5) ! multiply by gamma5
     !call gammamultshort(bo(:,:,:,ieo,ibleo), temp, 5) ! multiply by gamma5
     !call mulfor2(gtemp,u,temp,mu,gblclr,ieo,ibleo,vecbl)
     call mulfor2(temp,u,bo,mu,gblclr,ieo,ibleo,vecbl)!BS 12/3/2015
    ! temp(:,:,:)=-temp(:,:,:)!BS
     call gammamultshort(temp, gtemp, 5) ! multiply by gamma5
     call vgv(xe,gtemp,Jo(:,:,mu),gblclr,ieo,ibleo,mu,bc,vecbl,vecblinv, &
              lbd,ldiv,ib,lvbc,nms,myid,nn,iblv,MRT)

! Je_mu(x)
! = Je_mu(x) - bobar(x+mu) (1+gamma_mu) U_mu^dagger(x) xe(x)
! = Je_mu(x) - [ xe^dagger(x) (1+gamma_mu) U_mu(x) bobar^dagger(x+mu) ]^*
! (The complex conjugation will be done at the end of this subroutine.)
     gblclr = 1
     !call mulbac2(temp,u,bo,mu,bc,gblclr,ieo,ibleo,vecbl,vecblinv,myid,nn, &
     !             ldiv,nms,lvbc,ib,lbd,iblv,MRT)
     !call gammamultshort(temp, gtemp, 5) ! multiply by gamma5
     !call gammamultshort(bo(:,:,:,ieo,ibleo), temp, 5) ! multiply by gamma5
     !call mulbac2(gtemp,u,temp,mu,bc,gblclr,ieo,ibleo,vecbl,vecblinv,myid,nn, &
     !             ldiv,nms,lvbc,ib,lbd,iblv,MRT)
     !bs call mulbac3(temp,u,bo,mu,bc,gblclr,ieo,ibleo,vecbl,vecblinv,myid,nn, &
     call mulbac2(temp,u,bo,mu,bc,gblclr,ieo,ibleo,vecbl,vecblinv,myid,nn, &
                  ldiv,nms,lvbc,ib,lbd,iblv,MRT)
     call gammamultshort(temp, gtemp, 5) ! multiply by gamma5
     call vv(xe(:,:,:,ieo,ibleo),gtemp,Je(:,:,mu))
 
! Jo_mu(x)
! = Jo_mu(x) - bebar(x+mu) (1+gamma_mu) U_mu^dagger(x) xo(x)
! = Jo_mu(x) - [ xo^dagger(x) (1+gamma_mu) U_mu(x) bebar^dagger(x+mu) ]^*
! (The complex conjugation will be done at the end of this subroutine.)
     gblclr = 2
     !call mulbac2(temp,u,be,mu,bc,gblclr,ieo,ibleo,vecbl,vecblinv,myid,nn, &
     !             ldiv,nms,lvbc,ib,lbd,iblv,MRT)
     !call gammamultshort(temp, gtemp, 5) ! multiply by gamma5
     !call gammamultshort(be(:,:,:,ieo,ibleo), temp, 5) ! multiply by gamma5
     !call mulbac2(gtemp,u,temp,mu,bc,gblclr,ieo,ibleo,vecbl,vecblinv,myid,nn, &
     !             ldiv,nms,lvbc,ib,lbd,iblv,MRT)
     call mulbac2(temp,u,be,mu,bc,gblclr,ieo,ibleo,vecbl,vecblinv,myid,nn, &
                  ldiv,nms,lvbc,ib,lbd,iblv,MRT)
     call gammamultshort(temp, gtemp, 5) ! multiply by gamma5
     call vv(xo(:,:,:,ieo,ibleo),gtemp,Jo(:,:,mu))
 
    enddo ! mu

!*Local scalar operator [ psibar(x) psi(x) ] :
    call gammamultshort(be(:,:,:,ieo,ibleo), gtemp, 5) ! multiply by gamma5
    call vv(xe(:,:,:,ieo,ibleo),gtemp,Je(:,:,5))
    call gammamultshort(bo(:,:,:,ieo,ibleo), gtemp, 5) ! multiply by gamma5
    call vv(xo(:,:,:,ieo,ibleo),gtemp,Jo(:,:,5))
    if (ksub==1) then
     Je(1,:,5) = Je(1,:,5) + 12.0_KR ! 12 = 3 colours times 4 Dirac indices.
     Jo(1,:,5) = Jo(1,:,5) + 12.0_KR ! 12 = 3 colours times 4 Dirac indices.
    endif

! Local vector operator [ psibar(x) gamma_mu psi(x) ]:
   do mu = 1,4
     !call gammamultshort(be(:,:,:,ieo,ibleo), gtemp, mu) ! multiply by gamma_mu
     !call gammamultshort(gtemp, temp, 5)                 ! multiply by gamma5
     !call gammamultshort(be(:,:,:,ieo,ibleo),temp, 5)                 ! multiply by gamma5
     !call gammamultshort(temp, gtemp, mu) ! multiply by gamma_mu
     call gammamultshort(be(:,:,:,ieo,ibleo), gtemp, mu) ! multiply by gamma_mu
     call gammamultshort(gtemp, temp, 5)                 ! multiply by gamma5
     !temp(:,:,:)=-temp(:,:,:)!BS 12/4/2015
     call vv(xe(:,:,:,ieo,ibleo),temp,Je(:,:,mu+5))
     !call gammamultshort(bo(:,:,:,ieo,ibleo), gtemp, mu) ! multiply by gamma_mu
     !call gammamultshort(gtemp, temp, 5)                 ! multiply by gamma5
     !call gammamultshort(bo(:,:,:,ieo,ibleo), temp, 5)                 ! multiply by gamma5
     !call gammamultshort(temp, gtemp, mu)                 ! multiply by gamma_mu
     call gammamultshort(bo(:,:,:,ieo,ibleo), gtemp, mu) ! multiply by gamma_mu
     call gammamultshort(gtemp, temp, 5)                 ! multiply by gamma5
    ! temp(:,:,:)=-temp(:,:,:)!BS
     call vv(xo(:,:,:,ieo,ibleo),temp,Jo(:,:,mu+5))
   enddo ! mu


!*Now do the complex conjugation.
    Je(2,:,:) = -Je(2,:,:)
    Jo(2,:,:) = -Jo(2,:,:)

 end subroutine eigloopops

!!-------------------------------------------------------------------------------------------------------



! THIS FUNCTION WORKS EXACTLY LIKE eigloopops() EXCEPT THAT I USE IT FOR 
! EXPERIEMNTING TO KEEP VIC'S CODE INTACT -Abdou
!
! Changes w.r.t. eigloopops:
!    - the noise vector need to be multiplied by gamma5 first before applying the operator

 subroutine eigloopops_abdou(ksub,xe,xo,u,be,bo,Je,Jo,ieo,ibleo,bc,vecbl,vecblinv, &
                    myid,nn,ldiv,nms,lvbc,ib,lbd,iblv,MRT)
! Construct the operators to be included in disconnected loop contributions.

! INPUT:
!   ksub = 0 if the propagator has NOT been subtracted.
!   ksub = 1 if the propagator HAS been subtracted.
!   be() contains the source vector on even (gblclr=1) lattice sites.
!   bo() contains the source vector on odd (gblclr=2) lattice sites.
!        expected size: be(6,ntotal,4,2,8), bo(6,ntotal,4,2,8)
!                       where the first entry is real/imaginary and colour
!                       and the 3rd entry is the Dirac index
!                       and the other three entries give the lattice site.
!   xe() contains the sink vector on even (gblclr=1) lattice sites.
!   xo() contains the sink vector on odd (gblclr=2) lattice sites.
!        expected size: xe(6,ntotal,4,2,8), xo(6,ntotal,4,2,8)
!                       where the first entry is real/imaginary and colour
!                       and the 3rd entry is the Dirac index
!                       and the other three entries give the lattice site.
! OUTPUT:
!   Je() contains the operators on even (gblclr=1) lattice sites.
!   Jo() contains the operators on odd (gblclr=2) lattice sites.
!        expected size: Je(2,nvhalf,9), Jo(2,nvhalf,9)
!                       where the first entry is 1=real or 2=imaginary
!                       and the middle entry gives the lattice site
!                       and the last entry denotes one of the five operators.

    integer(kind=KI), intent(in)                 :: ksub, ieo, ibleo, myid, MRT
    real(kind=KR),    intent(inout), dimension(:,:,:,:,:) :: xe, xo
    real(kind=KR),    intent(in),    dimension(:,:,:,:,:) :: u
    real(kind=KR),    intent(inout), dimension(:,:,:,:,:) :: be, bo
    real(kind=KR),    intent(out),   dimension(:,:,:)     :: Je, Jo
    integer(kind=KI), intent(in),    dimension(:)         :: bc, nms
    integer(kind=KI), intent(in),    dimension(:,:)       :: vecbl, vecblinv
    integer(kind=KI), intent(in),    dimension(:,:)       :: nn, iblv
    logical,          intent(in),    dimension(:)         :: ldiv
    integer(kind=KI), intent(in),    dimension(:,:,:)     :: lvbc
    integer(kind=KI), intent(in),    dimension(:,:,:,:)   :: ib
    logical,          intent(in),    dimension(:,:)       :: lbd

    integer(kind=KI)                     :: gblclr, mu, ierr
    real(kind=KR), dimension(6,ntotal,4,2,8) :: temp,gbe,gbo
    real(kind=KR), dimension(6,nvhalf,4) :: gtemp
    integer(kind=KI) :: id, isite, ibl, jeo, jbleo;


!*Initializations.
    Je = 0.0_KR
    Jo = 0.0_KR

!*Multiply the souce by gamma5
    do jeo=1,2
     do jbleo=1,8
       call gammamult(be, bo, gbe, gbo, 5,jeo,jbleo)
     enddo
    enddo


!*Point-split vector operators:

   do mu = 1,4

! Je_mu(x)
! = Je_mu(x) + bebar(x) (1-gamma_mu) U_mu(x) xo(x+mu)
! = Je_mu(x) + [ xo^dagger(x+mu) (1-gamma_mu) U_mu^dagger(x) bebar^dagger(x) ]^*
! (The complex conjugation will be done at the end of this subroutine.)
     gblclr = 1
     ! using gbe and gbo --Abdou
     !call mulfor2(temp(:,:,:,ieo,ibleo),u,be,mu,gblclr,ieo,ibleo,vecbl)
     !call gammamultshort(temp(:,:,:,ieo,ibleo), gtemp, 5) ! multiply by gamma5
     !call vgv(xo,gtemp,Je(:,:,mu),gblclr,ieo,ibleo,mu,bc,vecbl,vecblinv, &
     !         lbd,ldiv,ib,lvbc,nms,myid,nn,iblv,MRT)
     call mulfor2(temp(:,:,:,ieo,ibleo),u,gbe,mu,gblclr,ieo,ibleo,vecbl)
     call vgv(xo,temp(:,:,:,ieo,ibleo),Je(:,:,mu),gblclr,ieo,ibleo,mu,bc,vecbl,vecblinv, &
              lbd,ldiv,ib,lvbc,nms,myid,nn,iblv,MRT)
 
! Jo_mu(x)
! = Jo_mu(x) + bobar(x) (1-gamma_mu) U_mu(x) xe(x+mu)
! = Jo_mu(x) + [ xe^dagger(x+mu) (1-gamma_mu) U_mu^dagger(x) bobar^dagger(x) ]^*
! (The complex conjugation will be done at the end of this subroutine.)
     gblclr = 2
     !call mulfor2(temp(:,:,:,ieo,ibleo),u,bo,mu,gblclr,ieo,ibleo,vecbl)
     !call gammamultshort(temp(:,:,:,ieo,ibleo), gtemp, 5) ! multiply by gamma5
     !call vgv(xe,gtemp,Jo(:,:,mu),gblclr,ieo,ibleo,mu,bc,vecbl,vecblinv, &
     !         lbd,ldiv,ib,lvbc,nms,myid,nn,iblv,MRT)
     call mulfor2(temp(:,:,:,ieo,ibleo),u,gbo,mu,gblclr,ieo,ibleo,vecbl)
     call vgv(xe,temp(:,:,:,ieo,ibleo),Jo(:,:,mu),gblclr,ieo,ibleo,mu,bc,vecbl,vecblinv, &
              lbd,ldiv,ib,lvbc,nms,myid,nn,iblv,MRT)

! Je_mu(x)
! = Je_mu(x) - bobar(x+mu) (1+gamma_mu) U_mu^dagger(x) xe(x)
! = Je_mu(x) - [ xe^dagger(x) (1+gamma_mu) U_mu(x) bobar^dagger(x+mu) ]^*
! (The complex conjugation will be done at the end of this subroutine.)
     gblclr = 1
     !call mulbac2(temp(:,:,:,ieo,ibleo),u,bo,mu,bc,gblclr,ieo,ibleo,vecbl,vecblinv,myid,nn, &
     !             ldiv,nms,lvbc,ib,lbd,iblv,MRT)
     !call gammamultshort(temp(:,:,:,ieo,ibleo), gtemp, 5) ! multiply by gamma5
     !call vv(xe(:,:,:,ieo,ibleo),gtemp,Je(:,:,mu))
     call mulbac2(temp(:,:,:,ieo,ibleo),u,gbo,mu,bc,gblclr,ieo,ibleo,vecbl,vecblinv,myid,nn, &
                  ldiv,nms,lvbc,ib,lbd,iblv,MRT)
     call vv(xe(:,:,:,ieo,ibleo),temp(:,:,:,ieo,ibleo),Je(:,:,mu))
 
! Jo_mu(x)
! = Jo_mu(x) - bebar(x+mu) (1+gamma_mu) U_mu^dagger(x) xo(x)
! = Jo_mu(x) - [ xo^dagger(x) (1+gamma_mu) U_mu(x) bebar^dagger(x+mu) ]^*
! (The complex conjugation will be done at the end of this subroutine.)
     gblclr = 2
     !call mulbac2(temp(:,:,:,ieo,ibleo),u,be,mu,bc,gblclr,ieo,ibleo,vecbl,vecblinv,myid,nn, &
     !             ldiv,nms,lvbc,ib,lbd,iblv,MRT)
     !call gammamultshort(temp(:,:,:,ieo,ibleo), gtemp, 5) ! multiply by gamma5
     !call vv(xo(:,:,:,ieo,ibleo),gtemp,Jo(:,:,mu))
     call mulbac2(temp(:,:,:,ieo,ibleo),u,gbe,mu,bc,gblclr,ieo,ibleo,vecbl,vecblinv,myid,nn, &
                  ldiv,nms,lvbc,ib,lbd,iblv,MRT)
     call vv(xo(:,:,:,ieo,ibleo),temp(:,:,:,ieo,ibleo),Jo(:,:,mu))
 
    enddo ! mu

!*Local scalar operator [ psibar(x) psi(x) ] :
    call gammamultshort(be(:,:,:,ieo,ibleo), gtemp, 5) ! multiply by gamma5
    call vv(xe(:,:,:,ieo,ibleo),gtemp,Je(:,:,5))
    call gammamultshort(bo(:,:,:,ieo,ibleo), gtemp, 5) ! multiply by gamma5
    call vv(xo(:,:,:,ieo,ibleo),gtemp,Jo(:,:,5))
    !call gammamultshort(be(:,:,:,ieo,ibleo), gtemp, 5) ! multiply by gamma5
    !call vv(xe(:,:,:,ieo,ibleo),gtemp,Je(:,:,5))
    !call gammamultshort(bo(:,:,:,ieo,ibleo), gtemp, 5) ! multiply by gamma5
    !call vv(xo(:,:,:,ieo,ibleo),gtemp,Jo(:,:,5))
    if (ksub==1) then
     Je(1,:,5) = Je(1,:,5) + 12.0_KR ! 12 = 3 colours times 4 Dirac indices.
     Jo(1,:,5) = Jo(1,:,5) + 12.0_KR ! 12 = 3 colours times 4 Dirac indices.
    endif

! Local vector operator [ psibar(x) gamma_mu psi(x) ]:
   do mu = 1,4
     !call gammamultshort(be(:,:,:,ieo,ibleo), gtemp, mu) ! multiply by gamma_mu
     !call gammamultshort(gtemp, temp(:,:,:,ieo,ibleo), 5)                 ! multiply by gamma5
     !call vv(xe(:,:,:,ieo,ibleo),temp(:,:,:,ieo,ibleo),Je(:,:,mu+5))
     call gammamultshort(gbe(:,:,:,ieo,ibleo), gtemp, mu) ! multiply by gamma_mu
     call vv(xe(:,:,:,ieo,ibleo),gtemp,Je(:,:,mu+5))
     !call gammamultshort(bo(:,:,:,ieo,ibleo), gtemp, mu) ! multiply by gamma_mu
     !call gammamultshort(gtemp, temp(:,:,:,ieo,ibleo), 5)                 ! multiply by gamma5
     !call vv(xo(:,:,:,ieo,ibleo),temp(:,:,:,ieo,ibleo),Jo(:,:,mu+5))
     call gammamultshort(gbo(:,:,:,ieo,ibleo), gtemp, mu) ! multiply by gamma_mu
     call vv(xo(:,:,:,ieo,ibleo),gtemp,Jo(:,:,mu+5))
   enddo ! mu


!*Now do the complex conjugation.
    Je(2,:,:) = -Je(2,:,:)
    Jo(2,:,:) = -Jo(2,:,:)

 end subroutine eigloopops_abdou
!-------------------------------------------------------------------------------------------------------
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 subroutine mulfor2(h,u,g,mu,gblclr,ieo,ibleo,vecbl)
! Compute the forward propagation term from the point-split vector operator:
! h(x) = (1 - gamma_mu) U_mu^dagger(x) g(x)
! INPUT:
!   u() contains the gauge fields for this sublattice.
!   g() contains the Dirac spinor for this sublattice.
!   mu is the direction of the gauge field (and Dirac gamma).
!   bc(mu) = 1,-1,0 for periodic,antiperiodic,fixed boundary conditions
!            respectively.
!   gblclr is the "global checkerboard colour" of x,
!          and thus of h(x) and g(x).
!          NOTE: This is NOT ieo, because of the lattice's blocked structure.
!                gblclr=1 means sites on blocks 1,4,6,7,10,11,13,16.
!                gblclr=2 means sites on blocks 2,3,5,8,9,12,14,15.
!   vecbl() defines the checkerboarding of the entire lattice.
!           In particular, the blocks vecbl(1,:), including ieo=1 and ieo=2,
!           form a checkerboard pattern.  So do the complement: vecbl(2,:).
! OUTPUT:
!   h(x)
!   expected size: h(a,b,c),
!                  where a=1..6 is the real/imaginary and colour index.
!                        b=1,2,3,...,nvhalf combines colour and position.
!                        c=1,2,3,4 is the Dirac spinor index.

    real(kind=KR),    intent(out), dimension(:,:,:)     :: h
    real(kind=KR),    intent(in),  dimension(:,:,:,:,:) :: u, g
    integer(kind=KI), intent(in)                      :: mu, gblclr, ieo, ibleo
    integer(kind=KI), intent(in),  dimension(:,:)       :: vecbl

    integer(kind=KI)                     :: ibl, id, i, ir, icount
    real(kind=KR), dimension(6,nvhalf,4) :: wx
    integer(kind=KI)                     :: isite, myid, ierr

! Perform the matrix-vector multiplication: h(x) = U_mu^dagger(x) * g(x).
    icount = nvhalf
    ibl = vecbl(gblclr,ibleo)

    do id = 1,4
       call mdv(icount,u(:,:,mu,ieo,ibl),g(:,:,id,ieo,ibleo),wx(:,:,id))
    enddo ! id

! Construct h() = (1 - gamma_mu) * wx().
    select case(mu)
     case(1)
      do i = 1,nvhalf
       do ir = 1,6 ! 6=nri*nc
        h(ir,i,1) = wx(ir,i,1) - wx(ir,i,4)
        h(ir,i,2) = wx(ir,i,2) - wx(ir,i,3)
        h(ir,i,3) = wx(ir,i,3) - wx(ir,i,2)
        h(ir,i,4) = wx(ir,i,4) - wx(ir,i,1)
       enddo ! ir
      enddo ! i
     case(2)
      do i = 1,nvhalf
       do ir = 1,5,2 ! 6=nri*nc
        h(ir  ,i,1) = wx(ir  ,i,1) - wx(ir+1,i,4)
        h(ir+1,i,1) = wx(ir+1,i,1) + wx(ir  ,i,4)
        h(ir  ,i,2) = wx(ir  ,i,2) + wx(ir+1,i,3)
        h(ir+1,i,2) = wx(ir+1,i,2) - wx(ir  ,i,3)
        h(ir  ,i,3) = wx(ir  ,i,3) - wx(ir+1,i,2)
        h(ir+1,i,3) = wx(ir+1,i,3) + wx(ir  ,i,2)
        h(ir  ,i,4) = wx(ir  ,i,4) + wx(ir+1,i,1)
        h(ir+1,i,4) = wx(ir+1,i,4) - wx(ir  ,i,1)
       enddo ! ir
      enddo ! i
     case(3)
      do i = 1,nvhalf
       do ir = 1,6 ! 6=nri*nc
        h(ir,i,1) = wx(ir,i,1) - wx(ir,i,3)
        h(ir,i,2) = wx(ir,i,2) + wx(ir,i,4)
        h(ir,i,3) = wx(ir,i,3) - wx(ir,i,1)
        h(ir,i,4) = wx(ir,i,4) + wx(ir,i,2)
       enddo ! ir
      enddo ! i
     case(4)
      do i = 1,nvhalf
       do ir = 1,6 ! 6=nri*nc
        h(ir,i,1) = 0.0_KR
        h(ir,i,2) = 0.0_KR
        h(ir,i,3) = 2.0_KR*wx(ir,i,3)
        h(ir,i,4) = 2.0_KR*wx(ir,i,4)
       enddo ! ir
      enddo ! i
     case default
      open(unit=8,file="DISCONLOOPS.ERROR",action="write",status="replace", &
           form="formatted")
       write(unit=8,fmt=*) "subroutine mulfor2: mu =", mu
      close(unit=8,status="keep")
      stop
    end select

 end subroutine mulfor2







 subroutine mulfor3(h,u,g,mu,gblclr,ieo,ibleo,vecbl)
! Compute the forward propagation term from the point-split vector operator:
! h(x) = (1 + gamma_mu) U_mu^dagger(x) g(x)
! INPUT:
!   u() contains the gauge fields for this sublattice.
!   g() contains the Dirac spinor for this sublattice.
!   mu is the direction of the gauge field (and Dirac gamma).
!   bc(mu) = 1,-1,0 for periodic,antiperiodic,fixed boundary conditions
!            respectively.
!   gblclr is the "global checkerboard colour" of x,
!          and thus of h(x) and g(x).
!          NOTE: This is NOT ieo, because of the lattice's blocked structure.
!                gblclr=1 means sites on blocks 1,4,6,7,10,11,13,16.
!                gblclr=2 means sites on blocks 2,3,5,8,9,12,14,15.
!   vecbl() defines the checkerboarding of the entire lattice.
!           In particular, the blocks vecbl(1,:), including ieo=1 and ieo=2,
!           form a checkerboard pattern.  So do the complement: vecbl(2,:).
! OUTPUT:
!   h(x)
!   expected size: h(a,b,c),
!                  where a=1..6 is the real/imaginary and colour index.
!                        b=1,2,3,...,nvhalf combines colour and position.
!                        c=1,2,3,4 is the Dirac spinor index.

    real(kind=KR),    intent(out), dimension(:,:,:)     :: h
    real(kind=KR),    intent(in),  dimension(:,:,:,:,:) :: u, g
    integer(kind=KI), intent(in)                      :: mu, gblclr, ieo, ibleo
    integer(kind=KI), intent(in),  dimension(:,:)       :: vecbl

    integer(kind=KI)                     :: ibl, id, i, ir, icount
    real(kind=KR), dimension(6,nvhalf,4) :: wx
    integer(kind=KI)                     :: isite, myid, ierr

! Perform the matrix-vector multiplication: h(x) = U_mu^dagger(x) * g(x).
    icount = nvhalf
    ibl = vecbl(gblclr,ibleo)

    do id = 1,4
       call mdv(icount,u(:,:,mu,ieo,ibl),g(:,:,id,ieo,ibleo),wx(:,:,id))
    enddo ! id

! Construct h() = (1 + gamma_mu) * wx().
    select case(mu)
     case(1)
      do i = 1,nvhalf
       do ir = 1,6 ! 6=nri*nc
        h(ir,i,1) = wx(ir,i,1) + wx(ir,i,4)
        h(ir,i,2) = wx(ir,i,2) + wx(ir,i,3)
        h(ir,i,3) = wx(ir,i,3) + wx(ir,i,2)
        h(ir,i,4) = wx(ir,i,4) + wx(ir,i,1)
       enddo ! ir
      enddo ! i
     case(2)
      do i = 1,nvhalf
       do ir = 1,5,2 ! 6=nri*nc
        h(ir  ,i,1) = wx(ir  ,i,1) + wx(ir+1,i,4)
        h(ir+1,i,1) = wx(ir+1,i,1) - wx(ir  ,i,4)
        h(ir  ,i,2) = wx(ir  ,i,2) - wx(ir+1,i,3)
        h(ir+1,i,2) = wx(ir+1,i,2) + wx(ir  ,i,3)
        h(ir  ,i,3) = wx(ir  ,i,3) + wx(ir+1,i,2)
        h(ir+1,i,3) = wx(ir+1,i,3) - wx(ir  ,i,2)
        h(ir  ,i,4) = wx(ir  ,i,4) - wx(ir+1,i,1)
        h(ir+1,i,4) = wx(ir+1,i,4) + wx(ir  ,i,1)
       enddo ! ir
      enddo ! i
     case(3)
      do i = 1,nvhalf
       do ir = 1,6 ! 6=nri*nc
        h(ir,i,1) = wx(ir,i,1) + wx(ir,i,3)
        h(ir,i,2) = wx(ir,i,2) - wx(ir,i,4)
        h(ir,i,3) = wx(ir,i,3) + wx(ir,i,1)
        h(ir,i,4) = wx(ir,i,4) - wx(ir,i,2)
       enddo ! ir
      enddo ! i
     case(4)
      do i = 1,nvhalf
       do ir = 1,6 ! 6=nri*nc
        h(ir,i,1) = 2.0_KR*wx(ir,i,1)
        h(ir,i,2) = 2.0_KR*wx(ir,i,2)
        h(ir,i,3) = 0.0_KR
        h(ir,i,4) = 0.0_KR
       enddo ! ir
      enddo ! i
     case default
      open(unit=8,file="DISCONLOOPS.ERROR",action="write",status="replace", &
           form="formatted")
       write(unit=8,fmt=*) "subroutine mulfor3: mu =", mu
      close(unit=8,status="keep")
      stop
    end select
 end subroutine mulfor3




! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 




 subroutine mulbac2(h,u,g,mu,bc,gblclr,ieo,ibleo,vecbl,vecblinv,myid,nn,ldiv, &
                    nms,lvbc,ib,lbd,iblv,MRT)
! Comute the backward propagation term from the action:
! h(x) = -(1 + gamma_mu) U_mu(x) g(x+mu)
! INPUT:
!   u() contains the gauge fields for this sublattice.
!   g() contains the Dirac spinor for this sublattice.
!   mu is the direction of the gauge field (and Dirac gamma).
!   bc(mu) = 1,-1,0 for periodic,antiperiodic,fixed boundary conditions
!            respectively.
!   gblclr is the "global checkerboard colour" of x.
!          x+mu has the opposite checkerboard colour.
!          NOTE: gblclr=1 means sites on blocks 1,4,6,7,10,11,13,16.
!                gblclr=2 means sites on blocks 2,3,5,8,9,12,14,15.
!   vecbl() defines the checkerboarding of the entire lattice.
!           In particular, the blocks vecbl(1,:), including ieo=1 and ieo=2,
!           form a checkerboard pattern.  So do the complement: vecbl(2,:).
! OUTPUT:
!   h(x)
!   expected size: h(a,b,c),
!                  where a=1..6 is the real/imaginary and colour index.
!                        b=1,2,3,...,nvhalf combines colour and position.
!                        c=1,2,3,4 is the Dirac spinor index.
 
    real(kind=KR),    intent(out),   dimension(:,:,:)     :: h
    real(kind=KR),    intent(inout), dimension(:,:,:,:,:) :: g
    real(kind=KR),    intent(in),    dimension(:,:,:,:,:) :: u
    integer(kind=KI), intent(in)           :: mu, gblclr, ieo, ibleo, myid, MRT
    integer(kind=KI), intent(in),    dimension(:)         :: bc, nms
    integer(kind=KI), intent(in),    dimension(:,:)       :: vecbl, vecblinv
    integer(kind=KI), intent(in),    dimension(:,:)       :: nn, iblv
    logical,          intent(in),    dimension(:)         :: ldiv
    integer(kind=KI), intent(in),    dimension(:,:,:)     :: lvbc
    integer(kind=KI), intent(in),    dimension(:,:,:,:)   :: ib
    logical,          intent(in),    dimension(:,:)       :: lbd
 
    integer(kind=KI)                     :: jeo, ibl, jbl, jbleo, id, i, ir
    real(kind=KR), dimension(6,nvhalf,4) :: wx
 
! Perform the matrix-vector multiplication: h(x) = U_mu(x) * g(x+mu).
    jeo = 3 - ieo
    ibl = vecbl(gblclr,ibleo)
    jbl = iblv(ibl,mu)
    jbleo = vecblinv(2,jbl)
    do id = 1,4
     call setmvf(mu,bc,lbd(ibl,mu),g(:,:,id,ieo,jbleo),         &
                 ldiv(mu),ib(:,mu,1,jeo),u(:,:,mu,ieo,ibl),     &
                 g(:,:,id,jeo,jbleo),wx(:,:,id),lvbc(:,mu,ieo), &
                 1,nms,myid,nn,MRT)
    enddo ! id
 
! Construct h() = -(1 + gamma_mu) * wx().
    select case(mu)
     case(1)
      do i = 1,nvhalf
       do ir = 1,6 ! 6=nri*nc
        h(ir,i,1) = -wx(ir,i,1) - wx(ir,i,4)
        h(ir,i,2) = -wx(ir,i,2) - wx(ir,i,3)
        h(ir,i,3) = -wx(ir,i,3) - wx(ir,i,2)
        h(ir,i,4) = -wx(ir,i,4) - wx(ir,i,1)
       enddo ! ir
      enddo ! i
     case(2)
      do i = 1,nvhalf
       do ir = 1,5,2 ! 6=nri*nc
        h(ir  ,i,1) = -wx(ir  ,i,1) - wx(ir+1,i,4)
        h(ir+1,i,1) = -wx(ir+1,i,1) + wx(ir  ,i,4)
        h(ir  ,i,2) = -wx(ir  ,i,2) + wx(ir+1,i,3)
        h(ir+1,i,2) = -wx(ir+1,i,2) - wx(ir  ,i,3)
        h(ir  ,i,3) = -wx(ir  ,i,3) - wx(ir+1,i,2)
        h(ir+1,i,3) = -wx(ir+1,i,3) + wx(ir  ,i,2)
        h(ir  ,i,4) = -wx(ir  ,i,4) + wx(ir+1,i,1)
        h(ir+1,i,4) = -wx(ir+1,i,4) - wx(ir  ,i,1)
       enddo ! ir
      enddo ! i
     case(3)
      do i = 1,nvhalf
       do ir = 1,6 ! 6=nri*nc
        h(ir,i,1) = -wx(ir,i,1) - wx(ir,i,3)
        h(ir,i,2) = -wx(ir,i,2) + wx(ir,i,4)
        h(ir,i,3) = -wx(ir,i,3) - wx(ir,i,1)
        h(ir,i,4) = -wx(ir,i,4) + wx(ir,i,2)
       enddo ! ir
      enddo ! i
     case(4)
      do i = 1,nvhalf
       do ir = 1,6 ! 6=nri*nc
        h(ir,i,1) = -2.0_KR*wx(ir,i,1)
        h(ir,i,2) = -2.0_KR*wx(ir,i,2)
        h(ir,i,3) = 0.0_KR
        h(ir,i,4) = 0.0_KR
       enddo ! ir
      enddo ! i
     case default
      open(unit=8,file="DISCONLOOPS.ERROR",action="write",status="replace", &
           form="formatted")
       write(unit=8,fmt=*) "subroutine mulbac2: mu =", mu
      close(unit=8,status="keep")
      stop
    end select
 
 end subroutine mulbac2



 subroutine mulbac3(h,u,g,mu,bc,gblclr,ieo,ibleo,vecbl,vecblinv,myid,nn,ldiv, &
                    nms,lvbc,ib,lbd,iblv,MRT)
! Comute the backward propagation term from the action:
! h(x) = (1 - gamma_mu) U_mu(x) g(x+mu)
! INPUT:
!   u() contains the gauge fields for this sublattice.
!   g() contains the Dirac spinor for this sublattice.
!   mu is the direction of the gauge field (and Dirac gamma).
!   bc(mu) = 1,-1,0 for periodic,antiperiodic,fixed boundary conditions
!            respectively.
!   gblclr is the "global checkerboard colour" of x.
!          x+mu has the opposite checkerboard colour.
!          NOTE: gblclr=1 means sites on blocks 1,4,6,7,10,11,13,16.
!                gblclr=2 means sites on blocks 2,3,5,8,9,12,14,15.
!   vecbl() defines the checkerboarding of the entire lattice.
!           In particular, the blocks vecbl(1,:), including ieo=1 and ieo=2,
!           form a checkerboard pattern.  So do the complement: vecbl(2,:).
! OUTPUT:
!   h(x)
!   expected size: h(a,b,c),
!                  where a=1..6 is the real/imaginary and colour index.
!                        b=1,2,3,...,nvhalf combines colour and position.
!                        c=1,2,3,4 is the Dirac spinor index.
 
    real(kind=KR),    intent(out),   dimension(:,:,:)     :: h
    real(kind=KR),    intent(inout), dimension(:,:,:,:,:) :: g
    real(kind=KR),    intent(in),    dimension(:,:,:,:,:) :: u
    integer(kind=KI), intent(in)           :: mu, gblclr, ieo, ibleo, myid, MRT
    integer(kind=KI), intent(in),    dimension(:)         :: bc, nms
    integer(kind=KI), intent(in),    dimension(:,:)       :: vecbl, vecblinv
    integer(kind=KI), intent(in),    dimension(:,:)       :: nn, iblv
    logical,          intent(in),    dimension(:)         :: ldiv
    integer(kind=KI), intent(in),    dimension(:,:,:)     :: lvbc
    integer(kind=KI), intent(in),    dimension(:,:,:,:)   :: ib
    logical,          intent(in),    dimension(:,:)       :: lbd
 
    integer(kind=KI)                     :: jeo, ibl, jbl, jbleo, id, i, ir
    real(kind=KR), dimension(6,nvhalf,4) :: wx
 
! Perform the matrix-vector multiplication: h(x) = U_mu(x) * g(x+mu).
    jeo = 3 - ieo
    ibl = vecbl(gblclr,ibleo)
    jbl = iblv(ibl,mu)
    jbleo = vecblinv(2,jbl)
    do id = 1,4
     call setmvf(mu,bc,lbd(ibl,mu),g(:,:,id,ieo,jbleo),         &
                 ldiv(mu),ib(:,mu,1,jeo),u(:,:,mu,ieo,ibl),     &
                 g(:,:,id,jeo,jbleo),wx(:,:,id),lvbc(:,mu,ieo), &
                 1,nms,myid,nn,MRT)
    enddo ! id
 
! Construct h() = (1 - gamma_mu) * wx().!test
    select case(mu)
     case(1)
      do i = 1,nvhalf
       do ir = 1,6 ! 6=nri*nc
        h(ir,i,1) = wx(ir,i,1) - wx(ir,i,4)
        h(ir,i,2) = wx(ir,i,2) - wx(ir,i,3)
        h(ir,i,3) = wx(ir,i,3) - wx(ir,i,2)
        h(ir,i,4) = wx(ir,i,4) - wx(ir,i,1)
       enddo ! ir
      enddo ! i
     case(2)
      do i = 1,nvhalf
       do ir = 1,5,2 ! 6=nri*nc
        h(ir  ,i,1) = wx(ir  ,i,1) - wx(ir+1,i,4)
        h(ir+1,i,1) = wx(ir+1,i,1) + wx(ir  ,i,4)
        h(ir  ,i,2) = wx(ir  ,i,2) + wx(ir+1,i,3)
        h(ir+1,i,2) = wx(ir+1,i,2) - wx(ir  ,i,3)
        h(ir  ,i,3) = wx(ir  ,i,3) - wx(ir+1,i,2)
        h(ir+1,i,3) = wx(ir+1,i,3) + wx(ir  ,i,2)
        h(ir  ,i,4) = wx(ir  ,i,4) + wx(ir+1,i,1)
        h(ir+1,i,4) = wx(ir+1,i,4) - wx(ir  ,i,1)
       enddo ! ir
      enddo ! i
     case(3)
      do i = 1,nvhalf
       do ir = 1,6 ! 6=nri*nc
        h(ir,i,1) = wx(ir,i,1) - wx(ir,i,3)
        h(ir,i,2) = wx(ir,i,2) + wx(ir,i,4)
        h(ir,i,3) = wx(ir,i,3) - wx(ir,i,1)
        h(ir,i,4) = wx(ir,i,4) + wx(ir,i,2)
       enddo ! ir
      enddo ! i
     case(4)
      do i = 1,nvhalf
       do ir = 1,6 ! 6=nri*nc
        h(ir,i,1) = 0.0_KR
        h(ir,i,2) = 0.0_KR
        h(ir,i,3) = 2.0_KR*wx(ir,i,3)
        h(ir,i,4) = 2.0_KR*wx(ir,i,4)
       enddo ! ir
      enddo ! i
     case default
      open(unit=8,file="DISCONLOOPS.ERROR",action="write",status="replace", &
           form="formatted")
       write(unit=8,fmt=*) "subroutine mulbac3: mu =", mu
      close(unit=8,status="keep")
      stop
    end select
 
 end subroutine mulbac3





! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 
 subroutine vv(v1,v2,w)
! Vector*vector multiplier: calculate
! w(i) = w(i) + v1^dagger(i) delta_{ic,colour} delta_{id,Dirac} v2(i)
! Note that there is no factor of gamma_4 between v1^dagger and v2, because
! v1^dagger and v2 are understood to be two ends of a fermion propagator
! and, as usual, should really be thought of as psibar (implicit gamma_4)
! and psi.

! INPUT:
!   v1() and v2() are vectors of the form     / v(1)+I*v(2) \
!                                         v = | v(3)+I*v(4) |
!                                             \ v(5)+I*v(6) /
!   in colour space.
!   There is one such vector for each of the 4 Dirac components.
!   expected size: v1(6,nvhalf,4), v2(6,nvhalf,4)
!   nvhalf is the number of v1 vectors = the number of v2 vectors.
! OUTPUT:
!   w(1,i) = w(1,i) + Re[ v1^dagger(i) v2(i)]
!   w(2,i) = w(2,i) + Im[ v1^dagger(i) v2(i)]
!   expected size: w(2,nvhalf) or w(2,ntotal)

    real(kind=KR),    intent(in),    dimension(:,:,:) :: v1, v2
    real(kind=KR),    intent(inout), dimension(:,:)   :: w

    integer(kind=KI) :: id, i

    do id = 1,4
     do i = 1,nvhalf
      w(1,i) = w(1,i) + v1(1,i,id)*v2(1,i,id) + v1(2,i,id)*v2(2,i,id) &
                      + v1(3,i,id)*v2(3,i,id) + v1(4,i,id)*v2(4,i,id) &
                      + v1(5,i,id)*v2(5,i,id) + v1(6,i,id)*v2(6,i,id)
      w(2,i) = w(2,i) + v1(1,i,id)*v2(2,i,id) - v1(2,i,id)*v2(1,i,id) &
                      + v1(3,i,id)*v2(4,i,id) - v1(4,i,id)*v2(3,i,id) &
                      + v1(5,i,id)*v2(6,i,id) - v1(6,i,id)*v2(5,i,id)
     enddo ! i
    enddo ! id

 end subroutine vv

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 subroutine vgv(v1,v2,w,gblclr,ieo,ibleo,mu,bc,vecbl,vecblinv,lbd,ldiv,ib, &
                lvbc,nms,myid,nn,iblv,MRT)
! Vector*vector multiplier: calculate
! w(i) = w(i) + v1^dagger(j) delta_{ic,colour} delta_{id,Dirac} v2(i)
! where j is the nearest neighbour of i in the +mu direction.
! Note that there is no factor of gamma_4 between v1^dagger and v2, because
! v1^dagger and v2 are understood to be two ends of a fermion propagator
! and, as usual, should really be thought of as psibar (implicit gamma_4)
! and psi.

! INPUT:
!   v1() and v2() are vectors of the form     / v(1)+I*v(2) \
!                                         v = | v(3)+I*v(4) |
!                                             \ v(5)+I*v(6) /
!   in colour space.
!   There is one such vector for each of the 4 Dirac components.
!   expected size: v1(6,ntotal,4,2,8), v2(6,nvhalf,4,2,8)
!   nvhalf is the number of v2 vectors.
!   ntotal>nvhalf since it includes the buffer region for process boundaries.
! OUTPUT:
!   w(1,i) = w(1,i) + Re[ v1^dagger(j) v2(i)]
!   w(2,i) = w(2,i) + Im[ v1^dagger(j) v2(i)]
!   expected size: w(2,nvhalf,4)

    real(kind=KR),    intent(inout), dimension(:,:,:,:,:) :: v1
    real(kind=KR),    intent(in),    dimension(:,:,:)     :: v2
    real(kind=KR),    intent(inout), dimension(:,:)       :: w
    integer(kind=KI), intent(in)           :: gblclr, ieo, ibleo, mu, myid, MRT
    integer(kind=KI), intent(in),    dimension(:)         :: bc, nms
    integer(kind=KI), intent(in),    dimension(:,:)       :: vecbl, vecblinv
    logical,          intent(in),    dimension(:,:)       :: lbd
    logical,          intent(in),    dimension(:)         :: ldiv
    integer(kind=KI), intent(in),    dimension(:,:,:,:)   :: ib
    integer(kind=KI), intent(in),    dimension(:,:,:)     :: lvbc
    integer(kind=KI), intent(in),    dimension(:,:)       :: nn, iblv

    integer(kind=KI) :: jeo, ibl, jbl, jbleo, id, i, j
    integer(kind=KI) :: ic
 
    jeo = 3 - ieo
    ibl = vecbl(gblclr,ibleo)
    jbl = iblv(ibl,mu)
    jbleo = vecblinv(2,jbl)

!-If this is the first of the two blocks in this direction, then no vectors
! are needed from a neighbouring process.
    if (lbd(ibl,mu)) then
     do id = 1,4
      do i = 1,nvhalf
       w(1,i) = w(1,i) + v1(1,i,id,ieo,jbleo)*v2(1,i,id) &
                       + v1(2,i,id,ieo,jbleo)*v2(2,i,id) &
                       + v1(3,i,id,ieo,jbleo)*v2(3,i,id) &
                       + v1(4,i,id,ieo,jbleo)*v2(4,i,id) &
                       + v1(5,i,id,ieo,jbleo)*v2(5,i,id) &
                       + v1(6,i,id,ieo,jbleo)*v2(6,i,id)
       w(2,i) = w(2,i) + v1(1,i,id,ieo,jbleo)*v2(2,i,id) &
                       - v1(2,i,id,ieo,jbleo)*v2(1,i,id) &
                       + v1(3,i,id,ieo,jbleo)*v2(4,i,id) &
                       - v1(4,i,id,ieo,jbleo)*v2(3,i,id) &
                       + v1(5,i,id,ieo,jbleo)*v2(6,i,id) &
                       - v1(6,i,id,ieo,jbleo)*v2(5,i,id)
      enddo ! i
     enddo ! id
!-If this is the second of the two blocks in this direction, then get the
! vectors from a neighbouring process if...
! ...there are multiple processes in the mu direction [ldiv(mu)=.true.],
! OR
! ...there is a single process in the mu direction but with non-periodic
!    boundary conditions [(.not.ldiv(mu)).and.bc(mu)/=1].
    else
     do id = 1,4
      if ( ldiv(mu) .or. ((.not.ldiv(mu)).and.bc(mu)/=1) )                &
       call setvector(v1(:,:,id,jeo,jbleo),ib(:,mu,1,jeo),mu,bc,nms,myid, &
                      nn,ldiv(mu),MRT)
      do i = 1,nvhalf
       j = lvbc(i,mu,ieo)
       w(1,i) = w(1,i) + v1(1,j,id,jeo,jbleo)*v2(1,i,id) &
                       + v1(2,j,id,jeo,jbleo)*v2(2,i,id) &
                       + v1(3,j,id,jeo,jbleo)*v2(3,i,id) &
                       + v1(4,j,id,jeo,jbleo)*v2(4,i,id) &
                       + v1(5,j,id,jeo,jbleo)*v2(5,i,id) &
                       + v1(6,j,id,jeo,jbleo)*v2(6,i,id)
       w(2,i) = w(2,i) + v1(1,j,id,jeo,jbleo)*v2(2,i,id) &
                       - v1(2,j,id,jeo,jbleo)*v2(1,i,id) &
                       + v1(3,j,id,jeo,jbleo)*v2(4,i,id) &
                       - v1(4,j,id,jeo,jbleo)*v2(3,i,id) &
                       + v1(5,j,id,jeo,jbleo)*v2(6,i,id) &
                       - v1(6,j,id,jeo,jbleo)*v2(5,i,id)
      enddo ! i
     enddo ! id
    endif

 end subroutine vgv

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 subroutine momfacs(io,momfac,myid,iblv)
! Define the functions that are the momentum dependences for disconnected loops.
! INPUT:
!   io(mu) is the location of the origin on the global (multi-process)
!          lattice in the mu direction.
!          expected size: io(1)=iox, io(2)=ioy, io(3)=ioz.  No 4th component.
! OUTPUT:
!   momfac(isite,ieo,ibl,imom,iop),
!          where the first 3 entries are the lattice site
!          and the 4th entry is the momentum
!          and the last entry denotes one of the five operators.

    integer(kind=KI), intent(in),  dimension(:)         :: io
    real(kind=KR),    intent(out), dimension(:,:,:,:,:) :: momfac
    integer(kind=KI), intent(in)                        :: myid
    integer(kind=KI), intent(in),  dimension(:,:)       :: iblv

    integer(kind=KI)               :: nx1, ny1, nz1, nt1, ixbit, iybit, izbit,&
                                      itbit, ix, iy, iz, it, intx, inty, intz,&
                                      intt, isite, ieo, ibl
    integer(kind=KI), dimension(4) :: np, ip
    real(kind=KR)                  :: dx, dy, dz
    real(kind=KR2),   parameter    :: twoPi=6.283185307179586

    nx1 = nx/npx
    ny1 = ny/npy
    nz1 = nz/npz
    nt1 = nt/npt
    np(1) = npx
    np(2) = npy
    np(3) = npz
    np(4) = npt
    call atoc(myid,np,ip)
    do itbit = 1,nt1
     it = itbit + ip(4)*nt1
     do izbit = 1,nz1
      iz = izbit + ip(3)*nz1
      dz = twoPi*real(iz-io(3),KR)/real(nz,KR2)
      do iybit = 1,ny1
       iy = iybit + ip(2)*ny1
       dy = twoPi*real(iy-io(2),KR)/real(ny,KR2)
       do ixbit = 1,nx1
        ix = ixbit + ip(1)*nx1
        dx = twoPi*real(ix-io(1),KR)/real(nx,KR2)
        intx = (ixbit-1)/4
        inty = (iybit-1)/2
        intz = (izbit-1)/2
        intt = (itbit-1)/2
        isite = 1 + intx + inty*nx1/4 + intz*nx1*ny1/8 + intt*nx1*ny1*nz1/16
        ieo = 1
        if (modulo(ix,4)==3.or.modulo(ix,4)==0) ieo = 3 - ieo
        if (modulo(iy,4)==3.or.modulo(iy,4)==0) ieo = 3 - ieo
        if (modulo(iz,4)==3.or.modulo(iz,4)==0) ieo = 3 - ieo
        if (modulo(it,4)==3.or.modulo(it,4)==0) ieo = 3 - ieo
        ibl = 1
        if (modulo(ix,2)==0) ibl = iblv(ibl,1)
        if (modulo(iy,2)==0) ibl = iblv(ibl,2)
        if (modulo(iz,2)==0) ibl = iblv(ibl,3)
        if (modulo(it,2)==0) ibl = iblv(ibl,4)
        momfac(isite,ieo,ibl,1,1) = sin(dx)
        momfac(isite,ieo,ibl,1,2) = sin(dy)
        momfac(isite,ieo,ibl,1,3) = sin(dz)
        momfac(isite,ieo,ibl,1,4) = (cos(dx)+cos(dy)+cos(dz))/3.0_KR
        momfac(isite,ieo,ibl,2,1) = sin(dx)*(cos(dy)+cos(dz))/2.0_KR
        momfac(isite,ieo,ibl,2,2) = sin(dy)*(cos(dx)+cos(dz))/2.0_KR
        momfac(isite,ieo,ibl,2,3) = sin(dz)*(cos(dx)+cos(dy))/2.0_KR
        momfac(isite,ieo,ibl,2,4) = (cos(dx)*cos(dy)+cos(dx)*cos(dz) &
                                    +cos(dy)*cos(dz))/3.0_KR
        momfac(isite,ieo,ibl,3,1) = sin(dx)*cos(dy)*cos(dz)
        momfac(isite,ieo,ibl,3,2) = sin(dy)*cos(dx)*cos(dz)
        momfac(isite,ieo,ibl,3,3) = sin(dz)*cos(dx)*cos(dy)
        momfac(isite,ieo,ibl,3,4) = cos(dx)*cos(dy)*cos(dz)
        momfac(isite,ieo,ibl,4,1) = sin(2.0_KR*dx)
        momfac(isite,ieo,ibl,4,2) = sin(2.0_KR*dy)
        momfac(isite,ieo,ibl,4,3) = sin(2.0_KR*dz)
        momfac(isite,ieo,ibl,4,4) = (cos(2.0_KR*dx)+cos(2.0_KR*dy) &
                                    +cos(2.0_KR*dz))/3.0_KR
       enddo ! ixbit
      enddo ! iybit
     enddo ! izbit
    enddo ! itbit
    momfac(:,:,:,:,5) = momfac(:,:,:,:,4)

 end subroutine momfacs

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 
 subroutine spacesum(opsout,opsine,opsino,ieo,ibleo,momfac,nmom,nop,myid, &
                     vecbl,vecblinv,MRT)
! Compute the observables per timeslice.
 
    real(kind=KR),    intent(inout), dimension(:,:,:,:)   :: opsout
    real(kind=KR),    intent(in),    dimension(:,:,:)     :: opsine, opsino
    integer(kind=KI), intent(in)                          :: ieo, ibleo, nmom, &
                                                             nop, myid, MRT
    real(kind=KR),    intent(in),    dimension(:,:,:,:,:) :: momfac
    integer(kind=KI), intent(in),    dimension(:,:)       :: vecbl,vecblinv
 
    integer(kind=KI) :: isite, ibl, jbl, jbleo, nsite, sitemin, &
                        sitemax, it, itbit, ierr, iop, imom, n, iri, itconstant
    integer(kind=KI), dimension(4)             :: np, ip
    real(kind=KR2),   dimension(2,nt,nmom,nop) :: opsbit
    real(kind=KR),    dimension(2,nt,nmom,nop) :: opscorr, opssum
 
! NOTE: it=1 is uniquely defined by          1<=isite<=  factor and 1<=ibl<= 8
!       it=2 is uniquely defined by          1<=isite<=  factor and 9<=ibl<=16
!       it=3 is uniquely defined by   factor+1<=isite<=2*factor and 1<=ibl<= 8
!       it=4 is uniquely defined by   factor+1<=isite<=2*factor and 9<=ibl<=16
!       it=5 is uniquely defined by 2*factor+1<=isite<=3*factor and 1<=ibl<= 8
!       it=6 is uniquely defined by 2*factor+1<=isite<=3*factor and 9<=ibl<=16
!       etc, where factor = 2*nvhalf*npt/nt.

! Sum the data according to timestep.
    ibl = vecbl(1,ibleo)
    jbl = vecbl(2,ibleo)
    jbleo = vecblinv(2,jbl)
    opsbit = 0.0_KR
    np(1) = npx
    np(2) = npy
    np(3) = npz
    np(4) = npt
    call atoc(myid,np,ip)
    itconstant = ip(4)*nt/npt - 1 + (ibleo-1)/4
    nsite = 2*nvhalf*npt/nt
    sitemin = 1
    sitemax = nsite
    do itbit = 2,nt/npt,2
     it = itbit + itconstant
     do isite = sitemin,sitemax
      do iop = 1,nop
       do iri = 1,2
        opsbit(iri,it,1,iop) = opsbit(iri,it,1,iop) &
                             + real(opsine(iri,isite,iop),KR2) &
                             + real(opsino(iri,isite,iop),KR2)
       enddo ! iri
      enddo ! iop
      do imom = 2,nmom
       opsbit(1,it,imom,1) = opsbit(1,it,imom,1)                         &
            + real(momfac(isite,ieo,ibl,imom-1,3)*opsine(2,isite,1),KR2) &
            + real(momfac(isite,ieo,jbl,imom-1,3)*opsino(2,isite,1),KR2)
       opsbit(2,it,imom,1) = opsbit(2,it,imom,1)                         &
            + real(momfac(isite,ieo,ibl,imom-1,2)*opsine(2,isite,1),KR2) &
            + real(momfac(isite,ieo,jbl,imom-1,2)*opsino(2,isite,1),KR2)
       opsbit(1,it,imom,2) = opsbit(1,it,imom,2)                         &
            + real(momfac(isite,ieo,ibl,imom-1,1)*opsine(2,isite,2),KR2) &
            + real(momfac(isite,ieo,jbl,imom-1,1)*opsino(2,isite,2),KR2)
       opsbit(2,it,imom,2) = opsbit(2,it,imom,2)                         &
            + real(momfac(isite,ieo,ibl,imom-1,3)*opsine(2,isite,2),KR2) &
            + real(momfac(isite,ieo,jbl,imom-1,3)*opsino(2,isite,2),KR2)
       opsbit(1,it,imom,3) = opsbit(1,it,imom,3)                         &
            + real(momfac(isite,ieo,ibl,imom-1,2)*opsine(2,isite,3),KR2) &
            + real(momfac(isite,ieo,jbl,imom-1,2)*opsino(2,isite,3),KR2)
       opsbit(2,it,imom,3) = opsbit(2,it,imom,3)                         &
            + real(momfac(isite,ieo,ibl,imom-1,1)*opsine(2,isite,3),KR2) &
            + real(momfac(isite,ieo,jbl,imom-1,1)*opsino(2,isite,3),KR2)
       do iri = 1,2
        opsbit(iri,it,imom,4) = opsbit(iri,it,imom,4)                      &
            + real(momfac(isite,ieo,ibl,imom-1,4)*opsine(iri,isite,4),KR2) &
            + real(momfac(isite,ieo,jbl,imom-1,4)*opsino(iri,isite,4),KR2)
        opsbit(iri,it,imom,5) = opsbit(iri,it,imom,5)                      &
            + real(momfac(isite,ieo,ibl,imom-1,5)*opsine(iri,isite,5),KR2) &
            + real(momfac(isite,ieo,jbl,imom-1,5)*opsino(iri,isite,5),KR2)
       enddo ! iri
      enddo ! imom
     enddo ! isite
     sitemin = sitemin + nsite
     sitemax = sitemax + nsite
    enddo ! itbit

! Collect the final results into opsout.
    opscorr = real(opsbit,KR)
    if (nps==1) then
     do iop = 1,nop
      do imom = 1,nmom
       do it = 1,nt
        do iri = 1,2
         opsout(iri,it,imom,iop) = opsout(iri,it,imom,iop) &
                                + opscorr(iri,it,imom,iop)
        enddo ! iri
       enddo ! it
      enddo ! imom
     enddo ! iop
    else
     n = 2*nt*nmom*nop
     call MPI_REDUCE(opscorr(1,1,1,1),opssum(1,1,1,1),n,MRT,MPI_SUM,0, &
                     MPI_COMM_WORLD,ierr)
     call MPI_BCAST(opssum(1,1,1,1),n,MRT,0,MPI_COMM_WORLD,ierr)
     do iop = 1,nop
      do imom = 1,nmom
       do it = 1,nt
        do iri = 1,2
         opsout(iri,it,imom,iop) = opsout(iri,it,imom,iop) &
                                 + opssum(iri,it,imom,iop)
        enddo ! iri
       enddo ! it
      enddo ! imom
     enddo ! iop
    endif
 
 end subroutine spacesum

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 subroutine checknonzero(xe,xo)

    real(kind=KR), intent(in), dimension(6,ntotal,4,2,8)     :: xe, xo
    integer(kind=KI) :: i,j,k,l,m

      print *, "Checking for ones"
      do i =1,5,2
         do j =1,nvhalf
            do k=1,4
               do l=1,2
                  do m=1,8
                     if(abs(xe(i,j,k,l,m)) /= 0.0_KR) then
                       print "(a32,5i3,1es17.10)", "xe(i,j,k,l,m)=", i,j,k,l,m,xe(i,j,k,l,m)
                     endif
!                    if(abs(xo(i,j,k,l,m)) /= 0.0_KR) then
!                      print *, "bo"
!                      print *, "icsrc,idsrc,j,k=", icsrc,idsrc,j,k
!                    endif
                  end do ! m
               end do ! l
           enddo ! k
        enddo ! j
      enddo ! i
   end subroutine checknonzero

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 
 subroutine fakegauge(uout,myid,rwdir,MRT)
   real(kind=KR),    intent(out),       dimension(18,ntotal,4,2,16)      :: uout
   integer(kind=KI), intent(in)                                          :: myid,MRT
   character(len=*), intent(in),        dimension(:)                     :: rwdir

   real(kind=KR),                       dimension(18,ntotal,4,2,16)      :: udebug
   integer(kind=KI)                                                      :: iblock,ieo,ixyz,j,i,&
                                                                            kc1,kc2,isite,site,inps
   integer(kind=KI)                                                      :: kc,proc,count,ierr
   integer(kind=KI),                    dimension(2)                     :: ix
   integer(kind=KI)                                                      :: iy,iz,it
   integer(kind=KI) :: ieo1,ieo2,itbit,izbit,iybit,ixbit,itbit2,&
                       izbit2,iybit2,ixbit2,ixbit3,iblbit
   integer(kind=KI), dimension(4)                                        :: ip,np
   real(kind=KR),    dimension(nxyzt)                                    :: gcos,gsin

   real(kind=KR)     :: twopi
   integer(kind=KI)  :: icount,ixx,iyy,izz,itt
   logical           :: nondiag, true,false


     uout   = 0.0_KR
     udebug = 0.0_KR
     twopi=8.0_KR*atan(1.0_KR)
    
! nondiag = (.true.) allows  non-diagonal elements for debugging.
  nondiag = .false.

     icount = 0
     do ixx =1,nx
        do iyy=1,ny
           do izz=1,nz
              do itt=1,nt
                 icount=icount+1
                 gcos(icount) = cos((twopi*(ixx+iyy+izz+itt)**1.5)/4.0_KR)
                 gsin(icount) = sin((twopi*(ixx+iyy+izz+itt)**1.5)/4.0_KR)
              enddo ! it
           enddo ! iz
        enddo ! iy
     enddo !ix
     if (myid==0) print *, "icount=",icount

     np(1) = npx
     np(2) = npy
     np(3) = npz
     np(4) = npt
     call atoc(myid,np,ip)

! Begin main loop, constructing ix,iy,iz,it from isite,ieo,ibl.
     isite = 0
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
              isite = isite + 1
              do ieo = 1,2
                 do iblock = 1,16
                    if (iblock>8) then
                        it = itbit2
                    else
                        it = itbit2 - 1
                    endif
                    iblbit = 1 + modulo(iblock-1,8)
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
                    if (modulo(iblock,2)==1) then
                        ixbit3 = ixbit2 - 1
                    else
                        ixbit3 = ixbit2
                    endif
                    ix(ieo1) = ixbit3 - 2
                    ix(ieo2) = ixbit3

                    site = ix(ieo) + (iy-1)*nx + (iz-1)*nx*ny + (it-1)*nx*ny*nz
!                   print "(i4,i10,i4,i4)", site,isite,ieo,iblock

                    do ixyz = 1,3
                       if (nondiag .eqv. .true.) then
                           udebug(3,isite,ixyz,ieo,iblock)    = 1.0_KR
                           udebug(5,isite,ixyz,ieo,iblock)    = 1.0_KR
                           udebug(7,isite,ixyz,ieo,iblock)    = 1.0_KR
                           udebug(11,isite,ixyz,ieo,iblock)   = 1.0_KR
                           udebug(13,isite,ixyz,ieo,iblock)   = 1.0_KR
                           udebug(15,isite,ixyz,ieo,iblock)   = 1.0_KR
                       endif ! nondiag

                       udebug(1,isite,ixyz,ieo,iblock)  = gcos(site)
                       if (ixyz==1 .and. ieo==1. .and. iblock==1) then
                          print *, "FOO: u=", udebug(1,isite,ixyz,ieo,iblock)
                       endif ! stuff
                       udebug(9,isite,ixyz,ieo,iblock)  = gcos(site)
                       udebug(17,isite,ixyz,ieo,iblock) = gcos(site)

                       udebug(2,isite,ixyz,ieo,iblock)  = gsin(site)
                       udebug(10,isite,ixyz,ieo,iblock) = gsin(site)
                       udebug(18,isite,ixyz,ieo,iblock) = gsin(site)

                    enddo ! ixyz
! Do the time componet with one's
                    if (nondiag .eqv. .true.) then
                        udebug(3,isite,4,ieo,iblock)    = 1.0_KR
                        udebug(5,isite,4,ieo,iblock)    = 1.0_KR
                        udebug(7,isite,4,ieo,iblock)    = 1.0_KR
                        udebug(11,isite,4,ieo,iblock)   = 1.0_KR
                        udebug(13,isite,4,ieo,iblock)   = 1.0_KR
                        udebug(15,isite,4,ieo,iblock)   = 1.0_KR
                    endif ! nondiag

                    udebug(1,isite,4,ieo,iblock)  = 1.0_KR
                    udebug(9,isite,4,ieo,iblock)  = 1.0_KR
                    udebug(17,isite,4,ieo,iblock) = 1.0_KR

                  enddo ! iblock
               enddo ! ieo
            enddo ! ixbit
          enddo ! iybit
        enddo ! izbit
      enddo ! itbit


      uout = udebug

 end subroutine fakegauge

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 end module disconloops

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
















