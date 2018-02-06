!------------------------------------------------------------------------------ 
 subroutine ppaverage1(Jtime,nsub,nmom,nop,momfac,u,be,bo,xe,xo,kappa,dobndry, &
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
                        ibleo, ieo, idag, i, j, p, iflag
    real(kind=KR)                               :: xk, fac1, fac2
    real(kind=KR),    dimension(6,ntotal,4,2,8) :: sub1e, sub1o, sub2e, sub2o,&
                                                   vecsrc, vecsrc1, temp, &
                                                   sub3e, sub3o
                                                                                    real(kind=KR),    dimension(6,nvhalf,4,2,8) :: sube, subo, GeeGooinv
    integer(kind=KI), parameter                 :: itstep=2*nvhalf*npt/nt
    real(kind=KR),    dimension(18,itstep,2,8)  :: ubndry
    !real(kind=KR),    dimension(2,nvhalf,5)     :: Je, Jo !5 should be replaced by nop
    real(kind=KR),    dimension(2,nvhalf,nop)     :: Je, Jo !5 should be replaced by nop
    integer(kind=KI) :: ierr, it, iop
    integer(kind=KI) :: icolor, idirac, mu, ibl1, gblclr1, jbl, jbleo, jeo
    real(kind=KR2),   dimension(2)                          :: beta
    ! This local gaugelink is used for debugging
    real(kind=KR),  dimension(18,ntotal,4,2,16) :: uout
    integer(kind=KI) :: exists!BS 3/2/2016
    character(len=128) :: startingvector !BS 3/2/2016
    ! Initialization.
    p = 12
    idag = 0
    Jtime = 0.0_KR
    GeeGooinv = 0.0_KR
    iflag = -1
    vprime = 0.0_KR2
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
      gblclr = 1
      idag = 0
      call Hsingle(temp,u,bo,idag,coact,bc,gblclr,vecbl, &
                   vecblinv,myid,nn,ldiv,nms,lvbc,ib,lbd,iblv,MRT)

      fac1 = 1.0_KR/kappa**2
      fac2 = 1.0_KR/kappa
      do isite = 1,nvhalf
       vecsrc(:,isite,:,:,:) = fac1*be(:,isite,:,:,:) + fac2*temp(:,isite,:,:,:)
      enddo ! isite
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!define the right hand side for xo
      gblclr = 2
      idag = 0
      call Hsingle(temp,u,be,idag,coact,bc,gblclr,vecbl, &
                   vecblinv,myid,nn,ldiv,nms,lvbc,ib,lbd,iblv,MRT)

      fac1 = 1.0_KR/kappa**2
      fac2 = 1.0_KR/kappa
      do isite = 1,nvhalf
       vecsrc1(:,isite,:,:,:)=fac1*bo(:,isite,:,:,:) + fac2*temp(:,isite,:,:,:)
      enddo ! isite
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

write(startingvector,"(a6,i3)")"startingvector"

if (myid==0) then
        inquire(file=trim(rwdir(myid+1))//trim(startingvector)//".LOG", exist=exists)
  if (.not. exists) then

         print *, "File does not exist. Creating it inside noise file."
open(unit=80,file=trim(rwdir(myid+1))//trim(startingvector)//".LOG",status="new",&
        action="write",form="formatted")
   close(unit=80,status="keep")
end if

end if

 call printferm(vecsrc,vecsrc1,rwdir,startingvector,myid,iblv,vecblinv,MRT2)!BS 3/2/2016
!Determine the polynomial for xe.!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    do isite = 1,nvhalf
     vprime(:,isite,:,:,:,1) = vecsrc(:,isite,:,:,:)
    enddo !k
    do i = 1,p
     call Hdouble(vprime(:,:,:,:,:,i+1),u,vprime(:,:,:,:,:,i),idag, &
                 coact,kappa,bc,vecbl,vecblinv,myid,nn, &
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
     call vecdot(vprime(:,:,:,:,:,i),vecsrc(:,:,:,:,:),beta,MRT2)
     cls(:,i-1,1) = beta(:)
!     print *, "i,cls(:,i)=", i-1, cls(:,i-1,1)
   enddo!i
    
    call linearsolver(p,1,lsmat,ipiv2,cls)
    co(:,:) = cls(:,:,1)    
!    co = 0.0_KR2    
!    co(1,1) = 4
   if(myid==0) then
    print * , "coefficients  following  Should be same as matlab"
    do i=1,p
     print *, "i,result(:,i)=", i, co(:,i)
    enddo!i  
   endif!myid
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Determine the polynomial for xo.!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    vprime = 0.0_KR2
    do isite = 1,nvhalf
     vprime(:,isite,:,:,:,1) = vecsrc1(:,isite,:,:,:)
    enddo !k
    do i = 1,p
     call Hdouble1(vprime(:,:,:,:,:,i+1),u,vprime(:,:,:,:,:,i),idag, &
                 coact,kappa,bc,vecbl,vecblinv,myid,nn, &
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
     call vecdot(vprime(:,:,:,:,:,i),vecsrc1(:,:,:,:,:),beta,MRT2)
     cls(:,i-1,1) = beta(:)
!     print *, "i,cls(:,i)=", i-1, cls(:,i-1,1)
   enddo!i
    
    call linearsolver(p,1,lsmat,ipiv2,cls)
    co1(:,:) = cls(:,:,1)    
!    co = 0.0_KR2    
!    co(1,1) = 4
   if(myid==0) then
    do i=1,p
     print *, "i,result1(:,i)=", i, co1(:,i)
    enddo!i  
   endif!myid
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    ! Compute the operators for each level of subtraction.Andy is here.In the pp
    ! average,we use P(M) as the M^(-1)_pert.The P(M) is generated when we do th    ! ppgmresdr invertion.The even and odd part are calculated seperately.

    do isub = 1,nsub
      select case(isub)
        case(1) ! no subtraction
          ksub = 0
          sube = 0.0_KR
          subo = 0.0_KR
        case(2) ! subtraction of O(kappa), O(kappa^2) and O(kappa^3)
        !even part    
          ksub = 1
          do icri=1,5,2
           do isite=1,nvhalf
            sube(icri,isite,:,:,:) = co(1,1)*vecsrc(icri,isite,:,:,:) &
                                      -co(2,1)*vecsrc(icri+1,isite,:,:,:)
            sube(icri+1,isite,:,:,:) = co(1,1)*vecsrc(icri+1,isite,:,:,:) &
                                        +co(2,1)*vecsrc(icri,isite,:,:,:)
           enddo!isite
          enddo!icri
          call Hdouble(sub1e,u,vecsrc,idag,coact,kappa,bc, &
                       vecbl,vecblinv,myid,nn,ldiv,nms,lvbc,ib,lbd,iblv,MRT)
        do icri=1,5,2
         do isite=1,nvhalf
          sube(icri  ,isite,:,:,:) = sube(icri ,isite,:,:,:) &
                                     +co(1,2)*sub1e(icri,isite,:,:,:) &
                                     -co(2,2)*sub1e(icri+1,isite,:,:,:)
          sube(icri+1,isite,:,:,:) = sube(icri+1,isite,:,:,:) &
                                     +co(1,2)*sub1e(icri+1,isite,:,:,:) &
                                     +co(2,2)*sub1e(icri,isite,:,:,:)
         enddo!isite
        enddo!icri             

        !odd part
          do icri=1,5,2
           do isite=1,nvhalf
            subo(icri,isite,:,:,:) = co1(1,1)*vecsrc1(icri,isite,:,:,:) &
                                      -co1(2,1)*vecsrc1(icri+1,isite,:,:,:)
            subo(icri+1,isite,:,:,:) = co1(1,1)*vecsrc1(icri+1,isite,:,:,:) &
                                        +co1(2,1)*vecsrc1(icri,isite,:,:,:)
           enddo!isite
          enddo!icri
          call Hdouble1(sub1o,u,vecsrc1,idag,coact,kappa,bc, &
                       vecbl,vecblinv,myid,nn,ldiv,nms,lvbc,ib,lbd,iblv,MRT)
        do icri=1,5,2
         do isite=1,nvhalf
          subo(icri  ,isite,:,:,:) = subo(icri ,isite,:,:,:) &
                                     +co1(1,2)*sub1o(icri,isite,:,:,:) &
                                     -co1(2,2)*sub1o(icri+1,isite,:,:,:)
          subo(icri+1,isite,:,:,:) = subo(icri+1,isite,:,:,:) &
                                     +co1(1,2)*sub1o(icri+1,isite,:,:,:) &
                                     +co1(2,2)*sub1o(icri,isite,:,:,:)
         enddo!isite
        enddo!icri             
!!!!!!!!!!!!!!!!!!!!!!!order1!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !even part
        call Hdouble(sub2e,u,sub1e,idag,coact,kappa,bc, &
                   vecbl,vecblinv,myid,nn,ldiv,nms,lvbc,ib,lbd,iblv,MRT)
        do icri=1,5,2
         do isite=1,nvhalf
          sube(icri  ,isite,:,:,:) = sube(icri ,isite,:,:,:) &
                                     +co(1,3)*sub2e(icri,isite,:,:,:) &
                                     -co(2,3)*sub2e(icri+1,isite,:,:,:)
          sube(icri+1,isite,:,:,:) = sube(icri+1,isite,:,:,:) &
                                     +co(1,3)*sub2e(icri+1,isite,:,:,:) &
                                     +co(2,3)*sub2e(icri,isite,:,:,:)
         enddo!isite
        enddo!icri             
         !odd part
         call Hdouble1(sub2o,u,sub1o,idag,coact,kappa,bc, &
                   vecbl,vecblinv,myid,nn,ldiv,nms,lvbc,ib,lbd,iblv,MRT)
        do icri=1,5,2
         do isite=1,nvhalf
          subo(icri  ,isite,:,:,:) = subo(icri ,isite,:,:,:) &
                                     +co1(1,3)*sub2o(icri,isite,:,:,:) &
                                     -co1(2,3)*sub2o(icri+1,isite,:,:,:)
          subo(icri+1,isite,:,:,:) = subo(icri+1,isite,:,:,:) &
                                     +co1(1,3)*sub2o(icri+1,isite,:,:,:) &
                                     +co1(2,3)*sub2o(icri,isite,:,:,:)
         enddo!isite
        enddo!icri             
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!order2!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !even part           
         call Hdouble(sub1e,u,sub2e,idag,coact,kappa,bc, &
                   vecbl,vecblinv,myid,nn,ldiv,nms,lvbc,ib,lbd,iblv,MRT)
        do icri=1,5,2
         do isite=1,nvhalf
          sube(icri  ,isite,:,:,:) = sube(icri ,isite,:,:,:) &
                                     +co(1,4)*sub1e(icri,isite,:,:,:) &
                                     -co(2,4)*sub1e(icri+1,isite,:,:,:)
          sube(icri+1,isite,:,:,:) = sube(icri+1,isite,:,:,:) &
                                     +co(1,4)*sub1e(icri+1,isite,:,:,:) &
                                     +co(2,4)*sub1e(icri,isite,:,:,:)
         enddo!isite
        enddo!icri             
         !odd part
        call Hdouble1(sub1o,u,sub2o,idag,coact,kappa,bc, &
                   vecbl,vecblinv,myid,nn,ldiv,nms,lvbc,ib,lbd,iblv,MRT)
        do icri=1,5,2
         do isite=1,nvhalf
          subo(icri  ,isite,:,:,:) = subo(icri ,isite,:,:,:) &
                                     +co1(1,4)*sub1o(icri,isite,:,:,:) &
                                     -co1(2,4)*sub1o(icri+1,isite,:,:,:)
          subo(icri+1,isite,:,:,:) = subo(icri+1,isite,:,:,:) &
                                     +co1(1,4)*sub1o(icri+1,isite,:,:,:) &
                                     +co1(2,4)*sub1o(icri,isite,:,:,:)
         enddo!isite
        enddo!icri             
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!order3!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!                                
        case(3) !order4,order5,order6 
         !even part
         call Hdouble(sub2e,u,sub1e,idag,coact,kappa,bc, &
                   vecbl,vecblinv,myid,nn,ldiv,nms,lvbc,ib,lbd,iblv,MRT)
        do icri=1,5,2
         do isite=1,nvhalf
          sube(icri  ,isite,:,:,:) = sube(icri ,isite,:,:,:) &
                                     +co(1,5)*sub2e(icri,isite,:,:,:) &
                                     -co(2,5)*sub2e(icri+1,isite,:,:,:)
          sube(icri+1,isite,:,:,:) = sube(icri+1,isite,:,:,:) &
                                     +co(1,5)*sub2e(icri+1,isite,:,:,:) &
                                     +co(2,5)*sub2e(icri,isite,:,:,:)
         enddo!isite
        enddo!icri             
         call Hdouble(sub1e,u,sub2e,idag,coact,kappa,bc, &
                   vecbl,vecblinv,myid,nn,ldiv,nms,lvbc,ib,lbd,iblv,MRT)
        do icri=1,5,2
         do isite=1,nvhalf
          sube(icri  ,isite,:,:,:) = sube(icri ,isite,:,:,:) &
                                     +co(1,6)*sub1e(icri,isite,:,:,:) &
                                     -co(2,6)*sub1e(icri+1,isite,:,:,:)
          sube(icri+1,isite,:,:,:) = sube(icri+1,isite,:,:,:) &
                                     +co(1,6)*sub1e(icri+1,isite,:,:,:) &
                                     +co(2,6)*sub1e(icri,isite,:,:,:)
         enddo!isite
        enddo!icri             

         call Hdouble(sub2e,u,sub1e,idag,coact,kappa,bc, &
                   vecbl,vecblinv,myid,nn,ldiv,nms,lvbc,ib,lbd,iblv,MRT)
        do icri=1,5,2
         do isite=1,nvhalf
          sube(icri  ,isite,:,:,:) = sube(icri ,isite,:,:,:) &
                                     +co(1,7)*sub2e(icri,isite,:,:,:) &
                                     -co(2,7)*sub2e(icri+1,isite,:,:,:)
          sube(icri+1,isite,:,:,:) = sube(icri+1,isite,:,:,:) &
                                     +co(1,7)*sub2e(icri+1,isite,:,:,:) &
                                     +co(2,7)*sub2e(icri,isite,:,:,:)
         enddo!isite
        enddo!icri             

          !odd part
         call Hdouble1(sub2o,u,sub1o,idag,coact,kappa,bc, &
                   vecbl,vecblinv,myid,nn,ldiv,nms,lvbc,ib,lbd,iblv,MRT)
        do icri=1,5,2
         do isite=1,nvhalf
          subo(icri  ,isite,:,:,:) = subo(icri ,isite,:,:,:) &
                                     +co1(1,5)*sub2o(icri,isite,:,:,:) &
                                     -co1(2,5)*sub2o(icri+1,isite,:,:,:)
          subo(icri+1,isite,:,:,:) = subo(icri+1,isite,:,:,:) &
                                     +co1(1,5)*sub2o(icri+1,isite,:,:,:) &
                                     +co1(2,5)*sub2o(icri,isite,:,:,:)
         enddo!isite
        enddo!icri             

         call Hdouble1(sub1o,u,sub2o,idag,coact,kappa,bc, &
                   vecbl,vecblinv,myid,nn,ldiv,nms,lvbc,ib,lbd,iblv,MRT)
        do icri=1,5,2
         do isite=1,nvhalf
          subo(icri  ,isite,:,:,:) = subo(icri ,isite,:,:,:) &
                                     +co1(1,6)*sub1o(icri,isite,:,:,:) &
                                     -co1(2,6)*sub1o(icri+1,isite,:,:,:)
          subo(icri+1,isite,:,:,:) = subo(icri+1,isite,:,:,:) &
                                     +co1(1,6)*sub1o(icri+1,isite,:,:,:) &
                                     +co1(2,6)*sub1o(icri,isite,:,:,:)
         enddo!isite
        enddo!icri             

         call Hdouble1(sub2o,u,sub1o,idag,coact,kappa,bc, &
                   vecbl,vecblinv,myid,nn,ldiv,nms,lvbc,ib,lbd,iblv,MRT)
        do icri=1,5,2
         do isite=1,nvhalf
          subo(icri  ,isite,:,:,:) = subo(icri ,isite,:,:,:) &
                                     +co1(1,7)*sub2o(icri,isite,:,:,:) &
                                     -co1(2,7)*sub2o(icri+1,isite,:,:,:)
          subo(icri+1,isite,:,:,:) = subo(icri+1,isite,:,:,:) &
                                     +co1(1,7)*sub2o(icri+1,isite,:,:,:) &
                                     +co1(2,7)*sub2o(icri,isite,:,:,:)
         enddo!isite
        enddo!icri             

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!order4,5,6!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        case(4)!subtraction of order7,order8,order9 
          !even part
         call Hdouble(sub1e,u,sub2e,idag,coact,kappa,bc, &
                   vecbl,vecblinv,myid,nn,ldiv,nms,lvbc,ib,lbd,iblv,MRT)
        do icri=1,5,2
         do isite=1,nvhalf
          sube(icri  ,isite,:,:,:) = sube(icri ,isite,:,:,:) &
                                     +co(1,8)*sub1e(icri,isite,:,:,:) &
                                     -co(2,8)*sub1e(icri+1,isite,:,:,:)
          sube(icri+1,isite,:,:,:) = sube(icri+1,isite,:,:,:) &
                                     +co(1,8)*sub1e(icri+1,isite,:,:,:) &
                                     +co(2,8)*sub1e(icri,isite,:,:,:)
         enddo!isite
        enddo!icri             
         call Hdouble(sub2e,u,sub1e,idag,coact,kappa,bc, &
                   vecbl,vecblinv,myid,nn,ldiv,nms,lvbc,ib,lbd,iblv,MRT)
        do icri=1,5,2
         do isite=1,nvhalf
          sube(icri  ,isite,:,:,:) = sube(icri ,isite,:,:,:) &
                                     +co(1,9)*sub2e(icri,isite,:,:,:) &
                                     -co(2,9)*sub2e(icri+1,isite,:,:,:)
          sube(icri+1,isite,:,:,:) = sube(icri+1,isite,:,:,:) &
                                     +co(1,9)*sub2e(icri+1,isite,:,:,:) &
                                     +co(2,9)*sub2e(icri,isite,:,:,:)
         enddo!isite
        enddo!icri             

         call Hdouble(sub1e,u,sub2e,idag,coact,kappa,bc, &
                   vecbl,vecblinv,myid,nn,ldiv,nms,lvbc,ib,lbd,iblv,MRT)
        do icri=1,5,2
         do isite=1,nvhalf
          sube(icri  ,isite,:,:,:) = sube(icri ,isite,:,:,:) &
                                     +co(1,10)*sub1e(icri,isite,:,:,:) &
                                     -co(2,10)*sub1e(icri+1,isite,:,:,:)
          sube(icri+1,isite,:,:,:) = sube(icri+1,isite,:,:,:) &
                                     +co(1,10)*sub1e(icri+1,isite,:,:,:) &
                                     +co(2,10)*sub1e(icri,isite,:,:,:)
         enddo!isite
        enddo!icri             

         !odd part
         call Hdouble1(sub1o,u,sub2o,idag,coact,kappa,bc, &
                   vecbl,vecblinv,myid,nn,ldiv,nms,lvbc,ib,lbd,iblv,MRT)
        do icri=1,5,2
         do isite=1,nvhalf
          subo(icri  ,isite,:,:,:) = subo(icri ,isite,:,:,:) &
                                     +co1(1,8)*sub1o(icri,isite,:,:,:) &
                                     -co1(2,8)*sub1o(icri+1,isite,:,:,:)
          subo(icri+1,isite,:,:,:) = subo(icri+1,isite,:,:,:) &
                                     +co1(1,8)*sub1o(icri+1,isite,:,:,:) &
                                     +co1(2,8)*sub1o(icri,isite,:,:,:)
         enddo!isite
        enddo!icri             
         call Hdouble1(sub2o,u,sub1o,idag,coact,kappa,bc, &
                   vecbl,vecblinv,myid,nn,ldiv,nms,lvbc,ib,lbd,iblv,MRT)
        do icri=1,5,2
         do isite=1,nvhalf
          subo(icri  ,isite,:,:,:) = subo(icri ,isite,:,:,:) &
                                     +co1(1,9)*sub2o(icri,isite,:,:,:) &
                                     -co1(2,9)*sub2o(icri+1,isite,:,:,:)
          subo(icri+1,isite,:,:,:) = subo(icri+1,isite,:,:,:) &
                                     +co1(1,9)*sub2o(icri+1,isite,:,:,:) &
                                     +co1(2,9)*sub2o(icri,isite,:,:,:)
         enddo!isite
        enddo!icri             

         call Hdouble1(sub1o,u,sub2o,idag,coact,kappa,bc, &
                   vecbl,vecblinv,myid,nn,ldiv,nms,lvbc,ib,lbd,iblv,MRT)
        do icri=1,5,2
         do isite=1,nvhalf
          subo(icri  ,isite,:,:,:) = subo(icri ,isite,:,:,:) &
                                     +co1(1,10)*sub1o(icri,isite,:,:,:) &
                                     -co1(2,10)*sub1o(icri+1,isite,:,:,:)
          subo(icri+1,isite,:,:,:) = subo(icri+1,isite,:,:,:) &
                                     +co1(1,10)*sub1o(icri+1,isite,:,:,:) &
                                     +co1(2,10)*sub1o(icri,isite,:,:,:)
         enddo!isite
        enddo!icri             

   
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!order7,8,9!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        case(5) 
          !even part
         call Hdouble(sub2e,u,sub1e,idag,coact,kappa,bc, &
                   vecbl,vecblinv,myid,nn,ldiv,nms,lvbc,ib,lbd,iblv,MRT)
        do icri=1,5,2
         do isite=1,nvhalf
          sube(icri  ,isite,:,:,:) = sube(icri ,isite,:,:,:) &
                                     +co(1,11)*sub2e(icri,isite,:,:,:) &
                                     -co(2,11)*sub2e(icri+1,isite,:,:,:)
          sube(icri+1,isite,:,:,:) = sube(icri+1,isite,:,:,:) &
                                     +co(1,11)*sub2e(icri+1,isite,:,:,:) &
                                     +co(2,11)*sub2e(icri,isite,:,:,:)
         enddo!isite
        enddo!icri             
         !odd part
         call Hdouble1(sub2o,u,sub1o,idag,coact,kappa,bc, &
                   vecbl,vecblinv,myid,nn,ldiv,nms,lvbc,ib,lbd,iblv,MRT)
        do icri=1,5,2
         do isite=1,nvhalf
          subo(icri  ,isite,:,:,:) = subo(icri ,isite,:,:,:) &
                                     +co1(1,11)*sub2o(icri,isite,:,:,:) &
                                     -co1(2,11)*sub2o(icri+1,isite,:,:,:)
          sube(icri+1,isite,:,:,:) = sube(icri+1,isite,:,:,:) &
                                     +co1(1,11)*sub2o(icri+1,isite,:,:,:) &
                                     +co1(2,11)*sub2o(icri,isite,:,:,:)
         enddo!isite
        enddo!icri             


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!order10!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        case(6) ! subtraction of order 11
                  !even part
         call Hdouble(sub1e,u,sub2e,idag,coact,kappa,bc, &
                   vecbl,vecblinv,myid,nn,ldiv,nms,lvbc,ib,lbd,iblv,MRT)
        do icri=1,5,2
         do isite=1,nvhalf
          sube(icri  ,isite,:,:,:) = sube(icri ,isite,:,:,:) &
                                     +co(1,12)*sub1e(icri,isite,:,:,:) &
                                     -co(2,12)*sub1e(icri+1,isite,:,:,:)
          sube(icri+1,isite,:,:,:) = sube(icri+1,isite,:,:,:) &
                                     +co(1,12)*sub1e(icri+1,isite,:,:,:) &
                                     +co(2,12)*sub1e(icri,isite,:,:,:)
         enddo!isite
        enddo!icri             
         !odd part
         call Hdouble1(sub1o,u,sub2o,idag,coact,kappa,bc, &
                   vecbl,vecblinv,myid,nn,ldiv,nms,lvbc,ib,lbd,iblv,MRT)
        do icri=1,5,2
         do isite=1,nvhalf
          subo(icri  ,isite,:,:,:) = subo(icri ,isite,:,:,:) &
                                     +co1(1,12)*sub1o(icri,isite,:,:,:) &
                                     -co1(2,12)*sub1o(icri+1,isite,:,:,:)
          subo(icri+1,isite,:,:,:) = subo(icri+1,isite,:,:,:) &
                                     +co1(1,12)*sub1o(icri+1,isite,:,:,:) &
                                     +co1(2,12)*sub1o(icri,isite,:,:,:)
         enddo!isite
        enddo!icri             


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!order11!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         
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
    end subroutine ppaverage1
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

