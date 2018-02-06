subroutine xipolygenerator(xiinvPOLY,numEV,eigstart,eigstep,u,kappa, &
                        evecReven,evecRodd,evecLeven,evecLodd,z2e,z2o, &
                        coact,bc,vecbl,vecblinv,myid,nn,ldiv,nms,lvbc,ib,lbd,
&
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
 real(kind=KR2), dimension(6,ntotal,4,2,8) :: tmpE,tmpO,vprime,vprime1,tue,tuo
 real(kind=KR2), dimension(6,ntotal,4,2,8) ::sube,subo,tempsube,tempsubo
 real(kind=KR2), dimension(6,ntotal,4,2,8) :: sub1e,sub1o,temptempsube
 real(kind=KR2), dimension(6,ntotal,4,2,8) :: sub2e,sub2o,temptempsubo
 integer(kind=KI)   :: isite, iri, ibl, itbit, gblcrl
 integer(kind=KI)   :: isub, ksub, id,nsub, idag
 real(kind=KR2)     :: xk
 integer(kind=KI) :: exists
 integer(kind=KI) ::  idag, p, i, j
 real(kind=KR2),   dimension(2)                          :: beta, beta1
 real(kind=KR2),   dimension(2,8,8)        ::  lsmat,lsmat1
 real(kind=KR2),   dimension(2,8,1)        ::  cls1
 real(kind=KR2),   dimension(2,8)          ::  co2


 idag = 0
 nsub = 6



!copy from ppaverage1 begins now


    p = 8
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

  enddo




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
!    co = 0.0_KR2
!    co(1,1) = 4
   if(myid==0) then
    do i=1,p
     print *, "i,result(:,i)=", i, co2(:,i)
    enddo!i
   endif!myid



 do eig = 1, numEV

  !BS   if (gammaid == 1_KI) gamma5 not done here so its same for both
  !gammaid
      be = evecReven(:,:,:,:,:,eig)
      bo = evecRodd(:,:,:,:,:,eig)








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
                    sube(icri,isite,id,ieo,ibleo) = co2(:,1)*be(icri,isite,id,ieo,ibleo)
&
                                         +co2(:,2)* kappa*sub1e(icri,isite,id,ieo,ibleo)
                    subo(icri,isite,id,ieo,ibleo) =co2(:,1)*  bo(icri,isite,id,ieo,ibleo)
&
                                         +co2(:,2)* kappa*sub1o(icri,isite,id,ieo,ibleo)
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
                    sube(icri,isite,id,ieo,ibleo) =sube(icri,isite,id,ieo,ibleo) &
                                              +co2(:,3)* xk*sub2e(icri,isite,id,ieo,ibleo)
                    subo(icri,isite,id,ieo,ibleo) =subo(icri,isite,id,ieo,ibleo) &
                                              +co2(:,3)* xk*sub2o(icri,isite,id,ieo,ibleo)
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
                    sube(icri,isite,id,ieo,ibleo) =sube(icri,isite,id,ieo,ibleo) &
                                              +co2(:,4)* xk*sub1e(icri,isite,id,ieo,ibleo)
                    subo(icri,isite,id,ieo,ibleo) =subo(icri,isite,id,ieo,ibleo) &
                                              +co2(:,4)* xk*sub1o(icri,isite,id,ieo,ibleo)
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
                    sube(icri,isite,id,ieo,ibleo) =sube(icri,isite,id,ieo,ibleo) &
                                              +co2(:,5)* xk*sub2e(icri,isite,id,ieo,ibleo)
                    subo(icri,isite,id,ieo,ibleo) =subo(icri,isite,id,ieo,ibleo) &
                                              +co2(:,5)* xk*sub2o(icri,isite,id,ieo,ibleo)
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
                    sube(icri,isite,id,ieo,ibleo) =sube(icri,isite,id,ieo,ibleo) &
                                              +co2(:,6)* xk*sub1e(icri,isite,id,ieo,ibleo)
                    subo(icri,isite,id,ieo,ibleo) =subo(icri,isite,id,ieo,ibleo) &
                                              +co2(:,6)* xk*sub1o(icri,isite,id,ieo,ibleo)
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
                  do icri = 1,6
                    sube(icri,isite,id,ieo,ibleo) =sube(icri,isite,id,ieo,ibleo)&
                                              +co2(:,7)* xk*sub2e(icri,isite,id,ieo,ibleo)
                    subo(icri,isite,id,ieo,ibleo) =subo(icri,isite,id,ieo,ibleo)&
                                              +co2(:,7)* xk*sub2o(icri,isite,id,ieo,ibleo)
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
                  do icri = 1,6
                    sube(icri,isite,id,ieo,ibleo) =sube(icri,isite,id,ieo,ibleo) &
                                              +co2(:,8)* xk*sub1e(icri,isite,id,ieo,ibleo)
                    subo(icri,isite,id,ieo,ibleo) =subo(icri,isite,id,ieo,ibleo) &
                                              +co2(:,8)* xk*sub1o(icri,isite,id,ieo,ibleo)
                  enddo ! icri
                enddo ! isite
              enddo ! id
            enddo ! ieo
          enddo ! ibleo


       case default





          open(unit=8,file="DISCONLOOPS.ERROR",action="write",status="replace",
&
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
         if (myid==0) then
            print *,'gammaid',gammaid
         endif
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
      xiinvPOLY(:2,isub, eig) = tmp1(:2) + tmp2(:2)


   sube=temptempsube
   subo=temptempsubo
    enddo
 enddo


 end subroutine xigenerator

