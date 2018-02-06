!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
! inverters.f90, Randy Lewis, randy.lewis@uregina.ca
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
! This module performs the fermion matrix inversion.
! Useful references: 
!
! "Accelerating Wilson fermion matrix inversions by means of the
! stabilized biconjugate gradient algorithm"
! Frommer, Hannemann, Nockel, Lippert & Schilling, Int.J.Mod.Phys.C5,1073(1994)
!
! "Progress on lattice QCD algorithms"
! Ph. de Forcrand, Nucl.Phys.Proc.Suppl.47:228(1996)
!
! "Krylov space solvers for shifted linear systems"
! B. Jegerlehner, hep-lat/9612014.
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -     

 module working
!   use MPI
    use kinds
    use latdims
    use basics
    use diracops
    use pseudolapack
    use gmresrhs
    use inverters
    implicit none
    private

! Use the following line if the MPI module is not available.
    include 'mpif.h'

! Define access to subroutines.
    public  ::  gmresdrEIGBS
    private :: qrfactorizationLAPACK, leastsquaresLAPACK,eigmodesLAPACK,matrixmultiplylike1 

 contains

!!!!!!!




 subroutine gmresdrEIGBS(rwdir,b,x,GMRES,rtol,itermin,itercount,u,GeeGooinv, &
                    iflag,idag,kappa,coact,bc,vecbl,vecblinv,myid,nn,ldiv,nms, &
                    lvbc,ib,lbd,iblv,MRT,MRT2)
    use shift

    character(len=*), intent(in),    dimension(:)         :: rwdir
    real(kind=KR),    intent(in),    dimension(:,:,:,:,:) :: b
    real(kind=KR),    intent(inout), dimension(:,:,:,:,:) :: x
    integer(kind=KI), intent(in),    dimension(:)         :: GMRES, bc, nms
    integer(kind=KI), intent(in)                          :: itermin, iflag, &
                                                             myid, MRT, MRT2, idag
    real(kind=KR),    intent(in)                          :: rtol
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

    integer(kind=KI) :: icycle, i, j,  jp1, jj, ir, irm1, is, ivb, ii, &
                        idis, icri, m, k,kk, nrhs, ilo, ihi, ischur, &
                        id, ieo, ibleo
    logical,          dimension(nmaxGMRES)                  :: myselect
    integer(kind=KI), dimension(nmaxGMRES)                  :: ipiv,ind
    real(kind=KR2)                                          :: const, tval, &
                                                               amags, con2, rv, &
                                                               normnum
    real(kind=KR2),   dimension(2)                          :: beta, tv1, tv2, &
                                                               tv,rninit,rn,vn, vnf,xrn,rina
    real(kind=KR2),   dimension(nmaxGMRES)                  :: mag,dabbs
    real(kind=KR2),   dimension(2,nmaxGMRES)                :: ss, gs, gc, &
                                                               tau, w
    real(kind=KR),       dimension(2,201,201)                     :: jamuna!BS
    real(kind=KR2),   dimension(2,nmaxGMRES+1)              :: c, c2, srv,d
    real(kind=KR2),   dimension(2,nmaxGMRES,1)              :: ff, punty !Chris
    real(kind=KR2),   dimension(2,nmaxGMRES)                :: dd,th,tmpVec,resi
    real(kind=KR2),   dimension(2,kmaxGMRES)                :: rho
    real(kind=KR2),   dimension(kmaxGMRES)                  :: rna,sita
    real(kind=KR2),   dimension(2,nmaxGMRES+1,nmaxGMRES+1)  :: z,gon,rr,gondag
    real(kind=KR2),   dimension(2,nmaxGMRES,nmaxGMRES+1)    :: hcht
    real(kind=KR2),   dimension(2,nmaxGMRES+1,nmaxGMRES)    :: h, h2, h3, hprint, hh,g,gg,greal, &
                                                               tmpmat,hnew
    real(kind=KR2),   dimension(2,nmaxGMRES+1,kmaxGMRES)    :: ws
    real(kind=KR2),   dimension(2,kmaxGMRES+1,kmaxGMRES)    :: gca, gsa
    real(kind=KR2),   dimension(6,ntotal,4,2,8)             :: r, xt
!   real(kind=KR2),   dimension(6,ntotal,4,2,8)             :: xb
    !real(kind=KR2),   dimension(6,ntotal,4,2,8,nmaxGMRES+1) :: v
    !real(kind=KR2),   dimension(6,nvhalf,4,2,8,nmaxGMRES+1) :: vtemp
    !real(kind=KR2),   dimension(6,ntotal,4,2,8)             :: xb
    !real(kind=KR2),   dimension(2,nmaxGMRES+1,nmaxGMRES)    :: hcnew
    real(kind=KR2),   dimension(6,ntotal,4,2,8,nmaxGMRES+1) :: work
    real(kind=KR2),   dimension(6,ntotal,4,2,8) :: tmpE, tmpE2
    real(kind=KR2),   dimension(6,ntotal,4,2,8) :: tmpO, tmpO2
    real(kind=KR), dimension(6,ntotal,4,2,8,nmaxGMRES+1)  ::muna,tina,matmul,matmully
    real(kind=KR), dimension(6,ntotal,4,2,8)  ::getemp,wv,f,xtmp,mina 
   ! real(kind=KR2),   dimension(6,nmaxGMRES+1)              :: vt
    real(kind=KR2),   dimension(6)              :: vt
    real(kind=KR2), dimension(2) :: tmp1,tmp2,tmp3,tmp4,broda
    real(kind=KR2)   :: temp11,temp22,temp33,normmatmully
    real(kind=KR), dimension(6,ntotal,4,2,8)   :: htemp
    integer(kind=KI) :: iblock, isite, idirac,icolorir, site, icolorr
    integer(kind=KI) :: didmaindo,pickle
    integer(kind=KI) :: exists,vena  ! initialize variables
  m = GMRES(1)  ! size of Krylov Subspace
  k = GMRES(2)  ! number of eigenvector/values to deflate
 ! print *,'k=',k
  !idag = 0      ! flag to do M*x NOT Mdag*x

  ! ignore initial guess and use 0
  x(:6,:ntotal,:4,:2,:8) = 0.0_KR2

  icycle = 1    ! initialize cycle counter
  j = 1         ! intiialize iteration counter
  h = 0.0_KR2   ! initialize h matrix

  ! calculate inital residual by ignoring intial guess and use x=0
  ! r = Ax - b = A*0 - b = b
  r = b
  call vecdot(r,r,rninit,MRT2)
  rninit(1) = sqrt(rninit(1))
  rninit(2) = 0.0_KR2
  rn = rninit
  vn = rn

  ! first vector
  v(:6,:ntotal,:4,:2,:8,1) = (1.0_KR2 / vn(1)) * r(:6,:ntotal,:4,:2,:8)

  ! right hand side to future least square problem
  c = 0.0_KR2
  c(1,1) = vn(1)

  if (myid==0) then
     print *,"gmresdrprinttest1"
  endif
  ! begin cycles
  do while(((rn(1)/rninit(1)) > 1e-160_KR2) .AND. (icycle <= 150))!becomes zero in 29 cycles
 ! do while(icycle <= 39)!BS 1/22/2016 to make sure enough evalue pass tolerance
    ! begin gmres(m) iterations
    do while(j <= m)
      call Hdbletm(wv,u,GeeGooinv,v(:,:,:,:,:,j),idag,coact,kappa,iflag,bc,vecbl,vecblinv, &
                   myid,nn,ldiv,nms,lvbc,ib,lbd,iblv,MRT)

      call vecdot(wv,wv,vnf,MRT2)
      vnf(1) = sqrt(vnf(1))
      vnf(2) = 0.0_KR2

      do i=1,j
        call vecdot(v(:,:,:,:,:,i),wv,h(:,i,j),MRT2)
        do icri = 1,5,2 ! 6=nri*nc
          do jj = 1,nvhalf
            wv(icri  ,jj,:,:,:) = wv(icri  ,jj,:,:,:) &
                                  - h(1,i,j)*v(icri  ,jj,:,:,:,i) &
                                  + h(2,i,j)*v(icri+1,jj,:,:,:,i)
            wv(icri+1,jj,:,:,:) = wv(icri+1,jj,:,:,:) &
                                  - h(2,i,j)*v(icri  ,jj,:,:,:,i) &
                                  - h(1,i,j)*v(icri+1,jj,:,:,:,i)
          enddo
        enddo
      enddo

      call vecdot(wv,wv,vn,MRT2)
      vn(1) = sqrt(vn(1))
      vn(2) = 0.0_KR2

      ! --- reorthogonalization section ---
      if (vn(1) < (1.1_KR * vnf(1))) then
        do i=1,j
          call vecdot(v(:,:,:,:,:,i),wv,tmp1,MRT2)
          do icri = 1,5,2 ! 6=nri*nc
            do jj = 1,nvhalf
              wv(icri  ,jj,:,:,:) = wv(icri  ,jj,:,:,:) &
                                  - tmp1(1)*v(icri  ,jj,:,:,:,i) &
                                  + tmp1(2)*v(icri+1,jj,:,:,:,i)
              wv(icri+1,jj,:,:,:) = wv(icri+1,jj,:,:,:) &
                                  - tmp1(2)*v(icri  ,jj,:,:,:,i) &
                                  - tmp1(1)*v(icri+1,jj,:,:,:,i)
            enddo
          enddo
          h(:,i,j) = h(:,i,j) + tmp1(:)
        enddo
        call vecdot(wv,wv,vn,MRT2)
        vn(1) = sqrt(vn(1))
        vn(2) = 0.0_KR2
      endif
      ! --- --- --- --- --- --- --- --- ---

      h(:,j+1,j) = vn
      v(:6,:ntotal,:4,:2,:8,j+1) = (1.0_KR2 / h(1,j+1,j)) * wv(:6,:ntotal,:4,:2,:8)

      j = j + 1
    enddo


!------------------************************----------------------------

!BS 5/12/2016 print jamuna to check orthogonality condition for all
!icycles,creates orthonormal.dat !look at scratch file,remember to set nps=1
!though

!   call vecdothigherrank (muna,muna,jamuna,m+1,MRT2)! vecdothigherrank was
!   supposed to give jamuna replacing the algorithms below
 
if (.false.) then
  do i=1,m+1
   do j=1,m+1

      broda=0.0_KR2

     call vecdot(v(:,:,:,:,:,i),v(:,:,:,:,:,j),broda,MRT2)

     jamuna(1,i,j)=broda(1)
     jamuna(2,i,j)=broda(2)

   enddo!j
  enddo!i


        if (myid==0) then
           inquire(file=trim(rwdir(myid+1))//"orthonormal.dat", exist=exists)
            if (.not. exists) then

               open(unit=47,file=trim(rwdir(myid+1))//"orthonormal.dat",status="new",&
               action="write",form="formatted")
               close(unit=47,status="keep")
            endif
       endif


     do i=1,m+1
       do j=1,m+1
            if (myid==0) then
                 open(unit=47,file=trim(rwdir(myid+1))//"orthonormal.dat",status="old",action="write",&
                 form="formatted",position="append")
                 write(unit=47,fmt="(a9,i7,a10,i3,a2,i3,a3,es22.12,es22.12)") "icycle=",icycle,"V'V(",i,",",j,")=",jamuna(1,i,j),jamuna(2,i,j)
                 close(unit=47,status="keep")
            endif
      enddo
    enddo




endif!true or false



!BS 05/12/016 print of jamuna completed succesfully....change (.false.) to (.true.) if want to create
!orthonormal.dat

!------------------------*******************************---------------------------















!-------------------********************************************--------------------------

!BS 05/12/2016 creates orthogonali.dat to see ||AV(1:m)-V(1:(m+1)*h(m+1,m)||

   if (.false.) then
 !    if (icycle>=28) then

          muna(:6,:nvhalf,:4,:2,:8,:(m+1)) =  v(:6,:nvhalf,:4,:2,:8,:(m+1))
          do i=1,m
             call Hdbletm(wv,u,GeeGooinv,v(:,:,:,:,:,i),idag,coact,kappa,iflag,bc,vecbl,vecblinv,&
                          myid,nn,ldiv,nms,lvbc,ib,lbd,iblv,MRT)
             tina(:,:,:,:,:,i) = wv(:,:,:,:,:)
          enddo

          call matrixmultiplylike1(muna,h,m,matmul,MRT2)




if(.false.) then
if (myid ==0) then
     print *,"tina(6,128,4,2,8,5) =", tina(6,128,4,2,8,5)
     print *,"tina(6,128,4,2,7,5) =", tina(6,128,4,2,7,5)
     print *,"muna(6,128,4,2,6,5) =", muna(6,128,4,2,6,5)
     print *,"v(6,128,4,2,5,5) =", v(6,128,4,2,5,5)
    ! print *,"tina(6,128,4,2,4,5) =", tina(6,128,4,2,4,5)
    ! print *,"tina(6,128,4,2,3,5) =", tina(6,128,4,2,3,5)
    ! print *,"tina(6,18,4,2,3,5) =",  tina(6,18,4,2,3,5)
     print *,"matmul(6,128,4,2,8,5) =", matmul(6,128,4,2,8,5)
     print *,"matmul(6,128,4,2,7,5) =", matmul(6,128,4,2,7,5)
    ! print *,"matmul(6,128,4,2,6,5) =", matmul(6,128,4,2,6,5)
    ! print *,"matmul(6,128,4,2,5,5) =", matmul(6,128,4,2,5,5)
    ! print *,"matmul(6,128,4,2,4,5) =", matmul(6,128,4,2,4,5)
    ! print *,"matmul(6,128,4,2,3,5) =", matmul(6,128,4,2,3,5)
    ! print *,"matmul(6,18,4,2,3,5) =",  matmul(6,18,4,2,3,5)
endif
endif !true or false


         normmatmully=0.0_KR2

          do ibleo = 1,8
           do ieo = 1,2
            do id = 1,4
             do isite = 1,nvhalf
              do icri =1,6
               do j = 1,m

           
               matmully(icri,isite,id,ieo,ibleo,j) = matmul(icri,isite,id,ieo,ibleo,j)  - tina(icri,isite,id,ieo,ibleo,j)
               normmatmully = normmatmully + matmully(icri,isite,id,ieo,ibleo,j)**2
               enddo!j
              enddo!icri
             enddo!isite
            enddo!id
           enddo!ieo
          enddo!ibleo

        if (myid==0) then
           print *,"normmatmully = ",normmatmully
        endif   

        if (myid==0) then
            inquire(file=trim(rwdir(myid+1))//"orthogonali.dat", exist=exists)
            if (.not. exists) then

               open(unit=48,file=trim(rwdir(myid+1))//"orthogonali.dat",status="new",&
               action="write",form="formatted")
               close(unit=48,status="keep")
            endif
        endif



        if (myid==0) then
            open(unit=48,file=trim(rwdir(myid+1))//"orthogonali.dat",status="old",action="write",&
                 form="formatted",position="append")
                do ibleo = 1,8
                 do ieo = 1,2
                  do id = 1,4
                   do isite = 1,nvhalf
                    do icri =1,5,2
                     do j = 1,m

                        write(unit=48,fmt="(i6,i6,i6,i6,i6,i6,i6,es22.12,es22.12)") icycle,&
                              icri,isite,id,ieo,ibleo,j,matmully(icri,isite,id,ieo,ibleo,m),&
                               matmully(icri+1,isite,id,ieo,ibleo,m)
                     enddo!j
                    enddo!icri
                   enddo!isite
                  enddo!id
                 enddo!ieo
                enddo!ibleo
            close(unit=48,status="keep")
        endif! myid


!     endif!if cycle==28
   endif! true or false




!BS 05/12/016 end creating orthogonali.dat

!----------------------------*************************-------------------------







    call leastsquaresLAPACK(h,m+1,m,c,d,myid)
    call matvecmult(h,m+1,m,d,m,srv)
    srv(:,:m+1) = c(:,:m+1) - srv(:,:m+1)

    ! x(:) = x(:) + v(:,1:m)*d(1:m)
    do i=1,m
      do icri = 1,5,2 ! 6=nri*nc
        do jj = 1,nvhalf
          x(icri  ,jj,:,:,:) = x(icri  ,jj,:,:,:) &
                                + d(1,i)*v(icri  ,jj,:,:,:,i) &
                                - d(2,i)*v(icri+1,jj,:,:,:,i)
          x(icri+1,jj,:,:,:) = x(icri+1,jj,:,:,:) &
                                + d(2,i)*v(icri  ,jj,:,:,:,i) &
                                + d(1,i)*v(icri+1,jj,:,:,:,i)
        enddo
      enddo
    enddo





!BS added this part for  linear equation residual 




      call Hdbletm(wv,u,GeeGooinv,x(:,:,:,:,:),idag,coact,kappa,iflag,bc,vecbl,vecblinv,&
                   myid,nn,ldiv,nms,lvbc,ib,lbd,iblv,MRT)


      xtmp = wv - b

     
      call vecdot(xtmp,xtmp,xrn,MRT2)
            xrn(1)=sqrt(xrn(1))
            xrn(2)=0.0_KR2

      


!BS residual for linear equation ends here











    ! r = v(:,1:m+1)*srv(1:m+1)
    r = 0.0_KR2
    do icri = 1,5,2 ! 6=nri*nc
      do jj = 1,nvhalf
        do i=1,m+1
          r(icri  ,jj,:,:,:) = r(icri,jj,:,:,:)   &
                              + srv(1,i)*v(icri  ,jj,:,:,:,i) &
                              - srv(2,i)*v(icri+1,jj,:,:,:,i)
          r(icri+1,jj,:,:,:) = r(icri+1,jj,:,:,:)  &
                              + srv(2,i)*v(icri  ,jj,:,:,:,i) &
                              + srv(1,i)*v(icri+1,jj,:,:,:,i)
        enddo
      enddo
    enddo

    call vecdot(r,r,rn,MRT2)
    rn(1) = sqrt(rn(1))
    rn(2) = 0.0_KR2

    ! prepare for next cycle
    hh(:2,1:m,1:m) = h(:2,1:m,1:m)
    do i=1,m
      do jj=1,m
        hcht(1,i,jj) = h(1,jj,i)
        hcht(2,i,jj) = -1.0_KR2 * h(2,jj,i)
      enddo
    enddo

    ff=0.0_KR2
    ff(1,m,1)=1.0_KR2

    call linearsolver(m,1,hcht,ipiv,ff)

    do i=1,m
      hh(1,i,m) = hh(1,i,m)-h(2,m+1,m)**2*ff(1,i,1)-2.0_KR2*h(2,m+1,m)*h(1,m+1,m)*ff(2,i,1)+h(1,m+1,m)**2*ff(1,i,1)
      hh(2,i,m) = hh(2,i,m)-h(2,m+1,m)**2*ff(2,i,1)+2.0_KR2*h(2,m+1,m)*h(1,m+1,m)*ff(1,i,1)+h(1,m+1,m)**2*ff(2,i,1)
    enddo

    ! sorted from smallest to biggest eigenpairs [th,g]
    call eigencalc(hh,m,1,dd,gg)!BS 1/19/2016
    !call eigmodesLAPACK(hh,m,1,dd,gg)

    do i=1,m
      dabbs(i) = dd(1,i)**2+dd(2,i)**2
    enddo

    call sort(dabbs,m,ind)


    do i=1,k
      th(:,i) = dd(:,ind(i))
      g(:,:,i) = gg(:,:,ind(i))
    enddo


   ! Compute Residual Norm of Eigenvectors of hh matrix
    do i=1,k
      call matvecmult(h,m,m,g(:,:,i),m,tmpVec)
      call vecdagvec(g(:,:,i),m,tmpVec,m,rho(:,i))

      do jj=1,m
        call cmpxmult(rho(:,i),g(:,jj,i),tmp1)
        tmpVec(:,jj) = tmpVec(:,jj) - tmp1
      enddo
      call vecdagvec(tmpVec,m,tmpVec,m,tmp1)
      tmp1(1) = sqrt(tmp1(1))
      tmp1(2) = 0.0_KR2

      rna(i) = sqrt((tmp1(1)*tmp1(1)) + (((h(1,m+1,m)*h(1,m+1,m)) + (h(2,m+1,m)*h(2,m+1,m))) &
                * ((g(1,m,i)*g(1,m,i)) + (g(2,m,i)*g(2,m,i)))))!BS changed + to *



    enddo!i



    call sort(rna,k,ind)

    do i=1,k
       sita(i)=rna(ind(i))!BS 5/4/2016
    enddo



    do i=1,k
       rna(i)=sita(i)!BS 5/4/2016
    enddo



if (.false.) then

    do i=1,k  !BS 5/4/2016

        if (myid==0) then
           inquire(file=trim(rwdir(myid+1))//"residual.dat", exist=exists)
            if (.not. exists) then

               open(unit=33,file=trim(rwdir(myid+1))//"residual.dat",status="new",&
               action="write",form="formatted")
               close(unit=33,status="keep")
            endif
       endif




             if (myid==0) then
            open(unit=33,file=trim(rwdir(myid+1))//"residual.dat",status="old",action="write",&
            form="formatted",position="append")
                 write(unit=33,fmt="(i7,a6,i7,a6,es19.12)") icycle,"   ",i," ",rna(i)
           close(unit=33,status="keep")
          endif

   enddo !i

endif !true or false















    do i=1,k
      gg(:,:,i) = g(:,:,ind(i))
    enddo

    do i=1,k
      gg(:,m+1,i) = 0.0_KR2
    enddo

! Chris

    beta(1) = h(1,m+1,m)
    beta(2) = h(2,m+1,m)

    do j = 1,m
        punty(1,j,1) = -beta(1)*ff(1,j,1)+beta(2)*ff(2,j,1)
        punty(2,j,1) = -beta(2)*ff(1,j,1)-beta(1)*ff(2,j,1)
    enddo

    do j = 1,m
        gg(:,j,k+1) = punty(:,j,1)
    enddo
    gg(1,m+1,k+1) = 1.0_KR2

!end Chris

!    gg(:,:,k+1) = srv !-----Chris

    call qrfactorizationLAPACK(gg,m+1,k+1,gon,rr,myid)

    do i=1,m+1
      do jj=1,m+1
        gondag(1,i,jj) = gon(1,jj,i)
        gondag(2,i,jj) =-1.0_KR2 * gon(2,jj,i)
      enddo
    enddo

    ! hcnew = gon'*h*gon(1:m,1:k) 
    call matmult(gondag,k+1,m+1,h,m+1,m,tmpmat)
    call matmult(tmpmat,k+1,m,gon,m,k,hcnew)

    h(:,k+1,:m) = 0.0_KR2

    i = 1
    do while (rna(i) < 1E-12)
      hnew(:,i+1:k+1,i) = 0.0_KR2
      i = i + 1
    enddo

    ! form right eigenvectors; evector is in shift module
    do j=1,k
      do ibleo = 1,8
        do ieo = 1,2
          do id = 1,4
            !do i = 1,ntotal
            do i = 1,nvhalf
              vt(:) = 0.0_KR2
              do kk=1,m
                do icri = 1,5,2
                  vt(icri)   = vt(icri) &
                        + v(icri  ,i,id,ieo,ibleo,kk)*gg(1,kk,j) &
                        - v(icri+1,i,id,ieo,ibleo,kk)*gg(2,kk,j) 
                  vt(icri+1) = vt(icri+1) &
                        + v(icri  ,i,id,ieo,ibleo,kk)*gg(2,kk,j) &
                        + v(icri+1,i,id,ieo,ibleo,kk)*gg(1,kk,j) 
                enddo
              enddo
              evector(:,i,id,ieo,ibleo,j) = vt(:)
            enddo
          enddo
        enddo
      enddo
      evalue(:,j) = th(:,ind(j))
! if (myid==0) then
 !              write(*,*) 'evalueEIG:',j,'j:',evalue
  !          endif

    enddo
       !BS do while(j <= m)
          !BS if (myid==0) then
            !BS   write(*,*) 'evalueEIG:',j,'j:',evalue
         !BS  endif
       !BS enddo

    do i=1,k
      do jj=1,k+1
        h(:,jj,i) = hcnew(:,jj,i)
      enddo
    enddo

    call matvecmult(gondag,k+1,m+1,srv,m+1,c)

    do i=k+2,m+1
      c(:,i) = 0.0_KR2
    enddo

    do jj=1,k+1
      do ibleo = 1,8
        do ieo = 1,2
          do id = 1,4
            !do i = 1,ntotal
            do i = 1,nvhalf
              vt(:) = 0.0_KR2
              do kk=1,m+1
                do icri = 1,5,2
                  vt(icri)   = vt(icri) &
                      + v(icri  ,i,id,ieo,ibleo,kk)*gon(1,kk,jj) &
                      - v(icri+1,i,id,ieo,ibleo,kk)*gon(2,kk,jj) 
                  vt(icri+1) = vt(icri+1) &
                      + v(icri  ,i,id,ieo,ibleo,kk)*gon(2,kk,jj) &
                      + v(icri+1,i,id,ieo,ibleo,kk)*gon(1,kk,jj) 
                enddo
              enddo
              work(:,i,id,ieo,ibleo,jj) = vt(:)
            enddo
          enddo
        enddo
      enddo
    enddo
    do i=1,k+1
     ! v(:6,:nvhalf,:4,:2,:8,i) = work(:6,:nvhalf,:4,:2,:8,i)
      v(:,:,:,:,:,i) = work(:,:,:,:,:,i)
    enddo

    do i=1,k
      call vecdot(v(:,:,:,:,:,i),v(:,:,:,:,:,k+1),tmp1,MRT2)

      do icri = 1,5,2 ! 6=nri*nc
        !do jj = 1,ntotal
        do jj = 1,nvhalf
          v(icri  ,jj,:,:,:,k+1) = v(icri  ,jj,:,:,:,k+1) &
                             - tmp1(1)*v(icri  ,jj,:,:,:,i) &
                             + tmp1(2)*v(icri+1,jj,:,:,:,i)
          v(icri+1,jj,:,:,:,k+1) = v(icri+1,jj,:,:,:,k+1) &
                             - tmp1(2)*v(icri  ,jj,:,:,:,i) &
                             - tmp1(1)*v(icri+1,jj,:,:,:,i)
        enddo
      enddo
    enddo

    call vecdot(v(:,:,:,:,:,k+1),v(:,:,:,:,:,k+1),tmp1,MRT2)
    tmp1(1) = sqrt(tmp1(1))

    !do jj = 1,ntotal
    do jj = 1,nvhalf
      v(:,jj,:,:,:,k+1) = (1.0_KR2 / tmp1(1))*v(:,jj,:,:,:,k+1)
    enddo

!BS8888888888888888888*****************
if (.false.) then
!muna =0.0_KR2
!tina=0.0_KR2
 muna(:6,:ntotal,:4,:2,:8,:(k+1)) =  v(:6,:ntotal,:4,:2,:8,:(k+1))
          do i=1,k
             call Hdbletm(wv,u,GeeGooinv,v(:,:,:,:,:,i),idag,coact,kappa,iflag,bc,vecbl,vecblinv,&
                          myid,nn,ldiv,nms,lvbc,ib,lbd,iblv,MRT)
             tina(:,:,:,:,:,i) = wv(:,:,:,:,:)
          enddo
if (myid == 0) then
print *, "tina(6,123,4,1,8,1)=",tina(6,123,4,1,8,1)
print *, "muna(6,123,4,1,8,1)=",muna(6,123,4,1,8,1)
endif
          call matrixmultiplylike1(muna,h,k,matmul,MRT2)
if (myid == 0) then
print *, "matmul(6,123,4,1,8,1)=",matmul(6,123,4,1,8,1)
endif

        normmatmully=0.0_KR2

          do ibleo = 1,8
           do ieo = 1,2
            do id = 1,4
             do isite = 1,nvhalf
              do icri =1,6
               do j = 1,k


               matmully(icri,isite,id,ieo,ibleo,j) = matmul(icri,isite,id,ieo,ibleo,j)  - tina(icri,isite,id,ieo,ibleo,j)
               normmatmully = normmatmully + matmully(icri,isite,id,ieo,ibleo,j)**2


               enddo!j
              enddo!icri
             enddo!isite
            enddo!id
           enddo!ieo
          enddo!ibleo

if (myid==0) then
print *,matmully(6,123,4,2,8,34)

endif
        if (myid==0) then
           print *,"normmatmully2 = ",normmatmully
        endif

endif !true or false

















!if (myid==0) then
!  write(*,*) 'cycle:',icycle,'resnorm:',rn(1)/rninit(1)
!endif

      if (myid==0) then
       open(unit=8,file=trim(rwdir(myid+1))//"CFGSPROPS.LOG",action="write", form="formatted",status="old",position="append")
!      BS  write(unit=8,fmt="(a12,i9,es17.10)") "gmresdr",itercount,beta(1)
       !BS write(unit=8,fmt="(a9,i4,a3,a5,es10.6,es10.6,es17.10)") "gmresdr",icycle,"rn","rnit",rn(1),rninit(1),rn(1)/rninit(1)
                write(unit=8,fmt="(a9,i4,a4,es11.4,a6,es11.4,es17.10)") "gmrsEIG",icycle,"rn=",rn(1),"rnit=",rninit(1),rn(1)/rninit(1)
                    write(unit=8,fmt="(a9,i4,a6,es11.4,a6,es11.4,es17.10)")"gmrsEIG",icycle,"xrn=",xrn(1),"xrnit=",rninit(1),xrn(1)/rninit(1)

     !  write(unit=8,fmt="(a12,i9,es17.10)") "gmresdr",icycle,rn(1)/rninit(1)
       close(unit=8,status="keep")
      endif



do vena=1,k


!BS added this part for true residual



      call Hdbletm(wv,u,GeeGooinv,evector(:,:,:,:,:,vena),idag,coact,kappa,iflag,bc,vecbl,vecblinv,&
                   myid,nn,ldiv,nms,lvbc,ib,lbd,iblv,MRT)

do ibleo = 1,8
        do ieo = 1,2
          do id = 1,4
            !do i = 1,ntotal
            do i = 1,nvhalf
              do icri =1,5,2 


    mina(icri,i,id,ieo,ibleo) = evalue(1,vena)*evector(icri,i,id,ieo,ibleo,vena)&
                                   - evalue(2,vena)*evector(icri+1,i,id,ieo,ibleo,vena)
    mina(icri+1,i,id,ieo,ibleo) = evalue(1,vena)*evector(icri+1,i,id,ieo,ibleo,vena)&
                                   + evalue(2,vena)*evector(icri,i,id,ieo,ibleo,vena)

              enddo
            enddo
         enddo
     enddo
enddo



      xtmp = wv - mina


      call vecdot(xtmp,xtmp,xrn,MRT2)
            rina(1)=sqrt(xrn(1))
            rina(2)=0.0_KR2


 if (.false.)then


        if (myid==0) then
            inquire(file=trim(rwdir(myid+1))//"trueresidual.dat", exist=exists)
            if (.not. exists) then

               print *, "File does not exist. Creating it."
               open(unit=46,file=trim(rwdir(myid+1))//"trueresidual.dat",status="new",&
               action="write",form="formatted")
               close(unit=46,status="keep")
            endif
       endif


             if (myid==0) then
            open(unit=46,file=trim(rwdir(myid+1))//"trueresidual.dat",status="old",action="write",&
            form="formatted",position="append")
                 write(unit=46,fmt="(i5,a5,i5,a5,es20.12)") icycle,"   ",vena," ",rina(1)
           close(unit=46,status="keep")
          endif

 endif!true or false 
enddo






    j = k+1
    icycle = icycle + 1
  enddo
        
        !  if (myid==0) then !BS
              
              !BS   write(*,fmt="(a12,f19.11)") "evalueEIG:",evalue
          !    print *,'evalueEIG:',evalue
         ! endif
     if (myid==0) then
    open(unit=8,file=trim(rwdir(myid+1))//"CFGSPROPS.LOG",action="write", &
         form="formatted",status="old",position="append")
   !BS  write(unit=8,fmt="(a12,i9,es17.10)")
   !"--gmresdr",icycle-1,rn(1)/rninit(1)
                write(unit=8,fmt="(a9,i4,a4,es11.4,a6,es11.4,es17.10)") "gmrsEIG",icycle-1,"rn=",rn(1),"rnit=",rninit(1),rn(1)/rninit(1)

    close(unit=8,status="keep")
  endif
  


  if (myid==0) then
    do i =1,k 
    open(unit=8,file=trim(rwdir(myid+1))//"CFGSPROPS.LOG",action="write", &
         form="formatted",status="old",position="append")
   !BS  write(unit=8,fmt="(a12,i9,es17.10)") "--gmresdr",icycle-1,rn(1)/rninit(1)
                write(unit=8,fmt="(i9,F17.12,F17.12)") i,evalue(1,i),evalue(2,i)

    close(unit=8,status="keep")
    enddo  
  endif
end subroutine gmresdrEIGBS



 subroutine qrfactorizationLAPACK(a,m,n,qfact,rfact,myid)
 integer(kind=KI), intent(in)                    :: m, n,myid
 real(kind=KR2),   intent(in),  dimension(:,:,:) :: a
 real(kind=KR2),   intent(out), dimension(:,:,:) :: qfact
 real(kind=KR2),   intent(out), dimension(:,:,:) :: rfact

 complex*16, allocatable, dimension(:,:) :: a_lap
 complex*16, allocatable, dimension(:)   :: tau
 complex*16, allocatable, dimension(:)   :: work
 integer(kind=KI) :: lda, info, lwork

 integer(kind=KI) :: i,j,k

 allocate(a_lap(m,n))
 allocate(tau(min(m,n)))
 allocate(work(10*n))

 lda = m
 lwork = 10*n
 k = min(m,n)

 do i=1,m
   do j=1,n
     a_lap(i,j) = cmplx(a(1,i,j),a(2,i,j),KR2)
   enddo
 enddo

 call zgeqrf(m,n,a_lap,lda,tau,work,lwork,info)
 rfact = 0.0_KR2
 do i=1,m
   do j=i,n
     rfact(1,i,j) = real(a_lap(i,j),KR2)
     rfact(2,i,j) = aimag(a_lap(i,j))
   enddo
 enddo

 call zungqr(m,n,k,a_lap,lda,tau,work,lwork,info)

 qfact = 0.0_KR2
 do i=1,m
   do j=1,n
     qfact(1,i,j) = real(a_lap(i,j),KR2)
     qfact(2,i,j) = aimag(a_lap(i,j))
   enddo
 enddo

 deallocate(a_lap)
 deallocate(tau)
 deallocate(work)
 end subroutine qrfactorizationLAPACK


 subroutine leastsquaresLAPACK(a,m,n,b,x,myid)
 real(kind=KR2),   intent(in), dimension(:,:,:) :: a
 real(kind=KR2),   intent(in), dimension(:,:)   :: b
 real(kind=KR2),   intent(out), dimension(:,:)  :: x
 integer(kind=KI), intent(in)                   :: m,n,myid

 complex*16, allocatable, dimension(:,:) :: a_lap
 complex*16, allocatable, dimension(:) :: work
 complex*16, allocatable, dimension(:,:) :: b_lap

! complex*16, dimension(m,n) :: a_lap
! complex*16, dimension(10*m*n) :: work
! complex*16, dimension(m,1) :: b_lap

 CHARACTER*1 trans
 integer(kind=KI) info, lda, ldb, lwork, nrhs

 integer(kind=KI) i,j,k



 allocate(a_lap(m,n))
 allocate(work(10*m*n))
 allocate(b_lap(m,1))


 trans = 'N'
 nrhs = 1
 lda = m
 ldb = m
 lwork = 10*m*n
 do i=1,m
   do j=1,n
     a_lap(i,j) = cmplx(a(1,i,j),a(2,i,j),KR2)
   enddo
 enddo

 do i=1,m
   b_lap(i,1) = cmplx(b(1,i),b(2,i),KR2)
 enddo

 call zgels(trans,m,n,nrhs,a_lap,lda,b_lap,ldb,work,lwork,info)

 do i=1,n
   x(1,i) = real(b_lap(i,1),KR2)
   x(2,i) = aimag(b_lap(i,1))
 enddo

 deallocate(a_lap)
 deallocate(work)
 deallocate(b_lap)
 end subroutine leastsquaresLAPACK

subroutine eigmodesLAPACK(mat,m,mode,eval,evec)
 real(kind=KR2),   intent(in), dimension(:,:,:)  :: mat
 integer(kind=KI), intent(in)                    :: m
 integer(kind=KI), intent(in)                    :: mode
 real(kind=KR2),   intent(out), dimension(:,:)   :: eval
 real(kind=KR2),   intent(out), dimension(:,:,:) :: evec

 ! local LAPACK variables
 complex*16, allocatable, dimension(:,:) :: a
 complex*16, allocatable, dimension(:,:) :: vl
 complex*16, allocatable, dimension(:,:) :: vr
 complex*16, allocatable, dimension(:)   :: w
 complex*16, allocatable, dimension(:)   :: work
!complex(kind=KR2), allocatable, dimension(:,:) :: a!BS changed comples to KR2
!complex(kind=KR2), allocatable, dimension(:,:) :: vl
!complex(kind=KR2), allocatable, dimension(:,:) :: vr
!complex(kind=KR2), allocatable, dimension(:)   :: w
! complex(kind=KR2), allocatable, dimension(:)   :: work
 ! real*16,    allocatable, dimension(:)   :: rwork
 real(kind=KR2),    allocatable, dimension(:)   :: rwork

 CHARACTER*1 jobvl
 CHARACTER*1 jobvr

 integer(kind=KI) info, lda, ldvl, ldvr, lwork, n

 real(kind=KR2), parameter :: DOUBLEPRECISIONLIMIT = 0 !1e-14_KR2

 ! loop variables
 integer(kind=KI) i,j,k
 ! set LAPACK integer values
 n     = m
 lda   = m
 ldvl  = m
 ldvr  = m
 lwork = 10*n ! notes give as >= max(1,2*n)
 jobvl = 'N'
 if (mode == 1) then
   jobvr = 'V'! else
   jobvr = 'N'
 endif

 ! allocate appropriate space
 allocate(a(lda,n))
 allocate(vl(ldvl,n))
 allocate(vr(ldvr,n))
 allocate(w(n))
 allocate(work(lwork))
 allocate(rwork(2*n))
 ! copy into LAPACK type arrays
 do i=1,lda
   do j=1,n
       a(i,j) = cmplx(mat(1,i,j),mat(2,i,j),KR2)
          if (abs(mat(1,i,j)) <= DOUBLEPRECISIONLIMIT .and. abs(mat(2,i,j)) >=DOUBLEPRECISIONLIMIT) then
       a(i,j) = cmplx(0.0,mat(2,i,j),KR2)
     else if (abs(mat(1,i,j)) >= DOUBLEPRECISIONLIMIT .and. abs(mat(2,i,j)) <= DOUBLEPRECISIONLIMIT) then
       a(i,j) = cmplx(mat(1,i,j),0.0,KR2)
     else if (abs(mat(1,i,j)) >= DOUBLEPRECISIONLIMIT .and. abs(mat(2,i,j)) >= DOUBLEPRECISIONLIMIT) then
       a(i,j) = cmplx(mat(1,i,j),mat(2,i,j),KR2)
     else if (abs(mat(1,i,j)) <= DOUBLEPRECISIONLIMIT .and. abs(mat(2,i,j)) <= DOUBLEPRECISIONLIMIT) then
       a(i,j) = cmplx(0.0,0.0,KR2)
     else
       write(*,*) 'AHHHHHHHHHHHHH!!!!!!!!!!!',DOUBLEPRECISIONLIMIT
     endif
   enddo
 enddo

    print * ,'about to call zgeev','  mode=',mode
 ! call LAPACK routine
 call zgeev(jobvl,jobvr,n,a,lda,w,vl,ldvl,vr,ldvr,work,lwork,rwork,info)

    print * ,'just  called zgeev','mode=',mode
 !   print * , 'w=',w,'vr =',vr,'info',info
 ! error checking
 if (info > 0) then
   write(*,*) 'ERROR: QR algorithm failed - zgeev'
   stop
 endif

 ! move evalues into final storage position
 do i=1,n
   eval(1,i) = real(w(i),KR2)
   eval(2,i) = aimag(w(i))
 enddo

 ! move evectors into final storage position
 do i=1,ldvr
   do j=1,n
     evec(1,i,j) = real(vr(i,j),KR2)
     evec(2,i,j) = aimag(vr(i,j))
   enddo
 enddo

 deallocate(a)
 deallocate(vl)
 deallocate(vr)
 deallocate(w)
 deallocate(work)
 deallocate(rwork)
 end subroutine eigmodesLAPACK










  subroutine matrixmultiplylike1(a,h,m,matmul,MRT2)
! Calculate the dot product of two vectors, a and b, where the dot
! product
! is understood to mean    sum_i a^dagger(i)*b(i) .
! Define the result to be  adb(1)+i*abd(2), where i=sqrt(-1).
! MRT2 is MPIDBLETYPE.


    real(kind=KR),    intent(in),  dimension(:,:,:,:,:,:) :: a
    real(kind=KR),    intent(in),  dimension(:,:,:)       :: h
    real(kind=KR),    intent(out),  dimension(:,:,:,:,:,:) :: matmul
    integer(kind=KI), intent(in)                          :: MRT2,m

    integer(kind=KI)             :: isite, id, ieo, ibleo,ierr,j,Aham,icri
    real(kind=KR2), dimension(2) :: adbbit
    real(kind=KR2)               :: temp1,temp2

    do ibleo = 1,8
     do ieo = 1,2
      do id = 1,4
       do isite = 1,nvhalf
        do icri =1,5,2
         do j = 1,m
         temp1=0.0_KR2
         temp2=0.0_KR2
          do Aham = 1,m+1

            matmul(icri,isite,id,ieo,ibleo,j)   =temp1             &
                                                +a(icri,isite,id,ieo,ibleo,Aham)*h(1,Aham,j) &
                                                -a(icri+1,isite,id,ieo,ibleo,Aham)*h(2,Aham,j)
            matmul(icri+1,isite,id,ieo,ibleo,j) =temp2             &
                                                + a(icri,isite,id,ieo,ibleo,Aham)*h(2,Aham,j)   &
                                                + a(icri+1,isite,id,ieo,ibleo,Aham)*h(1,Aham,j)

            temp1=matmul(icri,isite,id,ieo,ibleo,j)
            temp2=matmul(icri+1,isite,id,ieo,ibleo,j)


          enddo !Aham
           ! print *,matmul(icri,isite,id,ieo,ibleo,j)
           ! print *,matmul(icri+1,isite,id,ieo,ibleo,j)

         enddo !j
        enddo !icri
       enddo  !isite
      enddo !id
     enddo  !ieo
    enddo   !ibleo


            print *,"inside matrixmultiply matmul(6,128,4,2,8,5)",matmul(6,128,4,2,8,5)
            print *,"inside matrixmultiply matmul(6,128,4,2,7,5)",matmul(6,128,4,2,7,5)

! Sum the contributions from all processes.
    if (nps==1) then
    print*,"matrixmultiplylike1 is working"
    else
print *,"Error in subroutine matrrixmultiplylike , nps not equal to one"

!     call
!     MPI_REDUCE(adbbit(1),adb(1),2,MRT2,MPI_SUM,0,MPI_COMM_WORLD,ierr)
 !    call MPI_BCAST(adb(1),2,MRT2,0,MPI_COMM_WORLD,ierr)
    endif

 end subroutine matrixmultiplylike1










































 end module working 
