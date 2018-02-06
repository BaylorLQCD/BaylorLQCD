!-------------------------------------------------------------------------
! Utility subroutines to print the Dirac operator matrices or the hopping
! matrices and also sources to be used in MATLAB tetsing.
!
! Added by Abdou: January 27th, 2014
!--------------------------------------------------------------------------


! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
! printops.f90, Abdou Abdel-Rehim, amabdelrehim@gmail.com
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 module printops

!    use MPI
    use kinds
    use latdims
    use basics
    use lattice
    use gaugetools
    use diracops
    use quark
    use inverters

    implicit none
    private

! Use the following line if the MPI module is not available.
   include 'mpif.h'

   public    :: printH,printferm

 contains

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 subroutine printH(rwdir,u,coact,bc,vecbl,vecblinv,myid,nn,ldiv,nms,lvbc,ib,lbd,iblv,MRT,MRT2)

    character(len=*), intent(in),    dimension(:)         :: rwdir
    real(kind=KR),    intent(in),    dimension(:,:,:,:,:) :: u
    integer(kind=KI), intent(in),    dimension(:)         :: bc, nms
    integer(kind=KI), intent(in)                          :: myid, MRT, MRT2
    real(kind=KR),    intent(in),    dimension(:,:,:)     :: coact
    integer(kind=KI), intent(in),    dimension(:,:)       :: vecbl, vecblinv
    integer(kind=KI), intent(in),    dimension(:,:)       :: nn, iblv
    logical,          intent(in),    dimension(:)         :: ldiv
    integer(kind=KI), intent(in),    dimension(:,:,:)     :: lvbc
    integer(kind=KI), intent(in),    dimension(:,:,:,:)   :: ib
    logical,          intent(in),    dimension(:,:)       :: lbd

    real(kind=KR), dimension(6,ntotal,4,2,8)  :: bel,bol,ber,bor,betmp,botmp !right and left sources
    integer(kind=KI), dimension(4)            :: srcl,srcr
    integer(kind=KI)                          :: it,iz,iy,ix,ic,id,is
    integer(kind=KI)                          :: jt,jz,jy,jx,jc,jd,js
    integer(kind=KI)  :: nc,ns
    integer(kind=KI)  ::ieo,jeo
    integer(kind=KI)  ::idag,gblclr
    integer(kind=KI)  :: i,j,k
    real(kind=KR), dimension(2) :: matelem,tmp1,tmp2
    real(kind=KR)  :: normval
    integer(kind=KI), dimension(4) :: latsize
    integer(kind=KI) :: exists


    !store the lattice size in a local array
    latsize(1)=nx
    latsize(2)=ny
    latsize(3)=nz
    latsize(4)=nt

    idag=0 !applying the matrix itself not the hermitian conjugate

    nc=3 !number of colors
    ns=4 !number of spins


     if (myid == 0) then
 inquire(file=trim(rwdir(myid+1))//"hopping-matrix.LOG", exist=exists)
  if (.not. exists) then
         print *, "File does not exist. Creating it."

open(unit=500,file=trim(rwdir(myid+1))//"hopping-matrix.LOG",status="new",     &
        action="write",form="formatted")
   close(unit=500,status="keep")
end if
end if



    if(myid==0) then ! open the file
       open(unit=500,file=trim(rwdir(myid+1))//"hopping-matrix.LOG",action="write",&
                      form="formatted",status="old",position="append")
    endif

    do it=1,nt
     do iz=1,nz
      do iy=1,ny
       do ix=1,nx
        !source location for the row
        srcl(1)=ix
        srcl(2)=iy
        srcl(3)=iz
        srcl(4)=it
      
        if(modulo(ix+iy+iz+it,2)==0) then
          ieo=1   !even row
        else
          ieo=2   !odd row
        endif

        do id= 1,ns   !spin
         do ic=1,nc   !color

            !single row index
            is=ic+(id-1)*nc+(ix-1)*nc*ns+(iy-1)*nc*ns*nx+(iz-1)*nc*ns*nx*ny+(it-1)*nc*ns*nx*ny*nz
            call pointsource(bel,bol,srcl,ic,id,myid,iblv,vecblinv)

            !compute all nonzero column elemnts taking into account that only nearest neighbours contribute
            do i=1,4
               srcr(i)=srcl(i)
            enddo
            
            !neighbours in the +mu direction taking inot account the boundary conditions
            do i=1,4

               !take into account the boundary conditions      
               if( (bc(i) .ne. 0) .or. ( (bc(i)==0) .and. (srcl(i) .lt. latsize(i)) ) ) then !other than these cases the hopping matrix should be zero
                  srcr(i) = srcl(i) + 1
                  if(srcr(i) == (latsize(i)+1) ) then
                     srcr(i) = 1
                  endif

                  do jd=1,ns
                     do jc=1,nc

                        !single column index
                        js=jc+(jd-1)*nc+(srcr(1)-1)*nc*ns+(srcr(2)-1)*nc*ns*nx+(srcr(3)-1)*nc*ns*nx*ny+(srcr(4)-1)*nc*ns*nx*ny*nz
                        call pointsource(ber,bor,srcr,jc,jd,myid,iblv,vecblinv)
             
                        !apply the operator
                        betmp= 0.0_KR
                        botmp= 0.0_KR

                        if (ieo==1) then ! row is even and column is odd
                           gblclr=1
                           call Hsingle(betmp,u,bor,idag,coact,bc,gblclr,vecbl,vecblinv,myid,nn, &
                                    ldiv,nms,lvbc,ib,lbd,iblv,MRT)
                       else
                           gblclr=2  !row is odd and column is even
                           call Hsingle(botmp,u,ber,idag,coact,bc,gblclr,vecbl,vecblinv,myid,nn, &
                                       ldiv,nms,lvbc,ib,lbd,iblv,MRT)
                       endif

               
                       call vecdot(bel,betmp,tmp1,MRT2)
                       call vecdot(bol,botmp,tmp2,MRT2)
                       matelem(:)=tmp1(:)+tmp2(:)
                       normval = sqrt(matelem(1)*matelem(1)+matelem(2)*matelem(2))
                       if(normval > 1e-16) then
                         if(myid==0) then   ! write out this matrix element of the hopping matrix
                            write(unit=500,fmt="(2i20,2es20.10)")is,js,matelem(1),matelem(2)
                         endif
                       endif
                     enddo !jc
                  enddo !jd
                  srcr(i) = srcl(i) ! this brings us back to the original point
               endif ! if( (bc(i)/= 0) .......)
            enddo !do i=1,4
            
            !neighbours in the -mu directions
            do i=1,4
               !take into account the boundary conditions      
               if( (bc(i) /= 0) .or. ( (bc(i)==0) .and. (srcl(i) .gt. 1) ) )then !other than these cases the hopping matrix should be zero
                  srcr(i) = srcr(i)-1
                  if(srcr(i) == 0) then
                    srcr(i) = latsize(i)
                  endif
                                   
               
                  do jd=1,ns
                     do jc=1,nc

                        !single column index
                        js=jc+(jd-1)*nc+(srcr(1)-1)*nc*ns+(srcr(2)-1)*nc*ns*nx+(srcr(3)-1)*nc*ns*nx*ny+(srcr(4)-1)*nc*ns*nx*ny*nz
                        call pointsource(ber,bor,srcr,jc,jd,myid,iblv,vecblinv)
             
                        !apply the operator
                        betmp= 0.0_KR
                        botmp= 0.0_KR

                        if (ieo==1) then ! row is even and column is odd
                           gblclr=1
                           call Hsingle(betmp,u,bor,idag,coact,bc,gblclr,vecbl,vecblinv,myid,nn, &
                                    ldiv,nms,lvbc,ib,lbd,iblv,MRT)
                       else
                          gblclr=2  !row is odd and column is even
                          call Hsingle(botmp,u,ber,idag,coact,bc,gblclr,vecbl,vecblinv,myid,nn, &
                                    ldiv,nms,lvbc,ib,lbd,iblv,MRT)
                       endif

               
                       call vecdot(bel,betmp,tmp1,MRT2)
                       call vecdot(bol,botmp,tmp2,MRT2)
                       matelem(:)=tmp1(:)+tmp2(:)
                       normval = sqrt(matelem(1)*matelem(1)+matelem(2)*matelem(2))

                       if(normval > 1e-16) then
                          if(myid==0) then ! write out this matrix element of the hopping matrix
                                    write(unit=500,fmt="(2i20,2es20.10)")is,js,matelem(1),matelem(2)
                          endif
                       endif
                    enddo !jc
                  enddo !jd
                  srcr(i) = srcl(i) ! this brings us back to the original point
               endif ! if( (bc(i) /= 0) .or. ( (bc(i)==0) .and. (srcl(i) .gt. 1 ) )
            enddo !i=1,4 i.e. looping over all sites

       enddo !ic
      enddo !id
    enddo !ix
   enddo !iy
  enddo !iz
 enddo !it

 if(myid==0) then ! close the file
    close(unit=500,status="keep")
 endif
     
 end subroutine printH
!---------------------------------------------------------------------------------------------------------------

!print a fermion vector 
subroutine printferm(be,bo,rwdir,filename,myid,iblv,vecblinv,MRT2)

    character(len=*), intent(in),    dimension(:)         :: rwdir
    character(len=*), intent(in)                          :: filename
    integer(kind=KI), intent(in)                          :: myid, MRT2
    integer(kind=KI), intent(in),    dimension(:,:)       :: iblv, vecblinv
    real(kind=KR), intent(in), dimension(:,:,:,:,:)  :: be,bo !input source vector

    real(kind=KR), dimension(6,ntotal,4,2,8)  :: betmp,botmp  !auxilary point source for projecting components of the input vector
    integer(kind=KI), dimension(4)            :: src 
    integer(kind=KI)                          :: it,iz,iy,ix,ic,id,is
    integer(kind=KI)  :: nc,ns
    integer(kind=KI)  :: i,j,k
    real(kind=KR), dimension(2) :: vecelem,tmp1,tmp2
    integer(kind=KI) :: exists


    
    nc=3 !number of colors
    ns=4 !number of spins

    if(myid==0) then   ! open file for writing
      open(unit=80,file=trim(rwdir(myid+1))//trim(filename)//".LOG",action="write",&
           form="formatted",status="old",position="append")
    endif

    do it=1,nt
     do iz=1,nz
      do iy=1,ny
       do ix=1,nx
        !source location
        src(1)=ix
        src(2)=iy
        src(3)=iz
        src(4)=it
        do id= 1,ns   !spin
         do ic=1,nc   !color

            !single row index
            is=ic+(id-1)*nc+(ix-1)*nc*ns+(iy-1)*nc*ns*nx+(iz-1)*nc*ns*nx*ny+(it-1)*nc*ns*nx*ny*nz
            call pointsource(betmp,botmp,src,ic,id,myid,iblv,vecblinv)
            call vecdot(betmp,be,tmp1,MRT2)
            call vecdot(botmp,bo,tmp2,MRT2)
            vecelem(:)=tmp1(:)+tmp2(:)

            if(myid==0) then   ! write out this element of the vector
                write(unit=80,fmt="(i20,2es20.10)")is,vecelem(1),vecelem(2)
            endif

        enddo !ic
       enddo !id
      enddo !ix
     enddo !iy
    enddo !iz
   enddo !it

if(myid==0) then   ! close the file
  close(unit=80,status="keep")
endif
     
 end subroutine printferm

 end module printops
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

