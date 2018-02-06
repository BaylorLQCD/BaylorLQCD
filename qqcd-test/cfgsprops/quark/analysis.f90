!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
! analysis.f90, Dean Darnell, Dean_Darnell@baylor.edu
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
! "Data Reduction and Error Analysis for the Physical Sciences"
!  P. Bevington, McGraw-Hill Book Co., 1969
!
!  Algorithums for analysis were taked from the above books.
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 module analysis

!   use MPI
    use kinds
    use latdims
    use basics
!   use vevbleft
    implicit none
    private

! Define access to subroutines.
    public  :: stdev, scalarsignal, scalarsignal2

!   private :: 

 contains

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 subroutine stdev(a,numcount,subtokeep,momtokeep,optokeep, &
                     nmom,nop,nsub,sigma,realchoose,maxtime)

! This subroutine calculates the time sliced standard deviation for the vector
! a.

! For exeriments where the uncertainties are all equal, delta(q_i) = delta(q_j)
!     the uncertainty in the measurement can be estimated by 
!
!    sigma = sqrt{1/(N-1) * (sum(x_i^2)/N - xbar^2)}
!
!    where xbar = (1/N)*sum(x_i) = mean
!
!

   real(kind=KR),    intent(in), dimension(:,:,:)                           :: a
   integer(kind=KI), intent(in)                                             :: nsub, nmom, nop
   integer(kind=KI), intent(in)                                             :: subtokeep, momtokeep,&
                                                                               optokeep
   integer,          intent(in)                                             :: numcount, maxtime
   real(kind=KR),    intent(out),dimension(2,maxtime)                       :: sigma
   logical,          intent(in)                                             :: realchoose
   real(kind=KR),                dimension(2,maxtime)                       :: x2,xbar
   real(kind=KR)                                                            :: totalcount
   real(kind=KR)                                                            :: fac1, fac2 , fac3
   integer(kind=KI)                                                         :: iop,it,imom,isub,icount,&
                                                                               usub  


  sigma = 0.0_KR
  xbar  = 0.0_KR
  x2    = 0.0_KR

  totalcount = real(numcount,KR)
  fac1 = 1.0_KR/(totalcount) 
  fac2 = 1.0_KR/(totalcount-1.0_KR)
  
  do it = 1,maxtime
     do icount= 1,numcount
        x2(:,it)   = x2(:,it)   + a(:,it,icount)**2 
        xbar(:,it) = xbar(:,it) + a(:,it,icount)
     enddo ! icount

     sigma(:,it) = x2(:,it) - fac1*xbar(:,it)**2
     sigma(:,it) = fac2*sigma(:,it)
     sigma(:,it) = sqrt(sigma(:,it))

     if (realchoose) then
        print "(a27,4i3,1es17.10)", "it,isub,imom,iop, Re(sigma)", it,subtokeep,momtokeep,optokeep, sigma(1,it)
     else 
        print "(a27,4i3,1es17.10)", "it,isub,imom,iop, Im(sigma)", it,subtokeep,momtokeep,optokeep, sigma(2,it)
     endif ! real or imag
  enddo ! it

 end subroutine stdev

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 subroutine scalarsignal(a,noiseval,isigwrite,ntmqcd,nmom,nop,nsub,myid,rwdir)

 real(kind=KR),    intent(in),     dimension(:,:,:,:,:,:)   :: a
 character(len=*), intent(in),     dimension(:)           :: rwdir
 integer(kind=KI), intent(in)                             :: nmom,nop,nsub,ntmqcd
 integer(kind=KI), intent(in)                             :: myid
 integer(kind=KI), intent(in)                             :: isigwrite, noiseval
 
 real(kind=KR),                    dimension(3,nt,2,nmom,nop) :: sig
 real(kind=KR)                            :: fac,fac1,fac2
 character(len=128)                       :: opfile,tempfile
 character(len=2)                         :: trailer
 integer(kind=KI) :: iop,inoise,iri,it,imom
 
 integer(kind=KI) :: ltemp,utemp,ntt
 logical :: fileexists


  ltemp = 1
  utemp = nt
  ntt = utemp - ltemp + 1
   
 sig = 0.0_KR

! Need to normalize these values by the lattice size nxyzt=nx*ny*nz*nt
  fac1 = 1.0_KR/real(nxyz*ntt,KR)
  fac2 = 1.0_KR/real(nxyz,KR)

 if (myid==0) then
   if(ntmqcd > 0) then
      opfile = trim(rwdir(myid+1))//trim("Signal")//trim(".dat")//trim(".tmU")
   else
      opfile = trim(rwdir(myid+1))//trim("Signal")//trim(".dat")//trim(".tmD")
   endif ! ntmqcd
     
   tempfile = opfile
   trailer = ".x"

   if(nsub==6)  then

   do iop=1,nop
     write(unit=trailer(2:2),fmt="(i1.1)") iop
     opfile = trim(opfile)//trim(trailer)
     open(unit=8,file=opfile,action="write",form="formatted",status="new")
     do imom = 1,nmom
       do iri=1,2
         do it=ltemp,utemp
           do inoise=1,noiseval
      
             sig(1,it,iri,imom,iop) = sig(1,it,iri,imom,iop) + a(iri,it,1,imom,iop,inoise)
             sig(2,it,iri,imom,iop) = sig(2,it,iri,imom,iop) + a(iri,it,4,imom,iop,inoise)
             sig(3,it,iri,imom,iop) = sig(3,it,iri,imom,iop) + a(iri,it,6,imom,iop,inoise)

             if(mod(inoise,isigwrite)==0 .and. it==nt) then
             fac = 1.0_KR/real(inoise,KR)
                    write(unit=8, fmt="(3i3,3d24.10)") inoise,imom,iri, fac*fac1*sig(1,it,iri,imom,iop),&
                                                       fac*fac1*sig(2,it,iri,imom,iop), fac*fac1*sig(3,it,iri,imom,iop)
             endif ! mod

         enddo ! inoise
       enddo ! it
     enddo ! iri 
     enddo ! imom
   close(unit=8,status="keep")
   opfile = tempfile
   enddo ! iop

  else if (nsub==4) then

   do iop=1,nop
     write(unit=trailer(2:2),fmt="(i1.1)") iop
     opfile = trim(opfile)//trim(trailer)
     open(unit=8,file=opfile,action="write",form="formatted",status="new")
     do imom = 1,nmom
       do iri=1,2
         do it=ltemp,utemp
           do inoise=1,noiseval

             sig(1,it,iri,imom,iop) = sig(1,it,iri,imom,iop) + a(iri,it,1,imom,iop,inoise)
             sig(2,it,iri,imom,iop) = sig(2,it,iri,imom,iop) + a(iri,it,4,imom,iop,inoise)

             if(mod(inoise,isigwrite)==0 .and. it==nt) then
             fac = 1.0_KR/real(inoise,KR)
                    write(unit=8, fmt="(3i3,2d24.10)") inoise,imom,iri, fac*fac1*sig(1,it,iri,imom,iop),&
                                                       fac*fac1*sig(2,it,iri,imom,iop)
             endif ! mod

         enddo ! inoise
       enddo ! it
     enddo ! iri
     enddo ! imom
   close(unit=8,status="keep")
   opfile = tempfile
   enddo ! iop

   else if (nsub==0) then
!  do iop=1,nop
   iop=2
     write(unit=trailer(2:2),fmt="(i1.1)") iop
     opfile = trim(opfile)//trim(trailer)

     inquire(file=opfile, exist=fileexists)
     if (.not. fileexists) then
      open(unit=8,file=opfile,action="write",form="formatted",status="new")
      else
      open(unit=8,file=opfile,action="write",form="formatted",status="old",position="append")
     endif ! open
!       do inoise=1,noiseval

          write(unit=8,  fmt="(i4)") noiseval

           do it=ltemp,utemp

! This prints our zero momentum .....

             write(unit=8, fmt="(1i3,2d24.10)") it,fac2*a(1,it,1,1,2,noiseval),fac2*a(2,it,1,1,2,noiseval)

         enddo ! it
!      enddo ! inoise
   close(unit=8,status="keep")

   endif ! nsub
 endif ! myid 

!call printlog("Done with Signal.dat",myid,rwdir)

 end subroutine scalarsignal
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 subroutine scalarsignal2(a,noiseval,isigwrite,ntmqcd,nmom,nop,nsub,myid,rwdir)

!real(kind=KR),    intent(in),    dimension(numniose,2,nt,nsub,nmom,nop) :: a
 real(kind=KR),    intent(in),     dimension(:,:,:,:,:,:)   :: a
 character(len=*), intent(in),     dimension(:)           :: rwdir
 integer(kind=KI), intent(in)                             :: nmom,nop,nsub,ntmqcd
 integer(kind=KI), intent(in)                             :: myid
 integer(kind=KI), intent(in)                             :: isigwrite, noiseval

 real(kind=KR),                    dimension(3,nt,2,nmom,nop) :: sig
 real(kind=KR)                            :: fac,fac1,fac2
 character(len=128)                       :: opfile,tempfile
 character(len=2)                         :: trailer
 integer(kind=KI) :: iop,inoise,iri,it,imom

  integer(kind=KI) :: ltemp,utemp,ntt
 logical :: fileexists

  ltemp = 1
  utemp = nt
  ntt = utemp - ltemp + 1

 sig = 0.0_KR

! Need to normalize these values by the lattice size nxyzt=nx*ny*nz*nt
  fac1 = 1.0_KR/real(nx*ny*nz*ntt,KR)
  fac2 = 1.0_KR/real(nx*ny*nz,KR)

 if (myid==0) then
   if(ntmqcd > 0) then
      opfile = trim(rwdir(myid+1))//trim("Signal")//trim(".dat")//trim(".tmU")
   else
      opfile = trim(rwdir(myid+1))//trim("Signal")//trim(".dat")//trim(".tmD")
   endif ! ntmqcd

   tempfile = opfile
   trailer = ".x"

   if(nsub==6)  then

   do iop=1,nop
     write(unit=trailer(2:2),fmt="(i1.1)") iop
     opfile = trim(opfile)//trim(trailer)
     open(unit=8,file=opfile,action="write",form="formatted",status="new")
     do imom = 1,nmom
       do iri=1,2
         do it=ltemp,utemp
           do inoise=1,noiseval

             sig(1,it,iri,imom,iop) = sig(1,it,iri,imom,iop) + a(inoise,iri,it,1,imom,iop)
             sig(2,it,iri,imom,iop) = sig(2,it,iri,imom,iop) + a(inoise,iri,it,4,imom,iop)
             sig(3,it,iri,imom,iop) = sig(3,it,iri,imom,iop) + a(inoise,iri,it,6,imom,iop)

             if(mod(inoise,isigwrite)==0 .and. it==nt) then
             fac = 1.0_KR/real(inoise,KR)
                    write(unit=8, fmt="(3i3,3d24.10)") inoise,imom,iri, fac*fac1*sig(1,it,iri,imom,iop),&
                                                       fac*fac1*sig(2,it,iri,imom,iop), fac*fac1*sig(3,it,iri,imom,iop)
             endif ! mod

         enddo ! inoise
       enddo ! it
     enddo ! iri
     enddo ! imom
   close(unit=8,status="keep")
   opfile = tempfile
   enddo ! iop

  else if (nsub==4) then

   do iop=1,nop
     write(unit=trailer(2:2),fmt="(i1.1)") iop
     opfile = trim(opfile)//trim(trailer)
     open(unit=8,file=opfile,action="write",form="formatted",status="new")
     do imom = 1,nmom
       do iri=1,2
         do it=ltemp,utemp
           do inoise=1,noiseval

             sig(1,it,iri,imom,iop) = sig(1,it,iri,imom,iop) + a(inoise,iri,it,1,imom,iop)
             sig(2,it,iri,imom,iop) = sig(2,it,iri,imom,iop) + a(inoise,iri,it,4,imom,iop)

             if(mod(inoise,isigwrite)==0 .and. it==nt) then
             fac = 1.0_KR/real(inoise,KR)
                    write(unit=8, fmt="(3i3,2d24.10)") inoise,imom,iri, fac*fac1*sig(1,it,iri,imom,iop),&
                                                       fac*fac1*sig(2,it,iri,imom,iop)
             endif ! mod

         enddo ! inoise
       enddo ! it
     enddo ! iri
     enddo ! imom
   close(unit=8,status="keep")
   opfile = tempfile
   enddo ! iop

   else if (nsub==0) then
! This prints our zero momentum .....
!  do iop=1,nop
   iop=2
     write(unit=trailer(2:2),fmt="(i1.1)") iop
     opfile = trim(opfile)//trim(trailer)

     inquire(file=opfile, exist=fileexists)

     if (.not. fileexists) then
         open(unit=8,file=opfile,action="write",form="formatted",status="new")
     else
         open(unit=8,file=opfile,action="write",form="formatted",status="old",position="append")
         write(unit=8, fmt="(a32)") "Appending file" 
     endif ! open
     do inoise=1,noiseval
        write(unit=8,  fmt="(i4)") inoise
        do it=ltemp,utemp
           write(unit=8, fmt="(1i3,2d24.10)") it,fac2*a(1,it,1,1,2,inoise),fac2*a(2,it,1,1,2,inoise)
        enddo ! it
     enddo ! inoise
     write(unit=8,fmt=*) "fac2=",fac2
     close(unit=8,status="keep")

     endif ! nsub
 endif ! myid

! call printlog("Done with Signal.dat",myid,rwdir)

 end subroutine scalarsignal2
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 end module analysis
