! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
! basics.f90, Randy Lewis, randy.lewis@uregina.ca
!    
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
! nps is the total number of processes.
! nvhalf is the number of even (or odd) lattices sites per process per block.
! nbmax is the number of spacetime sites in the largest boundary
!       shared by two processes.  (The boundary is 3-D.)
! nwhole is the total number of even and odd lattice sites per process per block.
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 module basics
    use kinds
    use latdims
    implicit none
    private
    integer(kind=KI), public, parameter :: nxyzt=nx*ny*nz*nt
    integer(kind=KI), public, parameter :: nxyz =nx*ny*nz
    integer(kind=KI), public, parameter :: nps=npx*npy*npz*npt
    integer(kind=KI), public, parameter :: nvhalf=nxyzt/32/nps
!   integer(kind=KI), public, parameter :: nvwhole=2*nvhalf
    integer(kind=KI), public, parameter :: nbmax=max(2*nvhalf*npx/nx, &
                                                     2*nvhalf*npy/ny, &
                                                     2*nvhalf*npz/nz, &
                                                     2*nvhalf*npt/nt)
    integer(kind=KI), public, parameter :: ntotal=nvhalf+nbmax
    integer(kind=KI), public, parameter :: nvhalff=ntotal
!   integer(kind=KI), public, parameter :: nwhole=2*ntotal


    public:: gam5x,cmpxsqrt,cmpxmult,cmpxdiv,conj,printmat, &
             matvecmult, matmult, printvec, vecdotvec, & 
             vecdagvec, sort, oneover


    contains

subroutine npsvalue
print *,"npsvalue=",nps
end subroutine npsvalue

 ! Complex Division : 1/(c+di)
 ! Input: denom : denominator of return value
 ! Output: retVal : 1 / (c+di)
 subroutine oneover(denom, retVal)
 real(kind=KR2), intent(in),  dimension(2) :: denom
 real(kind=KR2), intent(out), dimension(2) :: retVal

 real(kind=KR2) :: tmp1,tmp2
 real(kind=KR2) :: c,d

 if (.false.) then ! the naive technique
   call conj(denom, retVal)
   retVal(1) = retVal(1) / (denom(1)*denom(1) + denom(2)*denom(2))
   retVal(2) = retVal(2) / (denom(1)*denom(1) + denom(2)*denom(2))
 else ! the Num. Rec. technique
   c = denom(1)
   d = denom(2)

   if (abs(c) .ge. abs(d)) then
     tmp1 = d / c ! (d/c)
     tmp2 = (tmp1 * d) + c

     retVal(1) = (1.0_KR2) / tmp2
     retVal(2) = -tmp1 / tmp2
   else ! abs(denom(1)) .lt. abs(denom(2))
     tmp1 = c / d ! (c/d)
     tmp2 = (tmp1 * c) + d

     retVal(1) = tmp1 / tmp2
     retVal(2) = (-1.0_KR2) / tmp2
   endif 
 endif
 end subroutine oneover



 subroutine sort(vec,m,indx)
 real(kind=KR2),   intent(in),  dimension(:) :: vec
 integer(kind=KI), intent(in)                  :: m
 integer(kind=KI), intent(out), dimension(:) :: indx

 real(kind=KR2), allocatable, dimension(:)   :: vectmp
 real(kind=KR2)  :: tmp
 integer(kind=KI) :: i,j
 integer(kind=KI) :: tmpint

 allocate(vectmp(m))

 vectmp(1:m) = vec(1:m)

 do i=1,m
   indx(i) = i
 enddo

 do j=1,m-1
   do i=1,m-j

     if (vectmp(i) > vectmp(i+1)) then
       tmpint = indx(i+1)
       indx(i+1) = indx(i)
       indx(i) = tmpint

       tmp = vectmp(i+1)
       vectmp(i+1) = vectmp(i)
       vectmp(i) = tmp
     endif
   enddo
 enddo

 deallocate(vectmp)
 end subroutine sort




    ! Matrix-Matrix multiplication
    ! Input: a   input matrix
    !        m1  leading dimension of a
    !        n1  second dimension of a
    !        b   input matrix
    !        m2  leading dimension of b
    !        n2  second dimension of b
    ! Output: c  Output matrix (solution of a*b); size m1,n2
    subroutine matmult(a,m1,n1,b,m2,n2,c)
    real(kind=KR2), intent(in), dimension(:,:,:)  :: a
    real(kind=KR2), intent(in), dimension(:,:,:)  :: b
    integer, intent(in)                :: m1,n1
    real(kind=KR2), intent(out), dimension(:,:,:) :: c
    integer, intent(in)                :: m2,n2

    integer   :: i,j,k
    real(kind=KR2), dimension(2)  :: tmp1

    if (n1 == m2) then
      c(:2,:m1,:n2) = 0.0_KR2
      do i=1,m1
        do j=1,n2
          do k=1,n1
            call cmpxmult(a(:2,i,k),b(:2,k,j),tmp1)
            c(:2,i,j) = c(:2,i,j) + tmp1
          enddo
        enddo
      enddo
    else
      write(*,*) 'ERROR in matrix size'
    endif
    end subroutine matmult

    ! Matrix-Vector multiplication
    ! Input: a   input matrix
    !        m1  leading dimension of a
    !        n1  second dimension of a
    !        b   input vector
    !        m2  leading dimension of b
    ! Output: c  Output vector (solution of a*b); size m1
    subroutine matvecmult(a,m1,n1,b,m2,c)
    real(kind=KR2), intent(in), dimension(:,:,:)  :: a
    real(kind=KR2), intent(in), dimension(:,:)  :: b
    integer, intent(in)                :: m1,n1
    real(kind=KR2), intent(out), dimension(:,:) :: c
    integer, intent(in)                :: m2

    integer   :: i,k
    real(kind=KR2), dimension(2)  :: tmp1

    if (n1 == m2) then
      c(:2,:m1) = 0.0_KR2
      do i=1,m1
        do k=1,n1
          call cmpxmult(a(:2,i,k),b(:2,k),tmp1)
          c(:2,i) = c(:2,i) + tmp1
        enddo
      enddo
    endif
    end subroutine matvecmult

    ! Vector Dot Vector (a^transpose * b)
    ! Input: a   input vector
    !        m  leading dimension of a
    !        b   input vector
    !        n  leading dimension of b
    ! Output: c  Output soluction
    subroutine vecdotvec(a,m,b,n,c)
    real(kind=KR2), intent(in), dimension(:,:) :: a
    real(kind=KR2), intent(in), dimension(:,:) :: b
    integer, intent(in)         :: m,n
    real(kind=KR2), intent(out), dimension(:)  :: c

    integer i

    if (m == n) then
      c(:2) = 0.0_KR2
      do i=1,n
        c(1) = c(1) + (a(1,i) * b(1,i)) - (a(2,i) * b(2,i))
        c(2) = c(2) + (a(1,i) * b(2,i)) + (a(2,i) * b(1,i))
      enddo
    else
      write(*,*) 'ERROR IN DIMENSIONS: vecdotvec'
    endif
    end subroutine vecdotvec

    ! Vector Dagger Vector (a^dagger * b)
    ! Input: a   input vector
    !        m  leading dimension of a
    !        b   input vector
    !        n  leading dimension of b
    ! Output: c  Output soluction
    subroutine vecdagvec(a,m,b,n,c)
    real(kind=KR2), intent(in), dimension(:,:) :: a
    real(kind=KR2), intent(in), dimension(:,:) :: b
    integer, intent(in)         :: m,n
    real(kind=KR2), intent(out), dimension(:)  :: c

    integer i

    if (m == n) then
      c(:2) = 0.0_KR2
      do i=1,n
        c(1) = c(1) + (a(1,i) * b(1,i)) + (a(2,i) * b(2,i))
        c(2) = c(2) + (a(1,i) * b(2,i)) - (a(2,i) * b(1,i))
      enddo
    else
      write(*,*) 'ERROR IN DIMENSIONS: vecdotvec'
    endif
    end subroutine vecdagvec



    ! Output Matrix to **Standard Output**
    ! Input: a  input matrix
    !        m  leading dimension of a
    !        n  second dimension of a
    subroutine printmat(a,m,n)
    integer, intent(in)    :: m,n
    real(kind=KR2), intent(in), dimension(:,:,:) :: a
    integer :: i,j
    do i=1,n
      do j=1,m
        write(*,*) a(:2,j,i)
      enddo
      write(*,*) '=========='
    enddo
    end subroutine printmat

    ! Output Vector to **Standard Output**
    ! Input: a  input vector
    !        m  leading dimension of a
    subroutine printvec(a,m)
    integer, intent(in)    :: m
    real(kind=KR2), intent(in), dimension(:,:) :: a
    integer :: j
    do j=1,m
      write(*,*) a(:2,j)
    enddo
    write(*,*) '=========='
    end subroutine printvec



    ! Complex Conjugate of Complex Number
    ! INPUT:
    !  z Complex number
    ! OUTPUT:
    !  zconj  complex conjugate of z
    subroutine conj(z,zconj)
    real(kind=KR2), intent(in), dimension(:) :: z
    real(kind=KR2), intent(out),dimension(:) ::zconj 
    zconj(1) =  z(1)
    zconj(2) = -z(2)
    end subroutine conj

    ! Square Root of Complex Number
    ! Algorithm:
    !   * let cpmxsqrt(a+bi)=p+qi
    !   * a+bi=(p+qi)**2
    !   * p=1/sqrt(2) * sqrt(a+sqrt(a**2+b**2))
    !   * q=sgn(b)/sqrt(2) * sqrt(sqrt(a**2+b**2)-a)
    ! INPUT:
    !  z  complex number (z(1)=a   z(2)=b)
    ! OUTPUT:
    !  zsqrt complex number (zsqrt(1)=p  zsqrt(2)=q)
    subroutine cmpxsqrt(z, zsqr)
    real(kind=KR2), intent(in), dimension(:) :: z
    real(kind=KR2), intent(out), dimension(:) :: zsqr

    real(kind=KR2) :: tmp1,tmp2
    real(kind=KR2) :: oneosqrt2
    real(kind=KR2) :: c,d,w


    if (.false.) then ! naive way
      tmp1 = sqrt(z(1)*z(1) + z(2)*z(2))
      zsqr(1) = sqrt((tmp1+z(1))/(2.0_KR2))
      zsqr(2) = sqrt((tmp1-z(1))/(2.0_KR2))

      if (z(2) < 0.0_KR2) then
        zsqr(2) = -1.0_KR2 * zsqr(2)
      endif
    else ! Num. Rec. way
      c = z(1)
      d = z(2)

      ! Compute "w"
      if ((c==0.0_KR2) .and. (d==0.0_KR2)) then
        w = 0.0_KR2
      else if (abs(c) .ge. abs(d)) then
        tmp1 = d/c
        tmp2 = tmp1 * tmp1
        tmp2 = 1.0_KR2 + sqrt(1.0_KR2 + tmp2)
        tmp2 = sqrt(tmp2 / 2.0_KR2)

        tmp1 = sqrt(abs(c))
        w = tmp1 * tmp2
      else ! abs(c) .lt. abs(d)
        tmp1 = c/d
        tmp2 = tmp1 * tmp1
        tmp2 = abs(tmp1) + sqrt(1.0_KR2 + tmp2)
        tmp2 = sqrt(tmp2 / 2.0_KR2)

        tmp1 = sqrt(abs(d))
        w = tmp1 * tmp2
      endif

      if (w == 0.0_KR2) then
        zsqr(1) = 0.0_KR2
        zsqr(2) = 0.0_KR2
      else if (c .ge. 0.0_KR2) then
        zsqr(1) = w
        zsqr(2) = d / (2.0_KR2 * w)
      else if ( (c .lt. 0.0_KR2) .and. (d .ge. 0.0_KR2)) then
        zsqr(1) = abs(d) / (2.0_KR2 * w)
        zsqr(2) = w
      else
        zsqr(1) = abs(d) / (2.0_KR2 * w)
        zsqr(2) = -w
      endif
    
    endif
    end subroutine cmpxsqrt


    ! Complex number multiplier a*b
    ! INPUT:
    !  a Complex Number
    !  b Complex Number
    ! OUTPUT:
    !  atimeb Complex Number
    subroutine cmpxmult(a,b,atimesb)
    real(kind=KR2), intent(in), dimension(:) :: a,b
    real(kind=KR2), intent(out), dimension(:) :: atimesb
 
    atimesb(1) = a(1)*b(1)-a(2)*b(2)
    atimesb(2) = a(2)*b(1)+a(1)*b(2)
    end subroutine cmpxmult


    ! Complex number division a/b
    ! INPUT:
    !  a Complex Number
    !  b Complex Number
    ! OUTPUT:
    !  aoverb Complex Number
    subroutine cmpxdiv(num,den,retVal)
    real(kind=KR2), intent(in), dimension(:)  :: num,den
    real(kind=KR2), intent(out), dimension(:) :: retVal

    real(kind=KR2) :: denom
    real(kind=KR2) :: a,b,c,d
    real(kind=KR2) :: tmp1, tmp2

    if (.false.) then ! naive way
      denom = den(1)**2+den(2)**2
      retVal(1) = num(1)*den(1)+num(2)*den(2)
      retVal(2) = num(2)*den(1)-num(1)*den(2)
      retVal(:2) = (1.0_KR2/denom) * retVal(:2)
    else ! Num Rec. way
      a = num(1)
      b = num(2)
      c = den(1)
      d = den(2)

      if (abs(c) .ge. abs(d)) then
        tmp1 = d / c ! (d/c)
        denom = c + (d * tmp1)

        retVal(1) = a + (b * tmp1)
        retVal(2) = b - (a * tmp1)
        retVal(:2) = retVal(:2) / denom
      else ! abs(c) .lt. abs(d)
        tmp1 = c / d ! (c/d)
        denom = d + (c * tmp1)

        retVal(1) = b + (a * tmp1)
        retVal(2) = (b * tmp1) - a
        retVal(:2) = retVal(:2) / denom
      endif
    endif
    end subroutine cmpxdiv





    subroutine gam5x(a)

      !calculates the product of a Dirac vector defined over globally odd
      !or globally even sites with gamma5. The convention for gamma5 is:
      !gamma5(1,3)=-i , gamma5(2,4)=-i , gamma5(3,1)=i , gamma5(4,2)=i
      !and the rest of elements are zero. 
      
      real(kind=KR), intent(inout)  ,  dimension(:,:,:,:,:) :: a
      real(kind=KR), dimension(6,nvhalf,4,2,8)                   :: b
      integer(kind=KI) :: isite,icri,idira,ieo,ibl


      b=0.0_KR
      do icri=1,6
       do isite=1,nvhalf
        do idira=1,4
         do ieo=1,2
          do ibl=1,8

           b(icri,isite,idira,ieo,ibl)=a(icri,isite,idira,ieo,ibl)

          enddo !ibl
         enddo !ieo
        enddo !idira
       enddo !isite
      enddo !icri
     
     do isite=1,nvhalf
      do icri=1,5,2
       do ieo=1,2
        do ibl=1,8

         a(icri  ,isite,1,:,:) =  b(icri+1,isite,3,:,:)
         a(icri+1,isite,1,:,:) = -b(icri  ,isite,3,:,:)
         a(icri  ,isite,2,:,:) =  b(icri+1,isite,4,:,:)
         a(icri+1,isite,2,:,:) = -b(icri  ,isite,4,:,:)
         a(icri  ,isite,3,:,:) = -b(icri+1,isite,1,:,:)
         a(icri+1,isite,3,:,:) =  b(icri  ,isite,1,:,:)
         a(icri  ,isite,4,:,:) = -b(icri+1,isite,2,:,:)
         a(icri+1,isite,4,:,:) =  b(icri  ,isite,2,:,:)

         enddo !ibl
        enddo !ieo
       enddo !icri
      enddo !isite
 
    end subroutine gam5x


!    subroutine gam5xgamm5(a)

      !calculates the product of a Dirac vector defined over globally odd
      !or globally even sites with gamma5 on left and right, i.e. gamma5*a*gamm5. 
      !The convention for gamma5 is:
      !gamma5(1,3)=-i , gamma5(2,4)=-i , gamma5(3,1)=i , gamma5(4,2)=i
      !and the rest of elements are zero. 
      
!      real(kind=KR), intent(inout)  ,  dimension(:,:,:,:,:) :: a
!      real(kind=KR), dimension(6,nvhalf,4,2,8)                   :: b
!      integer(kind=KI) :: isite,icri,idira,ieo,ibl


!      b=0.0_KR
!      do icri=1,6
!       do isite=1,nvhalf
!        do idira=1,4
!         do ieo=1,2
!          do ibl=1,8

!           b(icri,isite,idira,ieo,ibl)=a(icri,isite,idira,ieo,ibl)

!          enddo !ibl
!         enddo !ieo
!        enddo !idira
!       enddo !isite
!      enddo !icri
     
!     do isite=1,nvhalf
!      do icri=1,5,2
!       do ieo=1,2
!        do ibl=1,8

!         a(icri  ,isite,1,:,:) =  b(icri+1,isite,3,:,:)
!         a(icri+1,isite,1,:,:) = -b(icri  ,isite,3,:,:)
!         a(icri  ,isite,2,:,:) =  b(icri+1,isite,4,:,:)
!         a(icri+1,isite,2,:,:) = -b(icri  ,isite,4,:,:)
!         a(icri  ,isite,3,:,:) = -b(icri+1,isite,1,:,:)
!         a(icri+1,isite,3,:,:) =  b(icri  ,isite,1,:,:)
!         a(icri  ,isite,4,:,:) = -b(icri+1,isite,2,:,:)
!         a(icri+1,isite,4,:,:) =  b(icri  ,isite,2,:,:)

!         enddo !ibl
!        enddo !ieo
!       enddo !icri
!      enddo !isite
 
!    end subroutine gam5x

 end module basics

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
