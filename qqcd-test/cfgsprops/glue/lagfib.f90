! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
! lagfib.f90, Randy Lewis, randy.lewis@uregina.ca
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
! Random number algorithm of 
! *Jun Makino, Parallel Computing, 20 (1994) 1357-1367.
!  "Lagged-Fibonacci random number generators on parallel computers"
! and
! *Jun Makino, Tetsuya Takaishi and Osamu Miyamura,
!  Computer Physics Communications 70 (1992) 495.
!  "Generation of shift register random numbers on distributed memory
!   multiprocessors"
! Makino considers four options: xor, +, - and *.  Here xor is used.
! (In f90, xor can also be called ieor.  The F compiler only accepts 
!  ieor so it will be used herein.)
!
! This algorithm is also used in 
! *QCDMPI by Shinji Hioki, Parallel Computing 22 (1997) 1335.
! and
! *QCDimMPI by Astushi Nakamura and Shinji Hioki,
!           Nucl Phys B (Proc Suppl) 73 (1999) 895.
!
! The algorithm assumes a 32-bit architecture and therefore all integers
! are explicitly set to 32-bit size.
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
! HOW TO USE THESE ROUTINES:
!  STEP 1: DO IT ONCE: call initlagfib(myid,iseed) from the main program.
!          myid is expected to identify each process: 0, 1, 2, ...
!          iseed must be a NONNEGATIVE integer.
!          Each process initiates random numbers beginning with the same
!          initial string of random numbers, but the starting points
!          within this string are delayed by (iseed+myid)*2^n.
!          This subroutine chooses n=40, as does QCDMPI; Makino chose 261.
!  STEP 2: REPEAT AS REQUIRED:
!          call rnds(m,rn) to get a set of m random numbers stored in
!                          at the beginning of array rn(:)
!          Note: the rn(:) are in [0,1), so exact 1's never appear.
!                Subroutine rnds has a simple test for avoiding exact 0's
!                if desired (otherwise comment it out).
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 module lagfib
!   parameters and subroutines relevant to this set of random number routines.
!   (lagged-Fibonacci method)

    use kinds
    implicit none
    private

! The algorithm assumes a 32-bit architecture and therefore all integers
! must be set to 32-bit size.  I.e. KI=4 for f90; KI=3 for F.
! (This is done in module kinds.)

! The lagged-Fibonacci method generates random numbers recursively
! from x_k = x_(k-p) xor x_(k-q).
! Here p and q are called iplagfib and iqlagfib; their values are
! taken from Makino, page 1364.
    integer(kind=KI), private, parameter :: iplagfib=521, iqlagfib=32

! The following three variables are for local use only, and they
! should not be changed by the user.
! NOTE: iwlagfib, jrlagfib, krlagfib are public so the main program can
!       write them to disk or define then from a file.
    integer(kind=KI), public, dimension(0:iplagfib-1), save :: iwlagfib
    integer(kind=KI), public, save :: jrlagfib, krlagfib

! Define access to subroutines.
    private :: initbirep, bigdelay
    public  :: initlagfib, rnds

 contains

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 subroutine initlagfib(myid,iseed)
! Initiate the random number generator.  Call once from main program.

    integer(kind=KI), intent(in) :: myid, iseed
    integer(kind=KI)             :: ndelay

! Set the random seed on each node as ndelay.
    ndelay = myid + iseed
    call initbirep()
! Important: ndelay must be a positive integer
    call bigdelay(ndelay)

 end subroutine initlagfib

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 subroutine initbirep()
! Initiate the random number generator (identically on each node).
! An F77 version of this subroutine is in Makino, Takaishi and Miyamura.

! ishft(a,b) is an intrinsic fortran function that shifts all bits in 
!            a to the left by b, truncating on the left and adding
!            zeros on the right.
!            The arguments and the function itself are 32-bit.
!            Example: In binary, 9 is '1001', so ishft(9,-2)='10'=2
! ieor(i,j) is the intrinsic fortran function "exclusive or".
!           The arguments and the function itself are 32-bit.
!           It uses exclusive or (1 for unequal, 0 for equal) for each bit.
!           Example: ieor(7,5)=2, since 7='111', 5='101' and 2='10'.

    integer(kind=KI), dimension(0:iplagfib-1) :: ib
    integer(kind=KI)                          :: ix, i, j, iwork

! Set ib(i) to a 'random' sequence of 0's and 1's.
    ix = 1
    do i = 0,iplagfib-1
     ix = ix*69069
     ib(i) = ishft(ix,-31)
    enddo ! i

! Initialize jrlagfib and krlagfib.
    jrlagfib = 0
    krlagfib = iplagfib - iqlagfib

! The first time jrlagfib counts up to iplagfib-1, iwork acquires 
! the initial entry from ib(), then jrlagfib begins at 0 again with
! a new ib() chosen according to the lagged-Fibonacci recurrence
! prescription (equation 1 of Makino).
    do j = 0,iplagfib-1
     iwork = 0
     do i = 0,31
      iwork = 2*iwork + ib(jrlagfib)
      ib(jrlagfib) = ieor(ib(jrlagfib),ib(krlagfib))
      jrlagfib = jrlagfib + 1
      if (jrlagfib==iplagfib) jrlagfib=0
      krlagfib = krlagfib + 1
      if (krlagfib==iplagfib) krlagfib=0
     enddo ! i
! Set iwlagfib(j) to a 'random' sequence of 0's and 1's, which is the
! initial entry in the original sequence of random numbers;
! each process will derive its own initialization from this one.
! Note: The ishft by -1 enforces iwlagfib(j) < 2^(31), as chosen by
!       Makino (see parenthetical remark between equations 1 and 2).
     iwlagfib(j) = ishft(iwork,-1)
    enddo ! j

 end subroutine initbirep

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 subroutine bigdelay(ndelay)
! Initiate the random number generator (independently on each node).
! An F77 version of this subroutine is in Makino, Takaishi and Miyamura.
! See also page 1365 of Makino.  (Makino might call this a combination 
! of SUBROUTINE BIREPX and SUBROUTINE DELAYX.)

    integer(kind=KI), intent(in)                :: ndelay
    integer(kind=KI)                            :: mu, m, i, j, nb, iwork
    integer(kind=KI), dimension(0:2*iplagfib-2) :: iwk
    integer(kind=KI), dimension(0:2*iplagfib-1) :: c
    integer(kind=KI), dimension(0:iplagfib+31)  :: ib

! Makino, and Makino,Takaishi,Miyamura choose mu=261,
! but QCDMPI and QCDimMPI choose mu=40.
    mu = 40

! Define iwk() to be iwlagfib(), but extended to include more terms
! according to the lagged-Fibonacci recurrence relation (Makino equation 1).
    do i = 0,iplagfib-1
     iwk(i) = iwlagfib(i)
    enddo ! i
    do i = iplagfib,2*iplagfib-2
     iwk(i) = ieor(iwk(i-iplagfib),iwk(i-iqlagfib))
    enddo ! i

! Initialize nb, m and the first few entries of ib().
    do i = 0,mu-1
     ib(i) = 0
    enddo ! i
    m = ndelay
    nb = mu - 1

! Increase nb if ndelay is greater than iplagfib-1.
    loopy: do
     if (m<=iplagfib-1) exit loopy
     nb = nb + 1
     ib(nb) = modulo(m,2)
     m = m/2
    enddo loopy

! Define c(), which will distinguish random numbers on different processes.
    do i = 0,iplagfib-1
     c(i) = 0
    enddo ! i
    c(m) = 1
    do j = nb,0,-1
     do i = iplagfib-1,0,-1
      c(2*i+ib(j)) = c(i)
      c(2*i+1-ib(j)) = 0
     enddo ! i
     do i = 2*iplagfib-1,iplagfib,-1
      c(i-iplagfib)=ieor(c(i-iplagfib),c(i))
      c(i-iqlagfib)=ieor(c(i-iqlagfib),c(i))
     enddo ! i
    enddo ! j

! Define iwlagfib() for this process.
    do j = 0,iplagfib-1
     iwork = 0
     do i = 0,iplagfib-1
      iwork = ieor(iwork,c(i)*iwk(j+i))
     enddo ! i
     iwlagfib(j) = iwork
    enddo ! j

 end subroutine bigdelay

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

     subroutine rnds(m,rn)
! This subroutine returns 'm' real random numbers in the inverval [0,1),
! and stores them as the first 'm' entries of rn(:).

     integer(kind=KI), intent(in)                  :: m
     real(kind=KR),    intent(inout), dimension(:) :: rn
     real(kind=KR), parameter                      :: fnorm=0.465661e-9
     integer(kind=KI)                              :: i 

     !print *, "hello1"
 
     do i = 1,m
      !print '(2i10)', jrlagfib,krlagfib
      !print '(2i10)', jrlagfib,krlagfib
! The following do loop removes random numbers which are exactly zero.
! (Reason: my pseudo-heathbath application does not appreciate exact zeros.)
! If you want to keep exact zeros, just comment out these three lines...
! "nozeros: do", "if (rn(i)/=0.0_KR) exit nozeros", "enddo nozeros".
     nozeros: do
        iwlagfib(jrlagfib) = ieor(iwlagfib(jrlagfib),iwlagfib(krlagfib))
! fnorm = 2^(-31), so rn(i) is between 0 and 1.
      rn(i) = iwlagfib(jrlagfib)*fnorm
! Note: 2^(-31)=0.46566129...e-9 > fnorm, so rn(i)=1.0 should not be possible.
      if (rn(i)/=0.0_KR) exit nozeros
     enddo nozeros
     if (rn(i)>=1.0_KR) then
      open(unit=8,file="LAGFIB.ERROR",action="write",status="replace" &
          ,form="formatted")
       write(unit=8,fmt=*) "subroutine rnds: random number = ",rn(i)
      close(unit=8,status="keep")
      stop
     endif
      !print '(2i10)', jrlagfib,krlagfib
      jrlagfib = jrlagfib + 1
      if (jrlagfib>=iplagfib) jrlagfib=jrlagfib-iplagfib
      krlagfib = krlagfib + 1
      if (krlagfib>=iplagfib) krlagfib=krlagfib-iplagfib
     enddo ! i

     !print *,"hello2"
     end subroutine rnds

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 end module lagfib

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
