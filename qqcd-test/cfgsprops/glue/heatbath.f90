! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
! heatbath.f90, Randy Lewis, randy.lewis@uregina.ca
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
! Perform a pseudo-heatbath update for a single link.
!   
! parts of this code are based on: 
! *QCDMPI (version 5, 1996) by Shinji Hioki, Parallel Computing 22 (1997) 1335.
! and
! *QCDimMPI (version 0.9, 1998) by Astushi Nakamura and Shinji Hioki,
!           Nucl Phys B (Proc Suppl) 73 (1999) 895.
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 module heatbath

    use kinds
    use basics
    use lagfib
    implicit none
    private

! Define access to subroutines.
    public  :: phbup
    private :: phbsub1, phbsub2, phbsub3, su2up, sprd1

 contains

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 subroutine phbup(mu,ieo,ibl,c,u,nhit)
! Pseudo-heatbath update.
! INPUT:
!   mu is the spacetime direction of the links to be updated.
!   ieo is the part of the checkerboard containing the links to be updated.
!   ibl is the block of the sublattice containing the links of interest.
!   c() contains the staples.
!       expected size: c(18,nvhalf)
!   nhit is the number of updates desired for each lattice site.
! OUTPUT:
!   The links, u(), are updated.

    integer(kind=KI), intent(in)                          :: mu, ieo, ibl, nhit
    real(kind=KR),    intent(in),    dimension(:,:)       :: c
    real(kind=KR),    intent(inout), dimension(:,:,:,:,:) :: u

    real(kind=KR), dimension(4,nvhalf) :: z
! z() has twice as many entries in QCDMPI&QCDimMPI, but they seem to go unused.

!*First SU(2) subgroup
! determine SU(2) projection
    call phbsub1(mu,ieo,ibl,c,u,z)
! perform the pseudo-heatbath update for this SU(2) subgroup
    call su2up(z,nhit)
! record the updated link variables
    call sprd1(mu,ieo,ibl,z,u,1)
 
!*Second SU(2) subgroup
! determine SU(2) projection
    call phbsub2(mu,ieo,ibl,c,u,z)
! perform the pseudo-heatbath update for this SU(2) subgroup
    call su2up(z,nhit)
! record the updated link variables
    call sprd1(mu,ieo,ibl,z,u,2)
 
!*Third SU(2) subgroup
! determine SU(2) projection
    call phbsub3(mu,ieo,ibl,c,u,z)
! perform the pseudo-heatbath update for this SU(2) subgroup
    call su2up(z,nhit)
! record the updated link variables
    call sprd1(mu,ieo,ibl,z,u,3)
 
 end subroutine phbup

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 subroutine phbsub1(mu,ieo,ibl,c,u,z)
! Contribution to first SU(2) subgroup pseudo-heatbath update.
! INPUT:
!   mu is the spacetime direction of the links to be updated.
!   ieo is the part of the checkerboard containing the links to be updated.
!   ibl is the block of the sublattice containing the links of interest.
!   c(:,i) is a sum of staples (not an SU(3) matrix).
!          expected size: c(18,nvhalf)
! OUTPUT:
!   The 2x2 matrix represented by rows 1 and 2 and columns 1 and 2 of
!   u*(c^dagger) is projected onto the following SU(2) matrix,
!   z(1,:)*identity + i*z(2,:)*sigma(1) + i*z(3,:)*sigma(2) + i*z(4,:)*sigma(3)
!   where sigma(j) is a Pauli matrix and i=sqrt(-1).
!   expected size: z(4,nvhalf)

    integer(kind=KI), intent(in)                        :: mu, ieo, ibl
    real(kind=KR),    intent(in),  dimension(:,:)       :: c
    real(kind=KR),    intent(in),  dimension(:,:,:,:,:) :: u
    real(kind=KR),    intent(out), dimension(:,:)       :: z

    integer(kind=KI) :: i

    do i = 1,nvhalf
     z(1,i) =   u(1 ,i,mu,ieo,ibl) * c(1 ,i) + u(2 ,i,mu,ieo,ibl) * c(2 ,i) &
              + u(3 ,i,mu,ieo,ibl) * c(3 ,i) + u(4 ,i,mu,ieo,ibl) * c(4 ,i) &
              + u(5 ,i,mu,ieo,ibl) * c(5 ,i) + u(6 ,i,mu,ieo,ibl) * c(6 ,i) &
              + u(7 ,i,mu,ieo,ibl) * c(7 ,i) + u(8 ,i,mu,ieo,ibl) * c(8 ,i) &
              + u(9 ,i,mu,ieo,ibl) * c(9 ,i) + u(10,i,mu,ieo,ibl) * c(10,i) &
              + u(11,i,mu,ieo,ibl) * c(11,i) + u(12,i,mu,ieo,ibl) * c(12,i)
     z(2,i) = - u(1 ,i,mu,ieo,ibl) * c(8 ,i) + u(2 ,i,mu,ieo,ibl) * c(7 ,i) &
              - u(3 ,i,mu,ieo,ibl) * c(10,i) + u(4 ,i,mu,ieo,ibl) * c(9 ,i) &
              - u(5 ,i,mu,ieo,ibl) * c(12,i) + u(6 ,i,mu,ieo,ibl) * c(11,i) &
              - u(7 ,i,mu,ieo,ibl) * c(2 ,i) + u(8 ,i,mu,ieo,ibl) * c(1 ,i) &
              - u(9 ,i,mu,ieo,ibl) * c(4 ,i) + u(10,i,mu,ieo,ibl) * c(3 ,i) &
              - u(11,i,mu,ieo,ibl) * c(6 ,i) + u(12,i,mu,ieo,ibl) * c(5 ,i)
     z(3,i) =   u(1 ,i,mu,ieo,ibl) * c(7 ,i) + u(2 ,i,mu,ieo,ibl) * c(8 ,i) &
              + u(3 ,i,mu,ieo,ibl) * c(9 ,i) + u(4 ,i,mu,ieo,ibl) * c(10,i) &
              + u(5 ,i,mu,ieo,ibl) * c(11,i) + u(6 ,i,mu,ieo,ibl) * c(12,i) &
              - u(7 ,i,mu,ieo,ibl) * c(1 ,i) - u(8 ,i,mu,ieo,ibl) * c(2 ,i) &
              - u(9 ,i,mu,ieo,ibl) * c(3 ,i) - u(10,i,mu,ieo,ibl) * c(4 ,i) &
              - u(11,i,mu,ieo,ibl) * c(5 ,i) - u(12,i,mu,ieo,ibl) * c(6 ,i)
     z(4,i) = - u(1 ,i,mu,ieo,ibl) * c(2 ,i) + u(2 ,i,mu,ieo,ibl) * c(1 ,i) &
              - u(3 ,i,mu,ieo,ibl) * c(4 ,i) + u(4 ,i,mu,ieo,ibl) * c(3 ,i) &
              - u(5 ,i,mu,ieo,ibl) * c(6 ,i) + u(6 ,i,mu,ieo,ibl) * c(5 ,i) &
              + u(7 ,i,mu,ieo,ibl) * c(8 ,i) - u(8 ,i,mu,ieo,ibl) * c(7 ,i) &
              + u(9 ,i,mu,ieo,ibl) * c(10,i) - u(10,i,mu,ieo,ibl) * c(9 ,i) &
              + u(11,i,mu,ieo,ibl) * c(12,i) - u(12,i,mu,ieo,ibl) * c(11,i)
    enddo ! i
 
 end subroutine phbsub1

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 subroutine phbsub2(mu,ieo,ibl,c,u,z)
! Contribution to second SU(2) subgroup pseudo-heatbath update.
! INPUT:
!   mu is the spacetime direction of the links to be updated.
!   ieo is the part of the checkerboard containing the links to be updated.
!   ibl is the block of the sublattice containing the links of interest.
!   c(:,i) is a sum of staples (not an SU(3) matrix).
!          expected size: c(18,nvhalf)
! OUTPUT:
!   The 2x2 matrix represented by rows 2 and 3 and columns 2 and 3 of
!   u*(c^dagger) is projected onto the following SU(2) matrix,
!   z(1,:)*identity + i*z(2,:)*sigma(1) + i*z(3,:)*sigma(2) + i*z(4,:)*sigma(3)
!   where sigma(j) is a Pauli matrix and i=sqrt(-1).
!   expected size: z(4,nvhalf)

    integer(kind=KI), intent(in)                        :: mu, ieo, ibl
    real(kind=KR),    intent(in),  dimension(:,:)       :: c
    real(kind=KR),    intent(in),  dimension(:,:,:,:,:) :: u
    real(kind=KR),    intent(out), dimension(:,:)       :: z

    integer(kind=KI) :: i
 
    do i = 1,nvhalf
     z(1,i) =   u(7 ,i,mu,ieo,ibl) * c(7 ,i) + u(8 ,i,mu,ieo,ibl) * c(8 ,i) &
              + u(9 ,i,mu,ieo,ibl) * c(9 ,i) + u(10,i,mu,ieo,ibl) * c(10,i) &
              + u(11,i,mu,ieo,ibl) * c(11,i) + u(12,i,mu,ieo,ibl) * c(12,i) &
              + u(13,i,mu,ieo,ibl) * c(13,i) + u(14,i,mu,ieo,ibl) * c(14,i) &
              + u(15,i,mu,ieo,ibl) * c(15,i) + u(16,i,mu,ieo,ibl) * c(16,i) &
              + u(17,i,mu,ieo,ibl) * c(17,i) + u(18,i,mu,ieo,ibl) * c(18,i)
     z(2,i) = - u(7 ,i,mu,ieo,ibl) * c(14,i) + u(8 ,i,mu,ieo,ibl) * c(13,i) &
              - u(9 ,i,mu,ieo,ibl) * c(16,i) + u(10,i,mu,ieo,ibl) * c(15,i) &
              - u(11,i,mu,ieo,ibl) * c(18,i) + u(12,i,mu,ieo,ibl) * c(17,i) &
              - u(13,i,mu,ieo,ibl) * c(8 ,i) + u(14,i,mu,ieo,ibl) * c(7 ,i) &
              - u(15,i,mu,ieo,ibl) * c(10,i) + u(16,i,mu,ieo,ibl) * c(9 ,i) &
              - u(17,i,mu,ieo,ibl) * c(12,i) + u(18,i,mu,ieo,ibl) * c(11,i)
     z(3,i) = - u(7 ,i,mu,ieo,ibl) * c(13,i) - u(8 ,i,mu,ieo,ibl) * c(14,i) &
              - u(9 ,i,mu,ieo,ibl) * c(15,i) - u(10,i,mu,ieo,ibl) * c(16,i) &
              - u(11,i,mu,ieo,ibl) * c(17,i) - u(12,i,mu,ieo,ibl) * c(18,i) &
              + u(13,i,mu,ieo,ibl) * c(7 ,i) + u(14,i,mu,ieo,ibl) * c(8 ,i) &
              + u(15,i,mu,ieo,ibl) * c(9 ,i) + u(16,i,mu,ieo,ibl) * c(10,i) &
              + u(17,i,mu,ieo,ibl) * c(11,i) + u(18,i,mu,ieo,ibl) * c(12,i) 
     z(4,i) = - u(7 ,i,mu,ieo,ibl) * c(8 ,i) + u(8 ,i,mu,ieo,ibl) * c(7 ,i) &
              - u(9 ,i,mu,ieo,ibl) * c(10,i) + u(10,i,mu,ieo,ibl) * c(9 ,i) &
              - u(11,i,mu,ieo,ibl) * c(12,i) + u(12,i,mu,ieo,ibl) * c(11,i) &
              + u(13,i,mu,ieo,ibl) * c(14,i) - u(14,i,mu,ieo,ibl) * c(13,i) &
              + u(15,i,mu,ieo,ibl) * c(16,i) - u(16,i,mu,ieo,ibl) * c(15,i) &
              + u(17,i,mu,ieo,ibl) * c(18,i) - u(18,i,mu,ieo,ibl) * c(17,i)
    enddo ! i
 
 end subroutine phbsub2

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 subroutine phbsub3(mu,ieo,ibl,c,u,z)
! Contribution to third SU(2) subgroup pseudo-heatbath update.
! INPUT:
!   mu is the spacetime direction of the links to be updated.
!   ieo is the part of the checkerboard containing the links to be updated.
!   ibl is the block of the sublattice containing the links of interest.
!   c(:,i) is a sum of staples (not an SU(3) matrix).
!          expected size: c(18,nvhalf)
! OUTPUT:
!   The 2x2 matrix represented by rows 1 and 3 and columns 1 and 3 of
!   u*(c^dagger) is projected onto the following SU(2) matrix,
!   z(1,:)*identity + i*z(2,:)*sigma(1) + i*z(3,:)*sigma(2) + i*z(4,:)*sigma(3)
!   where sigma(j) is a Pauli matrix and i=sqrt(-1).
!   expected size: z(4,nvhalf)

    integer(kind=KI), intent(in)                        :: mu, ieo, ibl
    real(kind=KR),    intent(in),  dimension(:,:)       :: c
    real(kind=KR),    intent(in),  dimension(:,:,:,:,:) :: u
    real(kind=KR),    intent(out), dimension(:,:)       :: z

    integer(kind=KI) :: i
 
    do i = 1,nvhalf
     z(1,i) =   u(1 ,i,mu,ieo,ibl) * c(1 ,i) + u(2 ,i,mu,ieo,ibl) * c(2 ,i) &
              + u(3 ,i,mu,ieo,ibl) * c(3 ,i) + u(4 ,i,mu,ieo,ibl) * c(4 ,i) &
              + u(5 ,i,mu,ieo,ibl) * c(5 ,i) + u(6 ,i,mu,ieo,ibl) * c(6 ,i) &
              + u(13,i,mu,ieo,ibl) * c(13,i) + u(14,i,mu,ieo,ibl) * c(14,i) &
              + u(15,i,mu,ieo,ibl) * c(15,i) + u(16,i,mu,ieo,ibl) * c(16,i) &
              + u(17,i,mu,ieo,ibl) * c(17,i) + u(18,i,mu,ieo,ibl) * c(18,i)
     z(2,i) = - u(1 ,i,mu,ieo,ibl) * c(14,i) + u(2 ,i,mu,ieo,ibl) * c(13,i) &
              - u(3 ,i,mu,ieo,ibl) * c(16,i) + u(4 ,i,mu,ieo,ibl) * c(15,i) &
              - u(5 ,i,mu,ieo,ibl) * c(18,i) + u(6 ,i,mu,ieo,ibl) * c(17,i) &
              - u(13,i,mu,ieo,ibl) * c(2 ,i) + u(14,i,mu,ieo,ibl) * c(1 ,i) &
              - u(15,i,mu,ieo,ibl) * c(4 ,i) + u(16,i,mu,ieo,ibl) * c(3 ,i) &
              - u(17,i,mu,ieo,ibl) * c(6 ,i) + u(18,i,mu,ieo,ibl) * c(5 ,i)
     z(3,i) =   u(1 ,i,mu,ieo,ibl) * c(13,i) + u(2 ,i,mu,ieo,ibl) * c(14,i) &
              + u(3 ,i,mu,ieo,ibl) * c(15,i) + u(4 ,i,mu,ieo,ibl) * c(16,i) &
              + u(5 ,i,mu,ieo,ibl) * c(17,i) + u(6 ,i,mu,ieo,ibl) * c(18,i) &
              - u(13,i,mu,ieo,ibl) * c(1 ,i) - u(14,i,mu,ieo,ibl) * c(2 ,i) &
              - u(15,i,mu,ieo,ibl) * c(3 ,i) - u(16,i,mu,ieo,ibl) * c(4 ,i) &
              - u(17,i,mu,ieo,ibl) * c(5 ,i) - u(18,i,mu,ieo,ibl) * c(6 ,i)
     z(4,i) = - u(1 ,i,mu,ieo,ibl) * c(2 ,i) + u(2 ,i,mu,ieo,ibl) * c(1 ,i) &
              - u(3 ,i,mu,ieo,ibl) * c(4 ,i) + u(4 ,i,mu,ieo,ibl) * c(3 ,i) &
              - u(5 ,i,mu,ieo,ibl) * c(6 ,i) + u(6 ,i,mu,ieo,ibl) * c(5 ,i) &
              + u(13,i,mu,ieo,ibl) * c(14,i) - u(14,i,mu,ieo,ibl) * c(13,i) &
              + u(15,i,mu,ieo,ibl) * c(16,i) - u(16,i,mu,ieo,ibl) * c(15,i) &
              + u(17,i,mu,ieo,ibl) * c(18,i) - u(18,i,mu,ieo,ibl) * c(17,i)
    enddo ! i
 
 end subroutine phbsub3

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 subroutine su2up(z,nhit)
! Pseudo-heatbath update for an SU(2) subgroup.
! Algorithm is taken from Kennedy and Pendleton, Phys Lett 156B, 393 (1985).
! (also used by QCDMPI and QCDimMPI)
! INPUT:
!   z() is an SU(2) projection of the link to be updated times all staples.
!       expected size: z(4,nvhalf)
!   nhit is the number of pseudo-heatbath updates to be performed at each
!        lattice site.
! OUTPUT:
!   z() is updated and its value is output.

    real(kind=KR),    intent(inout), dimension(:,:) :: z
    integer(kind=KI), intent(in)                    :: nhit

    real(kind=KR), dimension(nvhalf)   :: a0, dk, id
    real(kind=KR), dimension(4*nvhalf) :: rn
    real(kind=KR)                      :: z1, z2, z3, z4, tt, d2, a1, &
                                          a2, a3, rad, rad2, theta
    integer(kind=KI)                   :: i, j, inotup, ihit, idum
    real(kind=KR), parameter           :: twoPi=6.283185307179586

!*Initialize id().
! id(i) will be set to 1 if the i'th link gets updated in this subroutine, 
!       else id(i) is set to 0.
    id = 0

!*Compute the inverse of the determinant of z().
    do i = 1,nvhalf
     dk(i) = 1.0_KR/sqrt(z(1,i)*z(1,i)+z(2,i)*z(2,i)+z(3,i)*z(3,i) &
                        +z(4,i)*z(4,i))
    enddo ! i

!*The Kennedy-Pendleton algorithm.
    do ihit = 1,nhit
     missnone: do
      inotup = nvhalf
      do i = 1,nvhalf
       inotup = inotup - id(i)
      enddo ! i
      if (inotup==0) then
       exit missnone
      else
       idum = 4*inotup
       call rnds(idum,rn)
       j = -3
       do i = 1,nvhalf
        if (id(i)/=1) then
         j = j + 4
         tt = cos(twoPi*rn(j+2))
         d2 = -(log(rn(j))+log(rn(j+1))*tt*tt)*dk(i)*3.0_KR
! (The above factor of 3 changes from SU(2) subgroup to SU(3) normalization.)
         if (rn(j+3)*rn(j+3)<=(1.0_KR-0.5_KR*d2)) then
          id(i) = 1
          a0(i) = 1.0_KR - d2
         endif
        endif
       enddo ! i
      endif
     enddo missnone
    enddo ! ihit

!*For each a0(i), generate (a1,a2,a3) uniformly on a hypersphere of radius 
! sqrt(1.-a0(i)*a0(i)).
    idum = 2*nvhalf
    call rnds(idum,rn)
    do i = 1,nvhalf
     j = i + i - 1
     rad = 1.0_KR - a0(i)*a0(i)
     a3 = sqrt(rad)*(2.0_KR*rn(j)-1.0_KR)
     rad2 = sqrt(abs(rad-a3*a3))
     theta = twoPi*rn(j+1)
     a1 = rad2*cos(theta)
     a2 = rad2*sin(theta)
     z1 = z(1,i)
     z2 = z(2,i)
     z3 = z(3,i)
     z4 = z(4,i)
     z(1,i) = dk(i)*( a0(i)*z1 + a1*z2 + a2*z3 + a3*z4)
     z(2,i) = dk(i)*(-a0(i)*z2 + a1*z1 + a2*z4 - a3*z3)
     z(3,i) = dk(i)*(-a0(i)*z3 - a1*z4 + a2*z1 + a3*z2)
     z(4,i) = dk(i)*(-a0(i)*z4 + a1*z3 - a2*z2 + a3*z1)
    enddo ! i

 end subroutine su2up

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 subroutine sprd1(mu,ieo,ibl,z,u,isu2)
! Replace the old links, u(), with the updated links.
! u[su(3)] = z[su(2)]*u[su(3)]
! INPUT:
!   mu is the spacetime direction of the links to be updated.
!   ieo is the part of the checkerboard containing the links to be updated.
!   ibl is the block of the sublattice containing the links of interest.
!   z() is the updated SU(2) projection of the links.
!       expected size: z(4,nvhalf)
!   isu2 = 1,2,3 identifies which of the 3 SU(2) subgroups contains z().
! OUTPUT:
!   The links, u(), are updated.

    integer(kind=KI), intent(in)                          :: mu, ieo, ibl, isu2
    real(kind=KR),    intent(in),    dimension(:,:)       :: z
    real(kind=KR),    intent(inout), dimension(:,:,:,:,:) :: u

    integer(kind=KI) :: i
    real(kind=KR)    :: y1, y2, y3, y4, y5, y6, y7, y8, y9, y10, &
                        y11, y12, y13, y14, y15, y16, y17, y18
 
    if (isu2==1) then
     do i = 1,nvhalf
      y1 =  z(1 ,i) * u(1 ,i,mu,ieo,ibl) - z(4 ,i) * u(2 ,i,mu,ieo,ibl) &
         +  z(3 ,i) * u(7 ,i,mu,ieo,ibl) - z(2 ,i) * u(8 ,i,mu,ieo,ibl) 
      y2 =  z(1 ,i) * u(2 ,i,mu,ieo,ibl) + z(4 ,i) * u(1 ,i,mu,ieo,ibl) &
         +  z(3 ,i) * u(8 ,i,mu,ieo,ibl) + z(2 ,i) * u(7 ,i,mu,ieo,ibl) 
      y3 =  z(1 ,i) * u(3 ,i,mu,ieo,ibl) - z(4 ,i) * u(4 ,i,mu,ieo,ibl) &
         +  z(3 ,i) * u(9 ,i,mu,ieo,ibl) - z(2 ,i) * u(10,i,mu,ieo,ibl) 
      y4 =  z(1 ,i) * u(4 ,i,mu,ieo,ibl) + z(4 ,i) * u(3 ,i,mu,ieo,ibl) &
         +  z(3 ,i) * u(10,i,mu,ieo,ibl) + z(2 ,i) * u(9 ,i,mu,ieo,ibl) 
      y5 =  z(1 ,i) * u(5 ,i,mu,ieo,ibl) - z(4 ,i) * u(6 ,i,mu,ieo,ibl) &
         +  z(3 ,i) * u(11,i,mu,ieo,ibl) - z(2 ,i) * u(12,i,mu,ieo,ibl) 
      y6 =  z(1 ,i) * u(6 ,i,mu,ieo,ibl) + z(4 ,i) * u(5 ,i,mu,ieo,ibl) &
         +  z(3 ,i) * u(12,i,mu,ieo,ibl) + z(2 ,i) * u(11,i,mu,ieo,ibl) 
      y7 = -z(3 ,i) * u(1 ,i,mu,ieo,ibl) - z(2 ,i) * u(2 ,i,mu,ieo,ibl) &
         +  z(1 ,i) * u(7 ,i,mu,ieo,ibl) + z(4 ,i) * u(8 ,i,mu,ieo,ibl) 
      y8 = -z(3 ,i) * u(2 ,i,mu,ieo,ibl) + z(2 ,i) * u(1 ,i,mu,ieo,ibl) &
         +  z(1 ,i) * u(8 ,i,mu,ieo,ibl) - z(4 ,i) * u(7 ,i,mu,ieo,ibl) 
      y9 = -z(3 ,i) * u(3 ,i,mu,ieo,ibl) - z(2 ,i) * u(4 ,i,mu,ieo,ibl) &
         +  z(1 ,i) * u(9 ,i,mu,ieo,ibl) + z(4 ,i) * u(10,i,mu,ieo,ibl) 
      y10= -z(3 ,i) * u(4 ,i,mu,ieo,ibl) + z(2 ,i) * u(3 ,i,mu,ieo,ibl) &
         +  z(1 ,i) * u(10,i,mu,ieo,ibl) - z(4 ,i) * u(9 ,i,mu,ieo,ibl) 
      y11= -z(3 ,i) * u(5 ,i,mu,ieo,ibl) - z(2 ,i) * u(6 ,i,mu,ieo,ibl) &
         +  z(1 ,i) * u(11,i,mu,ieo,ibl) + z(4 ,i) * u(12,i,mu,ieo,ibl) 
      y12= -z(3 ,i) * u(6 ,i,mu,ieo,ibl) + z(2 ,i) * u(5 ,i,mu,ieo,ibl) &
         +  z(1 ,i) * u(12,i,mu,ieo,ibl) - z(4 ,i) * u(11,i,mu,ieo,ibl) 
      u(1 ,i,mu,ieo,ibl) = y1
      u(2 ,i,mu,ieo,ibl) = y2
      u(3 ,i,mu,ieo,ibl) = y3
      u(4 ,i,mu,ieo,ibl) = y4
      u(5 ,i,mu,ieo,ibl) = y5
      u(6 ,i,mu,ieo,ibl) = y6
      u(7 ,i,mu,ieo,ibl) = y7
      u(8 ,i,mu,ieo,ibl) = y8
      u(9 ,i,mu,ieo,ibl) = y9
      u(10,i,mu,ieo,ibl) = y10
      u(11,i,mu,ieo,ibl) = y11
      u(12,i,mu,ieo,ibl) = y12
     enddo ! i
    elseif (isu2==2) then
     do i = 1,nvhalf
      y7 = z(1 ,i) * u(7 ,i,mu,ieo,ibl) - z(4 ,i) * u(8 ,i,mu,ieo,ibl) &
         - z(3 ,i) * u(13,i,mu,ieo,ibl) - z(2 ,i) * u(14,i,mu,ieo,ibl) 
      y8 = z(1 ,i) * u(8 ,i,mu,ieo,ibl) + z(4 ,i) * u(7 ,i,mu,ieo,ibl) &
         - z(3 ,i) * u(14,i,mu,ieo,ibl) + z(2 ,i) * u(13,i,mu,ieo,ibl) 
      y9 = z(1 ,i) * u(9 ,i,mu,ieo,ibl) - z(4 ,i) * u(10,i,mu,ieo,ibl) &
         - z(3 ,i) * u(15,i,mu,ieo,ibl) - z(2 ,i) * u(16,i,mu,ieo,ibl) 
      y10= z(1 ,i) * u(10,i,mu,ieo,ibl) + z(4 ,i) * u(9 ,i,mu,ieo,ibl) &
         - z(3 ,i) * u(16,i,mu,ieo,ibl) + z(2 ,i) * u(15,i,mu,ieo,ibl) 
      y11= z(1 ,i) * u(11,i,mu,ieo,ibl) - z(4 ,i) * u(12,i,mu,ieo,ibl) &
         - z(3 ,i) * u(17,i,mu,ieo,ibl) - z(2 ,i) * u(18,i,mu,ieo,ibl) 
      y12= z(1 ,i) * u(12,i,mu,ieo,ibl) + z(4 ,i) * u(11,i,mu,ieo,ibl) &
         - z(3 ,i) * u(18,i,mu,ieo,ibl) + z(2 ,i) * u(17,i,mu,ieo,ibl) 
      y13= z(3 ,i) * u(7 ,i,mu,ieo,ibl) - z(2 ,i) * u(8 ,i,mu,ieo,ibl) &
         + z(1 ,i) * u(13,i,mu,ieo,ibl) + z(4 ,i) * u(14,i,mu,ieo,ibl) 
      y14= z(3 ,i) * u(8 ,i,mu,ieo,ibl) + z(2 ,i) * u(7 ,i,mu,ieo,ibl) &
         + z(1 ,i) * u(14,i,mu,ieo,ibl) - z(4 ,i) * u(13,i,mu,ieo,ibl) 
      y15= z(3 ,i) * u(9 ,i,mu,ieo,ibl) - z(2 ,i) * u(10,i,mu,ieo,ibl) &
         + z(1 ,i) * u(15,i,mu,ieo,ibl) + z(4 ,i) * u(16,i,mu,ieo,ibl) 
      y16= z(3 ,i) * u(10,i,mu,ieo,ibl) + z(2 ,i) * u(9 ,i,mu,ieo,ibl) &
         + z(1 ,i) * u(16,i,mu,ieo,ibl) - z(4 ,i) * u(15,i,mu,ieo,ibl) 
      y17= z(3 ,i) * u(11,i,mu,ieo,ibl) - z(2 ,i) * u(12,i,mu,ieo,ibl) &
         + z(1 ,i) * u(17,i,mu,ieo,ibl) + z(4 ,i) * u(18,i,mu,ieo,ibl) 
      y18= z(3 ,i) * u(12,i,mu,ieo,ibl) + z(2 ,i) * u(11,i,mu,ieo,ibl) &
         + z(1 ,i) * u(18,i,mu,ieo,ibl) - z(4 ,i) * u(17,i,mu,ieo,ibl) 
      u(7 ,i,mu,ieo,ibl) = y7
      u(8 ,i,mu,ieo,ibl) = y8
      u(9 ,i,mu,ieo,ibl) = y9
      u(10,i,mu,ieo,ibl) = y10
      u(11,i,mu,ieo,ibl) = y11
      u(12,i,mu,ieo,ibl) = y12
      u(13,i,mu,ieo,ibl) = y13
      u(14,i,mu,ieo,ibl) = y14
      u(15,i,mu,ieo,ibl) = y15
      u(16,i,mu,ieo,ibl) = y16
      u(17,i,mu,ieo,ibl) = y17
      u(18,i,mu,ieo,ibl) = y18
     enddo ! i
    elseif (isu2==3) then
     do i = 1,nvhalf
      y1 =  z(1 ,i) * u(1 ,i,mu,ieo,ibl) - z(4 ,i) * u(2 ,i,mu,ieo,ibl) &
         +  z(3 ,i) * u(13,i,mu,ieo,ibl) - z(2 ,i) * u(14,i,mu,ieo,ibl) 
      y2 =  z(1 ,i) * u(2 ,i,mu,ieo,ibl) + z(4 ,i) * u(1 ,i,mu,ieo,ibl) &
         +  z(3 ,i) * u(14,i,mu,ieo,ibl) + z(2 ,i) * u(13,i,mu,ieo,ibl) 
      y3 =  z(1 ,i) * u(3 ,i,mu,ieo,ibl) - z(4 ,i) * u(4 ,i,mu,ieo,ibl) &
         +  z(3 ,i) * u(15,i,mu,ieo,ibl) - z(2 ,i) * u(16,i,mu,ieo,ibl) 
      y4 =  z(1 ,i) * u(4 ,i,mu,ieo,ibl) + z(4 ,i) * u(3 ,i,mu,ieo,ibl) &
         +  z(3 ,i) * u(16,i,mu,ieo,ibl) + z(2 ,i) * u(15,i,mu,ieo,ibl) 
      y5 =  z(1 ,i) * u(5 ,i,mu,ieo,ibl) - z(4 ,i) * u(6 ,i,mu,ieo,ibl) &
         +  z(3 ,i) * u(17,i,mu,ieo,ibl) - z(2 ,i) * u(18,i,mu,ieo,ibl) 
      y6 =  z(1 ,i) * u(6 ,i,mu,ieo,ibl) + z(4 ,i) * u(5 ,i,mu,ieo,ibl) &
         +  z(3 ,i) * u(18,i,mu,ieo,ibl) + z(2 ,i) * u(17,i,mu,ieo,ibl) 
      y13= -z(3 ,i) * u(1 ,i,mu,ieo,ibl) - z(2 ,i) * u(2 ,i,mu,ieo,ibl) &
         +  z(1 ,i) * u(13,i,mu,ieo,ibl) + z(4 ,i) * u(14,i,mu,ieo,ibl) 
      y14= -z(3 ,i) * u(2 ,i,mu,ieo,ibl) + z(2 ,i) * u(1 ,i,mu,ieo,ibl) &
         +  z(1 ,i) * u(14,i,mu,ieo,ibl) - z(4 ,i) * u(13,i,mu,ieo,ibl) 
      y15= -z(3 ,i) * u(3 ,i,mu,ieo,ibl) - z(2 ,i) * u(4 ,i,mu,ieo,ibl) &
         +  z(1 ,i) * u(15,i,mu,ieo,ibl) + z(4 ,i) * u(16,i,mu,ieo,ibl) 
      y16= -z(3 ,i) * u(4 ,i,mu,ieo,ibl) + z(2 ,i) * u(3 ,i,mu,ieo,ibl) &
         +  z(1 ,i) * u(16,i,mu,ieo,ibl) - z(4 ,i) * u(15,i,mu,ieo,ibl) 
      y17= -z(3 ,i) * u(5 ,i,mu,ieo,ibl) - z(2 ,i) * u(6 ,i,mu,ieo,ibl) &
         +  z(1 ,i) * u(17,i,mu,ieo,ibl) + z(4 ,i) * u(18,i,mu,ieo,ibl) 
      y18= -z(3 ,i) * u(6 ,i,mu,ieo,ibl) + z(2 ,i) * u(5 ,i,mu,ieo,ibl) &
         +  z(1 ,i) * u(18,i,mu,ieo,ibl) - z(4 ,i) * u(17,i,mu,ieo,ibl) 
      u(1 ,i,mu,ieo,ibl) = y1
      u(2 ,i,mu,ieo,ibl) = y2
      u(3 ,i,mu,ieo,ibl) = y3
      u(4 ,i,mu,ieo,ibl) = y4
      u(5 ,i,mu,ieo,ibl) = y5
      u(6 ,i,mu,ieo,ibl) = y6
      u(13,i,mu,ieo,ibl) = y13
      u(14,i,mu,ieo,ibl) = y14
      u(15,i,mu,ieo,ibl) = y15
      u(16,i,mu,ieo,ibl) = y16
      u(17,i,mu,ieo,ibl) = y17
      u(18,i,mu,ieo,ibl) = y18
     enddo ! i
    endif

 end subroutine sprd1

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 end module heatbath

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
