! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
! bigjob.f90, Randy Lewis, randy.lewis@uregina.ca
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 program bigjob
!   A management tool for monitoring and halting parallel jobs.

    implicit none
    integer                :: i, numnodes, newnumnodes, myint, system, status
    integer, dimension(13) :: node, newnode
    character(64)          :: jobname, mystring

! Welcome message.
    write(*,*)
    write(*,*) "openSPACE management tool for parallel codes"
    write(*,*)

    maindo : do

! Define the entire theory cluster.
     numnodes = 0
     do i = 1,13
      numnodes = numnodes + 1
      node(i) = 51 + i
     enddo ! i

! Identify the job of interest.
     write(*,*) 'Enter the program name to be monitored (or "quit"):'
     read(*,*) jobname
     if (jobname=="quit") exit maindo

     progdo : do

! Determine which nodes have this job running.
! (status=256 means that the job is not running on that node.)
      mystring = 'ssh space53 "ps -ef | grep '//trim(jobname)//'| grep -v grep"'
      newnumnodes = 0
      newnode = 0
      do i = 1,numnodes
       write(unit=mystring(10:11),fmt="(i2.2)") node(i)
       write(*,"(a10,i2,a1)") "On space",node(i),":"
       status = system(mystring)
       if (status/=256) then
        newnumnodes = newnumnodes + 1
        newnode(newnumnodes) = node(i)
       endif
      enddo ! i
      numnodes = newnumnodes
      node = newnode

! Kill a process.
      killdo : do
       write(*,*) "Make a kill or refresh the list? [2=kill, 1=refresh, 0=no]"
       read(*,*) myint
       if (myint==0) exit progdo
       if (myint==1) exit killdo
       mystring = 'ssh space53 "kill       "'
       write(*,*) "Choose the computer on which to make the kill"
       read(*,*) myint
       write(unit=mystring(10:11),fmt="(i2.2)") myint
       write(*,*) "Enter a PID to kill [none=any negative number]"
       read(*,*) myint
       write(unit=mystring(19:24),fmt="(i6.6)") myint
       if (myint>0) status = system(mystring)
      enddo killdo

     enddo progdo

    enddo maindo

 end program bigjob

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
