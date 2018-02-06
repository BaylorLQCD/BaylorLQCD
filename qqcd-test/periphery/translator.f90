! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
! translator.f90, Randy Lewis, randy.lewis@uregina.ca
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 program translator
!   A management tool for conversion of files from one format to another.

    use kinds
    implicit none

    integer(kind=KI), parameter                     :: nxyztmax=100000
    real(kind=KR),    dimension(12,nxyztmax,4,2,16) :: u
    character(len=64)                               :: rwdir, cfgfile
    character(len=4)                                :: trailer
    real(kind=KR)                                   :: beta, kappa
    integer(kind=KI), dimension(4)                  :: npin, npout, ip, gaction
    real(kind=KR),    dimension(4)                  :: alat, tad
    integer(kind=KI), dimension(2,16)               :: vecblinv
    integer(kind=KI), dimension(16,4)               :: iblv
    integer(kind=KI) :: mainchoice, nx, ny, nz, nt, npsin, npsout, nxin, nyin, &
                        nzin, ntin, nxout, nyout, nzout, ntout, icfg, cfgmin, &
                        cfgmax, myid, ipoffset, icri, ix, iy, iz, it, ibl, &
                        isite, icfgsave, myidin, mu, icfgin, id, icsrcin, &
                        idsrcin, icsrc, idsrc, ikappa, ivbl, ibleo, jbleo, &
                        ieo, intx, inty, intz, intt

! A definition (copied from subroutine gblchker in module lattice).
    vecblinv(1,1:16) = (/ 1, 2, 2, 1, 2, 1, 1, 2, 2, 1, 1, 2, 1, 2, 2, 1 /)
    vecblinv(2,1:16) = (/ 1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7, 8, 8 /)

! Welcome message.
    write(unit=*,fmt=*)
    write(unit=*,fmt=*) "*** management tool for file conversion ***"
    write(unit=*,fmt=*)
    write(unit=*,fmt=*) "No files will be erased by this code."
    write(unit=*,fmt=*)
    write(unit=*,fmt=*) "Choose the operation to be performed:"
    write(unit=*,fmt=*)"   translate gauge fields to new number of processors.1"
    write(unit=*,fmt=*)"   translate propagators to new number of processors..2"
    write(unit=*,fmt=*)"   translate gauge fields from qqcd to rwww format....3"
    read(unit=*,fmt=*) mainchoice

! User inputs.
    write(unit=*,fmt=*) &
          "Enter the directory for all files (include final slash):"
    read(unit=*,fmt=*) rwdir
    write(unit=*,fmt=*) "Filenames must take the form AxxxxxyyyBCDE..."
    write(unit=*,fmt=*) "where xxxxx is the configuration number"
    write(unit=*,fmt=*) "and yyy is the process number."
    select case(mainchoice)
     case(1,3) ! gauge field translation
      write(unit=*,fmt=*) "Enter a template for the existing files:"
     case(2) ! propagator translation
      trailer = ".xxk"
      write(unit=*,fmt=*) "All filenames must end with the trailer .cdk ."
      write(unit=*,fmt=*) "Enter c,d,k as 3 numbers separated by spaces:"
      read(unit=*,fmt=*) icsrc, idsrc, ikappa
      trailer = ".xxk"
      write(unit=trailer(2:2),fmt="(i1.1)") icsrc
      write(unit=trailer(3:3),fmt="(i1.1)") idsrc
      write(unit=trailer(4:4),fmt="(i1.1)") ikappa
      write(unit=*,fmt=*) "Enter a template for existing files (no trailer):"
     case default
      open(unit=8,file="TRANSLATOR.ERROR",action="write",status="replace", &
           form="formatted")
       write(unit=8,fmt=*) "mainchoice =", mainchoice
      close(unit=8,status="keep")
      stop
    end select
    read(unit=*,fmt=*) cfgfile
    select case(mainchoice)
     case(1,3) ! gauge field translation
      write(unit=*,fmt=*) &
        "NOTE: The output filename template is ",trim(cfgfile),".trans"
     case(2) ! propagator translation
      write(unit=*,fmt=*) &
        "NOTE: The output filename template is ",trim(cfgfile),".trans"//trailer
     case default
      open(unit=8,file="TRANSLATOR.ERROR",action="write",status="replace", &
           form="formatted")
       write(unit=8,fmt=*) "mainchoice =", mainchoice
      close(unit=8,status="keep")
      stop
    end select
    write(unit=*,fmt=*) "Enter the minimum configuration number:"
    read(unit=*,fmt=*) cfgmin
    write(unit=*,fmt=*) "Enter the maximum configuration number:"
    read(unit=*,fmt=*) cfgmax
    write(unit=*,fmt=*) "Enter nx,ny,nz,nt as 4 numbers separated by spaces:"
    read(unit=*,fmt=*) nx, ny, nz, nt
    if (nx*ny*nz*nt>32*nxyztmax) then
     write(unit=*,fmt=*)
     write(unit=*,fmt=*) "ERROR: Edit translator.f90 to increase nxyztmax"
     write(unit=*,fmt=*)
     stop
    endif
    write(unit=*,fmt=*) &
          "Enter npx,npy,npz,npt (incoming) as 4 numbers separated by spaces:"
    read(unit=*,fmt=*) npin
    npsin = npin(1)*npin(2)*npin(3)*npin(4)
    if (mainchoice<3) then
     write(unit=*,fmt=*) &
           "Enter npx,npy,npz,npt (outgoing) as 4 numbers separated by spaces:"
     read(unit=*,fmt=*) npout
    else
     npout(1) = 1
     npout(2) = 1
     npout(3) = 1
     npout(4) = 1
    endif
    npsout = npout(1)*npout(2)*npout(3)*npout(4)

! Set up the blocking structure (copied directly from subroutine seblk).
    if (mainchoice==3) then
     ibl = 0
     do it = 1,2
      do iz = 1,2
       do iy = 1,2
        do ix = 1,2
         ibl = ibl + 1
         iblv(ibl,1) = ibl +    3-2*ix
         iblv(ibl,2) = ibl + 2*(3-2*iy)
         iblv(ibl,3) = ibl + 4*(3-2*iz)
         iblv(ibl,4) = ibl + 8*(3-2*it)
        enddo ! ix
       enddo ! iy
      enddo ! iz
     enddo ! it
    endif

! PART ONE: GAUGE FIELDS.
    select case(mainchoice)
     case(1,3) ! gauge field translation
      do icfg = cfgmin,cfgmax
       write(unit=cfgfile(2:6),fmt="(i5.5)") icfg
       nxin = nx/(4*npin(1))
       nyin = ny/(2*npin(2))
       nzin = nz/(2*npin(3))
       ntin = nt/(2*npin(4))
       nxout = nx/(4*npout(1))
       nyout = ny/(2*npout(2))
       nzout = nz/(2*npout(3))
       ntout = nt/(2*npout(4))
       do myid = 0,npsin-1
        write(unit=cfgfile(7:9),fmt="(i3.3)") myid
        ip(1) = modulo(myid,npin(1))
        ip(2) = modulo(myid/npin(1),npin(2))
        ip(3) = modulo(myid/(npin(1)*npin(2)),npin(3))
        ip(4) = modulo(myid/(npin(1)*npin(2)*npin(3)),npin(4))
        ipoffset = nxin*ip(1) + nyin*ip(2)*nx/4 + nzin*ip(3)*nx*ny/8 &
                 + ntin*ip(4)*nx*ny*nz/16
        open(unit=9,file=trim(rwdir)//trim(cfgfile),action="read", &
             position="rewind",status="old",form="unformatted")
         do icri = 1,12
          do it = 1,ntin
           do iz = 1,nzin
            do iy = 1,nyin
             do ix = 1,nxin
              isite = ipoffset + ix + (iy-1)*nx/4 + (iz-1)*nx*ny/8 &
                    + (it-1)*nx*ny*nz/16
              do ibl = 1,16
               read(unit=9)(u(icri,isite,mu,1,ibl),u(icri,isite,mu,2,ibl), &
                            mu=1,4)
              enddo ! ibl
             enddo ! ix
            enddo ! iy
           enddo ! iz
          enddo ! it
         enddo ! icri
         read(unit=9) gaction, beta, alat, tad, icfgsave, myidin
        close(unit=9,status="keep")
       enddo ! myid
       if (mainchoice==1) then
        do myid = 0,npsout-1
         write(unit=cfgfile(7:9),fmt="(i3.3)") myid
         ip(1) = modulo(myid,npout(1))
         ip(2) = modulo(myid/npout(1),npout(2))
         ip(3) = modulo(myid/(npout(1)*npout(2)),npout(3))
         ip(4) = modulo(myid/(npout(1)*npout(2)*npout(3)),npout(4))
         ipoffset = nxout*ip(1) + nyout*ip(2)*nx/4 + nzout*ip(3)*nx*ny/8 &
                  + ntout*ip(4)*nx*ny*nz/16
         open(unit=9,file=trim(rwdir)//trim(cfgfile)//".trans",action="write", &
              status="new",form="unformatted")
          do icri = 1,12
           do it = 1,ntout
            do iz = 1,nzout
             do iy = 1,nyout
              do ix = 1,nxout
               isite = ipoffset + ix + (iy-1)*nx/4 + (iz-1)*nx*ny/8 &
                     + (it-1)*nx*ny*nz/16
               do ibl = 1,16
                write(unit=9)(u(icri,isite,mu,1,ibl),u(icri,isite,mu,2,ibl), &
                              mu=1,4)
               enddo ! ibl
              enddo ! ix
             enddo ! iy
            enddo ! iz
           enddo ! it
          enddo ! icri
          write(unit=9) gaction, beta, alat, tad, icfgsave, myid
         close(unit=9,status="keep")
        enddo ! myid
       elseif (mainchoice==3) then
        open(unit=9,file=trim(rwdir)//trim(cfgfile)//".trans",action="write", &
             status="new",form="formatted")
         do mu = 1,4
          do it = 1,nt
           intt = (it-1)/2
           do iz = 1,nz
            intz = (iz-1)/2
            do iy = 1,ny
             inty = (iy-1)/2
             do ix = 1,nx
              intx = (ix-1)/4
              isite = 1 + intx + inty*nx/4 + intz*nx*ny/8 + intt*nx*ny*nz/16
              ieo = 1
              ibl = 1
              if (modulo(ix,4)==3.or.modulo(ix,4)==0) ieo = 3 - ieo
              if (modulo(iy,4)==3.or.modulo(iy,4)==0) ieo = 3 - ieo
              if (modulo(iz,4)==3.or.modulo(iz,4)==0) ieo = 3 - ieo
              if (modulo(it,4)==3.or.modulo(it,4)==0) ieo = 3 - ieo
              if (modulo(ix,2)==0) ibl = iblv(ibl,1)
              if (modulo(iy,2)==0) ibl = iblv(ibl,2)
              if (modulo(iz,2)==0) ibl = iblv(ibl,3)
              if (modulo(it,2)==0) ibl = iblv(ibl,4)
              write(unit=9,fmt="(12es12.5)") &
                    u(1,isite,mu,ieo,ibl),u( 7,isite,mu,ieo,ibl), &
                    u(2,isite,mu,ieo,ibl),u( 8,isite,mu,ieo,ibl), &
                    u(3,isite,mu,ieo,ibl),u( 9,isite,mu,ieo,ibl), &
                    u(4,isite,mu,ieo,ibl),u(10,isite,mu,ieo,ibl), &
                    u(5,isite,mu,ieo,ibl),u(11,isite,mu,ieo,ibl), &
                    u(6,isite,mu,ieo,ibl),u(12,isite,mu,ieo,ibl) 
             enddo ! ix
            enddo ! iy
           enddo ! iz
          enddo ! it
         enddo ! mu
         write(unit=9,fmt="(i6)") myid
         write(unit=9,fmt="(2es12.5)") tad(1), tad(4)
        close(unit=9,status="keep")
       endif
      enddo ! icfg
! PART TWO: QUARK PROPAGATORS.
     case(2) ! propagator translation
      do icfg = cfgmin,cfgmax
       write(unit=cfgfile(2:6),fmt="(i5.5)") icfg
       nxin = nx/(4*npin(1))
       nyin = ny/(2*npin(2))
       nzin = nz/(2*npin(3))
       ntin = nt/(2*npin(4))
       nxout = nx/(4*npout(1))
       nyout = ny/(2*npout(2))
       nzout = nz/(2*npout(3))
       ntout = nt/(2*npout(4))
       do myid = 0,npsin-1
        write(unit=cfgfile(7:9),fmt="(i3.3)") myid
        ip(1) = modulo(myid,npin(1))
        ip(2) = modulo(myid/npin(1),npin(2))
        ip(3) = modulo(myid/(npin(1)*npin(2)),npin(3))
        ip(4) = modulo(myid/(npin(1)*npin(2)*npin(3)),npin(4))
        ipoffset = nxin*ip(1) + nyin*ip(2)*nx/4 + nzin*ip(3)*nx*ny/8 &
                 + ntin*ip(4)*nx*ny*nz/16
        open(unit=9,file=trim(rwdir)//trim(cfgfile)//trailer,action="read", &
             position="rewind",status="old",form="unformatted")
         do it = 1,ntin
          do iz = 1,nzin
           do iy = 1,nyin
            do ix = 1,nxin
             isite = ipoffset + ix + (iy-1)*nx/4 + (iz-1)*nx*ny/8 &
                   + (it-1)*nx*ny*nz/16
             do ibl = 1,16
              ivbl = vecblinv(1,ibl)
              ibleo = vecblinv(2,ibl)
              jbleo = ibleo + 8
              do icri = 1,6
               if (ivbl==1) then
                read(unit=9) (u(icri,isite,id,1,ibleo), &
                              u(icri,isite,id,2,ibleo),id=1,4)
               else
                read(unit=9) (u(icri,isite,id,1,jbleo), &
                              u(icri,isite,id,2,jbleo),id=1,4)
               endif
              enddo ! icri
             enddo ! ibl
            enddo ! ix
           enddo ! iy
          enddo ! iz
         enddo ! it
         read(unit=9) kappa, icsrcin, idsrcin, icfgin, myidin
        close(unit=9,status="keep")
       enddo ! myid
       do myid = 0,npsout-1
        write(unit=cfgfile(7:9),fmt="(i3.3)") myid
        ip(1) = modulo(myid,npout(1))
        ip(2) = modulo(myid/npout(1),npout(2))
        ip(3) = modulo(myid/(npout(1)*npout(2)),npout(3))
        ip(4) = modulo(myid/(npout(1)*npout(2)*npout(3)),npout(4))
        ipoffset = nxout*ip(1) + nyout*ip(2)*nx/4 + nzout*ip(3)*nx*ny/8 &
                 + ntout*ip(4)*nx*ny*nz/16
        open(unit=9,file=trim(rwdir)//trim(cfgfile)//".trans"//trailer, &
             action="write",status="new",form="unformatted")
         do it = 1,ntout
          do iz = 1,nzout
           do iy = 1,nyout
            do ix = 1,nxout
             isite = ipoffset + ix + (iy-1)*nx/4 + (iz-1)*nx*ny/8 &
                   + (it-1)*nx*ny*nz/16
             do ibl = 1,16
              ivbl = vecblinv(1,ibl)
              ibleo = vecblinv(2,ibl)
              jbleo = ibleo + 8
              do icri = 1,6
               if (ivbl==1) then
                write(unit=9) (u(icri,isite,id,1,ibleo), &
                               u(icri,isite,id,2,ibleo),id=1,4)
               else
                write(unit=9) (u(icri,isite,id,1,jbleo), &
                               u(icri,isite,id,2,jbleo),id=1,4)
               endif
              enddo ! icri
             enddo ! ibl
            enddo ! ix
           enddo ! iy
          enddo ! iz
         enddo ! it
         write(unit=9) kappa, icsrcin, idsrcin, icfgin, myid
        close(unit=9,status="keep")
       enddo ! myid
      enddo ! icfg
     case default
      open(unit=8,file="TRANSLATOR.ERROR",action="write",status="replace", &
           form="formatted")
       write(unit=8,fmt=*) "mainchoice =", mainchoice
      close(unit=8,status="keep")
      stop
    end select

 end program translator

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
