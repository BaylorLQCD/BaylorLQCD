 module gaugefilewrite

   implicit none

subroutine getfilename(filename, config,configwrite)
   character(len=*), intent(out)   :: filename
   integer, intent(in)             :: config, configwrite

   integer :: whichset

   ! whichset will be incremented every 10th configuration, i.e.,
   ! configs 1-10 have whichset 1, configs, 11-20 have whichset 2, etc.

   whichset = ((config - 1) / configwrite) + 1

   filename = "cfgsprops-xxxx.log"
   write(unit=filename(11:14), fmt="(i3.3)") whichset

end subroutine getfilename

subroutine gettimestamp(str)

   character(len=*), intent(out) :: str

   integer(kind=4), dimension(8) :: dt
   character(len=10)             :: igndate, igntime, ignzone


   call date_and_time(igndate, igntime, ignzone, dt)

   write(unit=str,fmt="(i2.2,'-',i2.2,'-',i4,'-',i2.2,':',i2.2,':',i2.2)") dt(2), dt(3), dt(1), dt(5), dt(6), dt(7)


 end module gaugefilewrite
