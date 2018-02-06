# BS removed qqcddata = /home/barals/andyfiles/scratch/cfgfiles/compileout/test
qqcdhome = /home/YOURKODIAKUSER/BaylorLQCD/qqcd-test
qqcddata = /home/YOURKODIAKUSER/BaylorLQCD/qqcd-test/sbin/objs
##BSqqcddata = /data/barals/qqcd/scratch/cfgfiles/compileout/test
#cc = /home/qcd2/qqcd-test/sbin/mpif90
#cc = mpifort   # for scali mpi 
cc = mpif90     #for openmpi
#cxxflags = -132 -w95 -O3 -shared-intel -march=core2 -falign-functions -mcmodel=large -unroll
#BScxxflags =-warn -132 -w95 -O3 -shared-intel -falign-functions -mcmodel=medium -unroll -check -g  -traceback 
#cxxflags =-warn -132 -w95 -O3 -shared-intel -falign-functions -mcmodel=medium -unroll -check bounds -g  -traceback 
#cxxflags = -132 -w95 -O3 -shared-intel -falign-functions -mcmodel=medium -unroll -check bounds -g  -traceback 
cxxflags = -132 -w95 -O3 -shared-intel -falign-functions -mcmodel=medium -unroll -g  -traceback
##cxxflags =-warn -132 -w95  -shared-intel -falign-functions -mcmodel=medium -unroll -g  -traceback

# BS previous line
#cxxflags = -Wall -ffree-line-length-none -O3  -falign-functions -mcmodel=medium -unroll -g -fbacktrace

target = qqcd-testout
#target = qqcd-testout#BS

objects = $(qqcddata)/kinds.o \
	$(qqcddata)/latdims.o \
	$(qqcddata)/basics.o        $(qqcddata)/lagfib.o \
	$(qqcddata)/lattice.o  $(qqcddata)/gaugetools.o \
	$(qqcddata)/heatbath.o \
	$(qqcddata)/commonkind.o \
	$(qqcddata)/gaugelinks.o \
	$(qqcddata)/input1.o \
	$(qqcddata)/input2.o \
	$(qqcddata)/input5a.o \
	$(qqcddata)/input5b.o \
	$(qqcddata)/operator.o \
	$(qqcddata)/seed.o \
	$(qqcddata)/sub.o \
	$(qqcddata)/multstorage.o \
	$(qqcddata)/diracops.o \
	$(qqcddata)/debug.o  $(qqcddata)/pseudolapack.o \
	$(qqcddata)/shift.o  $(qqcddata)/gmresrhs.o \
	$(qqcddata)/inverters.o \
	$(qqcddata)/working.o \
	$(qqcddata)/quark.o \
	$(qqcddata)/vevbleft.o \
	$(qqcddata)/analysis.o \
	$(qqcddata)/printops.o \
	$(qqcddata)/disconloops.o \
	$(qqcddata)/cfgsprops.o  \
	$(qqcddata)/cfgspropsmain.o
#$(target): $(objects)
#	$(cc) $(cxxflags) -o $(target) $(objects) /home/abdou/lapack-3.2.2/lapack_LINUX.a   /home/abdou/lapack-3.2.2/blas_LINUX.a
#	mv $(target) ../bin

$(target): $(objects)
	$(cc) $(cxxflags) -o $(target) $(objects) -lblas -llapack
	mv $(target) ../bin

$(qqcddata)/kinds.o: $(qqcdhome)/user/kinds.f90
	$(cc) $(cxxflags) -c $(qqcdhome)/user/kinds.f90 -o $(qqcddata)/kinds.o 

$(qqcddata)/latdims.o: $(qqcdhome)/user/latdims.f90
	$(cc) $(cxxflags) -c $(qqcdhome)/user/latdims.f90 -o $(qqcddata)/latdims.o

$(qqcddata)/basics.o: $(qqcdhome)/cfgsprops/basics.f90
	$(cc) $(cxxflags) -c $(qqcdhome)/cfgsprops/basics.f90 -o $(qqcddata)/basics.o

$(qqcddata)/lagfib.o: $(qqcdhome)/cfgsprops/glue/lagfib.f90
	$(cc) $(cxxflags) -c $(qqcdhome)/cfgsprops/glue/lagfib.f90 -o $(qqcddata)/lagfib.o

$(qqcddata)/lattice.o: $(qqcdhome)/cfgsprops/glue/lattice.f90
	$(cc) $(cxxflags) -c $(qqcdhome)/cfgsprops/glue/lattice.f90 -o $(qqcddata)/lattice.o

$(qqcddata)/gaugetools.o: $(qqcdhome)/cfgsprops/glue/gaugetools.f90
	$(cc) $(cxxflags) -c $(qqcdhome)/cfgsprops/glue/gaugetools.f90 -o $(qqcddata)/gaugetools.o

$(qqcddata)/heatbath.o: $(qqcdhome)/cfgsprops/glue/heatbath.f90
	$(cc) $(cxxflags) -c $(qqcdhome)/cfgsprops/glue/heatbath.f90 -o $(qqcddata)/heatbath.o

$(qqcddata)/gaugelinks.o: $(qqcdhome)/cfgsprops/quark/common/gaugelinks.f90
	$(cc) $(cxxflags) -c $(qqcdhome)/cfgsprops/quark/common/gaugelinks.f90 -o $(qqcddata)/gaugelinks.o

$(qqcddata)/input1.o: $(qqcdhome)/cfgsprops/quark/common/input1.f90
	$(cc) $(cxxflags) -c $(qqcdhome)/cfgsprops/quark/common/input1.f90 -o $(qqcddata)/input1.o

$(qqcddata)/input2.o: $(qqcdhome)/cfgsprops/quark/common/input2.f90
	$(cc) $(cxxflags) -c $(qqcdhome)/cfgsprops/quark/common/input2.f90 -o $(qqcddata)/input2.o

$(qqcddata)/input5a.o: $(qqcdhome)/cfgsprops/quark/common/input5a.f90
	$(cc) $(cxxflags) -c $(qqcdhome)/cfgsprops/quark/common/input5a.f90 -o $(qqcddata)/input5a.o

$(qqcddata)/input5b.o: $(qqcdhome)/cfgsprops/quark/common/input5b.f90
	$(cc) $(cxxflags) -c $(qqcdhome)/cfgsprops/quark/common/input5b.f90 -o $(qqcddata)/input5b.o

$(qqcddata)/operator.o: $(qqcdhome)/cfgsprops/quark/common/operator.f90
	$(cc) $(cxxflags) -c $(qqcdhome)/cfgsprops/quark/common/operator.f90 -o $(qqcddata)/operator.o

$(qqcddata)/seed.o: $(qqcdhome)/cfgsprops/quark/common/seed.f90
	$(cc) $(cxxflags) -c $(qqcdhome)/cfgsprops/quark/common/seed.f90 -o $(qqcddata)/seed.o

$(qqcddata)/sub.o: $(qqcdhome)/cfgsprops/quark/common/sub.f90
	$(cc) $(cxxflags) -c $(qqcdhome)/cfgsprops/quark/common/sub.f90 -o $(qqcddata)/sub.o

$(qqcddata)/multstorage.o: $(qqcdhome)/cfgsprops/quark/common/multstorage.f90
	$(cc) $(cxxflags) -c $(qqcdhome)/cfgsprops/quark/common/multstorage.f90 -o $(qqcddata)/multstorage.o

$(qqcddata)/commonkind.o: $(qqcdhome)/cfgsprops/quark/common/commonkind.f90
	$(cc) $(cxxflags) -c $(qqcdhome)/cfgsprops/quark/common/commonkind.f90 -o $(qqcddata)/commonkind.o

$(qqcddata)/debug.o: $(qqcdhome)/cfgsprops/debug/debug.f90
	$(cc) $(cxxflags) -c $(qqcdhome)/cfgsprops/debug/debug.f90 -o $(qqcddata)/debug.o

$(qqcddata)/pseudolapack.o: $(qqcdhome)/cfgsprops/quark/pseudolapack.f90
	$(cc) $(cxxflags) -c $(qqcdhome)/cfgsprops/quark/pseudolapack.f90 -o $(qqcddata)/pseudolapack.o

$(qqcddata)/shift.o: $(qqcdhome)/cfgsprops/quark/shift.f90
	$(cc) $(cxxflags) -c $(qqcdhome)/cfgsprops/quark/shift.f90 -o $(qqcddata)/shift.o

$(qqcddata)/gmresrhs.o: $(qqcdhome)/cfgsprops/quark/gmresrhs.f90
	$(cc) $(cxxflags) -c $(qqcdhome)/cfgsprops/quark/gmresrhs.f90 -o $(qqcddata)/gmresrhs.o

$(qqcddata)/mylapack.o: $(qqcdhome)/cfgsprops/quark/mylapack.f
	$(cc) $(cxxflags) -c $(qqcdhome)/cfgsprops/quark/mylapack.f -o $(qqcddata)/mylapack.o

$(qqcddata)/inverters.o: $(qqcdhome)/cfgsprops/quark/inverters.f90
	$(cc) $(cxxflags) -c $(qqcdhome)/cfgsprops/quark/inverters.f90 -o $(qqcddata)/inverters.o

$(qqcddata)/working.o: $(qqcdhome)/cfgsprops/quark/working.f90
	$(cc) $(cxxflags) -c $(qqcdhome)/cfgsprops/quark/working.f90 -o $(qqcddata)/working.o

$(qqcddata)/vevbleft.o: $(qqcdhome)/cfgsprops/quark/vevbleft.f90
	$(cc) $(cxxflags) -c $(qqcdhome)/cfgsprops/quark/vevbleft.f90 -o $(qqcddata)/vevbleft.o

$(qqcddata)/printops.o: $(qqcdhome)/cfgsprops/quark/printops.f90
	$(cc) $(cxxflags) -c $(qqcdhome)/cfgsprops/quark/printops.f90 -o $(qqcddata)/printops.o

$(qqcddata)/disconloops.o: $(qqcdhome)/cfgsprops/quark/disconloops.f90
	$(cc) $(cxxflags) -c $(qqcdhome)/cfgsprops/quark/disconloops.f90 -o $(qqcddata)/disconloops.o

$(qqcddata)/quark.o: $(qqcdhome)/cfgsprops/quark/quark.f90
	$(cc) $(cxxflags) -c $(qqcdhome)/cfgsprops/quark/quark.f90 -o $(qqcddata)/quark.o

$(qqcddata)/analysis.o: $(qqcdhome)/cfgsprops/quark/analysis.f90
	$(cc) $(cxxflags) -c $(qqcdhome)/cfgsprops/quark/analysis.f90 -o $(qqcddata)/analysis.o

$(qqcddata)/diracops.o: $(qqcdhome)/cfgsprops/quark/diracops.f90
	$(cc) $(cxxflags) -c $(qqcdhome)/cfgsprops/quark/diracops.f90 -o $(qqcddata)/diracops.o

$(qqcddata)/cfgsprops.o: $(qqcdhome)/cfgsprops/cfgsprops.f90
	$(cc) $(cxxflags) -c $(qqcdhome)/cfgsprops/cfgsprops.f90 -o $(qqcddata)/cfgsprops.o

$(qqcddata)/cfgspropsmain.o: $(qqcdhome)/user/cfgspropsmain.f90
	$(cc) $(cxxflags) -c $(qqcdhome)/user/cfgspropsmain.f90 -o $(qqcddata)/cfgspropsmain.o

clean:
	rm -f $(qqcddata)/*
	rm -f $(qqcdhome)/bin/$(target)
	rm -f *.mod
