#!/bin/sh
#PBS -l nodes=36:ppn=4
###PBS -l nodes=2:ppn=4,walltime=00:30:00

#PBS -m ea -M YOUR_EMAIL@baylor.edu

module add mvapich2-intel/2.2-2

cd $PBS_O_WORKDIR

num=`cat $PBS_NODEFILE | wc -l`
echo "Requested processors: $num"
echo "Node(s):"
uniq $PBS_NODEFILE
echo

echo "Cleaning scratch files"
rm /data/YOURKODIAKUSER/qqcd/scratch/DUBLIN.OUT
rm /data/YOURKODIAKUSER/qqcd/scratch/CFGSPROPS.LOG
rm /data/YOURKODIAKUSER/qqcd/scratch/EIGEN_VALS.LOG
rm /data/YOURKODIAKUSER/qqcd/scratch/trueresidual.dat
rm /data/YOURKODIAKUSER/qqcd/scratch/residual.dat
rm /data/YOURKODIAKUSER/qqcd/scratch/eigresidual.dat
rm /data/YOURKODIAKUSER/qqcd/scratch/wresidual.dat
rm /data/YOURKODIAKUSER/qqcd/scratch/linresidual.dat
rm /data/YOURKODIAKUSER/qqcd/scratch/rnherm.dat


echo

echo "Job starting at `date`"
echo

START=`date '+%s'`
export MV2_ENABLE_AFFINITY=0
mpiexec -n $num -machinefile $PBS_NODEFILE ./qqcd-testout >quick.out
OMP_NUM_THREADS=1
#/qqcd-testout1
END=`date '+%s'`

echo
echo "Job finished at `date`"
echo
echo 'Total Execution time: '`expr $END - $START`' seconds'
