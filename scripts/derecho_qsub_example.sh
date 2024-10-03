#!/bin/bash
#PBS -N AOS
#PBS -W block=False
#PBS -A NEOL0011
#PBS -l walltime=1:30:00
#PBS -q main@desched1
#PBS -l job_priority=premium
#PBS -r n
#PBS -j n
#PBS -o /glade/derecho/scratch/kmanning/AOS/KWM_Supercell_100m/d.out
#PBS -e /glade/derecho/scratch/kmanning/AOS/KWM_Supercell_100m/d.err
#PBS -l select=8:ncpus=4:mpiprocs=4
#PBS -l place=scatter:exclhost
module --force purge
module load ncarenv
module load intel
module load cray-mpich
module load mkl
module load netcdf-mpi
module load ncarcompilers

export TMPDIR=/glade/derecho/scratch/kmanning/temp
mkdir -p ${TMPDIR}
ulimit -s unlimited
cd /glade/derecho/scratch/kmanning/AOS/KWM_Supercell_100m/

# mpiexec --cpu-bind depth -n 4 /glade/derecho/scratch/kmanning/AOS/APAR-Observing-Simulator/code/embed-crsim/a.out
mpirun -n 32 /glade/derecho/scratch/kmanning/AOS/APAR-Observing-Simulator/code/embed-crsim/a.out namelist
