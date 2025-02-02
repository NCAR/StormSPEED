# StormSpeed - A CAM port of E3SM SE non-hydrostatic dycores

CAM E3SM NH development code is stored in this repository on the stormspeed branch.

This repository has been forked from ESCOMP cam and contains the following branches:
* main - contains this readme and the Code of Conduct information
* stormspeed - contains the development branch for the StormSPEED project which is based off CAM6_3_085

## How to checkout and use StormSPEED CAM on the NCAR Izumi cluster:

Begin by cloning the StormSpeed repository and checking out the CAM code and externals
```
cd /scratch/cluster/$USER
git clone https://github.com/NCAR/StormSPEED
cd StormSPEED
git checkout stormspeed
./manage_externals/checkout_externals
```
Setup the Izumi environment for configuring and compiling StormSPEED.  For now, this will need to be done for each new terminal session.  We are working on making this work out of the box for supported hardware.

```
module load compiler/intel/20.0.1 lang/python/3.11.5
setenv Kokkos_ROOT "/project/amp/jet/install_izumi/spack/opt/spack/linux-almalinux8-skylake_avx512/gcc-9.3.0/kokkos-3.0.00-f6s53cvwkjph5gbxt7cadkq4psqnph22"
setenv GCC_ROOT "/project/amp/jet/install_izumi/spack/opt/spack/linux-almalinux8-skylake_avx512/gcc-9.3.0/gcc-9.3.0-w7rjyw3sdfhyxaih5huqx7q6tpxhnxmr"
setenv USER_SLIBS "$Kokkos_ROOT/lib64/libkokkoscontainers.so.3.0.0 $Kokkos_ROOT/lib64/libkokkoscore.so.3.0.0 $GCC_ROOT/lib64/libstdc++.so.6.0.28"
setenv LD_LIBRARY_PATH "${LD_LIBRARY_PATH}:${Kokkos_ROOT}/lib64:${GCC_ROOT}/lib64"
```
Create a new case, set the dycore target to theta-l, configure run, setup, build, and run a Kessler Simple physics test on Izumi.

```
cd /scratch/cluster/$USER
set CESMDIR=/scratch/cluster/$USER/StormSPEED
$CESMDIR/cime/scripts/create_newcase --compset FKESSLER  --res ne16_g37 --compiler intel --case FKESSLER_StormSPEED.01 --driver mct --walltime 24:00:00 --queue long --run-unsupported
cd FKESSLER_StormSPEED.01
./xmlchange CAM_TARGET=theta-l
./xmlchange STOP_OPTION=ndays,STOP_N=12,RESUBMIT=0
./xmlchange PROJECT=none
./xmlchange JOB_WALLCLOCK_TIME=24:00:00
./xmlchange DOUT_S='FALSE'
./xmlchange BFBFLAG=TRUE
./case.setup
./case.build
set RUNDIR=`./xmlquery -value RUNDIR`
mkdir -p $RUNDIR/timing/checkpoints
./case.submit
```

### To view the release branches in Github, go to the "Branch:main" pulldown menu and select the appropriate branch.

## NOTE: This is **unsupported** development code and is subject to the [CESM developer's agreement](https://www.cgd.ucar.edu/sections/cseg/policies).

### CAM Documentation - https://ncar.github.io/CAM/doc/build/html/index.html

### CAM6 namelist settings - https://docs.cesm.ucar.edu/models/cesm2/settings/current/cam_nml.html

