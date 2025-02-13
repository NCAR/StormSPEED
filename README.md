# StormSpeed - A CAM port of E3SM SE non-hydrostatic dycores

CAM E3SM NH development code is stored in this repository on the stormspeed branch.

This repository has been forked from ESCOMP cam and contains the following branches:
* main - contains this readme and the Code of Conduct information
* stormspeed - contains the development branch for the StormSPEED project which is based off CAM6_3_085

## How to checkout and use StormSPEED CAM on the Derecho:

Begin by cloning the StormSpeed repository and checking out the CAM code and externals on derecho
```
cd /glade/derecho/scratch/$USER
git clone https://github.com/NCAR/StormSPEED
cd StormSPEED
git checkout stormspeed
./manage_externals/checkout_externals
```
Create a new case, set the dycore target to theta-l, configure run, setup, build, and run a Kessler Simple physics test.

```
cd /glade/derecho/scratch/$USER
set CESMDIR=/glader/derecho/scratch/$USER/StormSPEED
$CESMDIR/cime/scripts/create_newcase --compset FKESSLER  --res ne16_g37 --compiler intel --case FKESSLER_StormSPEED.01 --driver mct --run-unsupported
cd FKESSLER_StormSPEED.01
./xmlchange CAM_TARGET=theta-l
./xmlchange STOP_OPTION=ndays,STOP_N=12,RESUBMIT=0
./xmlchange PROJECT=**ENTER USER_PROJECT_NUMBER HERE**
./xmlchange JOB_WALLCLOCK_TIME=00:40:00
./xmlchange DOUT_S='FALSE'
./xmlchange BFBFLAG=TRUE
./case.setup
cat >> user_nl_cam << EOF
 interpolate_output = .true.
 interpolate_nlat = 96
 interpolate_nlon = 144
EOF
./case.build
./case.submit
```

### To view the release branches in Github, go to the "Branch:main" pulldown menu and select the stormspeed branch.

## NOTE: This is **unsupported** development code and is subject to the [CESM developer's agreement](https://www.cgd.ucar.edu/sections/cseg/policies).

### CAM Documentation - https://ncar.github.io/CAM/doc/build/html/index.html

### CAM6 namelist settings - https://docs.cesm.ucar.edu/models/cesm2/settings/current/cam_nml.html

