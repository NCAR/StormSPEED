# versions and sizes from here: https://hub.docker.com/r/jiansunncar/ubuntu24.04_gcc12.4.0_stormspeed
FROM jiansunncar/ubuntu24.04_gcc12.4.0_stormspeed:v0

# Copy the source code from a pull request into the container
COPY . /tmp/StormSPEED

# Set up some git configurations; it can be any user name and email address
RUN git config --global user.email "example_user@example.com" \
    && git config --global user.name "example_user"

# # Check out individual components of StormSPEED
# RUN cd /tmp/StormSPEED \
#     && ./bin/git-fleximod update

# Set up some case-independent environment variables
ENV TMP=tmp
ENV USER=example_user

# build an FADIAB test at 1-degree resolution
RUN . /home/spack/share/spack/setup-env.sh \
    && spack find \
    && spack load gcc@12.4.0/jwf2kjs \
    && spack load libxml2 cmake openmpi netcdf-c netcdf-fortran parallel-netcdf parallelio esmf
    # && export MPI_ROOT=$(spack location -i openmpi@5.0.8 ^gcc@12.4.0) \
    # && export NETCDF_C_PATH=$(spack location -i netcdf-c@4.9.2) \
    # && export NETCDF_FORTRAN_PATH=$(spack location -i netcdf-fortran@4.6.1) \
    # && export PNETCDF=$(spack location -i parallel-netcdf@1.14.0) \
    # && export PIO=$(spack location -i parallelio@2.6.5) \
    # && export ESMFMKFILE=$(spack location -i esmf@8.8.1)/lib/esmf.mk \
    # && export LAPACK=$(spack location -i netlib-lapack@3.12.1) \
    # && export PIO_VERSION_MAJOR=2 \
    # && export PIO_TYPENAME_VALID_VALUES="netcdf, pnetcdf, netcdf4c, netcdf4p" \
    # && cd /tmp/StormSPEED/cime/scripts \
    # && ./create_newcase --case /$TMP/ci_test --machine cirrus --compset FADIAB --res ne30_ne30_mg17 --compiler gnu --run-unsupported \
    # && cd /$TMP/ci_test \
    # && ./case.setup \
    # && ./case.build

WORKDIR /$TMP/ci_test