# StormSPEED - The 'Storm-resolving SPEctral Element Dycore' for CESM3

[StormSPEED](https://sites.google.com/umich.edu/nsf-stormspeed/home), funded by the U.S. National Science Foundation (NSF), will integrate the DoE's nonhydrostatic (NH) version of the Spectral Element (SE) dynamical core and performant Semi-Lagrangian tracer advection scheme into NCAR's Community Atmosphere Model (CAM). These new features will enable nonhydrostatic, computationally-efficient, kilometer-scale Earth system simulations with the Community Earth System Model (CESM) developed at NCAR.

## Project Goals

- Provide a significant enhancement to CESM which will enable kilometer-scale Earth system simulations
- Innovate CESM's computational design
- Educate the next generation of Earth system scientists

## Key Features

- Highly optimized and efficient dycore design
- CPU/GPU implementation making use of the Kokkos library.
- Reuse of the infrastructure and dynamics/physics interface coupling layer.
- Full integration into the CIME framework to provide seamless access to the new functionality

## Documentation and Support

For information on installation, usage, development guidelines, and troubleshooting, please visit the [StormSPEED Wiki](https://github.com/NCAR/StormSPEED/wiki).

[CAM Documentation](https://ncar.github.io/CAM/doc/build/html/index.html) - https://ncar.github.io/CAM/doc/build/html/index.html

[CAM6 namelist settings](https://docs.cesm.ucar.edu/models/cesm2/settings/current/cam_nml.html)  - https://docs.cesm.ucar.edu/models/cesm2/settings/current/cam_nml.html

[StormSPEED webpage](https://sites.google.com/umich.edu/nsf-stormspeed)

## License

StormSPEED is released under the [MIT License](LICENSE).

---

For all other details, including instructions for obtaining and running the code, please see the [StormSPEED Wiki](https://github.com/NCAR/StormSPEED/wiki).

**NOTE: This is **unsupported** development code and is subject to the [CESM developer's agreement](https://www.cgd.ucar.edu/sections/cseg/policies).**
