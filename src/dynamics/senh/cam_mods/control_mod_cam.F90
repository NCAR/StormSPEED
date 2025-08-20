! This module contains additional CAM-specific settings beside the dycore/components/homme/src/share/control_mod.F90 module.
!
module control_mod_cam
  use control_mod
  implicit none
  logical, public :: refined_mesh
  integer, public :: multilevel
  integer, public :: tasknum
  integer, public :: remapfreq      ! remap frequency of synopsis of system state (steps)
end module control_mod_cam
