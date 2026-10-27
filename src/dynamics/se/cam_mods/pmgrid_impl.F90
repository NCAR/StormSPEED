module pmgrid_impl
  use edgetype_mod, only : EdgeBuffer_t
  use dimensions_mod_cam, only : nnodes,npart,nmpi_per_node
  use reduction_mod, only : ReductionBuffer_ordered_1d_t
  implicit none
  save
  integer, parameter :: plon   = 1
  integer, parameter :: plev   = PLEV
  integer, parameter :: plat   = 1
  integer, parameter :: plevp  = plev + 1
  type (EdgeBuffer_t) :: edge_3
  type (EdgeBuffer_t) :: edge_1lp1
  type (EdgeBuffer_t) :: edge_3lp1
  type (ReductionBuffer_ordered_1d_t) :: red
  integer, parameter :: nxpt   = 0
  integer, parameter :: jintmx = 0
  integer, parameter :: i1     = 1
  integer, parameter :: j1     = 1
  integer, parameter :: numbnd = 0
  integer :: beglat, endlat, beglatex, endlatex, numlats
  logical :: dyndecomp_set = .false.
  parameter (beglat = 1, endlat = plat, numlats = plat)
  integer, parameter :: iam = 0
  logical, parameter :: masterproc = .true.
end module pmgrid_impl
