! This module contains additional CAM-specific settings beside the dycore/components/homme/src/share/element_mod.F90 module.
!
module element_mod_cam
  use element_mod

  implicit none
  private

  type, public :: element_t
     integer(kind=int_kind) :: LocalId       ! element numbering on each MPI task
     integer(kind=int_kind) :: GlobalId      ! global element numbering independent of MPI decomposition

     ! Coordinate values of element points
     type (spherical_polar_t) :: spherep(np,np)                       ! Spherical or physical coords of GLL points
     ! EVENTUALLY HERE spherical_polar_t should become cartesian3D_t
     ! However, this breaks a bunch of code that expects spherep to be spherical_polar_t
     ! So instead we will let lat=y, lon=x, r=z=0
     ! This is an ugly hack that we can fix once the "planar" geometry code can handle arbitrary geometries/topologies

     ! Equ-angular gnomonic projection coordinates
     type (cartesian2D_t)     :: cartp(np,np)                         ! gnomonic or reference coords of GLL points
     type (cartesian2D_t)     :: corners(4)                           ! gnomonic or reference coords of element corners
     real (kind=real_kind)    :: u2qmap(4,2)                          ! bilinear map from ref element to quad in cubedsphere coordinates

     ! 3D cartesian coordinates
     type (cartesian3D_t)     :: corners3D(4)                         ! Physical coords of corners


     ! Element diagnostics
     real (kind=real_kind)    :: area                                 ! Area of element
     real (kind=real_kind)    :: normDinv                             ! some type of norm of Dinv used for CFL
     real (kind=real_kind)    :: dx_short                             ! short length scale in km
     real (kind=real_kind)    :: dx_long                              ! long length scale in km

     real (kind=real_kind)    :: variable_hyperviscosity(np,np)       ! hyperviscosity based on above
     real (kind=real_kind)    :: hv_courant                           ! hyperviscosity courant number
     real (kind=real_kind)    :: tensorVisc(np,np,2,2)                !og, matrix V for tensor viscosity

     type (GridVertex_t)      :: vertex                               ! element grid vertex information
     type (EdgeDescriptor_t)  :: desc

     type (elem_state_t)      :: state
     type (derived_state_t)   :: derived
     type (elem_accum_t)       :: accum

     ! Metric terms
     real (kind=real_kind)    :: met(np,np,2,2)                       ! metric tensor on velocity and pressure grid
     real (kind=real_kind)    :: metinv(np,np,2,2)                    ! metric tensor on velocity and pressure grid
     real (kind=real_kind)    :: metdet(np,np)                        ! g = SQRT(det(g_ij)) on velocity and pressure grid
     real (kind=real_kind)    :: rmetdet(np,np)                       ! 1/metdet on velocity pressure grid
     real (kind=real_kind)    :: D(np,np,2,2)                         ! Map covariant field on cube or reference plane to vector field on the sphere or physical plane
     real (kind=real_kind)    :: Dinv(np,np,2,2)                      ! Map vector field on the sphere or physical plane to covariant v on cube or reference plane


     ! Mass flux across the sides of each sub-element.
     ! The storage is redundent since the mass across shared sides
     ! must be equal in magnitude and opposite in sign.
     ! The layout is like:
     !   --------------------------------------------------------------
     ! ^|    (1,4,3)     |                |              |    (4,4,3) |
     ! ||                |                |              |            |
     ! ||(1,4,4)         |                |              |(4,4,4)     |
     ! ||         (1,4,2)|                |              |     (4,4,2)|
     ! ||                |                |              |            |
     ! ||   (1,4,1)      |                |              |  (4,4,1)   |
     ! |---------------------------------------------------------------
     ! S|                |                |              |            |
     ! e|                |                |              |            |
     ! c|                |                |              |            |
     ! o|                |                |              |            |
     ! n|                |                |              |            |
     ! d|                |                |              |            |
     !  ---------------------------------------------------------------
     ! C|                |                |              |            |
     ! o|                |                |              |            |
     ! o|                |                |              |            |
     ! r|                |                |              |            |
     ! d|                |                |              |            |
     ! i|                |                |              |            |
     ! n---------------------------------------------------------------
     ! a|    (1,1,3)     |                |              |    (4,1,3) |
     ! t|                |                |              |(4,1,4)     |
     ! e|(1,1,4)         |                |              |            |
     !  |         (1,1,2)|                |              |     (4,1,2)|
     !  |                |                |              |            |
     !  |    (1,1,1)     |                |              |  (4,1,1)   |
     !  ---------------------------------------------------------------
     !          First Coordinate ------->


     ! Convert vector fields from spherical to rectangular components
     ! The transpose of this operation is its pseudoinverse.
     ! This is just "identity" for plane

     real (kind=real_kind)    :: vec_sphere2cart(np,np,3,2)

     ! Mass matrix terms for an element on a cube or reference face
     real (kind=real_kind)    :: mp(np,np)                            ! mass matrix on v and p grid
     real (kind=real_kind)    :: rmp(np,np)                           ! inverse mass matrix on v and p grid

     ! Mass matrix terms for an element on the sphere
     ! This mass matrix is used when solving the equations in weak form
     ! with the natural (surface area of the sphere or plane) inner product
     real (kind=real_kind)    :: spheremp(np,np)                      ! mass matrix on v and p grid
     real (kind=real_kind)    :: rspheremp(np,np)                     ! inverse mass matrix on v and p grid

     integer(kind=long_kind)  :: gdofP(np,np)                         ! global degree of freedom (P-grid)

     real (kind=real_kind)    :: fcor(np,np)                          ! Coriolis term

     type (index_t) :: idxP
     integer :: FaceNum

     ! force element_t to be a multiple of 8 bytes.
     ! on BGP, code will crash (signal 7, or signal 15) if 8 byte alignment is off
     ! check core file for:
     ! core.63:Generated by interrupt..(Alignment Exception DEAR=0xa1ef671c ESR=0x01800000 CCR0=0x4800a002)
     !integer :: dummy
  end type element_t

end module element_mod_cam
