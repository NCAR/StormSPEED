module dyn_grid
!-------------------------------------------------------------------------------
!
! Define SE computational grids on the dynamics decomposition.
!

! The grid used by the SE dynamics is called the GLL grid.  It is
! decomposed into elements which correspond to "blocks" in the
! physics/dynamics coupler terminology.  The columns in this grid are
! located at the Gauss-Lobatto-Legendre (GLL) quadrature points.  The GLL
! grid will also be used by the physics if the CSLAM advection is not used.
! If CSLAM is used for tracer advection then it uses an FVM grid and the
! physics will either use the same FVM grid or an FVM grid with a different
! number of equal area subcells.  The FVM grid used by the physics is
! referred to as the "physgrid".
!
! Module responsibilities:
!
! . Provide the physics/dynamics coupler (in module phys_grid) with data for the
!   physics grid on the dynamics decomposition.
!
! . Create CAM grid objects that are used by the I/O functionality to read
!   data from an unstructured grid format to the dynamics data structures, and
!   to write from the dynamics data structures to unstructured grid format.  The
!   global column ordering for the unstructured grid is determined by the SE dycore.
!
!-------------------------------------------------------------------------------

use shr_kind_mod,           only: r8 => shr_kind_r8, shr_kind_cl
use spmd_utils,             only: masterproc, iam, mpicom, mstrid=>masterprocid
use spmd_utils,             only: npes, mpi_integer, mpi_real8, mpi_success, mpi_max
use constituents,           only: pcnst
use physconst,              only: pi
use cam_initfiles,          only: initial_file_get_id
use cam_grid_support,       only: iMap
use physics_column_type,    only: physics_column_t

use cam_logfile,            only: iulog
use cam_abortutils,         only: endrun

use pio,                    only: pio_inq_dimid, pio_seterrorhandling
use pio,                    only: PIO_BCAST_ERROR, PIO_NOERR, file_desc_t
use pio,                    only: pio_internal_error, pio_inq_dimlen

use domain_mod,             only: domain1d_t

use dimensions_mod_cam,     only: globaluniquecols, nelem, nelemd, nelemdmax, &
                                  ne, np, npsq, fv_nphys, nlev, ntrac
use element_mod,            only: element_t
use hybvcoord_mod,          only: hvcoord_t, set_layer_locations
use prim_driver_mod,        only: prim_init1
use time_mod,               only: TimeLevel_t

implicit none
private
save

integer, parameter :: dyn_decomp = 101 ! The SE dynamics grid
!jtinteger, parameter :: fvm_decomp = 102 ! The FVM (CSLAM) grid
integer, parameter :: physgrid_d = 103 ! physics grid on dynamics decomp
integer, parameter :: ini_decomp = 104 ! alternate dynamics grid for reading initial file

character(len=3), protected :: ini_grid_name

! Name of horizontal grid dimension in initial file.
character(len=6), protected :: ini_grid_hdim_name = ''


integer, parameter :: ptimelevels = 2

type (TimeLevel_t)         :: TimeLevel     ! main time level struct (used by tracers)
type (domain1d_t), pointer :: dom_mt(:) => null()
type (hvcoord_t)           :: hvcoord
type(element_t),   pointer :: elem(:) => null()
real(r8),          pointer :: w(:) => null()        ! weights

! FV physics grid resolution (physics on GLL grid if NPG=0)
public ::  dyn_decomp
public ::  ini_grid_name
public ::  ini_grid_hdim_name
public ::  ptimelevels
public ::  timelevel
public ::  dom_mt
public ::  hvcoord
public ::  elem
public ::  w

public :: dyn_grid_init
public :: get_dyn_grid_info    ! Return physics grid column information
public :: physgrid_copy_attributes_d

public :: get_horiz_grid_d
public :: get_horiz_grid_dim_d
public :: get_dyn_grid_parm
public :: get_dyn_grid_parm_real1d
public :: dyn_grid_get_elem_coords
public :: dyn_grid_get_colndx
public :: fv_physgrid_init, fv_physgrid_final

! number of global dynamics columns. Set by SE dycore init.
integer :: ngcols_d = 0     ! number of dynamics columns
! number of global elements. Set by SE dycore init.
integer :: nelem_d = 0

real(r8), parameter :: rad2deg = 180.0_r8 / pi

type(physics_column_t), allocatable, target :: local_dyn_columns(:)

!=========================================================================================
contains
!=========================================================================================

subroutine dyn_grid_init()

   ! Initialize SE grid, and decomposition.

   use hycoef,         only: hycoef_init, hypi, hypm, nprlev, &
                             hyam, hybm, hyai, hybi, ps0
   use ref_pres,       only: ref_pres_init
   use time_manager,   only: get_nstep, get_step_size

   use thread_mod_cam,   only: nthreads
   use parallel_mod_cam, only: par
   use control_mod_cam,  only: qsplit, rsplit, dt_remap_factor, dt_tracer_factor, &
                             timestep_make_eam_parameters_consistent
   use time_mod,       only: tstep, nsplit

   ! Local variables
   type(file_desc_t), pointer :: fh_ini

#ifdef _OPENMP
   integer :: omp_get_num_threads
#endif

   integer :: k
   integer :: ncolid
   integer :: ncollen
   integer :: ierr
   integer :: neltmp(3)
   integer :: dtime
   integer :: nstep_factor

   character(len=*), parameter :: sub = 'dyn_grid_init'
   !----------------------------------------------------------------------------

   dtime = get_step_size()

   ! Get file handle for initial file and first consistency check
   fh_ini => initial_file_get_id()

   call pio_seterrorhandling(fh_ini, pio_bcast_error)
   ierr = pio_inq_dimid(fh_ini, 'ncol', ncolid)
   call pio_seterrorhandling(fh_ini, pio_internal_error)

   if (ierr /= pio_noerr) then
      call endrun(sub//': ERROR: initial dataset not on unstructured grid')
   else
      ierr = pio_inq_dimlen(fh_ini, ncolid, ncollen)
   end if

   ! Initialize hybrid coordinate arrays
   call hycoef_init(fh_ini)

   hvcoord%hyam = hyam
   hvcoord%hyai = hyai
   hvcoord%hybm = hybm
   hvcoord%hybi = hybi
   hvcoord%ps0  = ps0
   do k=1,nlev
      hvcoord%hybd(k) = hvcoord%hybi(k+1) - hvcoord%hybi(k)
   end do
   call set_layer_locations(hvcoord, .false., par%masterproc)
   ! Initialize reference pressures
   call ref_pres_init(hypi, hypm, nprlev)

#ifdef _OPENMP
!   Set by driver
!$omp parallel
   nthreads = omp_get_num_threads()
!$omp end parallel
   if (masterproc) then
      write(iulog,*) sub//": INFO: number of OpenMP threads = ", nthreads
   end if
#if defined (COLUMN_OPENMP)
   if (masterproc) then
      write(iulog,*) sub//": INFO: using OpenMP within element instead of across elements"
   end if
#endif
#else
   nthreads = 1
   if (masterproc) then
      write(iulog,*) sub//": INFO: openmp not activated"
   end if
#endif

   if (iam < par%nprocs) then

      call prim_init1(elem,par,dom_mt,TimeLevel)

      ! globaluniquecols set by call to prim_init1
      if (ncollen /= GlobalUniqueCols) then
         write(iulog,*) sub//': ERROR: model parameters do not match initial dataset parameters'
         write(iulog,*)'  Model Parameters:    globaluniquecols = ', globaluniquecols
         write(iulog,*)'  Dataset Parameters:  ncol             = ', ncollen
         call endrun(sub//': ERROR: model parameters do not match initial dataset parameters')
      end if

      neltmp(1) = nelemdmax
      neltmp(2) = nelem
      neltmp(3) = globaluniquecols
   else
      globaluniquecols = 0
      nelemd = 0
      neltmp(1) = 0
      neltmp(2) = 0
      neltmp(3) = 0
   endif

   ! nelemdmax is computed on the dycore comm, we need it globally.
   ngcols_d = nelemdmax
   call MPI_Allreduce(ngcols_d, nelemdmax, 1, MPI_INTEGER, MPI_MAX, mpicom, ierr)
   ! All pes might not have the correct global grid size
   call MPI_Allreduce(globaluniquecols, ngcols_d, 1, MPI_INTEGER, MPI_MAX, mpicom, ierr)
   ! All pes might not have the correct number of elements
   call MPI_Allreduce(nelem, nelem_d, 1, MPI_INTEGER, MPI_MAX, mpicom, ierr)


   if (par%nprocs .lt. npes) then

      ! Broadcast quantities to auxiliary processes
      call mpi_bcast(neltmp, 3, mpi_integer, mstrid, mpicom, ierr)
      if (ierr /= mpi_success) then
         call endrun(sub//': FATAL: mpi_bcast: neltmp')
      end if

      if (iam .ge. par%nprocs) then
         nelemdmax = neltmp(1)
         nelem     = neltmp(2)
      end if
   end if

!jt run by peter
!!$   ! CAM code for Dynamics timestep
!!$   ! Dynamics timestep
!!$   !
!!$   !  Note: dtime = progress made in one timestep.  value in namelist
!!$   !        dtime = the frequency at which physics is called
!!$   !        tstep = the dynamics timestep:
!!$   dtime = get_step_size()
!!$   if (rsplit==0) then
!!$      ! non-lagrangian code
!!$      tstep = dtime/real(nsplit*qsplit,r8)
!!$      TimeLevel%nstep = get_nstep()*nsplit*qsplit
!!$   else
!!$      ! lagrangian code
!!$      tstep = dtime/real(nsplit*qsplit*rsplit,r8)
!!$      TimeLevel%nstep = get_nstep()*nsplit*qsplit*rsplit
!!$   endif

   ! Dynamics timestep    !
   !  Note: dtime = progress made in one timestep.  value in namelist
   !        dtime = the frequency at which physics is called
   !        tstep = the dynamics timestep:
   ! Ignore ierr, as on error, timestep_make_eam_parameters_consistent defaults
   ! to printing an error and then aborting.

   ierr = timestep_make_eam_parameters_consistent(par, dt_remap_factor, dt_tracer_factor, &
        nsplit, nstep_factor, tstep, dtime)
   tstep = dtime/real(nstep_factor,r8)
   TimeLevel%nstep = get_nstep()*nstep_factor

   ! initial SE (subcycled) nstep
   TimeLevel%nstep0 = 0

   ! Initialize FV physics grid variables
   if (fv_nphys > 0) then
!jt      call fv_physgrid_init(check=2)
      call fv_physgrid_init(check=2)
   end if

#ifdef HAVE_MOAB
   call create_moab_meshes(par, elem)
#endif

   ! determine whether initial file uses 'ncol' or 'ncol_d'
   call get_hdim_name(fh_ini, ini_grid_hdim_name)

   ! Define the CAM grids.
   ! Physics-grid will be defined later by phys_grid_init
   call define_cam_grids()

end subroutine dyn_grid_init

!==============================================================================

subroutine get_dyn_grid_info(hdim1_d, hdim2_d, num_lev,                       &
     index_model_top_layer, index_surface_layer, unstructured, dyn_columns)
   !------------------------------------------------------------
   !
   ! get_dyn_grid_info returns physics grid column information
   !
   !------------------------------------------------------------
   use shr_const_mod,          only: SHR_CONST_PI
   use cam_abortutils,         only: endrun
   use spmd_utils,             only: iam
   use gllfvremap_mod,         only: gfr_f_get_latlon, gfr_f_get_area

   ! Dummy arguments
   integer,          intent(out)   :: hdim1_d ! # longitudes or grid size
   integer,          intent(out)   :: hdim2_d ! # latitudes or 1
   integer,          intent(out)   :: num_lev ! # levels
   integer,          intent(out)   :: index_model_top_layer
   integer,          intent(out)   :: index_surface_layer
   logical,          intent(out)   :: unstructured
   ! dyn_columns will contain a copy of the physics column info local to this
   ! dynamics task
   type(physics_column_t), allocatable, intent(out) :: dyn_columns(:)
   ! Local variables
   integer                         :: lindex
   integer                         :: gindex
   integer                         :: elem_ind, col_ind, ii, jj
   integer                         :: num_local_cols
   real(r8)                        :: dcoord
   real(r8),         parameter     :: radtodeg = 180.0_r8 / SHR_CONST_PI
   real(r8),         parameter     :: degtorad = SHR_CONST_PI / 180.0_r8
   character(len=*), parameter     :: subname = 'get_dyn_grid_info'

   unstructured = .true. ! SE is an unstructured dycore
   if (fv_nphys > 0) then ! physics uses an FVM grid
      num_local_cols = nelemd * fv_nphys * fv_nphys
   else
      num_local_cols = 0
      do elem_ind = 1, nelemd
         num_local_cols = num_local_cols + elem(elem_ind)%idxP%NumUniquePts
      end do
   end if
   if (allocated(local_dyn_columns)) then
      ! Check for correct number of columns
      if (size(local_dyn_columns) /= num_local_cols) then
         call endrun(subname//': called with inconsistent column numbers')
      end if
   else
      allocate(local_dyn_columns(num_local_cols))
      if (fv_nphys > 0) then ! physics uses an FVM grid
         hdim1_d = nelem * fv_nphys * fv_nphys
      else
         hdim1_d = ngcols_d
      end if
      hdim2_d = 1
      num_lev = nlev
      index_model_top_layer = 1
      index_surface_layer = nlev
      lindex = 0
      do elem_ind = 1, nelemd
         if (fv_nphys > 0) then ! physics uses an FVM grid
            do col_ind = 0, (fv_nphys * fv_nphys) - 1
               lindex = lindex + 1
               ii = MOD(col_ind, fv_nphys) + 1
               jj = (col_ind / fv_nphys) + 1
               call gfr_f_get_latlon(elem_ind, ii, jj, &
                    local_dyn_columns(lindex)%lat_rad, &
                    local_dyn_columns(lindex)%lon_rad)
               dcoord = local_dyn_columns(lindex)%lat_rad * radtodeg
               local_dyn_columns(lindex)%lat_deg = dcoord
               dcoord = local_dyn_columns(lindex)%lon_rad * radtodeg
               local_dyn_columns(lindex)%lon_deg = dcoord
               local_dyn_columns(lindex)%area = gfr_f_get_area(elem_ind, ii, jj)
               local_dyn_columns(lindex)%weight =                             &
                    local_dyn_columns(lindex)%area
               ! File decomposition
               gindex = ((elem(elem_ind)%GlobalId-1) * fv_nphys * fv_nphys) + &
                    col_ind + 1
               local_dyn_columns(lindex)%global_col_num = gindex
               ! Note, coord_indices not used for unstructured dycores
               ! Dynamics decomposition
               local_dyn_columns(lindex)%dyn_task = iam
               local_dyn_columns(lindex)%local_dyn_block = elem_ind
               local_dyn_columns(lindex)%global_dyn_block =                   &
                    elem(elem_ind)%GlobalId
               allocate(local_dyn_columns(lindex)%dyn_block_index(1))
               local_dyn_columns(lindex)%dyn_block_index(1) = col_ind + 1
            end do
         else
            do col_ind = 1, elem(elem_ind)%idxP%NumUniquePts
               lindex = lindex + 1
               ii = elem(elem_ind)%idxP%ia(col_ind)
               jj = elem(elem_ind)%idxP%ja(col_ind)

               dcoord = elem(elem_ind)%spherep(ii,jj)%lat
               local_dyn_columns(lindex)%lat_rad = dcoord
               dcoord = local_dyn_columns(lindex)%lat_rad * radtodeg
               local_dyn_columns(lindex)%lat_deg = dcoord
               dcoord = elem(elem_ind)%spherep(ii,jj)%lon
               local_dyn_columns(lindex)%lon_rad = dcoord
               dcoord = local_dyn_columns(lindex)%lon_rad * radtodeg
               local_dyn_columns(lindex)%lon_deg = dcoord
               local_dyn_columns(lindex)%area =                               &
                    1.0_r8 / elem(elem_ind)%rspheremp(ii,jj)
               local_dyn_columns(lindex)%weight = local_dyn_columns(lindex)%area
               ! File decomposition
               gindex = elem(elem_ind)%idxP%UniquePtoffset + col_ind - 1
               local_dyn_columns(lindex)%global_col_num = gindex
               ! Note, coord_indices not used for unstructured dycores
               ! Dynamics decomposition
               local_dyn_columns(lindex)%dyn_task = iam
               local_dyn_columns(lindex)%local_dyn_block = elem_ind
               local_dyn_columns(lindex)%global_dyn_block =                   &
                    elem(elem_ind)%GlobalId
               allocate(local_dyn_columns(lindex)%dyn_block_index(1))
               local_dyn_columns(lindex)%dyn_block_index(1) = col_ind
            end do
         end if
      end do
   end if
   ! Copy the information to the output array
   if (allocated(dyn_columns)) then
      deallocate(dyn_columns)
   end if
   allocate(dyn_columns(lindex))
   do lindex = 1, num_local_cols
      dyn_columns(lindex) = local_dyn_columns(lindex)
   end do

   end subroutine get_dyn_grid_info

  !=================================================================================================
  !
  subroutine get_horiz_grid_dim_d(hdim1_d,hdim2_d)
    !---------------------------------------------------------------------------
    ! Purpose: Returns declared horizontal dimensions of computational grid.
    !          For non-lon/lat grids, declare grid to be one-dimensional,
    !          i.e., (ngcols_d x 1)
    !
    ! Author: Patrick Worley
    !------------------------------Arguments------------------------------------
    integer, intent(out)           :: hdim1_d  ! first horizontal dimension
    integer, intent(out), optional :: hdim2_d  ! second horizontal dimension
    !---------------------------------------------------------------------------
    if (fv_nphys > 0) then
      hdim1_d = fv_nphys*fv_nphys*nelem_d
    else
      hdim1_d = ngcols_d
    end if ! fv_nphys > 0

    if (present(hdim2_d)) hdim2_d = 1

    return
  end subroutine get_horiz_grid_dim_d

!=========================================================================================

  subroutine physgrid_copy_attributes_d(gridname, grid_attribute_names)
    use cam_grid_support, only: max_hcoordname_len
    !------------------------------Arguments------------------------------------
    character(len=max_hcoordname_len),          intent(out) :: gridname
    character(len=max_hcoordname_len), pointer, intent(out) :: grid_attribute_names(:)
    !---------------------------------------------------------------------------
    if (fv_nphys > 0) then
      gridname = 'physgrid_d'
      allocate(grid_attribute_names(2))
      grid_attribute_names(1) = 'fv_nphys'
      grid_attribute_names(2) = 'ne'
    else
      gridname = 'GLL'
      allocate(grid_attribute_names(2))
      grid_attribute_names(1) = 'np'
      grid_attribute_names(2) = 'ne'
    end if ! fv_nphys > 0

  end subroutine physgrid_copy_attributes_d

  !=========================================================================================

integer function get_dyn_grid_parm(name) result(ival)

   ! This function is in the process of being deprecated, but is still needed
   ! as a dummy interface to satisfy external references from some chemistry routines.

   use pmgrid,          only: plat, plev

   character(len=*), intent(in) :: name
   !----------------------------------------------------------------------------

   if (name.eq.'plat') then
      ival = plat
   else if(name.eq.'plon') then
      if (fv_nphys>0) then
         ival = fv_nphys*fv_nphys*nelem_d
      else
         ival = ngcols_d
      end if
   else if(name.eq.'plev') then
      ival = plev

   else
      ival = -1
   end if

end function get_dyn_grid_parm

  !=========================================================================================

  subroutine dyn_grid_get_colndx(igcol, ncols, owners, col, lbk)

    ! For each global column index return the owning task.  If the column is owned
    ! by this task, then also return the local block number and column index in that
    ! block.
    !
    ! NOTE: this routine needs to be updated for the physgrid

    integer, intent(in)  :: ncols
    integer, intent(in)  :: igcol(ncols)
    integer, intent(out) :: owners(ncols)
    integer, intent(out) :: col(ncols)
    integer, intent(out) :: lbk(ncols)

    !----------------------------------------------------------------------------

    owners = (igcol * 0) -1 ! Kill compiler warnings
    col    = -1             ! Kill compiler warnings
    lbk    = -1             ! Kill compiler warnings
    call endrun('dyn_grid_get_colndx: not implemented for unstructured grids')

  end subroutine dyn_grid_get_colndx

!=========================================================================================

  subroutine dyn_grid_get_elem_coords( ie, rlon, rlat, cdex )
    !---------------------------------------------------------------------------
    ! Purpose: this returns coordinates of a specified block element of the dyn grid
    !---------------------------------------------------------------------------
    use dof_mod, only: UniqueCoords
    !------------------------------Arguments------------------------------------
    integer,           intent(in ) :: ie      ! block element index
    real(r8),optional, intent(out) :: rlon(:) ! longitudes of the columns in the element
    real(r8),optional, intent(out) :: rlat(:) ! latitudes of the columns in the element
    integer, optional, intent(out) :: cdex(:) ! global column index
    !----------------------------Local-Variables--------------------------------
    integer :: sb,eb, ii, i,j, icol, igcol
    real(r8), allocatable :: clat(:)
    real(r8), allocatable :: clon(:)
    !---------------------------------------------------------------------------
    if (fv_nphys > 0) then
      call endrun('dyn_grid_get_elem_coords: not implemented for the FV physics grid')
    end if

    sb = elem(ie)%idxp%UniquePtOffset
    eb = sb + elem(ie)%idxp%NumUniquePts-1

    allocate( clat(sb:eb), clon(sb:eb) )
    call UniqueCoords( elem(ie)%idxP, elem(ie)%spherep, clat(sb:eb), clon(sb:eb) )

    if (present(cdex)) cdex(:) = -1
    if (present(rlat)) rlat(:) = -999._r8
    if (present(rlon)) rlon(:) = -999._r8

    do ii=1,elem(ie)%idxp%NumUniquePts
      i=elem(ie)%idxp%ia(ii)
      j=elem(ie)%idxp%ja(ii)
      icol = i+(j-1)*np
      igcol = elem(ie)%idxp%UniquePtoffset+ii-1
      if (present(cdex)) cdex(icol) = igcol
      if (present(rlat)) rlat(icol) = clat( igcol )
      if (present(rlon)) rlon(icol) = clon( igcol )
    end do

    deallocate( clat, clon )

  end subroutine dyn_grid_get_elem_coords

  !=========================================================================================
  ! Private routines.
  !=========================================================================================

 subroutine get_hdim_name(fh_ptr, grid_hdim_name)
   use pio, only: pio_inq_dimid, pio_seterrorhandling
   use pio, only: PIO_BCAST_ERROR, PIO_NOERR

   ! Determine whether the supplied file uses 'ncol' or 'ncol_d' horizontal
   ! dimension in the unstructured grid.  It is also possible when using
   ! analytic initial conditions that the file only contains
   ! vertical coordinates.
   ! Return 'none' if that is the case.

   ! Arguments
   type(file_desc_t),   pointer  :: fh_ptr
   character(len=6), intent(out) :: grid_hdim_name ! horizontal dimension name

   ! local variables
   integer  :: ierr, pio_errtype
   integer  :: ncol_did

   character(len=*), parameter :: sub = 'get_hdim_name'
   !----------------------------------------------------------------------------

   ! Set PIO to return error flags.
   call pio_seterrorhandling(fh_ptr, PIO_BCAST_ERROR, pio_errtype)

   ! Check for ncol_d first just in case the file also contains fields on
   ! the physics grid.
   ierr = pio_inq_dimid(fh_ptr, 'ncol_d', ncol_did)
   if (ierr == PIO_NOERR) then

      grid_hdim_name = 'ncol_d'

   else

      ! if 'ncol_d' not in file, check for 'ncol'
      ierr = pio_inq_dimid(fh_ptr, 'ncol', ncol_did)

      if (ierr == PIO_NOERR) then

         grid_hdim_name = 'ncol'

      else

         grid_hdim_name = 'none'

      end if
   end if

   ! Return PIO to previous error handling.
   call pio_seterrorhandling(fh_ptr, pio_errtype)

 end subroutine get_hdim_name

  !
  !=================================================================================================
  !
  subroutine define_cam_grids()
    use cam_grid_support, only: horiz_coord_t, horiz_coord_create
    use cam_grid_support, only: cam_grid_register, cam_grid_attribute_register
    use gllfvremap_mod,   only: gfr_f_get_area, gfr_f_get_latlon
    !----------------------------Local-Variables--------------------------------
    character(len=16)            :: gridname, latname, lonname, ncolname, areaname
    type(horiz_coord_t), pointer :: lat_coord
    type(horiz_coord_t), pointer :: lon_coord
    integer(iMap),       pointer :: grid_map(:,:)
    integer                      :: ie, i, j, k, mapind ! Loop variables
    real(r8)                     :: area_scm(1), lat, lon
    integer                      :: ncols_physgrid      ! local column count
    integer                      :: ngcols_physgrid     ! global column count
    integer(iMap),       pointer :: physgrid_map(:)
    real(r8),            pointer :: physgrid_area(:)
    real(r8),            pointer :: physgrid_areawt(:)
    real(r8),        allocatable :: physgrid_lat(:)
    real(r8),        allocatable :: physgrid_lon(:)
    real(r8),        allocatable :: pelat_deg(:)  ! pe-local latitudes (degrees)
    real(r8),        allocatable :: pelon_deg(:)  ! pe-local longitudes (degrees)
    real(r8),        pointer     :: pearea(:) => null()  ! pe-local areas
    real(r8),        pointer     :: pearea_wt(:) => null()  ! pe-local areas
    integer                      :: ii,jj
    integer(iMap)                :: fdofP_local(npsq,nelemd) ! pe-local map for dynamics decomp
    integer(iMap),   allocatable :: pemap(:)                 ! pe-local map for PIO decomp
    real(r8), parameter          :: rarea_unit_sphere = 1.0_r8 / (4.0_r8*PI)

    nullify(pearea_wt)
    nullify(pearea)
    nullify(physgrid_map)
    nullify(physgrid_area)
    nullify(physgrid_areawt)
    nullify(grid_map)

    !-----------------------
    ! Create GLL grid object
    !-----------------------

    ! Calculate the mapping between element GLL points and file order
    fdofp_local = 0_iMap
    do ie = 1, nelemd
       do ii = 1, elem(ie)%idxP%NumUniquePts
          i = elem(ie)%idxP%ia(ii)
          j = elem(ie)%idxP%ja(ii)
          fdofp_local((np*(j-1))+i,ie) = elem(ie)%idxP%UniquePtoffset + ii - 1
       end do
    end do


    allocate(pelat_deg(np*np*nelemd))
    allocate(pelon_deg(np*np*nelemd))
    allocate(pearea(np*np*nelemd))
    allocate(pearea_wt(np*np*nelemd))
    allocate(pemap(np*np*nelemd))

    pemap = 0_iMap
    ii = 1
    do ie = 1, nelemd
       pemap(ii:ii+npsq-1) = fdofp_local(:,ie)
       do j = 1, np
          do i = 1, np
             pearea(ii) = elem(ie)%mp(i,j)*elem(ie)%metdet(i,j)
             pearea_wt(ii) = pearea(ii)*rarea_unit_sphere
             pelat_deg(ii) = elem(ie)%spherep(i,j)%lat * rad2deg
             pelon_deg(ii) = elem(ie)%spherep(i,j)%lon * rad2deg
             ii = ii + 1
          end do
       end do
    end do

   ! If using the physics grid then the GLL grid will use the names with
   ! '_d' suffixes and the physics grid will use the unadorned names.
   ! This allows fields on both the GLL and physics grids to be written to history
   ! output files.
   if (trim(ini_grid_hdim_name) == 'ncol_d') then
      latname  = 'lat_d'
      lonname  = 'lon_d'
      ncolname = 'ncol_d'
      areaname = 'area_d'
   else
      latname  = 'lat'
      lonname  = 'lon'
      ncolname = 'ncol'
      areaname = 'area'
   end if
   lat_coord => horiz_coord_create('lat_d', 'ncol_d', ngcols_d,  &
         'latitude', 'degrees_north', 1, size(pelat_deg), pelat_deg, map=pemap)
   lon_coord => horiz_coord_create('lon_d', 'ncol_d', ngcols_d,  &
         'longitude', 'degrees_east', 1, size(pelon_deg), pelon_deg, map=pemap)

   ! Map for GLL grid
   allocate(grid_map(3,npsq*nelemd))
   grid_map = 0_iMap
   mapind = 1
   do j = 1, nelemd
      do i = 1, npsq
         grid_map(1, mapind) = i
         grid_map(2, mapind) = j
         grid_map(3, mapind) = pemap(mapind)
         mapind = mapind + 1
      end do
   end do

   ! The native SE GLL grid
   call cam_grid_register('GLL', dyn_decomp, lat_coord, lon_coord,           &
         grid_map, block_indexed=.false., unstruct=.true.)
   call cam_grid_attribute_register('GLL', 'area_d', 'gll grid areas', &
         'ncol_d', pearea, map=pemap)
   call cam_grid_attribute_register('GLL', 'area_weight_gll', 'gll grid area weights', &
         'ncol_d', pearea_wt, map=pemap)
   call cam_grid_attribute_register('GLL', 'np', '', np)
   call cam_grid_attribute_register('GLL', 'ne', '', ne)

   ! If dim name is 'ncol', create INI grid
   ! We will read from INI grid, but use GLL grid for all output
   if (trim(ini_grid_hdim_name) == 'ncol') then

      lat_coord => horiz_coord_create('lat', 'ncol', ngcols_d,  &
         'latitude', 'degrees_north', 1, size(pelat_deg), pelat_deg, map=pemap)
      lon_coord => horiz_coord_create('lon', 'ncol', ngcols_d,  &
         'longitude', 'degrees_east', 1, size(pelon_deg), pelon_deg, map=pemap)

      call cam_grid_register('INI', ini_decomp, lat_coord, lon_coord,         &
         grid_map, block_indexed=.false., unstruct=.true.)
      call cam_grid_attribute_register('INI', 'area', 'ini grid areas', &
               'ncol', pearea, map=pemap)
      call cam_grid_attribute_register('INI', 'area_weight_ini', 'ini grid area weights', &
           'ncol', pearea_wt, map=pemap)

      ini_grid_name = 'INI'
   else
      ! The dyn_decomp grid can be used to read the initial file.
      ini_grid_name = 'GLL'
   end if

   ! Coordinate values and maps are copied into the coordinate and attribute objects.
   ! Locally allocated storage is no longer needed.
   deallocate(pelat_deg)
   deallocate(pelon_deg)
   deallocate(pemap)

   ! pearea cannot be deallocated as the attribute object is just pointing
   ! to that memory.  It can be nullified since the attribute object has
   ! the reference.
   nullify(pearea)
   nullify(pearea_wt)

   ! grid_map cannot be deallocated as the cam_filemap_t object just points
   ! to it.  It can be nullified.
   nullify(grid_map)

    !---------------------------------------------------------------------------
    ! Create grid object for physics grid on the dynamics decomposition
    !---------------------------------------------------------------------------
    ! Following CAM6-SE the 'physgrid_d' grid is created and used to load the
    ! PHIS field (i.e. surface topo) on the FV physics grid and then
    ! interpolated to the GLL grid. This ensures consistent treatment of SGH.
    if (fv_nphys>0) then

      gridname = 'physgrid_d'
      latname  = 'lat'
      lonname  = 'lon'
      ncolname = 'ncol'
      areaname = 'area'

      ncols_physgrid = fv_nphys * fv_nphys * nelemd
      ngcols_physgrid = fv_nphys * fv_nphys * nelem_d
      allocate(physgrid_map(ncols_physgrid))
      allocate(physgrid_lat(ncols_physgrid))
      allocate(physgrid_lon(ncols_physgrid))
      allocate(physgrid_area(ncols_physgrid))
      allocate(physgrid_areawt(ncols_physgrid))

      ! copy local grid properties
      do ie = 1,nelemd
        k = 1
        do j = 1,fv_nphys
          do i = 1,fv_nphys
            mapind = k + (ie-1)*fv_nphys*fv_nphys
            physgrid_map(mapind) = k + (elem(ie)%GlobalId-1)*fv_nphys*fv_nphys
            physgrid_area(mapind)= gfr_f_get_area(ie, i, j)
            physgrid_areawt(mapind) = physgrid_area(mapind)*rarea_unit_sphere
            call gfr_f_get_latlon(ie, i, j, lat, lon)
            lat = lat*rad2deg
            lon = lon*rad2deg
            physgrid_lat(mapind) = lat
            physgrid_lon(mapind) = lon
            k = k + 1
          end do ! i
        end do ! j
      end do ! ie

      lat_coord => horiz_coord_create('lat', 'ncol', ngcols_physgrid,          &
                                      'latitude', 'degrees_north', 1,          &
                                      size(physgrid_lat),physgrid_lat,map=physgrid_map)

      lon_coord => horiz_coord_create('lon', 'ncol', ngcols_physgrid,          &
                                      'longitude', 'degrees_east', 1,          &
                                      size(physgrid_lon),physgrid_lon,map=physgrid_map)

      ! Map for physics grid
      allocate(grid_map(3, ncols_physgrid))
      grid_map = 0_iMap
      mapind = 1
      do j = 1, nelemd
        do i = 1, fv_nphys*fv_nphys
          grid_map(1,mapind) = i
          grid_map(2,mapind) = j
          grid_map(3,mapind) = physgrid_map(mapind)
          mapind = mapind + 1
        end do ! i
      end do ! j

      ! create physics grid object
      call cam_grid_register('physgrid_d', physgrid_d, lat_coord, lon_coord, &
                             grid_map, block_indexed=.false., unstruct=.true.)
      call cam_grid_attribute_register('physgrid_d', 'area_physgrid',         &
                                       'physics grid areas', 'ncol',   &
                                       physgrid_area, map=physgrid_map)
      call cam_grid_attribute_register('physgrid_d', 'area_weight_physgrid', 'physics grid area weight',   &
           'ncol', physgrid_areawt, map=physgrid_map)
      call cam_grid_attribute_register('physgrid_d', 'fv_nphys', '', fv_nphys)
      call cam_grid_attribute_register('physgrid_d', 'ne',       '', ne)

      deallocate(physgrid_lat)
      deallocate(physgrid_lon)
      deallocate(physgrid_map)
      nullify(physgrid_area)
      nullify(physgrid_areawt)
      nullify(grid_map)

    end if ! fv_nphys>0
    !---------------------------------------------------------------------------
    !---------------------------------------------------------------------------
    nullify(lat_coord) ! Belongs to grid
    nullify(lon_coord) ! Belongs to grid
  end subroutine define_cam_grids

  !========================================================================================

  subroutine write_grid_mapping(par, elem)

    use parallel_mod,  only: parallel_t
    use cam_pio_utils, only: cam_pio_createfile, pio_subsystem
    use pio,           only: pio_def_dim, var_desc_t, pio_int, pio_def_var, &
         pio_enddef, pio_closefile, pio_initdecomp, io_desc_t, &
         pio_write_darray, pio_freedecomp
    use dof_mod,       only: createmetadata

    ! arguments
    type(parallel_t), intent(in) :: par
    type(element_t),  intent(in) :: elem(:)

    ! local variables
    integer, parameter :: npm12 = (np-1)*(np-1)

    type(file_desc_t) :: nc
    type(var_desc_t)  :: vid
    type(io_desc_t)   :: iodesc
    integer :: dim1, dim2, ierr, i, j, ie, cc, base, ii, jj
    integer :: subelement_corners(npm12*nelemd,4)
    integer :: dof(npm12*nelemd*4)
    !----------------------------------------------------------------------------

    ! Create a CS grid mapping file for postprocessing tools

    ! write meta data for physics on GLL nodes
    call cam_pio_createfile(nc, 'SEMapping.nc', 0)

    ierr = pio_def_dim(nc, 'ncenters', npm12*nelem_d, dim1)
    ierr = pio_def_dim(nc, 'ncorners', 4, dim2)
    ierr = pio_def_var(nc, 'element_corners', PIO_INT, (/dim1,dim2/), vid)

    ierr = pio_enddef(nc)
    call createmetadata(par, elem, subelement_corners)

    jj=0
    do cc = 0, 3
       do ie = 1, nelemd
          base = ((elem(ie)%globalid-1)+cc*nelem_d)*npm12
          ii=0
          do j = 1, np-1
             do i = 1, np-1
                ii=ii+1
                jj=jj+1
                dof(jj) = base+ii
             end do
          end do
       end do
    end do

    call pio_initdecomp(pio_subsystem, pio_int, (/nelem_d*npm12,4/), dof, iodesc)

    call pio_write_darray(nc, vid, iodesc, &
         reshape(subelement_corners, (/nelemd*npm12*4/)), ierr)

    call pio_freedecomp(nc, iodesc)

    call pio_closefile(nc)

  end subroutine write_grid_mapping

!=============================================================================
!==
!!!!!! DUMMY INTERFACE TO TEST WEAK SCALING FIX, THIS SHOULD GO AWAY
!==
!=============================================================================
!=========================================================================================

subroutine get_horiz_grid_d(nxy,clat_d_out,clon_d_out,area_d_out, &
                              wght_d_out,lat_d_out,lon_d_out,cost_d_out)
    !---------------------------------------------------------------------------
    ! Purpose: Return latitude and longitude (in radians), column surface
    !          area (in radians squared) and surface integration weights
    !          for global column indices that will be passed to/from
    !          physics. Optionally also return estimated physics
    !          computational cost per global column for use in load
    !          balancing.
    !
    ! Author: Jim Edwards
    !------------------------------Arguments------------------------------------
    integer,  intent(in )                   :: nxy           ! array sizes
    real(r8), intent(out),         optional :: clat_d_out(:) ! column latitudes
    real(r8), intent(out),         optional :: clon_d_out(:) ! column longitudes
    real(r8), intent(out), target, optional :: area_d_out(:) ! column surface area
    real(r8), intent(out), target, optional :: wght_d_out(:) ! column integration weight
    real(r8), intent(out),         optional :: lat_d_out(:)  ! column degree latitudes
    real(r8), intent(out),         optional :: lon_d_out(:)  ! column degree longitudes
    real(r8), intent(out),        optional :: cost_d_out(:) ! column cost
    character(len=*), parameter     :: subname = 'get_horiz_grid_d'


   call endrun(subname//': NOT SUPPORTED WITH WEAK SCALING FIX')

 end subroutine get_horiz_grid_d

  function get_dyn_grid_parm_real1d(name) result(rval)
    !------------------------------Arguments------------------------------------
    character(len=*), intent(in) :: name
    real(r8),            pointer :: rval(:)
    !---------------------------------------------------------------------------
    if(name.eq.'w') then
      if (.not. associated(w)) then
        call endrun('get_dyn_grid_parm_real1d: w not defined')
      end if
      rval => w
    else if(name.eq.'clat') then
      call endrun('get_dyn_grid_parm_real1d: clat not supported, use get_horiz_grid_d')
    else if(name.eq.'latdeg') then
      ! This is never set so I'm calling bogus and stomping calls --goldy
      call endrun('get_dyn_grid_parm_real1d: latdeg not defined')
    else
      nullify(rval)
    end if
  end function get_dyn_grid_parm_real1d

  !===========================================================================

  subroutine fv_physgrid_init(check)
    use gllfvremap_mod,         only: gfr_init
    use parallel_mod,           only: par
    integer, intent(in), optional :: check
    call gfr_init(par, elem, fv_nphys, check=check)
  end subroutine fv_physgrid_init
  !
  !=================================================================================================
  !
  subroutine fv_physgrid_final()
    use gllfvremap_mod,   only: gfr_finish

    call gfr_finish()
  end subroutine fv_physgrid_final
  !
  !=================================================================================================
  !
end module dyn_grid
