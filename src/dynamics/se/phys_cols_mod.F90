module phys_cols_mod
   !----------------------------------------------------------------------------
   ! Physics-column arrays owned by the SE dycore backend (se_backend_t%phys).
   !
   ! The backend fills these from the dynamics state on the physics grid; each
   ! host model (CAM, CAM-SIMA) copies them into its own physics state.  The
   ! arrays use only intrinsic types so the module carries no host dependency.
   !
   ! Layout (CCPP order): column index first, then level, then constituent.
   ! Columns follow the host's physics-column order on the task
   ! (col_ind = 1..ncol, as in phys_grid's get_dyn_col_p / get_chunk_info_p).
   !----------------------------------------------------------------------------
   use shr_kind_mod, only: r8 => shr_kind_r8

   implicit none
   private

   type, public :: phys_cols_t
      integer :: ncol  = 0           ! physics columns on this task
      integer :: nlev  = 0           ! layers
      integer :: ncnst = 0           ! constituents
      logical :: has_front = .false. ! frontgf/frontga allocated
      logical :: has_vort  = .false. ! vort4gw allocated

      ! Dynamics -> physics
      real(r8), allocatable :: ps(:)          ! surface pressure (Pa)
      real(r8), allocatable :: phis(:)        ! surface geopotential (m2/s2)
      real(r8), allocatable :: t(:,:)         ! temperature (K)
      real(r8), allocatable :: u(:,:)         ! zonal wind (m/s)
      real(r8), allocatable :: v(:,:)         ! meridional wind (m/s)
      real(r8), allocatable :: omega(:,:)     ! vertical pressure velocity
      real(r8), allocatable :: pdel(:,:)      ! layer thickness, when the dycore provides it (Pa)
      real(r8), allocatable :: q(:,:,:)       ! constituent mixing ratios
      real(r8), allocatable :: frontgf(:,:)   ! frontogenesis function
      real(r8), allocatable :: frontga(:,:)   ! frontogenesis angle
      real(r8), allocatable :: vort4gw(:,:)   ! relative vorticity for gravity waves
   contains
      procedure :: alloc   => phys_cols_alloc
      procedure :: dealloc => phys_cols_dealloc
   end type phys_cols_t

contains

   subroutine phys_cols_alloc(this, ncol, nlev, ncnst, with_front, with_vort)
      ! Allocate (or keep) the arrays for this task's columns.  Safe to call
      ! every step: it only reallocates if a dimension or option changed.
      class(phys_cols_t), intent(inout) :: this
      integer,            intent(in)    :: ncol, nlev, ncnst
      logical,            intent(in)    :: with_front, with_vort

      if (allocated(this%t)) then
         if (this%ncol == ncol .and. this%nlev == nlev .and. this%ncnst == ncnst .and. &
             (this%has_front .eqv. with_front) .and. (this%has_vort .eqv. with_vort)) return
         call this%dealloc()
      end if

      this%ncol      = ncol
      this%nlev      = nlev
      this%ncnst     = ncnst
      this%has_front = with_front
      this%has_vort  = with_vort

      allocate(this%ps(ncol), this%phis(ncol))
      allocate(this%t(ncol,nlev), this%u(ncol,nlev), this%v(ncol,nlev))
      allocate(this%omega(ncol,nlev), this%pdel(ncol,nlev))
      allocate(this%q(ncol,nlev,ncnst))
      if (with_front) allocate(this%frontgf(ncol,nlev), this%frontga(ncol,nlev))
      if (with_vort)  allocate(this%vort4gw(ncol,nlev))
   end subroutine phys_cols_alloc

   subroutine phys_cols_dealloc(this)
      class(phys_cols_t), intent(inout) :: this
      if (allocated(this%ps))      deallocate(this%ps)
      if (allocated(this%phis))    deallocate(this%phis)
      if (allocated(this%t))       deallocate(this%t)
      if (allocated(this%u))       deallocate(this%u)
      if (allocated(this%v))       deallocate(this%v)
      if (allocated(this%omega))   deallocate(this%omega)
      if (allocated(this%pdel))    deallocate(this%pdel)
      if (allocated(this%q))       deallocate(this%q)
      if (allocated(this%frontgf)) deallocate(this%frontgf)
      if (allocated(this%frontga)) deallocate(this%frontga)
      if (allocated(this%vort4gw)) deallocate(this%vort4gw)
      this%ncol = 0; this%nlev = 0; this%ncnst = 0
      this%has_front = .false.; this%has_vort = .false.
   end subroutine phys_cols_dealloc

end module phys_cols_mod
