module se_backend_mod
   !----------------------------------------------------------------------------
   ! VPATH-selected F2003 Object-Oriented Strategy Backend for 'cam' target
   ! (Dry-pressure vertical coordinate, CSLAM FVM tracer transport & physgrid)
   !----------------------------------------------------------------------------
   use shr_kind_mod,           only: r8 => shr_kind_r8
   use ppgrid,                 only: begchunk, endchunk
   use physics_types,          only: physics_state, physics_tend
   use physics_buffer,         only: physics_buffer_desc
   use camsrfexch,             only: cam_out_t
   use pio,                    only: file_desc_t
   use physics_column_type,    only: physics_column_t
   use cam_grid_support,       only: max_hcoordname_len
   use cam_history_support,    only: interp_info_t

   use dimensions_mod,         only: np, npsq, ne, nelem, nelemd, nelemdmax, &
                                     nlev, fv_nphys, qsize, ntrac, use_cslam, &
                                     cnst_name_gll, cnst_longname_gll
   use control_mod,            only: ftype, qsplit, rsplit, statefreq, runtype
   use parallel_mod,           only: par, parallel_t
   use thread_mod,             only: max_num_threads, horz_num_threads, &
                                     vert_num_threads, tracer_num_threads
   use se_dyn_time_mod,        only: TimeLevel_t, TimeLevel_Qdp, tstep, nsplit
   use element_mod,            only: element_t, elem_state_t
   use fvm_control_volume_mod, only: fvm_struct
   use hybvcoord_mod,          only: hvcoord_t

   implicit none
   private

   ! Describes the shared element data (dyn_in/dyn_out point at the same elem),
   ! not one backend instance; kept at module level so read-only (intent(in))
   ! paths such as write_restart may still sync and clear it.
   logical, save :: host_state_dirty = .false.

   public :: np, npsq, ne, nelem, nelemd, nelemdmax, nlev, fv_nphys, qsize, ntrac, use_cslam
   public :: cnst_name_gll, cnst_longname_gll
   public :: ftype, qsplit, rsplit, statefreq, runtype
   public :: par, parallel_t
   public :: max_num_threads, horz_num_threads, vert_num_threads, tracer_num_threads
   public :: TimeLevel_t, TimeLevel_Qdp, tstep, nsplit
   public :: element_t, elem_state_t, fvm_struct, hvcoord_t

   type, public :: se_backend_t
      type(element_t),  pointer :: elem(:) => null()
      type(fvm_struct), pointer :: fvm(:)  => null()
      logical                   :: gpu_resident     = .false.
   contains
      procedure, nopass :: readnl         => backend_readnl
      procedure, nopass :: register_pbuf  => backend_register_pbuf
      procedure         :: init_model     => backend_init_model
      procedure         :: run_subcycle   => backend_run_subcycle
      procedure, nopass :: finalize       => backend_finalize

      procedure         :: sync_to_device => backend_sync_to_device
      procedure         :: sync_to_host   => backend_sync_to_host

      procedure         :: d_p_coupling   => backend_d_p_coupling
      procedure         :: p_d_coupling   => backend_p_d_coupling

      procedure         :: stepon_init    => backend_stepon_init
      procedure         :: stepon_run1    => backend_stepon_run1
      procedure         :: stepon_run2    => backend_stepon_run2
      procedure         :: stepon_run3    => backend_stepon_run3
      procedure         :: stepon_final   => backend_stepon_final

      procedure         :: init_restart   => backend_init_restart
      procedure         :: write_restart  => backend_write_restart
      procedure         :: read_restart   => backend_read_restart
   end type se_backend_t

contains

   subroutine backend_readnl(NLFileName)
      use dyn_comp_impl, only: impl_readnl => dyn_readnl
      character(len=*), intent(in) :: NLFileName
      call impl_readnl(NLFileName)
   end subroutine backend_readnl

   subroutine backend_register_pbuf()
      use dyn_comp_impl, only: impl_reg => dyn_register
      call impl_reg()
   end subroutine backend_register_pbuf

   subroutine backend_init_model(this, dyn_out)
      use dyn_comp_impl, only: impl_init => dyn_init, impl_imp_t => dyn_import_t, impl_exp_t => dyn_export_t
      class(se_backend_t), intent(inout) :: this
      class(se_backend_t), intent(inout) :: dyn_out
      type(impl_imp_t) :: imp
      type(impl_exp_t) :: exp
      call impl_init(imp, exp)
      this%elem    => imp%elem
      this%fvm     => imp%fvm
      dyn_out%elem => exp%elem
      dyn_out%fvm  => exp%fvm
   end subroutine backend_init_model

   subroutine backend_run_subcycle(this)
      use dyn_comp_impl, only: impl_run => dyn_run, impl_exp_t => dyn_export_t
      class(se_backend_t), intent(inout) :: this
      type(impl_exp_t) :: exp
      exp%elem => this%elem
      exp%fvm  => this%fvm
      call impl_run(exp)
   end subroutine backend_run_subcycle

   subroutine backend_finalize(DYN_STATE, RESTART_FILE)
      use dyn_comp_impl, only: impl_final => dyn_final
      type(elem_state_t), target     :: DYN_STATE
      character(LEN=*),   intent(in) :: RESTART_FILE
      call impl_final(DYN_STATE, RESTART_FILE)
   end subroutine backend_finalize

   subroutine backend_sync_to_device(this)
      class(se_backend_t), intent(inout) :: this
   end subroutine backend_sync_to_device

   subroutine backend_sync_to_host(this)
      class(se_backend_t), intent(in)    :: this
   end subroutine backend_sync_to_host

   subroutine backend_d_p_coupling(this, phys_state, phys_tend, pbuf2d)
      use dp_coupling_impl, only: impl_dp => d_p_coupling
      use dyn_comp_impl,    only: impl_exp_t => dyn_export_t
      class(se_backend_t),               intent(inout) :: this
      type(physics_state),               intent(inout) :: phys_state(begchunk:endchunk)
      type(physics_tend),                intent(inout) :: phys_tend(begchunk:endchunk)
      type(physics_buffer_desc), pointer               :: pbuf2d(:,:)
      type(impl_exp_t) :: exp
      exp%elem => this%elem
      exp%fvm  => this%fvm
      call impl_dp(phys_state, phys_tend, pbuf2d, exp)
   end subroutine backend_d_p_coupling

   subroutine backend_p_d_coupling(this, phys_state, phys_tend)
      use dp_coupling_impl, only: impl_pd => p_d_coupling
      use dyn_comp_impl,    only: impl_imp_t => dyn_import_t
      use dyn_grid_impl,    only: TimeLevel
      class(se_backend_t), intent(inout) :: this
      type(physics_state), intent(inout) :: phys_state(begchunk:endchunk)
      type(physics_tend),  intent(inout) :: phys_tend(begchunk:endchunk)
      type(impl_imp_t) :: imp
      integer          :: tl_f, tl_fQdp
      tl_f = TimeLevel%n0
      call TimeLevel_Qdp(TimeLevel, qsplit, tl_fQdp)
      imp%elem => this%elem
      imp%fvm  => this%fvm
      call impl_pd(phys_state, phys_tend, imp, tl_f, tl_fQdp)
   end subroutine backend_p_d_coupling

   subroutine backend_stepon_init(this, dyn_out)
      use stepon_impl,   only: impl_sinit => stepon_init
      use dyn_comp_impl, only: impl_imp_t => dyn_import_t, impl_exp_t => dyn_export_t
      class(se_backend_t), intent(inout) :: this
      class(se_backend_t), intent(inout) :: dyn_out
      type(impl_imp_t) :: imp
      type(impl_exp_t) :: exp
      imp%elem => this%elem; imp%fvm => this%fvm
      exp%elem => dyn_out%elem; exp%fvm => dyn_out%fvm
      call impl_sinit(imp, exp)
   end subroutine backend_stepon_init

   subroutine backend_stepon_run1(this, dtime_out, phys_state, phys_tend, pbuf2d, dyn_out)
      use stepon_impl,   only: impl_srun1 => stepon_run1
      use dyn_comp_impl, only: impl_imp_t => dyn_import_t, impl_exp_t => dyn_export_t
      class(se_backend_t),               intent(inout) :: this
      real(r8),                          intent(out)   :: dtime_out
      type(physics_state),               intent(inout) :: phys_state(begchunk:endchunk)
      type(physics_tend),                intent(inout) :: phys_tend(begchunk:endchunk)
      type(physics_buffer_desc), pointer               :: pbuf2d(:,:)
      class(se_backend_t),               intent(inout) :: dyn_out
      type(impl_imp_t) :: imp
      type(impl_exp_t) :: exp
      imp%elem => this%elem; imp%fvm => this%fvm
      exp%elem => dyn_out%elem; exp%fvm => dyn_out%fvm
      call impl_srun1(dtime_out, phys_state, phys_tend, pbuf2d, imp, exp)
   end subroutine backend_stepon_run1

   subroutine backend_stepon_run2(this, phys_state, phys_tend, dyn_out)
      use stepon_impl,   only: impl_srun2 => stepon_run2
      use dyn_comp_impl, only: impl_imp_t => dyn_import_t, impl_exp_t => dyn_export_t
      class(se_backend_t), intent(inout) :: this
      type(physics_state), intent(inout) :: phys_state(begchunk:endchunk)
      type(physics_tend),  intent(inout) :: phys_tend(begchunk:endchunk)
      class(se_backend_t), intent(inout) :: dyn_out
      type(impl_imp_t) :: imp
      type(impl_exp_t) :: exp
      imp%elem => this%elem; imp%fvm => this%fvm
      exp%elem => dyn_out%elem; exp%fvm => dyn_out%fvm
      call impl_srun2(phys_state, phys_tend, imp, exp)
   end subroutine backend_stepon_run2

   subroutine backend_stepon_run3(this, dtime, cam_out, phys_state, dyn_out)
      use stepon_impl,   only: impl_srun3 => stepon_run3
      use dyn_comp_impl, only: impl_imp_t => dyn_import_t, impl_exp_t => dyn_export_t
      class(se_backend_t), intent(inout) :: this
      real(r8),            intent(in)    :: dtime
      type(cam_out_t),     intent(inout) :: cam_out(:)
      type(physics_state), intent(inout) :: phys_state(begchunk:endchunk)
      class(se_backend_t), intent(inout) :: dyn_out
      type(impl_imp_t) :: imp
      type(impl_exp_t) :: exp
      imp%elem => this%elem; imp%fvm => this%fvm
      exp%elem => dyn_out%elem; exp%fvm => dyn_out%fvm
      call impl_srun3(dtime, cam_out, phys_state, imp, exp)
   end subroutine backend_stepon_run3

   subroutine backend_stepon_final(this, dyn_out)
      use stepon_impl,   only: impl_sfinal => stepon_final
      use dyn_comp_impl, only: impl_imp_t => dyn_import_t, impl_exp_t => dyn_export_t
      class(se_backend_t), intent(inout) :: this
      class(se_backend_t), intent(inout) :: dyn_out
      type(impl_imp_t) :: imp
      type(impl_exp_t) :: exp
      imp%elem => this%elem; imp%fvm => this%fvm
      exp%elem => dyn_out%elem; exp%fvm => dyn_out%fvm
      call impl_sfinal(imp, exp)
   end subroutine backend_stepon_final

   subroutine backend_init_restart(this, File)
      use restart_dynamics_impl, only: impl_irest => init_restart_dynamics
      use dyn_comp_impl,         only: impl_exp_t => dyn_export_t
      class(se_backend_t), intent(in)    :: this
      type(file_desc_t),   intent(inout) :: File
      type(impl_exp_t) :: exp
      exp%elem => this%elem; exp%fvm => this%fvm
      call impl_irest(File, exp)
   end subroutine backend_init_restart

   subroutine backend_write_restart(this, File)
      use restart_dynamics_impl, only: impl_wrest => write_restart_dynamics
      use dyn_comp_impl,         only: impl_exp_t => dyn_export_t
      class(se_backend_t), intent(in)    :: this
      type(file_desc_t),   intent(inout) :: File
      type(impl_exp_t) :: exp
      exp%elem => this%elem; exp%fvm => this%fvm
      call impl_wrest(File, exp)
   end subroutine backend_write_restart

   subroutine backend_read_restart(this, File, dyn_out)
      use restart_dynamics_impl, only: impl_rrest => read_restart_dynamics
      use dyn_comp_impl,         only: impl_imp_t => dyn_import_t, impl_exp_t => dyn_export_t
      class(se_backend_t), intent(inout) :: this
      type(file_desc_t),   intent(inout) :: File
      class(se_backend_t), intent(inout) :: dyn_out
      type(impl_imp_t) :: imp
      type(impl_exp_t) :: exp
      call impl_rrest(File, imp, exp)
      this%elem    => imp%elem; this%fvm    => imp%fvm
      dyn_out%elem => exp%elem; dyn_out%fvm => exp%fvm
   end subroutine backend_read_restart

end module se_backend_mod
