module dyn_comp
   !----------------------------------------------------------------------------
   ! Unified CAM-SE Dynamical Core Component Interface
   ! Supports: 'cam', 'preqx', 'theta-l', and 'theta-l_kokkos' via se_backend_t
   !----------------------------------------------------------------------------
   use se_backend_mod, only: se_backend_t, elem_state_t
   use dyn_comp_impl,  only: frontgf_idx, frontga_idx, vort4gw_idx, write_restart_unstruct

   implicit none
   private
   save

   public :: dyn_import_t, dyn_export_t
   public :: dyn_readnl, dyn_register, dyn_init, dyn_run, dyn_final
   public :: frontgf_idx, frontga_idx, vort4gw_idx, write_restart_unstruct

   type, extends(se_backend_t) :: dyn_import_t
   end type dyn_import_t

   type, extends(se_backend_t) :: dyn_export_t
   end type dyn_export_t

contains

   subroutine dyn_readnl(NLFileName)
      character(len=*), intent(in) :: NLFileName
      type(dyn_import_t) :: backend
      call backend%readnl(NLFileName)
   end subroutine dyn_readnl

   subroutine dyn_register()
      type(dyn_import_t) :: backend
      call backend%register_pbuf()
   end subroutine dyn_register

   subroutine dyn_init(dyn_in, dyn_out)
      type(dyn_import_t), intent(out) :: dyn_in
      type(dyn_export_t), intent(out) :: dyn_out
      call dyn_in%init_model(dyn_out)
   end subroutine dyn_init

   subroutine dyn_run(dyn_state)
      type(dyn_export_t), intent(inout) :: dyn_state
      call dyn_state%run_subcycle()
   end subroutine dyn_run

   subroutine dyn_final(DYN_STATE, RESTART_FILE)
      type(elem_state_t), target     :: DYN_STATE
      character(LEN=*),   intent(in) :: RESTART_FILE
      type(dyn_export_t) :: backend
      call backend%finalize(DYN_STATE, RESTART_FILE)
   end subroutine dyn_final

end module dyn_comp
