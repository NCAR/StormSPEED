module restart_dynamics
   !----------------------------------------------------------------------------
   ! Unified SE Restart Interface
   ! Dispatches via Type-Bound Procedures on dyn_export_t / dyn_import_t
   !----------------------------------------------------------------------------
   use pio,      only: file_desc_t
   use dyn_comp, only: dyn_import_t, dyn_export_t

   implicit none
   private
   save

   public :: init_restart_dynamics, write_restart_dynamics, read_restart_dynamics

contains

   subroutine init_restart_dynamics(File, dyn_out)
      type(file_desc_t),  intent(inout) :: File
      type(dyn_export_t), intent(in)    :: dyn_out
      call dyn_out%init_restart(File)
   end subroutine init_restart_dynamics

   subroutine write_restart_dynamics(File, dyn_out)
      type(file_desc_t),  intent(inout) :: File
      type(dyn_export_t), intent(in)    :: dyn_out
      call dyn_out%write_restart(File)
   end subroutine write_restart_dynamics

   subroutine read_restart_dynamics(File, dyn_in, dyn_out)
      type(file_desc_t),  intent(inout) :: File
      type(dyn_import_t), intent(out)   :: dyn_in
      type(dyn_export_t), intent(out)   :: dyn_out
      call dyn_in%read_restart(File, dyn_out)
   end subroutine read_restart_dynamics

end module restart_dynamics
