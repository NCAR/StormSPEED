module stepon
   !----------------------------------------------------------------------------
   ! Unified SE Dynamics Time-Stepping Interface
   !----------------------------------------------------------------------------
   use shr_kind_mod,   only: r8 => shr_kind_r8
   use ppgrid,         only: begchunk, endchunk
   use physics_types,  only: physics_state, physics_tend
   use physics_buffer, only: physics_buffer_desc
   use camsrfexch,     only: cam_out_t
   use dyn_comp,       only: dyn_import_t, dyn_export_t

   implicit none
   private
   save

   public :: stepon_init, stepon_run1, stepon_run2, stepon_run3, stepon_final

contains

   subroutine stepon_init(dyn_in, dyn_out)
      type(dyn_import_t), intent(inout) :: dyn_in
      type(dyn_export_t), intent(inout) :: dyn_out
      call dyn_in%stepon_init(dyn_out)
   end subroutine stepon_init

   subroutine stepon_run1(dtime_out, phys_state, phys_tend, pbuf2d, dyn_in, dyn_out)
      real(r8),                  intent(out)   :: dtime_out
      type(physics_state),       intent(inout) :: phys_state(begchunk:endchunk)
      type(physics_tend),        intent(inout) :: phys_tend(begchunk:endchunk)
      type(physics_buffer_desc), pointer       :: pbuf2d(:,:)
      type(dyn_import_t),        intent(inout) :: dyn_in
      type(dyn_export_t),        intent(inout) :: dyn_out
      call dyn_in%stepon_run1(dtime_out, phys_state, phys_tend, pbuf2d, dyn_out)
   end subroutine stepon_run1

   subroutine stepon_run2(phys_state, phys_tend, dyn_in, dyn_out)
      type(physics_state), intent(inout) :: phys_state(begchunk:endchunk)
      type(physics_tend),  intent(inout) :: phys_tend(begchunk:endchunk)
      type(dyn_import_t),  intent(inout) :: dyn_in
      type(dyn_export_t),  intent(inout) :: dyn_out
      call dyn_in%stepon_run2(phys_state, phys_tend, dyn_out)
   end subroutine stepon_run2

   subroutine stepon_run3(dtime, cam_out, phys_state, dyn_in, dyn_out)
      real(r8),            intent(in)    :: dtime
      type(cam_out_t),     intent(inout) :: cam_out(:)
      type(physics_state), intent(inout) :: phys_state(begchunk:endchunk)
      type(dyn_import_t),  intent(inout) :: dyn_in
      type(dyn_export_t),  intent(inout) :: dyn_out
      call dyn_in%stepon_run3(dtime, cam_out, phys_state, dyn_out)
   end subroutine stepon_run3

   subroutine stepon_final(dyn_in, dyn_out)
      type(dyn_import_t), intent(inout) :: dyn_in
      type(dyn_export_t), intent(inout) :: dyn_out
      call dyn_in%stepon_final(dyn_out)
   end subroutine stepon_final

end module stepon
