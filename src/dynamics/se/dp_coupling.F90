module dp_coupling
   !----------------------------------------------------------------------------
   ! Unified Physics-Dynamics Coupling Interface
   ! Dispatches via Type-Bound Procedures on dyn_export_t / dyn_import_t
   !----------------------------------------------------------------------------
   use ppgrid,         only: begchunk, endchunk
   use physics_types,  only: physics_state, physics_tend
   use physics_buffer, only: physics_buffer_desc
   use dyn_comp,       only: dyn_export_t, dyn_import_t

   implicit none
   private
   save

   public :: d_p_coupling, p_d_coupling

contains

   subroutine d_p_coupling(phys_state, phys_tend, pbuf2d, dyn_out)
      type(physics_state), intent(inout), dimension(begchunk:endchunk) :: phys_state
      type(physics_tend ), intent(inout), dimension(begchunk:endchunk) :: phys_tend
      type(physics_buffer_desc), pointer                               :: pbuf2d(:,:)
      type(dyn_export_t),  intent(inout)                               :: dyn_out

      call dyn_out%d_p_coupling(phys_state, phys_tend, pbuf2d)
   end subroutine d_p_coupling

   subroutine p_d_coupling(phys_state, phys_tend, dyn_in)
      type(physics_state), intent(inout), dimension(begchunk:endchunk) :: phys_state
      type(physics_tend),  intent(inout), dimension(begchunk:endchunk) :: phys_tend
      type(dyn_import_t),  intent(inout)                               :: dyn_in

      call dyn_in%p_d_coupling(phys_state, phys_tend)
   end subroutine p_d_coupling

end module dp_coupling
