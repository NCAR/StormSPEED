module dp_coupling
   !----------------------------------------------------------------------------
   ! Unified Physics-Dynamics Coupling Interface
   ! Dispatches via Type-Bound Procedures on dyn_export_t / dyn_import_t
   !----------------------------------------------------------------------------
   use shr_kind_mod,   only: r8 => shr_kind_r8
   use ppgrid,         only: begchunk, endchunk
   use physics_types,  only: physics_state, physics_tend
   use physics_buffer, only: physics_buffer_desc
   use dyn_comp,       only: dyn_export_t, dyn_import_t
   use phys_cols_mod,  only: phys_cols_t
   use dp_coupling_cam, only: phys_state_from_backend, phys_tend_to_backend

   implicit none
   private
   save

   public :: d_p_coupling, p_d_coupling
   public :: d_p_coupling_cam_finish  ! CAM-side part of d_p_coupling (also used by stepon_run1)

contains

   subroutine d_p_coupling(phys_state, phys_tend, pbuf2d, dyn_out)
      type(physics_state), intent(inout), dimension(begchunk:endchunk) :: phys_state
      type(physics_tend ), intent(inout), dimension(begchunk:endchunk) :: phys_tend
      type(physics_buffer_desc), pointer                               :: pbuf2d(:,:)
      type(dyn_export_t),  intent(inout)                               :: dyn_out

      call dyn_out%d_p_coupling()
      call phys_state_from_backend(phys_state, phys_tend, dyn_out)
      call d_p_coupling_cam_finish(phys_state, phys_tend, pbuf2d, dyn_out%phys)
   end subroutine d_p_coupling

   subroutine p_d_coupling(phys_state, phys_tend, dyn_in)
      type(physics_state), intent(inout), dimension(begchunk:endchunk) :: phys_state
      type(physics_tend),  intent(inout), dimension(begchunk:endchunk) :: phys_tend
      type(dyn_import_t),  intent(inout)                               :: dyn_in

      call phys_tend_to_backend(phys_state, phys_tend, dyn_in)
      call dyn_in%p_d_coupling()
   end subroutine p_d_coupling

   subroutine d_p_coupling_cam_finish(phys_state, phys_tend, pbuf2d, cols)
      !-------------------------------------------------------------------------
      ! CAM's part of dynamics-to-physics coupling, done after the backend has
      ! filled the physics-column arrays (cols) and CAM's physics state:
      !   - gravity-wave fields (FRONTGF, FRONTGA, VORT4GW) into CAM's pbuf
      !   - energy and water integrals of the input state (check_energy)
      ! Keeps CAM's physics buffer out of the SE dycore backend.
      !-------------------------------------------------------------------------
      use physics_buffer, only: pbuf_get_chunk, pbuf_get_field, pbuf_get_index
      use phys_grid,      only: get_chunk_info_p, phys_columns_on_task
      use check_energy,   only: check_energy_timestep_init

      type(physics_state),       intent(inout) :: phys_state(begchunk:endchunk)
      type(physics_tend ),       intent(inout) :: phys_tend(begchunk:endchunk)
      type(physics_buffer_desc), pointer       :: pbuf2d(:,:)
      type(phys_cols_t),         intent(in)    :: cols

      type(physics_buffer_desc), pointer :: pbuf_chnk(:)
      real(r8),                  pointer :: fld1(:,:), fld2(:,:)
      integer                            :: col_ind, lchnk, icol, k
      integer, save                      :: frontgf_idx = -1, frontga_idx = -1
      integer, save                      :: vort4gw_idx = -1

      if (cols%has_front) then
         if (frontgf_idx < 0) frontgf_idx = pbuf_get_index('FRONTGF')
         if (frontga_idx < 0) frontga_idx = pbuf_get_index('FRONTGA')
         !$omp parallel do private (col_ind, lchnk, icol, k, pbuf_chnk, fld1, fld2)
         do col_ind = 1, phys_columns_on_task
            call get_chunk_info_p(col_ind, lchnk, icol)
            pbuf_chnk => pbuf_get_chunk(pbuf2d, lchnk)
            call pbuf_get_field(pbuf_chnk, frontgf_idx, fld1)
            call pbuf_get_field(pbuf_chnk, frontga_idx, fld2)
            do k = 1, cols%nlev
               fld1(icol, k) = cols%frontgf(col_ind, k)
               fld2(icol, k) = cols%frontga(col_ind, k)
            end do
         end do
      end if

      if (cols%has_vort) then
         if (vort4gw_idx < 0) vort4gw_idx = pbuf_get_index('VORT4GW')
         !$omp parallel do private (col_ind, lchnk, icol, k, pbuf_chnk, fld1)
         do col_ind = 1, phys_columns_on_task
            call get_chunk_info_p(col_ind, lchnk, icol)
            pbuf_chnk => pbuf_get_chunk(pbuf2d, lchnk)
            call pbuf_get_field(pbuf_chnk, vort4gw_idx, fld1)
            do k = 1, cols%nlev
               fld1(icol, k) = cols%vort4gw(col_ind, k)
            end do
         end do
      end if

      ! Energy and water integrals of the input state
      do lchnk = begchunk, endchunk
         pbuf_chnk => pbuf_get_chunk(pbuf2d, lchnk)
         call check_energy_timestep_init(phys_state(lchnk), phys_tend(lchnk), pbuf_chnk)
      end do
   end subroutine d_p_coupling_cam_finish

end module dp_coupling
