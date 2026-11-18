module dp_coupling_impl

!-------------------------------------------------------------------------------
! dynamics - physics coupling module
!-------------------------------------------------------------------------------

use shr_kind_mod,   only: r8=>shr_kind_r8
use ppgrid,         only: begchunk, endchunk, pcols, pver, pverp
use constituents,   only: pcnst, cnst_type

use spmd_dyn,       only: local_dp_map
use spmd_utils,     only: iam
use dyn_grid_impl,       only: TimeLevel, edgebuf
use dyn_comp_impl,       only: dyn_export_t, dyn_import_t

use physics_types,  only: physics_state, physics_tend, physics_cnst_limit
use phys_grid,      only: get_ncols_p
use phys_grid,      only: get_dyn_col_p, columns_on_task, get_chunk_info_p, phys_columns_on_task

use dp_mapping,     only: nphys_pts

use perf_mod,       only: t_startf, t_stopf
use cam_abortutils, only: endrun

use parallel_mod,   only: par
use thread_mod,     only: horz_num_threads, max_num_threads
use phys_cols_mod,  only: phys_cols_t
use hybrid_mod,     only: config_thread_region, get_loop_ranges, hybrid_t
use dimensions_mod, only: np, nelemd, nlev, qsize, ntrac, fv_nphys

use dof_mod,        only: UniquePoints, PutUniquePoints
use element_mod,    only: element_t

implicit none
private
save

public :: d_p_coupling, p_d_coupling
public :: debug_overwrite_tendencies   ! -Ddebug_coupling only, via se_dynamics_t%debug_phys_tend

real (kind=r8), allocatable :: q_prev(:,:,:,:) ! Previous Q for computing tendencies

!=========================================================================================
CONTAINS
!=========================================================================================

subroutine d_p_coupling(dyn_out)

   ! Convert the dynamics output state into the physics input state.
   ! Note that all pressures and tracer mixing ratios coming from the dycore are based on
   ! dry air mass.

   use gravity_waves_sources_impl,  only: gws_src_fnct,gws_src_vort
   use phys_control,           only: use_gw_front, use_gw_front_igw,  use_gw_movmtn_pbl
   use hycoef,                 only: hyai, ps0
   use fvm_mapping,            only: dyn2phys_vector, dyn2phys_all_vars
   use se_dyn_time_mod,        only: timelevel_qdp
   use control_mod,            only: qsplit
   use test_fvm_mapping_impl,       only: test_mapping_overwrite_dyn_state
   use prim_advance_mod,       only: tot_energy_dyn
   use phys_cols_mod,          only: phys_cols_t
   use dyn_tests_utils,        only: vc_dry_pressure
   ! arguments
   type(dyn_export_t),  intent(inout)                               :: dyn_out             ! dynamics export


   ! LOCAL VARIABLES
   type(element_t), pointer     :: elem(:)             ! pointer to dyn_out element array
   integer                      :: ie                  ! indices over elements
   integer                      :: lchnk, icol, ilyr   ! indices over chunks, columns, layers

   real (kind=r8),  allocatable :: ps_tmp(:,:)         ! temp array to hold ps
   real (kind=r8),  allocatable :: dp3d_tmp(:,:,:)     ! temp array to hold dp3d
   real (kind=r8),  allocatable :: dp3d_tmp_tmp(:,:)
   real (kind=r8),  allocatable :: phis_tmp(:,:)       ! temp array to hold phis
   real (kind=r8),  allocatable :: T_tmp(:,:,:)        ! temp array to hold T
   real (kind=r8),  allocatable :: uv_tmp(:,:,:,:)     ! temp array to hold u and v
   real (kind=r8),  allocatable :: q_tmp(:,:,:,:)      ! temp to hold advected constituents
   real (kind=r8),  allocatable :: omega_tmp(:,:,:)    ! temp array to hold omega

   ! Frontogenesis
   real (kind=r8),  allocatable :: frontgf(:,:,:)      ! temp arrays to hold frontogenesis
   real (kind=r8),  allocatable :: frontga(:,:,:)      ! function (frontgf) and angle (frontga)

   ! Vorticity
   real (kind=r8),  allocatable :: vort4gw(:,:,:)      ! temp arrays to hold vorticity


   type(phys_cols_t), pointer   :: cols

   integer                      :: ncols, ierr
   integer                      :: col_ind, blk_ind(1), m
   integer                      :: nphys

   real (kind=r8),  allocatable :: qgll(:,:,:,:)
   real (kind=r8)               :: inv_dp3d(np,np,nlev)
   integer                      :: tl_f, tl_qdp_np0, tl_qdp_np1

   !----------------------------------------------------------------------------

   if (.not. local_dp_map) then
      call endrun('d_p_coupling: Weak scaling does not support load balancing')
   end if

   elem => dyn_out%elem
   tl_f = TimeLevel%n0
   call TimeLevel_Qdp(TimeLevel, qsplit, tl_qdp_np0,tl_qdp_np1)




   if (fv_nphys > 0) then
      nphys = fv_nphys
   else
     allocate(qgll(np,np,nlev,pcnst))
     nphys = np
   end if

   ! Allocate temporary arrays to hold data for physics decomposition
   allocate(ps_tmp(nphys_pts,nelemd))
   allocate(dp3d_tmp(nphys_pts,pver,nelemd))
   allocate(dp3d_tmp_tmp(nphys_pts,pver))
   allocate(phis_tmp(nphys_pts,nelemd))
   allocate(T_tmp(nphys_pts,pver,nelemd))
   allocate(uv_tmp(nphys_pts,2,pver,nelemd))
   allocate(q_tmp(nphys_pts,pver,pcnst,nelemd))
   allocate(omega_tmp(nphys_pts,pver,nelemd))

   call tot_energy_dyn(elem,dyn_out%fvm, 1, nelemd,tl_f , tl_qdp_np0,'dBF')

   if (use_gw_front .or. use_gw_front_igw) then
      allocate(frontgf(nphys_pts,pver,nelemd), stat=ierr)
      if (ierr /= 0) call endrun("dp_coupling: Allocate of frontgf failed.")
      allocate(frontga(nphys_pts,pver,nelemd), stat=ierr)
      if (ierr /= 0) call endrun("dp_coupling: Allocate of frontga failed.")
   end if
   if (use_gw_movmtn_pbl) then
      allocate(vort4gw(nphys_pts,pver,nelemd), stat=ierr)
      if (ierr /= 0) call endrun("dp_coupling: Allocate of vort4gw failed.")
   end if

   if (iam < par%nprocs) then
      if (use_gw_front .or. use_gw_front_igw ) then
         call gws_src_fnct(elem, tl_f, tl_qdp_np0, frontgf, frontga, nphys)
      end if
      if (use_gw_movmtn_pbl ) then
         call gws_src_vort(elem, tl_f, tl_qdp_np0, vort4gw, nphys)
      end if

      if (fv_nphys > 0) then
         call test_mapping_overwrite_dyn_state(elem,dyn_out%fvm)
         !******************************************************************
         ! physics runs on an FVM grid: map GLL vars to physics grid
         !******************************************************************
         call t_startf('dyn2phys')
         ! note that the fvm halo has been filled in prim_run_subcycle
         ! if physics grid resolution is not equal to fvm resolution
         call dyn2phys_all_vars(1,nelemd,elem, dyn_out%fvm,&
              pcnst,hyai(1)*ps0,tl_f,                      &
              ! output
              dp3d_tmp, ps_tmp, q_tmp, T_tmp,              &
              omega_tmp, phis_tmp                          &
              )
         do ie = 1, nelemd
            uv_tmp(:,:,:,ie) = &
               dyn2phys_vector(elem(ie)%state%v(:,:,:,:,tl_f),elem(ie))
         end do
         call t_stopf('dyn2phys')
      else

         !******************************************************************
         ! Physics runs on GLL grid: collect unique points before mapping to
         ! physics decomposition
         !******************************************************************

         if (qsize < pcnst) then
            call endrun('d_p_coupling: Fewer GLL tracers advected than required')
         end if

         call t_startf('UniquePoints')
         do ie = 1, nelemd
           inv_dp3d(:,:,:) = 1.0_r8/elem(ie)%state%dp3d(:,:,:,tl_f)
           do m=1,pcnst
             qgll(:,:,:,m) = elem(ie)%state%Qdp(:,:,:,m,tl_qdp_np0)*inv_dp3d(:,:,:)
           end do
            ncols = elem(ie)%idxP%NumUniquePts
            call UniquePoints(elem(ie)%idxP, elem(ie)%state%psdry(:,:), ps_tmp(1:ncols,ie))
            call UniquePoints(elem(ie)%idxP, nlev, elem(ie)%state%dp3d(:,:,:,tl_f), dp3d_tmp(1:ncols,:,ie))
            call UniquePoints(elem(ie)%idxP, nlev, elem(ie)%state%T(:,:,:,tl_f), T_tmp(1:ncols,:,ie))
            call UniquePoints(elem(ie)%idxV, 2, nlev, elem(ie)%state%V(:,:,:,:,tl_f), uv_tmp(1:ncols,:,:,ie))
            call UniquePoints(elem(ie)%idxV, nlev, elem(ie)%derived%omega, omega_tmp(1:ncols,:,ie))

            call UniquePoints(elem(ie)%idxP, elem(ie)%state%phis, phis_tmp(1:ncols,ie))
            call UniquePoints(elem(ie)%idxP, nlev, pcnst, qgll,Q_tmp(1:ncols,:,:,ie))
         end do
         call t_stopf('UniquePoints')

      end if ! if fv_nphys>0

   else

      ps_tmp(:,:)      = 0._r8
      T_tmp(:,:,:)     = 0._r8
      uv_tmp(:,:,:,:)  = 0._r8
      omega_tmp(:,:,:) = 0._r8
      phis_tmp(:,:)    = 0._r8
      Q_tmp(:,:,:,:)   = 0._r8

      if (use_gw_front .or. use_gw_front_igw) then
         frontgf(:,:,:) = 0._r8
         frontga(:,:,:) = 0._r8
      end if
      if (use_gw_movmtn_pbl) then
         vort4gw(:,:,:) = 0._r8
      end if

   endif ! iam < par%nprocs

   if (fv_nphys < 1) then
      deallocate(qgll)
   end if

   ! q_prev is for saving the tracer fields for calculating tendencies
   if (.not. allocated(q_prev)) then
      allocate(q_prev(pcols,pver,pcnst,begchunk:endchunk))
   end if
   q_prev = 0.0_R8

   call t_startf('dpcopy')
   ! Fill the backend's physics-column arrays (host independent); the copy
   ! into CAM's chunked physics state below reads only from them.
   cols => dyn_out%phys
   if (.not. associated(cols)) call endrun('d_p_coupling: dyn_out%phys not associated')
   call cols%alloc(phys_columns_on_task, pver, pcnst, &
        use_gw_front .or. use_gw_front_igw, use_gw_movmtn_pbl)
   cols%vcoord      = vc_dry_pressure
   cols%q_is_wet(:) = .false.   ! all mixing ratios are dry here
   !$omp parallel do num_threads(max_num_threads) private (col_ind, ie, blk_ind, ilyr, m)
   do col_ind = 1, phys_columns_on_task
      call get_dyn_col_p(col_ind, ie, blk_ind)
      cols%ps(col_ind)   = ps_tmp(blk_ind(1), ie)
      cols%phis(col_ind) = phis_tmp(blk_ind(1), ie)
      do ilyr = 1, pver
         cols%pdel(col_ind, ilyr)  = dp3d_tmp(blk_ind(1), ilyr, ie)
         cols%t(col_ind, ilyr)     = T_tmp(blk_ind(1), ilyr, ie)
         cols%u(col_ind, ilyr)     = uv_tmp(blk_ind(1), 1, ilyr, ie)
         cols%v(col_ind, ilyr)     = uv_tmp(blk_ind(1), 2, ilyr, ie)
         cols%omega(col_ind, ilyr) = omega_tmp(blk_ind(1), ilyr, ie)
         if (cols%has_front) then
            cols%frontgf(col_ind, ilyr) = frontgf(blk_ind(1), ilyr, ie)
            cols%frontga(col_ind, ilyr) = frontga(blk_ind(1), ilyr, ie)
         end if
         if (cols%has_vort) cols%vort4gw(col_ind, ilyr) = vort4gw(blk_ind(1), ilyr, ie)
      end do
      do m = 1, pcnst
         do ilyr = 1, pver
            cols%q(col_ind, ilyr, m) = Q_tmp(blk_ind(1), ilyr, m, ie)
         end do
      end do
   end do


   call t_stopf('dpcopy')

   ! Save the tracer fields input to physics package for calculating tendencies
   ! The mixing ratios are all dry at this point.
   !$omp parallel do num_threads(max_num_threads) private (col_ind, lchnk, icol, ilyr, m)
   do col_ind = 1, phys_columns_on_task
      call get_chunk_info_p(col_ind, lchnk, icol)
      do m = 1, pcnst
         do ilyr = 1, pver
            q_prev(icol, ilyr, m, lchnk) = cols%q(col_ind, ilyr, m)
         end do
      end do
   end do
   ! The debug_coupling output (test_mapping_output_phys_state) now runs on
   ! CAM's side after the copy, through se_dynamics_t%debug_phys_state.

   ! Deallocate the temporary arrays
   deallocate(ps_tmp)
   deallocate(dp3d_tmp)
   deallocate(phis_tmp)
   deallocate(T_tmp)
   deallocate(uv_tmp)
   deallocate(q_tmp)
   deallocate(omega_tmp)

end subroutine d_p_coupling

!=========================================================================================

subroutine debug_overwrite_tendencies(phys_state, phys_tend, fvm)
   ! CAM-SE coupling test (does nothing unless cpp macro debug_coupling is
   ! defined): overwrite the physics tendencies with test functions.  Called
   ! by CAM's side before it fills the return arrays.
   use test_fvm_mapping_impl, only: test_mapping_overwrite_tendencies
   use fvm_control_volume_mod, only: fvm_struct
   type(physics_state), intent(inout) :: phys_state(begchunk:endchunk)
   type(physics_tend),  intent(inout) :: phys_tend(begchunk:endchunk)
   type(fvm_struct),    intent(inout) :: fvm(:)
   integer :: lchnk, ncols
   do lchnk = begchunk, endchunk
      ncols = get_ncols_p(lchnk)
      call test_mapping_overwrite_tendencies(phys_state(lchnk), phys_tend(lchnk), ncols, lchnk, &
           q_prev(1:ncols,:,:,lchnk), fvm)
   end do
end subroutine debug_overwrite_tendencies

!=========================================================================================

subroutine p_d_coupling(dyn_in, tl_f, tl_qdp)

   ! Convert the physics output state into the dynamics input state.

   use phys_grid,        only: get_dyn_col_p, columns_on_task, get_chunk_info_p, phys_columns_on_task
   use bndry_mod,        only: bndry_exchange
   use edge_mod,         only: edgeVpack, edgeVunpack
   use fvm_mapping,      only: phys2dyn_forcings_fvm
   use test_fvm_mapping_impl, only: test_mapping_output_mapped_tendencies
   use dimensions_mod,   only: use_cslam
   use control_mod,      only: gll_advect_q
   ! arguments
   integer,             intent(in)                                  :: tl_qdp, tl_f
   type(dyn_import_t),  intent(inout)                               :: dyn_in
   type(hybrid_t)                                                   :: hybrid

   ! LOCAL VARIABLES
   integer                      :: ncols             ! index
   type(element_t), pointer     :: elem(:)           ! pointer to dyn_in element array
   integer                      :: ie                ! index for elements
   integer                      :: col_ind           ! index over columns
   integer                      :: blk_ind(1)        ! element offset
   integer                      :: lchnk, icol, ilyr ! indices for chunk, column, layer

   real (kind=r8),  allocatable :: dp_phys(:,:,:)    ! temp array to hold dp on physics grid
   real (kind=r8),  allocatable :: T_tmp(:,:,:)      ! temp array to hold T
   real (kind=r8),  allocatable :: dq_tmp(:,:,:,:)   ! temp array to hold q
   real (kind=r8),  allocatable :: uv_tmp(:,:,:,:)   ! temp array to hold uv
   integer                      :: m, i, j, k

   type(phys_cols_t),   pointer :: cols
   integer                      :: num_trac
   integer                      :: nets, nete
   integer                      :: kptr, ii
   !----------------------------------------------------------------------------

   if (.not. local_dp_map) then
      call endrun('p_d_coupling: Weak scaling does not support load balancing')
   end if

   if (iam < par%nprocs) then
      elem => dyn_in%elem
   else
      nullify(elem)
   end if

   allocate(T_tmp(nphys_pts,pver,nelemd))
   allocate(uv_tmp(nphys_pts,2,pver,nelemd))
   allocate(dq_tmp(nphys_pts,pver,pcnst,nelemd))
   allocate(dp_phys(nphys_pts,pver,nelemd))

   T_tmp  = 0.0_r8
   uv_tmp = 0.0_r8
   dq_tmp = 0.0_r8

   if (.not. allocated(q_prev)) then
      call endrun('p_d_coupling: q_prev not allocated')
   end if


   ! Wet-to-dry conversion and the debug_coupling tendency test are done on
   ! CAM's side (dp_coupling_cam), which fills cols%dtdt/dudt/dvdt/q_out/pdeldry_out.
   cols => dyn_in%phys
   call t_startf('pd_copy')
   !$omp parallel do num_threads(max_num_threads) private (col_ind, lchnk, icol, ie, blk_ind, ilyr, m)
   do col_ind = 1, phys_columns_on_task
      call get_dyn_col_p(col_ind, ie, blk_ind)
      call get_chunk_info_p(col_ind, lchnk, icol)


      do ilyr = 1, pver
         dp_phys(blk_ind(1),ilyr,ie)  = cols%pdeldry_out(col_ind,ilyr)
         T_tmp(blk_ind(1),ilyr,ie)    = cols%dtdt(col_ind,ilyr)
         uv_tmp(blk_ind(1),1,ilyr,ie) = cols%dudt(col_ind,ilyr)
         uv_tmp(blk_ind(1),2,ilyr,ie) = cols%dvdt(col_ind,ilyr)
         do m = 1, pcnst
            dq_tmp(blk_ind(1),ilyr,m,ie) =                                    &
                 (cols%q_out(col_ind,ilyr,m) - q_prev(icol,ilyr,m,lchnk))
         end do
      end do
   end do
   call t_stopf('pd_copy')

   if (iam < par%nprocs) then

      if (fv_nphys > 0) then

         ! put forcings into fvm structure
         num_trac = max(qsize,ntrac)
         do ie = 1, nelemd
            do j = 1, fv_nphys
               do i = 1, fv_nphys
                  ii = i + (j-1)*fv_nphys
                  dyn_in%fvm(ie)%ft(i,j,1:pver)                 = T_tmp(ii,1:pver,ie)
                  dyn_in%fvm(ie)%fm(i,j,1:2,1:pver)             = uv_tmp(ii,1:2,1:pver,ie)
                  dyn_in%fvm(ie)%fc_phys(i,j,1:pver,1:num_trac) = dq_tmp(ii,1:pver,1:num_trac,ie)
                  dyn_in%fvm(ie)%dp_phys(i,j,1:pver)            = dp_phys(ii,1:pver,ie)
               end do
            end do
         end do

         !JMD $OMP PARALLEL NUM_THREADS(horz_num_threads), DEFAULT(SHARED), PRIVATE(hybrid,nets,nete,n)
         !JMD        hybrid = config_thread_region(par,'horizontal')
         hybrid = config_thread_region(par,'serial')
         call get_loop_ranges(hybrid,ibeg=nets,iend=nete)
         !
         ! high-order mapping of ft and fm using fvm technology
         !
         call t_startf('phys2dyn')
         call phys2dyn_forcings_fvm(elem, dyn_in%fvm, hybrid,nets,nete,ntrac==0, tl_f, tl_qdp)
         call t_stopf('phys2dyn')
      else

         call t_startf('putUniquePoints')

         !$omp parallel do num_threads(max_num_threads) private(ie,ncols)
         do ie = 1, nelemd
            ncols = elem(ie)%idxP%NumUniquePts
            call putUniquePoints(elem(ie)%idxP, nlev, T_tmp(1:ncols,:,ie),       &
               elem(ie)%derived%fT(:,:,:))
            call putUniquePoints(elem(ie)%idxV, 2, nlev, uv_tmp(1:ncols,:,:,ie), &
               elem(ie)%derived%fM(:,:,:,:))
            call putUniquePoints(elem(ie)%idxV, nlev,pcnst, dq_tmp(1:ncols,:,:,ie), &
               elem(ie)%derived%fQ(:,:,:,:))
         end do
         call t_stopf('putUniquePoints')
      end if
   end if

   deallocate(T_tmp)
   deallocate(uv_tmp)
   deallocate(dq_tmp)

   ! Boundary exchange for physics forcing terms.
   ! For physics on GLL grid, for points with duplicate degrees of freedom,
   ! putuniquepoints() set one of the element values and set the others to zero,
   ! so do a simple sum (boundary exchange with no weights).
   ! For physics grid, we interpolated into all points, so do weighted average.

   call t_startf('p_d_coupling:bndry_exchange')

   do ie = 1, nelemd
      if (fv_nphys > 0) then
         do k = 1, nlev
            dyn_in%elem(ie)%derived%FM(:,:,1,k) =                          &
                 dyn_in%elem(ie)%derived%FM(:,:,1,k) *                     &
                 dyn_in%elem(ie)%spheremp(:,:)
            dyn_in%elem(ie)%derived%FM(:,:,2,k) =                          &
                 dyn_in%elem(ie)%derived%FM(:,:,2,k) *                     &
                 dyn_in%elem(ie)%spheremp(:,:)
            dyn_in%elem(ie)%derived%FT(:,:,k) =                            &
                 dyn_in%elem(ie)%derived%FT(:,:,k) *                       &
                 dyn_in%elem(ie)%spheremp(:,:)
            if (gll_advect_q) then
               do m = 1, qsize
                  dyn_in%elem(ie)%derived%FQ(:,:,k,m) =                    &
                       dyn_in%elem(ie)%derived%FQ(:,:,k,m) *               &
                       dyn_in%elem(ie)%spheremp(:,:)
               end do
            end if
         end do
      end if
      kptr = 0
      call edgeVpack(edgebuf, dyn_in%elem(ie)%derived%FM(:,:,:,:), 2*nlev, kptr, ie)
      kptr = kptr + 2*nlev
      call edgeVpack(edgebuf, dyn_in%elem(ie)%derived%FT(:,:,:), nlev, kptr, ie)
      if (.not. use_cslam .or. gll_advect_q) then
         !
         ! if using CSLAM qdp is being overwritten with CSLAM values in the dynamics
         ! so no need to do boundary exchange of tracer tendency on GLL grid here
         ! (unless gll_advect_q, in which case the GLL tracers are forced via FQ too)
         !
         kptr = kptr + nlev
         call edgeVpack(edgebuf, dyn_in%elem(ie)%derived%FQ(:,:,:,:), nlev*qsize, kptr, ie)
      end if
   end do

   if (iam < par%nprocs) then
     call bndry_exchange(par, edgebuf, location='p_d_coupling')
   end if

   do ie = 1, nelemd
      kptr = 0
      call edgeVunpack(edgebuf, dyn_in%elem(ie)%derived%FM(:,:,:,:), 2*nlev, kptr, ie)
      kptr = kptr + 2*nlev
      call edgeVunpack(edgebuf, dyn_in%elem(ie)%derived%FT(:,:,:), nlev, kptr, ie)
      kptr = kptr + nlev
      if (.not. use_cslam .or. gll_advect_q) then
         call edgeVunpack(edgebuf, dyn_in%elem(ie)%derived%FQ(:,:,:,:), nlev*qsize, kptr, ie)
      end if
      if (fv_nphys > 0) then
         do k = 1, nlev
            dyn_in%elem(ie)%derived%FM(:,:,1,k) =                             &
                 dyn_in%elem(ie)%derived%FM(:,:,1,k) *                        &
                 dyn_in%elem(ie)%rspheremp(:,:)
            dyn_in%elem(ie)%derived%FM(:,:,2,k) =                             &
                 dyn_in%elem(ie)%derived%FM(:,:,2,k) *                        &
                 dyn_in%elem(ie)%rspheremp(:,:)
            dyn_in%elem(ie)%derived%FT(:,:,k) =                               &
                 dyn_in%elem(ie)%derived%FT(:,:,k) *                          &
                 dyn_in%elem(ie)%rspheremp(:,:)
            if (gll_advect_q) then
               do m = 1, qsize
                  dyn_in%elem(ie)%derived%FQ(:,:,k,m) =                       &
                       dyn_in%elem(ie)%derived%FQ(:,:,k,m) *                  &
                       dyn_in%elem(ie)%rspheremp(:,:)
               end do
            end if
         end do
      end if
   end do
   call t_stopf('p_d_coupling:bndry_exchange')

   if (iam < par%nprocs .and. fv_nphys > 0) then
      call test_mapping_output_mapped_tendencies(dyn_in%fvm(1:nelemd), elem(1:nelemd), &
                                                 1, nelemd, tl_f, tl_qdp)
   end if
end subroutine p_d_coupling

!=========================================================================================
end module dp_coupling_impl
