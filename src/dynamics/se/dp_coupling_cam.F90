module dp_coupling_cam
   !----------------------------------------------------------------------------
   ! CAM's side of dynamics-to-physics coupling for the SE dycores.
   !
   ! The backend (se/dycore or se/cam_mods) fills its physics-column arrays
   ! (se_backend_t%phys, a phys_cols_t) and returns.  This module copies them
   ! into CAM's chunked physics_state and computes CAM's derived quantities.
   ! Which derivation applies is decided by what the columns hold
   ! (phys%vcoord), not by which backend was compiled:
   !   vc_dry_pressure   -> derived_phys_dry   (formerly in se/dycore)
   !   vc_moist_pressure -> derived_phys_moist (formerly derived_phys in se/cam_mods)
   ! Both routines are moved unchanged (nlev -> pver).
   !----------------------------------------------------------------------------
   use shr_kind_mod,    only: r8 => shr_kind_r8
   use ppgrid,          only: begchunk, endchunk, pcols, pver, pverp
   use constituents,    only: pcnst, cnst_type
   use physics_types,   only: physics_state, physics_tend, physics_cnst_limit
   use phys_grid,       only: get_ncols_p, get_chunk_info_p, phys_columns_on_task
   use air_composition, only: rairv, cpairv
   use qneg_module,     only: qneg3
   use perf_mod,        only: t_startf, t_stopf
   use cam_abortutils,  only: endrun
   use cam_logfile,     only: iulog
   use phys_cols_mod,   only: phys_cols_t
   use se_backend_mod,  only: se_backend_t

   implicit none
   private
   save

   public :: phys_state_from_backend   ! copy + derived quantities (d_p_coupling, stepon_run1)
   public :: aerosol_state_from_phys   ! aerosol state update from phys_state%q (stepon_run1)

contains

   subroutine phys_state_from_backend(phys_state, phys_tend, dyn_out)
      use dyn_tests_utils, only: vc_dry_pressure, vc_moist_pressure
      type(physics_state), intent(inout) :: phys_state(begchunk:endchunk)
      type(physics_tend ), intent(inout) :: phys_tend(begchunk:endchunk)
      class(se_backend_t), intent(inout) :: dyn_out

      type(phys_cols_t), pointer :: cols
      integer :: col_ind, lchnk, icol, ilyr, m, ncols
      logical :: copy_pdel

      cols => dyn_out%phys
      if (.not. associated(cols)) call endrun('phys_state_from_backend: dyn_out%phys not associated')
      if (cols%ncnst /= pcnst .or. cols%nlev /= pver) &
           call endrun('phys_state_from_backend: phys columns do not match pcnst/pver')
      ! The dry-pressure dycore provides pdel; the moist one leaves it to derived_phys_moist
      copy_pdel = (cols%vcoord == vc_dry_pressure)

      call t_startf('dpcopy_cam')
      !$omp parallel do private (col_ind, lchnk, icol, ilyr, m)
      do col_ind = 1, phys_columns_on_task
         call get_chunk_info_p(col_ind, lchnk, icol)
         phys_state(lchnk)%ps(icol)   = cols%ps(col_ind)
         phys_state(lchnk)%phis(icol) = cols%phis(col_ind)
         do ilyr = 1, pver
            if (copy_pdel) phys_state(lchnk)%pdel(icol, ilyr) = cols%pdel(col_ind, ilyr)
            phys_state(lchnk)%t(icol, ilyr)     = cols%t(col_ind, ilyr)
            phys_state(lchnk)%u(icol, ilyr)     = cols%u(col_ind, ilyr)
            phys_state(lchnk)%v(icol, ilyr)     = cols%v(col_ind, ilyr)
            phys_state(lchnk)%omega(icol, ilyr) = cols%omega(col_ind, ilyr)
         end do
         do m = 1, pcnst
            do ilyr = 1, pver
               phys_state(lchnk)%q(icol, ilyr, m) = cols%q(col_ind, ilyr, m)
            end do
         end do
      end do
      call t_stopf('dpcopy_cam')

      ! CAM-SE's coupling test output (-Ddebug_coupling); a no-op otherwise
      call dyn_out%debug_phys_state(phys_state)

      call t_startf('derived_phys')
      if (cols%vcoord == vc_dry_pressure) then
         call derived_phys_dry(phys_state, phys_tend)
      else if (cols%vcoord == vc_moist_pressure) then
         call derived_phys_moist(phys_state, phys_tend)
      else
         call endrun('phys_state_from_backend: unknown vertical coordinate in phys%vcoord')
      end if
      call t_stopf('derived_phys')

      ! preqx delivers omega/p; physics wants omega
      if (cols%omega_over_p) then
         !$omp parallel do private (lchnk, ncols, ilyr, icol)
         do lchnk = begchunk, endchunk
            ncols = get_ncols_p(lchnk)
            do ilyr = 1, pver
               do icol = 1, ncols
                  phys_state(lchnk)%omega(icol,ilyr) = phys_state(lchnk)%omega(icol,ilyr) &
                                                      *phys_state(lchnk)%pmid(icol,ilyr)
               end do
            end do
         end do
      end if

      ! Unused columns of each chunk
      !$omp parallel do private (lchnk, ncols)
      do lchnk = begchunk, endchunk
         ncols = get_ncols_p(lchnk)
         if (pcols > ncols) then
            phys_state(lchnk)%phis(ncols+1:) = 0.0_r8
         end if
      end do
   end subroutine phys_state_from_backend

   subroutine aerosol_state_from_phys(phys_state)
      ! Update the aerosol state objects from the physics state constituents
      ! (moved from the backends' stepon_run1).
      use aerosol_properties_mod, only: aerosol_properties
      use aerosol_state_mod,      only: aerosol_state
      use microp_aero,            only: aerosol_state_object, aerosol_properties_object
      type(physics_state), intent(in) :: phys_state(begchunk:endchunk)

      class(aerosol_properties), pointer :: aero_props_obj
      class(aerosol_state),      pointer :: aero_state_obj
      logical, save :: checked = .false., aerosols_transported = .false.
      integer :: c

      if (.not. checked) then
         aero_props_obj => aerosol_properties_object()
         if (associated(aero_props_obj)) then
            aerosols_transported = aero_props_obj%number_transported() > 0
         end if
         checked = .true.
      end if
      if (aerosols_transported) then
         do c = begchunk, endchunk
            aero_state_obj => aerosol_state_object(c)
            call aero_state_obj%set_transported(phys_state(c)%q)
         end do
      end if
   end subroutine aerosol_state_from_phys

!=========================================================================================
! Moved unchanged from se/dycore/dp_coupling_impl.F90 (nlev -> pver)
!=========================================================================================

subroutine derived_phys_dry(phys_state, phys_tend)

   ! The ps, pdel, and q components of phys_state are all dry on input.
   ! On output the psdry and pdeldry components are initialized; ps and pdel are
   ! updated to contain contribution from water vapor only; the 'wet' constituent
   ! mixing ratios are converted to a wet basis.  Initialize geopotential heights.
   ! Finally compute energy and water column integrals of the physics input state.

   use constituents,    only: qmin
   use physconst,       only: gravit, zvir
   use cam_thermo,      only: cam_thermo_dry_air_update, cam_thermo_water_update
   use air_composition, only: thermodynamic_active_species_num
   use air_composition, only: thermodynamic_active_species_idx
   use air_composition, only: cpairv, rairv, cappav, dry_air_species_num
   use shr_const_mod,   only: shr_const_rwv
   use phys_control,    only: waccmx_is
   use geopotential,    only: geopotential_t
   use static_energy,   only: update_dry_static_energy_run
   use hycoef,          only: hyai, ps0
   use shr_vmath_mod,   only: shr_vmath_log
   use qneg_module,     only: qneg3
   use dyn_tests_utils, only: vc_dry_pressure
   use shr_kind_mod,    only: shr_kind_cx
   ! arguments
   type(physics_state), intent(inout), dimension(begchunk:endchunk) :: phys_state
   type(physics_tend ), intent(inout), dimension(begchunk:endchunk) :: phys_tend

   ! local variables
   integer :: lchnk

   real(r8) :: zvirv(pcols,pver)    ! Local zvir array pointer
   real(r8) :: factor_array(pcols,pver)

   integer :: m, i, k, ncol, m_cnst

   !Needed for "update_dry_static_energy" CCPP scheme
   integer :: errflg
   character(len=shr_kind_cx) :: errmsg
   !----------------------------------------------------------------------------

   ! Evaluate derived quantities

   !!$omp parallel do num_threads(horz_num_threads) private (lchnk, ncol, k, i, m , zvirv, factor_array)
   do lchnk = begchunk,endchunk

      ncol = get_ncols_p(lchnk)

      ! dry pressure variables

      do i = 1, ncol
         phys_state(lchnk)%psdry(i) = hyai(1)*ps0 + sum(phys_state(lchnk)%pdel(i,:))
      end do
      do i = 1, ncol
         phys_state(lchnk)%pintdry(i,1) = hyai(1)*ps0
      end do
      call shr_vmath_log(phys_state(lchnk)%pintdry(1:ncol,1), &
                         phys_state(lchnk)%lnpintdry(1:ncol,1),ncol)
      do k = 1, pver
         do i = 1, ncol
            phys_state(lchnk)%pintdry(i,k+1) = phys_state(lchnk)%pintdry(i,k) + &
                                               phys_state(lchnk)%pdel(i,k)
         end do
         call shr_vmath_log(phys_state(lchnk)%pintdry(1:ncol,k+1),&
                            phys_state(lchnk)%lnpintdry(1:ncol,k+1),ncol)
      end do

      do k=1,pver
         do i=1,ncol
            phys_state(lchnk)%pdeldry (i,k) = phys_state(lchnk)%pdel(i,k)
            phys_state(lchnk)%rpdeldry(i,k) = 1._r8/phys_state(lchnk)%pdeldry(i,k)
            phys_state(lchnk)%pmiddry (i,k) = 0.5D0*(phys_state(lchnk)%pintdry(i,k+1) + &
                                                     phys_state(lchnk)%pintdry(i,k))
         end do
         call shr_vmath_log(phys_state(lchnk)%pmiddry(1:ncol,k), &
                            phys_state(lchnk)%lnpmiddry(1:ncol,k),ncol)
      end do

      ! wet pressure variables (should be removed from physics!)
      factor_array(:,:) = 1.0_r8
      do m_cnst=dry_air_species_num + 1,thermodynamic_active_species_num
        m = thermodynamic_active_species_idx(m_cnst)
        do k=1,pver
          do i=1,ncol
            ! at this point all q's are dry
            factor_array(i,k) = factor_array(i,k)+phys_state(lchnk)%q(i,k,m)
          end do
        end do
      end do

      do k=1,pver
         do i=1,ncol
            phys_state(lchnk)%pdel (i,k) = phys_state(lchnk)%pdeldry(i,k)*factor_array(i,k)
         end do
      end do

      ! initialize vertical loop - model top pressure

      do i=1,ncol
         phys_state(lchnk)%ps(i)     = phys_state(lchnk)%pintdry(i,1)
         phys_state(lchnk)%pint(i,1) = phys_state(lchnk)%pintdry(i,1)
      end do
      do k = 1, pver
         do i=1,ncol
            phys_state(lchnk)%pint(i,k+1) =  phys_state(lchnk)%pint(i,k)+phys_state(lchnk)%pdel(i,k)
            phys_state(lchnk)%pmid(i,k)   = (phys_state(lchnk)%pint(i,k+1)+phys_state(lchnk)%pint(i,k))/2._r8
            phys_state(lchnk)%ps  (i)     =  phys_state(lchnk)%ps(i) + phys_state(lchnk)%pdel(i,k)
         end do
         call shr_vmath_log(phys_state(lchnk)%pint(1:ncol,k),phys_state(lchnk)%lnpint(1:ncol,k),ncol)
         call shr_vmath_log(phys_state(lchnk)%pmid(1:ncol,k),phys_state(lchnk)%lnpmid(1:ncol,k),ncol)
      end do
      call shr_vmath_log(phys_state(lchnk)%pint(1:ncol,pverp),phys_state(lchnk)%lnpint(1:ncol,pverp),ncol)

      do k = 1, pver
         do i = 1, ncol
            phys_state(lchnk)%rpdel(i,k)  = 1._r8/phys_state(lchnk)%pdel(i,k)
         end do
      end do

      !------------------------------------------------------------
      ! Apply limiters to mixing ratios of major species (waccmx)
      !------------------------------------------------------------
      if (dry_air_species_num>0) then
        call physics_cnst_limit( phys_state(lchnk) )
        !-----------------------------------------------------------------------------
        ! Call cam_thermo_dry_air_update to compute cpairv, rairv, mbarv, and cappav as
        ! constituent dependent variables.
        ! Compute molecular viscosity(kmvis) and conductivity(kmcnd).
        ! Fill local zvirv variable; calculated for WACCM-X.
        !-----------------------------------------------------------------------------
        call cam_thermo_dry_air_update(phys_state(lchnk)%q, phys_state(lchnk)%t, lchnk, ncol)
        zvirv(:,:) = shr_const_rwv / rairv(:,:,lchnk) -1._r8
      else
        zvirv(:,:) = zvir
      end if
      !
      ! update cp_dycore in module air_composition.
      ! (note: at this point q is dry)
      !
      call cam_thermo_water_update(phys_state(lchnk)%q(1:ncol,:,:), lchnk, ncol, vc_dry_pressure)
      do k = 1, pver
         do i = 1, ncol
            phys_state(lchnk)%exner(i,k) = (phys_state(lchnk)%pint(i,pver+1) &
                                            / phys_state(lchnk)%pmid(i,k))**cappav(i,k,lchnk)
         end do
      end do
      !
      ! CAM physics: water tracers are moist; the rest dry
      !
      factor_array(1:ncol,1:pver) = 1._r8/factor_array(1:ncol,1:pver)
      do m = 1,pcnst
         if (cnst_type(m) == 'wet') then
            do k = 1, pver
               do i = 1, ncol
                  phys_state(lchnk)%q(i,k,m) = factor_array(i,k)*phys_state(lchnk)%q(i,k,m)
               end do
            end do
         end if
      end do

      ! Ensure tracers are all positive
      call qneg3('D_P_COUPLING',lchnk  ,ncol    ,pcols   ,pver    , &
           1, pcnst, qmin  ,phys_state(lchnk)%q)

      ! Compute initial geopotential heights - based on full pressure
      call geopotential_t(phys_state(lchnk)%lnpint, phys_state(lchnk)%lnpmid  , phys_state(lchnk)%pint, &
         phys_state(lchnk)%pmid  , phys_state(lchnk)%pdel    , phys_state(lchnk)%rpdel                , &
         phys_state(lchnk)%t     , phys_state(lchnk)%q(:,:,:), rairv(:,:,lchnk), gravit, zvirv        , &
         phys_state(lchnk)%zi    , phys_state(lchnk)%zm      , ncol)
      ! Compute initial dry static energy, include surface geopotential
      call update_dry_static_energy_run(pver, gravit, phys_state(lchnk)%t(1:ncol,:),  &
                                        phys_state(lchnk)%zm(1:ncol,:),               &
                                        phys_state(lchnk)%phis(1:ncol),               &
                                        phys_state(lchnk)%s(1:ncol,:),                &
                                        cpairv(1:ncol,:,lchnk), errflg, errmsg)
      ! Energy and water integrals of the input state: CAM side (d_p_coupling_cam_finish)


   end do  ! lchnk

end subroutine derived_phys_dry

!=========================================================================================
! Moved unchanged from se/cam_mods/dp_coupling_impl.F90, where it was derived_phys
! (nlev -> pver)
!=========================================================================================

  subroutine derived_phys_moist(phys_state, phys_tend)
    use constituents,  only: qmin
    use physconst,     only: gravit, zvir, cappa
    use spmd_utils,    only: masterproc
    use ppgrid,        only: pver
    use geopotential,  only: geopotential_t
    use static_energy, only: update_dry_static_energy_run
    use physics_types, only: set_state_pdry, set_wet_to_dry
    use hycoef,   only : hyam, hybm, hyai, hybi, ps0
    use shr_vmath_mod, only: shr_vmath_log
    use gmean_mod,     only: gmean
    use shr_kind_mod,  only: shr_kind_cx


    implicit none
    type(physics_state), intent(inout), dimension(begchunk:endchunk) :: phys_state
    type(physics_tend ), intent(inout), dimension(begchunk:endchunk) :: phys_tend

    integer :: lchnk
    real(r8) :: qbot                 ! bottom level q before change
    real(r8) :: qbotm1               ! bottom-1 level q before change
    real(r8) :: dqreq                ! q change at pver-1 required to remove q<qmin at pver
    real(r8) :: qmavl                ! available q at level pver-1

    real(r8) :: ke(pcols,begchunk:endchunk)
    real(r8) :: se(pcols,begchunk:endchunk)
    real(r8) :: ke_glob(1),se_glob(1)
    real(r8) :: zvirv(pcols,pver)    ! Local zvir array pointer

    integer :: m, i, k, ncol

    ! Needed for "update_dry_static_energy" CCPP scheme
    integer :: errflg
    character(len=shr_kind_cx) :: errmsg


!
! Evaluate derived quantities
!
!$omp parallel do private (lchnk, ncol, k, i, zvirv, errflg, errmsg)
    do lchnk = begchunk,endchunk
       ncol = get_ncols_p(lchnk)
       do k=1,pver
          do i=1,ncol
             phys_state(lchnk)%pint(i,k)=hyai(k)*ps0+hybi(k)*phys_state(lchnk)%ps(i)
             phys_state(lchnk)%pmid(i,k)=hyam(k)*ps0+hybm(k)*phys_state(lchnk)%ps(i)
          end do
        call shr_vmath_log(phys_state(lchnk)%pint(1:ncol,k), &
                           phys_state(lchnk)%lnpint(1:ncol,k),ncol)
        call shr_vmath_log(phys_state(lchnk)%pmid(1:ncol,k), &
                           phys_state(lchnk)%lnpmid(1:ncol,k),ncol)
       end do
       do i=1,ncol
          phys_state(lchnk)%pint(i,pverp)=hyai(pverp)*ps0+hybi(pverp)*phys_state(lchnk)%ps(i)
       end do
      call shr_vmath_log(phys_state(lchnk)%pint(1:ncol,pverp), &
                         phys_state(lchnk)%lnpint(1:ncol,pverp),ncol)

       do k=1,pver
          do i=1,ncol
          phys_state(lchnk)%pdel (i,k)  = phys_state(lchnk)%pint(i,k+1) &
                                         -phys_state(lchnk)%pint(i,k)
             phys_state(lchnk)%rpdel(i,k) = 1._r8/phys_state(lchnk)%pdel(i,k)
             phys_state(lchnk)%exner (i,k) = (phys_state(lchnk)%pint(i,pver+1) &
                                             / phys_state(lchnk)%pmid(i,k))**cappa
          end do
       end do

!-----------------------------------------------------------------------------------
!  Need to fill zvirv 2D variables to be compatible with geopotential_t interface
!-----------------------------------------------------------------------------------
       zvirv(:,:) = zvir
!
! Compute initial geopotential heights
!

! NOTE:  if a tracer is marked "dry", that means physics wants it dry
!        if dycore advects it wet, it should be converted here
!        FV dycore does this, and in physics/cam/tphysac.F90 it will
!        be converted back to wet, BUT ONLY FOR FV dycore
!
!        EUL: advects all tracers (except q1) as dry.  so it never
!        calls this.
!
!        SE:  we should follow FV and advect all tracers wet (especially
!        since we will be switching to conservation form of advection).
!        So this is broken since dry tracers will never get converted
!        back to wet. (in APE, all tracers are wet, so it is ok for now)
!
! Convert dry type constituents from moist to dry mixing ratio
       call set_state_pdry(phys_state(lchnk), .true.)	 ! First get dry pressure to use for this timestep
       call set_wet_to_dry(phys_state(lchnk), convert_cnst_type='dry')    ! Dynamics had moist, physics wants dry.
!
! Ensure tracers are all positive
!
!jt       Could cause differences if checking BFB
       call qneg3('D_P_COUPLING',lchnk  ,ncol    ,pcols   ,pver    , &
                  1, pcnst, qmin  ,phys_state(lchnk)%q)

      ! Compute initial geopotential heights
      call geopotential_t(phys_state(lchnk)%lnpint, phys_state(lchnk)%lnpmid  ,&
                          phys_state(lchnk)%pint  , phys_state(lchnk)%pmid    ,&
                          phys_state(lchnk)%pdel  , phys_state(lchnk)%rpdel   ,&
                          phys_state(lchnk)%t     , phys_state(lchnk)%q       ,&
                          rairv(:,:,lchnk)        , gravit, zvirv             ,&
                          phys_state(lchnk)%zi    , phys_state(lchnk)%zm      ,&
                          ncol)

       ! Compute initial dry static energy s = g*z + c_p*T
       call update_dry_static_energy_run(pver, gravit, phys_state(lchnk)%t(1:ncol,:),  &
                                         phys_state(lchnk)%zm(1:ncol,:),               &
                                         phys_state(lchnk)%phis(1:ncol),               &
                                         phys_state(lchnk)%s(1:ncol,:),                &
                                         cpairv(1:ncol,:,lchnk), errflg, errmsg)


! Energy and water integrals of the input state: CAM side (d_p_coupling_cam_finish)


#if 0
       ke(:,lchnk) = 0._r8
       se(:,lchnk) = 0._r8
!       wv = 0._r8
!       wl = 0._r8
!       wi = 0._r8
       do k = 1, pver
          do i = 1, ncol
             ke(i,lchnk) = ke(i,lchnk) + ( 0.5_r8*(phys_state(lchnk)%u(i,k)**2 + &
                  phys_state(lchnk)%v(i,k)**2)*phys_state(lchnk)%pdel(i,k) )/gravit
             se(i,lchnk) = se(i,lchnk) + phys_state(lchnk)%s(i,k         )*phys_state(lchnk)%pdel(i,k)/gravit
!             wv = wv + phys_state(lchnk)%q(i,k,1       )*phys_state(lchnk)%pdel(i,k)
!             wl = wl + phys_state(lchnk)%q(i,k,ixcldliq)*phys_state(lchnk)%pdel(i,k)
!             wi = wi + phys_state(lchnk)%q(i,k,ixcldice)*phys_state(lchnk)%pdel(i,k)
          end do
       end do
#endif
    end do

#if 0
! This wont match SE exactly.  SE computes KE at half levels
! SE includes cp_star (physics SE uses cp )
! CAM stdout of total energy also includes latent energy of Q,Q1,Q2
! making it a little larger
    call gmean(ke,ke_glob,1)
    call gmean(se,se_glob,1)
    if (masterproc) then
       write(iulog,'(a,e20.8,a,e20.8)') 'KE = ',ke_glob(1),' SE = ',se_glob(1)
       write(iulog,'(a,e20.8)') 'TOTE = ',ke_glob(1)+se_glob(1)
    endif
#endif

  end subroutine derived_phys_moist

end module dp_coupling_cam
