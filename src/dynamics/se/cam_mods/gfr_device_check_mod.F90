#if (defined(SE_GFR_DEVICE_CHECK) || defined(SE_GFR_DEVICE_COUPLING)) && defined(KOKKOS_TARGET)
#define SE_GFR_DEVICE
#endif
module gfr_device_check_mod
  !-----------------------------------------------------------------------------
  ! Shadow check for the device GLL -> physics-grid remap (stormspeed_sima,
  ! GPU work step 1).  Compiled in only with -DSE_GFR_DEVICE_CHECK (and the
  ! Kokkos target); the model's answers are unchanged.
  !
  ! Each d_p_coupling: run Hommexx's GllFvRemap on the device into arrays owned
  ! here, gather them on the device into physics-column order (the layout of
  ! phys_cols_t), and compare with the Fortran remap that filled phys.  The
  ! differences (C++ vs Fortran remap: roundoff expected) go to the atm log.
  !-----------------------------------------------------------------------------
  use shr_kind_mod,   only: r8 => shr_kind_r8
  use phys_cols_mod,  only: phys_cols_t
  implicit none
  private
  public :: gfr_device_check_init, gfr_device_check, gfr_device_check_return
  public :: gfr_device_fill_cols        ! -DSE_GFR_DEVICE_COUPLING: device remap fills phys
  public :: gfr_device_forcing_from_cols ! -DSE_GFR_DEVICE_COUPLING: device remap gives the forcing
  public :: gfr_device_forcing_on_device ! -DSE_GFR_DEVICE_RESIDENT: forcing stays on the device

#if defined(SE_GFR_DEVICE)
  real(r8), allocatable :: ps_c(:,:), phis_c(:,:)        ! (npg, nelemd)
  real(r8), allocatable :: t_c(:,:,:), om_c(:,:,:)       ! (nlev, npg, nelemd)
  real(r8), allocatable :: uv_c(:,:,:,:)                 ! (nlev, 2, npg, nelemd)
  real(r8), allocatable :: q_c(:,:,:,:)                  ! (nlev, nq, npg, nelemd)
  real(r8), allocatable :: g2(:,:,:)                     ! gathered (ncol, nlev, 4): t, u, v, omega
  real(r8), allocatable :: g1(:,:)                       ! gathered (ncol, 2): ps, phis
  real(r8), allocatable :: gq(:,:,:)                     ! gathered (ncol, nlev, nq)
  integer,  allocatable :: map_ie(:), map_col(:)         ! physics column -> (ie, col)
  integer :: npg = 0, nq = 0, ncall = 0, ncall_r = 0
  ! return path (physics -> dynamics)
  real(r8), allocatable :: rt_c(:,:,:), ruv_c(:,:,:,:), rq_c(:,:,:,:)   ! GllFvRemap layout
  real(r8), allocatable :: r2(:,:,:), rq(:,:,:)                         ! (ncol,nlev,3), (ncol,nlev,nq)

  interface
     subroutine gfr_device_init_c() bind(c, name='gfr_device_init_c')
     end subroutine gfr_device_init_c
     subroutine gfr_device_phys_to_dyn_c(t, uv, q, nelem, npg, nlev, nq) &
          bind(c, name='gfr_device_phys_to_dyn_c')
       use iso_c_binding, only: c_double, c_int
       real(c_double) :: t(*), uv(*), q(*)                          ! device addresses
       integer(c_int), value :: nelem, npg, nlev, nq
     end subroutine gfr_device_phys_to_dyn_c
     subroutine gfr_device_compare_forcing_c(ft, fm, fq, nq, out) &
          bind(c, name='gfr_device_compare_forcing_c')
       use iso_c_binding, only: c_double, c_int
       real(c_double) :: ft(*), fm(*), fq(*), out(6)                ! host arrays
       integer(c_int), value :: nq
     end subroutine gfr_device_compare_forcing_c
     subroutine gfr_device_forcing_to_f90_c(ft, fm, fq, nq) bind(c, name='gfr_device_forcing_to_f90_c')
       use iso_c_binding, only: c_double, c_int
       real(c_double) :: ft(*), fm(*), fq(*)                        ! host arrays
       integer(c_int), value :: nq
     end subroutine gfr_device_forcing_to_f90_c
     subroutine gfr_device_dyn_to_phys_c(ps, phis, t, om, uv, q, nelem, npg, nlev, nq) &
          bind(c, name='gfr_device_dyn_to_phys_c')
       use iso_c_binding, only: c_double, c_int
       real(c_double) :: ps(*), phis(*), t(*), om(*), uv(*), q(*)   ! device addresses
       integer(c_int), value :: nelem, npg, nlev, nq
     end subroutine gfr_device_dyn_to_phys_c
  end interface
#endif

contains

  subroutine gfr_device_check_init()
#if defined(SE_GFR_DEVICE)
    use gllfvremap_mod,  only: gfr_init_hxx
    use cam_logfile,     only: iulog
    use spmd_utils,      only: masterproc
    call gfr_device_init_c()
    call gfr_init_hxx()
    if (masterproc) write(iulog,*) 'gfr_device_check: Hommexx GllFvRemap set up for the device shadow check'
#endif
  end subroutine gfr_device_check_init

  subroutine gfr_device_check(cols)
#if defined(SE_GFR_DEVICE)
    use dimensions_mod_cam, only: nelemd, nlev, fv_nphys, qsize
    use phys_grid,          only: get_dyn_col_p
    use spmd_utils,         only: masterproc, mpicom, mpi_max, mpi_real8
    use cam_logfile,        only: iulog
#endif
    type(phys_cols_t), intent(in) :: cols
#if defined(SE_GFR_DEVICE)
    integer  :: col_ind, ie, k, m, ncol, ierr, blk(1), nqc
    real(r8) :: dloc(2*7), dglob(2*7)

    ncol = cols%ncol
    call device_dyn_to_phys(ncol)

    ! 3. Compare with the Fortran remap (already in cols): max |diff| and max relative diff
    dloc = 0._r8
    nqc = min(nq, cols%ncnst)
    dloc(1) = maxval(abs(g1(:,1) - cols%ps));             dloc(8)  = maxval(abs(g1(:,1) - cols%ps)  /max(abs(cols%ps),  tiny(1._r8)))
    dloc(2) = maxval(abs(g1(:,2) - cols%phis));           dloc(9)  = maxval(abs(g1(:,2) - cols%phis)/max(abs(cols%phis),tiny(1._r8)))
    dloc(3) = maxval(abs(g2(:,:,1) - cols%t));            dloc(10) = maxval(abs(g2(:,:,1) - cols%t)    /max(abs(cols%t),    tiny(1._r8)))
    dloc(4) = maxval(abs(g2(:,:,2) - cols%u));            dloc(11) = maxval(abs(g2(:,:,2) - cols%u)    /max(abs(cols%u),    tiny(1._r8)))
    dloc(5) = maxval(abs(g2(:,:,3) - cols%v));            dloc(12) = maxval(abs(g2(:,:,3) - cols%v)    /max(abs(cols%v),    tiny(1._r8)))
    dloc(6) = maxval(abs(g2(:,:,4) - cols%omega));        dloc(13) = maxval(abs(g2(:,:,4) - cols%omega)/max(abs(cols%omega),tiny(1._r8)))
    dloc(7) = maxval(abs(gq(:,:,1:nqc) - cols%q(:,:,1:nqc)))
    dloc(14) = maxval(abs(gq(:,:,1:nqc) - cols%q(:,:,1:nqc))/max(abs(cols%q(:,:,1:nqc)), tiny(1._r8)))
    call mpi_allreduce(dloc, dglob, size(dloc), mpi_real8, mpi_max, mpicom, ierr)
    ncall = ncall + 1
    if (masterproc) then
       write(iulog,'(a,i4,a)') 'gfr_device_check call ', ncall, ' (device C++ remap vs Fortran remap; max abs / max rel diff):'
       write(iulog,'(a,7(1x,a,es9.2,"/",es9.2))') '  ', &
            'ps',  dglob(1), dglob(8),  'phis', dglob(2), dglob(9),  'T', dglob(3), dglob(10), &
            'u',   dglob(4), dglob(11), 'v',    dglob(5), dglob(12), 'omega', dglob(6), dglob(13), &
            'q',   dglob(7), dglob(14)
    end if

#endif
  end subroutine gfr_device_check

  subroutine gfr_device_fill_cols(cols)
    ! -DSE_GFR_DEVICE_COUPLING: fill phys's dynamics->physics arrays from the
    ! device remap (Hommexx GllFvRemap) instead of the Fortran remap.
#if defined(SE_GFR_DEVICE)
    use cam_abortutils, only: endrun
#endif
    type(phys_cols_t), intent(inout) :: cols
#if defined(SE_GFR_DEVICE)
    integer :: nqc
    call device_dyn_to_phys(cols%ncol)
    nqc = min(nq, cols%ncnst)
    cols%ps    = g1(:,1)
    cols%phis  = g1(:,2)
    cols%t     = g2(:,:,1)
    cols%u     = g2(:,:,2)
    cols%v     = g2(:,:,3)
    cols%omega = g2(:,:,4)
    cols%q(:,:,1:nqc) = gq(:,:,1:nqc)
    if (nqc < cols%ncnst) cols%q(:,:,nqc+1:) = 0._r8
#endif
  end subroutine gfr_device_fill_cols

  subroutine gfr_device_forcing_from_cols(cols)
    ! -DSE_GFR_DEVICE_COUPLING: physics->dynamics through the device remap.
    ! Hommexx's forcing views are filled on the device; then copied into the
    ! Fortran forcing arrays, so CAM mode's push_forcing_to_c hands Hommexx
    ! the same values (no HOMME change needed for this step).
#if defined(SE_GFR_DEVICE)
    use dimensions_mod_cam, only: nelemd, nlev
    use element_state,      only: elem_derived_FT, elem_derived_FM, elem_derived_FQ
#endif
    type(phys_cols_t), intent(in) :: cols
#if defined(SE_GFR_DEVICE)
    call device_phys_to_dyn(cols)
    call gfr_device_forcing_to_f90_c(elem_derived_FT, elem_derived_FM, elem_derived_FQ, nq)
#endif
  end subroutine gfr_device_forcing_from_cols

  subroutine gfr_device_forcing_on_device(cols)
    ! -DSE_GFR_DEVICE_RESIDENT: physics->dynamics through the device remap only;
    ! Hommexx's forcing views are the forcing (HOMME's push is switched off).
    type(phys_cols_t), intent(in) :: cols
#if defined(SE_GFR_DEVICE)
    call device_phys_to_dyn(cols)
#endif
  end subroutine gfr_device_forcing_on_device

  subroutine gfr_device_check_return(cols)
    ! Return path: after the Fortran physics->dynamics remap and its DSS have
    ! filled elem_derived_FT/FM/FQ, scatter the physics tendencies (phys) into
    ! GllFvRemap's layout on the device, run Hommexx's remap + DSS into its
    ! forcing views, and compare those with the Fortran forcing.  Hommexx's
    ! forcing is overwritten by push_forcing_to_c before it is used, so the
    ! model's answers do not change.
#if defined(SE_GFR_DEVICE)
    use dimensions_mod_cam, only: nelemd, nlev
    use element_state,      only: elem_derived_FT, elem_derived_FM, elem_derived_FQ
    use spmd_utils,         only: masterproc, mpicom, mpi_max, mpi_real8
    use cam_logfile,        only: iulog
#endif
    type(phys_cols_t), intent(in) :: cols
#if defined(SE_GFR_DEVICE)
    integer  :: col_ind, k, m, ncol, ierr
    real(r8) :: dloc(6), dglob(6)

    if (.not. allocated(t_c)) return       ! dynamics->physics check not set up yet
    call device_phys_to_dyn(cols)
    call gfr_device_compare_forcing_c(elem_derived_FT, elem_derived_FM, elem_derived_FQ, nq, dloc)
    call mpi_allreduce(dloc, dglob, size(dloc), mpi_real8, mpi_max, mpicom, ierr)
    ncall_r = ncall_r + 1
    if (masterproc) then
       write(iulog,'(a,i4,a)') 'gfr_device_check_return call ', ncall_r, &
            ' (device C++ phys->dyn remap+DSS vs Fortran; max abs / max abs normalized by max|Fortran|):'
       write(iulog,'(a,3(1x,a,es9.2,"/",es9.2))') '  ', &
            'FT', dglob(1), dglob(4), 'FM', dglob(2), dglob(5), 'FQ', dglob(3), dglob(6)
    end if
#endif
  end subroutine gfr_device_check_return

#if defined(SE_GFR_DEVICE)
  subroutine device_dyn_to_phys(ncol)
    ! The remap on the device (C++ GllFvRemap) into this module's device arrays,
    ! then a device gather into physics-column order (the phys_cols_t layout);
    ! results copied to the host arrays g1, g2, gq.
    use dimensions_mod_cam, only: nelemd, nlev, fv_nphys, qsize
    use phys_grid,          only: get_dyn_col_p
    integer, intent(in) :: ncol
    integer :: col_ind, ie, k, m, blk(1)
    if (.not. allocated(t_c)) then
       npg = fv_nphys*fv_nphys
       nq  = qsize
       allocate(ps_c(npg,nelemd), phis_c(npg,nelemd), t_c(nlev,npg,nelemd), om_c(nlev,npg,nelemd))
       allocate(uv_c(nlev,2,npg,nelemd), q_c(nlev,nq,npg,nelemd))
       allocate(g2(ncol,nlev,4), g1(ncol,2), gq(ncol,nlev,nq), map_ie(ncol), map_col(ncol))
       do col_ind = 1, ncol
          call get_dyn_col_p(col_ind, ie, blk)
          map_ie(col_ind)  = ie
          map_col(col_ind) = blk(1)
       end do
       !$acc enter data create(ps_c, phis_c, t_c, om_c, uv_c, q_c, g2, g1, gq) copyin(map_ie, map_col)
    end if

    ! 1. The remap on the device (C++ GllFvRemap), into this module's device arrays
    call run_remap(ps_c, phis_c, t_c, om_c, uv_c, q_c, nelemd, npg, nlev, nq)

    ! 2. Gather on the device: GllFvRemap layout (lev fastest, by element) ->
    !    physics-column order, columns fastest (the phys_cols_t layout)
    !$acc parallel loop collapse(2) present(t_c, om_c, uv_c, ps_c, phis_c, g2, g1, map_ie, map_col)
    do k = 1, nlev
       do col_ind = 1, ncol
          g2(col_ind,k,1) = t_c (k,   map_col(col_ind), map_ie(col_ind))
          g2(col_ind,k,2) = uv_c(k,1, map_col(col_ind), map_ie(col_ind))
          g2(col_ind,k,3) = uv_c(k,2, map_col(col_ind), map_ie(col_ind))
          g2(col_ind,k,4) = om_c(k,   map_col(col_ind), map_ie(col_ind))
          if (k == 1) then
             g1(col_ind,1) = ps_c  (map_col(col_ind), map_ie(col_ind))
             g1(col_ind,2) = phis_c(map_col(col_ind), map_ie(col_ind))
          end if
       end do
    end do
    !$acc parallel loop collapse(3) present(q_c, gq, map_ie, map_col)
    do m = 1, nq
       do k = 1, nlev
          do col_ind = 1, ncol
             gq(col_ind,k,m) = q_c(k, m, map_col(col_ind), map_ie(col_ind))
          end do
       end do
    end do
    !$acc update self(g2, g1, gq)
  end subroutine device_dyn_to_phys

  subroutine device_phys_to_dyn(cols)
    ! Copy phys's return arrays to the device, scatter them on the device into
    ! GllFvRemap's layout, and run the C++ remap + DSS into Hommexx's forcing.
    use dimensions_mod_cam, only: nelemd, nlev
    use cam_abortutils,     only: endrun
    type(phys_cols_t), intent(in) :: cols
    integer :: col_ind, k, m, ncol
    if (.not. allocated(t_c)) call endrun('gfr_device: dynamics->physics device remap not set up')
    ncol = cols%ncol
    if (.not. allocated(rt_c)) then
       allocate(rt_c(nlev,npg,nelemd), ruv_c(nlev,2,npg,nelemd), rq_c(nlev,nq,npg,nelemd))
       allocate(r2(ncol,nlev,3), rq(ncol,nlev,nq))
       !$acc enter data create(rt_c, ruv_c, rq_c, r2, rq)
    end if
    r2(:,:,1) = cols%dtdt
    r2(:,:,2) = cols%dudt
    r2(:,:,3) = cols%dvdt
    rq(:,:,1:nq) = cols%q_out(:,:,1:nq)
    !$acc update device(r2, rq)

    ! Scatter on the device: physics-column order -> GllFvRemap layout
    !$acc parallel loop collapse(2) present(r2, rt_c, ruv_c, map_ie, map_col)
    do k = 1, nlev
       do col_ind = 1, ncol
          rt_c (k,   map_col(col_ind), map_ie(col_ind)) = r2(col_ind,k,1)
          ruv_c(k,1, map_col(col_ind), map_ie(col_ind)) = r2(col_ind,k,2)
          ruv_c(k,2, map_col(col_ind), map_ie(col_ind)) = r2(col_ind,k,3)
       end do
    end do
    !$acc parallel loop collapse(3) present(rq, rq_c, map_ie, map_col)
    do m = 1, nq
       do k = 1, nlev
          do col_ind = 1, ncol
             rq_c(k, m, map_col(col_ind), map_ie(col_ind)) = rq(col_ind,k,m)
          end do
       end do
    end do

    call run_remap_return(rt_c, ruv_c, rq_c, nelemd, npg, nlev, nq)
  end subroutine device_phys_to_dyn

  subroutine run_remap(ps, phis, t, om, uv, q, nelem, npg_, nlev, nq_)
    ! Explicit-shape dummies: their device addresses go to the C++ remap
    integer,  intent(in)    :: nelem, npg_, nlev, nq_
    real(r8), intent(inout) :: ps(npg_,nelem), phis(npg_,nelem)
    real(r8), intent(inout) :: t(nlev,npg_,nelem), om(nlev,npg_,nelem)
    real(r8), intent(inout) :: uv(nlev,2,npg_,nelem), q(nlev,nq_,npg_,nelem)
    !$acc host_data use_device(ps, phis, t, om, uv, q)
    call gfr_device_dyn_to_phys_c(ps, phis, t, om, uv, q, nelem, npg_, nlev, nq_)
    !$acc end host_data
  end subroutine run_remap

  subroutine run_remap_return(t, uv, q, nelem, npg_, nlev, nq_)
    integer,  intent(in)    :: nelem, npg_, nlev, nq_
    real(r8), intent(inout) :: t(nlev,npg_,nelem), uv(nlev,2,npg_,nelem), q(nlev,nq_,npg_,nelem)
    !$acc host_data use_device(t, uv, q)
    call gfr_device_phys_to_dyn_c(t, uv, q, nelem, npg_, nlev, nq_)
    !$acc end host_data
  end subroutine run_remap_return
#endif

end module gfr_device_check_mod
