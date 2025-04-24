! This module contains additional CAM-specific settings beside the dycore/components/homme/src/share/dimensions_mod.F90 module.
!
module dimensions_mod_cam
  use dimensions_mod
  use shr_kind_mod, only : r8=>shr_kind_r8
#ifdef FVM_TRACERS
  use constituents, only : ntrac_d=>pcnst ! _EXTERNAL
#else
  use constituents, only : qsize_d=>pcnst ! _EXTERNAL
#endif

  implicit none

#ifdef FVM_TRACERS
  integer, parameter                      :: qsize_d =10 ! SE tracers (currently SE supports 10 condensate loading tracers)
#else
  integer, parameter                      :: ntrac_d = 0 ! No fvm tracers if CSLAM is off
#endif
  character(len=16),  allocatable, public :: cnst_name_gll(:)     ! constituent names for SE tracers
  character(len=128), allocatable, public :: cnst_longname_gll(:) ! long name of SE tracers
  logical                        , public :: lcp_moist = .true.

  integer, parameter             , public :: nc = 3               ! cslam resolution
  integer                        , public :: fv_nphys             ! physics-grid resolution - the "MAX" is so that the code compiles with NC=0
  integer                        , public :: ntrac = 0            ! ntrac is set in dyn_comp
  logical                        , public :: use_cslam = .false.  ! logical for CSLAM
  !
  ! fvm dimensions:
  logical, public                         :: lprint               ! for debugging
  integer, parameter,              public :: ngpc=3               ! number of Gausspoints for the fvm integral approximation   !phl change from 4
  integer, parameter,              public :: irecons_tracer=6     ! =1 is PCoM, =3 is PLM, =6 is PPM for tracer reconstruction
  integer,                         public :: irecons_tracer_lev(PLEV)
  integer, parameter,              public :: nhe=1                ! Max. Courant number
  integer, parameter,              public :: nhr=2                ! halo width needed for reconstruction - phl
  integer, parameter,              public :: nht=nhe+nhr          ! total halo width where reconstruction is needed (nht<=nc) - phl
  integer, parameter,              public :: ns=3                 ! quadratic halo interpolation - recommended setting for nc=3
  !nhc determines width of halo exchanged with neighboring elements
  integer, parameter,              public :: nhc = nhr+(nhe-1)+(ns-MOD(ns,2))/2 ! (different from halo needed for elements on edges and corners
  integer, parameter,              public :: lbc = 1-nhc
  integer, parameter,              public :: ubc = nc+nhc
  logical,                         public :: large_Courant_incr

  integer,                         public :: kmin_jet,kmax_jet    ! min and max level index for the jet
  integer,                         public :: fvm_supercycling    
  integer,                         public :: fvm_supercycling_jet

  integer, allocatable,            public :: kord_tr(:), kord_tr_cslam(:)
  
  real(r8),                        public :: nu_scale_top(PLEV)        ! scaling of del2 viscosity in sopnge layer (initialized in dyn_comp)
  real(r8),                        public :: nu_lev(PLEV)              ! level dependent del4 (u,v) damping
  real(r8),                        public :: nu_t_lev(PLEV)            ! level depedendet del4 T damping
  integer,                         public :: ksponge_end               ! sponge is active k=1,ksponge_end
  real(r8),                        public :: nu_div_lev(PLEV) = 1.0_r8 ! scaling of viscosity in sponge layer
                                                                       ! (set in prim_state; if applicable)
  real(r8),                        public :: kmvis_ref(PLEV)           ! reference profiles for molecular diffusion 
  real(r8),                        public :: kmcnd_ref(PLEV)           ! reference profiles for molecular diffusion  
  real(r8),                        public :: rho_ref(PLEV)             ! reference profiles for rho
  real(r8),                        public :: km_sponge_factor(PLEV)    ! scaling for molecular diffusion (when used as sponge)

  integer,                         public :: nhc_phys 
  integer,                         public :: nhe_phys 
  integer,                         public :: nhr_phys 
  integer,                         public :: ns_phys  

  integer,                         public :: npdg = 0  ! dg degree for hybrid cg/dg element  0=disabled 

  integer,                         public :: s_nv = 2*max_elements_attached_to_node

  integer,                         public :: nPhysProc ! This is the number of physics processors/ per dynamics processor

  public :: ntrac_d

end module dimensions_mod_cam