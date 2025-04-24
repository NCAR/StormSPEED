! This module contains additional CAM-specific settings beside the dycore/components/homme/src/share/hybvcoord_mod.F90 module.
!
module hybvcoord_mod_cam

   use hybvcoord_mod

   implicit none
   private

   !----------------------------------------------------------------------- 
   ! hvcoord_t: Hybrid level definitions: p = a*p0 + b*ps
   !            interfaces   p(k) = hyai(k)*ps0 + hybi(k)*ps
   !            midpoints    p(k) = hyam(k)*ps0 + hybm(k)*ps
   !-----------------------------------------------------------------------
   type, public :: hvcoord_t
   real(r8) ps0          ! base state surface-pressure for level definitions
   real(r8) hyai(plevp)  ! ps0 component of hybrid coordinate - interfaces
   real(r8) hyam(plev)   ! ps0 component of hybrid coordinate - midpoints
   real(r8) hybi(plevp)  ! ps  component of hybrid coordinate - interfaces
   real(r8) hybm(plev)   ! ps  component of hybrid coordinate - midpoints
   real(r8) etam(plev)   ! eta-levels at midpoints
   real(r8) etai(plevp)  ! eta-levels at interfaces
   real(r8) dp0(plev)      ! average layer thickness
   real(r8) hybd(plev)   ! difference in b (hybi) across layers
   real(r8) prsfac       ! log pressure extrapolation factor (time, space independent)
   end type

end module hybvcoord_mod_cam
