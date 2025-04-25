! This module contains additional CAM-specific settings beside the dycore/components/homme/src/share/interpolate_mod.F90 module.
!
module interpolate_mod_cam
   use interpolate_mod
   use element_mod,            only : element_t
   use kinds,                  only : real_kind
   use dimensions_mod_cam,     only : np
   use control_mod_cam,        only : cubed_sphere_map
   use coordinate_systems_mod, only : cartesian2d_t
   use cube_mod_cam,           only : dmap_cam 
   use string_utils,           only : int2str
   use cam_abortutils,         only : endrun

   implicit none
   private
   save

   interface interpolate_scalar
      module procedure interpolate_scalar2d
      module procedure interpolate_scalar3d
   end interface

   interface interpolate_vector
      module procedure interpolate_vector3d
   end interface

   public :: vec_latlon_to_contra

contains

   !===============================
   !(Nair) Bilinear interpolation for every GLL grid cell
   !===============================

   function interpol_bilinear(cart, f, xoy, imin, imax, fillvalue) result(fxy)
      integer, intent(in)               :: imin,imax
      type (cartesian2D_t), intent(in)  :: cart
      real (kind=real_kind), intent(in) :: f(imin:imax,imin:imax)
      real (kind=real_kind)             :: xoy(imin:imax)
      real (kind=real_kind)             :: fxy     ! value of f interpolated to (x,y)
      real (kind=real_kind), intent(in), optional :: fillvalue
      ! local variables

      real (kind=real_kind) :: p,q,xp,yp ,y4(4)
      integer        :: l,j,k, ii, jj, na,nb,nm

      xp = cart%x
      yp = cart%y

      ! Search index along "x"  (bisection method)

      na = imin
      nb = imax
      do
         if  ((nb-na) <=  1)  exit
         nm = (nb + na)/2
         if (xp  >  xoy(nm)) then
            na = nm
         else
            nb = nm
         endif
      enddo
      ii = na

      ! Search index along "y"

      na = imin
      nb = imax
      do
         if  ((nb-na) <=  1)  exit
         nm = (nb + na)/2
         if (yp  >  xoy(nm)) then
            na = nm
         else
            nb = nm
         endif
      enddo
      jj = na

      ! GLL cell containing (xp,yp)

      y4(1) = f(ii,jj)
      y4(2) = f(ii+1,jj)
      y4(3) = f(ii+1,jj+1)
      y4(4) = f(ii,jj+1)

      if(present(fillvalue)) then
         if (any(y4==fillvalue)) then
            fxy = fillvalue
            return
         endif
      endif

      p = (xp - xoy(ii))/(xoy(ii+1) - xoy(ii))
      q = (yp - xoy(jj))/(xoy(jj+1) - xoy(jj))

      fxy = (1.0_real_kind - p)*(1.0_real_kind - q)* y4(1) + p*(1.0_real_kind - q) * y4(2)   &
            + p*q* y4(3) + (1.0_real_kind - p)*q * y4(4)
   end function interpol_bilinear

   ! =======================================
   ! interpolate_scalar
   !
   ! Interpolate a scalar field given in an element (fld_cube) to the points in
   ! interpdata%interp_xy(i), i=1 .. interpdata%n_interp.
   !
   ! Note that it is possible the given element contains none of the interpolation points
   ! =======================================
   subroutine interpolate_scalar2d(interpdata,fld_cube,nsize,nhalo,fld, fillvalue)
      use dimensions_mod_cam, only : npsq, fv_nphys,nc
      integer,             intent(in) ::  nsize,nhalo
      real (kind=real_kind),      intent(in) ::  fld_cube(1-nhalo:nsize+nhalo,1-nhalo:nsize+nhalo) ! cube field
      real (kind=real_kind),      intent(out)::  fld(:)          ! field at new grid lat,lon coordinates
      type (interpdata_t), intent(in) ::  interpdata
      real (kind=real_kind),      intent(in), optional :: fillvalue
      ! Local variables
      type (interpolate_t), pointer  ::  interp          ! interpolation structure

      integer :: i,imin,imax,ne
      real (kind=real_kind):: xoy(1-nhalo:nsize+nhalo),dx
      type (cartesian2D_t) :: cart

      if (nsize==np.and.nhalo==0) then
         !
         ! GLL grid
         !
         interp => interp_p
         xoy = interp%glp(:)
         imin = 1
         imax = np
      else if (nhalo>0.and.(nsize==fv_nphys.or.nsize==nc)) then
         !
         ! finite-volume grid
         !
         if (itype.ne.1) then
         call endrun('itype must be 1 for latlon output from finite-volume (non-GLL) grids')
         end if
         imin = 1-nhalo
         imax = nsize+nhalo
         !
         ! create normalized coordinates
         !
         dx      = 2.0_real_kind/REAL(nsize,KIND=real_kind)
         do i=imin,imax
         xoy(i) = -1.0_real_kind+(i-0.5_real_kind)*dx
         end do
      else
         call endrun('interpolate_scalar2d: resolution not supported')
      endif

      ! Choice for Native (high-order) or Bilinear interpolations
      if (itype == 0) then
         do i=1,interpdata%n_interp
         fld(i)=interpolate_2d(interpdata%interp_xy(i),fld_cube,interp,nsize,fillvalue)
         end do
      else if (itype == 1) then
         do i=1,interpdata%n_interp
         fld(i)=interpol_bilinear(interpdata%interp_xy(i),fld_cube,xoy,imin,imax,fillvalue)
         end do
      else
         call endrun("interpolate_scalar2d: wrong interpolation type: "//int2str(itype))
      end if

   end subroutine interpolate_scalar2d

   subroutine interpolate_scalar3d(interpdata,fld_cube,nsize,nhalo,nlev,fld, fillvalue)
      use dimensions_mod_cam, only : npsq, fv_nphys,nc
      integer ,            intent(in)  ::  nsize, nhalo, nlev
      real (kind=real_kind),      intent(in)  ::  fld_cube(1-nhalo:nsize+nhalo,1-nhalo:nsize+nhalo,nlev) ! cube field
      real (kind=real_kind),      intent(out) ::  fld(:,:)          ! field at new grid lat,lon coordinates
      type (interpdata_t), intent(in)  ::  interpdata
      real (kind=real_kind), intent(in), optional :: fillvalue
      ! Local variables
      type (interpolate_t), pointer  ::  interp          ! interpolation structure

      integer :: ne

      integer :: i, k, imin, imax
      real (kind=real_kind) :: xoy(1-nhalo:nsize+nhalo),dx

      type (cartesian2D_t) :: cart

      if (nsize==np.and.nhalo==0) then
         !
         ! GLL grid
         !
         interp => interp_p
         xoy = interp%glp(:)
         imin = 1
         imax = np
      else if (nhalo>0.and.(nsize==fv_nphys.or.nsize==nc)) then
         !
         ! finite-volume grid
         !
         if (itype.ne.1) then
         call endrun('itype must be 1 for latlon output from finite-volume (non-GLL) grids')
         end if
         imin = 1-nhalo
         imax = nsize+nhalo
         !
         ! create normalized coordinates
         !
         dx      = 2.0_real_kind/REAL(nsize,KIND=real_kind)
         do i=imin,imax
         xoy(i) = -1.0_real_kind+(i-0.5_real_kind)*dx
         end do
      else
         call endrun('interpolate_scalar3d: resolution not supported')
      endif

      ! Choice for Native (high-order) or Bilinear interpolations
      if (itype == 0) then
         do k=1,nlev
         do i=1,interpdata%n_interp
            fld(i,k)=interpolate_2d(interpdata%interp_xy(i),fld_cube(:,:,k),interp,nsize,fillvalue)
         end do
         end do
      elseif (itype == 1) then
         do k=1,nlev
         do i=1,interpdata%n_interp
            fld(i,k)=interpol_bilinear(interpdata%interp_xy(i),fld_cube(:,:,k),xoy,imin,imax,fillvalue)
         end do
         end do
      else
         call endrun("interpolate_scalar3d: wrong interpolation type: "//int2str(itype))
      endif
   end subroutine interpolate_scalar3d

   ! =======================================
   ! interpolate_vector
   !
   ! Interpolate a vector field given in an element (fld_cube)
   ! to the points in interpdata%interp_xy(i), i=1 .. interpdata%n_interp.
   !
   ! input_coords = 0    fld_cube given in lat-lon
   ! input_coords = 1    fld_cube given in contravariant
   !
   ! Note that it is possible the given element contains none of the interpolation points
   ! =======================================
   subroutine interpolate_vector3d(interpdata,elem,fld_cube,npts,nlev,fld,input_coords, fillvalue)
      implicit none
      type (interpdata_t),intent(in)       ::  interpdata
      type (element_t), intent(in)         ::  elem
      integer, intent(in)                  ::  npts, nlev
      real (kind=real_kind), intent(in)    ::  fld_cube(npts,npts,2,nlev) ! vector field
      real (kind=real_kind), intent(out)   ::  fld(:,:,:)          ! field at new grid lat,lon coordinates
      real (kind=real_kind), intent(in),optional :: fillvalue
      integer, intent(in)                  ::  input_coords

      ! Local variables
      real (kind=real_kind)    ::  fld_contra(npts,npts,2,nlev) ! vector field
      type (interpolate_t), pointer  ::  interp          ! interpolation structure

      real (kind=real_kind)    ::  v1,v2
      real (kind=real_kind)    ::  D(2,2)   ! derivative of gnomonic mapping
      real (kind=real_kind)    ::  JJ(2,2), tmpD(2,2)   ! derivative of gnomonic mapping

      integer :: i,j,k

      type (cartesian2D_t) :: cart
      if(present(fillvalue)) then
         if (any(fld_cube==fillvalue)) then
            fld = fillvalue
            return
         end if
      end if
      if (input_coords==0 ) then
         ! convert to contra
         do k=1,nlev
            do j=1,npts
               do i=1,npts
                  ! latlon->contra
                  fld_contra(i,j,1,k) = elem%Dinv(i,j,1,1)*fld_cube(i,j,1,k) + elem%Dinv(i,j,1,2)*fld_cube(i,j,2,k)
                  fld_contra(i,j,2,k) = elem%Dinv(i,j,2,1)*fld_cube(i,j,1,k) + elem%Dinv(i,j,2,2)*fld_cube(i,j,2,k)
               enddo
            enddo
         end do
      else
         fld_contra=fld_cube
      endif

      if (npts==np) then
         interp => interp_p
      else if (npts==np) then
         call endrun('interpolate_vector3d: Error in interpolate_vector(): input must be on velocity grid')
      endif

      ! Choice for Native (high-order) or Bilinear interpolations

      if (itype == 0) then
         do k=1,nlev
            do i=1,interpdata%n_interp
               fld(i,k,1)=interpolate_2d(interpdata%interp_xy(i),fld_contra(:,:,1,k),interp,npts)
               fld(i,k,2)=interpolate_2d(interpdata%interp_xy(i),fld_contra(:,:,2,k),interp,npts)
            end do
         end do
      elseif (itype == 1) then
         do k=1,nlev
            do i=1,interpdata%n_interp
               fld(i,k,1)=interpol_bilinear(interpdata%interp_xy(i),fld_contra(:,:,1,k),interp%glp(:),1,np)
               fld(i,k,2)=interpol_bilinear(interpdata%interp_xy(i),fld_contra(:,:,2,k),interp%glp(:),1,np)
            end do
         end do
      else
         call endrun("interpolate_vector3d: wrong interpolation type: "//int2str(itype))
      endif

      do i=1,interpdata%n_interp
         ! compute D(:,:) at the point elem%interp_cube(i)
         call dmap_cam(D,interpdata%interp_xy(i)%x,interpdata%interp_xy(i)%y,&
               elem%corners3D,cubed_sphere_map,elem%corners,elem%u2qmap,elem%facenum)
         do k=1,nlev
            ! convert fld from contra->latlon
            v1 = fld(i,k,1)
            v2 = fld(i,k,2)

            fld(i,k,1)=D(1,1)*v1 + D(1,2)*v2
            fld(i,k,2)=D(2,1)*v1 + D(2,2)*v2
         end do
      end do
   end subroutine interpolate_vector3d

   subroutine vec_latlon_to_contra(elem,nphys,nhalo,nlev,fld,fvm)
      use fvm_control_volume_mod, only : fvm_struct
      use dimensions_mod_cam,     only : fv_nphys
      integer      , intent(in)   :: nphys,nhalo,nlev
      real(kind=real_kind), intent(inout):: fld(1-nhalo:nphys+nhalo,1-nhalo:nphys+nhalo,2,nlev)
      type (element_t), intent(in)           :: elem
      type(fvm_struct), intent(in), optional :: fvm
      !
      ! local variables
      !
      integer :: i,j,k
      real(real_kind):: v1,v2

      if (nhalo==0.and.nphys==np) then
         do k=1,nlev
         do j=1,nphys
            do i=1,nphys
               ! latlon->contra
               v1 = fld(i,j,1,k)
               v2 = fld(i,j,2,k)
               fld(i,j,1,k) = elem%Dinv(i,j,1,1)*v1 + elem%Dinv(i,j,1,2)*v2
               fld(i,j,2,k) = elem%Dinv(i,j,2,1)*v1 + elem%Dinv(i,j,2,2)*v2
            enddo
         enddo
         end do
      else if (nphys==fv_nphys.and.nhalo.le.fv_nphys) then
         do k=1,nlev
         do j=1-nhalo,nphys+nhalo
            do i=1-nhalo,nphys+nhalo
               ! latlon->contra
               v1 = fld(i,j,1,k)
               v2 = fld(i,j,2,k)
               fld(i,j,1,k) = fvm%Dinv_physgrid(i,j,1,1)*v1 + fvm%Dinv_physgrid(i,j,1,2)*v2
               fld(i,j,2,k) = fvm%Dinv_physgrid(i,j,2,1)*v1 + fvm%Dinv_physgrid(i,j,2,2)*v2
            enddo
         enddo
         end do
      else
         call endrun('ERROR: vec_latlon_to_contra - grid not supported or halo too large')
      end if
   end subroutine vec_latlon_to_contra

end module interpolate_mod_cam
