
! This module contains additional CAM-specific settings beside the dycore/components/homme/src/share/cube_mod.F90 module.
!
module cube_mod_cam
   use cube_mod
   public :: dmap_cam

     ! ========================================================
  ! Dmap:
  !
  ! Initialize mapping that tranforms contravariant
  ! vector fields on the reference element onto vector fields on
  ! the sphere.
  ! ========================================================
  subroutine Dmap_cam(D, a,b, corners3D, ref_map, corners, u2qmap, facenum)
    real (kind=real_kind), intent(out)  :: D(2,2)
    real (kind=real_kind), intent(in)     :: a,b
    type (cartesian3D_t)   :: corners3D(4)  !x,y,z coords of element corners
    integer :: ref_map
    ! only needed for ref_map=0,1
    type (cartesian2D_t),optional   :: corners(4)    ! gnomonic coords of element corners
    real (kind=real_kind),optional  :: u2qmap(4,2)
    integer,optional  :: facenum



    if (ref_map==0) then
       if (.not. present ( corners ) ) &
            call abortmp('Dmap(): missing arguments for equiangular map')
       call dmap_equiangular_cam(D,a,b,corners,u2qmap,facenum)
    else if (ref_map==1) then
       call abortmp('equi-distance gnomonic map not yet implemented')
    else if (ref_map==2) then
       call dmap_elementlocal(D,a,b,corners3D)
    else
       call abortmp('bad value of ref_map')
    endif
  end subroutine Dmap_cam

  ! ========================================================
  ! Dmap:
  !
  ! Equiangular Gnomonic Projection
  ! Composition of equiangular Gnomonic projection to cubed-sphere face,
  ! followd by bilinear map to reference element
  ! ========================================================
  subroutine dmap_equiangular_cam(D, a,b, corners,u2qmap,facenum )
    use dimensions_mod, only : np
    real (kind=real_kind), intent(out)  :: D(2,2)
    real (kind=real_kind), intent(in)     :: a,b
    real (kind=real_kind)    :: u2qmap(4,2)
    type (cartesian2D_t)     :: corners(4)                          ! gnomonic coords of element corners
    integer :: facenum
    ! local
    real (kind=real_kind)  :: tmpD(2,2), Jp(2,2),x1,x2,pi,pj,qi,qj
    real (kind=real_kind), dimension(4,2) :: unif2quadmap

    ! input (a,b) shold be a point in the reference element [-1,1]
    ! compute Jp(a,b)
    Jp(1,1) = u2qmap(2,1) + u2qmap(4,1)*b
    Jp(1,2) = u2qmap(3,1) + u2qmap(4,1)*a
    Jp(2,1) = u2qmap(2,2) + u2qmap(4,2)*b
    Jp(2,2) = u2qmap(3,2) + u2qmap(4,2)*a

    ! map (a,b) to the [-pi/2,pi/2] equi angular cube face:  x1,x2
    ! a = gp%points(i)
    ! b = gp%points(j)
    pi = (1-a)/2
    pj = (1-b)/2
    qi = (1+a)/2
    qj = (1+b)/2
    x1 = pi*pj*corners(1)%x &
         + qi*pj*corners(2)%x &
         + qi*qj*corners(3)%x &
         + pi*qj*corners(4)%x
    x2 = pi*pj*corners(1)%y &
         + qi*pj*corners(2)%y &
         + qi*qj*corners(3)%y &
         + pi*qj*corners(4)%y

    call vmap(tmpD,x1,x2,facenum)

    ! Include map from element -> ref element in D
    D(1,1) = tmpD(1,1)*Jp(1,1) + tmpD(1,2)*Jp(2,1)
    D(1,2) = tmpD(1,1)*Jp(1,2) + tmpD(1,2)*Jp(2,2)
    D(2,1) = tmpD(2,1)*Jp(1,1) + tmpD(2,2)*Jp(2,1)
    D(2,2) = tmpD(2,1)*Jp(1,2) + tmpD(2,2)*Jp(2,2)
  end subroutine dmap_equiangular_cam

end module cube_mod_cam