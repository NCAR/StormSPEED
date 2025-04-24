! This module contains additional CAM-specific settings beside the dycore/components/homme/src/share/parallel_mod.F90 module.
!
module parallel_mod_cam
  use parallel_mod

  implicit none
  public 

  public :: copy_par

  interface assignment ( = )
    module procedure copy_par
  end interface

contains

! ================================================
!   copy_par: copy constructor for parallel_t type
!
!
!   Overload assignment operator for parallel_t
! ================================================

  subroutine copy_par(par2,par1)
    type(parallel_t), intent(out) :: par2
    type(parallel_t), intent(in)  :: par1

    par2%rank       = par1%rank
    par2%root       = par1%root
    par2%nprocs     = par1%nprocs
    par2%comm       = par1%comm
    par2%masterproc = par1%masterproc
    par2%dynproc    = par1%dynproc

  end subroutine copy_par

end module parallel_mod_cam
