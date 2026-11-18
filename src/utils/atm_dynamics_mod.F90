module atm_dynamics_mod
   !----------------------------------------------------------------------------
   ! The dynamics process: the interface every dynamical core offers its host
   ! model for physics coupling and time stepping (after EAMxx's "atmosphere
   ! process").  A dycore provides a type that extends atm_dynamics_t and
   ! implements the deferred procedures; the host calls only these.
   !
   !   dyn_to_phys  fill this%phys from the dynamics state (on the physics columns)
   !   phys_to_dyn  read the physics tendencies from this%phys into the dynamics
   !   run          advance the dynamics one physics time step
   !
   ! The exchange goes only through this%phys (phys_cols_t: plain arrays plus
   ! labels for the vertical coordinate and wet/dry tracers), so the type has
   ! no host-model dependency.  Implemented today by se_dynamics_t (SE: CAM-SE
   ! and the E3SM HOMME dycores); FV, FV3 and MPAS can implement it the same way.
   !----------------------------------------------------------------------------
   use phys_cols_mod, only: phys_cols_t

   implicit none
   private

   type, abstract, public :: atm_dynamics_t
      ! Physics-column arrays exchanged with the host (shared by the dycore's
      ! import and export objects)
      type(phys_cols_t), pointer :: phys => null()
   contains
      procedure(dyn_process_iface), deferred :: dyn_to_phys
      procedure(dyn_process_iface), deferred :: phys_to_dyn
      procedure(dyn_process_iface), deferred :: run
   end type atm_dynamics_t

   abstract interface
      subroutine dyn_process_iface(this)
         import :: atm_dynamics_t
         class(atm_dynamics_t), intent(inout) :: this
      end subroutine dyn_process_iface
   end interface

end module atm_dynamics_mod
