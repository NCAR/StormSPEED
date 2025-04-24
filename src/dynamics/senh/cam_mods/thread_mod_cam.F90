! This module contains additional CAM-specific settings beside the dycore/components/homme/src/share/thread_mod.F90 module.
!
module thread_mod_cam
  use thread_mod
  use cam_logfile,            only: iulog
  use spmd_utils,             only: masterproc

  implicit none
  private

  integer, public,  TARGET :: max_num_threads ! maximum number of OpenMP threads
  integer, public          :: tracer_num_threads
  integer, public,  TARGET :: horz_num_threads , vert_num_threads

  integer, public, pointer :: NThreads   ! total number of threads
                                         ! standalone HOMME: from namelist
                                         ! in CAM: set by driver
  integer, public, pointer :: hthreads   ! computed based on nthreads, vthreads,nelemd
  integer, public, pointer :: vthreads   ! not used unless set in namelist

  public :: initomp

contains

#ifndef _OPENMP
  subroutine initomp
    max_num_threads = 1
    NThreads=>max_num_threads
    hthreads=>horz_num_threads
    vthreads => vert_num_threads
    if (masterproc) then
      write(iulog,*) "INITOMP: INFO: openmp not activated"
    end if
  end subroutine initomp
#else
  subroutine initomp
    !$OMP PARALLEL
    max_num_threads = omp_get_num_threads()
    !$OMP END PARALLEL
    NThreads=>max_num_threads
    hthreads=>horz_num_threads
    vthreads => vert_num_threads
    if (masterproc) then
      write(iulog,*) "INITOMP: INFO: number of OpenMP threads = ", max_num_threads
    end if
  end subroutine initomp
#endif

end module thread_mod
