// Bridge between the SE HOMME backend (Fortran) and Hommexx's GllFvRemap (C++),
// so the GLL <-> physics-grid remap can run on the device on arrays owned by
// the Fortran side (device addresses passed in from OpenACC).
//
// GllFvRemap is normally set up by EAMxx.  In CAM, Hommexx's functors are
// already initialized (prim_init2) before we get here, and their shared buffer
// cannot be re-allocated, so GllFvRemap gets its own FunctorsBuffersManager.
#ifdef KOKKOS_TARGET

#include "Context.hpp"
#include "GllFvRemap.hpp"
#include "FunctorsBuffersManager.hpp"
#include "SimulationParams.hpp"
#include "TimeLevel.hpp"

namespace {
Homme::FunctorsBuffersManager* g_gfr_fbm = nullptr;
}

extern "C" {

// Create and set up GllFvRemap after prim_init2.  Call gfr_init_hxx (Fortran)
// afterwards to pass it the remap operators computed by gfr_init.
void gfr_device_init_c () {
  using namespace Homme;
  auto& c = Context::singleton();
  auto& gfr = c.create_if_not_there<GllFvRemap>();
  gfr.setup();
  if (g_gfr_fbm == nullptr) {
    g_gfr_fbm = new FunctorsBuffersManager();
    g_gfr_fbm->request_size(gfr.requested_buffer_size());
    g_gfr_fbm->allocate();
  }
  gfr.init_buffers(*g_gfr_fbm);
  gfr.reset(c.get<SimulationParams>());
  gfr.init_boundary_exchanges();
}

// Dynamics -> physics grid on the device.  All pointers are DEVICE addresses of
// arrays in GllFvRemap's layout (LayoutRight, i.e. Fortran order reversed):
//   ps, phis          (ie, col)            Fortran: (npg, nelem)
//   T, omega          (ie, col, lev)       Fortran: (nlev, npg, nelem)
//   uv                (ie, col, 2, lev)    Fortran: (nlev, 2, npg, nelem)
//   q                 (ie, col, iq, lev)   Fortran: (nlev, nq, npg, nelem)
void gfr_device_dyn_to_phys_c (double* ps, double* phis, double* T, double* omega,
                               double* uv, double* q,
                               const int nelem, const int npg, const int nlev, const int nq) {
  using namespace Homme;
  auto& c = Context::singleton();
  auto& gfr = c.get<GllFvRemap>();
  const int n0 = c.get<TimeLevel>().n0;
  const GllFvRemap::Phys1T ps_v(ps, nelem, npg), phis_v(phis, nelem, npg);
  const GllFvRemap::Phys2T T_v(T, nelem, npg, nlev), omega_v(omega, nelem, npg, nlev);
  const GllFvRemap::Phys3T uv_v(uv, nelem, npg, 2, nlev), q_v(q, nelem, npg, nq, nlev);
  gfr.run_dyn_to_fv_phys(n0, ps_v, phis_v, T_v, omega_v, nullptr, nullptr, uv_v, q_v);
  Kokkos::fence();
}

} // extern "C"

#endif // KOKKOS_TARGET
