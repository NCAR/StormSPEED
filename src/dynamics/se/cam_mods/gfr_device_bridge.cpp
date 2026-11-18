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
#include "ElementsForcing.hpp"
#include "Tracers.hpp"
#include "Types.hpp"

#include <cmath>
#include <algorithm>

namespace {
Homme::FunctorsBuffersManager* g_gfr_fbm = nullptr;
}

extern "C" {

// Create and set up GllFvRemap after prim_init2.  Call gfr_init_hxx (Fortran)
// afterwards to pass it the remap operators computed by gfr_init.
void gfr_device_init_c () {
  using namespace Homme;
  auto& c = Context::singleton();
  // CAM mode allocates the tracer forcing lazily (first push_forcing_to_c);
  // GllFvRemap keeps a copy of the Tracers views at setup and registers fq for
  // its DSS, so allocate it now, exactly as push_forcing_to_c would.
  auto& tracers = c.get<Tracers>();
  if (tracers.fq.data() == nullptr) {
    tracers.fq = decltype(tracers.fq)("fq", tracers.num_elems(), tracers.num_tracers());
  }
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

// Physics grid -> dynamics on the device.  T, uv: tendencies; q: constituents
// after physics.  Device addresses in GllFvRemap's layout (see above).  Writes
// Hommexx's forcing views (m_ft, m_fm, tracers.fq) and applies their DSS.
void gfr_device_phys_to_dyn_c (const double* T, const double* uv, const double* q,
                               const int nelem, const int npg, const int nlev, const int nq) {
  using namespace Homme;
  auto& c = Context::singleton();
  auto& gfr = c.get<GllFvRemap>();
  const int n0 = c.get<TimeLevel>().n0;
  const GllFvRemap::CPhys2T T_v(T, nelem, npg, nlev);
  const GllFvRemap::CPhys3T uv_v(uv, nelem, npg, 2, nlev), q_v(q, nelem, npg, nq, nlev);
  gfr.run_fv_phys_to_dyn(n0, T_v, uv_v, q_v);
  Kokkos::fence();
  gfr.run_fv_phys_to_dyn_dss();
  Kokkos::fence();
}

// Shadow check only: compare Hommexx's forcing views with the Fortran forcing
// arrays (theta-l_kokkos element_state: elem_derived_FT(np,np,nlev,nelemd),
// elem_derived_FM(np,np,3,nlev,nelemd), elem_derived_FQ(np,np,nlev,qsize_d,nelemd)),
// using the same index mapping as push_forcing_to_c.  out(1:6) = max|diff| and
// max|diff|/max|Fortran| for FT, FM (u and v), FQ (first nq tracers).
void gfr_device_compare_forcing_c (const double* ft_f, const double* fm_f, const double* fq_f,
                                   const int nq, double* out) {
  using namespace Homme;
  auto& forcing = Context::singleton().get<ElementsForcing>();
  auto& tracers = Context::singleton().get<Tracers>();
  const int ne = forcing.num_elems();
  HostViewUnmanaged<const Real*[NUM_PHYSICAL_LEV][NP][NP]>          ft90(ft_f, ne);
  HostViewUnmanaged<const Real*[NUM_PHYSICAL_LEV][3][NP][NP]>       fm90(fm_f, ne);
  HostViewUnmanaged<const Real*[QSIZE_D][NUM_PHYSICAL_LEV][NP][NP]> fq90(fq_f, ne);
  auto ft = Kokkos::create_mirror_view(forcing.m_ft); Kokkos::deep_copy(ft, forcing.m_ft);
  auto fm = Kokkos::create_mirror_view(forcing.m_fm); Kokkos::deep_copy(fm, forcing.m_fm);
  auto fq = Kokkos::create_mirror_view(tracers.fq);  Kokkos::deep_copy(fq, tracers.fq);
  const int nqc = std::min(nq, tracers.num_tracers());
  double d[3] = {0,0,0}, m[3] = {0,0,0};
  for (int ie = 0; ie < ne; ++ie)
    for (int k = 0; k < NUM_PHYSICAL_LEV; ++k) {
      const int ilev = k / VECTOR_SIZE, ivec = k % VECTOR_SIZE;
      for (int i = 0; i < NP; ++i)
        for (int j = 0; j < NP; ++j) {
          d[0] = std::max(d[0], std::abs(ft(ie,i,j,ilev)[ivec] - ft90(ie,k,i,j)));
          m[0] = std::max(m[0], std::abs(ft90(ie,k,i,j)));
          for (int dd = 0; dd < 2; ++dd) {
            d[1] = std::max(d[1], std::abs(fm(ie,dd,i,j,ilev)[ivec] - fm90(ie,k,dd,i,j)));
            m[1] = std::max(m[1], std::abs(fm90(ie,k,dd,i,j)));
          }
          for (int q = 0; q < nqc; ++q) {
            d[2] = std::max(d[2], std::abs(fq(ie,q,i,j,ilev)[ivec] - fq90(ie,q,k,i,j)));
            m[2] = std::max(m[2], std::abs(fq90(ie,q,k,i,j)));
          }
        }
    }
  for (int f = 0; f < 3; ++f) {
    out[f]   = d[f];
    out[3+f] = d[f] / std::max(m[f], 1e-300);
  }
}

} // extern "C"

#endif // KOKKOS_TARGET
