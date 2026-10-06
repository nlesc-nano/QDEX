#include "basis_to_k.hpp"
#include <cmath>
#include <stdexcept>
#include <thread>
#include <Eigen/Dense>

#ifdef _OPENMP
  #include <omp.h>
#endif

namespace licpp {

/* =================== Parallel launcher ========================== */
template <typename Lambda>
inline void parallel_do_k(Lambda&& body, int nthreads) {
#ifdef _OPENMP
#pragma omp parallel num_threads(nthreads)
  {
    body(omp_get_thread_num());
  }
#else
  std::vector<std::thread> workers;
  for (int id = 1; id < nthreads; ++id)
    workers.emplace_back(body, id);
  body(0);
  for (auto& t : workers) t.join();
#endif
}

/* -----------------------------------------------------------------------
 * Real solid harmonics in Libint SO order (m = -l..l), up to l = 3, on the
 * unit sphere. Libint normalizes a contracted shell so that x^l exp(-a r^2)
 * has unit norm; the prefactors below make every m component unit-normalized
 * with those coefficients, as in libint's own integrals:
 *   l=2: sqrt3 xy, sqrt3 yz, (3z^2-r^2)/2, sqrt3 zx, sqrt3/2 (x^2-y^2)
 *   l=3: sqrt(5/8) y(3x^2-y^2), sqrt15 xyz, sqrt(3/8) y(5z^2-r^2),
 *        z(5z^2-3r^2)/2, sqrt(3/8) x(5z^2-r^2), sqrt15/2 z(x^2-y^2),
 *        sqrt(5/8) x(x^2-3y^2)
 * --------------------------------------------------------------------- */
constexpr double kSqrt3 = 1.7320508075688772;
constexpr double kSqrt15 = 3.872983346207417;
constexpr double kSqrt5_8 = 0.7905694150420949;
constexpr double kSqrt3_8 = 0.6123724356957945;

inline void rsh_array_l0(double, double, double, double* o) { o[0] = 1.0; }

inline void rsh_array_l1(double x,double y,double z,double* o) {
  o[0]=y;  o[1]=z;  o[2]=x;
}

inline void rsh_array_l2(double x,double y,double z,double* o) {
  const double r2 = 1.0;
  o[0]=kSqrt3*x*y;                // m=-2
  o[1]=kSqrt3*y*z;                // m=-1
  o[2]=0.5*(3.0*z*z-r2);          // m= 0
  o[3]=kSqrt3*z*x;                // m=+1
  o[4]=0.5*kSqrt3*(x*x-y*y);      // m=+2
}

inline void rsh_array_l3(double x,double y,double z,double* o) {
  o[0] = kSqrt5_8*y*(3.0*x*x - y*y);        // m = -3
  o[1] = kSqrt15*x*y*z;                     // m = -2
  o[2] = kSqrt3_8*y*(5.0*z*z - 1.0);        // m = -1   (since r^2=1 on unit sphere)
  o[3] = 0.5*z*(5.0*z*z - 3.0);             // m =  0
  o[4] = kSqrt3_8*x*(5.0*z*z - 1.0);        // m = +1
  o[5] = 0.5*kSqrt15*z*(x*x - y*y);         // m = +2
  o[6] = kSqrt5_8*x*(x*x - 3.0*y*y);        // m = +3
}

inline void rsh_array(int l,double x,double y,double z,double* o){
  switch(l){
    case 0: rsh_array_l0(x,y,z,o); break;
    case 1: rsh_array_l1(x,y,z,o); break;
    case 2: rsh_array_l2(x,y,z,o); break;
    case 3: rsh_array_l3(x,y,z,o); break;
    default:
      throw std::runtime_error("rsh_array: l>3 not implemented");
  }
}

/* ---- analytic FT of one Libint shell at a single k-point --------------- */
static std::vector<std::complex<double>>
shell_ft_complex(const libint2::Shell& sh,
                 const Eigen::Vector3d& k) {

  const double klen = k.norm();
  const int l = sh.contr[0].l;
  const size_t nfunc = sh.size();
  std::vector<std::complex<double>> vals(nfunc, std::complex<double>(0.0, 0.0));

  // (-i)^l overall phase for momentum-space solid harmonics
  const std::complex<double> il = std::pow(std::complex<double>(0.0, -1.0), l);

  if (klen < 1e-15) {
    double ang[9]; 
    rsh_array(l, 0.0, 0.0, 1.0, ang); 
    for (size_t p = 0; p < sh.nprim(); ++p) {
      const double alpha = sh.alpha[p];
      const double coeff = sh.contr[0].coeff[p];

      const std::complex<double> pref =
          il * (coeff)
             * std::pow(M_PI/alpha, 1.5)
             * std::pow(klen/(2.0*alpha), l)
             * std::exp(-klen*klen/(4.0*alpha));

      for (size_t f = 0; f < nfunc; ++f)
        vals[f] += pref * ang[f];
    }
    return vals;
  }

  const double x = k(0)/klen, y = k(1)/klen, z = k(2)/klen;
  double ang[9]; rsh_array(l, x, y, z, ang);

  for (size_t p = 0; p < sh.nprim(); ++p) {
    const double alpha = sh.alpha[p];
    const double coeff = sh.contr[0].coeff[p];

    const std::complex<double> pref =
        il * (coeff)
           * std::pow(M_PI/alpha, 1.5)
           * std::pow(klen/(2.0*alpha), l)
           * std::exp(-klen*klen/(4.0*alpha));

    for (size_t f = 0; f < nfunc; ++f)
      vals[f] += pref * ang[f];
  }

  // center phase e^{-i k·R}
  const Eigen::Vector3d R{sh.O[0], sh.O[1], sh.O[2]};
  const std::complex<double> phase = std::exp(std::complex<double>(0.0, -k.dot(R)));
  for (auto& v : vals) v *= phase;

  return vals;
}

/* AO-FT matrix: (n_ao × n_k) complex-valued */
Eigen::MatrixXcd ao_ft_complex(const std::vector<libint2::Shell>& shells,
                               const std::vector<KPoint>& kpts,
                               int nthreads) {
  const size_t n_ao = nbasis(shells);
  const size_t n_k  = kpts.size();
  Eigen::MatrixXcd F(n_ao, n_k); F.setZero();

  parallel_do_k([&](int tid){
    for (size_t kp = 0; kp < n_k; ++kp) {
      if (kp % nthreads != (size_t)tid) continue;
      Eigen::Vector3d kvec = kpts[kp].vec();
      size_t ao0 = 0;
      for (size_t s = 0; s < shells.size(); ++s) {
        auto vec = shell_ft_complex(shells[s], kvec);
        for (size_t f = 0; f < vec.size(); ++f)
          F(ao0 + f, kp) = vec[f];
        ao0 += vec.size();
      }
    }
  }, nthreads);
  return F;
}

/* ---- Concatenate AO values for all shells on many points ---------------- */
Eigen::MatrixXd ao_values_at_points(
    const std::vector<libint2::Shell>& shells,
    const std::vector<Eigen::Vector3d>& points,
    int nthreads) {

  const size_t npts = points.size();
  const size_t nao  = nbasis(shells);
  Eigen::MatrixXd A(nao, npts);
  A.setZero();

  Eigen::VectorXd Px(npts), Py(npts), Pz(npts);
  for (size_t i = 0; i < npts; ++i) {
    Px[i] = points[i][0];
    Py[i] = points[i][1];
    Pz[i] = points[i][2];
  }

  const auto s2bf = map_shell_to_basis_function(shells);

  #pragma omp parallel for num_threads(nthreads) schedule(dynamic)
  for (ptrdiff_t s = 0; s < (ptrdiff_t)shells.size(); ++s) {
    const auto& sh = shells[(size_t)s];
    const int l = sh.contr[0].l;
    const size_t nfunc = sh.size();
    const double Cx = sh.O[0], Cy = sh.O[1], Cz = sh.O[2];

    Eigen::ArrayXd x = Px.array() - Cx;
    Eigen::ArrayXd y = Py.array() - Cy;
    Eigen::ArrayXd z = Pz.array() - Cz;
    Eigen::ArrayXd xx = x*x, yy = y*y, zz = z*z;
    Eigen::ArrayXd r2 = xx + yy + zz;

    Eigen::ArrayXd radial = Eigen::ArrayXd::Zero(npts);
    for (size_t p = 0; p < sh.nprim(); ++p) {
      const double a = sh.alpha[p];
      const double c = sh.contr[0].coeff[p];
      radial += c * (-a * r2).exp();
    }

    // Same normalized solid harmonics as rsh_array (libint conventions)
    Eigen::MatrixXd block(nfunc, npts);
    if (l == 0) {
      block.row(0) = radial.matrix();
    } else if (l == 1) {
      block.row(0) = (radial * y).matrix();
      block.row(1) = (radial * z).matrix();
      block.row(2) = (radial * x).matrix();
    } else if (l == 2) {
      block.row(0) = (radial * (kSqrt3*x*y)).matrix();
      block.row(1) = (radial * (kSqrt3*y*z)).matrix();
      block.row(2) = (radial * (0.5*(3.0*zz - r2))).matrix();
      block.row(3) = (radial * (kSqrt3*z*x)).matrix();
      block.row(4) = (radial * (0.5*kSqrt3*(xx - yy))).matrix();
    } else if (l == 3) {
      block.row(0) = (radial * (kSqrt5_8 * y * (3.0*xx - yy))).matrix();
      block.row(1) = (radial * (kSqrt15*x*y*z)).matrix();
      block.row(2) = (radial * (kSqrt3_8 * y * (5.0*zz - r2))).matrix();
      block.row(3) = (radial * (0.5 * z * (5.0*zz - 3.0*r2))).matrix();
      block.row(4) = (radial * (kSqrt3_8 * x * (5.0*zz - r2))).matrix();
      block.row(5) = (radial * (0.5*kSqrt15 * z * (xx - yy))).matrix();
      block.row(6) = (radial * (kSqrt5_8 * x * (xx - 3.0*yy))).matrix();
    } else {
      continue;
    }
    
    const size_t bf0 = s2bf[(size_t)s];
    A.middleRows((Eigen::Index)bf0, (Eigen::Index)nfunc) = block;
  }

  return A;
}

} // namespace licpp

