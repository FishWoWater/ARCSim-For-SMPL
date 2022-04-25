#ifndef _MATHEMATICADEFINITIONS_H_
#define _MATHEMATICADEFINITIONS_H_

#include <math.h>

namespace hylc {
namespace mathematica {

typedef double Real;
inline Real Power(Real A, Real B) { return std::pow(A, B); }
// inline Real Power(const Real &A, const int &B) { return std::pow(A, B); }
inline Real Power(Real A, int B) {
  assert(abs(B) < 8 && "require higher pow, implement using std.");
  Real out = 1.0;
  if (B >= 0) {
    for (int i = 0; i < B; ++i)
      out *= A;
  } else {
    Real invA = 1.0 / A;
    for (int i = 0; i < -B; ++i)
      out *= invA;
  }
  return out;
}
inline Real Abs(Real A) { return std::abs(A); }
inline Real Sqrt(Real A) { return std::sqrt(A); }
inline Real Exp(Real A) { return std::exp(A); }
inline Real Sin(Real A) { return std::sin(A); }
inline Real Cos(Real A) { return std::cos(A); }
inline Real Log(Real A) { return std::log(A); }
inline Real ArcCos(Real A) { return std::acos(A); }
inline Real ArcTan(Real X, const Real Y) { return std::atan2(Y, X); }
inline Real Tan(Real X) { return std::tan(X); }
inline Real Sec(Real X) { return 1.0/std::cos(X); }
inline Real RealAbs(Real X) { return std::abs(X); }
inline Real Sign(Real X) { return (X > 0) - (X < 0); }

} // namespace Mathematica
} // namespace HYLC

#endif /* end of include guard: _MATHEMATICADEFINITIONS_H_ */
