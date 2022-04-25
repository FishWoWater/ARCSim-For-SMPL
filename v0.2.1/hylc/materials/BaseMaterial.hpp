#ifndef _BASEMATERIAL_HPP_
#define _BASEMATERIAL_HPP_

#include "../../vectors.hpp"
#include "../MathematicaDefinitions.h"
#include <utility>

namespace hylc {
typedef Mat<6, 6> Mat6x6;
typedef Vec<6> Vec6;

// Base class for energy densities
// depending on the stretching and curvature tensors
class BaseMaterial {
public:
  // energy density value
  virtual double psi(const Vec6 &strain) = 0;
  // pair<Hessian,Gradient>
  virtual std::pair<Mat6x6, Vec6> psi_drv(const Vec6 &strain) = 0;
  // gradient
  virtual Vec6 psi_grad(const Vec6 &strain) = 0;

  double density=1;
};
} // namespace hylc

#endif /* end of include guard: _BASEMATERIAL_HPP_ */
