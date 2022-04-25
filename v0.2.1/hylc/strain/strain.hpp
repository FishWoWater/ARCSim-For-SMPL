#ifndef _STRAIN_HPP_
#define _STRAIN_HPP_
#define hylc_strain_II

#include "../../vectors.hpp"
#include "../MathematicaDefinitions.h"

#include <tuple>
#include <utility>
#include <vector>

namespace hylc {
typedef Mat<18, 18> Mat18x18;
typedef Vec<18> Vec18;
typedef Vec<6> Vec6;
typedef Mat<6, 18> Mat6x18;
namespace mathematica {
Vec6 strain(const Vec18 &xloc, const Mat2x2 &invDm, const Vec3& niexists);

std::tuple<std::vector<Mat18x18>, Mat6x18, Vec6>
strain_valdrv(const Vec18 &xloc, const Mat2x2 &invDm, const Vec3& niexists);

// for explicit integration without hessians
std::tuple<Mat6x18, Vec6>
strain_valgrad(const Vec18 &xloc, const Mat2x2 &invDm, const Vec3& niexists);
} // namespace mathematica
} // namespace hylc

#endif /* end of include guard: _STRAIN_HPP_ */
