#ifndef __SEAMFORCE__HPP__
#define __SEAMFORCE__HPP__

#include "../vectors.hpp"
#include <utility>

namespace hylc {
double seamforce_val(const Vec<3> &x0, const Vec<3> &x1, double L);
Vec<6> seamforce_grad(const Vec<3> &x0, const Vec<3> &x1, double L);
std::pair<Mat<6, 6>, Vec<6>> seamforce_drv(const Vec<3> &x0, const Vec<3> &x1,
                                           double L);
}  // namespace hylc

#endif  // __SEAMFORCE__HPP__