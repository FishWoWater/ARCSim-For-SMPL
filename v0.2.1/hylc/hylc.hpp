#ifndef _HYLC_H_
#define _HYLC_H_

#include <Eigen/Dense>
#include <Eigen/SVD>
#include <chrono> // DEBUG timing
#include <memory>
#include <utility>
#include <vector>

#include "../cloth.hpp"
#include "../geometry.hpp"
#include "../sparse.hpp"

#include "hylc_conf.hpp"
#include "strain/strain.hpp"

#include <iostream>
#include <string>

#define USE_SPARSE3

namespace hylc {

bool hylc_enabled();
std::shared_ptr<BaseMaterial> get_material();
double get_density();

template <Space s> double hylc_internal_energy(const Cloth &cloth);

template <Space s>
void hylc_add_internal_forces(const Cloth &cloth, SpMat<Mat3x3> &A,
                              std::vector<Vec3> &b, double dt);
template <Space s>
void hylc_add_internal_forces(const Cloth &cloth, std::vector<Vec3> &b,
                              double dt, double stretchdamping);

void hylc_write_strains(const std::string & filename, const Cloth &cloth);
} // namespace hylc
#endif /* end of include guard: _HYLC_H_ */
