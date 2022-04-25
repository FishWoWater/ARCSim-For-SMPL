#include "hylc.hpp"

#include "../blockvectors.hpp"
#include "materials/SplineMaterial.hpp"
#include "seamforce.hpp"

namespace hylc {

bool hylc_enabled() { return config.enabled; }

// NOTE bad design - should store somewhere locally..
std::shared_ptr<BaseMaterial> get_material() { return config.material; }

double get_density() {
  if (config.material == nullptr)
    return 0;
  else
    return config.material->density *config.weight_mult;
}

Node *find_across(const Edge *edge, const Node *oppv1) {
  Node *found = nullptr;
  Vert *v     = edge_opp_vert(edge, 0);
  if (v)
    if (v->node == oppv1)  // wrong side
      v = edge_opp_vert(edge, 1);
  if (v)
    found = v->node;
  return found;
}

template <Space s>
void compute_local_information(const Face *face, Vec18 &xlocal,
                               bool &nn0_exists, bool &nn1_exists,
                               bool &nn2_exists) {
  xlocal         = Vec18(0);
  const Node *n0 = face->v[0]->node, *n1 = face->v[1]->node,
             *n2 = face->v[2]->node;
  Vec3 x0 = pos<s>(n0), x1 = pos<s>(n1), x2 = pos<s>(n2);
  for (int j = 0; j < 3; ++j) {
    xlocal[0 + j] = x0[j];
    xlocal[3 + j] = x1[j];
    xlocal[6 + j] = x2[j];
  }

  // potential neighbors
  Node *nn0 = nullptr, *nn1 = nullptr, *nn2 = nullptr;
  {
    Edge *e = get_edge(n0, n1);
    nn0     = find_across(e, n2);  // vtx opposite of n2 across edge n0n1
    if (nn0) {
      Vec3 x = pos<s>(nn0);
      for (int j = 0; j < 3; ++j) xlocal[9 + j] = x[j];
    }
  }
  {
    Edge *e = get_edge(n1, n2);
    nn1     = find_across(e, n0);
    if (nn1) {
      Vec3 x = pos<s>(nn1);
      for (int j = 0; j < 3; ++j) xlocal[12 + j] = x[j];
    }
  }
  {
    Edge *e = get_edge(n2, n0);
    nn2     = find_across(e, n1);
    if (nn2) {
      Vec3 x = pos<s>(nn2);
      for (int j = 0; j < 3; ++j) xlocal[15 + j] = x[j];
    }
  }

  nn0_exists = nn0 ? true : false;
  nn1_exists = nn1 ? true : false;
  nn2_exists = nn2 ? true : false;
}

template <Space s>
double hylc_local_energy(const Face *face) {
  if (config.material == nullptr)
    return 0;

  namespace mm = hylc::mathematica;
  // 0. compute local primitive values
  Vec18 xlocal;
  // double l0, l1, l2;
  // double theta_rest0, theta_rest1, theta_rest2;
  bool nn0_exists, nn1_exists, nn2_exists;
  compute_local_information<s>(face, xlocal, nn0_exists, nn1_exists,
                               nn2_exists);

  double A     = face->a;      // material space / reference config area
  Mat2x2 invDm = face->invDm;  // shape matrix

  // 1. mmcpp compute epsilon(x),kappa(x) as vec6
  Vec6 strain;
  strain = mm::strain(xlocal, invDm, Vec3(nn0_exists, nn1_exists, nn2_exists));
  strain[0] -= 1;
  strain[2] -= 1;
  for (int i = 3; i < 6; i++) strain(i) *= config.bend_scale;

  // 2. mmcpp compute psi(epskappa)
  double E = get_material()->psi(strain);
  // psi is energy density, constant over triangle
  E *= A;

  // SEAM FORCES
  // spring for each boundary/seam-edge in a face
  // (seams shared by two faces will get two forces!)
  if (config.seam_stiffness > 0)
    for (int i = 0; i < 3; i++) {
      Edge *e = face->adje[i];
      if (!is_seam_or_boundary(e))
        continue;
      Vec3 x0  = e->n[0]->x;
      Vec3 x1  = e->n[1]->x;
      double L = e->l;
      E += config.seam_stiffness * seamforce_val(x0, x1, L);
    }

  return E * config.stiffness_mult;

  // NOTE: probably not optimal recomputing normals within mmcpp of strains
  // same goes for dihedral angle theta
}

template <Space s>
std::pair<Mat18x18, Vec18> hylc_local_forces(/* const*/ Face *face) {
  if (config.material == nullptr)
    return std::make_pair(Mat18x18(), Vec18());

  namespace mm = hylc::mathematica;
  using namespace hylc::mathematica;
  // 0. compute local primitive values
  Vec18 xlocal;
  // double l0, l1, l2;
  // double theta_rest0, theta_rest1, theta_rest2;
  bool nn0_exists, nn1_exists, nn2_exists;
  compute_local_information<s>(face, xlocal, nn0_exists, nn1_exists,
                               nn2_exists);

  double A     = face->a;      // material space / reference config area
  Mat2x2 invDm = face->invDm;  // shape matrix

  // 1. mmcpp compute epsilon, kappa as vec6 strain, and simulatenous grad and
  // hess
  std::tuple<std::vector<Mat18x18>, Mat6x18, Vec6> strain_hgv;
  strain_hgv = mm::strain_valdrv(xlocal, invDm,
                                 Vec3(nn0_exists, nn1_exists, nn2_exists));

  Vec6 &strain = std::get<2>(strain_hgv);
  strain[0] -= 1;
  strain[2] -= 1;
  Mat6x18 &strain_grad               = std::get<1>(strain_hgv);
  std::vector<Mat18x18> &strain_hess = std::get<0>(strain_hgv);
  for (int i = 3; i < 6; i++) strain(i) *= config.bend_scale;

  // 2. mmcpp compute dpsidek(ek),d2psidekdek(ek) at the same time
  std::pair<Mat6x6, Vec6> psidrv = get_material()->psi_drv(strain);
  // 4. assemble local grad: dekdx.T dpsidek
  Vec18 g = transpose(strain_grad) * psidrv.second;
  // 5. assemble local hess: dpsidek.T d2ekdxdx + dekdx.T d2psidekdek dekdx
  // second part 18x6 * 6x6 * 6x18
  Mat18x18 H = transpose(strain_grad) * psidrv.first * strain_grad;
  // first part 1x6 * 6x18x18
  for (int i = 0; i < 6; ++i) H += psidrv.second[i] * strain_hess[i];

  // DEBUG color by strain
  if (debug.debug_color != 0) {
    std::vector<int> shifts   = {1, 0, 1, 0, 0};
    std::vector<double> norms = {1.0, 1.0, 1.0, 200.0, 200.0};
    int i                     = debug.debug_color - 1;
    double x;
    if (i < 3)
      x = (strain(i) - shifts[i]) / norms[i];
    else {
      double lam1 =
          (strain(3) + strain(5)) * 0.5 +
          std::sqrt((strain(3) - strain(5)) * (strain(3) - strain(5)) * 0.25 +
                    strain(4) * strain(4));
      double lam2 =
          (strain(3) + strain(5)) * 0.5 -
          std::sqrt((strain(3) - strain(5)) * (strain(3) - strain(5)) * 0.25 +
                    strain(4) * strain(4));
      double l1, l2;
      if (std::abs(lam1) > std::abs(lam2)) {
        l1 = lam1;
        l2 = lam2;
      } else {
        l1 = lam2;
        l2 = lam1;
      }
      x = std::vector<double>{l1, l2}[i - 3] / norms[3];
    }

    face->hylc_in_range(0) = x;
  }

  g = g * A;
  H = H * A;
  // SEAM FORCES
  // spring for each boundary/seam-edge in a face
  // (seams shared by two faces will get two forces!)
  if (config.seam_stiffness > 0)
    for (int i = 0; i < 3; i++) {
      Edge *e = face->adje[i];
      if (!is_seam_or_boundary(e))
        continue;
      Vec3 x0  = e->n[0]->x;
      Vec3 x1  = e->n[1]->x;
      double L = e->l;
      auto drv = seamforce_drv(x0, x1, L);

      // GET CORRECT INDEX, there might be a simpler prior structure
      auto getIndexInFace = [&](Node *n) {
        return n == face->v[0]->node ? 0 : (n == face->v[1]->node ? 1 : 2);
      };
      int i0 = getIndexInFace(e->n[0]);
      int i1 = getIndexInFace(e->n[1]);

      for (int j = 0; j < 3; j++) {
        g(i0 * 3 + j) += config.seam_stiffness * drv.second(j);
        g(i1 * 3 + j) += config.seam_stiffness * drv.second(3 + j);
        for (int k = 0; k < 3; k++) {
          H(i0 * 3 + j, i0 * 3 + k) += config.seam_stiffness * drv.first(j, k);
          H(i1 * 3 + j, i0 * 3 + k) +=
              config.seam_stiffness * drv.first(3 + j, k);
          H(i0 * 3 + j, i1 * 3 + k) +=
              config.seam_stiffness * drv.first(j, 3 + k);
          H(i1 * 3 + j, i1 * 3 + k) +=
              config.seam_stiffness * drv.first(3 + j, 3 + k);
        }
      }
    }

  // SVD clamp eigenvalues of local matrix

  typedef Eigen::Matrix<double, 18, 18> emat18;
  typedef Eigen::Matrix<double, 18, 1> evec18;
  emat18 eH;
  for (int i = 0; i < 18; ++i) {
    for (int j = 0; j < 18; ++j) {
      eH(i, j) = H(i, j);
    }
  }

  Eigen::SelfAdjointEigenSolver<emat18> slv;
  slv.compute(eH);

  double mineig = 1e-8;
  if (slv.eigenvalues()(0) < mineig) {
    evec18 D = slv.eigenvalues().cwiseMax(mineig);

    eH = slv.eigenvectors() * D.asDiagonal() * slv.eigenvectors().transpose();

    for (int i = 0; i < 18; ++i) {
      for (int j = 0; j < 18; ++j) {
        H(i, j) = eH(i, j);
      }
    }
  }

  return std::make_pair(H * config.stiffness_mult, g * config.stiffness_mult);
}

template <Space s>
Vec18 hylc_local_forces_nojac(const Face *face) {
  if (config.material == nullptr)
    return Vec18();

  namespace mm = hylc::mathematica;
  using namespace hylc::mathematica;
  // 0. compute local primitive values
  Vec18 xlocal;
  // double l0, l1, l2;
  // double theta_rest0, theta_rest1, theta_rest2;
  bool nn0_exists, nn1_exists, nn2_exists;
  compute_local_information<s>(face, xlocal, nn0_exists, nn1_exists,
                               nn2_exists);

  double A     = face->a;      // material space / reference config area
  Mat2x2 invDm = face->invDm;  // shape matrix

  // 1. mmcpp compute epsilon, kappa as vec6 ek, and simulatenous grad and hess
  std::tuple<Mat6x18, Vec6> strain_gv;
  strain_gv = mm::strain_valgrad(xlocal, invDm,
                                 Vec3(nn0_exists, nn1_exists, nn2_exists));

  Vec6 &strain = std::get<1>(strain_gv);
  strain[0] -= 1;
  strain[2] -= 1;
  Mat6x18 &strain_grad = std::get<0>(strain_gv);
  for (int i = 3; i < 6; i++) strain(i) *= config.bend_scale;

  Vec6 psi_grad = get_material()->psi_grad(strain);
  Vec18 g       = transpose(strain_grad) * psi_grad;

  g = g * A;

  // SEAM FORCES
  // spring for each boundary/seam-edge in a face
  // (seams shared by two faces will get two forces!)
  if (config.seam_stiffness > 0)
    for (int i = 0; i < 3; i++) {
      Edge *e = face->adje[i];
      if (!is_seam_or_boundary(e))
        continue;
      Vec3 x0  = e->n[0]->x;
      Vec3 x1  = e->n[1]->x;
      double L = e->l;
      Vec6 drv = seamforce_grad(x0, x1, L);

      // GET CORRECT INDEX, there might be a simpler prior structure
      auto getIndexInFace = [&](Node *n) {
        return n == face->v[0]->node ? 0 : (n == face->v[1]->node ? 1 : 2);
      };
      int i0 = getIndexInFace(e->n[0]);
      int i1 = getIndexInFace(e->n[1]);

      for (int j = 0; j < 3; j++) {
        g(i0 * 3 + j) += config.seam_stiffness * drv(j);
        g(i1 * 3 + j) += config.seam_stiffness * drv(3 + j);
      }
    }

  return g * config.stiffness_mult;
}

template <int m, int n>
Mat<3, 3> submat3(const Mat<m, n> &A, int i, int j) {
  Mat3x3 Asub;
  for (int k = 0; k < 3; k++)
    for (int l = 0; l < 3; l++) Asub(k, l) = A(i * 3 + k, j * 3 + l);
  return Asub;
}

template <int n>
Vec<3> subvec3(const Vec<n> &b, int i) {
  Vec3 bsub;
  for (int k = 0; k < 3; k++) bsub[k] = b[i * 3 + k];
  return bsub;
}

template <int m>
void add_submat(const Mat<m * 3, m * 3> &Asub, const Vec<m, int> &ix,
                SpMat<Mat3x3> &A) {
  for (int i = 0; i < m; i++)
    if (ix[i] >= 0)
      for (int j = 0; j < m; j++)
        if (ix[j] >= 0)
          A(ix[i], ix[j]) += submat3(Asub, i, j);
}

template <int m>
void add_subvec(const Vec<m * 3> &bsub, const Vec<m, int> &ix,
                std::vector<Vec3> &b) {
  for (int i = 0; i < m; i++)
    if (ix[i] >= 0)
      b[ix[i]] += subvec3(bsub, i);
}

template <int n>
Mat<1, n> rowmat(const Vec<n> &v) {
  Mat<1, n> A;
  for (int i = 0; i < n; i++) A(0, i) = v[i];
  return A;
}

template <int m, int n, int p, int q>
Mat<m * p, n * q> kronecker(const Mat<m, n> &A, const Mat<p, q> &B) {
  Mat<m * p, n * q> C;
  for (int i = 0; i < m; i++)
    for (int j = 0; j < n; j++)
      for (int k = 0; k < p; k++)
        for (int l = 0; l < q; l++) C(i * p + k, j * q + l) = A(i, j) * B(k, l);
  return C;
}

Vec<6, int> indices(const Node *n0, const Node *n1, const Node *n2,
                    const Node *nn0, const Node *nn1, const Node *nn2) {
  Vec<6, int> ix;
  ix[0] = n0->index;
  ix[1] = n1->index;
  ix[2] = n2->index;
  ix[3] = (nn0) ? nn0->index : -1;
  ix[4] = (nn1) ? nn1->index : -1;
  ix[5] = (nn2) ? nn2->index : -1;
  return ix;
}

Vec<3, int> indices(const Node *n0, const Node *n1, const Node *n2) {
  Vec<3, int> ix;
  ix[0] = n0->index;
  ix[1] = n1->index;
  ix[2] = n2->index;
  return ix;
}

}  // namespace hylc

template <Space s>
double hylc::hylc_internal_energy(const Cloth &cloth) {
  // local, parallel per triangle
  const Mesh &mesh = cloth.mesh;
  int n_triangles  = mesh.faces.size();
  std::vector<double> local_E(n_triangles);
#pragma omp parallel for
  for (int i = 0; i < n_triangles; ++i) {
    local_E[i] = hylc_local_energy<s>(mesh.faces[i]);
  }
  // global sequential / omp reduction
  double E = 0;
#pragma omp parallel for reduction(+ : E)
  for (size_t i = 0; i < local_E.size(); i++) {
    E += local_E[i];
  }

  return E;
}

template <Space s>
void hylc::hylc_add_internal_forces(const Cloth &cloth, SpMat<Mat3x3> &A,
                                    std::vector<Vec3> &b, double dt) {
  // local, parallel per triangle
  const Mesh &mesh = cloth.mesh;
  int n_triangles  = mesh.faces.size();
  std::vector<std::pair<Mat18x18, Vec18>> local_Hg(n_triangles);

#pragma omp parallel for
  for (int i = 0; i < n_triangles; ++i) {
    local_Hg[i] = hylc_local_forces<s>(mesh.faces[i]);
  }

  // global sequential assembly
  // sum force and hess with correct indexing
  for (int i = 0; i < n_triangles; ++i) {
    auto Hg = local_Hg[i];

    Mat18x18 &H = Hg.first;
    Vec18 &g    = Hg.second;

    const Face *face = mesh.faces[i];
    const Node *n0 = face->v[0]->node, *n1 = face->v[1]->node,
               *n2 = face->v[2]->node;

    // try finding neighbors
    const Node *nn0 = find_across(get_edge(n0, n1), n2),
               *nn1 = find_across(get_edge(n1, n2), n0),
               *nn2 = find_across(get_edge(n2, n0), n1);

    // get velocities (missing neighbors set to 0)
    Vec18 vs;
    for (int j = 0; j < 3; ++j) {
      vs[3 * 0 + j] = n0->v[j];
      vs[3 * 1 + j] = n1->v[j];
      vs[3 * 2 + j] = n2->v[j];
      vs[3 * 3 + j] = (nn0) ? nn0->v[j] : 0;
      vs[3 * 4 + j] = (nn1) ? nn1->v[j] : 0;
      vs[3 * 5 + j] = (nn2) ? nn2->v[j] : 0;
    }
    // get indices (missing neighbors set to -1)
    auto ixs = hylc::indices(n0, n1, n2, nn0, nn1, nn2);
    // add to global matrix (missing neighbors with ix -1 ignored)
    if (dt == 0) {
      hylc::add_submat(H, ixs, A);   // -J
      hylc::add_subvec(-g, ixs, b);  // F
    } else {
      double damping = cloth.materials[face->label]->damping;  // or EDGEbased ?
      //-dt * (dt + damping) * J
      hylc::add_submat(dt * (dt + damping) * H, ixs, A);
      // dt * (F + (dt + damping) * J * vs)
      hylc::add_subvec(-dt * (g + (dt + damping) * H * vs), ixs, b);
    }
  }
}

template <Space s>
void hylc::hylc_add_internal_forces(const Cloth &cloth, std::vector<Vec3> &b,
                                    double dt, double stretchdamping) {
  typedef Mat<9, 9> Mat9x9;
  typedef Vec<9> Vec9;

  // local, parallel per triangle
  const Mesh &mesh = cloth.mesh;
  int n_triangles  = mesh.faces.size();
  std::vector<Vec18> local_g(n_triangles);
  std::vector<Vec9> local_Jfake_vs(n_triangles);

#pragma omp parallel for
  for (int i = 0; i < n_triangles; ++i) {
    local_g[i] = hylc_local_forces_nojac<s>(mesh.faces[i]);
  }

#pragma omp parallel for
  for (int i = 0; i < n_triangles; ++i) {
    auto &face = mesh.faces[i];
    Mat3x2 F   = derivative(pos<s>(face->v[0]->node), pos<s>(face->v[1]->node),
                          pos<s>(face->v[2]->node), face);
    Mat2x2 G   = (F.t() * F - Mat2x2(1)) / 2.;
    // NOTE dependency on non-hylc material params removed!
    // Vec4 k = stretching_stiffness(G,
    // (*::materials)[face->label]->stretching); double weakening =
    // (*::materials)[face->label]->weakening; k *= 1 / (1 + weakening *
    // face->damage);
    Mat2x3 D = derivative(face);
    Vec3 du = D.row(0), dv = D.row(1);
    Mat<3, 9> Du   = kronecker(rowmat(du), Mat3x3(1)),
              Dv   = kronecker(rowmat(dv), Mat3x3(1));
    const Vec3 &xu = F.col(0), &xv = F.col(1);  // should equal Du*mat_to_vec(X)
    Vec9 fuu = Du.t() * xu, fvv = Dv.t() * xv,
         fuv      = (Du.t() * xv + Dv.t() * xu) / 2.;
    Mat9x9 hess_e = (outer(fuu, fuu) + std::max(G(0, 0), 0.) * Du.t() * Du) +
                    (outer(fvv, fvv) + std::max(G(1, 1), 0.) * Dv.t() * Dv) +
                    (outer(fuu, fvv) + std::max(G(0, 0), 0.) * Dv.t() * Dv +
                     outer(fvv, fuu) + std::max(G(1, 1), 0.) * Du.t() * Du) +
                    2. * (outer(fuv, fuv));

    const Node *n0 = face->v[0]->node, *n1 = face->v[1]->node,
               *n2    = face->v[2]->node;
    Vec9 vs           = mat_to_vec(Mat3x3(n0->v, n1->v, n2->v));
    local_Jfake_vs[i] = -face->a * hess_e * vs;
  }

  // global sequential assembly
  // sum force and hess with correct indexing
  for (int i = 0; i < n_triangles; ++i) {
    Vec18 &g = local_g[i];

    const Face *face = mesh.faces[i];
    const Node *n0 = face->v[0]->node, *n1 = face->v[1]->node,
               *n2 = face->v[2]->node;

    // try finding neighbors
    const Node *nn0 = find_across(get_edge(n0, n1), n2),
               *nn1 = find_across(get_edge(n1, n2), n0),
               *nn2 = find_across(get_edge(n2, n0), n1);

    // get indices (missing neighbors set to -1)
    auto ixs = hylc::indices(n0, n1, n2, nn0, nn1, nn2);
    // add to global matrix (missing neighbors with ix -1 ignored)
    if (dt == 0) {
      hylc::add_subvec(-g, ixs, b);  // F
    } else {
      hylc::add_subvec(-dt * g, ixs, b);  // dt * F
      // dt*(dt+damp)*J*vs, only with "fake" stretching J (from default arcsim)
      hylc::add_subvec(dt * ((dt + stretchdamping) * local_Jfake_vs[i]),
                       hylc::indices(n0, n1, n2), b);
    }
  }
}

template double hylc::hylc_internal_energy<WS>(const Cloth &);
template double hylc::hylc_internal_energy<PS>(const Cloth &);
template void hylc::hylc_add_internal_forces<WS>(const Cloth &, SpMat<Mat3x3> &,
                                                 std::vector<Vec3> &, double);
template void hylc::hylc_add_internal_forces<PS>(const Cloth &, SpMat<Mat3x3> &,
                                                 std::vector<Vec3> &, double);
template void hylc::hylc_add_internal_forces<WS>(const Cloth &,
                                                 std::vector<Vec3> &, double,
                                                 double);
template void hylc::hylc_add_internal_forces<PS>(const Cloth &,
                                                 std::vector<Vec3> &, double,
                                                 double);

#include <fstream>
void hylc::hylc_write_strains(const std::string &filename, const Cloth &cloth) {
  namespace mm = hylc::mathematica;

  const Mesh &mesh = cloth.mesh;
  int n_triangles  = mesh.faces.size();
  std::vector<Vec6> strains(n_triangles);
  std::vector<Vec2> uvs(n_triangles);
#pragma omp parallel for
  for (int i = 0; i < n_triangles; ++i) {
    Face *face = mesh.faces[i];
    // 0. compute local primitive values
    Vec18 xlocal;
    // double l0, l1, l2;
    // double theta_rest0, theta_rest1, theta_rest2;
    bool nn0_exists, nn1_exists, nn2_exists;
    compute_local_information<WS>(face, xlocal, nn0_exists, nn1_exists,
                                  nn2_exists);

    // double A     = face->a;      // material space / reference config area
    Mat2x2 invDm = face->invDm;  // shape matrix
    // Vec2 t0 = face->t0, t1 = face->t1, t2 = face->t2;

    // t0 *= (double)nn0_exists;
    // t1 *= (double)nn1_exists;
    // t2 *= (double)nn2_exists;

    // 1. mmcpp compute epsilon(x),kappa(x) as vec6
    strains[i] =
        mm::strain(xlocal, invDm, Vec3(nn0_exists, nn1_exists, nn2_exists));
    strains[i][0] -= 1;
    strains[i][2] -= 1;
    for (int j = 3; j < 6; j++) strains[i](j) *= config.bend_scale;

    // centroid in uv for filtering
    const Vert *v0 = face->v[0], *v1 = face->v[1], *v2 = face->v[2];
    uvs[i] = 1.0 / 3.0 * (v0->u + v1->u + v2->u);
  }

  std::ofstream myfile;
  myfile.open(filename);
  for (int i = 0; i < n_triangles; ++i) {
    for (int j = 0; j < 6; ++j) myfile << strains[i][j] << " ";
    for (int j = 0; j < 2; ++j) myfile << uvs[i][j] << " ";
    myfile << "\n";
  }
  myfile.close();
}