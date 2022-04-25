/*
  Copyright �2013 The Regents of the University of California
  (Regents). All Rights Reserved. Permission to use, copy, modify, and
  distribute this software and its documentation for educational,
  research, and not-for-profit purposes, without fee and without a
  signed licensing agreement, is hereby granted, provided that the
  above copyright notice, this paragraph and the following two
  paragraphs appear in all copies, modifications, and
  distributions. Contact The Office of Technology Licensing, UC
  Berkeley, 2150 Shattuck Avenue, Suite 510, Berkeley, CA 94720-1620,
  (510) 643-7201, for commercial licensing opportunities.

  IN NO EVENT SHALL REGENTS BE LIABLE TO ANY PARTY FOR DIRECT,
  INDIRECT, SPECIAL, INCIDENTAL, OR CONSEQUENTIAL DAMAGES, INCLUDING
  LOST PROFITS, ARISING OUT OF THE USE OF THIS SOFTWARE AND ITS
  DOCUMENTATION, EVEN IF REGENTS HAS BEEN ADVISED OF THE POSSIBILITY
  OF SUCH DAMAGE.

  REGENTS SPECIFICALLY DISCLAIMS ANY WARRANTIES, INCLUDING, BUT NOT
  LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
  FOR A PARTICULAR PURPOSE. THE SOFTWARE AND ACCOMPANYING
  DOCUMENTATION, IF ANY, PROVIDED HEREUNDER IS PROVIDED "AS
  IS". REGENTS HAS NO OBLIGATION TO PROVIDE MAINTENANCE, SUPPORT,
  UPDATES, ENHANCEMENTS, OR MODIFICATIONS.
*/

#include "display.hpp"
#include "solvers.h"
#include "sparse_solver.hpp"
#include "timer.hpp"
#include <cstdlib>
#include <iostream>
#include <numeric>

#include <Eigen/Sparse>

using namespace std;
using namespace alglib;
using namespace Eigen;

namespace {

struct Comparator {
  Comparator(std::vector<int> const &vec_) : vec(&vec_) {}
  bool operator()(int a, int b) { return (*vec)[a] < (*vec)[b]; }
  std::vector<int> const *vec;
};

std::vector<int> sort_permutation(std::vector<int> const &vec) {
  std::vector<int> p(vec.size());
  for (int i = 0; i < p.size(); ++i)
    p[i] = i;
  Comparator comp(vec);
  std::sort(p.begin(), p.end(), comp);
  return p;
}

template <typename T>
std::vector<T> apply_permutation(std::vector<T> const &vec,
                                 std::vector<int> const &p) {
  std::vector<T> sorted_vec(p.size());
  for (int i = 0; i < p.size(); ++i)
    sorted_vec[i] = vec[p[i]];
  return sorted_vec;
}
} // namespace

vector<Node *> *debug_nodes = 0;

// ostream &operator<< (ostream &out, taucs_ccs_matrix *A) {
//    out << "n: " << A->n << endl;
//    out << "m: " << A->m << endl;
//    out << "flags: " << A->flags << endl;
//    out << "colptr: ";
//    for (int i = 0; i <= A->n; i++)
//        out << (i==0?"":", ") << A->colptr[i];
//    out << endl;
//    out << "rowind: ";
//    for (int j = 0; j <= A->colptr[A->n]; j++)
//        out << (j==0?"":", ") << A->rowind[j];
//    out << endl;
//    out << "values.d: ";
//    for (int j = 0; j <= A->colptr[A->n]; j++)
//        out << (j==0?"":", ") << A->values.d[j];
//    out << endl;
//    return out;
//}

Eigen::SparseMatrix<double> sparse_to_eigen(const SpMat<double> &As_) {

  SpMat<double> As = As_;
  for (int i = 0; i < As.n; ++i) {
    std::vector<int> p = sort_permutation(As.rows[i].indices);
    As.rows[i].indices = apply_permutation(As.rows[i].indices, p);
    As.rows[i].entries = apply_permutation(As.rows[i].entries, p);
  }

  // assumption: A is square and symmetric
  int n = As.n;
  int nnz = 0;
  for (int i = 0; i < n; i++) {
    for (int k = 0; k < (int)As.rows[i].indices.size(); k++) {
      int j = As.rows[i].indices[k];
      if (j < i)
        continue;
      nnz++;
    }
  }
  // taucs_ccs_matrix *At = taucs_ccs_create(n,n, nnz, TAUCS_DOUBLE |
  // TAUCS_SYMMETRIC | TAUCS_LOWER);
  SparseMatrix<double> At(n, n);
  At.resizeNonZeros(nnz);

  int pos = 0;
  for (int i = 0; i < n; i++) {
    At.outerIndexPtr()[i] = pos;
    for (int k = 0; k < (int)As.rows[i].indices.size(); k++) {
      int j = As.rows[i].indices[k];
      if (j < i)
        continue;
      At.innerIndexPtr()[pos] = j;
      At.valuePtr()[pos] = As.rows[i].entries[k];
      pos++;
    }
  }
  At.outerIndexPtr()[n] = pos;
  return At;
}

template <int m>
SparseMatrix<double> sparse_to_eigen(const SpMat<Mat<m, m>> &As_) {

  SpMat<Mat<m, m>> As = As_;
  for (int i = 0; i < As.n; ++i) {
    // auto p = sort_permutation(As.rows[i].indices, std::less<int>());
    std::vector<int> p = sort_permutation(As.rows[i].indices);
    As.rows[i].indices = apply_permutation(As.rows[i].indices, p);
    As.rows[i].entries = apply_permutation(As.rows[i].entries, p);
  }

  // assumption: A is square and symmetric
  int n = As.n;
  int nnz = 0;
  for (int i = 0; i < n; i++) {
    for (int jj = 0; jj < (int)As.rows[i].indices.size(); jj++) {
      int j = As.rows[i].indices[jj];
      if (j < i)
        continue;
      nnz += (j == i) ? m * (m + 1) / 2 : m * m;
    }
  }
  // taucs_ccs_matrix *At = taucs_ccs_create
  //    (n*m,n*m, nnz, TAUCS_DOUBLE | TAUCS_SYMMETRIC | TAUCS_LOWER);
  SparseMatrix<double> At(m * n, m * n);
  At.resizeNonZeros(nnz);
  int pos = 0;
  for (int i = 0; i < n; i++) {
    for (int k = 0; k < m; k++) {
      At.outerIndexPtr()[i * m + k] = pos;
      for (int jj = 0; jj < (int)As.rows[i].indices.size(); jj++) {
        int j = As.rows[i].indices[jj];
        if (j < i)
          continue;
        const Mat<m, m> &Aij = As.rows[i].entries[jj];
        for (int l = (i == j) ? k : 0; l < m; l++) {
          At.innerIndexPtr()[pos] = j * m + l;
          At.valuePtr()[pos] = Aij(k, l);
          pos++;
        }
      }
    }
  }
  At.outerIndexPtr()[n * m] = pos;
  return At;
}

vector<double> alglib_linear_solve(const SpMat<double> &A,
                                   const vector<double> &b) {
  const int n = b.size();
  real_2d_array M;
  real_1d_array x, c;
  M.setlength(n, n);
  x.setlength(n);
  c.setcontent(n, &b[0]);
  for (int i = 0; i < A.m; i++) {
    const SpVec<double> &row = A.rows[i];
    for (int j = 0; j < A.n; j++)
      M(i, j) = 0;
    for (size_t jj = 0; jj < row.indices.size(); jj++) {
      int j = row.indices[jj];
      M(i, j) = row.entries[jj];
    }
  }

  ae_int_t info;
  densesolverreport rep;
  rmatrixsolve(M, n, c, info, rep, x);

  vector<double> ret(n);
  for (int i = 0; i < n; i++)
    ret[i] = x[i];

  return ret;
}

template <int C>
vector<Vec<C>> alglib_linear_solve_vec(const SpMat<Mat<C, C>> &A,
                                       const vector<Vec<C>> &b) {
  const int n = b.size() * C;
  real_2d_array M;
  real_1d_array x, c;
  M.setlength(n, n);
  x.setlength(n);
  c.setlength(n);
  for (int i = 0; i < n; i++)
    c[i] = b[i / C][i % C];
  for (int i = 0; i < A.m; i++) {
    const SpVec<Mat<C, C>> &row = A.rows[i];
    for (size_t jj = 0; jj < row.indices.size(); jj++) {
      int j = row.indices[jj];
      for (int si = 0; si < C; si++)
        for (int sj = 0; sj < C; sj++)
          M(i * C + si, j * C + sj) = row.entries[jj](si, sj);
    }
  }

  ae_int_t info;
  densesolverreport rep;
  rmatrixsolve(M, n, c, info, rep, x);

  vector<Vec<C>> ret(n);
  for (int i = 0; i < n; i++)
    ret[i / C][i % C] = x[i];

  return ret;
}

#include <unsupported/Eigen/IterativeSolvers> // for minres
// NOTE: unified solver function
void eigen_linear_solve(SparseMatrix<double> &A, const Map<VectorXd const> &b,
                        Map<VectorXd> &x) {

  // SimplicialLDLT<SparseMatrix<double>, Lower> solver(A);
  // VectorXd D = solver.vectorD();
  // double lnegcheck = D.minCoeff();
  // double lmin = D.cwiseAbs().minCoeff();
  // double lmax = D.cwiseAbs().maxCoeff();
  // printf("lmin: %.2e, lmax: %.2e,\n cond: %.2e\n",lmin,lmax,lmax/lmin);
  // if (lnegcheck < 0)
  //   printf("  matrix is not positive definite!\n");
  // if (lnegcheck = 0)
  //   printf("  matrix is only positive semi-definite!\n");
  // printf("\n");
  // x = solver.solve(b);

  SimplicialLDLT<SparseMatrix<double>, Lower | Upper> solver(A);
  // VectorXd D = solver.vectorD();
  // double lnegcheck = D.minCoeff();
  // double lmin = D.cwiseAbs().minCoeff();
  // double lmax = D.cwiseAbs().maxCoeff();
  // printf("lmin: %.2e, lmax: %.2e,\n cond: %.2e\n",lmin,lmax,lmax/lmin);
  // if (lnegcheck < 0)
  //   printf("  matrix is not positive definite!\n");
  // if (lnegcheck == 0)
  //   printf("  matrix is only positive semi-definite!\n");

  // SparseLU< SparseMatrix<double> > solver(A);
  // SparseQR< SparseMatrix<double>, COLAMDOrdering<int> > solver(A);
  // BiCGSTAB<SparseMatrix<double>/*, IncompleteLUT<double>*/> solver(A);
  // LeastSquaresConjugateGradient<SparseMatrix<double>> solver(A);
  // solver.setMaxIterations(10000);

  if (solver.info() != Success) {
    printf("Eigen: Factorization failed\n");
  }

  x = solver.solve(b);

  if (solver.info() != Success) {
    printf("Eigen: Solve failed\n");
  }

  // printf("xTb=xTAx=%.2e,  res=%.2e\n", x.dot(b), (A * x - b).norm() / b.norm());
}

vector<double> eigen_linear_solve(const SpMat<double> &A,
                                  const vector<double> &b) {
  if (b.size() < 20)
    return alglib_linear_solve(A, b);
  SparseMatrix<double> Aeigen = sparse_to_eigen(A);
  Map<VectorXd const> b_(b.data(), b.size());

  vector<double> x(b.size());
  Map<VectorXd> x_(x.data(), x.size());

  eigen_linear_solve(Aeigen, b_, x_);

  return x;
}

template <int m>
vector<Vec<m>> eigen_linear_solve(const SpMat<Mat<m, m>> &A,
                                  const vector<Vec<m>> &b) {
  if (b.size() < 6)
    return alglib_linear_solve_vec(A, b);
  SparseMatrix<double> Aeigen = sparse_to_eigen(A);
  Map<VectorXd const> b_(&b[0][0], m * b.size());

  vector<Vec<m>> x(b.size());
  Map<VectorXd> x_(&x[0][0], m * x.size());

  eigen_linear_solve(Aeigen, b_, x_);

  return x;
}

template vector<Vec3> eigen_linear_solve(const SpMat<Mat3x3> &A,
                                         const vector<Vec3> &b);
