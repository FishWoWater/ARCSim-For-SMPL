#ifndef _HSPLINES2_H_
#define _HSPLINES2_H_

#include <algorithm>
#include <vector>
#include "../../vectors.hpp"

struct HermiteSpline2D {
  int k0, k1;
  std::vector<double> tu, tv, p, mu, mv, muv;
  int ext = 0;

  typedef Mat<4, 4> Mat4x4;
  typedef Vec<4> Vec4;

  // M = np.array([
  //         [2,-2,1,1],
  //         [-3,3,-2,-1],
  //         [0,0,1,0],
  //         [1,0,0,0]
  //     ])
  static Mat4x4 M, MT, MextL, MextR, MextLT, MextRT;

  double eval(double x, double y, int dx = 0, int dy = 0) {
    if (tu.size() == 0 || tv.size() == 0)
      return 0;

    auto itu = std::upper_bound(tu.begin(), tu.end(),
                                x);  // index to the right of element
    auto itv = std::upper_bound(tv.begin(), tv.end(),
                                y);  // index to the right of element
    int i = itv - tv.begin();
    int j = itu - tu.begin();

    bool extrapol_u = false;
    bool extrapol_v = false;
    bool extrapol_u_right = false;
    bool extrapol_v_right = false;

    int ntv = (int)tv.size();
    int ntu = (int)tu.size();

    if (i <= 0 || i >= ntv) {
      extrapol_v = true;
      extrapol_v_right = i >= ntv;
      i          = std::max(1, std::min(ntv - 1, i));
    }
    if (j <= 0 || j >= ntu) {
      extrapol_u = true;
      extrapol_u_right = j >= ntu;
      j          = std::max(1, std::min(ntu - 1, j));
    }
    double deltau = (tu[j] - tu[j - 1]);
    double deltav = (tv[i] - tv[i - 1]);
    double u,v;
    if (extrapol_u_right)
     u      = (x - tu[j]) / deltau;
     else
     u      = (x - tu[j - 1]) / deltau;
    if (extrapol_v_right)
     v      = (y - tv[i]) / deltav;
    else
     v      = (y - tv[i - 1]) / deltav;

    i = i - 1; // shift ixs to be left of x,y
    j = j - 1;
    double p00, p01, p10, p11;
    double mu00, mu01, mu10, mu11;
    double mv00, mv01, mv10, mv11;
    double muv00, muv01, muv10, muv11;
    p00   = p[i * ntu + j];
    p01   = p[i * ntu + j + 1];
    p10   = p[(i + 1) * ntu + j];
    p11   = p[(i + 1) * ntu + j + 1];
    mu00  = mu[i * ntu + j] * deltau;
    mu01  = mu[i * ntu + j + 1] * deltau;
    mu10  = mu[(i + 1) * ntu + j] * deltau;
    mu11  = mu[(i + 1) * ntu + j + 1] * deltau;
    mv00  = mv[i * ntu + j] * deltav;
    mv01  = mv[i * ntu + j + 1] * deltav;
    mv10  = mv[(i + 1) * ntu + j] * deltav;
    mv11  = mv[(i + 1) * ntu + j + 1] * deltav;
    muv00 = muv[i * ntu + j] * deltav * deltau;
    muv01 = muv[i * ntu + j + 1] * deltav * deltau;
    muv10 = muv[(i + 1) * ntu + j] * deltav * deltau;
    muv11 = muv[(i + 1) * ntu + j + 1] * deltav * deltau;

    double uu = u * u, vv = v * v;
    Vec4 U, V;
    double deltas = 1.0;  // accumulate chain rule derivatives
    assert(dx >= 0 && dx < 3 && dy >= 0 && dy < 3);
    if (dx == 0)
      U = Vec4{uu * u, uu, u, 1.0};
    else if (dx == 1) {
      U = Vec4{3 * uu, 2 * u, 1.0, 0.0};
      deltas *= 1.0 / deltau;
    } else if (dx == 2) {
      U = Vec4{6 * u, 2, 0.0, 0.0};
      deltas *= 1.0 / (deltau * deltau);
    }
    if (dy == 0)
      V = Vec4{vv * v, vv, v, 1.0};
    else if (dy == 1) {
      V = Vec4{3 * vv, 2 * v, 1.0, 0.0};
      deltas *= 1.0 / deltav;
    } else if (dy == 2) {
      V = Vec4{6 * v, 2, 0.0, 0.0};
      deltas *= 1.0 / (deltav * deltav);
    }

    // get Mu and Mv as ptr based on ext, remove u3 u2
    // if const remove u1 0, if clamp set all U to 0

    Mat4x4 *Mu = &M, *MvT = &MT;
    // Mat4x4 *Mu = &M, *MvT = &M;

//     Mat4x4 Mu_dbg, MvT_dbg;
//     if (!(k0 == 0 && k1 == 2)) { // DEBUG
//       // linear else
//  Mu_dbg = Mat4x4 {Vec4{0,0,-1,1},Vec4{0,0,1,0},Vec4{0,0,0,0},Vec4{0,0,0,0}};
//  MvT_dbg = Mu_dbg.t();
//       Mu = &Mu_dbg;
//       MvT = &MvT_dbg;
//     }


    assert(ext >= 0 && ext <= 2);
    if (extrapol_u) {
      if (!extrapol_u_right)
        Mu = &MextL;
      else
        Mu = &MextR;
      U[0] = 0; // u^3 -> 0 // shouldnt matter because rows of M are 0 anyway
      U[1] = 0; // u^2 -> 0 // same
      if (ext > 0)
        U[2] = 0; // u -> 0
      if (ext > 1)
        U[3] = 0; // u^0 -> 0
    }
    if (extrapol_v) {
      if (!extrapol_v_right)
        // MvT = &MextL;
        MvT = &MextLT;
      else
        // MvT = &MextR;
        MvT = &MextRT;
      V[0]  = 0;
      V[1]  = 0;
      if (ext > 0)
        V[2] = 0;
      if (ext > 1)
        V[3] = 0;
    }

    // B = np.array([
    //     [p00,p10, mv00, mv10],
    //     [p01,p11, mv01, mv11],
    //     [mu00,mu10, muv00, muv10],
    //     [mu01,mu11, muv01, muv11]
    // ])
    // B = B * np.array([
    //     [1,1,deltav,deltav],
    //     [1,1,deltav,deltav],
    //     [deltau,deltau,deltau*deltav,deltau*deltav],
    //     [deltau,deltau,deltau*deltav,deltau*deltav],
    // ])
    // delta stuff already done above
    Mat4x4 B{Vec4{p00, p01, mu00, mu01}, Vec4{p10, p11, mu10, mu11},
             Vec4{mv00, mv01, muv00, muv01}, Vec4{mv10, mv11, muv10, muv11}};

    // derivatives are 1/deltau^n d^nU.T/du^n M B M.T d^mV/dv^m 1/deltav^m

    // # U.T M B M.T V
    // Y[k] = np.einsum("i,ij,jk,lk,l->",U,M,B,M,V)
    return deltas * dot(U, ((*Mu) * (B * ((*MvT) * V))));
    // return deltas * dot(U, ((*Mu) * (B * ((*MvT).t() * V))));
  }

  inline double dx(double x, double y) { return eval(x, y, 1, 0); }
  inline double dy(double x, double y) { return eval(x, y, 0, 1); }
  inline double dxdx(double x, double y) { return eval(x, y, 2, 0); }
  inline double dxdy(double x, double y) { return eval(x, y, 1, 1); }
  inline double dydy(double x, double y) { return eval(x, y, 0, 2); }
};

#endif  // _HSPLINES2_H_