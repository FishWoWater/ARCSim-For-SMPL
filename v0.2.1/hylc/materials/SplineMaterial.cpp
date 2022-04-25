#include "SplineMaterial.hpp"
#include <memory>
#include <vector>

using namespace hylc;
// using namespace fitpackpp;

// NOTE: various debugging definitions for playing with energy model
// #define BENDBARRIER
// #define BENDPOWER 1e0 // 1e3 too stiff, 1e2 also died with limit 400
// #define BENDLIMIT 200
// #define BENDPERCENT 0.25
// #define BENDCLAMPING
// #define BENDCLAMPINGPC
// #define BENDCLAMPLIMIT 230
// #define DOUBLECURVEPENALTY 1e-6
// #define NO1D

int sgn(double x) {
  return (x > 0) - (x < 0);
}


bool select_2D(int k0, int k1) {
  // return k0 == 0 && k1 == 2; // only poisson
  return true; // use all
  return false; // use none
}

SplineMaterial::SplineMaterial() {
  // s_kx = 1e-2;
  // s_ky = 1e-2;
  // s_kx = 1e1;
  // s_ky = 1e1;
  s_kx = 1.0;
  s_ky = 1.0;
}

double SplineMaterial::psi(const Vec6 &strain) {
  double val = 0;
  if (!initialized)
    return val;

  Vec3 sbend(0);

#ifndef BENDCLAMPING
  for (int i = 0; i < 3; i++)
    sbend[i] = strain[3+i];
#else
  for (int i = 0; i < 3; i++) {
    sbend[i] = strain[3+i];
    if (std::abs(sbend[i]) > BENDCLAMPLIMIT)
      sbend[i] = sgn(sbend[i]) * BENDCLAMPLIMIT;
  }
#endif //BENDCLAMPING

  // PC is a vector of (k1,k2,c^2),
  // where k1,k2 are principal curvatures and c^2 is the squared cos of the
  // angle between x axis and first eigenvector
  // double ss=0.25;
  // Vec3 PC_val = pc_val(ss*strain[3], ss*strain[4], ss*strain[5]);
  Vec3 PC_val = pc_val(sbend[0],sbend[1],sbend[2]);
  // Vec3 PC_val = pc_val(strain[3], strain[4], strain[5]);
  double l1   = PC_val(0);
  double l2   = PC_val(1);
  double cc   = PC_val(2);

#ifdef BENDCLAMPINGPC
  if (std::abs(l1) > BENDCLAMPLIMIT)
    l1 = sgn(l1) * BENDCLAMPLIMIT;
  if (std::abs(l2) > BENDCLAMPLIMIT)
    l2 = sgn(l2) * BENDCLAMPLIMIT;
#endif //BENDCLAMPINGPC

  // only normalize inplane components, bc bending normalization depends on
  // PC
  Vec3 S;
  for (int i = 0; i < 3; i++)
    S(i) = (strain(i)) / this->strainscale[i];

  // const
  val += C0;

  // 1D
  #ifndef NO1D
  for (auto &s : hsplines_1d) {
    if (s.k < 3) {  // in-plane
      double x = S(s.k);
      val += s.eval(x);
    } else if (s.k == 3) {  // x-bending
      val += s_kx * cc * s.eval(l1 / this->strainscale[s.k]) +
             (1 - cc) * s.eval(l2 / this->strainscale[s.k]);
    } else {  // y-bending
      val += s_ky * cc * s.eval(l2 / this->strainscale[s.k]) +
             (1 - cc) * s.eval(l1 / this->strainscale[s.k]);
    }
  }
  #endif

  // 2D
  for (auto &s : hsplines_2d) {
    if (!select_2D(s.k0,s.k1))
      continue;

    // assume sorted, k1 > k0
    assert(s.k1 > s.k0 && s.k0 < 3);
    if (s.k1 < 3) {  // in plane
      double x = S(s.k0);
      double y = S(s.k1);
      val += s.eval(x, y);
    } else {  // since we dont have double curve, k0 has to be in plane
      double x = S(s.k0);
      if (s.k1 == 3) {
        val += s_kx * (cc * s.eval(x, l1 / this->strainscale[s.k1]) +
                       (1 - cc) * s.eval(x, l2 / this->strainscale[s.k1]));
      } else {
        val += s_ky * (cc * s.eval(x, l2 / this->strainscale[s.k1]) +
                       (1 - cc) * s.eval(x, l1 / this->strainscale[s.k1]));
      }
    }
  }

#ifdef BENDBARRIER
  double mult = BENDPOWER; 
  double T = BENDLIMIT;
  if (std::abs(l1) > T) {
    double r = ((l1 - sgn(l1) * T) / (BENDPERCENT*T));
    val += mult * r * r;
  }
  if (std::abs(l2) > T) {
    double r = ((l2 - sgn(l2) * T) / (BENDPERCENT*T));
    val += mult * r * r;
  }
#endif  // BENDBARRIER


#ifdef DOUBLECURVEPENALTY
  double k = (sbend[0]*sbend[2] - sbend[1]*sbend[1]);
  val += DOUBLECURVEPENALTY*k*k;
#endif
  return val;
}

Vec6 SplineMaterial::psi_grad(const Vec6 &strain) {
  Vec6 grad(0);
  if (!initialized)
    return grad;

  Vec3 gradpc(0);
  Vec3 sbend(0);

#ifndef BENDCLAMPING
  for (int i = 0; i < 3; i++)
    sbend[i] = strain[3+i];
#else
  for (int i = 0; i < 3; i++) {
    sbend[i] = strain[3+i];
    if (std::abs(sbend[i]) > BENDCLAMPLIMIT)
      sbend[i] = sgn(sbend[i]) * BENDCLAMPLIMIT;
  }
#endif //BENDCLAMPING

  // auto PC         = pc_valgrad(strain[3], strain[4], strain[5]);
  auto PC         = pc_valgrad(sbend[0],sbend[1],sbend[2]);
  Mat3x3 &PC_grad = std::get<0>(PC);
  Vec3 &PC_val    = std::get<1>(PC);
  double l1       = PC_val(0);
  double l2       = PC_val(1);
  double cc       = PC_val(2);

#ifdef BENDCLAMPINGPC
  if (std::abs(l1) > BENDCLAMPLIMIT)
    l1 = sgn(l1) * BENDCLAMPLIMIT;
  if (std::abs(l2) > BENDCLAMPLIMIT)
    l2 = sgn(l2) * BENDCLAMPLIMIT;
#endif //BENDCLAMPINGPC

  Vec3 S;
  for (int i = 0; i < 3; i++)
    S(i) = (strain(i)) / this->strainscale[i];

  // 1D
  #ifndef NO1D
  for (auto &s : hsplines_1d) {
    double invsc = 1.0 / this->strainscale[s.k];
    if (s.k < 3) {  // in-plane
      double x = S(s.k);
      grad(s.k) += s.dx(x) * invsc;
    } else if (s.k == 3) {  // x-bending
      // d(l1,l2,cc)
      gradpc(0) += s_kx * cc * s.dx(l1 * invsc) * invsc;
      gradpc(1) += s_kx * (1 - cc) * s.dx(l2 * invsc) * invsc;
      gradpc(2) += s_kx * (s.eval(l1 * invsc) - s.eval(l2 * invsc));

    } else {  // y-bending
      // d(l1,l2,cc)
      gradpc(0) += s_ky * (1 - cc) * s.dx(l1 * invsc) * invsc;
      gradpc(1) += s_ky * cc * s.dx(l2 * invsc) * invsc;
      gradpc(2) += s_ky * (s.eval(l2 * invsc) - s.eval(l1 * invsc));
    }
  }
  #endif

  // 2D
  for (auto &s : hsplines_2d) {
    if (!select_2D(s.k0,s.k1))
      continue;

    double invsc0 = 1.0 / this->strainscale[s.k0];
    double invsc1 = 1.0 / this->strainscale[s.k1];
    // assume sorted, k1 > k0
    assert(s.k1 > s.k0 && s.k0 < 3);
    if (s.k1 < 3) {  // in plane
      double x = S(s.k0);
      double y = S(s.k1);
      grad(s.k0) += s.dx(x, y) * invsc0;
      grad(s.k1) += s.dy(x, y) * invsc1;
    } else {  // since we dont have double curve, k0 has to be in plane
      double x     = S(s.k0);
      double dxl1 = s.dx(x, l1 * invsc1);
      double dxl2 = s.dx(x, l2 * invsc1);
      double dyl1 = s.dy(x, l1 * invsc1);
      double dyl2 = s.dy(x, l2 * invsc1);
      if (s.k1 == 3) {
        grad(s.k0) += s_kx * (cc * dxl1 + (1 - cc) * dxl2) * invsc0;
        gradpc(0) += s_kx * cc * dyl1 * invsc1;
        gradpc(1) += s_kx * (1 - cc) * dyl2 * invsc1;
        gradpc(2) += s_kx * (s.eval(x, l1 * invsc1) - s.eval(x, l2 * invsc1));
      } else {
        grad(s.k0) += s_ky * (cc * dxl2 + (1 - cc) * dxl1) * invsc0;
        gradpc(0) += s_ky * (1 - cc) * dyl1 * invsc1;
        gradpc(1) += s_ky * cc * dyl2 * invsc1;
        gradpc(2) += s_ky * (s.eval(x, l2 * invsc1) - s.eval(x, l1 * invsc1));
      }
    }
  }


#ifdef BENDBARRIER
  double mult = BENDPOWER; 
  double T = BENDLIMIT;
  if (std::abs(l1) > T) {
    double r = ((l1 - sgn(l1) * T) / (BENDPERCENT*T));
    double dr = 1/(BENDPERCENT*T);
    gradpc(0) += mult * 2 * r * dr;
  }
  if (std::abs(l2) > T) {
    double r = ((l2 - sgn(l2) * T) / (BENDPERCENT*T));
    double dr = 1/(BENDPERCENT*T);
    gradpc(1) += mult * 2 * r * dr;
  }
#endif  // BENDBARRIER

  // grad_pc Psi -> grad_pc(k) Psi . dpc(k)/dlam(i)
  for (int i = 0; i < 3; ++i) {
    for (int k = 0; k < 3; ++k) {
      grad(3 + i) += gradpc(k) * PC_grad(k, i);
    }
  }


#ifdef DOUBLECURVEPENALTY
  double k = (sbend[0]*sbend[2] - sbend[1]*sbend[1]);
  // val += DOUBLECURVEPENALTY*k*k;
  grad(3) += DOUBLECURVEPENALTY * 2 * k * sbend[2];
  grad(4) += DOUBLECURVEPENALTY * 2 * k * -2*sbend[1];
  grad(5) += DOUBLECURVEPENALTY * 2 * k * sbend[0] + 10;
#endif

  return grad;
}

std::pair<Mat6x6, Vec6> SplineMaterial::psi_drv(const Vec6 &strain) {
  Vec6 grad(0);
  Mat6x6 hess(0);
  if (!initialized)
    return std::make_pair(hess, grad);

  Vec3 gradpc(0);
  Mat3x3 hesspc(0);     // d2 / dpc dpc
  Mat3x3 hessmixed(0);  // d2 / ds dpc
  Vec3 sbend(0);

#ifndef BENDCLAMPING
  for (int i = 0; i < 3; i++)
    sbend[i] = strain[3+i];
#else
  for (int i = 0; i < 3; i++) {
    sbend[i] = strain[3+i];
    if (std::abs(sbend[i]) > BENDCLAMPLIMIT)
      sbend[i] = sgn(sbend[i]) * BENDCLAMPLIMIT;
  }
#endif //BENDCLAMPING

  auto PC = pc_valdrv(sbend[0],sbend[1],sbend[2]);

  std::vector<Mat3x3> &PC_hess = std::get<0>(PC);
  Mat3x3 &PC_grad              = std::get<1>(PC);
  Vec3 &PC_val                 = std::get<2>(PC);
  double l1                    = PC_val(0);
  double l2                    = PC_val(1);
  double cc                    = PC_val(2);


#ifdef BENDCLAMPINGPC
  if (std::abs(l1) > BENDCLAMPLIMIT)
    l1 = sgn(l1) * BENDCLAMPLIMIT;
  if (std::abs(l2) > BENDCLAMPLIMIT)
    l2 = sgn(l2) * BENDCLAMPLIMIT;
#endif //BENDCLAMPINGPC

  Vec3 S;
  for (int i = 0; i < 3; i++)
    S(i) = (strain(i)) / this->strainscale[i];

  // 1D
  #ifndef NO1D
  for (auto &s : hsplines_1d) {
    double invsc = 1.0 / this->strainscale[s.k];
    if (s.k < 3) {  // in-plane
      double x = S(s.k);
      grad(s.k) += s.dx(x) * invsc;
      hess(s.k, s.k) += s.dxdx(x) * invsc * invsc;
    } else {
      double dxl1 = s.dx(l1 * invsc);
      double dxl2 = s.dx(l2 * invsc);

      if (s.k == 3) {  // x-bending
        // d(l1,l2,cc)
        gradpc(0) += s_kx * cc * dxl1 * invsc;
        gradpc(1) += s_kx * (1 - cc) * dxl2 * invsc;
        gradpc(2) += s_kx * (s.eval(l1 * invsc) - s.eval(l2 * invsc));
        // only filling nonzero upper triangular
        hesspc(0, 0) += s_kx * cc * s.dxdx(l1 * invsc) * invsc * invsc;
        // hesspc(0, 1) += 0;
        hesspc(0, 2) += s_kx * dxl1 * invsc;
        hesspc(1, 1) += s_kx * (1 - cc) * s.dxdx(l2 * invsc) * invsc * invsc;
        hesspc(1, 2) += s_kx * -dxl2 * invsc;
        // hesspc(2, 2) += 0;
      } else {  // y-bending
        // d(l1,l2,cc)
        gradpc(0) += s_ky * (1 - cc) * dxl1 * invsc;
        gradpc(1) += s_ky * cc * dxl2 * invsc;
        gradpc(2) += s_ky * (s.eval(l2 * invsc) - s.eval(l1 * invsc));
        hesspc(0, 0) += s_ky * (1 - cc) * s.dxdx(l1 * invsc) * invsc * invsc;
        hesspc(0, 2) += s_ky * -dxl1 * invsc;
        hesspc(1, 1) += s_ky * cc * s.dxdx(l2 * invsc) * invsc * invsc;
        hesspc(1, 2) += s_ky * dxl2 * invsc;
      }
    }
  }
  #endif

  // 2D
  for (auto &s : hsplines_2d) {
    if (!select_2D(s.k0,s.k1))
      continue;

    double invsc0 = 1.0 / this->strainscale[s.k0];
    double invsc1 = 1.0 / this->strainscale[s.k1];
    // assume sorted, k1 > k0
    assert(s.k1 > s.k0 && s.k0 < 3);
    if (s.k1 < 3) {  // in plane
      double x = S(s.k0);
      double y = S(s.k1);
      grad(s.k0) += s.dx(x, y) * invsc0;
      grad(s.k1) += s.dy(x, y) * invsc1;
      hess(s.k0, s.k0) += s.dxdx(x, y) * invsc0 * invsc0;
      double dxdy = s.dxdy(x, y) * invsc0 * invsc1;
      hess(s.k0, s.k1) += dxdy;
      hess(s.k1, s.k0) += dxdy;
      hess(s.k1, s.k1) += s.dydy(x, y) * invsc1 * invsc1;
    } else {  // since we dont have double curve, k0 has to be in plane

      double x    = S(s.k0);

      double dxl1 = s.dx(x, l1 * invsc1);
      double dxl2 = s.dx(x, l2 * invsc1);
      double dyl1 = s.dy(x, l1 * invsc1);
      double dyl2 = s.dy(x, l2 * invsc1);
      if (s.k1 == 3) {
        //   val += cc * s.eval(x, l1 / this->strainscale[s.k1]) +
        //          (1 - cc) * s.eval(x, l2 / this->strainscale[s.k1]);

        grad(s.k0) += s_kx * (cc * dxl1 + (1 - cc) * dxl2) * invsc0;
        gradpc(0) += s_kx * cc * dyl1 * invsc1;
        gradpc(1) += s_kx * (1 - cc) * dyl2 * invsc1;
        gradpc(2) += s_kx * (s.eval(x, l1 * invsc1) - s.eval(x, l2 * invsc1));

        hess(s.k0, s.k0) +=
            s_kx *
            (cc * s.dxdx(x, l1 * invsc1) + (1 - cc) * s.dxdx(x, l2 * invsc1)) *
            invsc0 * invsc0;

        hessmixed(s.k0, 0) += s_kx * cc * s.dxdy(x, l1 * invsc1) * invsc1 * invsc0;
        hessmixed(s.k0, 1) += s_kx * (1 - cc) * s.dxdy(x, l2 * invsc1) * invsc1 * invsc0;
        hessmixed(s.k0, 2) += s_kx * (dxl1 - dxl2) * invsc0;

        hesspc(0, 0) += s_kx * cc * s.dydy(x, l1 * invsc1) * invsc1 * invsc1;
        hesspc(0, 2) += s_kx * dyl1 * invsc1;
        hesspc(1, 1) +=
            s_kx * (1 - cc) * s.dydy(x, l2 * invsc1) * invsc1 * invsc1;
        hesspc(1, 2) += s_kx * -dyl2 * invsc1;
      } else {
        //   val += cc * s.eval(x, l2 / this->strainscale[s.k1]) +
        //          (1 - cc) * s.eval(x, l1 / this->strainscale[s.k1]);
        grad(s.k0) += s_ky * (cc * dxl2 + (1 - cc) * dxl1) * invsc0;
        gradpc(0) += s_ky * (1 - cc) * dyl1 * invsc1;  // this adds jitter..?
        gradpc(1) += s_ky * cc * dyl2 * invsc1;
        gradpc(2) += s_ky * (s.eval(x, l2 * invsc1) - s.eval(x, l1 * invsc1));

        hess(s.k0, s.k0) +=
            s_ky *
            (cc * s.dxdx(x, l2 * invsc1) + (1 - cc) * s.dxdx(x, l1 * invsc1)) *
            invsc0 * invsc0;

        hessmixed(s.k0, 0) += s_ky * (1 - cc) * s.dxdy(x, l1 * invsc1) * invsc0 * invsc1;
        hessmixed(s.k0, 1) += s_ky * cc * s.dxdy(x, l2 * invsc1) * invsc0 * invsc1;
        hessmixed(s.k0, 2) += s_ky * (dxl2 - dxl1) * invsc0;

        hesspc(0, 0) +=
            s_ky * (1 - cc) * s.dydy(x, l1 * invsc1) * invsc1 * invsc1;
        hesspc(0, 2) += s_ky * -dyl1 * invsc1;
        hesspc(1, 1) += s_ky * cc * s.dydy(x, l2 * invsc1) * invsc1 * invsc1;
        hesspc(1, 2) += s_ky * dyl2 * invsc1;
      }
    }
  }

#ifdef BENDBARRIER
  double mult = BENDPOWER; 
  double T = BENDLIMIT;
  if (std::abs(l1) > T) {
    double r = ((l1 - sgn(l1) * T) / (BENDPERCENT*T));
    double dr = 1/(BENDPERCENT*T);
    gradpc(0) += mult * 2 * r * dr;
    hesspc(0,0) += mult * 2 * dr * dr;
  }
  if (std::abs(l2) > T) {
    double r = ((l2 - sgn(l2) * T) / (BENDPERCENT*T));
    double dr = 1/(BENDPERCENT*T);
    gradpc(1) += mult * 2 * r * dr;
    hesspc(1,1) += mult * 2 * dr * dr;
  }
#endif  // BENDBARRIER


  // grad_pc Psi -> grad_pc(k) Psi . dPC(k)/dL(i)
  for (int i = 0; i < 3; ++i) {
    for (int k = 0; k < 3; ++k) {
      grad(3 + i) += gradpc(k) * PC_grad(k, i);
    }
  }

  // by now we have topleft hess as symmetric dS dS
  // hessmixed as uppertri dS dPC
  // hesspc as uppertri dPC dPC

  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      for (int k = 0; k < 3; ++k) {

        // dS(i)dL(j) = dS(i)dPC(k) . dPC(k)/dL(j)
        hess(i, 3 + j) += hessmixed(i, k) * PC_grad(k, j);

        if (i > j) // upper tri, symmetrize later
          continue;

        // dL(i)dL(j) = dPC(k) . dPC(k)/dL(i)dL(j)
        //           + dPC(k)dPC(l) . dPC(k)/dL(i) . dPC(l)/dL(j)
        hess(3 + i, 3 + j) += gradpc(k) * PC_hess[k](i, j);
        for (int l = 0; l < 3; ++l) {
          double dPCkdPCl = (l > k) ? hesspc(k, l) : hesspc(l, k);  // upper tri
          hess(3 + i, 3 + j) += dPCkdPCl * PC_grad(k, i) * PC_grad(l, j);
        }
      }
    }
  }



#ifdef DOUBLECURVEPENALTY
  double k = (sbend[0]*sbend[2] - sbend[1]*sbend[1]);
  // val += DOUBLECURVEPENALTY*k*k;
  grad(3) += DOUBLECURVEPENALTY * 2 * k * sbend[2];
  grad(4) += DOUBLECURVEPENALTY * 2 * k * -2*sbend[1];
  grad(5) += DOUBLECURVEPENALTY * 2 * k * sbend[0];
  hess(3,3) += DOUBLECURVEPENALTY * 2 * sbend[2]*sbend[2];
  hess(3,4) += DOUBLECURVEPENALTY * -4 * sbend[1]*sbend[2];
  hess(3,5) += DOUBLECURVEPENALTY * -2 * (sbend[1] * sbend[1] - 2*sbend[0]*sbend[2]);
  hess(4,4) += DOUBLECURVEPENALTY * 4 * (3 * sbend[1] * sbend[1] - sbend[0]*sbend[2]);
  hess(4,5) += DOUBLECURVEPENALTY * -4 * sbend[1]*sbend[0];
  hess(5,5) += DOUBLECURVEPENALTY * 2 * sbend[0]*sbend[0];
#endif

  // symmetrize
  for (int i = 0; i < 6; ++i)
    for (int j = i + 1; j < 6; ++j)
      hess(j,i) = hess(i,j);

  // for (int i = 0; i < 6; ++i)
  //   hess(i,i)+= (i<3)? 1e-3 : 1e-3;

  return std::make_pair(hess, grad);
}

using namespace hylc::mathematica;

Vec3 SplineMaterial::pc_val(double lxx, double lxy, double lyy) {
  Vec3 out;
  Real sgn               = lxx < lyy ? -1 : 1;
  static const Real eps2 = 1e-15;

  Real copt1  = Power(lxy, 2);
  Real copt2  = -lyy;
  Real copt3  = copt2 + lxx;
  Real copt4  = Power(copt3, 2);
  Real copt5  = copt4 / 4.;
  Real copt6  = copt1 + copt5 + eps2;
  Real copt7  = Sqrt(copt6);
  Real copt8  = lxx + lyy;
  Real copt9  = copt8 / 2.;
  Real copt13 = copt3 / 2.;
  Real copt14 = copt7 * sgn;
  Real copt15 = copt13 + copt14;
  Real copt16 = Power(copt15, 2);
  Real copt17 = copt1 + copt16;
  Real copt18 = 1 / copt17;
  out(0)      = copt7 + copt9;
  out(1)      = -copt7 + copt9;
  out(2)      = 0.5 - (-0.5 + copt1 * copt18) * sgn;

  return out;
}

std::tuple<Mat3x3, Vec3> SplineMaterial::pc_valgrad(double lxx, double lxy,
                                                    double lyy) {
  Mat3x3 grad(0);
  Vec3 val(0);
  auto out1              = [&](int i) -> Real & { return val[i]; };
  auto out2              = [&](int i, int j) -> Real & { return grad(i, j); };
  Real sgn               = lxx < lyy ? -1 : 1;
  static const Real eps2 = 1e-15;

  Real copt1  = Power(lxy, 2);
  Real copt2  = -lyy;
  Real copt3  = copt2 + lxx;
  Real copt4  = Power(copt3, 2);
  Real copt5  = copt4 / 4.;
  Real copt6  = copt1 + copt5 + eps2;
  Real copt7  = Sqrt(copt6);
  Real copt8  = lxx + lyy;
  Real copt9  = copt8 / 2.;
  Real copt13 = copt3 / 2.;
  Real copt14 = copt7 * sgn;
  Real copt15 = copt13 + copt14;
  Real copt16 = Power(copt15, 2);
  Real copt17 = copt1 + copt16;
  Real copt18 = 1 / copt17;
  Real copt23 = 1 / copt7;
  Real copt27 = -lxx;
  Real copt28 = copt27 + lyy;
  Real copt29 = (copt23 * copt28) / 4.;
  Real copt30 = 0.5 + copt29;
  Real copt24 = (copt23 * copt3) / 4.;
  Real copt25 = 0.5 + copt24;
  Real copt34 = Power(copt17, 2);
  Real copt35 = 1 / copt34;
  Real copt37 = copt1 * lxy;
  Real copt38 = 4 * eps2;
  Real copt39 = Power(lxx, 2);
  Real copt40 = 4 * copt1;
  Real copt41 = -2 * lxx * lyy;
  Real copt42 = Power(lyy, 2);
  Real copt43 = copt38 + copt39 + copt40 + copt41 + copt42;
  Real copt44 = Sqrt(copt43);
  Real copt45 = 1 / copt44;
  out1(0)     = copt7 + copt9;
  out1(1)     = -copt7 + copt9;
  out1(2)     = 0.5 - (-0.5 + copt1 * copt18) * sgn;
  out2(0, 0)  = copt25;
  out2(0, 1)  = copt23 * lxy;
  out2(0, 2)  = copt30;
  out2(1, 0)  = copt30;
  out2(1, 1)  = -(copt23 * lxy);
  out2(1, 2)  = copt25;
  out2(2, 0) =
      2 * copt1 * copt15 * copt35 * sgn * (0.5 + (copt23 * copt3 * sgn) / 4.);
  out2(2, 1) = -2 * copt18 * lxy * sgn +
               2 * copt35 * copt37 * sgn * (2 + copt3 * copt45 * sgn);
  out2(2, 2) =
      2 * copt1 * copt15 * copt35 * sgn * (-0.5 + (copt23 * copt28 * sgn) / 4.);

  return std::make_tuple(grad, val);
}

std::tuple<std::vector<Mat3x3>, Mat3x3, Vec3> SplineMaterial::pc_valdrv(
    double lxx, double lxy, double lyy) {
  std::vector<Mat3x3> hess(3);  // 6x18x18
  Mat3x3 grad(0);
  Vec3 val(0);
  auto out1 = [&](int i) -> Real & { return val[i]; };
  auto out2 = [&](int i, int j) -> Real & { return grad(i, j); };
  auto out3 = [&](int i, int j, int k) -> Real & { return hess[i](j, k); };
  Real sgn  = lxx < lyy ? -1 : 1;
  static const Real eps2 = 1e-15;

  Real copt1   = Power(lxy, 2);
  Real copt2   = -lyy;
  Real copt3   = copt2 + lxx;
  Real copt4   = Power(copt3, 2);
  Real copt5   = copt4 / 4.;
  Real copt6   = copt1 + copt5 + eps2;
  Real copt7   = Sqrt(copt6);
  Real copt8   = lxx + lyy;
  Real copt9   = copt8 / 2.;
  Real copt13  = copt3 / 2.;
  Real copt14  = copt7 * sgn;
  Real copt15  = copt13 + copt14;
  Real copt16  = Power(copt15, 2);
  Real copt17  = copt1 + copt16;
  Real copt18  = 1 / copt17;
  Real copt23  = 1 / copt7;
  Real copt27  = -lxx;
  Real copt28  = copt27 + lyy;
  Real copt29  = (copt23 * copt28) / 4.;
  Real copt30  = 0.5 + copt29;
  Real copt24  = (copt23 * copt3) / 4.;
  Real copt25  = 0.5 + copt24;
  Real copt34  = Power(copt17, 2);
  Real copt35  = 1 / copt34;
  Real copt37  = copt1 * lxy;
  Real copt38  = 4 * eps2;
  Real copt39  = Power(lxx, 2);
  Real copt40  = 4 * copt1;
  Real copt41  = -2 * lxx * lyy;
  Real copt42  = Power(lyy, 2);
  Real copt43  = copt38 + copt39 + copt40 + copt41 + copt42;
  Real copt44  = Sqrt(copt43);
  Real copt45  = 1 / copt44;
  Real copt55  = copt43 * copt44;
  Real copt56  = 1 / copt55;
  Real copt58  = copt6 * copt7;
  Real copt59  = 1 / copt58;
  Real copt54  = copt1 + eps2;
  Real copt60  = -(copt3 * copt59 * lxy) / 4.;
  Real copt61  = -2 * copt54 * copt56;
  Real copt64  = 2 * copt3 * copt56 * lxy;
  Real copt57  = 2 * copt54 * copt56;
  Real copt62  = copt38 + copt4;
  Real copt70  = 2 * eps2;
  Real copt71  = copt44 * lxx * sgn;
  Real copt72  = -(copt44 * lyy * sgn);
  Real copt73  = copt39 + copt40 + copt41 + copt42 + copt70 + copt71 + copt72;
  Real copt74  = Power(copt73, 2);
  Real copt75  = 1 / copt74;
  Real copt77  = 2 * copt1;
  Real copt78  = copt39 + copt41 + copt42 + copt70 + copt77;
  Real copt79  = copt45 * copt78 * sgn;
  Real copt80  = copt2 + copt79 + lxx;
  Real copt81  = Power(copt80, 2);
  Real copt82  = copt17 * copt34;
  Real copt83  = 1 / copt82;
  Real copt32  = (copt23 * copt3 * sgn) / 4.;
  Real copt33  = 0.5 + copt32;
  Real copt85  = Power(copt33, 2);
  Real copt88  = copt73 * copt74;
  Real copt89  = 1 / copt88;
  Real copt90  = Power(eps2, 2);
  Real copt91  = copt90 * eps2;
  Real copt94  = Power(copt39, 2);
  Real copt96  = copt39 * lxx;
  Real copt97  = Power(copt96, 2);
  Real copt101 = Power(copt1, 2);
  Real copt104 = Power(copt37, 2);
  Real copt108 = copt94 * lxx;
  Real copt119 = copt42 * lyy;
  Real copt123 = Power(copt42, 2);
  Real copt127 = copt123 * lyy;
  Real copt129 = Power(copt119, 2);
  Real copt46  = copt3 * copt45 * sgn;
  Real copt47  = 2 + copt46;
  Real copt186 = Power(copt47, 2);
  Real copt51  = (copt23 * copt28 * sgn) / 4.;
  Real copt52  = -0.5 + copt51;
  Real copt149 = 18 * copt90 * lxx;
  Real copt150 = 16 * copt96 * eps2;
  Real copt151 = 3 * copt108;
  Real copt152 = 30 * copt1 * eps2 * lxx;
  Real copt153 = 15 * copt1 * copt96;
  Real copt154 = 12 * copt101 * lxx;
  Real copt155 = -18 * copt90 * lyy;
  Real copt156 = -48 * copt39 * eps2 * lyy;
  Real copt157 = -15 * copt94 * lyy;
  Real copt158 = -30 * copt1 * eps2 * lyy;
  Real copt159 = -45 * copt1 * copt39 * lyy;
  Real copt160 = -12 * copt101 * lyy;
  Real copt161 = 48 * copt42 * eps2 * lxx;
  Real copt162 = 30 * copt42 * copt96;
  Real copt163 = 45 * copt1 * copt42 * lxx;
  Real copt164 = -16 * copt119 * eps2;
  Real copt165 = -30 * copt119 * copt39;
  Real copt166 = -15 * copt1 * copt119;
  Real copt167 = 15 * copt123 * lxx;
  Real copt168 = -3 * copt127;
  Real copt169 = 4 * copt90;
  Real copt170 = 3 * copt1;
  Real copt171 = copt170 + copt39 + copt41 + copt42;
  Real copt172 = 3 * copt171 * copt4;
  Real copt173 = 5 * copt39;
  Real copt174 = -10 * lxx * lyy;
  Real copt175 = 5 * copt42;
  Real copt176 = copt173 + copt174 + copt175 + copt77;
  Real copt177 = 2 * copt176 * eps2;
  Real copt178 = copt169 + copt172 + copt177;
  Real copt179 = copt178 * copt44 * sgn;
  Real copt180 = copt149 + copt150 + copt151 + copt152 + copt153 + copt154 +
                 copt155 + copt156 + copt157 + copt158 + copt159 + copt160 +
                 copt161 + copt162 + copt163 + copt164 + copt165 + copt166 +
                 copt167 + copt168 + copt179;
  Real copt181 = (copt1 * copt180 * copt59 * copt83) / 8.;
  Real copt92  = 8 * copt91;
  Real copt93  = 18 * copt39 * copt90;
  Real copt95  = 8 * copt94 * eps2;
  Real copt98  = 4 * copt1 * copt90;
  Real copt99  = 12 * copt1 * copt39 * eps2;
  Real copt100 = 3 * copt1 * copt94;
  Real copt102 = -12 * copt101 * eps2;
  Real copt103 = -6 * copt101 * copt39;
  Real copt105 = -8 * copt104;
  Real copt106 = -36 * copt90 * lxx * lyy;
  Real copt107 = -32 * copt96 * eps2 * lyy;
  Real copt109 = -6 * copt108 * lyy;
  Real copt110 = -24 * copt1 * eps2 * lxx * lyy;
  Real copt111 = -12 * copt1 * copt96 * lyy;
  Real copt112 = 12 * copt101 * lxx * lyy;
  Real copt113 = 18 * copt42 * copt90;
  Real copt114 = 48 * copt39 * copt42 * eps2;
  Real copt115 = 15 * copt42 * copt94;
  Real copt116 = 12 * copt1 * copt42 * eps2;
  Real copt117 = 18 * copt1 * copt39 * copt42;
  Real copt118 = -6 * copt101 * copt42;
  Real copt120 = -32 * copt119 * eps2 * lxx;
  Real copt121 = -20 * copt119 * copt96;
  Real copt122 = -12 * copt1 * copt119 * lxx;
  Real copt124 = 8 * copt123 * eps2;
  Real copt125 = 15 * copt123 * copt39;
  Real copt126 = 3 * copt1 * copt123;
  Real copt128 = -6 * copt127 * lxx;
  Real copt130 = 8 * copt90;
  Real copt131 = -6 * copt101;
  Real copt132 = -4 * copt96 * lyy;
  Real copt133 = copt1 * copt42;
  Real copt134 = 2 * copt42;
  Real copt135 = copt1 + copt134;
  Real copt136 = -2 * copt135 * lxx * lyy;
  Real copt137 = 3 * copt39;
  Real copt138 = -6 * lxx * lyy;
  Real copt139 = 3 * copt42;
  Real copt140 = copt1 + copt137 + copt138 + copt139;
  Real copt141 = 2 * copt140 * eps2;
  Real copt142 = 6 * copt42;
  Real copt143 = copt1 + copt142;
  Real copt144 = copt143 * copt39;
  Real copt145 = copt123 + copt130 + copt131 + copt132 + copt133 + copt136 +
                 copt141 + copt144 + copt94;
  Real copt146 = copt145 * copt3 * copt44 * sgn;
  Real copt147 = copt100 + copt102 + copt103 + copt105 + copt106 + copt107 +
                 copt109 + copt110 + copt111 + copt112 + copt113 + copt114 +
                 copt115 + copt116 + copt117 + copt118 + copt120 + copt121 +
                 copt122 + copt124 + copt125 + copt126 + copt128 + copt129 +
                 copt146 + copt92 + copt93 + copt95 + copt97 + copt98 + copt99;
  Real copt68   = copt44 * sgn;
  Real copt69   = copt2 + copt68 + lxx;
  Real copt76   = 8 * copt1 * copt54 * copt56 * copt69 * copt75;
  Real copt84   = -2 * copt1 * copt81 * copt83 * sgn;
  Real copt197  = copt3 * copt45 * lxy * sgn;
  Real copt198  = copt197 + lxy;
  Real copt199  = Power(copt198, 2);
  out1(0)       = copt7 + copt9;
  out1(1)       = -copt7 + copt9;
  out1(2)       = 0.5 - (-0.5 + copt1 * copt18) * sgn;
  out2(0, 0)    = copt25;
  out2(0, 1)    = copt23 * lxy;
  out2(0, 2)    = copt30;
  out2(1, 0)    = copt30;
  out2(1, 1)    = -(copt23 * lxy);
  out2(1, 2)    = copt25;
  out2(2, 0)    = 2 * copt1 * copt15 * copt33 * copt35 * sgn;
  out2(2, 1)    = 2 * copt35 * copt37 * copt47 * sgn - 2 * copt18 * lxy * sgn;
  out2(2, 2)    = 2 * copt1 * copt15 * copt35 * copt52 * sgn;
  out3(0, 0, 0) = copt57;
  out3(0, 0, 1) = copt60;
  out3(0, 0, 2) = copt61;
  out3(0, 1, 0) = copt60;
  out3(0, 1, 1) = 2 * copt56 * copt62;
  out3(0, 1, 2) = copt64;
  out3(0, 2, 0) = copt61;
  out3(0, 2, 1) = copt64;
  out3(0, 2, 2) = copt57;
  out3(1, 0, 0) = copt61;
  out3(1, 0, 1) = copt64;
  out3(1, 0, 2) = copt57;
  out3(1, 1, 0) = (copt3 * copt59 * lxy) / 4.;
  out3(1, 1, 1) = -2 * copt56 * copt62;
  out3(1, 1, 2) = (copt28 * copt59 * lxy) / 4.;
  out3(1, 2, 0) = copt57;
  out3(1, 2, 1) = copt60;
  out3(1, 2, 2) = copt61;
  out3(2, 0, 0) = copt76 + copt84 + 2 * copt1 * copt35 * copt85 * sgn;
  out3(2, 0, 1) = 16 * copt147 * copt56 * copt89 * lxy;
  out3(2, 0, 2) = copt181;
  out3(2, 1, 0) = 32 * copt37 * copt54 * copt56 * copt75 -
                  8 * copt15 * copt33 * copt37 * copt47 * copt83 * sgn +
                  4 * copt15 * copt33 * copt35 * lxy * sgn;
  out3(2, 1, 1) = -8 * copt101 * copt3 * copt35 * copt56 - 2 * copt18 * sgn +
                  10 * copt1 * copt35 * copt47 * sgn -
                  8 * copt101 * copt186 * copt83 * sgn;
  out3(2, 1, 2) = -32 * copt37 * copt54 * copt56 * copt75 -
                  8 * copt15 * copt37 * copt47 * copt52 * copt83 * sgn +
                  4 * copt15 * copt35 * copt52 * lxy * sgn;
  out3(2, 2, 0) = copt181;
  out3(2, 2, 1) = -16 * copt147 * copt56 * copt89 * lxy;
  out3(2, 2, 2) = copt76 + copt84 + (copt199 * copt35 * sgn) / 2.;

  return std::make_tuple(hess, grad, val);
}