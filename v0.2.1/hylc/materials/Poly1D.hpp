#ifndef _POLY1D_H_
#define _POLY1D_H_

#include <vector>

struct Poly1D {
  int k;
  std::vector<double> c;
  double xmin, xmax;
  bool clampx;
  bool compr;

  inline void clamp(double &x) {
    if (clampx)
      x = std::min(std::max(x, xmin), xmax);
  }

  // pows = [1,2,3,4]
  double eval(double x) {
    double xx = x*x;
    double xxx = xx * x;
    double val = 0;

    val = c[0]*x + c[1]*xx + c[2] *xxx + c[3] * xxx*x;
    if (compr) {
      val += c[4] * 1.0 / (x + c[5]);
    }
    return val;
  }

  double dx(double x) {
    double xx = x*x;
    double xxx = xx * x;
    double val = 0;

    val = c[0] + 2*c[1]*x + 3*c[2] *xx + 4 * c[3] * xxx;
    if (compr) {
      double tmp = (x + c[5]);
      val += -c[4] * 1.0 / (tmp*tmp);
    }
    return val;
  }

  double dxdx(double x) {
    double xx = x*x;
    double val = 0;
    val = 2*c[1] + 6*c[2] *x + 12 * c[3] * xx;
    if (compr) {
      double tmp = (x + c[5]);
      val += c[4] * 2.0 / (tmp*tmp*tmp);
    }
    return val;
  }
};

#endif  // _POLY1D_H_