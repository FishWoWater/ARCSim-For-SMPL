#ifndef _POLY2D_H_
#define _POLY2D_H_

#include <vector>

struct Poly2D {
  int k0, k1;
  std::vector<double> c;
  double xmin, xmax;
  bool clampx;
  double ymin, ymax;
  bool clampy;

  inline void clamp(double &x, double &y) {
    if (clampx)
      x = std::min(std::max(x, xmin), xmax);
    if (clampy)
      y = std::min(std::max(y, ymin), ymax);
  }

  // pows = [(1,1),(1,2),(2,1),(1,3),(2,2),(3,1)]
  double eval(double x, double y) {
    double xx  = x * x;
    double xxx = xx * x;
    double yy  = y * y;
    double yyy = yy * y;
    return c[0] * x * y + c[1] * x * yy + c[2] * xx * y + c[3] * x * yyy +
           c[4] * xx * yy + c[5] * xxx * y;
  }

  double dx(double x, double y) {
    double xx  = x * x;
    double yy  = y * y;
    double yyy = yy * y;
    return c[0] * y + c[1] * yy + c[2] * 2 * x * y + c[3] * yyy +
           c[4] * 2 * x * yy + c[5] * 3 * xx * y;
  }

  double dy(double x, double y) {
    double xx  = x * x;
    double xxx = xx * x;
    double yy  = y * y;
    return c[0] * x + c[1] * x * 2 * y + c[2] * xx + c[3] * x * 3 * yy +
           c[4] * xx * 2 * y + c[5] * xxx;
  }

  double dxdx(double x, double y) {
    return c[2] * 2 * y + c[4] * 2 * y * y + c[5] * 6 * x * y;
  }

  double dxdy(double x, double y) {
    return c[0] + c[1] * 2 * y + c[2] * 2 * x + c[3] * 3 * y * y +
           c[4] * x * 4 * y + c[5] * 3 * x * x;
  }

  double dydy(double x, double y) {
    return c[1] * x * 2 + c[3] * x * 6 * y + c[4] * x * x * 2;
  }
};

#endif  // _POLY2D_H_