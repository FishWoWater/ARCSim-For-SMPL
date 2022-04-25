#ifndef _HSPLINES1_H_
#define _HSPLINES1_H_

#include <algorithm>
#include <vector>

struct HermiteSpline1D {
  int k;
  std::vector<double> t, p, m;
  int ext = 0;

  double eval(double x, int der = 0) {
    if (t.size() == 0)
      return 0;

    auto it = std::upper_bound(t.begin(), t.end(),
                               x);  // index to the right of element
    int i   = it - t.begin();

    if (i == 0 || i == (int)t.size()) {  // outside domain
      double tt, pp, mm;
      if (i == 0) {  // left
        tt = t[0];
        pp = p[0];
        mm = m[0];
      } else {  // right
        tt = t[t.size() - 1];
        pp = p[t.size() - 1];
        mm = m[t.size() - 1];
      }
      if (ext == 1)  // const
        mm = 0;
      else if (ext == 2) {  // 0
        pp = 0;
        mm = 0;
      }
      if (der == 0)
        return pp + mm * (x - tt);
      else if (der == 1)
        return mm;
      else
        return 0;
    } else {
      // p(t) = (2t^3 - 3t^2 + 1) p0 + (t^3 - 2t^2 + t) m0 + (-2t^3 + 3t^2)p1 +
      // (t^3 - t^2) m1
      double t0 = t[i - 1], t1 = t[i];
      double p0 = p[i - 1], p1 = p[i];
      double m0 = m[i - 1], m1 = m[i];

      double invdt = 1.0 / (t1 - t0);
      double tt    = (x - t0) * invdt;
      double tt2   = tt * tt;
      double tt3   = tt2 * tt;
      double h00, h10, h01, h11;

      double dtdx;
      if (der == 0) {
        dtdx = 1.0;
        h00  = (2 * tt3 - 3 * tt2 + 1);
        h10  = (tt3 - 2 * tt2 + tt);
        h01  = (-2 * tt3 + 3 * tt2);
        h11  = (tt3 - tt2);
      } else if (der == 1) {
        dtdx = invdt;
        h00  = (6 * tt2 - 6 * tt);
        h10  = (3 * tt2 - 4 * tt + 1);
        h01  = (-6 * tt2 + 6 * tt);
        h11  = (3 * tt2 - 2 * tt);
      } else if (der == 2) {
        dtdx = invdt * invdt;
        h00  = (12 * tt - 6);
        h10  = (6 * tt - 4);
        h01  = (-12 * tt + 6);
        h11  = (6 * tt - 2);
      } else if (der == 3) {
        dtdx = invdt * invdt * invdt;
        h00  = 12;
        h10  = 6;
        h01  = -12;
        h11  = 6;
      } else {
        return 0;
      }
      return dtdx * (h00 * p0 + h10 * (t1 - t0) * m0 + h01 * p1 +
                     h11 * (t1 - t0) * m1);
    }
  }

  double dx(double x) { return eval(x, 1); }
  double dxdx(double x) { return eval(x, 2); }
};

#endif  // _HSPLINES1_H_