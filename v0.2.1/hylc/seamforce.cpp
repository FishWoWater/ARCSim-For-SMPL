
#include "seamforce.hpp"
#include "MathematicaDefinitions.h"

using namespace hylc;
using namespace hylc::mathematica;

namespace hylc {

double seamforce_val(const Vec<3> &x0, const Vec<3> &x1, double L) {
  Real copt1  = 1 / L;
  Real copt2  = x0(0);
  Real copt3  = -copt2;
  Real copt4  = x1(0);
  Real copt5  = copt3 + copt4;
  Real copt6  = Power(copt5, 2);
  Real copt7  = x0(1);
  Real copt8  = -copt7;
  Real copt9  = x1(1);
  Real copt10 = copt8 + copt9;
  Real copt11 = Power(copt10, 2);
  Real copt12 = x0(2);
  Real copt13 = -copt12;
  Real copt14 = x1(2);
  Real copt15 = copt13 + copt14;
  Real copt16 = Power(copt15, 2);
  Real copt17 = copt11 + copt16 + copt6;  // + eps2;
  Real copt18 = Sqrt(copt17);
  Real copt19 = copt1 * copt18;
  Real copt20 = -1 + copt19;
  Real copt21 = Power(copt20, 2);
  return copt21;
}
Vec<6> seamforce_grad(const Vec<3> &x0, const Vec<3> &x1, double L) {
  double eps2 = 1e-10;
  Vec<6> out(0);

  Real copt1  = 1 / L;
  Real copt2  = x0(0);
  Real copt3  = -copt2;
  Real copt4  = x1(0);
  Real copt5  = copt3 + copt4;
  Real copt6  = Power(copt5, 2);
  Real copt7  = x0(1);
  Real copt8  = -copt7;
  Real copt9  = x1(1);
  Real copt10 = copt8 + copt9;
  Real copt11 = Power(copt10, 2);
  Real copt12 = x0(2);
  Real copt13 = -copt12;
  Real copt14 = x1(2);
  Real copt15 = copt13 + copt14;
  Real copt16 = Power(copt15, 2);
  Real copt17 = copt11 + copt16 + copt6 + eps2;
  Real copt18 = Sqrt(copt17);
  Real copt19 = 1 / copt18;
  Real copt20 = copt1 * copt18;
  Real copt21 = -1 + copt20;
  out(0)      = -2 * copt1 * copt19 * copt21 * copt5;
  out(1)      = -2 * copt1 * copt10 * copt19 * copt21;
  out(2)      = -2 * copt1 * copt15 * copt19 * copt21;
  out(3)      = 2 * copt1 * copt19 * copt21 * copt5;
  out(4)      = 2 * copt1 * copt10 * copt19 * copt21;
  out(5)      = 2 * copt1 * copt15 * copt19 * copt21;

  return out;
}

std::pair<Mat<6, 6>, Vec<6>> seamforce_drv(const Vec<3> &x0, const Vec<3> &x1,
                                           double L) {
  double eps2 = 1e-10;
  Mat<6, 6> out2(0);
  Vec<6> out1(0);

  Real copt1  = 1 / L;
  Real copt2  = x0(0);
  Real copt3  = -copt2;
  Real copt4  = x1(0);
  Real copt5  = copt3 + copt4;
  Real copt6  = Power(copt5, 2);
  Real copt7  = x0(1);
  Real copt8  = -copt7;
  Real copt9  = x1(1);
  Real copt10 = copt8 + copt9;
  Real copt11 = Power(copt10, 2);
  Real copt12 = x0(2);
  Real copt13 = -copt12;
  Real copt14 = x1(2);
  Real copt15 = copt13 + copt14;
  Real copt16 = Power(copt15, 2);
  Real copt17 = copt11 + copt16 + copt6 + eps2;
  Real copt18 = Sqrt(copt17);
  Real copt19 = 1 / copt18;
  Real copt20 = copt1 * copt18;
  Real copt21 = -1 + copt20;
  Real copt28 = Power(L, 2);
  Real copt29 = 1 / copt28;
  Real copt30 = 1 / copt17;
  Real copt32 = copt17 * copt18;
  Real copt33 = 1 / copt32;
  Real copt37 = 2 * copt10 * copt29 * copt30 * copt5;
  Real copt38 = -2 * copt1 * copt10 * copt21 * copt33 * copt5;
  Real copt39 = copt37 + copt38;
  Real copt35 = 2 * copt1 * copt19 * copt21;
  Real copt47 = -2 * copt10 * copt29 * copt30 * copt5;
  Real copt48 = 2 * copt1 * copt10 * copt21 * copt33 * copt5;
  Real copt49 = copt47 + copt48;
  Real copt45 = -2 * copt1 * copt19 * copt21;
  Real copt40 = 2 * copt15 * copt29 * copt30 * copt5;
  Real copt41 = -2 * copt1 * copt15 * copt21 * copt33 * copt5;
  Real copt42 = copt40 + copt41;
  Real copt56 = 2 * copt10 * copt15 * copt29 * copt30;
  Real copt57 = -2 * copt1 * copt10 * copt15 * copt21 * copt33;
  Real copt58 = copt56 + copt57;
  Real copt50 = -2 * copt15 * copt29 * copt30 * copt5;
  Real copt51 = 2 * copt1 * copt15 * copt21 * copt33 * copt5;
  Real copt52 = copt50 + copt51;
  Real copt62 = -2 * copt10 * copt15 * copt29 * copt30;
  Real copt63 = 2 * copt1 * copt10 * copt15 * copt21 * copt33;
  Real copt64 = copt62 + copt63;
  Real copt43 = -2 * copt29 * copt30 * copt6;
  Real copt44 = 2 * copt1 * copt21 * copt33 * copt6;
  Real copt46 = copt43 + copt44 + copt45;
  Real copt31 = 2 * copt29 * copt30 * copt6;
  Real copt34 = -2 * copt1 * copt21 * copt33 * copt6;
  Real copt36 = copt31 + copt34 + copt35;
  Real copt59 = -2 * copt11 * copt29 * copt30;
  Real copt60 = 2 * copt1 * copt11 * copt21 * copt33;
  Real copt61 = copt45 + copt59 + copt60;
  Real copt53 = 2 * copt11 * copt29 * copt30;
  Real copt54 = -2 * copt1 * copt11 * copt21 * copt33;
  Real copt55 = copt35 + copt53 + copt54;
  Real copt68 = -2 * copt16 * copt29 * copt30;
  Real copt69 = 2 * copt1 * copt16 * copt21 * copt33;
  Real copt70 = copt45 + copt68 + copt69;
  Real copt65 = 2 * copt16 * copt29 * copt30;
  Real copt66 = -2 * copt1 * copt16 * copt21 * copt33;
  Real copt67 = copt35 + copt65 + copt66;
  out1(0)     = -2 * copt1 * copt19 * copt21 * copt5;
  out1(1)     = -2 * copt1 * copt10 * copt19 * copt21;
  out1(2)     = -2 * copt1 * copt15 * copt19 * copt21;
  out1(3)     = 2 * copt1 * copt19 * copt21 * copt5;
  out1(4)     = 2 * copt1 * copt10 * copt19 * copt21;
  out1(5)     = 2 * copt1 * copt15 * copt19 * copt21;
  out2(0, 0)  = copt36;
  out2(0, 1)  = copt39;
  out2(0, 2)  = copt42;
  out2(0, 3)  = copt46;
  out2(0, 4)  = copt49;
  out2(0, 5)  = copt52;
  out2(1, 0)  = copt39;
  out2(1, 1)  = copt55;
  out2(1, 2)  = copt58;
  out2(1, 3)  = copt49;
  out2(1, 4)  = copt61;
  out2(1, 5)  = copt64;
  out2(2, 0)  = copt42;
  out2(2, 1)  = copt58;
  out2(2, 2)  = copt67;
  out2(2, 3)  = copt52;
  out2(2, 4)  = copt64;
  out2(2, 5)  = copt70;
  out2(3, 0)  = copt46;
  out2(3, 1)  = copt49;
  out2(3, 2)  = copt52;
  out2(3, 3)  = copt36;
  out2(3, 4)  = copt39;
  out2(3, 5)  = copt42;
  out2(4, 0)  = copt49;
  out2(4, 1)  = copt61;
  out2(4, 2)  = copt64;
  out2(4, 3)  = copt39;
  out2(4, 4)  = copt55;
  out2(4, 5)  = copt58;
  out2(5, 0)  = copt52;
  out2(5, 1)  = copt64;
  out2(5, 2)  = copt70;
  out2(5, 3)  = copt42;
  out2(5, 4)  = copt58;
  out2(5, 5)  = copt67;

  return std::make_pair(out2, out1);
}

}  // namespace hylc