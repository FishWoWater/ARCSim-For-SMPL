#include "HermiteSpline2D.hpp"

typedef HermiteSpline2D::Mat4x4 Mat4x4;
typedef HermiteSpline2D::Vec4 Vec4;

auto makeM = [](){
    // M = np.array([
    //         [2,-2,1,1],
    //         [-3,3,-2,-1],
    //         [0,0,1,0],
    //         [1,0,0,0]
    //     ])
    // return Mat4x4 {Vec4{0,0,-1,1},Vec4{0,0,1,0},Vec4{0,0,0,0},Vec4{0,0,0,0}};
    return Mat4x4 {Vec4{2,-3,0,1},Vec4{-2,3,0,0},Vec4{1,-2,1,0},Vec4{1,-1,0,0}};
};
auto makeMextL = [](){
    // [0,0,0,0],
    // [0,0,0,0],
    // [0,0,lin,0],
    // [1,0,0,0]
    return (Mat4x4 {Vec4{0,0,0,1},Vec4{0,0,0,0},Vec4{0,0,1,0},Vec4{0,0,0,0}});
};
auto makeMextR = [](){
    // [0,0,0,0],
    // [0,0,0,0],
    // [0,0,0,lin],
    // [0,1,0,0]
    return (Mat4x4 {Vec4{0,0,0,0},Vec4{0,0,0,1},Vec4{0,0,0,0},Vec4{0,0,1,0}});
};

Mat4x4 HermiteSpline2D::M = makeM();
Mat4x4 HermiteSpline2D::MT = makeM().t();
Mat4x4 HermiteSpline2D::MextL = makeMextL();
Mat4x4 HermiteSpline2D::MextLT = makeMextL().t();
Mat4x4 HermiteSpline2D::MextR = makeMextR();
Mat4x4 HermiteSpline2D::MextRT = makeMextR().t();