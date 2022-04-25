#include "hylc_conf.hpp"
#include <fstream>
#include "materials/SplineMaterial.hpp"

using namespace hylc;

hylc::Config hylc::config{};
hylc::Debug hylc::debug{};

extern void complain(const Json::Value &json, const std::string &expected);
extern void parse(double &x, const Json::Value &json);
extern void parse(int &x, const Json::Value &json);
extern void parse(bool &x, const Json::Value &json);
extern void parse(std::string &x, const Json::Value &json);

template <typename T>
void parse(T &x, const Json::Value &json, const T &x0) {
  if (json.isNull())
    x = x0;
  else
    parse(x, json);
}

// void parse(SplineMaterial::Spline1D &spline1d, const Json::Value &json);
// void parse(SplineMaterial::Spline2D &spline2d, const Json::Value &json);
// void parse(SplineMaterial::HSpline2Das1D &, const Json::Value &);

template <int n>
void parse(Vec<n> &v, const Json::Value &json) {
  if (!json.isArray())
    complain(json, "array");
  assert(json.size() == n);
  for (int i = 0; i < n; i++) parse(v[i], json[i]);
}
template <typename T>
void parse(std::vector<T> &v, const Json::Value &json) {
  if (!json.isArray())
    complain(json, "array");
  v.resize(json.size());
  for (int i = 0; i < json.size(); i++) parse(v[i], json[i]);
}

template <typename T>
void parse_optional(T &x, const Json::Value &json) {
  if (json.isNull())
    return;  // leave x as its default value
  else
    parse(x, json);
}

// void parse(SplineMaterial::Spline1D &spline1d, const Json::Value &json) {
//   std::vector<double> tx, c;
//   int degree = 3;

//   parse(spline1d.k, json["k"]);
//   parse(tx, json["tx"]);
//   parse(c, json["c"]);

//   spline1d.spline = std::make_shared<fitpackpp::BSplineCurve>(tx, c, degree);

//   parse(spline1d.clampx, json["clampx"]);
//   if (spline1d.clampx) {
//     parse(spline1d.xmin, json["xmin"]);
//     parse(spline1d.xmax, json["xmax"]);
//   }

// }

// void parse(SplineMaterial::Spline2D &spline2d, const Json::Value &json) {
//   std::vector<double> tx, ty, c;
//   int degree = 3;

//   parse(spline2d.k0, json["k0"]);
//   parse(spline2d.k1, json["k1"]);
//   parse(tx, json["tx"]);
//   parse(ty, json["ty"]);
//   parse(c, json["c"]);

//   spline2d.spline =
//       std::make_shared<fitpackpp::BSplineSurface>(tx, ty, c, degree);
// }

// void parse(SplineMaterial::HSpline2Das1D &spline, const Json::Value &json) {
//   parse(spline.k0, json["k0"]);
//   parse(spline.k1, json["k1"]);
//   parse(spline.fun.t, json["t"]);
//   parse(spline.fun.p, json["p"]);
//   parse(spline.fun.m, json["m"]);
//   parse(spline.fun.ext, json["ext"], 0);
// }

// void parse(Poly1D &poly, const Json::Value &json) {
//   parse(poly.c, json["k"]);

//   parse(poly.compr, json["compr"]);
//   parse(poly.clampx, json["clampx"]);
//   if (poly.clampx) {
//     parse(poly.xmin, json["xmin"]);
//     parse(poly.xmax, json["xmax"]);
//   }
// }

// void parse(Poly2D &poly, const Json::Value &json) {
//   parse(poly.k0, json["k0"]);
//   parse(poly.k1, json["k1"]);
//   parse(poly.c, json["C"]);

//   parse(poly.clampx, json["clampx"]);
//   if (poly.clampx) {
//     parse(poly.xmin, json["xmin"]);
//     parse(poly.xmax, json["xmax"]);
//   }
//   parse(poly.clampy, json["clampy"]);
//   if (poly.clampy) {
//     parse(poly.ymin, json["ymin"]);
//     parse(poly.ymax, json["ymax"]);
//   }
// }

void parse(HermiteSpline1D &spline, const Json::Value &json) {
  parse(spline.k, json["k"]);
  parse(spline.t, json["t"]);
  parse(spline.p, json["p"]);
  parse(spline.m, json["m"]);
  parse(spline.ext, json["ext"], 0);
}
void parse(HermiteSpline2D &spline, const Json::Value &json) {
  parse(spline.k0, json["k0"]);
  parse(spline.k1, json["k1"]);
  parse(spline.tu, json["tu"]);
  parse(spline.tv, json["tv"]);
  parse(spline.p, json["p"]);
  parse(spline.mu, json["mu"]);
  parse(spline.mv, json["mv"]);
  if (!json["muv"].isNull())
    parse(spline.muv, json["muv"]);
  else
    spline.muv.resize(spline.mu.size(), 0);
  parse(spline.ext, json["ext"], 0);
}

std::shared_ptr<SplineMaterial> load_material(const std::string &filename,
                                              bool only1D = false) {
  Json::Value json;
  Json::Reader reader;
  std::ifstream file(filename.c_str());
  bool parsingSuccessful = reader.parse(file, json);
  if (!parsingSuccessful) {
    fprintf(stderr, "Error reading file: %s\n", filename.c_str());
    fprintf(stderr, "%s", reader.getFormatedErrorMessages().c_str());
    abort();
  }
  file.close();

  auto material = std::make_shared<SplineMaterial>();

  parse(material->density, json["density"]);

  parse(material->strainscale, json["strain scale"]);

  Json::Value jsoncoeff = json["coeffs"];
  parse(material->C0, jsoncoeff["const"]);
  // parse(material->polys_1d, jsoncoeff["1D"]);
  parse(material->hsplines_1d, jsoncoeff["1D"]);
  // parse(material->splines_1d, jsoncoeff["1D"]);
  // parse(material->splines_2d, jsoncoeff["2D"]);
  // parse(material->polys_2d, jsoncoeff["2D"]);
  if (!only1D)
    parse(material->hsplines_2d, jsoncoeff["2D"]);
  // parse(material->hsplines_2d1d, jsoncoeff["2D"]);

  material->initialized = true;

  return material;
}

void parse(hylc::Config &params, const Json::Value &json) {
  if (json.isNull())
    return;
  parse_optional(params.enabled, json["enabled"]);

  std::string filename;
  parse(filename, json["material"]);
  parse(params.stiffness_mult, json["stiffness_mult"], 1.0);
  parse(params.bend_scale, json["bend_scale"], 1.0);
  parse(params.weight_mult, json["weight_mult"], 1.0);
  // parse(params.center_grav, json["center_grav"], 0.0);
  bool only1D;
  parse(only1D, json["only1D"], false);
  params.material = load_material(filename, only1D);
  parse(params.seam_stiffness, json["seam_stiffness"], 0.0);
}

void hylc::parse_hylc(const Json::Value &json) { parse(hylc::config, json); }
