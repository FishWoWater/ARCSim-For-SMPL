#include "strain.hpp"
#ifndef hylc_strain_II

using namespace hylc;
using namespace hylc::mathematica;

std::tuple<Mat6x18, Vec6> hylc::mathematica::strain_valgrad(
    const Vec18 &xloc, const Mat2x2 &invDm, const Vec3 &niexists) {
  // define output
  Mat6x18 grad(0);
  Vec6 val(0);
  auto out1 = [&](int i) -> Real & { return val[i]; };
  auto out2 = [&](int i, int j) -> Real & { return grad(i, j); };

  Real copt2   = invDm(0, 0);
  Real copt3   = xloc(0);
  Real copt4   = -copt3;
  Real copt5   = xloc(3);
  Real copt6   = copt4 + copt5;
  Real copt8   = copt2 * copt6;
  Real copt9   = invDm(1, 0);
  Real copt10  = xloc(6);
  Real copt12  = copt10 + copt4;
  Real copt13  = copt12 * copt9;
  Real copt14  = copt13 + copt8;
  Real copt15  = Power(copt14, 2);
  Real copt16  = xloc(1);
  Real copt17  = -copt16;
  Real copt18  = xloc(4);
  Real copt19  = copt17 + copt18;
  Real copt20  = copt19 * copt2;
  Real copt22  = xloc(7);
  Real copt23  = copt17 + copt22;
  Real copt24  = copt23 * copt9;
  Real copt25  = copt20 + copt24;
  Real copt26  = Power(copt25, 2);
  Real copt27  = xloc(2);
  Real copt28  = -copt27;
  Real copt29  = xloc(5);
  Real copt30  = copt28 + copt29;
  Real copt32  = copt2 * copt30;
  Real copt33  = xloc(8);
  Real copt36  = copt28 + copt33;
  Real copt37  = copt36 * copt9;
  Real copt38  = copt32 + copt37;
  Real copt39  = Power(copt38, 2);
  Real copt43  = copt15 + copt26 + copt39;
  Real copt44  = Sqrt(copt43);
  Real copt48  = 1 / copt44;
  Real copt51  = invDm(0, 1);
  Real copt56  = invDm(1, 1);
  Real copt55  = copt51 * copt6;
  Real copt58  = copt12 * copt56;
  Real copt59  = copt55 + copt58;
  Real copt84  = Power(copt59, 2);
  Real copt67  = copt19 * copt51;
  Real copt69  = copt23 * copt56;
  Real copt71  = copt67 + copt69;
  Real copt85  = Power(copt71, 2);
  Real copt76  = copt30 * copt51;
  Real copt79  = copt36 * copt56;
  Real copt80  = copt76 + copt79;
  Real copt88  = Power(copt80, 2);
  Real copt89  = copt84 + copt85 + copt88;
  Real copt90  = Sqrt(copt89);
  Real copt91  = 1 / copt90;
  Real copt95  = -copt5;
  Real copt98  = copt3 + copt95;
  Real copt99  = Power(copt98, 2);
  Real copt100 = -copt18;
  Real copt101 = copt100 + copt16;
  Real copt102 = Power(copt101, 2);
  Real copt103 = -copt29;
  Real copt104 = copt103 + copt27;
  Real copt105 = Power(copt104, 2);
  Real copt106 = copt102 + copt105 + copt99;
  Real copt108 = Sqrt(copt106);
  Real copt109 = 1 / copt108;
  Real copt110 = -copt10;
  Real copt111 = copt110 + copt3;
  Real copt112 = Power(copt111, 2);
  Real copt113 = -copt22;
  Real copt114 = copt113 + copt16;
  Real copt115 = Power(copt114, 2);
  Real copt116 = -copt33;
  Real copt117 = copt116 + copt27;
  Real copt118 = Power(copt117, 2);
  Real copt119 = copt112 + copt115 + copt118;
  Real copt121 = Sqrt(copt119);
  Real copt122 = 1 / copt121;
  Real copt123 = copt110 + copt5;
  Real copt124 = Power(copt123, 2);
  Real copt125 = copt113 + copt18;
  Real copt126 = Power(copt125, 2);
  Real copt128 = copt116 + copt29;
  Real copt129 = Power(copt128, 2);
  Real copt131 = copt124 + copt126 + copt129;
  Real copt132 = Sqrt(copt131);
  Real copt133 = 1 / copt132;
  Real copt134 = Power(copt2, 2);
  Real copt135 = Power(copt3, 2);
  Real copt136 = Power(copt16, 2);
  Real copt137 = Power(copt27, 2);
  Real copt138 = -2 * copt3 * copt5;
  Real copt139 = Power(copt5, 2);
  Real copt140 = -2 * copt16 * copt18;
  Real copt141 = Power(copt18, 2);
  Real copt142 = -2 * copt27 * copt29;
  Real copt143 = Power(copt29, 2);
  Real copt144 = copt135 + copt136 + copt137 + copt138 + copt139 + copt140 +
                 copt141 + copt142 + copt143;
  Real copt145 = copt134 * copt144;
  Real copt146 = -(copt27 * copt29);
  Real copt147 = copt10 * copt5;
  Real copt148 = copt10 + copt5;
  Real copt149 = -(copt148 * copt3);
  Real copt150 = copt18 * copt22;
  Real copt151 = copt18 + copt22;
  Real copt153 = -(copt151 * copt16);
  Real copt154 = -(copt27 * copt33);
  Real copt155 = copt29 * copt33;
  Real copt156 = copt135 + copt136 + copt137 + copt146 + copt147 + copt149 +
                 copt150 + copt153 + copt154 + copt155;
  Real copt157 = 2 * copt156 * copt2 * copt9;
  Real copt158 = Power(copt9, 2);
  Real copt159 = -2 * copt10 * copt3;
  Real copt160 = Power(copt10, 2);
  Real copt161 = -2 * copt16 * copt22;
  Real copt162 = Power(copt22, 2);
  Real copt163 = -2 * copt27 * copt33;
  Real copt164 = Power(copt33, 2);
  Real copt165 = copt135 + copt136 + copt137 + copt159 + copt160 + copt161 +
                 copt162 + copt163 + copt164;
  Real copt166 = copt158 * copt165;
  Real copt167 = copt145 + copt157 + copt166;
  Real copt168 = 1 / copt167;
  Real copt169 = copt123 * copt16;
  Real copt170 = copt10 * copt18;
  Real copt171 = -(copt22 * copt5);
  Real copt172 = copt100 + copt22;
  Real copt173 = copt172 * copt3;
  Real copt174 = copt169 + copt170 + copt171 + copt173;
  Real copt175 = Power(copt174, 2);
  Real copt176 = copt123 * copt27;
  Real copt177 = copt10 * copt29;
  Real copt178 = -(copt33 * copt5);
  Real copt179 = copt103 + copt33;
  Real copt180 = copt179 * copt3;
  Real copt181 = copt176 + copt177 + copt178 + copt180;
  Real copt182 = Power(copt181, 2);
  Real copt183 = copt125 * copt27;
  Real copt184 = copt22 * copt29;
  Real copt185 = -(copt18 * copt33);
  Real copt186 = copt16 * copt179;
  Real copt187 = copt183 + copt184 + copt185 + copt186;
  Real copt188 = Power(copt187, 2);
  Real copt190 = copt175 + copt182 + copt188;
  Real copt191 = Sqrt(copt190);
  Real copt192 = xloc(12);
  Real copt193 = -copt192;
  Real copt194 = copt193 + copt5;
  Real copt197 = copt110 + copt192;
  Real copt199 = copt10 + copt95;
  Real copt200 = xloc(13);
  Real copt214 = xloc(14);
  Real copt196 = copt194 * copt22;
  Real copt198 = copt18 * copt197;
  Real copt208 = copt199 * copt200;
  Real copt209 = copt196 + copt198 + copt208;
  Real copt252 = copt2 + copt9;
  Real copt253 = Power(copt252, 2);
  Real copt255 = Sqrt(copt144);
  Real copt256 = Sqrt(copt165);
  Real copt258 = -2 * copt10 * copt5;
  Real copt259 = -2 * copt18 * copt22;
  Real copt260 = -2 * copt29 * copt33;
  Real copt261 = copt139 + copt141 + copt143 + copt160 + copt162 + copt164 +
                 copt258 + copt259 + copt260;
  Real copt262 = Sqrt(copt261);
  Real copt263 = xloc(15);
  Real copt264 = -copt263;
  Real copt265 = copt10 + copt264;
  Real copt267 = copt263 + copt4;
  Real copt272 = xloc(16);
  Real copt278 = xloc(17);
  Real copt310 = -(copt10 * copt18);
  Real copt312 = copt125 * copt3;
  Real copt313 = copt22 * copt5;
  Real copt338 = xloc(9);
  Real copt339 = -copt338;
  Real copt345 = xloc(10);
  Real copt367 = xloc(11);
  Real copt227 = -(copt22 * copt29 * copt3);
  Real copt228 = copt18 * copt3 * copt33;
  Real copt340 = copt339 + copt5;
  Real copt237 = copt170 + copt171 + copt173;
  Real copt353 = copt338 + copt95;
  Real copt721 = copt190 * copt191;
  Real copt722 = copt190 * copt43;
  Real copt723 = Sqrt(copt722);
  Real copt724 = copt722 * copt723;
  Real copt725 = 1 / copt724;
  Real copt210 = copt174 * copt209;
  Real copt211 = copt194 * copt33;
  Real copt213 = copt197 * copt29;
  Real copt215 = copt199 * copt214;
  Real copt216 = copt211 + copt213 + copt215;
  Real copt217 = copt181 * copt216;
  Real copt218 = -copt200;
  Real copt219 = copt18 + copt218;
  Real copt220 = copt219 * copt33;
  Real copt221 = copt113 + copt200;
  Real copt222 = copt221 * copt29;
  Real copt223 = copt172 * copt214;
  Real copt224 = copt220 + copt222 + copt223;
  Real copt225 = copt187 * copt224;
  Real copt226 = copt210 + copt217 + copt225;
  Real copt229 = copt192 * copt22 * copt29;
  Real copt230 = -(copt18 * copt192 * copt33);
  Real copt231 = copt200 * copt29 * copt3;
  Real copt232 = -(copt10 * copt200 * copt29);
  Real copt233 = -(copt200 * copt3 * copt33);
  Real copt234 = copt200 * copt33 * copt5;
  Real copt235 = copt209 * copt27;
  Real copt238 = copt214 * copt237;
  Real copt239 = copt10 + copt193;
  Real copt240 = copt239 * copt29;
  Real copt242 = copt192 + copt95;
  Real copt243 = copt242 * copt33;
  Real copt245 = copt123 * copt214;
  Real copt246 = copt240 + copt243 + copt245;
  Real copt247 = copt16 * copt246;
  Real copt248 = copt227 + copt228 + copt229 + copt230 + copt231 + copt232 +
                 copt233 + copt234 + copt235 + copt238 + copt247;
  Real copt249 = copt132 * copt248;
  Real copt250 = ArcTan(copt226, copt249);
  Real copt254 = niexists(1);
  Real copt727 = -(copt10 * copt5);
  Real copt731 = -(copt18 * copt22);
  Real copt733 = -(copt29 * copt33);
  Real copt328 = copt16 * copt199;
  Real copt329 = copt310 + copt312 + copt313 + copt328;
  Real copt342 = copt16 * copt340;
  Real copt343 = copt338 + copt4;
  Real copt344 = copt18 * copt343;
  Real copt346 = copt345 * copt98;
  Real copt347 = copt342 + copt344 + copt346;
  Real copt348 = copt329 * copt347;
  Real copt351 = copt3 + copt339;
  Real copt352 = copt29 * copt351;
  Real copt360 = copt27 * copt353;
  Real copt370 = copt367 * copt6;
  Real copt371 = copt352 + copt360 + copt370;
  Real copt372 = copt181 * copt371;
  Real copt375 = -(copt22 * copt29);
  Real copt376 = copt172 * copt27;
  Real copt378 = copt128 * copt16;
  Real copt379 = copt18 * copt33;
  Real copt380 = copt375 + copt376 + copt378 + copt379;
  Real copt383 = -copt345;
  Real copt384 = copt18 + copt383;
  Real copt385 = copt27 * copt384;
  Real copt389 = copt17 + copt345;
  Real copt390 = copt29 * copt389;
  Real copt391 = copt101 * copt367;
  Real copt401 = copt385 + copt390 + copt391;
  Real copt402 = copt380 * copt401;
  Real copt405 = copt348 + copt372 + copt402;
  Real copt406 = copt22 * copt29 * copt338;
  Real copt411 = -(copt18 * copt33 * copt338);
  Real copt543 = copt29 * copt3 * copt345;
  Real copt562 = -(copt10 * copt29 * copt345);
  Real copt570 = -(copt3 * copt33 * copt345);
  Real copt604 = copt33 * copt345 * copt5;
  Real copt627 = copt22 * copt340;
  Real copt634 = copt110 + copt338;
  Real copt662 = copt18 * copt634;
  Real copt663 = copt199 * copt345;
  Real copt664 = copt627 + copt662 + copt663;
  Real copt665 = copt27 * copt664;
  Real copt666 = copt237 * copt367;
  Real copt680 = copt10 + copt339;
  Real copt685 = copt29 * copt680;
  Real copt686 = copt33 * copt353;
  Real copt687 = copt123 * copt367;
  Real copt695 = copt685 + copt686 + copt687;
  Real copt696 = copt16 * copt695;
  Real copt697 = copt227 + copt228 + copt406 + copt411 + copt543 + copt562 +
                 copt570 + copt604 + copt665 + copt666 + copt696;
  Real copt698 = copt108 * copt697;
  Real copt699 = ArcTan(copt405, copt698);
  Real copt708 = niexists(0);
  Real copt266 = copt16 * copt265;
  Real copt271 = copt22 * copt267;
  Real copt273 = copt111 * copt272;
  Real copt274 = copt266 + copt271 + copt273;
  Real copt275 = copt174 * copt274;
  Real copt276 = copt265 * copt27;
  Real copt277 = copt267 * copt33;
  Real copt279 = copt111 * copt278;
  Real copt280 = copt276 + copt277 + copt279;
  Real copt281 = copt181 * copt280;
  Real copt282 = -copt272;
  Real copt283 = copt22 + copt282;
  Real copt284 = copt27 * copt283;
  Real copt285 = copt17 + copt272;
  Real copt286 = copt285 * copt33;
  Real copt287 = copt114 * copt278;
  Real copt288 = copt284 + copt286 + copt287;
  Real copt289 = copt187 * copt288;
  Real copt290 = copt275 + copt281 + copt289;
  Real copt291 = copt22 * copt29 * copt3;
  Real copt292 = -(copt18 * copt3 * copt33);
  Real copt293 = -(copt22 * copt263 * copt29);
  Real copt294 = copt18 * copt263 * copt33;
  Real copt295 = -(copt272 * copt29 * copt3);
  Real copt296 = copt10 * copt272 * copt29;
  Real copt300 = copt272 * copt3 * copt33;
  Real copt301 = -(copt272 * copt33 * copt5);
  Real copt302 = copt18 * copt265;
  Real copt303 = copt263 + copt95;
  Real copt304 = copt22 * copt303;
  Real copt305 = copt123 * copt272;
  Real copt308 = copt302 + copt304 + copt305;
  Real copt309 = copt27 * copt308;
  Real copt314 = copt310 + copt312 + copt313;
  Real copt315 = copt278 * copt314;
  Real copt316 = copt264 + copt5;
  Real copt317 = copt316 * copt33;
  Real copt318 = copt110 + copt263;
  Real copt319 = copt29 * copt318;
  Real copt320 = copt199 * copt278;
  Real copt321 = copt317 + copt319 + copt320;
  Real copt322 = copt16 * copt321;
  Real copt323 = copt291 + copt292 + copt293 + copt294 + copt295 + copt296 +
                 copt300 + copt301 + copt309 + copt315 + copt322;
  Real copt324 = -(copt121 * copt323);
  Real copt325 = ArcTan(copt290, copt324);
  Real copt326 = niexists(2);
  Real copt771 = copt106 * copt134;
  Real copt772 = copt119 * copt158;
  Real copt773 = copt157 + copt771 + copt772;
  Real copt774 = 1 / copt773;
  Real copt775 = 1 / copt191;
  Real copt776 = copt119 * copt9;
  Real copt753 = copt156 * copt2;
  Real copt777 = copt753 + copt776;
  Real copt778 = Power(copt777, 2);
  Real copt780 = copt106 * copt2;
  Real copt749 = copt156 * copt9;
  Real copt781 = copt749 + copt780;
  Real copt782 = Power(copt781, 2);
  Real copt726 = -(copt16 * copt18);
  Real copt728 = copt199 * copt3;
  Real copt729 = copt16 * copt22;
  Real copt732 = copt27 * copt33;
  Real copt734 = copt139 + copt141 + copt143 + copt146 + copt726 + copt727 +
                 copt728 + copt729 + copt731 + copt732 + copt733;
  Real copt784 = copt2 * copt734;
  Real copt736 = copt27 * copt29;
  Real copt737 = copt123 * copt3;
  Real copt743 = copt125 * copt16;
  Real copt744 = copt154 + copt160 + copt162 + copt164 + copt727 + copt731 +
                 copt733 + copt736 + copt737 + copt743;
  Real copt785  = -(copt744 * copt9);
  Real copt786  = copt784 + copt785;
  Real copt787  = Power(copt786, 2);
  Real copt860  = copt43 * copt44;
  Real copt861  = 1 / copt860;
  Real copt862  = copt89 * copt90;
  Real copt863  = 1 / copt862;
  Real copt66   = copt14 * copt59;
  Real copt74   = copt25 * copt71;
  Real copt82   = copt38 * copt80;
  Real copt83   = copt66 + copt74 + copt82;
  Real copt864  = copt51 + copt56;
  Real copt796  = copt2 * copt98;
  Real copt799  = copt111 * copt9;
  Real copt800  = copt796 + copt799;
  Real copt818  = copt101 * copt2;
  Real copt820  = copt114 * copt9;
  Real copt821  = copt818 + copt820;
  Real copt841  = copt104 * copt2;
  Real copt845  = copt117 * copt9;
  Real copt849  = copt841 + copt845;
  Real copt902  = copt51 * copt98;
  Real copt903  = copt111 * copt56;
  Real copt912  = copt902 + copt903;
  Real copt952  = copt101 * copt51;
  Real copt954  = copt114 * copt56;
  Real copt956  = copt952 + copt954;
  Real copt974  = copt104 * copt51;
  Real copt975  = copt117 * copt56;
  Real copt976  = copt974 + copt975;
  Real copt257  = copt250 * copt253 * copt254 * copt255 * copt256;
  Real copt327  = copt134 * copt255 * copt325 * copt326;
  Real copt712  = copt158 * copt256 * copt699 * copt708;
  Real copt714  = copt327 + copt712;
  Real copt715  = copt262 * copt714;
  Real copt717  = copt257 + copt715;
  Real copt1134 = Power(copt167, 2);
  Real copt1135 = 1 / copt1134;
  Real copt1143 = copt119 * copt121;
  Real copt1144 = 1 / copt1143;
  Real copt1147 = copt106 * copt108;
  Real copt1148 = 1 / copt1147;
  Real copt1152 = 1 / copt256;
  Real copt1154 = 1 / copt255;
  Real copt1158 = copt135 * copt141;
  Real copt1159 = copt135 * copt143;
  Real copt1160 = -2 * copt10 * copt141 * copt3;
  Real copt1161 = -2 * copt10 * copt143 * copt3;
  Real copt1162 = copt141 * copt160;
  Real copt1163 = copt143 * copt160;
  Real copt1164 = -2 * copt135 * copt18 * copt22;
  Real copt1165 = 2 * copt18 * copt22 * copt3 * copt5;
  Real copt1166 = 2 * copt10 * copt18 * copt22 * copt3;
  Real copt1167 = -2 * copt10 * copt18 * copt22 * copt5;
  Real copt1168 = copt135 * copt162;
  Real copt1170 = -2 * copt162 * copt3 * copt5;
  Real copt1171 = copt139 * copt162;
  Real copt1172 = copt143 * copt162;
  Real copt1173 = copt139 + copt141 + copt160 + copt162 + copt258 + copt259;
  Real copt1174 = copt1173 * copt137;
  Real copt1175 = -2 * copt135 * copt29 * copt33;
  Real copt1176 = 2 * copt29 * copt3 * copt33 * copt5;
  Real copt1177 = 2 * copt10 * copt29 * copt3 * copt33;
  Real copt1178 = -2 * copt10 * copt29 * copt33 * copt5;
  Real copt1179 = -2 * copt18 * copt22 * copt29 * copt33;
  Real copt1180 = copt135 * copt164;
  Real copt1182 = -2 * copt164 * copt3 * copt5;
  Real copt1183 = copt139 * copt164;
  Real copt1184 = copt141 * copt164;
  Real copt1185 = copt139 + copt143 + copt160 + copt164 + copt258 + copt260;
  Real copt1186 = copt1185 * copt136;
  Real copt1187 = -(copt10 * copt18 * copt5);
  Real copt1188 = copt160 * copt18;
  Real copt1189 = copt123 * copt125 * copt3;
  Real copt1190 = copt139 * copt22;
  Real copt1191 = copt143 * copt22;
  Real copt1192 = -(copt10 * copt22 * copt5);
  Real copt1193 = copt125 * copt128 * copt27;
  Real copt1200 = -(copt18 * copt29 * copt33);
  Real copt1201 = -(copt22 * copt29 * copt33);
  Real copt1202 = copt164 * copt18;
  Real copt1210 = copt1187 + copt1188 + copt1189 + copt1190 + copt1191 +
                  copt1192 + copt1193 + copt1200 + copt1201 + copt1202;
  Real copt1211 = -2 * copt1210 * copt16;
  Real copt1212 = copt160 * copt29;
  Real copt1213 = -(copt18 * copt22 * copt29);
  Real copt1214 = copt162 * copt29;
  Real copt1223 = copt123 * copt128 * copt3;
  Real copt1224 = copt139 * copt33;
  Real copt1225 = copt141 * copt33;
  Real copt1226 = -(copt18 * copt22 * copt33);
  Real copt1227 = copt29 + copt33;
  Real copt1228 = -(copt10 * copt1227 * copt5);
  Real copt1229 = copt1212 + copt1213 + copt1214 + copt1223 + copt1224 +
                  copt1225 + copt1226 + copt1228;
  Real copt1230 = -2 * copt1229 * copt27;
  Real copt1231 = copt1158 + copt1159 + copt1160 + copt1161 + copt1162 +
                  copt1163 + copt1164 + copt1165 + copt1166 + copt1167 +
                  copt1168 + copt1170 + copt1171 + copt1172 + copt1174 +
                  copt1175 + copt1176 + copt1177 + copt1178 + copt1179 +
                  copt1180 + copt1182 + copt1183 + copt1184 + copt1186 +
                  copt1211 + copt1230;
  Real copt1232 = 1 / copt1231;
  Real copt1239 = Power(copt405, 2);
  Real copt1240 = Power(copt697, 2);
  Real copt1241 = copt106 * copt1240;
  Real copt1242 = copt1239 + copt1241;
  Real copt1243 = 1 / copt1242;
  Real copt1247 = -copt367;
  Real copt1271 = Power(copt323, 2);
  Real copt1272 = copt119 * copt1271;
  Real copt1273 = Power(copt290, 2);
  Real copt1274 = copt1272 + copt1273;
  Real copt1276 = 1 / copt1274;
  Real copt1281 = copt116 + copt278;
  Real copt1467 = copt131 * copt132;
  Real copt1468 = 1 / copt1467;
  Real copt1472 = 1 / copt262;
  Real copt1506 = Power(copt248, 2);
  Real copt1509 = copt131 * copt1506;
  Real copt1512 = Power(copt226, 2);
  Real copt1513 = copt1509 + copt1512;
  Real copt1514 = 1 / copt1513;
  Real copt1515 = copt218 + copt22;
  Real copt1560 = -2 * copt16;
  Real copt1561 = 2 * copt18;
  Real copt1563 = copt1560 + copt1561;
  Real copt1518 = -copt214;
  Real copt1519 = copt1518 + copt33;
  Real copt1586 = copt3 * copt33;
  Real copt1691 = -2 * copt27;
  Real copt1694 = 2 * copt29;
  Real copt1722 = copt1691 + copt1694;
  Real copt1819 = -(copt22 * copt3);
  Real copt2115 = -2 * copt3;
  Real copt2122 = 2 * copt10;
  Real copt2123 = copt2115 + copt2122;
  Real copt2162 = copt100 + copt200;
  Real copt1244 = copt100 + copt345;
  Real copt1248 = copt1247 + copt29;
  Real copt2336 = 2 * copt22;
  Real copt2337 = copt1560 + copt2336;
  Real copt2165 = copt103 + copt214;
  Real copt2374 = -(copt29 * copt3);
  Real copt2251 = -copt278;
  Real copt2268 = copt2251 + copt27;
  Real copt2486 = 2 * copt33;
  Real copt2487 = copt1691 + copt2486;
  Real copt2519 = copt18 * copt3;
  Real copt2445 = copt29 * copt3;
  Real copt1330 = -(copt10 * copt29);
  Real copt1331 = copt199 * copt27;
  Real copt1645 = -(copt3 * copt33);
  Real copt1335 = copt33 * copt5;
  Real copt2684 = copt1330 + copt1331 + copt1335 + copt1645 + copt2445;
  Real copt1127 = 2 * copt172 * copt174;
  Real copt1131 = 2 * copt179 * copt181;
  Real copt1132 = copt1127 + copt1131;
  Real copt735  = -(copt2 * copt734);
  Real copt745  = copt744 * copt9;
  Real copt746  = copt735 + copt745;
  Real copt747  = copt250 * copt252 * copt254 * copt255 * copt256 * copt746;
  Real copt748  = copt144 * copt2;
  Real copt750  = copt748 + copt749;
  Real copt752  = -(copt256 * copt699 * copt708 * copt750 * copt9);
  Real copt754  = copt165 * copt9;
  Real copt756  = copt753 + copt754;
  Real copt757  = copt2 * copt255 * copt325 * copt326 * copt756;
  Real copt758  = copt752 + copt757;
  Real copt759  = copt262 * copt758;
  Real copt760  = copt747 + copt759;
  Real copt2807 = Power(copt722, 2);
  Real copt2808 = copt2807 * copt723;
  Real copt2809 = 1 / copt2808;
  Real copt2800 = -copt2;
  Real copt2802 = -copt9;
  Real copt2803 = copt2800 + copt2802;
  Real copt2821 = 2 * copt3;
  Real copt2159 = -2 * copt5;
  Real copt2822 = -2 * copt10;
  Real copt2826 = copt2821 + copt2822;
  Real copt2828 = copt2159 + copt2821;
  Real copt2832 = copt110 + copt2821 + copt95;
  Real copt1245 = copt1244 * copt329;
  Real copt1246 = copt125 * copt347;
  Real copt1249 = copt1248 * copt181;
  Real copt1250 = copt179 * copt371;
  Real copt1251 = copt1245 + copt1246 + copt1249 + copt1250;
  Real copt1252 = -(copt108 * copt1251 * copt697);
  Real copt1253 = -(copt33 * copt345);
  Real copt1254 = copt113 + copt345;
  Real copt1255 = copt1254 * copt29;
  Real copt1256 = copt1247 + copt33;
  Real copt1257 = copt1256 * copt18;
  Real copt1258 = copt22 * copt367;
  Real copt1259 = copt1253 + copt1255 + copt1257 + copt1258;
  Real copt1260 = copt106 * copt1259;
  Real copt1264 = copt697 * copt98;
  Real copt1265 = copt1260 + copt1264;
  Real copt1266 = copt109 * copt1265 * copt405;
  Real copt1269 = copt1252 + copt1266;
  Real copt1277 = copt113 + copt272;
  Real copt1279 = copt1277 * copt174;
  Real copt1280 = copt172 * copt274;
  Real copt1282 = copt1281 * copt181;
  Real copt1289 = copt179 * copt280;
  Real copt1290 = copt1279 + copt1280 + copt1282 + copt1289;
  Real copt1291 = copt119 * copt1290 * copt323;
  Real copt1293 = copt283 * copt29;
  Real copt1294 = copt272 * copt33;
  Real copt1295 = -(copt22 * copt278);
  Real copt1297 = copt1281 * copt18;
  Real copt1300 = copt1293 + copt1294 + copt1295 + copt1297;
  Real copt1307 = copt119 * copt1300;
  Real copt1308 = copt111 * copt323;
  Real copt1309 = copt1307 + copt1308;
  Real copt1310 = -(copt1309 * copt290);
  Real copt1311 = copt1291 + copt1310;
  Real copt1319 = 2 * copt123 * copt174;
  Real copt1320 = 2 * copt179 * copt187;
  Real copt1321 = copt1319 + copt1320;
  Real copt1576 = -2 * copt22;
  Real copt2861 = 2 * copt16;
  Real copt2344 = -2 * copt18;
  Real copt1334 = copt128 * copt3;
  Real copt1336 = copt1330 + copt1331 + copt1334 + copt1335;
  Real copt2862 = copt1576 + copt2861;
  Real copt2864 = copt2344 + copt2861;
  Real copt2868 = copt100 + copt113 + copt2861;
  Real copt1340 = copt329 * copt340;
  Real copt1341 = copt199 * copt347;
  Real copt1342 = copt103 + copt367;
  Real copt1343 = copt1342 * copt380;
  Real copt1344 = copt128 * copt401;
  Real copt1345 = copt1340 + copt1341 + copt1343 + copt1344;
  Real copt1347 = -(copt108 * copt1345 * copt697);
  Real copt1350 = copt106 * copt695;
  Real copt1351 = copt101 * copt697;
  Real copt1352 = copt1350 + copt1351;
  Real copt1353 = copt109 * copt1352 * copt405;
  Real copt1354 = copt1347 + copt1353;
  Real copt1361 = copt174 * copt265;
  Real copt1363 = copt123 * copt274;
  Real copt1364 = copt1281 * copt187;
  Real copt1372 = copt179 * copt288;
  Real copt1374 = copt1361 + copt1363 + copt1364 + copt1372;
  Real copt1375 = copt119 * copt1374 * copt323;
  Real copt1376 = copt119 * copt321;
  Real copt1377 = copt114 * copt323;
  Real copt1378 = copt1376 + copt1377;
  Real copt1379 = -(copt1378 * copt290);
  Real copt1380 = copt1375 + copt1379;
  Real copt1387 = 2 * copt123 * copt181;
  Real copt1388 = 2 * copt125 * copt187;
  Real copt1389 = copt1387 + copt1388;
  Real copt1749 = -2 * copt33;
  Real copt2905 = 2 * copt27;
  Real copt2510 = -2 * copt29;
  Real copt2906 = copt1749 + copt2905;
  Real copt2908 = copt2510 + copt2905;
  Real copt2912 = copt103 + copt116 + copt2905;
  Real copt1403 = copt181 * copt353;
  Real copt1404 = copt380 * copt384;
  Real copt1405 = copt123 * copt371;
  Real copt1407 = copt172 * copt401;
  Real copt1408 = copt1403 + copt1404 + copt1405 + copt1407;
  Real copt1409 = -(copt108 * copt1408 * copt697);
  Real copt1410 = copt106 * copt664;
  Real copt1411 = copt104 * copt697;
  Real copt1412 = copt1410 + copt1411;
  Real copt1413 = copt109 * copt1412 * copt405;
  Real copt1414 = copt1409 + copt1413;
  Real copt1416 = copt181 * copt265;
  Real copt1417 = copt187 * copt283;
  Real copt1419 = copt123 * copt280;
  Real copt1420 = copt125 * copt288;
  Real copt1421 = copt1416 + copt1417 + copt1419 + copt1420;
  Real copt1422 = copt119 * copt1421 * copt323;
  Real copt1423 = copt119 * copt308;
  Real copt1424 = copt117 * copt323;
  Real copt1425 = copt1423 + copt1424;
  Real copt1426 = -(copt1425 * copt290);
  Real copt1433 = copt1422 + copt1426;
  Real copt1454 = 2 * copt114 * copt174;
  Real copt1456 = 2 * copt117 * copt181;
  Real copt1457 = copt1454 + copt1456;
  Real copt2944 = 2 * copt5;
  Real copt2949 = copt2115 + copt2944;
  Real copt1477 = copt16 + copt383;
  Real copt1478 = copt1477 * copt329;
  Real copt1479 = copt23 * copt347;
  Real copt1480 = copt28 + copt367;
  Real copt1481 = copt1480 * copt181;
  Real copt1482 = copt117 * copt371;
  Real copt1483 = copt1478 + copt1479 + copt1481 + copt1482;
  Real copt1484 = -(copt108 * copt1483 * copt697);
  Real copt1485 = copt22 + copt383;
  Real copt1486 = copt1485 * copt27;
  Real copt1487 = copt33 * copt345;
  Real copt1490 = -(copt22 * copt367);
  Real copt1491 = copt116 + copt367;
  Real copt1492 = copt1491 * copt16;
  Real copt1493 = copt1486 + copt1487 + copt1490 + copt1492;
  Real copt1495 = copt106 * copt1493;
  Real copt1496 = -(copt697 * copt98);
  Real copt1499 = copt1495 + copt1496;
  Real copt1500 = copt109 * copt1499 * copt405;
  Real copt1501 = copt1484 + copt1500;
  Real copt1516 = copt1515 * copt174;
  Real copt1517 = copt114 * copt209;
  Real copt1520 = copt1519 * copt181;
  Real copt1521 = copt117 * copt216;
  Real copt1522 = copt1516 + copt1517 + copt1520 + copt1521;
  Real copt1523 = -(copt132 * copt1522 * copt248);
  Real copt1524 = copt1515 * copt27;
  Real copt1525 = copt200 * copt33;
  Real copt1526 = -(copt214 * copt22);
  Real copt1535 = copt116 + copt214;
  Real copt1536 = copt1535 * copt16;
  Real copt1537 = copt1524 + copt1525 + copt1526 + copt1536;
  Real copt1539 = copt131 * copt1537;
  Real copt1540 = copt123 * copt248;
  Real copt1541 = copt1539 + copt1540;
  Real copt1543 = copt133 * copt1541 * copt226;
  Real copt1550 = copt1523 + copt1543;
  Real copt1556 = 2 * copt12 * copt174;
  Real copt1557 = 2 * copt117 * copt187;
  Real copt1558 = copt1556 + copt1557;
  Real copt1577 = copt1561 + copt1576;
  Real copt1579 = copt174 * copt197;
  Real copt1580 = copt12 * copt209;
  Real copt1581 = copt1519 * copt187;
  Real copt1582 = copt117 * copt224;
  Real copt1583 = copt1579 + copt1580 + copt1581 + copt1582;
  Real copt1584 = -(copt132 * copt1514 * copt1583 * copt248);
  Real copt1587 = -(copt192 * copt33);
  Real copt1588 = copt197 * copt27;
  Real copt1589 = copt12 * copt214;
  Real copt1590 = copt1586 + copt1587 + copt1588 + copt1589;
  Real copt1591 = copt132 * copt1590;
  Real copt1594 = copt125 * copt133 * copt248;
  Real copt1595 = copt1591 + copt1594;
  Real copt1596 = copt1514 * copt1595 * copt226;
  Real copt1597 = copt1584 + copt1596;
  Real copt1601 = copt329 * copt343;
  Real copt1602 = copt111 * copt347;
  Real copt1603 = copt1247 + copt27;
  Real copt1604 = copt1603 * copt380;
  Real copt1605 = copt36 * copt401;
  Real copt1607 = copt1601 + copt1602 + copt1604 + copt1605;
  Real copt1615 = -(copt108 * copt1243 * copt1607 * copt697);
  Real copt1616 = -(copt33 * copt338);
  Real copt1617 = copt27 * copt634;
  Real copt1619 = copt12 * copt367;
  Real copt1621 = copt1586 + copt1616 + copt1617 + copt1619;
  Real copt1622 = copt108 * copt1621;
  Real copt1630 = -(copt101 * copt109 * copt697);
  Real copt1631 = copt1622 + copt1630;
  Real copt1633 = copt1243 * copt1631 * copt405;
  Real copt1634 = copt1615 + copt1633;
  Real copt1638 = copt12 * copt274;
  Real copt1639 = copt117 * copt288;
  Real copt1643 = copt1638 + copt1639;
  Real copt1644 = copt121 * copt1276 * copt1643 * copt323;
  Real copt1658 = copt263 * copt33;
  Real copt1662 = copt1645 + copt1658 + copt276 + copt279;
  Real copt1663 = -(copt121 * copt1276 * copt1662 * copt290);
  Real copt1664 = copt1644 + copt1663;
  Real copt1680 = 2 * copt12 * copt181;
  Real copt1688 = 2 * copt187 * copt23;
  Real copt1689 = copt1680 + copt1688;
  Real copt1750 = copt1694 + copt1749;
  Real copt1758 = copt181 * copt197;
  Real copt1759 = copt187 * copt221;
  Real copt1773 = copt12 * copt216;
  Real copt1800 = copt224 * copt23;
  Real copt1817 = copt1758 + copt1759 + copt1773 + copt1800;
  Real copt1818 = -(copt132 * copt1514 * copt1817 * copt248);
  Real copt1820 = copt16 * copt239;
  Real copt1822 = copt192 * copt22;
  Real copt1825 = copt200 * copt3;
  Real copt1830 = -(copt10 * copt200);
  Real copt1831 = copt1819 + copt1820 + copt1822 + copt1825 + copt1830;
  Real copt1832 = copt132 * copt1831;
  Real copt1834 = copt128 * copt133 * copt248;
  Real copt1835 = copt1832 + copt1834;
  Real copt1836 = copt1514 * copt1835 * copt226;
  Real copt1856 = copt1818 + copt1836;
  Real copt1887 = copt181 * copt351;
  Real copt1891 = copt380 * copt389;
  Real copt1898 = copt12 * copt371;
  Real copt1899 = copt114 * copt401;
  Real copt1900 = copt1887 + copt1891 + copt1898 + copt1899;
  Real copt1906 = -(copt108 * copt1243 * copt1900 * copt697);
  Real copt1907 = copt16 * copt680;
  Real copt1911 = copt22 * copt338;
  Real copt1912 = copt3 * copt345;
  Real copt1917 = -(copt10 * copt345);
  Real copt1918 = copt1819 + copt1907 + copt1911 + copt1912 + copt1917;
  Real copt1920 = copt108 * copt1918;
  Real copt1921 = -(copt104 * copt109 * copt697);
  Real copt1922 = copt1920 + copt1921;
  Real copt1959 = copt1243 * copt1922 * copt405;
  Real copt1960 = copt1906 + copt1959;
  Real copt1970 = copt12 * copt280;
  Real copt1971 = copt23 * copt288;
  Real copt1972 = copt1970 + copt1971;
  Real copt1982 = copt121 * copt1276 * copt1972 * copt323;
  Real copt1983 = copt22 * copt3;
  Real copt1987 = -(copt22 * copt263);
  Real copt1993 = copt16 * copt318;
  Real copt1994 = -(copt272 * copt3);
  Real copt1995 = copt10 * copt272;
  Real copt2051 = copt1983 + copt1987 + copt1993 + copt1994 + copt1995;
  Real copt2052 = -(copt121 * copt1276 * copt2051 * copt290);
  Real copt2059 = copt1982 + copt2052;
  Real copt2100 = 2 * copt174 * copt19;
  Real copt2104 = 2 * copt181 * copt30;
  Real copt2105 = copt2100 + copt2104;
  Real copt2160 = copt2122 + copt2159;
  Real copt2163 = copt174 * copt2162;
  Real copt2164 = copt19 * copt209;
  Real copt2166 = copt181 * copt2165;
  Real copt2167 = copt216 * copt30;
  Real copt2168 = copt2163 + copt2164 + copt2166 + copt2167;
  Real copt2169 = -(copt132 * copt1514 * copt2168 * copt248);
  Real copt2173 = -(copt200 * copt29);
  Real copt2183 = copt2162 * copt27;
  Real copt2195 = copt1518 + copt29;
  Real copt2196 = copt16 * copt2195;
  Real copt2197 = copt18 * copt214;
  Real copt2198 = copt2173 + copt2183 + copt2196 + copt2197;
  Real copt2199 = copt132 * copt2198;
  Real copt2200 = -(copt123 * copt133 * copt248);
  Real copt2201 = copt2199 + copt2200;
  Real copt2202 = copt1514 * copt2201 * copt226;
  Real copt2203 = copt2169 + copt2202;
  Real copt2209 = -(copt29 * copt345);
  Real copt2215 = copt1244 * copt27;
  Real copt2222 = copt1248 * copt16;
  Real copt2235 = copt18 * copt367;
  Real copt2236 = copt2209 + copt2215 + copt2222 + copt2235;
  Real copt2237 = copt108 * copt1243 * copt2236 * copt405;
  Real copt2238 = copt101 * copt347;
  Real copt2239 = copt30 * copt371;
  Real copt2240 = copt2238 + copt2239;
  Real copt2241 = -(copt108 * copt1243 * copt2240 * copt697);
  Real copt2242 = copt2237 + copt2241;
  Real copt2244 = copt16 + copt282;
  Real copt2245 = copt174 * copt2244;
  Real copt2249 = copt19 * copt274;
  Real copt2269 = copt181 * copt2268;
  Real copt2270 = copt280 * copt30;
  Real copt2271 = copt2245 + copt2249 + copt2269 + copt2270;
  Real copt2272 = copt121 * copt1276 * copt2271 * copt323;
  Real copt2273 = copt18 + copt282;
  Real copt2274 = copt2273 * copt27;
  Real copt2275 = copt272 * copt29;
  Real copt2276 = -(copt18 * copt278);
  Real copt2277 = copt103 + copt278;
  Real copt2278 = copt16 * copt2277;
  Real copt2279 = copt2274 + copt2275 + copt2276 + copt2278;
  Real copt2283 = -(copt121 * copt2279);
  Real copt2284 = copt111 * copt122 * copt323;
  Real copt2302 = copt2283 + copt2284;
  Real copt2303 = copt1276 * copt2302 * copt290;
  Real copt2304 = copt2272 + copt2303;
  Real copt2311 = 2 * copt174 * copt98;
  Real copt2315 = 2 * copt187 * copt30;
  Real copt2316 = copt2311 + copt2315;
  Real copt2345 = copt2336 + copt2344;
  Real copt2350 = copt174 * copt194;
  Real copt2365 = copt209 * copt98;
  Real copt2367 = copt187 * copt2165;
  Real copt2368 = copt224 * copt30;
  Real copt2369 = copt2350 + copt2365 + copt2367 + copt2368;
  Real copt2370 = -(copt132 * copt1514 * copt2369 * copt248);
  Real copt2375 = copt194 * copt27;
  Real copt2376 = copt192 * copt29;
  Real copt2377 = copt214 * copt98;
  Real copt2378 = copt2374 + copt2375 + copt2376 + copt2377;
  Real copt2379 = copt132 * copt2378;
  Real copt2380 = -(copt125 * copt133 * copt248);
  Real copt2384 = copt2379 + copt2380;
  Real copt2392 = copt1514 * copt226 * copt2384;
  Real copt2393 = copt2370 + copt2392;
  Real copt2396 = copt27 * copt340;
  Real copt2397 = copt29 * copt338;
  Real copt2401 = copt367 * copt98;
  Real copt2410 = copt2374 + copt2396 + copt2397 + copt2401;
  Real copt2411 = copt108 * copt1243 * copt2410 * copt405;
  Real copt2412 = copt347 * copt6;
  Real copt2413 = copt104 * copt401;
  Real copt2414 = copt2412 + copt2413;
  Real copt2415 = -(copt108 * copt1243 * copt2414 * copt697);
  Real copt2419 = copt2411 + copt2415;
  Real copt2421 = copt174 * copt267;
  Real copt2422 = copt274 * copt98;
  Real copt2423 = copt187 * copt2268;
  Real copt2442 = copt288 * copt30;
  Real copt2443 = copt2421 + copt2422 + copt2423 + copt2442;
  Real copt2444 = copt121 * copt1276 * copt2443 * copt323;
  Real copt2446 = -(copt263 * copt29);
  Real copt2447 = copt27 * copt303;
  Real copt2448 = copt278 * copt6;
  Real copt2449 = copt2445 + copt2446 + copt2447 + copt2448;
  Real copt2450 = -(copt121 * copt2449);
  Real copt2451 = copt114 * copt122 * copt323;
  Real copt2452 = copt2450 + copt2451;
  Real copt2453 = copt1276 * copt2452 * copt290;
  Real copt2454 = copt2444 + copt2453;
  Real copt2480 = 2 * copt181 * copt98;
  Real copt2481 = 2 * copt101 * copt187;
  Real copt2483 = copt2480 + copt2481;
  Real copt2511 = copt2486 + copt2510;
  Real copt2513 = copt181 * copt194;
  Real copt2514 = copt187 * copt219;
  Real copt2515 = copt216 * copt98;
  Real copt2516 = copt101 * copt224;
  Real copt2517 = copt2513 + copt2514 + copt2515 + copt2516;
  Real copt2518 = -(copt132 * copt1514 * copt248 * copt2517);
  Real copt2520 = -(copt18 * copt192);
  Real copt2521 = copt16 * copt242;
  Real copt2526 = -(copt200 * copt3);
  Real copt2538 = copt200 * copt5;
  Real copt2539 = copt2519 + copt2520 + copt2521 + copt2526 + copt2538;
  Real copt2540 = copt132 * copt2539;
  Real copt2544 = -(copt128 * copt133 * copt248);
  Real copt2547 = copt2540 + copt2544;
  Real copt2548 = copt1514 * copt226 * copt2547;
  Real copt2549 = copt2518 + copt2548;
  Real copt2562 = -(copt18 * copt338);
  Real copt2563 = copt16 * copt353;
  Real copt2567 = -(copt3 * copt345);
  Real copt2571 = copt345 * copt5;
  Real copt2572 = copt2519 + copt2562 + copt2563 + copt2567 + copt2571;
  Real copt2573 = copt108 * copt1243 * copt2572 * copt405;
  Real copt2576 = copt371 * copt98;
  Real copt2586 = copt19 * copt401;
  Real copt2589 = copt2576 + copt2586;
  Real copt2593 = -(copt108 * copt1243 * copt2589 * copt697);
  Real copt2594 = copt2573 + copt2593;
  Real copt2598 = copt181 * copt267;
  Real copt2604 = copt187 * copt285;
  Real copt2605 = copt280 * copt98;
  Real copt2606 = copt101 * copt288;
  Real copt2609 = copt2598 + copt2604 + copt2605 + copt2606;
  Real copt2617 = copt121 * copt1276 * copt2609 * copt323;
  Real copt2618 = -(copt18 * copt3);
  Real copt2619 = copt16 * copt316;
  Real copt2622 = copt18 * copt263;
  Real copt2624 = copt272 * copt3;
  Real copt2625 = -(copt272 * copt5);
  Real copt2626 = copt2618 + copt2619 + copt2622 + copt2624 + copt2625;
  Real copt2629 = -(copt121 * copt2626);
  Real copt2640 = copt117 * copt122 * copt323;
  Real copt2641 = copt2629 + copt2640;
  Real copt2642 = copt1276 * copt2641 * copt290;
  Real copt2644 = copt2617 + copt2642;
  Real copt2665 = copt108 * copt1243 * copt187 * copt405;
  Real copt2666 = copt19 * copt329;
  Real copt2667 = copt104 * copt181;
  Real copt2668 = copt2666 + copt2667;
  Real copt2671 = -(copt108 * copt1243 * copt2668 * copt697);
  Real copt2682 = copt2665 + copt2671;
  Real copt2687 = copt108 * copt1243 * copt2684 * copt405;
  Real copt2689 = copt329 * copt98;
  Real copt2690 = copt30 * copt380;
  Real copt2691 = copt2689 + copt2690;
  Real copt2692 = -(copt108 * copt1243 * copt2691 * copt697);
  Real copt2693 = copt2687 + copt2692;
  Real copt2695 = copt108 * copt1243 * copt174 * copt405;
  Real copt2696 = copt101 * copt380;
  Real copt2697 = copt181 * copt6;
  Real copt2698 = copt2696 + copt2697;
  Real copt2699 = -(copt108 * copt1243 * copt2698 * copt697);
  Real copt2701 = copt2695 + copt2699;
  Real copt2703 = copt125 * copt174;
  Real copt2704 = copt128 * copt181;
  Real copt2705 = copt2703 + copt2704;
  Real copt2706 = -(copt132 * copt1514 * copt248 * copt2705);
  Real copt2707 = copt132 * copt1514 * copt187 * copt226;
  Real copt2715 = copt2706 + copt2707;
  Real copt2744 = copt174 * copt199;
  Real copt2747 = copt128 * copt187;
  Real copt2750 = copt2744 + copt2747;
  Real copt2753 = -(copt132 * copt1514 * copt248 * copt2750);
  Real copt2755 = copt132 * copt1514 * copt226 * copt2684;
  Real copt2759 = copt2753 + copt2755;
  Real copt2762 = copt181 * copt199;
  Real copt2763 = copt172 * copt187;
  Real copt2765 = copt2762 + copt2763;
  Real copt2767 = -(copt132 * copt1514 * copt248 * copt2765);
  Real copt2768 = copt132 * copt1514 * copt174 * copt226;
  Real copt2770 = copt2767 + copt2768;
  Real copt2772 = copt174 * copt23;
  Real copt2773 = copt181 * copt36;
  Real copt2774 = copt2772 + copt2773;
  Real copt2776 = copt121 * copt1276 * copt2774 * copt323;
  Real copt2777 = -(copt121 * copt1276 * copt290 * copt380);
  Real copt2780 = copt2776 + copt2777;
  Real copt2782 = copt111 * copt174;
  Real copt2783 = copt187 * copt36;
  Real copt2784 = copt2782 + copt2783;
  Real copt2785 = copt121 * copt1276 * copt2784 * copt323;
  Real copt2787 = copt1586 + copt176 + copt177 + copt178 + copt2374;
  Real copt2788 = -(copt121 * copt1276 * copt2787 * copt290);
  Real copt2789 = copt2785 + copt2788;
  Real copt2791 = copt111 * copt181;
  Real copt2792 = copt114 * copt187;
  Real copt2793 = copt2791 + copt2792;
  Real copt2795 = copt121 * copt1276 * copt2793 * copt323;
  Real copt2796 = -(copt121 * copt1276 * copt290 * copt329);
  Real copt2797 = copt2795 + copt2796;
  Real copt3113 = 1 / copt721;
  Real copt779  = -(copt122 * copt325 * copt326 * copt778);
  Real copt783  = -(copt109 * copt699 * copt708 * copt782);
  Real copt790  = -(copt133 * copt250 * copt254 * copt787);
  Real copt793  = copt779 + copt783 + copt790;
  Real copt3122 = Power(copt773, 2);
  Real copt3123 = 1 / copt3122;
  Real copt2833 = copt2 * copt2832;
  Real copt2837 = copt2832 * copt9;
  Real copt3142 = 1 / copt119;
  Real copt2869 = copt2 * copt2868;
  Real copt2873 = copt2868 * copt9;
  Real copt2913 = copt2 * copt2912;
  Real copt2917 = copt2912 * copt9;
  Real copt2945 = copt110 + copt2944 + copt4;
  Real copt1565 = 2 * copt2 * copt23 * copt9;
  Real copt2974 = copt113 + copt1561 + copt17;
  Real copt1738 = 2 * copt2 * copt36 * copt9;
  Real copt3001 = copt116 + copt1694 + copt28;
  Real copt2112 = 2 * copt2 * copt6 * copt9;
  Real copt3029 = copt2122 + copt4 + copt95;
  Real copt2335 = 2 * copt19 * copt2 * copt9;
  Real copt3057 = copt100 + copt17 + copt2336;
  Real copt2485 = 2 * copt2 * copt30 * copt9;
  Real copt3085 = copt103 + copt2486 + copt28;
  out1(0)       = copt44;
  out1(1)       = copt48 * copt83 * copt91;
  out1(2)       = copt90;
  out1(3)       = -(copt109 * copt122 * copt133 * copt168 * copt191 * copt717);
  out1(4) = copt109 * copt122 * copt133 * copt44 * copt721 * copt725 * copt760;
  out1(5) = copt774 * copt775 * copt793;
  out2(0, 0)  = copt252 * copt48 * copt800;
  out2(0, 1)  = copt252 * copt48 * copt821;
  out2(0, 2)  = copt252 * copt48 * copt849;
  out2(0, 3)  = copt14 * copt2 * copt48;
  out2(0, 4)  = copt2 * copt25 * copt48;
  out2(0, 5)  = copt2 * copt38 * copt48;
  out2(0, 6)  = copt14 * copt48 * copt9;
  out2(0, 7)  = copt25 * copt48 * copt9;
  out2(0, 8)  = copt38 * copt48 * copt9;
  out2(0, 9)  = 0;
  out2(0, 10) = 0;
  out2(0, 11) = 0;
  out2(0, 12) = 0;
  out2(0, 13) = 0;
  out2(0, 14) = 0;
  out2(0, 15) = 0;
  out2(0, 16) = 0;
  out2(0, 17) = 0;
  out2(1, 0) =
      copt861 * copt863 *
      (copt43 * copt59 * copt83 * copt864 + copt14 * copt252 * copt83 * copt89 +
       copt43 * copt89 * (copt800 * copt864 + copt252 * copt912));
  out2(1, 1) =
      copt861 * copt863 *
      (copt43 * copt71 * copt83 * copt864 + copt25 * copt252 * copt83 * copt89 +
       copt43 * copt89 * (copt821 * copt864 + copt252 * copt956));
  out2(1, 2) =
      copt861 * copt863 *
      (copt43 * copt80 * copt83 * copt864 + copt252 * copt38 * copt83 * copt89 +
       copt43 * copt89 * (copt849 * copt864 + copt252 * copt976));
  out2(1, 3) =
      copt861 * copt863 *
      (-(copt43 * copt51 * copt59 * copt83) - copt14 * copt2 * copt83 * copt89 +
       copt43 * copt89 *
           (copt12 * copt51 * copt9 + copt2 * (copt58 - 2 * copt51 * copt98)));
  out2(1, 4) =
      copt861 * copt863 *
      (-(copt43 * copt51 * copt71 * copt83) - copt2 * copt25 * copt83 * copt89 +
       copt43 * copt89 *
           (copt2 * (-2 * copt101 * copt51 + copt69) +
            copt23 * copt51 * copt9));
  out2(1, 5) =
      copt861 * copt863 *
      (-(copt43 * copt51 * copt80 * copt83) - copt2 * copt38 * copt83 * copt89 +
       copt43 * copt89 *
           (copt2 * (-2 * copt104 * copt51 + copt79) +
            copt36 * copt51 * copt9));
  out2(1, 6) =
      copt861 * copt863 *
      (-(copt43 * copt56 * copt59 * copt83) - copt14 * copt83 * copt89 * copt9 +
       copt43 * copt89 *
           (copt51 * copt6 * copt9 + copt56 * (copt8 + 2 * copt12 * copt9)));
  out2(1, 7) =
      copt861 * copt863 *
      (-(copt43 * copt56 * copt71 * copt83) - copt25 * copt83 * copt89 * copt9 +
       copt43 * copt89 *
           (copt19 * copt51 * copt9 + copt56 * (copt20 + 2 * copt23 * copt9)));
  out2(1, 8) =
      copt861 * copt863 *
      (-(copt43 * copt56 * copt80 * copt83) - copt38 * copt83 * copt89 * copt9 +
       copt43 * copt89 *
           (copt30 * copt51 * copt9 + copt56 * (copt32 + 2 * copt36 * copt9)));
  out2(1, 9)  = 0;
  out2(1, 10) = 0;
  out2(1, 11) = 0;
  out2(1, 12) = 0;
  out2(1, 13) = 0;
  out2(1, 14) = 0;
  out2(1, 15) = 0;
  out2(1, 16) = 0;
  out2(1, 17) = 0;
  out2(2, 0)  = copt864 * copt91 * copt912;
  out2(2, 1)  = copt864 * copt91 * copt956;
  out2(2, 2)  = copt864 * copt91 * copt976;
  out2(2, 3)  = copt51 * copt59 * copt91;
  out2(2, 4)  = copt51 * copt71 * copt91;
  out2(2, 5)  = copt51 * copt80 * copt91;
  out2(2, 6)  = copt56 * copt59 * copt91;
  out2(2, 7)  = copt56 * copt71 * copt91;
  out2(2, 8)  = copt56 * copt80 * copt91;
  out2(2, 9)  = 0;
  out2(2, 10) = 0;
  out2(2, 11) = 0;
  out2(2, 12) = 0;
  out2(2, 13) = 0;
  out2(2, 14) = 0;
  out2(2, 15) = 0;
  out2(2, 16) = 0;
  out2(2, 17) = 0;
  out2(3, 0) =
      copt109 * copt111 * copt1144 * copt133 * copt168 * copt191 * copt717 -
      (copt109 * copt1132 * copt122 * copt133 * copt168 * copt717 * copt775) /
          2. +
      2 * copt109 * copt1135 * copt122 * copt133 * copt191 * copt252 * copt717 *
          copt800 +
      copt1148 * copt122 * copt133 * copt168 * copt191 * copt717 * copt98 -
      copt109 * copt122 * copt133 * copt168 * copt191 *
          (copt111 * copt1152 * copt250 * copt253 * copt254 * copt255 -
           copt1232 * copt132 * copt253 * copt254 * copt255 * copt256 *
               copt380 +
           copt1154 * copt250 * copt253 * copt254 * copt256 * copt98 +
           copt262 *
               (copt122 * copt1276 * copt1311 * copt134 * copt255 * copt326 +
                copt1243 * copt1269 * copt158 * copt256 * copt708 +
                copt111 * copt1152 * copt158 * copt699 * copt708 +
                copt1154 * copt134 * copt325 * copt326 * copt98));
  out2(3, 1) =
      -(copt109 * copt122 * copt133 * copt168 * copt191 *
        (copt114 * copt1152 * copt250 * copt253 * copt254 * copt255 +
         copt101 * copt1154 * copt250 * copt253 * copt254 * copt256 +
         copt1232 * copt132 * copt1336 * copt253 * copt254 * copt255 * copt256 +
         copt262 *
             (copt122 * copt1276 * copt134 * copt1380 * copt255 * copt326 +
              copt101 * copt1154 * copt134 * copt325 * copt326 +
              copt1243 * copt1354 * copt158 * copt256 * copt708 +
              copt114 * copt1152 * copt158 * copt699 * copt708))) +
      copt109 * copt114 * copt1144 * copt133 * copt168 * copt191 * copt717 +
      copt101 * copt1148 * copt122 * copt133 * copt168 * copt191 * copt717 -
      (copt109 * copt122 * copt1321 * copt133 * copt168 * copt717 * copt775) /
          2. +
      2 * copt109 * copt1135 * copt122 * copt133 * copt191 * copt252 * copt717 *
          copt821;
  out2(3, 2) =
      -(copt109 * copt122 * copt133 * copt168 * copt191 *
        (copt1152 * copt117 * copt250 * copt253 * copt254 * copt255 +
         copt104 * copt1154 * copt250 * copt253 * copt254 * copt256 +
         copt1232 * copt132 * copt174 * copt253 * copt254 * copt255 * copt256 +
         copt262 *
             (copt122 * copt1276 * copt134 * copt1433 * copt255 * copt326 +
              copt104 * copt1154 * copt134 * copt325 * copt326 +
              copt1243 * copt1414 * copt158 * copt256 * copt708 +
              copt1152 * copt117 * copt158 * copt699 * copt708))) +
      copt109 * copt1144 * copt117 * copt133 * copt168 * copt191 * copt717 +
      copt104 * copt1148 * copt122 * copt133 * copt168 * copt191 * copt717 -
      (copt109 * copt122 * copt133 * copt1389 * copt168 * copt717 * copt775) /
          2. +
      2 * copt109 * copt1135 * copt122 * copt133 * copt191 * copt252 * copt717 *
          copt849;
  out2(3, 3) =
      copt109 * copt122 * copt123 * copt1468 * copt168 * copt191 * copt717 +
      2 * copt109 * copt1135 * copt122 * copt133 * copt14 * copt191 * copt2 *
          copt717 +
      copt1148 * copt122 * copt133 * copt168 * copt191 * copt6 * copt717 -
      (copt109 * copt122 * copt133 * copt1457 * copt168 * copt717 * copt775) /
          2. -
      copt109 * copt122 * copt133 * copt168 * copt191 *
          (copt1514 * copt1550 * copt253 * copt254 * copt255 * copt256 +
           copt123 * copt1472 * copt714 -
           copt1154 * copt250 * copt253 * copt254 * copt256 * copt98 +
           copt262 *
               (copt121 * copt1232 * copt134 * copt187 * copt255 * copt326 +
                copt1243 * copt1501 * copt158 * copt256 * copt708 -
                copt1154 * copt134 * copt325 * copt326 * copt98));
  out2(3, 4) =
      -(copt109 * copt122 * copt133 * copt168 * copt191 *
        ((copt1154 * copt1563 * copt250 * copt253 * copt254 * copt256) / 2. +
         copt1597 * copt253 * copt254 * copt255 * copt256 +
         copt262 * (copt134 * copt1664 * copt255 * copt326 +
                    (copt1154 * copt134 * copt1563 * copt325 * copt326) / 2. +
                    copt158 * copt1634 * copt256 * copt708) +
         (copt1472 * copt1577 * copt714) / 2.)) +
      copt109 * copt1135 * copt122 * copt133 * (copt134 * copt1563 + copt1565) *
          copt191 * copt717 -
      copt101 * copt1148 * copt122 * copt133 * copt168 * copt191 * copt717 +
      copt109 * copt122 * copt125 * copt1468 * copt168 * copt191 * copt717 -
      (copt109 * copt122 * copt133 * copt1558 * copt168 * copt717 * copt775) /
          2.;
  out2(3, 5) =
      -(copt109 * copt122 * copt133 * copt168 * copt191 *
        ((copt1154 * copt1722 * copt250 * copt253 * copt254 * copt256) / 2. +
         copt1856 * copt253 * copt254 * copt255 * copt256 +
         copt262 * (copt134 * copt2059 * copt255 * copt326 +
                    (copt1154 * copt134 * copt1722 * copt325 * copt326) / 2. +
                    copt158 * copt1960 * copt256 * copt708) +
         (copt1472 * copt1750 * copt714) / 2.)) -
      copt104 * copt1148 * copt122 * copt133 * copt168 * copt191 * copt717 +
      copt109 * copt122 * copt128 * copt1468 * copt168 * copt191 * copt717 +
      copt109 * copt1135 * copt122 * copt133 * (copt134 * copt1722 + copt1738) *
          copt191 * copt717 -
      (copt109 * copt122 * copt133 * copt168 * copt1689 * copt717 * copt775) /
          2.;
  out2(3, 6) =
      -(copt109 * copt122 * copt133 * copt168 * copt191 *
        ((copt1152 * copt2123 * copt250 * copt253 * copt254 * copt255) / 2. +
         copt2203 * copt253 * copt254 * copt255 * copt256 +
         copt262 * (copt134 * copt2304 * copt255 * copt326 +
                    copt158 * copt2242 * copt256 * copt708 +
                    (copt1152 * copt158 * copt2123 * copt699 * copt708) / 2.) +
         (copt1472 * copt2160 * copt714) / 2.)) -
      copt109 * copt111 * copt1144 * copt133 * copt168 * copt191 * copt717 -
      copt109 * copt122 * copt123 * copt1468 * copt168 * copt191 * copt717 +
      copt109 * copt1135 * copt122 * copt133 * copt191 *
          (copt2112 + copt158 * copt2123) * copt717 -
      (copt109 * copt122 * copt133 * copt168 * copt2105 * copt717 * copt775) /
          2.;
  out2(3, 7) =
      -(copt109 * copt122 * copt133 * copt168 * copt191 *
        ((copt1152 * copt2337 * copt250 * copt253 * copt254 * copt255) / 2. +
         copt2393 * copt253 * copt254 * copt255 * copt256 +
         copt262 * (copt134 * copt2454 * copt255 * copt326 +
                    copt158 * copt2419 * copt256 * copt708 +
                    (copt1152 * copt158 * copt2337 * copt699 * copt708) / 2.) +
         (copt1472 * copt2345 * copt714) / 2.)) -
      copt109 * copt114 * copt1144 * copt133 * copt168 * copt191 * copt717 -
      copt109 * copt122 * copt125 * copt1468 * copt168 * copt191 * copt717 +
      copt109 * copt1135 * copt122 * copt133 * copt191 *
          (copt2335 + copt158 * copt2337) * copt717 -
      (copt109 * copt122 * copt133 * copt168 * copt2316 * copt717 * copt775) /
          2.;
  out2(3, 8) =
      -(copt109 * copt122 * copt133 * copt168 * copt191 *
        ((copt1152 * copt2487 * copt250 * copt253 * copt254 * copt255) / 2. +
         copt253 * copt254 * copt2549 * copt255 * copt256 +
         copt262 * (copt134 * copt255 * copt2644 * copt326 +
                    copt158 * copt256 * copt2594 * copt708 +
                    (copt1152 * copt158 * copt2487 * copt699 * copt708) / 2.) +
         (copt1472 * copt2511 * copt714) / 2.)) -
      copt109 * copt1144 * copt117 * copt133 * copt168 * copt191 * copt717 -
      copt109 * copt122 * copt128 * copt1468 * copt168 * copt191 * copt717 +
      copt109 * copt1135 * copt122 * copt133 * copt191 *
          (copt2485 + copt158 * copt2487) * copt717 -
      (copt109 * copt122 * copt133 * copt168 * copt2483 * copt717 * copt775) /
          2.;
  out2(3, 9)  = -(copt109 * copt122 * copt133 * copt158 * copt168 * copt191 *
                 copt256 * copt262 * copt2682 * copt708);
  out2(3, 10) = -(copt109 * copt122 * copt133 * copt158 * copt168 * copt191 *
                  copt256 * copt262 * copt2693 * copt708);
  out2(3, 11) = -(copt109 * copt122 * copt133 * copt158 * copt168 * copt191 *
                  copt256 * copt262 * copt2701 * copt708);
  out2(3, 12) = -(copt109 * copt122 * copt133 * copt168 * copt191 * copt253 *
                  copt254 * copt255 * copt256 * copt2715);
  out2(3, 13) = -(copt109 * copt122 * copt133 * copt168 * copt191 * copt253 *
                  copt254 * copt255 * copt256 * copt2759);
  out2(3, 14) = -(copt109 * copt122 * copt133 * copt168 * copt191 * copt253 *
                  copt254 * copt255 * copt256 * copt2770);
  out2(3, 15) = -(copt109 * copt122 * copt133 * copt134 * copt168 * copt191 *
                  copt255 * copt262 * copt2780 * copt326);
  out2(3, 16) = -(copt109 * copt122 * copt133 * copt134 * copt168 * copt191 *
                  copt255 * copt262 * copt2789 * copt326);
  out2(3, 17) = -(copt109 * copt122 * copt133 * copt134 * copt168 * copt191 *
                  copt255 * copt262 * copt2797 * copt326);
  out2(4, 0) =
      (-3 * copt109 * copt122 * copt133 * copt2809 *
       (2 * copt14 * copt190 * copt2803 + copt1132 * copt43) * copt44 *
       copt721 * copt760) /
          2. +
      (3 * copt109 * copt1132 * copt122 * copt133 * copt191 * copt44 * copt725 *
       copt760) /
          2. -
      copt109 * copt111 * copt1144 * copt133 * copt44 * copt721 * copt725 *
          copt760 +
      copt109 * copt122 * copt133 * copt14 * copt2803 * copt48 * copt721 *
          copt725 * copt760 +
      copt109 * copt122 * copt133 * copt44 * copt721 * copt725 *
          ((copt1152 * copt250 * copt252 * copt254 * copt255 * copt2826 *
            copt746) /
               2. +
           (copt1154 * copt250 * copt252 * copt254 * copt256 * copt2828 *
            copt746) /
               2. -
           copt1232 * copt132 * copt252 * copt254 * copt255 * copt256 *
               copt380 * copt746 +
           copt250 * copt252 * copt254 * copt255 * copt256 *
               (-(copt199 * copt2) + copt123 * copt9) +
           copt262 *
               (copt122 * copt1276 * copt1311 * copt2 * copt255 * copt326 *
                    copt756 +
                (copt1154 * copt2 * copt2828 * copt325 * copt326 * copt756) /
                    2. -
                copt256 * (copt2 * copt2828 + copt2837) * copt699 * copt708 *
                    copt9 -
                copt1243 * copt1269 * copt256 * copt708 * copt750 * copt9 -
                (copt1152 * copt2826 * copt699 * copt708 * copt750 * copt9) /
                    2. +
                copt2 * copt255 * copt325 * copt326 *
                    (copt2833 + copt2826 * copt9))) -
      copt1148 * copt122 * copt133 * copt44 * copt721 * copt725 * copt760 *
          copt98;
  out2(4, 1) =
      (-3 * copt109 * copt122 * copt133 * copt2809 *
       (2 * copt190 * copt25 * copt2803 + copt1321 * copt43) * copt44 *
       copt721 * copt760) /
          2. +
      (3 * copt109 * copt122 * copt1321 * copt133 * copt191 * copt44 * copt725 *
       copt760) /
          2. -
      copt109 * copt114 * copt1144 * copt133 * copt44 * copt721 * copt725 *
          copt760 -
      copt101 * copt1148 * copt122 * copt133 * copt44 * copt721 * copt725 *
          copt760 +
      copt109 * copt122 * copt133 * copt25 * copt2803 * copt48 * copt721 *
          copt725 * copt760 +
      copt109 * copt122 * copt133 * copt44 * copt721 * copt725 *
          (copt1232 * copt132 * copt1336 * copt252 * copt254 * copt255 *
               copt256 * copt746 +
           (copt1152 * copt250 * copt252 * copt254 * copt255 * copt2862 *
            copt746) /
               2. +
           (copt1154 * copt250 * copt252 * copt254 * copt256 * copt2864 *
            copt746) /
               2. +
           copt250 * copt252 * copt254 * copt255 * copt256 *
               (-(copt172 * copt2) + copt125 * copt9) +
           copt262 *
               (copt122 * copt1276 * copt1380 * copt2 * copt255 * copt326 *
                    copt756 +
                (copt1154 * copt2 * copt2864 * copt325 * copt326 * copt756) /
                    2. -
                copt256 * (copt2 * copt2864 + copt2873) * copt699 * copt708 *
                    copt9 -
                copt1243 * copt1354 * copt256 * copt708 * copt750 * copt9 -
                (copt1152 * copt2862 * copt699 * copt708 * copt750 * copt9) /
                    2. +
                copt2 * copt255 * copt325 * copt326 *
                    (copt2869 + copt2862 * copt9)));
  out2(4, 2) =
      (-3 * copt109 * copt122 * copt133 * copt2809 *
       (2 * copt190 * copt2803 * copt38 + copt1389 * copt43) * copt44 *
       copt721 * copt760) /
          2. +
      (3 * copt109 * copt122 * copt133 * copt1389 * copt191 * copt44 * copt725 *
       copt760) /
          2. -
      copt109 * copt1144 * copt117 * copt133 * copt44 * copt721 * copt725 *
          copt760 -
      copt104 * copt1148 * copt122 * copt133 * copt44 * copt721 * copt725 *
          copt760 +
      copt109 * copt122 * copt133 * copt2803 * copt38 * copt48 * copt721 *
          copt725 * copt760 +
      copt109 * copt122 * copt133 * copt44 * copt721 * copt725 *
          (copt1232 * copt132 * copt174 * copt252 * copt254 * copt255 *
               copt256 * copt746 +
           (copt1152 * copt250 * copt252 * copt254 * copt255 * copt2906 *
            copt746) /
               2. +
           (copt1154 * copt250 * copt252 * copt254 * copt256 * copt2908 *
            copt746) /
               2. +
           copt250 * copt252 * copt254 * copt255 * copt256 *
               (-(copt179 * copt2) + copt128 * copt9) +
           copt262 *
               (copt122 * copt1276 * copt1433 * copt2 * copt255 * copt326 *
                    copt756 +
                (copt1154 * copt2 * copt2908 * copt325 * copt326 * copt756) /
                    2. -
                copt256 * (copt2 * copt2908 + copt2917) * copt699 * copt708 *
                    copt9 -
                copt1243 * copt1414 * copt256 * copt708 * copt750 * copt9 -
                (copt1152 * copt2906 * copt699 * copt708 * copt750 * copt9) /
                    2. +
                copt2 * copt255 * copt325 * copt326 *
                    (copt2913 + copt2906 * copt9)));
  out2(4, 3) =
      (-3 * copt109 * copt122 * copt133 * copt2809 *
       (2 * copt14 * copt190 * copt2 + copt1457 * copt43) * copt44 * copt721 *
       copt760) /
          2. +
      (3 * copt109 * copt122 * copt133 * copt1457 * copt191 * copt44 * copt725 *
       copt760) /
          2. -
      copt109 * copt122 * copt123 * copt1468 * copt44 * copt721 * copt725 *
          copt760 +
      copt109 * copt122 * copt133 * copt14 * copt2 * copt48 * copt721 *
          copt725 * copt760 +
      copt109 * copt122 * copt133 * copt44 * copt721 * copt725 *
          (copt1514 * copt1550 * copt252 * copt254 * copt255 * copt256 *
               copt746 +
           (copt1154 * copt250 * copt252 * copt254 * copt256 * copt2949 *
            copt746) /
               2. +
           (copt1472 * (copt2822 + copt2944) * copt758) / 2. +
           copt250 * copt252 * copt254 * copt255 * copt256 *
               (-(copt2 * copt2945) + copt799) +
           copt262 *
               (copt12 * copt134 * copt255 * copt325 * copt326 +
                copt121 * copt1232 * copt187 * copt2 * copt255 * copt326 *
                    copt756 +
                (copt1154 * copt2 * copt2949 * copt325 * copt326 * copt756) /
                    2. -
                copt256 * (copt13 + copt2 * copt2949) * copt699 * copt708 *
                    copt9 -
                copt1243 * copt1501 * copt256 * copt708 * copt750 * copt9)) +
      copt1148 * copt122 * copt133 * copt44 * copt721 * copt725 * copt760 *
          copt98;
  out2(4, 4) =
      (-3 * copt109 * copt122 * copt133 * copt2809 *
       (2 * copt190 * copt2 * copt25 + copt1558 * copt43) * copt44 * copt721 *
       copt760) /
          2. +
      (3 * copt109 * copt122 * copt133 * copt1558 * copt191 * copt44 * copt725 *
       copt760) /
          2. +
      copt101 * copt1148 * copt122 * copt133 * copt44 * copt721 * copt725 *
          copt760 -
      copt109 * copt122 * copt125 * copt1468 * copt44 * copt721 * copt725 *
          copt760 +
      copt109 * copt122 * copt133 * copt2 * copt25 * copt48 * copt721 *
          copt725 * copt760 +
      copt109 * copt122 * copt133 * copt44 * copt721 * copt725 *
          ((copt1154 * copt1563 * copt250 * copt252 * copt254 * copt256 *
            copt746) /
               2. +
           copt1597 * copt252 * copt254 * copt255 * copt256 * copt746 +
           (copt1472 * copt1577 * copt758) / 2. +
           copt250 * copt252 * copt254 * copt255 * copt256 *
               (-(copt2 * copt2974) + copt820) +
           copt262 *
               (copt134 * copt23 * copt255 * copt325 * copt326 +
                copt1664 * copt2 * copt255 * copt326 * copt756 +
                (copt1154 * copt1563 * copt2 * copt325 * copt326 * copt756) /
                    2. -
                (copt1563 * copt2 + copt24) * copt256 * copt699 * copt708 *
                    copt9 -
                copt1634 * copt256 * copt708 * copt750 * copt9));
  out2(4, 5) =
      (-3 * copt109 * copt122 * copt133 * copt2809 *
       (2 * copt190 * copt2 * copt38 + copt1689 * copt43) * copt44 * copt721 *
       copt760) /
          2. +
      (3 * copt109 * copt122 * copt133 * copt1689 * copt191 * copt44 * copt725 *
       copt760) /
          2. +
      copt104 * copt1148 * copt122 * copt133 * copt44 * copt721 * copt725 *
          copt760 -
      copt109 * copt122 * copt128 * copt1468 * copt44 * copt721 * copt725 *
          copt760 +
      copt109 * copt122 * copt133 * copt2 * copt38 * copt48 * copt721 *
          copt725 * copt760 +
      copt109 * copt122 * copt133 * copt44 * copt721 * copt725 *
          ((copt1154 * copt1722 * copt250 * copt252 * copt254 * copt256 *
            copt746) /
               2. +
           copt1856 * copt252 * copt254 * copt255 * copt256 * copt746 +
           (copt1472 * copt1750 * copt758) / 2. +
           copt250 * copt252 * copt254 * copt255 * copt256 *
               (-(copt2 * copt3001) + copt845) +
           copt262 *
               (copt134 * copt255 * copt325 * copt326 * copt36 +
                copt2 * copt2059 * copt255 * copt326 * copt756 +
                (copt1154 * copt1722 * copt2 * copt325 * copt326 * copt756) /
                    2. -
                copt256 * (copt1722 * copt2 + copt37) * copt699 * copt708 *
                    copt9 -
                copt1960 * copt256 * copt708 * copt750 * copt9));
  out2(4, 6) =
      (3 * copt109 * copt122 * copt133 * copt191 * copt2105 * copt44 * copt725 *
       copt760) /
          2. +
      copt109 * copt111 * copt1144 * copt133 * copt44 * copt721 * copt725 *
          copt760 +
      copt109 * copt122 * copt123 * copt1468 * copt44 * copt721 * copt725 *
          copt760 +
      copt109 * copt122 * copt133 * copt14 * copt48 * copt721 * copt725 *
          copt760 * copt9 -
      (3 * copt109 * copt122 * copt133 * copt2809 * copt44 * copt721 * copt760 *
       (copt2105 * copt43 + 2 * copt14 * copt190 * copt9)) /
          2. +
      copt109 * copt122 * copt133 * copt44 * copt721 * copt725 *
          ((copt1152 * copt2123 * copt250 * copt252 * copt254 * copt255 *
            copt746) /
               2. +
           copt2203 * copt252 * copt254 * copt255 * copt256 * copt746 +
           (copt1472 * copt2160 * copt758) / 2. +
           copt262 *
               (-(copt158 * copt256 * copt6 * copt699 * copt708) +
                copt2 * copt2304 * copt255 * copt326 * copt756 -
                copt2242 * copt256 * copt708 * copt750 * copt9 -
                (copt1152 * copt2123 * copt699 * copt708 * copt750 * copt9) /
                    2. +
                copt2 * copt255 * copt325 * copt326 *
                    (copt8 + copt2123 * copt9)) +
           copt250 * copt252 * copt254 * copt255 * copt256 *
               (copt3029 * copt9 - copt2 * copt98));
  out2(4, 7) =
      (3 * copt109 * copt122 * copt133 * copt191 * copt2316 * copt44 * copt725 *
       copt760) /
          2. +
      copt109 * copt114 * copt1144 * copt133 * copt44 * copt721 * copt725 *
          copt760 +
      copt109 * copt122 * copt125 * copt1468 * copt44 * copt721 * copt725 *
          copt760 +
      copt109 * copt122 * copt133 * copt25 * copt48 * copt721 * copt725 *
          copt760 * copt9 -
      (3 * copt109 * copt122 * copt133 * copt2809 * copt44 * copt721 * copt760 *
       (copt2316 * copt43 + 2 * copt190 * copt25 * copt9)) /
          2. +
      copt109 * copt122 * copt133 * copt44 * copt721 * copt725 *
          ((copt1152 * copt2337 * copt250 * copt252 * copt254 * copt255 *
            copt746) /
               2. +
           copt2393 * copt252 * copt254 * copt255 * copt256 * copt746 +
           (copt1472 * copt2345 * copt758) / 2. +
           copt250 * copt252 * copt254 * copt255 * copt256 *
               (-(copt101 * copt2) + copt3057 * copt9) +
           copt262 *
               (-(copt158 * copt19 * copt256 * copt699 * copt708) +
                copt2 * copt2454 * copt255 * copt326 * copt756 -
                copt2419 * copt256 * copt708 * copt750 * copt9 -
                (copt1152 * copt2337 * copt699 * copt708 * copt750 * copt9) /
                    2. +
                copt2 * copt255 * copt325 * copt326 *
                    (copt20 + copt2337 * copt9)));
  out2(4, 8) =
      (3 * copt109 * copt122 * copt133 * copt191 * copt2483 * copt44 * copt725 *
       copt760) /
          2. +
      copt109 * copt1144 * copt117 * copt133 * copt44 * copt721 * copt725 *
          copt760 +
      copt109 * copt122 * copt128 * copt1468 * copt44 * copt721 * copt725 *
          copt760 +
      copt109 * copt122 * copt133 * copt38 * copt48 * copt721 * copt725 *
          copt760 * copt9 -
      (3 * copt109 * copt122 * copt133 * copt2809 * copt44 * copt721 * copt760 *
       (copt2483 * copt43 + 2 * copt190 * copt38 * copt9)) /
          2. +
      copt109 * copt122 * copt133 * copt44 * copt721 * copt725 *
          ((copt1152 * copt2487 * copt250 * copt252 * copt254 * copt255 *
            copt746) /
               2. +
           copt252 * copt254 * copt2549 * copt255 * copt256 * copt746 +
           (copt1472 * copt2511 * copt758) / 2. +
           copt250 * copt252 * copt254 * copt255 * copt256 *
               (-(copt104 * copt2) + copt3085 * copt9) +
           copt262 *
               (-(copt158 * copt256 * copt30 * copt699 * copt708) +
                copt2 * copt255 * copt2644 * copt326 * copt756 -
                copt256 * copt2594 * copt708 * copt750 * copt9 -
                (copt1152 * copt2487 * copt699 * copt708 * copt750 * copt9) /
                    2. +
                copt2 * copt255 * copt325 * copt326 *
                    (copt32 + copt2487 * copt9)));
  out2(4, 9)  = -(copt109 * copt122 * copt133 * copt256 * copt262 * copt2682 *
                 copt44 * copt708 * copt721 * copt725 * copt750 * copt9);
  out2(4, 10) = -(copt109 * copt122 * copt133 * copt256 * copt262 * copt2693 *
                  copt44 * copt708 * copt721 * copt725 * copt750 * copt9);
  out2(4, 11) = -(copt109 * copt122 * copt133 * copt256 * copt262 * copt2701 *
                  copt44 * copt708 * copt721 * copt725 * copt750 * copt9);
  out2(4, 12) = copt109 * copt122 * copt133 * copt252 * copt254 * copt255 *
                copt256 * copt2715 * copt44 * copt721 * copt725 * copt746;
  out2(4, 13) = copt109 * copt122 * copt133 * copt252 * copt254 * copt255 *
                copt256 * copt2759 * copt44 * copt721 * copt725 * copt746;
  out2(4, 14) = copt109 * copt122 * copt133 * copt252 * copt254 * copt255 *
                copt256 * copt2770 * copt44 * copt721 * copt725 * copt746;
  out2(4, 15) = copt109 * copt122 * copt133 * copt2 * copt255 * copt262 *
                copt2780 * copt326 * copt44 * copt721 * copt725 * copt756;
  out2(4, 16) = copt109 * copt122 * copt133 * copt2 * copt255 * copt262 *
                copt2789 * copt326 * copt44 * copt721 * copt725 * copt756;
  out2(4, 17) = copt109 * copt122 * copt133 * copt2 * copt255 * copt262 *
                copt2797 * copt326 * copt44 * copt721 * copt725 * copt756;
  out2(5, 0) = -(copt1132 * copt3113 * copt774 * copt793) / 2. -
               copt3123 * copt775 * copt793 *
                   (2 * copt111 * copt158 + 2 * copt2 * copt2832 * copt9 +
                    2 * copt134 * copt98) +
               copt774 * copt775 *
                   (-(copt1276 * copt1311 * copt3142 * copt326 * copt778) +
                    copt111 * copt1144 * copt325 * copt326 * copt778 -
                    copt109 * copt1243 * copt1269 * copt708 * copt782 +
                    copt1232 * copt254 * copt380 * copt787 -
                    2 * copt122 * copt325 * copt326 * copt777 *
                        (copt2833 + 2 * copt111 * copt9) -
                    2 * copt133 * copt250 * copt254 * copt786 *
                        (copt199 * copt2 - copt123 * copt9) +
                    copt1148 * copt699 * copt708 * copt782 * copt98 -
                    2 * copt109 * copt699 * copt708 * copt781 *
                        (copt2837 + 2 * copt2 * copt98));
  out2(5, 1) = -(copt1321 * copt3113 * copt774 * copt793) / 2. -
               copt3123 * copt775 * copt793 *
                   (2 * copt101 * copt134 + 2 * copt114 * copt158 +
                    2 * copt2 * copt2868 * copt9) +
               copt774 * copt775 *
                   (-(copt1276 * copt1380 * copt3142 * copt326 * copt778) +
                    copt114 * copt1144 * copt325 * copt326 * copt778 -
                    2 * copt109 * (2 * copt101 * copt2 + copt2873) * copt699 *
                        copt708 * copt781 -
                    copt109 * copt1243 * copt1354 * copt708 * copt782 +
                    copt101 * copt1148 * copt699 * copt708 * copt782 -
                    copt1232 * copt1336 * copt254 * copt787 -
                    2 * copt122 * copt325 * copt326 * copt777 *
                        (copt2869 + 2 * copt114 * copt9) -
                    2 * copt133 * copt250 * copt254 * copt786 *
                        (copt172 * copt2 - copt125 * copt9));
  out2(5, 2) = -(copt1389 * copt3113 * copt774 * copt793) / 2. -
               copt3123 * copt775 * copt793 *
                   (2 * copt104 * copt134 + 2 * copt117 * copt158 +
                    2 * copt2 * copt2912 * copt9) +
               copt774 * copt775 *
                   (-(copt1276 * copt1433 * copt3142 * copt326 * copt778) +
                    copt1144 * copt117 * copt325 * copt326 * copt778 -
                    2 * copt109 * (2 * copt104 * copt2 + copt2917) * copt699 *
                        copt708 * copt781 -
                    copt109 * copt1243 * copt1414 * copt708 * copt782 +
                    copt104 * copt1148 * copt699 * copt708 * copt782 -
                    copt1232 * copt174 * copt254 * copt787 -
                    2 * copt122 * copt325 * copt326 * copt777 *
                        (copt2913 + 2 * copt117 * copt9) -
                    2 * copt133 * copt250 * copt254 * copt786 *
                        (copt179 * copt2 - copt128 * copt9));
  out2(5, 3) =
      -(copt1457 * copt3113 * copt774 * copt793) / 2. -
      copt3123 * copt775 * copt793 *
          (2 * copt12 * copt2 * copt9 - 2 * copt134 * copt98) +
      copt774 * copt775 *
          (-2 * copt12 * copt122 * copt2 * copt325 * copt326 * copt777 -
           copt1232 * copt187 * copt326 * copt778 -
           copt109 * copt1243 * copt1501 * copt708 * copt782 -
           copt133 * copt1514 * copt1550 * copt254 * copt787 +
           copt123 * copt1468 * copt250 * copt254 * copt787 -
           2 * copt133 * copt250 * copt254 * copt786 *
               (copt2 * copt2945 - copt111 * copt9) -
           copt1148 * copt699 * copt708 * copt782 * copt98 -
           2 * copt109 * copt699 * copt708 * copt781 *
               (copt13 - 2 * copt2 * copt98));
  out2(5, 4) =
      -(copt1558 * copt3113 * copt774 * copt793) / 2. -
      (-2 * copt101 * copt134 + copt1565) * copt3123 * copt775 * copt793 +
      copt774 * copt775 *
          (-2 * copt122 * copt2 * copt23 * copt325 * copt326 * copt777 -
           copt122 * copt1664 * copt326 * copt778 -
           2 * copt109 * (-2 * copt101 * copt2 + copt24) * copt699 * copt708 *
               copt781 -
           copt109 * copt1634 * copt708 * copt782 -
           copt101 * copt1148 * copt699 * copt708 * copt782 -
           copt133 * copt1597 * copt254 * copt787 +
           copt125 * copt1468 * copt250 * copt254 * copt787 -
           2 * copt133 * copt250 * copt254 * copt786 *
               (copt2 * copt2974 - copt114 * copt9));
  out2(5, 5) =
      -(copt1689 * copt3113 * copt774 * copt793) / 2. -
      (-2 * copt104 * copt134 + copt1738) * copt3123 * copt775 * copt793 +
      copt774 * copt775 *
          (-2 * copt122 * copt2 * copt325 * copt326 * copt36 * copt777 -
           copt122 * copt2059 * copt326 * copt778 -
           2 * copt109 * (-2 * copt104 * copt2 + copt37) * copt699 * copt708 *
               copt781 -
           copt109 * copt1960 * copt708 * copt782 -
           copt104 * copt1148 * copt699 * copt708 * copt782 -
           copt133 * copt1856 * copt254 * copt787 +
           copt128 * copt1468 * copt250 * copt254 * copt787 -
           2 * copt133 * copt250 * copt254 * copt786 *
               (copt2 * copt3001 - copt117 * copt9));
  out2(5, 6) =
      -(copt2105 * copt3113 * copt774 * copt793) / 2. -
      (-2 * copt111 * copt158 + copt2112) * copt3123 * copt775 * copt793 +
      copt774 * copt775 *
          (-(copt122 * copt2304 * copt326 * copt778) -
           copt111 * copt1144 * copt325 * copt326 * copt778 -
           copt109 * copt2242 * copt708 * copt782 -
           copt133 * copt2203 * copt254 * copt787 -
           copt123 * copt1468 * copt250 * copt254 * copt787 -
           2 * copt109 * copt6 * copt699 * copt708 * copt781 * copt9 -
           2 * copt122 * copt325 * copt326 * copt777 *
               (copt8 - 2 * copt111 * copt9) -
           2 * copt133 * copt250 * copt254 * copt786 *
               (copt796 - copt3029 * copt9));
  out2(5, 7) =
      -(copt2316 * copt3113 * copt774 * copt793) / 2. -
      (-2 * copt114 * copt158 + copt2335) * copt3123 * copt775 * copt793 +
      copt774 * copt775 *
          (-(copt122 * copt2454 * copt326 * copt778) -
           copt114 * copt1144 * copt325 * copt326 * copt778 -
           copt109 * copt2419 * copt708 * copt782 -
           copt133 * copt2393 * copt254 * copt787 -
           copt125 * copt1468 * copt250 * copt254 * copt787 -
           2 * copt109 * copt19 * copt699 * copt708 * copt781 * copt9 -
           2 * copt122 * copt325 * copt326 * copt777 *
               (copt20 - 2 * copt114 * copt9) -
           2 * copt133 * copt250 * copt254 * copt786 *
               (copt818 - copt3057 * copt9));
  out2(5, 8) =
      -(copt2483 * copt3113 * copt774 * copt793) / 2. -
      (-2 * copt117 * copt158 + copt2485) * copt3123 * copt775 * copt793 +
      copt774 * copt775 *
          (-(copt122 * copt2644 * copt326 * copt778) -
           copt1144 * copt117 * copt325 * copt326 * copt778 -
           copt109 * copt2594 * copt708 * copt782 -
           copt128 * copt1468 * copt250 * copt254 * copt787 -
           copt133 * copt254 * copt2549 * copt787 -
           2 * copt109 * copt30 * copt699 * copt708 * copt781 * copt9 -
           2 * copt122 * copt325 * copt326 * copt777 *
               (copt32 - 2 * copt117 * copt9) -
           2 * copt133 * copt250 * copt254 * copt786 *
               (copt841 - copt3085 * copt9));
  out2(5, 9)  = -(copt109 * copt2682 * copt708 * copt774 * copt775 * copt782);
  out2(5, 10) = -(copt109 * copt2693 * copt708 * copt774 * copt775 * copt782);
  out2(5, 11) = -(copt109 * copt2701 * copt708 * copt774 * copt775 * copt782);
  out2(5, 12) = -(copt133 * copt254 * copt2715 * copt774 * copt775 * copt787);
  out2(5, 13) = -(copt133 * copt254 * copt2759 * copt774 * copt775 * copt787);
  out2(5, 14) = -(copt133 * copt254 * copt2770 * copt774 * copt775 * copt787);
  out2(5, 15) = -(copt122 * copt2780 * copt326 * copt774 * copt775 * copt778);
  out2(5, 16) = -(copt122 * copt2789 * copt326 * copt774 * copt775 * copt778);
  out2(5, 17) = -(copt122 * copt2797 * copt326 * copt774 * copt775 * copt778);

  return std::make_tuple(grad, val);
}

#endif  // hylc_strain_II
