#include "strain.hpp"
#ifdef hylc_strain_II

using namespace hylc;
using namespace hylc::mathematica;

std::tuple<Mat6x18, Vec6> hylc::mathematica::strain_valgrad(
    const Vec18 &xloc, const Mat2x2 &invDm, const Vec3 &niexists) {
  // define output
  Mat6x18 grad(0);
  Vec6 val(0);
  auto out1 = [&](int i) -> Real & { return val[i]; };
  auto out2 = [&](int i, int j) -> Real & { return grad(i, j); };

  Real copt41  = invDm(0, 0);
  Real copt45  = xloc(0);
  Real copt49  = -copt45;
  Real copt50  = xloc(3);
  Real copt57  = copt49 + copt50;
  Real copt158 = copt41 * copt57;
  Real copt236 = invDm(1, 0);
  Real copt297 = xloc(6);
  Real copt298 = copt297 + copt49;
  Real copt304 = copt236 * copt298;
  Real copt325 = copt158 + copt304;
  Real copt331 = Power(copt325, 2);
  Real copt332 = xloc(1);
  Real copt333 = -copt332;
  Real copt334 = xloc(4);
  Real copt335 = copt333 + copt334;
  Real copt336 = copt335 * copt41;
  Real copt338 = xloc(7);
  Real copt340 = copt333 + copt338;
  Real copt342 = copt236 * copt340;
  Real copt343 = copt336 + copt342;
  Real copt344 = Power(copt343, 2);
  Real copt345 = xloc(2);
  Real copt346 = -copt345;
  Real copt347 = xloc(5);
  Real copt348 = copt346 + copt347;
  Real copt349 = copt348 * copt41;
  Real copt350 = xloc(8);
  Real copt351 = copt346 + copt350;
  Real copt352 = copt236 * copt351;
  Real copt353 = copt349 + copt352;
  Real copt354 = Power(copt353, 2);
  Real copt355 = copt331 + copt344 + copt354;
  Real copt356 = Sqrt(copt355);
  Real copt357 = 1 / copt356;
  Real copt358 = invDm(0, 1);
  Real copt360 = invDm(1, 1);
  Real copt359 = copt358 * copt57;
  Real copt361 = copt298 * copt360;
  Real copt362 = copt359 + copt361;
  Real copt373 = Power(copt362, 2);
  Real copt364 = copt335 * copt358;
  Real copt365 = copt340 * copt360;
  Real copt366 = copt364 + copt365;
  Real copt374 = Power(copt366, 2);
  Real copt368 = copt348 * copt358;
  Real copt369 = copt351 * copt360;
  Real copt370 = copt368 + copt369;
  Real copt375 = Power(copt370, 2);
  Real copt376 = copt373 + copt374 + copt375;
  Real copt377 = Sqrt(copt376);
  Real copt378 = 1 / copt377;
  Real copt389 = xloc(9);
  Real copt394 = xloc(10);
  Real copt381 = -copt50;
  Real copt413 = xloc(11);
  Real copt396 = -copt334;
  Real copt406 = -copt347;
  Real copt438 = copt381 + copt45;
  Real copt439 = Power(copt438, 2);
  Real copt440 = copt332 + copt396;
  Real copt441 = Power(copt440, 2);
  Real copt442 = copt345 + copt406;
  Real copt443 = Power(copt442, 2);
  Real copt444 = copt439 + copt441 + copt443;
  Real copt445 = Sqrt(copt444);
  Real copt401 = -copt297;
  Real copt428 = -copt394;
  Real copt390 = -copt389;
  Real copt423 = -copt350;
  Real copt477 = Power(copt236, 2);
  Real copt479 = 1 / copt445;
  Real copt402 = copt401 + copt50;
  Real copt421 = copt338 + copt396;
  Real copt486 = xloc(15);
  Real copt384 = -copt338;
  Real copt491 = xloc(16);
  Real copt403 = copt345 * copt402;
  Real copt404 = copt297 * copt347;
  Real copt405 = -(copt350 * copt50);
  Real copt407 = copt350 + copt406;
  Real copt408 = copt407 * copt45;
  Real copt409 = copt403 + copt404 + copt405 + copt408;
  Real copt487 = -copt486;
  Real copt488 = copt297 + copt487;
  Real copt499 = xloc(17);
  Real copt385 = copt334 + copt384;
  Real copt501 = copt423 + copt499;
  Real copt519 = copt401 + copt45;
  Real copt520 = Power(copt519, 2);
  Real copt521 = copt332 + copt384;
  Real copt522 = Power(copt521, 2);
  Real copt523 = copt345 + copt423;
  Real copt524 = Power(copt523, 2);
  Real copt525 = copt520 + copt522 + copt524;
  Real copt526 = Sqrt(copt525);
  Real copt490 = copt338 * copt486;
  Real copt492 = -(copt297 * copt491);
  Real copt493 = copt384 + copt491;
  Real copt555 = Power(copt41, 2);
  Real copt557 = 1 / copt526;
  Real copt481 = copt332 * copt402;
  Real copt482 = copt297 * copt334;
  Real copt483 = -(copt338 * copt50);
  Real copt484 = copt421 * copt45;
  Real copt485 = copt481 + copt482 + copt483 + copt484;
  Real copt559 = xloc(12);
  Real copt563 = xloc(13);
  Real copt561 = copt401 + copt559;
  Real copt572 = xloc(14);
  Real copt505 = copt345 * copt385;
  Real copt506 = copt338 * copt347;
  Real copt507 = -(copt334 * copt350);
  Real copt508 = copt332 * copt407;
  Real copt509 = copt505 + copt506 + copt507 + copt508;
  Real copt573 = -copt572;
  Real copt574 = copt350 + copt573;
  Real copt587 = Power(copt402, 2);
  Real copt588 = Power(copt385, 2);
  Real copt424 = copt347 + copt423;
  Real copt589 = Power(copt424, 2);
  Real copt590 = copt587 + copt588 + copt589;
  Real copt591 = Sqrt(copt590);
  Real copt446 = -(copt338 * copt347 * copt45);
  Real copt447 = copt334 * copt350 * copt45;
  Real copt560 = -(copt338 * copt559);
  Real copt562 = copt334 * copt561;
  Real copt564 = -copt563;
  Real copt565 = copt338 + copt564;
  Real copt566 = copt50 * copt565;
  Real copt567 = copt297 * copt563;
  Real copt568 = copt560 + copt562 + copt566 + copt567;
  Real copt615 = copt236 + copt41;
  Real copt616 = Power(copt615, 2);
  Real copt618 = 1 / copt591;
  Real copt621 = Power(copt485, 2);
  Real copt622 = Power(copt409, 2);
  Real copt623 = Power(copt509, 2);
  Real copt624 = copt621 + copt622 + copt623;
  Real copt625 = Sqrt(copt624);
  Real copt569 = copt485 * copt568;
  Real copt570 = -(copt350 * copt559);
  Real copt571 = copt347 * copt561;
  Real copt575 = copt50 * copt574;
  Real copt576 = copt297 * copt572;
  Real copt577 = copt570 + copt571 + copt575 + copt576;
  Real copt578 = copt409 * copt577;
  Real copt579 = -(copt350 * copt563);
  Real copt580 = copt384 + copt563;
  Real copt581 = copt347 * copt580;
  Real copt582 = copt334 * copt574;
  Real copt583 = copt338 * copt572;
  Real copt584 = copt579 + copt581 + copt582 + copt583;
  Real copt585 = copt509 * copt584;
  Real copt586 = copt569 + copt578 + copt585;
  Real copt592 = copt338 * copt347 * copt559;
  Real copt593 = -(copt334 * copt350 * copt559);
  Real copt594 = copt347 * copt45 * copt563;
  Real copt595 = -(copt297 * copt347 * copt563);
  Real copt596 = -(copt350 * copt45 * copt563);
  Real copt597 = copt350 * copt50 * copt563;
  Real copt598 = copt345 * copt568;
  Real copt599 = -(copt334 * copt45 * copt572);
  Real copt600 = copt297 * copt334 * copt572;
  Real copt601 = copt338 * copt45 * copt572;
  Real copt602 = -(copt338 * copt50 * copt572);
  Real copt603 = -copt559;
  Real copt604 = copt297 + copt603;
  Real copt605 = copt347 * copt604;
  Real copt606 = copt350 * copt559;
  Real copt607 = -(copt297 * copt572);
  Real copt608 = copt423 + copt572;
  Real copt609 = copt50 * copt608;
  Real copt610 = copt605 + copt606 + copt607 + copt609;
  Real copt611 = copt332 * copt610;
  Real copt612 = copt446 + copt447 + copt592 + copt593 + copt594 + copt595 +
                 copt596 + copt597 + copt598 + copt599 + copt600 + copt601 +
                 copt602 + copt611;
  Real copt613 = copt591 * copt612;
  Real copt614 = ArcTan(copt586, copt613);
  Real copt617 = niexists(1);
  Real copt628 = Power(copt45, 2);
  Real copt629 = Power(copt332, 2);
  Real copt630 = Power(copt345, 2);
  Real copt631 = -2 * copt45 * copt50;
  Real copt632 = Power(copt50, 2);
  Real copt633 = -2 * copt332 * copt334;
  Real copt634 = Power(copt334, 2);
  Real copt635 = -2 * copt345 * copt347;
  Real copt636 = Power(copt347, 2);
  Real copt637 = copt628 + copt629 + copt630 + copt631 + copt632 + copt633 +
                 copt634 + copt635 + copt636;
  Real copt638 = Sqrt(copt637);
  Real copt639 = -2 * copt297 * copt45;
  Real copt640 = Power(copt297, 2);
  Real copt641 = -2 * copt332 * copt338;
  Real copt642 = Power(copt338, 2);
  Real copt643 = -2 * copt345 * copt350;
  Real copt644 = Power(copt350, 2);
  Real copt645 = copt628 + copt629 + copt630 + copt639 + copt640 + copt641 +
                 copt642 + copt643 + copt644;
  Real copt646 = Sqrt(copt645);
  Real copt648 = -2 * copt297 * copt50;
  Real copt649 = -2 * copt334 * copt338;
  Real copt650 = -2 * copt347 * copt350;
  Real copt651 = copt632 + copt634 + copt636 + copt640 + copt642 + copt644 +
                 copt648 + copt649 + copt650;
  Real copt652 = Sqrt(copt651);
  Real copt489 = copt332 * copt488;
  Real copt494 = copt45 * copt493;
  Real copt495 = copt489 + copt490 + copt492 + copt494;
  Real copt496 = copt485 * copt495;
  Real copt497 = copt345 * copt488;
  Real copt498 = copt350 * copt486;
  Real copt500 = -(copt297 * copt499);
  Real copt502 = copt45 * copt501;
  Real copt503 = copt497 + copt498 + copt500 + copt502;
  Real copt504 = copt409 * copt503;
  Real copt510 = -copt491;
  Real copt511 = copt338 + copt510;
  Real copt512 = copt345 * copt511;
  Real copt513 = copt350 * copt491;
  Real copt514 = -(copt338 * copt499);
  Real copt515 = copt332 * copt501;
  Real copt516 = copt512 + copt513 + copt514 + copt515;
  Real copt517 = copt509 * copt516;
  Real copt518 = copt496 + copt504 + copt517;
  Real copt527 = copt338 * copt347 * copt45;
  Real copt528 = -(copt334 * copt350 * copt45);
  Real copt529 = -(copt338 * copt347 * copt486);
  Real copt530 = copt334 * copt350 * copt486;
  Real copt531 = -(copt347 * copt45 * copt491);
  Real copt532 = copt297 * copt347 * copt491;
  Real copt533 = copt350 * copt45 * copt491;
  Real copt534 = -(copt350 * copt491 * copt50);
  Real copt535 = copt334 * copt488;
  Real copt536 = copt493 * copt50;
  Real copt537 = copt490 + copt492 + copt535 + copt536;
  Real copt538 = copt345 * copt537;
  Real copt539 = copt334 * copt45 * copt499;
  Real copt540 = -(copt297 * copt334 * copt499);
  Real copt541 = -(copt338 * copt45 * copt499);
  Real copt542 = copt338 * copt499 * copt50;
  Real copt543 = -(copt350 * copt486);
  Real copt544 = copt401 + copt486;
  Real copt545 = copt347 * copt544;
  Real copt546 = -copt499;
  Real copt547 = copt350 + copt546;
  Real copt548 = copt50 * copt547;
  Real copt549 = copt297 * copt499;
  Real copt550 = copt543 + copt545 + copt548 + copt549;
  Real copt551 = copt332 * copt550;
  Real copt552 = copt527 + copt528 + copt529 + copt530 + copt531 + copt532 +
                 copt533 + copt534 + copt538 + copt539 + copt540 + copt541 +
                 copt542 + copt551;
  Real copt553 = -(copt526 * copt552);
  Real copt554 = ArcTan(copt518, copt553);
  Real copt556 = niexists(2);
  Real copt380 = -(copt297 * copt334);
  Real copt382 = copt297 + copt381;
  Real copt383 = copt332 * copt382;
  Real copt386 = copt385 * copt45;
  Real copt387 = copt338 * copt50;
  Real copt388 = copt380 + copt383 + copt386 + copt387;
  Real copt391 = copt390 + copt50;
  Real copt392 = copt332 * copt391;
  Real copt393 = copt334 * copt389;
  Real copt395 = -(copt394 * copt50);
  Real copt397 = copt394 + copt396;
  Real copt398 = copt397 * copt45;
  Real copt399 = copt392 + copt393 + copt395 + copt398;
  Real copt400 = copt388 * copt399;
  Real copt410 = -(copt347 * copt389);
  Real copt411 = copt381 + copt389;
  Real copt412 = copt345 * copt411;
  Real copt414 = -copt413;
  Real copt415 = copt347 + copt414;
  Real copt416 = copt415 * copt45;
  Real copt417 = copt413 * copt50;
  Real copt418 = copt410 + copt412 + copt416 + copt417;
  Real copt419 = copt409 * copt418;
  Real copt420 = -(copt338 * copt347);
  Real copt422 = copt345 * copt421;
  Real copt425 = copt332 * copt424;
  Real copt426 = copt334 * copt350;
  Real copt427 = copt420 + copt422 + copt425 + copt426;
  Real copt429 = copt334 + copt428;
  Real copt430 = copt345 * copt429;
  Real copt431 = copt347 * copt394;
  Real copt432 = -(copt334 * copt413);
  Real copt433 = copt406 + copt413;
  Real copt434 = copt332 * copt433;
  Real copt435 = copt430 + copt431 + copt432 + copt434;
  Real copt436 = copt427 * copt435;
  Real copt437 = copt400 + copt419 + copt436;
  Real copt448 = copt338 * copt347 * copt389;
  Real copt449 = -(copt334 * copt350 * copt389);
  Real copt450 = copt347 * copt394 * copt45;
  Real copt451 = -(copt297 * copt347 * copt394);
  Real copt452 = -(copt350 * copt394 * copt45);
  Real copt453 = copt350 * copt394 * copt50;
  Real copt454 = -(copt338 * copt389);
  Real copt455 = copt389 + copt401;
  Real copt456 = copt334 * copt455;
  Real copt457 = copt338 + copt428;
  Real copt458 = copt457 * copt50;
  Real copt459 = copt297 * copt394;
  Real copt460 = copt454 + copt456 + copt458 + copt459;
  Real copt461 = copt345 * copt460;
  Real copt462 = -(copt334 * copt413 * copt45);
  Real copt463 = copt297 * copt334 * copt413;
  Real copt464 = copt338 * copt413 * copt45;
  Real copt465 = -(copt338 * copt413 * copt50);
  Real copt466 = copt297 + copt390;
  Real copt467 = copt347 * copt466;
  Real copt468 = copt350 * copt389;
  Real copt469 = -(copt297 * copt413);
  Real copt470 = copt413 + copt423;
  Real copt471 = copt470 * copt50;
  Real copt472 = copt467 + copt468 + copt469 + copt471;
  Real copt473 = copt332 * copt472;
  Real copt474 = copt446 + copt447 + copt448 + copt449 + copt450 + copt451 +
                 copt452 + copt453 + copt461 + copt462 + copt463 + copt464 +
                 copt465 + copt473;
  Real copt475  = copt445 * copt474;
  Real copt476  = ArcTan(copt437, copt475);
  Real copt478  = niexists(0);
  Real copt659  = Power(copt360, 2);
  Real copt661  = Power(copt358, 2);
  Real copt627  = copt358 + copt360;
  Real copt663  = Power(copt627, 2);
  Real copt685  = copt355 * copt356;
  Real copt686  = 1 / copt685;
  Real copt687  = copt376 * copt377;
  Real copt688  = 1 / copt687;
  Real copt363  = copt325 * copt362;
  Real copt367  = copt343 * copt366;
  Real copt371  = copt353 * copt370;
  Real copt372  = copt363 + copt367 + copt371;
  Real copt667  = copt41 * copt438;
  Real copt668  = copt236 * copt519;
  Real copt669  = copt667 + copt668;
  Real copt671  = copt41 * copt440;
  Real copt672  = copt236 * copt521;
  Real copt673  = copt671 + copt672;
  Real copt675  = copt41 * copt442;
  Real copt676  = copt236 * copt523;
  Real copt677  = copt675 + copt676;
  Real copt777  = -copt358;
  Real copt778  = -copt360;
  Real copt779  = copt777 + copt778;
  Real copt480  = -(copt476 * copt477 * copt478 * copt479);
  Real copt558  = -(copt554 * copt555 * copt556 * copt557);
  Real copt619  = -(copt614 * copt616 * copt617 * copt618);
  Real copt620  = copt480 + copt558 + copt619;
  Real copt792  = 1 / copt625;
  Real copt794  = copt444 * copt445;
  Real copt795  = 1 / copt794;
  Real copt797  = copt525 * copt526;
  Real copt798  = 1 / copt797;
  Real copt805  = Power(copt437, 2);
  Real copt806  = Power(copt474, 2);
  Real copt807  = copt444 * copt806;
  Real copt808  = copt805 + copt807;
  Real copt809  = 1 / copt808;
  Real copt823  = Power(copt586, 2);
  Real copt824  = Power(copt612, 2);
  Real copt825  = copt590 * copt824;
  Real copt826  = copt823 + copt825;
  Real copt827  = 1 / copt826;
  Real copt840  = Power(copt552, 2);
  Real copt841  = copt525 * copt840;
  Real copt842  = Power(copt518, 2);
  Real copt843  = copt841 + copt842;
  Real copt844  = 1 / copt843;
  Real copt943  = copt590 * copt591;
  Real copt944  = 1 / copt943;
  Real copt1011 = copt350 * copt45;
  Real copt1065 = -(copt338 * copt45);
  Real copt1128 = copt396 + copt563;
  Real copt1131 = copt406 + copt572;
  Real copt1178 = -(copt347 * copt45);
  Real copt1191 = copt50 + copt603;
  Real copt1151 = copt345 + copt546;
  Real copt1238 = copt334 * copt45;
  Real copt1209 = copt486 + copt49;
  Real copt1216 = copt347 * copt45;
  Real copt1042 = -(copt350 * copt45);
  Real copt1277 = -(copt334 * copt45);
  Real copt1098 = copt338 * copt45;
  Real copt1300 = -(copt297 * copt347);
  Real copt1301 = copt345 * copt382;
  Real copt1302 = copt350 * copt50;
  Real copt1303 = copt1042 + copt1216 + copt1300 + copt1301 + copt1302;
  Real copt1311 = copt1098 + copt1277 + copt481 + copt482 + copt483;
  Real copt789  = 2 * copt421 * copt485;
  Real copt790  = 2 * copt407 * copt409;
  Real copt791  = copt789 + copt790;
  Real copt647  = copt614 * copt615 * copt617 * copt627 * copt638 * copt646;
  Real copt653  = copt358 * copt41 * copt554 * copt556 * copt638;
  Real copt654  = copt236 * copt360 * copt476 * copt478 * copt646;
  Real copt655  = copt653 + copt654;
  Real copt656  = copt652 * copt655;
  Real copt657  = copt647 + copt656;
  Real copt1369 = 1 / copt646;
  Real copt1366 = 2 * copt45;
  Real copt1373 = 1 / copt638;
  Real copt820  = copt347 * copt563;
  Real copt821  = -(copt334 * copt572);
  Real copt822  = copt420 + copt426 + copt579 + copt583 + copt820 + copt821;
  Real copt828  = copt586 * copt591 * copt822 * copt827;
  Real copt829  = copt421 * copt568;
  Real copt830  = copt407 * copt577;
  Real copt831  = copt829 + copt830;
  Real copt832  = -(copt591 * copt612 * copt827 * copt831);
  Real copt833  = copt828 + copt832;
  Real copt1371 = -2 * copt50;
  Real copt1372 = copt1366 + copt1371;
  Real copt1367 = -2 * copt297;
  Real copt1368 = copt1366 + copt1367;
  Real copt800  = copt388 * copt397;
  Real copt801  = copt385 * copt399;
  Real copt802  = copt409 * copt415;
  Real copt803  = copt407 * copt418;
  Real copt804  = copt800 + copt801 + copt802 + copt803;
  Real copt810  = -(copt445 * copt474 * copt804 * copt809);
  Real copt811  = -(copt350 * copt394);
  Real copt812  = copt338 * copt413;
  Real copt813  = copt420 + copt426 + copt431 + copt432 + copt811 + copt812;
  Real copt814  = copt445 * copt813;
  Real copt815  = copt438 * copt474 * copt479;
  Real copt816  = copt814 + copt815;
  Real copt817  = copt437 * copt809 * copt816;
  Real copt818  = copt810 + copt817;
  Real copt835  = copt485 * copt493;
  Real copt836  = copt421 * copt495;
  Real copt837  = copt409 * copt501;
  Real copt838  = copt407 * copt503;
  Real copt839  = copt835 + copt836 + copt837 + copt838;
  Real copt845  = copt526 * copt552 * copt839 * copt844;
  Real copt846  = -(copt347 * copt491);
  Real copt847  = copt334 * copt499;
  Real copt848  = copt506 + copt507 + copt513 + copt514 + copt846 + copt847;
  Real copt849  = -(copt526 * copt848);
  Real copt850  = -(copt519 * copt552 * copt557);
  Real copt851  = copt849 + copt850;
  Real copt852  = copt518 * copt844 * copt851;
  Real copt853  = copt845 + copt852;
  Real copt858  = 2 * copt402 * copt485;
  Real copt859  = 2 * copt407 * copt509;
  Real copt860  = copt858 + copt859;
  Real copt1388 = 2 * copt332;
  Real copt876  = copt586 * copt591 * copt610 * copt827;
  Real copt877  = copt402 * copt568;
  Real copt878  = copt407 * copt584;
  Real copt879  = copt877 + copt878;
  Real copt880  = -(copt591 * copt612 * copt827 * copt879);
  Real copt881  = copt876 + copt880;
  Real copt1392 = -2 * copt334;
  Real copt1393 = copt1388 + copt1392;
  Real copt1389 = -2 * copt338;
  Real copt1390 = copt1388 + copt1389;
  Real copt864  = copt388 * copt391;
  Real copt865  = copt382 * copt399;
  Real copt866  = copt427 * copt433;
  Real copt867  = copt424 * copt435;
  Real copt868  = copt864 + copt865 + copt866 + copt867;
  Real copt869  = -(copt445 * copt474 * copt809 * copt868);
  Real copt870  = copt445 * copt472;
  Real copt871  = copt440 * copt474 * copt479;
  Real copt872  = copt870 + copt871;
  Real copt873  = copt437 * copt809 * copt872;
  Real copt874  = copt869 + copt873;
  Real copt883  = copt485 * copt488;
  Real copt884  = copt402 * copt495;
  Real copt885  = copt501 * copt509;
  Real copt886  = copt407 * copt516;
  Real copt887  = copt883 + copt884 + copt885 + copt886;
  Real copt888  = copt526 * copt552 * copt844 * copt887;
  Real copt889  = -(copt526 * copt550);
  Real copt890  = -(copt521 * copt552 * copt557);
  Real copt891  = copt889 + copt890;
  Real copt892  = copt518 * copt844 * copt891;
  Real copt893  = copt888 + copt892;
  Real copt898  = 2 * copt402 * copt409;
  Real copt899  = 2 * copt385 * copt509;
  Real copt900  = copt898 + copt899;
  Real copt1408 = 2 * copt345;
  Real copt916  = copt568 * copt586 * copt591 * copt827;
  Real copt917  = copt402 * copt577;
  Real copt918  = copt385 * copt584;
  Real copt919  = copt917 + copt918;
  Real copt920  = -(copt591 * copt612 * copt827 * copt919);
  Real copt921  = copt916 + copt920;
  Real copt1412 = -2 * copt347;
  Real copt1413 = copt1408 + copt1412;
  Real copt1409 = -2 * copt350;
  Real copt1410 = copt1408 + copt1409;
  Real copt904  = copt409 * copt411;
  Real copt905  = copt427 * copt429;
  Real copt906  = copt402 * copt418;
  Real copt907  = copt421 * copt435;
  Real copt908  = copt904 + copt905 + copt906 + copt907;
  Real copt909  = -(copt445 * copt474 * copt809 * copt908);
  Real copt910  = copt445 * copt460;
  Real copt911  = copt442 * copt474 * copt479;
  Real copt912  = copt910 + copt911;
  Real copt913  = copt437 * copt809 * copt912;
  Real copt914  = copt909 + copt913;
  Real copt923  = copt409 * copt488;
  Real copt924  = copt509 * copt511;
  Real copt925  = copt402 * copt503;
  Real copt926  = copt385 * copt516;
  Real copt927  = copt923 + copt924 + copt925 + copt926;
  Real copt928  = copt526 * copt552 * copt844 * copt927;
  Real copt929  = -(copt526 * copt537);
  Real copt930  = -(copt523 * copt552 * copt557);
  Real copt931  = copt929 + copt930;
  Real copt932  = copt518 * copt844 * copt931;
  Real copt933  = copt928 + copt932;
  Real copt938  = 2 * copt485 * copt521;
  Real copt939  = 2 * copt409 * copt523;
  Real copt940  = copt938 + copt939;
  Real copt1429 = 2 * copt50;
  Real copt1433 = 1 / copt652;
  Real copt965  = copt485 * copt565;
  Real copt966  = copt521 * copt568;
  Real copt967  = copt409 * copt574;
  Real copt968  = copt523 * copt577;
  Real copt969  = copt965 + copt966 + copt967 + copt968;
  Real copt970  = -(copt591 * copt612 * copt827 * copt969);
  Real copt971  = copt345 * copt565;
  Real copt972  = copt350 * copt563;
  Real copt973  = -(copt338 * copt572);
  Real copt974  = copt332 * copt608;
  Real copt975  = copt971 + copt972 + copt973 + copt974;
  Real copt976  = copt591 * copt975;
  Real copt977  = copt402 * copt612 * copt618;
  Real copt978  = copt976 + copt977;
  Real copt979  = copt586 * copt827 * copt978;
  Real copt980  = copt970 + copt979;
  Real copt1428 = -2 * copt45;
  Real copt1430 = copt1428 + copt1429;
  Real copt946  = copt332 + copt428;
  Real copt947  = copt388 * copt946;
  Real copt948  = copt340 * copt399;
  Real copt949  = copt346 + copt413;
  Real copt950  = copt409 * copt949;
  Real copt951  = copt418 * copt523;
  Real copt952  = copt947 + copt948 + copt950 + copt951;
  Real copt953  = -(copt445 * copt474 * copt809 * copt952);
  Real copt954  = copt345 * copt457;
  Real copt955  = copt350 * copt394;
  Real copt956  = -(copt338 * copt413);
  Real copt957  = copt332 * copt470;
  Real copt958  = copt954 + copt955 + copt956 + copt957;
  Real copt959  = copt445 * copt958;
  Real copt960  = -(copt438 * copt474 * copt479);
  Real copt961  = copt959 + copt960;
  Real copt962  = copt437 * copt809 * copt961;
  Real copt963  = copt953 + copt962;
  Real copt982  = copt495 * copt521;
  Real copt983  = copt503 * copt523;
  Real copt984  = copt982 + copt983;
  Real copt985  = copt526 * copt552 * copt844 * copt984;
  Real copt986  = -(copt350 * copt491);
  Real copt987  = copt345 * copt493;
  Real copt988  = copt332 * copt547;
  Real copt989  = copt338 * copt499;
  Real copt990  = copt986 + copt987 + copt988 + copt989;
  Real copt991  = -(copt518 * copt526 * copt844 * copt990);
  Real copt992  = copt985 + copt991;
  Real copt997  = 2 * copt298 * copt485;
  Real copt998  = 2 * copt509 * copt523;
  Real copt999  = copt997 + copt998;
  Real copt1448 = 2 * copt334;
  Real copt1023 = copt485 * copt561;
  Real copt1024 = copt298 * copt568;
  Real copt1025 = copt509 * copt574;
  Real copt1026 = copt523 * copt584;
  Real copt1027 = copt1023 + copt1024 + copt1025 + copt1026;
  Real copt1028 = -(copt1027 * copt591 * copt612 * copt827);
  Real copt1029 = copt345 * copt561;
  Real copt1030 = -(copt45 * copt572);
  Real copt1031 = copt1011 + copt1029 + copt1030 + copt570 + copt576;
  Real copt1032 = copt1031 * copt591;
  Real copt1033 = copt385 * copt612 * copt618;
  Real copt1034 = copt1032 + copt1033;
  Real copt1035 = copt1034 * copt586 * copt827;
  Real copt1036 = copt1028 + copt1035;
  Real copt1447 = -2 * copt332;
  Real copt1449 = copt1447 + copt1448;
  Real copt1003 = copt389 + copt49;
  Real copt1004 = copt1003 * copt388;
  Real copt1005 = copt399 * copt519;
  Real copt1006 = copt345 + copt414;
  Real copt1007 = copt1006 * copt427;
  Real copt1008 = copt351 * copt435;
  Real copt1009 = copt1004 + copt1005 + copt1007 + copt1008;
  Real copt1010 = -(copt1009 * copt445 * copt474 * copt809);
  Real copt1012 = -(copt350 * copt389);
  Real copt1013 = copt345 * copt455;
  Real copt1014 = -(copt413 * copt45);
  Real copt1015 = copt297 * copt413;
  Real copt1016 = copt1011 + copt1012 + copt1013 + copt1014 + copt1015;
  Real copt1017 = copt1016 * copt445;
  Real copt1018 = -(copt440 * copt474 * copt479);
  Real copt1019 = copt1017 + copt1018;
  Real copt1020 = copt1019 * copt437 * copt809;
  Real copt1021 = copt1010 + copt1020;
  Real copt1038 = copt298 * copt495;
  Real copt1039 = copt516 * copt523;
  Real copt1040 = copt1038 + copt1039;
  Real copt1041 = copt1040 * copt526 * copt552 * copt844;
  Real copt1043 = copt45 * copt499;
  Real copt1044 = copt1042 + copt1043 + copt497 + copt498 + copt500;
  Real copt1045 = -(copt1044 * copt518 * copt526 * copt844);
  Real copt1046 = copt1041 + copt1045;
  Real copt1051 = 2 * copt298 * copt409;
  Real copt1052 = 2 * copt340 * copt509;
  Real copt1053 = copt1051 + copt1052;
  Real copt1466 = 2 * copt347;
  Real copt1077 = copt409 * copt561;
  Real copt1078 = copt509 * copt580;
  Real copt1079 = copt298 * copt577;
  Real copt1080 = copt340 * copt584;
  Real copt1081 = copt1077 + copt1078 + copt1079 + copt1080;
  Real copt1082 = -(copt1081 * copt591 * copt612 * copt827);
  Real copt1083 = copt332 * copt604;
  Real copt1084 = copt338 * copt559;
  Real copt1085 = copt45 * copt563;
  Real copt1086 = -(copt297 * copt563);
  Real copt1087 = copt1065 + copt1083 + copt1084 + copt1085 + copt1086;
  Real copt1088 = copt1087 * copt591;
  Real copt1089 = copt424 * copt612 * copt618;
  Real copt1090 = copt1088 + copt1089;
  Real copt1091 = copt1090 * copt586 * copt827;
  Real copt1092 = copt1082 + copt1091;
  Real copt1465 = -2 * copt345;
  Real copt1467 = copt1465 + copt1466;
  Real copt1057 = copt390 + copt45;
  Real copt1058 = copt1057 * copt409;
  Real copt1059 = copt333 + copt394;
  Real copt1060 = copt1059 * copt427;
  Real copt1061 = copt298 * copt418;
  Real copt1062 = copt435 * copt521;
  Real copt1063 = copt1058 + copt1060 + copt1061 + copt1062;
  Real copt1064 = -(copt1063 * copt445 * copt474 * copt809);
  Real copt1066 = copt332 * copt466;
  Real copt1067 = copt338 * copt389;
  Real copt1068 = copt394 * copt45;
  Real copt1069 = -(copt297 * copt394);
  Real copt1070 = copt1065 + copt1066 + copt1067 + copt1068 + copt1069;
  Real copt1071 = copt1070 * copt445;
  Real copt1072 = -(copt442 * copt474 * copt479);
  Real copt1073 = copt1071 + copt1072;
  Real copt1074 = copt1073 * copt437 * copt809;
  Real copt1075 = copt1064 + copt1074;
  Real copt1094 = copt298 * copt503;
  Real copt1095 = copt340 * copt516;
  Real copt1096 = copt1094 + copt1095;
  Real copt1097 = copt1096 * copt526 * copt552 * copt844;
  Real copt1099 = -(copt338 * copt486);
  Real copt1100 = copt332 * copt544;
  Real copt1101 = -(copt45 * copt491);
  Real copt1102 = copt297 * copt491;
  Real copt1103 = copt1098 + copt1099 + copt1100 + copt1101 + copt1102;
  Real copt1104 = -(copt1103 * copt518 * copt526 * copt844);
  Real copt1105 = copt1097 + copt1104;
  Real copt1110 = 2 * copt335 * copt485;
  Real copt1111 = 2 * copt348 * copt409;
  Real copt1112 = copt1110 + copt1111;
  Real copt1483 = 2 * copt297;
  Real copt1129 = copt1128 * copt485;
  Real copt1130 = copt335 * copt568;
  Real copt1132 = copt1131 * copt409;
  Real copt1133 = copt348 * copt577;
  Real copt1134 = copt1129 + copt1130 + copt1132 + copt1133;
  Real copt1135 = -(copt1134 * copt591 * copt612 * copt827);
  Real copt1136 = -(copt347 * copt563);
  Real copt1137 = copt1128 * copt345;
  Real copt1138 = copt347 + copt573;
  Real copt1139 = copt1138 * copt332;
  Real copt1140 = copt334 * copt572;
  Real copt1141 = copt1136 + copt1137 + copt1139 + copt1140;
  Real copt1142 = copt1141 * copt591;
  Real copt1143 = -(copt402 * copt612 * copt618);
  Real copt1144 = copt1142 + copt1143;
  Real copt1145 = copt1144 * copt586 * copt827;
  Real copt1146 = copt1135 + copt1145;
  Real copt1484 = copt1428 + copt1483;
  Real copt1116 = -(copt347 * copt394);
  Real copt1117 = copt345 * copt397;
  Real copt1118 = copt332 * copt415;
  Real copt1119 = copt334 * copt413;
  Real copt1120 = copt1116 + copt1117 + copt1118 + copt1119;
  Real copt1121 = copt1120 * copt437 * copt445 * copt809;
  Real copt1122 = copt399 * copt440;
  Real copt1123 = copt348 * copt418;
  Real copt1124 = copt1122 + copt1123;
  Real copt1125 = -(copt1124 * copt445 * copt474 * copt809);
  Real copt1126 = copt1121 + copt1125;
  Real copt1148 = copt332 + copt510;
  Real copt1149 = copt1148 * copt485;
  Real copt1150 = copt335 * copt495;
  Real copt1152 = copt1151 * copt409;
  Real copt1153 = copt348 * copt503;
  Real copt1154 = copt1149 + copt1150 + copt1152 + copt1153;
  Real copt1155 = copt1154 * copt526 * copt552 * copt844;
  Real copt1156 = copt334 + copt510;
  Real copt1157 = copt1156 * copt345;
  Real copt1158 = copt347 * copt491;
  Real copt1159 = -(copt334 * copt499);
  Real copt1160 = copt406 + copt499;
  Real copt1161 = copt1160 * copt332;
  Real copt1162 = copt1157 + copt1158 + copt1159 + copt1161;
  Real copt1163 = -(copt1162 * copt526);
  Real copt1164 = copt519 * copt552 * copt557;
  Real copt1165 = copt1163 + copt1164;
  Real copt1166 = copt1165 * copt518 * copt844;
  Real copt1167 = copt1155 + copt1166;
  Real copt1172 = 2 * copt438 * copt485;
  Real copt1173 = 2 * copt348 * copt509;
  Real copt1174 = copt1172 + copt1173;
  Real copt1500 = 2 * copt338;
  Real copt1192 = copt1191 * copt485;
  Real copt1193 = copt438 * copt568;
  Real copt1194 = copt1131 * copt509;
  Real copt1195 = copt348 * copt584;
  Real copt1196 = copt1192 + copt1193 + copt1194 + copt1195;
  Real copt1197 = -(copt1196 * copt591 * copt612 * copt827);
  Real copt1198 = copt1191 * copt345;
  Real copt1199 = copt347 * copt559;
  Real copt1200 = copt45 * copt572;
  Real copt1201 = -(copt50 * copt572);
  Real copt1202 = copt1178 + copt1198 + copt1199 + copt1200 + copt1201;
  Real copt1203 = copt1202 * copt591;
  Real copt1204 = -(copt385 * copt612 * copt618);
  Real copt1205 = copt1203 + copt1204;
  Real copt1206 = copt1205 * copt586 * copt827;
  Real copt1207 = copt1197 + copt1206;
  Real copt1501 = copt1447 + copt1500;
  Real copt1179 = copt345 * copt391;
  Real copt1180 = copt347 * copt389;
  Real copt1181 = copt413 * copt45;
  Real copt1182 = -(copt413 * copt50);
  Real copt1183 = copt1178 + copt1179 + copt1180 + copt1181 + copt1182;
  Real copt1184 = copt1183 * copt437 * copt445 * copt809;
  Real copt1185 = copt399 * copt57;
  Real copt1186 = copt435 * copt442;
  Real copt1187 = copt1185 + copt1186;
  Real copt1188 = -(copt1187 * copt445 * copt474 * copt809);
  Real copt1189 = copt1184 + copt1188;
  Real copt1210 = copt1209 * copt485;
  Real copt1211 = copt438 * copt495;
  Real copt1212 = copt1151 * copt509;
  Real copt1213 = copt348 * copt516;
  Real copt1214 = copt1210 + copt1211 + copt1212 + copt1213;
  Real copt1215 = copt1214 * copt526 * copt552 * copt844;
  Real copt1217 = -(copt347 * copt486);
  Real copt1218 = copt381 + copt486;
  Real copt1219 = copt1218 * copt345;
  Real copt1220 = -(copt45 * copt499);
  Real copt1221 = copt499 * copt50;
  Real copt1222 = copt1216 + copt1217 + copt1219 + copt1220 + copt1221;
  Real copt1223 = -(copt1222 * copt526);
  Real copt1224 = copt521 * copt552 * copt557;
  Real copt1225 = copt1223 + copt1224;
  Real copt1226 = copt1225 * copt518 * copt844;
  Real copt1227 = copt1215 + copt1226;
  Real copt1232 = 2 * copt409 * copt438;
  Real copt1233 = 2 * copt440 * copt509;
  Real copt1234 = copt1232 + copt1233;
  Real copt1517 = 2 * copt350;
  Real copt1251 = copt1191 * copt409;
  Real copt1252 = copt334 + copt564;
  Real copt1253 = copt1252 * copt509;
  Real copt1254 = copt438 * copt577;
  Real copt1255 = copt440 * copt584;
  Real copt1256 = copt1251 + copt1253 + copt1254 + copt1255;
  Real copt1257 = -(copt1256 * copt591 * copt612 * copt827);
  Real copt1258 = -(copt334 * copt559);
  Real copt1259 = copt381 + copt559;
  Real copt1260 = copt1259 * copt332;
  Real copt1261 = -(copt45 * copt563);
  Real copt1262 = copt50 * copt563;
  Real copt1263 = copt1238 + copt1258 + copt1260 + copt1261 + copt1262;
  Real copt1264 = copt1263 * copt591;
  Real copt1265 = -(copt424 * copt612 * copt618);
  Real copt1266 = copt1264 + copt1265;
  Real copt1267 = copt1266 * copt586 * copt827;
  Real copt1268 = copt1257 + copt1267;
  Real copt1518 = copt1465 + copt1517;
  Real copt1239 = -(copt334 * copt389);
  Real copt1240 = copt332 * copt411;
  Real copt1241 = -(copt394 * copt45);
  Real copt1242 = copt394 * copt50;
  Real copt1243 = copt1238 + copt1239 + copt1240 + copt1241 + copt1242;
  Real copt1244 = copt1243 * copt437 * copt445 * copt809;
  Real copt1245 = copt418 * copt438;
  Real copt1246 = copt335 * copt435;
  Real copt1247 = copt1245 + copt1246;
  Real copt1248 = -(copt1247 * copt445 * copt474 * copt809);
  Real copt1249 = copt1244 + copt1248;
  Real copt1270 = copt1209 * copt409;
  Real copt1271 = copt333 + copt491;
  Real copt1272 = copt1271 * copt509;
  Real copt1273 = copt438 * copt503;
  Real copt1274 = copt440 * copt516;
  Real copt1275 = copt1270 + copt1272 + copt1273 + copt1274;
  Real copt1276 = copt1275 * copt526 * copt552 * copt844;
  Real copt1278 = copt487 + copt50;
  Real copt1279 = copt1278 * copt332;
  Real copt1280 = copt334 * copt486;
  Real copt1281 = copt45 * copt491;
  Real copt1282 = -(copt491 * copt50);
  Real copt1283 = copt1277 + copt1279 + copt1280 + copt1281 + copt1282;
  Real copt1284 = -(copt1283 * copt526);
  Real copt1285 = copt523 * copt552 * copt557;
  Real copt1286 = copt1284 + copt1285;
  Real copt1287 = copt1286 * copt518 * copt844;
  Real copt1288 = copt1276 + copt1287;
  Real copt1293 = copt437 * copt445 * copt509 * copt809;
  Real copt1294 = copt335 * copt388;
  Real copt1295 = copt409 * copt442;
  Real copt1296 = copt1294 + copt1295;
  Real copt1297 = -(copt1296 * copt445 * copt474 * copt809);
  Real copt1298 = copt1293 + copt1297;
  Real copt1304 = copt1303 * copt437 * copt445 * copt809;
  Real copt1305 = copt388 * copt438;
  Real copt1306 = copt348 * copt427;
  Real copt1307 = copt1305 + copt1306;
  Real copt1308 = -(copt1307 * copt445 * copt474 * copt809);
  Real copt1309 = copt1304 + copt1308;
  Real copt1312 = copt1311 * copt437 * copt445 * copt809;
  Real copt1313 = copt427 * copt440;
  Real copt1314 = copt409 * copt57;
  Real copt1315 = copt1313 + copt1314;
  Real copt1316 = -(copt1315 * copt445 * copt474 * copt809);
  Real copt1317 = copt1312 + copt1316;
  Real copt1319 = copt509 * copt586 * copt591 * copt827;
  Real copt1320 = copt385 * copt485;
  Real copt1321 = copt409 * copt424;
  Real copt1322 = copt1320 + copt1321;
  Real copt1323 = -(copt1322 * copt591 * copt612 * copt827);
  Real copt1324 = copt1319 + copt1323;
  Real copt1326 = copt1303 * copt586 * copt591 * copt827;
  Real copt1327 = copt382 * copt485;
  Real copt1328 = copt424 * copt509;
  Real copt1329 = copt1327 + copt1328;
  Real copt1330 = -(copt1329 * copt591 * copt612 * copt827);
  Real copt1331 = copt1326 + copt1330;
  Real copt1333 = copt1311 * copt586 * copt591 * copt827;
  Real copt1334 = copt382 * copt409;
  Real copt1335 = copt421 * copt509;
  Real copt1336 = copt1334 + copt1335;
  Real copt1337 = -(copt1336 * copt591 * copt612 * copt827);
  Real copt1338 = copt1333 + copt1337;
  Real copt1340 = copt340 * copt485;
  Real copt1341 = copt351 * copt409;
  Real copt1342 = copt1340 + copt1341;
  Real copt1343 = copt1342 * copt526 * copt552 * copt844;
  Real copt1344 = -(copt427 * copt518 * copt526 * copt844);
  Real copt1345 = copt1343 + copt1344;
  Real copt1347 = copt485 * copt519;
  Real copt1348 = copt351 * copt509;
  Real copt1349 = copt1347 + copt1348;
  Real copt1350 = copt1349 * copt526 * copt552 * copt844;
  Real copt1351 = copt1011 + copt1178 + copt403 + copt404 + copt405;
  Real copt1352 = -(copt1351 * copt518 * copt526 * copt844);
  Real copt1353 = copt1350 + copt1352;
  Real copt1355 = copt409 * copt519;
  Real copt1356 = copt509 * copt521;
  Real copt1357 = copt1355 + copt1356;
  Real copt1358 = copt1357 * copt526 * copt552 * copt844;
  Real copt1359 = copt1065 + copt1238 + copt380 + copt383 + copt387;
  Real copt1360 = -(copt1359 * copt518 * copt526 * copt844);
  Real copt1361 = copt1358 + copt1360;
  Real copt660  = -(copt476 * copt478 * copt479 * copt659);
  Real copt662  = -(copt554 * copt556 * copt557 * copt661);
  Real copt664  = -(copt614 * copt617 * copt618 * copt663);
  Real copt665  = copt660 + copt662 + copt664;
  out1(0)       = copt356;
  out1(1)       = copt357 * copt372 * copt378;
  out1(2)       = copt377;
  out1(3)       = copt620 * copt625;
  out1(4)       = -(copt479 * copt557 * copt618 * copt625 * copt657);
  out1(5)       = copt625 * copt665;
  out2(0, 0)    = copt357 * copt615 * copt669;
  out2(0, 1)    = copt357 * copt615 * copt673;
  out2(0, 2)    = copt357 * copt615 * copt677;
  out2(0, 3)    = copt325 * copt357 * copt41;
  out2(0, 4)    = copt343 * copt357 * copt41;
  out2(0, 5)    = copt353 * copt357 * copt41;
  out2(0, 6)    = copt236 * copt325 * copt357;
  out2(0, 7)    = copt236 * copt343 * copt357;
  out2(0, 8)    = copt236 * copt353 * copt357;
  out2(0, 9)    = 0;
  out2(0, 10)   = 0;
  out2(0, 11)   = 0;
  out2(0, 12)   = 0;
  out2(0, 13)   = 0;
  out2(0, 14)   = 0;
  out2(0, 15)   = 0;
  out2(0, 16)   = 0;
  out2(0, 17)   = 0;
  out2(1, 0)    = (copt325 * copt372 * copt376 * copt615 +
                copt355 * copt362 * copt372 * copt627 +
                copt355 * copt376 *
                    ((copt358 * copt438 + copt360 * copt519) * copt615 +
                     copt627 * copt669)) *
               copt686 * copt688;
  out2(1, 1) = (copt343 * copt372 * copt376 * copt615 +
                copt355 * copt366 * copt372 * copt627 +
                copt355 * copt376 *
                    ((copt358 * copt440 + copt360 * copt521) * copt615 +
                     copt627 * copt673)) *
               copt686 * copt688;
  out2(1, 2) = (copt353 * copt372 * copt376 * copt615 +
                copt355 * copt370 * copt372 * copt627 +
                copt355 * copt376 *
                    ((copt358 * copt442 + copt360 * copt523) * copt615 +
                     copt627 * copt677)) *
               copt686 * copt688;
  out2(1, 3) = (-(copt355 * copt358 * copt362 * copt372) -
                copt325 * copt372 * copt376 * copt41 +
                copt355 * copt376 *
                    (copt236 * copt298 * copt358 +
                     copt41 * (copt361 - 2 * copt358 * copt438))) *
               copt686 * copt688;
  out2(1, 4) = (-(copt355 * copt358 * copt366 * copt372) -
                copt343 * copt372 * copt376 * copt41 +
                copt355 * copt376 *
                    (copt236 * copt340 * copt358 +
                     copt41 * (copt365 - 2 * copt358 * copt440))) *
               copt686 * copt688;
  out2(1, 5) = (-(copt355 * copt358 * copt370 * copt372) -
                copt353 * copt372 * copt376 * copt41 +
                copt355 * copt376 *
                    (copt236 * copt351 * copt358 +
                     copt41 * (copt369 - 2 * copt358 * copt442))) *
               copt686 * copt688;
  out2(1, 6) = (-(copt355 * copt360 * copt362 * copt372) -
                copt236 * copt325 * copt372 * copt376 +
                copt355 * copt376 *
                    ((copt158 + 2 * copt236 * copt298) * copt360 +
                     copt236 * copt358 * copt57)) *
               copt686 * copt688;
  out2(1, 7) = (-(copt355 * copt360 * copt366 * copt372) +
                copt355 * (copt343 * copt360 + copt236 * copt366) * copt376 -
                copt236 * copt343 * copt372 * copt376) *
               copt686 * copt688;
  out2(1, 8) = copt357 * (copt353 * copt360 + copt236 * copt370) * copt378 -
               copt236 * copt353 * copt372 * copt378 * copt686 -
               copt357 * copt360 * copt370 * copt372 * copt688;
  out2(1, 9)  = 0;
  out2(1, 10) = 0;
  out2(1, 11) = 0;
  out2(1, 12) = 0;
  out2(1, 13) = 0;
  out2(1, 14) = 0;
  out2(1, 15) = 0;
  out2(1, 16) = 0;
  out2(1, 17) = 0;
  out2(2, 0)  = copt362 * copt378 * copt779;
  out2(2, 1)  = copt366 * copt378 * copt779;
  out2(2, 2)  = copt370 * copt378 * copt779;
  out2(2, 3)  = copt358 * copt362 * copt378;
  out2(2, 4)  = copt358 * copt366 * copt378;
  out2(2, 5)  = copt358 * copt370 * copt378;
  out2(2, 6)  = copt360 * copt362 * copt378;
  out2(2, 7)  = copt360 * copt366 * copt378;
  out2(2, 8)  = copt360 * copt370 * copt378;
  out2(2, 9)  = 0;
  out2(2, 10) = 0;
  out2(2, 11) = 0;
  out2(2, 12) = 0;
  out2(2, 13) = 0;
  out2(2, 14) = 0;
  out2(2, 15) = 0;
  out2(2, 16) = 0;
  out2(2, 17) = 0;
  out2(3, 0)  = (copt620 * copt791 * copt792) / 2. +
               copt625 * (copt438 * copt476 * copt477 * copt478 * copt795 +
                          copt519 * copt554 * copt555 * copt556 * copt798 -
                          copt477 * copt478 * copt479 * copt818 -
                          copt616 * copt617 * copt618 * copt833 -
                          copt555 * copt556 * copt557 * copt853);
  out2(3, 1) = (copt620 * copt792 * copt860) / 2. +
               copt625 * (copt440 * copt476 * copt477 * copt478 * copt795 +
                          copt521 * copt554 * copt555 * copt556 * copt798 -
                          copt477 * copt478 * copt479 * copt874 -
                          copt616 * copt617 * copt618 * copt881 -
                          copt555 * copt556 * copt557 * copt893);
  out2(3, 2) = (copt620 * copt792 * copt900) / 2. +
               copt625 * (copt442 * copt476 * copt477 * copt478 * copt795 +
                          copt523 * copt554 * copt555 * copt556 * copt798 -
                          copt477 * copt478 * copt479 * copt914 -
                          copt616 * copt617 * copt618 * copt921 -
                          copt555 * copt556 * copt557 * copt933);
  out2(3, 3) = (copt620 * copt792 * copt940) / 2. +
               copt625 * (-(copt438 * copt476 * copt477 * copt478 * copt795) +
                          copt402 * copt614 * copt616 * copt617 * copt944 -
                          copt477 * copt478 * copt479 * copt963 -
                          copt616 * copt617 * copt618 * copt980 -
                          copt555 * copt556 * copt557 * copt992);
  out2(3, 4) = copt625 * (-(copt1021 * copt477 * copt478 * copt479) -
                          copt1046 * copt555 * copt556 * copt557 -
                          copt1036 * copt616 * copt617 * copt618 -
                          copt440 * copt476 * copt477 * copt478 * copt795 +
                          copt385 * copt614 * copt616 * copt617 * copt944) +
               (copt620 * copt792 * copt999) / 2.;
  out2(3, 5) = (copt1053 * copt620 * copt792) / 2. +
               copt625 * (-(copt1075 * copt477 * copt478 * copt479) -
                          copt1105 * copt555 * copt556 * copt557 -
                          copt1092 * copt616 * copt617 * copt618 -
                          copt442 * copt476 * copt477 * copt478 * copt795 +
                          copt424 * copt614 * copt616 * copt617 * copt944);
  out2(3, 6) = (copt1112 * copt620 * copt792) / 2. +
               copt625 * (-(copt1126 * copt477 * copt478 * copt479) -
                          copt1167 * copt555 * copt556 * copt557 -
                          copt1146 * copt616 * copt617 * copt618 -
                          copt519 * copt554 * copt555 * copt556 * copt798 -
                          copt402 * copt614 * copt616 * copt617 * copt944);
  out2(3, 7) = (copt1174 * copt620 * copt792) / 2. +
               copt625 * (-(copt1189 * copt477 * copt478 * copt479) -
                          copt1227 * copt555 * copt556 * copt557 -
                          copt1207 * copt616 * copt617 * copt618 -
                          copt521 * copt554 * copt555 * copt556 * copt798 -
                          copt385 * copt614 * copt616 * copt617 * copt944);
  out2(3, 8) = (copt1234 * copt620 * copt792) / 2. +
               copt625 * (-(copt1249 * copt477 * copt478 * copt479) -
                          copt1288 * copt555 * copt556 * copt557 -
                          copt1268 * copt616 * copt617 * copt618 -
                          copt523 * copt554 * copt555 * copt556 * copt798 -
                          copt424 * copt614 * copt616 * copt617 * copt944);
  out2(3, 9)  = -(copt1298 * copt477 * copt478 * copt479 * copt625);
  out2(3, 10) = -(copt1309 * copt477 * copt478 * copt479 * copt625);
  out2(3, 11) = -(copt1317 * copt477 * copt478 * copt479 * copt625);
  out2(3, 12) = -(copt1324 * copt616 * copt617 * copt618 * copt625);
  out2(3, 13) = -(copt1331 * copt616 * copt617 * copt618 * copt625);
  out2(3, 14) = -(copt1338 * copt616 * copt617 * copt618 * copt625);
  out2(3, 15) = -(copt1345 * copt555 * copt556 * copt557 * copt625);
  out2(3, 16) = -(copt1353 * copt555 * copt556 * copt557 * copt625);
  out2(3, 17) = -(copt1361 * copt555 * copt556 * copt557 * copt625);
  out2(4, 0) =
      -(copt479 * copt557 * copt618 * copt657 * copt791 * copt792) / 2. +
      copt438 * copt557 * copt618 * copt625 * copt657 * copt795 +
      copt479 * copt519 * copt618 * copt625 * copt657 * copt798 -
      copt479 * copt557 * copt618 * copt625 *
          ((copt1368 * copt1369 * copt614 * copt615 * copt617 * copt627 *
            copt638) /
               2. +
           (copt1372 * copt1373 * copt614 * copt615 * copt617 * copt627 *
            copt646) /
               2. +
           copt615 * copt617 * copt627 * copt638 * copt646 * copt833 +
           copt652 *
               ((copt1368 * copt1369 * copt236 * copt360 * copt476 * copt478) /
                    2. +
                (copt1372 * copt1373 * copt358 * copt41 * copt554 * copt556) /
                    2. +
                copt236 * copt360 * copt478 * copt646 * copt818 +
                copt358 * copt41 * copt556 * copt638 * copt853));
  out2(4, 1) =
      copt440 * copt557 * copt618 * copt625 * copt657 * copt795 +
      copt479 * copt521 * copt618 * copt625 * copt657 * copt798 -
      (copt479 * copt557 * copt618 * copt657 * copt792 * copt860) / 2. -
      copt479 * copt557 * copt618 * copt625 *
          ((copt1369 * copt1390 * copt614 * copt615 * copt617 * copt627 *
            copt638) /
               2. +
           (copt1373 * copt1393 * copt614 * copt615 * copt617 * copt627 *
            copt646) /
               2. +
           copt615 * copt617 * copt627 * copt638 * copt646 * copt881 +
           copt652 *
               ((copt1369 * copt1390 * copt236 * copt360 * copt476 * copt478) /
                    2. +
                (copt1373 * copt1393 * copt358 * copt41 * copt554 * copt556) /
                    2. +
                copt236 * copt360 * copt478 * copt646 * copt874 +
                copt358 * copt41 * copt556 * copt638 * copt893));
  out2(4, 2) =
      copt442 * copt557 * copt618 * copt625 * copt657 * copt795 +
      copt479 * copt523 * copt618 * copt625 * copt657 * copt798 -
      (copt479 * copt557 * copt618 * copt657 * copt792 * copt900) / 2. -
      copt479 * copt557 * copt618 * copt625 *
          ((copt1369 * copt1410 * copt614 * copt615 * copt617 * copt627 *
            copt638) /
               2. +
           (copt1373 * copt1413 * copt614 * copt615 * copt617 * copt627 *
            copt646) /
               2. +
           copt615 * copt617 * copt627 * copt638 * copt646 * copt921 +
           copt652 *
               ((copt1369 * copt1410 * copt236 * copt360 * copt476 * copt478) /
                    2. +
                (copt1373 * copt1413 * copt358 * copt41 * copt554 * copt556) /
                    2. +
                copt236 * copt360 * copt478 * copt646 * copt914 +
                copt358 * copt41 * copt556 * copt638 * copt933));
  out2(4, 3) =
      -(copt438 * copt557 * copt618 * copt625 * copt657 * copt795) -
      (copt479 * copt557 * copt618 * copt657 * copt792 * copt940) / 2. +
      copt402 * copt479 * copt557 * copt625 * copt657 * copt944 -
      copt479 * copt557 * copt618 * copt625 *
          ((copt1373 * copt1430 * copt614 * copt615 * copt617 * copt627 *
            copt646) /
               2. +
           ((copt1367 + copt1429) * copt1433 * copt655) / 2. +
           copt615 * copt617 * copt627 * copt638 * copt646 * copt980 +
           copt652 *
               ((copt1373 * copt1430 * copt358 * copt41 * copt554 * copt556) /
                    2. +
                copt236 * copt360 * copt478 * copt646 * copt963 +
                copt358 * copt41 * copt556 * copt638 * copt992));
  out2(4, 4) =
      -(copt479 * copt557 * copt618 * copt625 *
        ((copt1373 * copt1449 * copt614 * copt615 * copt617 * copt627 *
          copt646) /
             2. +
         copt1036 * copt615 * copt617 * copt627 * copt638 * copt646 +
         ((copt1373 * copt1449 * copt358 * copt41 * copt554 * copt556) / 2. +
          copt1046 * copt358 * copt41 * copt556 * copt638 +
          copt1021 * copt236 * copt360 * copt478 * copt646) *
             copt652 +
         (copt1433 * (copt1389 + copt1448) * copt655) / 2.)) -
      copt440 * copt557 * copt618 * copt625 * copt657 * copt795 +
      copt385 * copt479 * copt557 * copt625 * copt657 * copt944 -
      (copt479 * copt557 * copt618 * copt657 * copt792 * copt999) / 2.;
  out2(4, 5) =
      -(copt479 * copt557 * copt618 * copt625 *
        ((copt1373 * copt1467 * copt614 * copt615 * copt617 * copt627 *
          copt646) /
             2. +
         copt1092 * copt615 * copt617 * copt627 * copt638 * copt646 +
         ((copt1373 * copt1467 * copt358 * copt41 * copt554 * copt556) / 2. +
          copt1105 * copt358 * copt41 * copt556 * copt638 +
          copt1075 * copt236 * copt360 * copt478 * copt646) *
             copt652 +
         (copt1433 * (copt1409 + copt1466) * copt655) / 2.)) -
      (copt1053 * copt479 * copt557 * copt618 * copt657 * copt792) / 2. -
      copt442 * copt557 * copt618 * copt625 * copt657 * copt795 +
      copt424 * copt479 * copt557 * copt625 * copt657 * copt944;
  out2(4, 6) =
      -(copt479 * copt557 * copt618 * copt625 *
        ((copt1369 * copt1484 * copt614 * copt615 * copt617 * copt627 *
          copt638) /
             2. +
         copt1146 * copt615 * copt617 * copt627 * copt638 * copt646 +
         ((copt1369 * copt1484 * copt236 * copt360 * copt476 * copt478) / 2. +
          copt1167 * copt358 * copt41 * copt556 * copt638 +
          copt1126 * copt236 * copt360 * copt478 * copt646) *
             copt652 +
         (copt1433 * (copt1371 + copt1483) * copt655) / 2.)) -
      (copt1112 * copt479 * copt557 * copt618 * copt657 * copt792) / 2. -
      copt479 * copt519 * copt618 * copt625 * copt657 * copt798 -
      copt402 * copt479 * copt557 * copt625 * copt657 * copt944;
  out2(4, 7) =
      -(copt479 * copt557 * copt618 * copt625 *
        ((copt1369 * copt1501 * copt614 * copt615 * copt617 * copt627 *
          copt638) /
             2. +
         copt1207 * copt615 * copt617 * copt627 * copt638 * copt646 +
         ((copt1369 * copt1501 * copt236 * copt360 * copt476 * copt478) / 2. +
          copt1227 * copt358 * copt41 * copt556 * copt638 +
          copt1189 * copt236 * copt360 * copt478 * copt646) *
             copt652 +
         (copt1433 * (copt1392 + copt1500) * copt655) / 2.)) -
      (copt1174 * copt479 * copt557 * copt618 * copt657 * copt792) / 2. -
      copt479 * copt521 * copt618 * copt625 * copt657 * copt798 -
      copt385 * copt479 * copt557 * copt625 * copt657 * copt944;
  out2(4, 8) =
      -(copt479 * copt557 * copt618 * copt625 *
        ((copt1369 * copt1518 * copt614 * copt615 * copt617 * copt627 *
          copt638) /
             2. +
         copt1268 * copt615 * copt617 * copt627 * copt638 * copt646 +
         ((copt1369 * copt1518 * copt236 * copt360 * copt476 * copt478) / 2. +
          copt1288 * copt358 * copt41 * copt556 * copt638 +
          copt1249 * copt236 * copt360 * copt478 * copt646) *
             copt652 +
         (copt1433 * (copt1412 + copt1517) * copt655) / 2.)) -
      (copt1234 * copt479 * copt557 * copt618 * copt657 * copt792) / 2. -
      copt479 * copt523 * copt618 * copt625 * copt657 * copt798 -
      copt424 * copt479 * copt557 * copt625 * copt657 * copt944;
  out2(4, 9)  = -(copt1298 * copt236 * copt360 * copt478 * copt479 * copt557 *
                 copt618 * copt625 * copt646 * copt652);
  out2(4, 10) = -(copt1309 * copt236 * copt360 * copt478 * copt479 * copt557 *
                  copt618 * copt625 * copt646 * copt652);
  out2(4, 11) = -(copt1317 * copt236 * copt360 * copt478 * copt479 * copt557 *
                  copt618 * copt625 * copt646 * copt652);
  out2(4, 12) = -(copt1324 * copt479 * copt557 * copt615 * copt617 * copt618 *
                  copt625 * copt627 * copt638 * copt646);
  out2(4, 13) = -(copt1331 * copt479 * copt557 * copt615 * copt617 * copt618 *
                  copt625 * copt627 * copt638 * copt646);
  out2(4, 14) = -(copt1338 * copt479 * copt557 * copt615 * copt617 * copt618 *
                  copt625 * copt627 * copt638 * copt646);
  out2(4, 15) = -(copt1345 * copt358 * copt41 * copt479 * copt556 * copt557 *
                  copt618 * copt625 * copt638 * copt652);
  out2(4, 16) = -(copt1353 * copt358 * copt41 * copt479 * copt556 * copt557 *
                  copt618 * copt625 * copt638 * copt652);
  out2(4, 17) = -(copt1361 * copt358 * copt41 * copt479 * copt556 * copt557 *
                  copt618 * copt625 * copt638 * copt652);
  out2(5, 0)  = (copt665 * copt791 * copt792) / 2. +
               copt625 * (copt438 * copt476 * copt478 * copt659 * copt795 +
                          copt519 * copt554 * copt556 * copt661 * copt798 -
                          copt478 * copt479 * copt659 * copt818 -
                          copt617 * copt618 * copt663 * copt833 -
                          copt556 * copt557 * copt661 * copt853);
  out2(5, 1) = (copt665 * copt792 * copt860) / 2. +
               copt625 * (copt440 * copt476 * copt478 * copt659 * copt795 +
                          copt521 * copt554 * copt556 * copt661 * copt798 -
                          copt478 * copt479 * copt659 * copt874 -
                          copt617 * copt618 * copt663 * copt881 -
                          copt556 * copt557 * copt661 * copt893);
  out2(5, 2) = (copt665 * copt792 * copt900) / 2. +
               copt625 * (copt442 * copt476 * copt478 * copt659 * copt795 +
                          copt523 * copt554 * copt556 * copt661 * copt798 -
                          copt478 * copt479 * copt659 * copt914 -
                          copt617 * copt618 * copt663 * copt921 -
                          copt556 * copt557 * copt661 * copt933);
  out2(5, 3) = (copt665 * copt792 * copt940) / 2. +
               copt625 * (-(copt438 * copt476 * copt478 * copt659 * copt795) +
                          copt402 * copt614 * copt617 * copt663 * copt944 -
                          copt478 * copt479 * copt659 * copt963 -
                          copt617 * copt618 * copt663 * copt980 -
                          copt556 * copt557 * copt661 * copt992);
  out2(5, 4) = copt625 * (-(copt1021 * copt478 * copt479 * copt659) -
                          copt1046 * copt556 * copt557 * copt661 -
                          copt1036 * copt617 * copt618 * copt663 -
                          copt440 * copt476 * copt478 * copt659 * copt795 +
                          copt385 * copt614 * copt617 * copt663 * copt944) +
               (copt665 * copt792 * copt999) / 2.;
  out2(5, 5) = (copt1053 * copt665 * copt792) / 2. +
               copt625 * (-(copt1075 * copt478 * copt479 * copt659) -
                          copt1105 * copt556 * copt557 * copt661 -
                          copt1092 * copt617 * copt618 * copt663 -
                          copt442 * copt476 * copt478 * copt659 * copt795 +
                          copt424 * copt614 * copt617 * copt663 * copt944);
  out2(5, 6) = (copt1112 * copt665 * copt792) / 2. +
               copt625 * (-(copt1126 * copt478 * copt479 * copt659) -
                          copt1167 * copt556 * copt557 * copt661 -
                          copt1146 * copt617 * copt618 * copt663 -
                          copt519 * copt554 * copt556 * copt661 * copt798 -
                          copt402 * copt614 * copt617 * copt663 * copt944);
  out2(5, 7) = (copt1174 * copt665 * copt792) / 2. +
               copt625 * (-(copt1189 * copt478 * copt479 * copt659) -
                          copt1227 * copt556 * copt557 * copt661 -
                          copt1207 * copt617 * copt618 * copt663 -
                          copt521 * copt554 * copt556 * copt661 * copt798 -
                          copt385 * copt614 * copt617 * copt663 * copt944);
  out2(5, 8) = (copt1234 * copt665 * copt792) / 2. +
               copt625 * (-(copt1249 * copt478 * copt479 * copt659) -
                          copt1288 * copt556 * copt557 * copt661 -
                          copt1268 * copt617 * copt618 * copt663 -
                          copt523 * copt554 * copt556 * copt661 * copt798 -
                          copt424 * copt614 * copt617 * copt663 * copt944);
  out2(5, 9)  = -(copt1298 * copt478 * copt479 * copt625 * copt659);
  out2(5, 10) = -(copt1309 * copt478 * copt479 * copt625 * copt659);
  out2(5, 11) = -(copt1317 * copt478 * copt479 * copt625 * copt659);
  out2(5, 12) = -(copt1324 * copt617 * copt618 * copt625 * copt663);
  out2(5, 13) = -(copt1331 * copt617 * copt618 * copt625 * copt663);
  out2(5, 14) = -(copt1338 * copt617 * copt618 * copt625 * copt663);
  out2(5, 15) = -(copt1345 * copt556 * copt557 * copt625 * copt661);
  out2(5, 16) = -(copt1353 * copt556 * copt557 * copt625 * copt661);
  out2(5, 17) = -(copt1361 * copt556 * copt557 * copt625 * copt661);
  return std::make_tuple(grad, val);
}

#endif  // hylc_strain_II
