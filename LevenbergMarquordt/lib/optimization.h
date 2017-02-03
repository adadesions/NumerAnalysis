#pragma once
#include <iostream>
#include <stdio.h>
#include <math.h>
#include <vector>
#include "./adalib.h"

using namespace std;
using namespace adalib;

adavector BFGS(adavector b, adavector y, adavector dx, double alpha){
  double yts, stbs;
  adavector yt(transpose(y));
  adavector s(scalarMulVector(alpha, dx));
  adavector st(transpose(s));
  adavector bs(multiplyVector(b, s));
  adavector btst(transpose(bs));
  adavector stb(multiplyVector(st, b));
  adavector result(b);
  adavector yyt(multiplyVector(y, yt));
  adavector bsbtst(multiplyVector(bs, btst));

  yts = multiplyVector(yt , s)[0][0];
  stbs = multiplyVector(stb, s)[0][0];
  yyt = scalarDivVector(yts, yyt);
  bsbtst = scalarDivVector(stbs, bsbtst);
  result = subtractVector(yyt, bsbtst);
  return addVector(result,b);
}

// adavector DFP(adavector b, adavector y, adavector dx, double alpha){}
