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

adavector DFP(adavector b, adavector y, adavector dx){
  // Init
  adavector I {{1,0},{0,1}};
  adavector yt(transpose(y));
  adavector dxt(transpose(dx));
  adavector first, second, third;
  adavector result(I);
  double ytdx(multiplyVector(yt, dx)[0][0]);

  //Computation
  //First Term
  first = scalarDivVector(ytdx, multiplyVector(y, dxt));
  first = subtractVector(I, first);
  //Second Term
  second = scalarDivVector(ytdx, multiplyVector(dx, yt));
  second  = subtractVector(I, second);
  //Third Term
  third = scalarDivVector(ytdx, multiplyVector(y, yt));

  //Merge for Result
  result = multiplyVector(multiplyVector(first, b), second);
  result = addVector(result, third);
  return result;
}
