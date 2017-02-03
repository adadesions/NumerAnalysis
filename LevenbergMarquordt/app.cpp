#include <iostream>
#include <stdio.h>
#include <math.h>
#include <vector>
#include "./lib/adavector.hpp"

using namespace std;

#define VNAME(x) cout << #x << " = " << endl;

typedef vector<vector<double>> adavector;

double f(adavector x){
  return pow(x[0][0]-9, 2) + pow(x[1][0]-6, 2) + 7;
}

vector<double> partialDiff(adavector x, int pos){
  const double alpha = 0.01;
  adavector temp(x);
  double fplus, fminus;
  // f(x+h/2,y) or f(x,y+h/2)
  temp[pos][0] = x[pos][0]+(double)alpha/2;
  fplus = f(temp);

  // f(x-h/2,y) or f(x,y-h/2)
  temp[pos][0] = x[pos][0]-(double)alpha/2;
  fminus = f(temp);
  return {(double)(fplus-fminus)/alpha};
}

adavector grad(adavector x){
  adavector result(x);
  result[0][0] = partialDiff(x, 0)[0];
  result[1][0] = partialDiff(x, 1)[0];
  return result;
}


adavector bfgsMethod(adavector b, adavector y, adavector dx, double alpha){
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

  //DEBUG ZONE
  // dispVector(yyt, "Debug yyt");
  // dispVector(bsbtst, "Debug bsbtst");
  // dispVector(result, "Before add B0");
  // dispVector(addVector(result,b), "added B0");
  //END DEBUG ZONE
  return addVector(result,b);
}

adavector cal_B(adavector b, adavector y, adavector dx, double alpha){
  return bfgsMethod(b, y, dx, alpha);
}

adavector LM(adavector b,adavector x1, double lambda){
  // Init
  double alpha = 1;
  adavector gradX1 = grad(x1);

  //Computation
  adavector dx = multiplyVector(inverseVector(b), grad(x1));
  adavector x2 = subtractVector(x1, scalarMulVector(alpha, dx));
  adavector y = subtractVector(grad(x2), grad(x1));
  adavector B(cal_B(b,y,dx,alpha)); // Bk+1
  adavector DB(diagonal(B)); // Diagonal of B
  DB = scalarMulVector(lambda, DB); // Multiply lambda with Diagonal B
  adavector IDB = inverseVector(addVector(B, DB)); // Inverse DB
  IDB = multiplyVector(IDB, grad(x1)); // Multiply with grad(f)
  adavector x_new = adavector(addVector(x1,IDB)); // Got new x

  //Halt case
  if( fabs(f(x_new)-f(x1)) < 0.01)
    return x_new;

  dispVector(x_new, "x* procrssing");

  if( f(x_new) >= f(x1) ){
    return LM(B, x_new, lambda*10);
  }
  return LM(B, x_new, lambda/10);
}

adavector QN(adavector b, adavector x1, double alpha){
  adavector gf1 = grad(x1);
  adavector dx = multiplyVector(inverseVector(b), gf1);
  adavector x_new = subtractVector(x1, scalarMulVector(alpha, dx));
  adavector gf2 = grad(x_new);
  adavector y = subtractVector(gf2, gf1);
  adavector b_new = cal_B(b, y, dx, alpha);
  printf("f_new, f_old : %.2f %.2f\n", f(x_new), f(x1));
  printf("diff = %.2f\n", f(x_new)-f(x1));
  printf("Dx = %f\n", norm(dx));
  if( norm(dx) <= 0.001 ){
    return x_new;
  }
  dispVector(x_new, "x* procrssing");
  if(f(x_new) > f(x1))
    return QN(b_new, x_new, alpha-0.5);
  return QN(b_new, x_new, alpha+0.5);
}

int main(){
  double lambda = 10;
  adavector x0 {{2}, {2}};
  adavector x00 {{5}, {8}};
  adavector b0 {{x0[0][0], 0}, {0, x0[1][0]}};
  // dispVector(QN(b0,x0,1), "QN - x*");
  dispVector(LM(b0, x00, lambda), "LM - x*");
  // dispVector(LM(b0, x0, lambda, 0.000001), "x*");
  return 0;
}
