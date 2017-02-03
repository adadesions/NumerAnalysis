#include <iostream>
#include <stdio.h>
#include <math.h>
#include <vector>
#include "./lib/adalib.h"
#include "./lib/optimization.h"

using namespace std;
using namespace adalib;

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

adavector cal_B(adavector b, adavector y, adavector dx, double alpha, int gramma){
  adavector fBFGS(scalarMulVector(1-gramma, BFGS(b, y, dx, alpha)));
  adavector fDFP(scalarMulVector(gramma, DFP(b, y, dx)));
  return addVector(fBFGS, fDFP);
}

adavector LM(adavector b,adavector x1, double lambda){
  // Init
  double alpha = 1;
  adavector gradX1 = grad(x1);

  //Computation
  adavector dx = multiplyVector(inverseVector(b), grad(x1));
  adavector x2 = subtractVector(x1, scalarMulVector(alpha, dx));
  adavector y = subtractVector(grad(x2), grad(x1));
  adavector B(cal_B(b,y,dx,alpha,(int)norm(y)%2)); // Bk+1
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
  adavector b_new = cal_B(b, y, dx, alpha, (int)norm(y)%2);
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
  double lambda = 0.001;
  adavector x0 {{2}, {1}};  
  adavector b0 {{x0[0][0], 0}, {0, x0[1][0]}};
  dispVector(LM(b0, x0, lambda), "LM - x*");
  return 0;
}
