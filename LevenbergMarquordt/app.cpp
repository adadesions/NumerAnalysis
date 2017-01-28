#include <iostream>
#include <stdio.h>
#include <math.h>
#include <vector>

using namespace std;

#define VNAME(x) cout << #x << " = " << endl;

typedef vector<vector<double>> adavector;

double f(adavector x){
  return pow(x[0][0]-9, 2) + pow(x[1][0]-6, 2) + 7;
}

vector<double> partialDiff(adavector x, int pos){
  const double alpha = 0.0001;
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

double det2(adavector v){
  return v[0][0]*v[1][1] - v[0][1]*v[1][0];
}

adavector inverseVector(adavector v){
  double det = det2(v);
  double temp = v[0][0];
  v[0][0] = v[1][1];
  v[1][1] = temp;
  v[0][1] = (-1)*v[0][1];
  v[1][0] = (-1)*v[1][0];
  for(int i=0; i<2; i++){
    for(int j=0; j<2; j++){
      v[i][j] = v[i][j]/det;
    }
  }
  return v;
}

adavector multiplyVector(adavector a, adavector b){
  adavector result;
  int dim1[2], dim2[2];
  dim1[0] = a.size();
  dim1[1] = a[0].size();
  dim2[0] = b.size();
  dim2[1] = b[0].size();

  if(dim1[1] == dim2[0]){
    result.resize(dim1[0]);
    for(int i=0; i<dim1[0]; i++)
      result[i].resize(dim2[1]);

    int rows, cols, inter;
    rows = result.size();
    cols = result[0].size();
    inter = dim1[1];

    for(int i=0; i<rows; i++){
      for(int j=0; j<cols; j++){
        double temp = 0;
        for(int k=0; k<inter; k++){
          temp += a[i][k]*b[k][j];
        }
        result[i][j] = temp;
      }
    }
    return result;
  }
  else{
    printf("Cannot multypy these two matrices because dimantion are not fit\n");
    return result;
  }
}

adavector operationVector(bool operation, adavector a, adavector b){
  adavector result(a);
  int dim1[2], dim2[2];
  dim1[0] = a.size();
  dim1[1] = a[0].size();
  dim2[0] = b.size();
  dim2[1] = b[0].size();
  bool isColsEqual = dim1[0] == dim2[0];
  bool isRowsEqual = dim1[1] == dim2[1];

  if(isRowsEqual && isColsEqual){
    for(int i=0; i<dim1[0]; i++){
      for(int j=0; j<dim1[1]; j++){
        if(operation)
          result[i][j] = a[i][j]+b[i][j];
        else
          result[i][j] = a[i][j]-b[i][j];
      }
    }
  }
  else{
    printf("Cannot add because dimantaion are not fit");
  }
  return result;
}

adavector addVector(adavector a, adavector b){
  return operationVector(true, a, b);
}

adavector subtractVector(adavector a, adavector b){
  return operationVector(false, a, b);
}

adavector scalarVector(bool operation, double s, adavector v){
  int rows, cols;
  rows = v.size();
  cols = v[0].size();
  for(int i=0; i<rows; i++){
    for(int j=0; j<cols; j++){
      if(operation)
        v[i][j] = s*v[i][j];
      else
        v[i][j] = v[i][j]/s;
    }
  }
  return v;
}

adavector scalarMulVector(double s, adavector v){
  return scalarVector(true, s, v);
}

adavector scalarDivVector(double s, adavector v){
  return scalarVector(false, s, v);
}

void dispVector(adavector v, string dispName){
  int rows, cols;
  rows = v.size();
  cols = v[0].size();
  cout << dispName << " = " << endl;
  for(int i=0; i<rows; i++){
    for(int j=0; j<cols; j++){
      printf("%.2f ", v[i][j]);
    }
    printf("\n");
  }
  printf("\n");
}

adavector transpose(adavector v){
  adavector t;
  int rows,cols;
  rows = v.size();
  cols = v[0].size();
  t.resize(cols);
  for(int i=0; i<rows; i++){
    t[i].resize(rows);
  }
  rows = t.size();
  cols = t[0].size();
  for (int i = 0; i < rows; i++) {
    for (int j = 0; j < cols; j++) {
      t[i][j] = v[j][i];
    }
  }
  return t;
}

adavector cal_B(adavector b, adavector y, adavector dx, double alpha){
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

  dispVector(yyt, "yyt");
  dispVector(bsbtst, "bsbtst");

  yts = multiplyVector(yt , s)[0][0];
  stbs = multiplyVector(stb, s)[0][0];
  yyt = scalarDivVector(yts, yyt);
  bsbtst = scalarDivVector(stbs, bsbtst);


  printf("yts = %.2f\n", yts);
  dispVector(yyt, "yyt/yts");

  printf("stbs = %.2f\n", stbs);
  dispVector(bsbtst, "bsbtst/stbs");

  result = subtractVector(yyt, bsbtst);
  b = addVector(b, result);
  dispVector(b, "New B");

  return b;
}

int main(){
  double alpha = 0.1;
  adavector x0 {{1}, {2}};
  adavector b0 {{x0[0][0], 0}, {0, x0[1][0]}};
  adavector ib0(inverseVector(b0));
  adavector gf(grad(x0));
  adavector dx0(multiplyVector(ib0, gf));
  adavector x1, y0;

  x1 = subtractVector(x0,scalarMulVector(alpha, dx0));
  y0 = subtractVector(grad(x1), gf);
  dispVector(x0, "x0");
  dispVector(scalarMulVector(alpha, dx0), "alpha*dx0");
  dispVector(x1, "x1");
  dispVector(y0, "y0");
  cal_B(b0, y0, dx0, alpha);

  return 0;
}
