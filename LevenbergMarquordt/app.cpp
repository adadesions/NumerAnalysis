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

    int rows, cols, curCol;
    rows = result.size();
    cols = result[0].size();

    for(int i=0; i<rows; i++){
      for(int j=0; j<cols; j++){
        double temp = 0;
        for(int k=0; k<rows; k++){
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

int main(){
  adavector x0 {{1}, {2}};
  adavector b0 {{x0[0][0], 0}, {0, x0[1][0]}};
  adavector ib0(inverseVector(b0));
  adavector gf(grad(x0));
  adavector x1;
  double alpha = 0.1;

  dispVector(ib0,"ib0");
  dispVector(gf, "gradient");
  dispVector(multiplyVector(ib0, gf), "delta x0");
  dispVector(subtractVector(b0, ib0), "Subtract 2b0");

  return 0;
}
