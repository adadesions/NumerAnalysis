#pragma once
#include <iostream>
#include <stdio.h>
#include <math.h>
#include <vector>

using namespace std;
typedef vector<vector<double>> adavector;

namespace adalib {
  double norm(adavector v){
    return sqrt(pow(v[0][0],2) + pow(v[1][0],2));
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
      printf("Cannot operate because dimantaion are not fit\n");
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

  adavector diagonal(adavector v){
    for(int i=0; i<v.size(); i++){
      for(int j=0; j<v[0].size(); j++){
        if( i != j)
          v[i][j] = 0;
      }
    }
    return v;
  }

}
