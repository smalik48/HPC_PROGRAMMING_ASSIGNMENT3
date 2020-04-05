/**
 * @file    jacobi.cpp
 * @author  Patrick Flick <patrick.flick@gmail.com>
 * @brief   Implements matrix vector multiplication and Jacobi's method.
 *
 * Copyright (c) 2014 Georgia Institute of Technology. All Rights Reserved.
 */
#include "jacobi.h"

/*
 * TODO: Implement your solutions here
 */

// my implementation:
#include <iostream>
#include <math.h>
#include <vector>
#include <algorithm>
#include <functional>

using namespace std;

void find_DandR(const int n, const double* A, vector<double>& D, vector<double>& R){
  //go over matrix A
  for(int i = 0; i < (n*n) ; i++){
    //make everything zero except diagonal Elements
    //diagonal elements
    if(i % (n + 1) == 0) D[i] = A[i];
    //non diagonal elements
    else R[i] = A[i];
  }
}

void invert_D(const vector<double>& D, vector<double>& D_inv){
  for(int i = 0; i < D.size(); i++){
    if(D[i]!=0) D_inv[i] = 1/D[i];
  }
}

// Calculates y = A*x for a square n-by-n matrix A, and n-dimensional vectors x
// and y
void matrix_vector_mult(const int n, const double* A, const double* x, double* y)
{
    //TODO
    for (int j = 0; j < n; j++){
      y[j] = 0.0;
      for(int i = 0; i< n; i++)
        y[j] += A[j*n+i] * x[i];
    }
}

// Calculates y = A*x for a n-by-m matrix A, a m-dimensional vector x
// and a n-dimensional vector y
void matrix_vector_mult(const int n, const int m, const double* A, const double* x, double* y)
{
    // TODO
    for (int j = 0; j < n; j++){
      y[j] = 0.0;
      for(int i = 0; i< m; i++)
        y[j] += A[j*m+i] * x[i];
    }
}

double l2_norm(vector<double> const& u) {
    double accum = 0.;
    for (int i = 0; i < u.size(); i++) {
        accum += u[i] * u[i];
    }
    return sqrt(accum);
}

double L2_norm(const int n, const double* A, const double* x, const double* b){
  //y_temp = A*x
  vector<double> y_temp(n);
  matrix_vector_mult(n, A, x, &y_temp[0]);

  //y_temp = y_temp - b
  transform(y_temp.begin(), y_temp.end(), b, y_temp.begin(), minus<double>());
  //norm of A*x - b
  return l2_norm(y_temp);
}

void update_x(const int n, const vector<double> &D_inv, const double* b, const vector<double>& R, double* x){
  //y_temp = R*x
  vector<double> y_temp(n);
  matrix_vector_mult(n, &R[0], x, &y_temp[0]);

  //y_temp = b-R*x
  transform(b, b + n, y_temp.begin(), y_temp.begin(), minus<double>());

  for(int i = 0; i < n; i++){
    x[i] = D_inv[i*(n+1)] * y_temp[i];
  }
}

// implements the sequential jacobi method
void jacobi(const int n, double* A, double* b, double* x, int max_iter, double l2_termination)
{
    // TODO
    //initialize x to zero
    for(int i = 0 ;i < n; i++){
      x[i] = 0;
    }

    //D = diag(A)
    vector<double> D(n*n, 0);
    //R = A - D
    vector<double> R(n*n, 0);
    find_DandR(n, A, D, R);

    vector<double> D_inv(D.size(), 0);
    invert_D(D, D_inv);

    int iter = 0;
    while (L2_norm(n, A, x, b) > l2_termination && iter < max_iter){
      update_x(n, D_inv, b, R, x);
      iter++;
    }
}
