// Copyright (C) 2014, HydroComplexity Group
// All rights reserved.
//
// GPU-BASED CONJUNCTIVE SURFACE-SUBSURFACE FLOW MODEL (GCSFlow)
// GCSFlow model is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published
// by the Free Software Foundation; either version 2.1 of the License, or
// (at your option) any later version.
//
// GCSFlow is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
// or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
// for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with GCSFlow; if not, see <http://www.gnu.org/licenses/>.
//
// Author: levuvietphong@gmail.com (Phong Le)


// --------------------------------------------------------------------
// ThomasAlgorithm()
//    Solves tridiagonal matrix (size N*N) system Af=d using Thomas
//    Algorithm - TDMA (or Gauss Elimination).
//          A - a tridiagonal matrix including:
//              a: sub-diagonal
//              b: main diagonal
//              c: sup-diagonal
//          d - vector of RHS
//          f - vector solution
//          c_star - modified values of c
//          d_star - modified values of d
//    Input: a(N-1), b(N), c(N-1), d(N)
//    Output: f(N)
// --------------------------------------------------------------------
__device__ void ThomasAlgorithm(double *a, double *b, double *c, double *d,
    double *f, double *c_star, double *d_star, int N) {
  double m;

  // Modify the first-row coefficients
  c_star[0] = c[0] / b[0];
  d_star[0] = d[0] / b[0];

  // Forward sweep and modification
  for (int i = 1; i < N; i++) {
    m = 1.0 / (b[i] - a[i] * c_star[i-1]);
    c_star[i] = c[i] * m;
    d_star[i] = (d[i] - a[i] * d_star[i-1]) * m;
  }

  // Backward substitution
  f[N-1] = d_star[N-1];
  for (int i = N-2; i >= 0; i--) {
    f[i] = d_star[i] - c_star[i] * f[i+1];
  }
}

// -----------------------------------------------------------------------------
// :::::::::::::::::::::::::::::::: END OF FILE ::::::::::::::::::::::::::::::::
// -----------------------------------------------------------------------------