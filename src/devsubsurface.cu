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

// This file provides device functions for subsurface flow model. Device 
// functions are linked and called from global subsurface flow functions 

#include "../include/constants.h"
#include "../include/constantextern.h"


// --------------------------------------------------------------------
// RichZMatrix()
// RichXMatrix()
// RichYMatrix()
//    These functions create the matrix A on LHS and vector RHS d for
//    tridiagonal systems in ADI method for sub-surface flow model.
//    RichZMatrix() for vertical z-direction calculation.
//    RichXMatrix() and RichYMatrix() for horizontal x and y directions.
//
//    Input...
//        * Cnp1m: Specific moisture capacity at time n+1, interation m [1/L]
//        * Knp1m: Hydraulic conductivity at time n+1, interation m [L/T]
//        * hnp1m: Pressure head at time n+1, interation m [L]
//        * thetanp1m: Soil moisture at time n+1, interation m [-]
//        * K_prev: Hydraulic conductivity at time n [L/T]
//        * h_prev: Pressure head at time n [L]
//        * thetan: Soil moisture at time n [L/T]
//        * Scalar indexes
//
//    Output...
//        Matrix A and vector RHS d
//        * a: sub-diagonal of A
//        * b: main diagonal of A
//        * c: sup-diagonal of A
//        * d: RHS vector
// --------------------------------------------------------------------
__device__ void RichZMatrix(double *a, double *b, double *c, double *d,
    double *Cnp1m, double *Knp1m, double *K_prev, double *hnp1m, double *h_prev,
    double *thetanp1m, double *thetan, int i, int j, int M, int N, int P) {
  int I;
  double Knp1m_up, Knp1m_down, Knp1m_east, Knp1m_west, Knp1m_north, Knp1m_south;

  for (int k = 0; k < P; k++) {
    I = k*M*N + j*N + i;    // Get 3D global index

    if (k == 0 || k == P-1) {
      // In z-direction
      // At top or bottom cells, no averaging
      Knp1m_up = Knp1m[k];
      Knp1m_down = Knp1m[k];
    } else {
      // Averaging for other cells
      Knp1m_up = 0.5 * (Knp1m[k] + Knp1m[k+1]);         // K(i,j,k+1/2)
      Knp1m_down = 0.5 * (Knp1m[k] + Knp1m[k-1]);       // K(i,j,k-1/2)
    }

    // In x-direction and y-direction
    Knp1m_east = 0.5 * (Knp1m[k] + K_prev[I+1]);        // K(i+1/2,j,k)
    Knp1m_west = 0.5 * (Knp1m[k] + K_prev[I-1]);        // K(i-1/2,j,k)
    Knp1m_north = 0.5 * (Knp1m[k] + K_prev[I+N]);       // K(i,j+1/2,k)
    Knp1m_south = 0.5 * (Knp1m[k] + K_prev[I-N]);       // K(i,j-1/2,k)

    if (k == 0 || k == P-1) {
      // At top and bottom boundary cells
      b[k] = (Cnp1m[k]+Ss*thetanp1m[k]/poros)/(1.0/3.0*dt);
      a[k] = 0.0;
      c[k] = 0.0;
      d[k] = (3.0/dt)*(Ss*thetanp1m[k]/poros)*(h_prev[I] - hnp1m[k])
            - (3.0/dt)*(thetanp1m[k] - thetan[k]);
    } else {
      // Not at top and bottom boundary cells
      b[k] = (Cnp1m[k]+Ss*thetanp1m[k]/poros)/(1.0/3.0*dt)
            + (Knp1m_up+Knp1m_down)/(dz*dz);
      a[k] = -Knp1m_down/(dz*dz);
      c[k] = -Knp1m_up/(dz*dz);
      d[k] =  (1.0/(dz*dz)) * (Knp1m_up*(hnp1m[k+1]-hnp1m[k])-Knp1m_down*(hnp1m[k]-hnp1m[k-1]))
          + (1.0/(dx*dx)) * (Knp1m_east*(h_prev[I+M*N]-hnp1m[k]) - Knp1m_west*(hnp1m[k]-h_prev[I-M*N]))
          + (1.0/(dy*dy)) * (Knp1m_north*(h_prev[I+N]-hnp1m[k]) - Knp1m_south*(hnp1m[k]-h_prev[I-N]))
          + (3.0/dt)*(Ss*thetanp1m[k]/poros)*(h_prev[I] - hnp1m[k])
          - (1.0/dz)*(Knp1m_up - Knp1m_down) - (3.0/dt)*(thetanp1m[k] - thetan[k]);
    }
  }
}

__device__ void RichXMatrix(double *a, double *b, double *c, double *d,
    double *Cnp1m, double *Knp1m, double *K_prev, double *hnp1m, double *h_prev,
    double *thetanp1m, double *thetan, int j, int k, int M, int N, int P) {
  int I, r;
  double Knp1m_east, Knp1m_west, Knp1m_north, Knp1m_south, Knp1m_up, Knp1m_down;
  double hnp1m_, Knp1m_;

  for (int i = 1; i < N-1; i++) {
    I = k*M*N + j*N + i;      // Get 3D global index
    r = i-1;

    Knp1m_ = Knp1m[i];
    hnp1m_ = hnp1m[i];
    Knp1m_east = 0.5 * (Knp1m_ + Knp1m[i+1]);
    Knp1m_west = 0.5 * (Knp1m_ + Knp1m[i-1]);
    Knp1m_north = 0.5 * (Knp1m_ + K_prev[I+N]);
    Knp1m_south = 0.5 * (Knp1m_ + K_prev[I-N]);
    Knp1m_up = 0.5 * (Knp1m_ + K_prev[I+M*N]);
    Knp1m_down = 0.5 * (Knp1m_ + K_prev[I-M*N]);

    b[r] = (Cnp1m[i]+Ss*thetanp1m[i]/poros)/(1.0/3.0*dt) + (Knp1m_east+Knp1m_west)/(dx*dx);
    a[r] = -Knp1m_west/(dx*dx);
    c[r] = -Knp1m_east/(dx*dx);
    d[r] =  (1.0/(dx*dx)) * (Knp1m_east*(hnp1m[i+1]-hnp1m_)-Knp1m_west*(hnp1m_-hnp1m[i-1]))
        + (1.0/(dy*dy)) * (Knp1m_north*(h_prev[I+N]-hnp1m_)-Knp1m_south*(hnp1m_-h_prev[I-N]))
        + (1.0/(dz*dz)) * (Knp1m_up*(h_prev[I+M*N]-hnp1m_)-Knp1m_down*(hnp1m_-h_prev[I-M*N]))
        + (3.0/dt)*(Ss*thetanp1m[i]/poros)*(h_prev[I] - hnp1m_)
        - (1.0/dz)*(Knp1m_up - Knp1m_down) 
        - (3.0/dt)*(thetanp1m[i] - thetan[i]);
  }
}

__device__ void RichYMatrix(double *a, double *b, double *c, double *d,
    double *Cnp1m, double *Knp1m, double *K_prev, double *hnp1m, double *h_prev,
    double *thetanp1m, double *thetan, int i, int k, int M, int N, int P) {
  int I, r;
  double Knp1m_east, Knp1m_west, Knp1m_north, Knp1m_south, Knp1m_up, Knp1m_down;
  double hnp1m_, Knp1m_;

  for (int j = 1; j < M-1; j++) {
    I = k*M*N + j*N + i;      // Get 3D global index
    r = j-1;
    Knp1m_ = Knp1m[j];
    hnp1m_ = hnp1m[j];
    Knp1m_east = 0.5 * (Knp1m_ + K_prev[I+1]);
    Knp1m_west = 0.5 * (Knp1m_ + K_prev[I-1]);
    Knp1m_north = 0.5 * (Knp1m_ + Knp1m[j+1]);
    Knp1m_south = 0.5 * (Knp1m_ + Knp1m[j-1]);
    Knp1m_up = 0.5 * (Knp1m_ + K_prev[I+M*N]);
    Knp1m_down = 0.5 * (Knp1m_ + K_prev[I-M*N]);

    b[r] = (Cnp1m[j]+Ss*thetanp1m[j]/poros)/(1.0/3.0*dt) + (Knp1m_south+Knp1m_north)/(dy*dy);
    a[r] = -Knp1m_south/(dy*dy);
    c[r] = -Knp1m_north/(dy*dy);
    d[r] =  (1.0/(dy*dy)) * (Knp1m_north*(hnp1m[j+1]-hnp1m_)-Knp1m_south*(hnp1m_-hnp1m[j-1]))
        + (1.0/(dx*dx)) * (Knp1m_east*(h_prev[I+1]-hnp1m_)-Knp1m_west*(hnp1m_-h_prev[I-1]))
        + (1.0/(dz*dz)) *(Knp1m_up*(h_prev[I+M*N]-hnp1m_)-Knp1m_down*(hnp1m_-h_prev[I-M*N]))
        + (3.0/dt)*(Ss*thetanp1m[j]/poros)*(h_prev[I] - hnp1m_)
        - (1.0/dz)*(Knp1m_up - Knp1m_down)
        - (3.0/(dt))*(thetanp1m[j] - thetan[j]);
  }
}
