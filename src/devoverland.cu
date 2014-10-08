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


// This file provides device functions for overland flow model.Device functions
// are linked and called from global overland flow functions 

#include "../include/constants.h"
#include "../include/constantextern.h"


// --------------------------------------------------------------------
// SweXMatrix()
// SweYMatrix()
//    Creates tridiagonal matrix A and vector RHS for SWE (See Lal 1998, JHE)
//    The ADI is implemented in x and y directions for half time step 0.5*dt
//    each. The order of x & y is reversed every time step.
//    SweXMatrix() for horizontal x-direction calculation.
//    SweyMatrix() for horizontal y-direction calculation.
//
//    Input...
//        * Hs: Water elevation [L]
//        * h: Water depth [L]
//        * mann: Manning's coefficient [L^{-1/3} T]
//        * eRF: effective rainfall [L/T]
//        * ET: Evaporation [L/T]
//        * IN: Infiltration [L/T]
//        * Scalar indexes
//
//    Output...
//        Matrix A and vector RHS d
//        * a: sub-diagonal of A
//        * b: main diagonal of A
//        * c: sup-diagonal of A
//        * d: RHS vector
// --------------------------------------------------------------------
__device__ void SweXMatrix(double *a_x, double *b_x, double *c_x, double *d_x,
                int N, int j, int t, double *Hs, double *h, double *mann,
                double *eRF, double *IN, double *ET, double *K2w, double *K2e) {
  double Keast, Kwest, Knorth, Ksouth, Sse, Ssw, Ssn, Sss;
  double INX, ETX, PPTX, dA, h_, h_p1, h_m1, h_pN, h_mN;
  double mann_, mann_p1, mann_m1, mann_pN, mann_mN;
  double Hs_, Hs_p1, Hs_pN, Hs_pNp1, Hs_pNm1, Hs_m1, Hs_mN, Hs_mNp1, Hs_mNm1;
  int I, k;

  for (int i = 1; i < N-1; i++) {
    k = i-1;
    I = j*N+i;

    // Transfer data from global to local memory for speed
    // p1: +1, m1: -1, pN: +N, mN: -N
    // water depth h
    h_ = h[I];
    h_p1 = h[I+1];
    h_pN = h[I+N];
    h_m1 = h[I-1];
    h_mN = h[I-N];

    // Water elevation Hs=h+z
    Hs_ = Hs[I];
    Hs_p1 = Hs[I+1];
    Hs_pN = Hs[I+N];
    Hs_pNp1 = Hs[I+N+1];
    Hs_pNm1 = Hs[I+N-1];
    Hs_m1 = Hs[I-1];
    Hs_mN = Hs[I-N];
    Hs_mNp1 = Hs[I-N+1];
    Hs_mNm1 = Hs[I-N-1];

    // Manning's coefficients
    mann_ = mann[I];
    mann_p1 = mann[I+1];
    mann_pN = mann[I+N];
    mann_m1 = mann[I-1];
    mann_mN = mann[I-N];

    // Calculate slope Ss east and K east - Eqn (11)
    Sse = sqrt(pow((Hs_p1 - Hs_)/dx, 2.0)
        + pow((Hs_pNp1+Hs_pN-Hs_mNp1-Hs_mN)/(4*dy), 2.0));

    // If Hs in right cell is greater
    if (Hs_p1 > Hs_) {
      // Check water depth condition
      if (h_p1 > hmin && (h_p1+h_)/2.0 > hmin && abs(Sse) > delta) {
        Keast = pow((h_p1+h_)/2.0, 5.0/3.0)/(0.5*(mann_p1+mann_)*sqrt(Sse));
      } else {
        Keast = K0;
      }

    } else {
      // Check water depth condition
      if (h_ > hmin && (h_p1+h_)/2.0 > hmin && abs(Sse) > delta) {
        Keast = pow((h_p1 + h_)/2.0, 5.0/3.0)/(0.5*(mann_p1+mann_)*sqrt(Sse));
      } else {
        Keast = K0;
      }
    }
    K2e[I] = Keast;    // return to global memory

    // Calculate slope Ss west and K west - Eqn (11)
    Ssw = sqrt(pow((Hs_ - Hs_m1)/dx, 2.0)
        + pow((Hs_pN+Hs_pNm1-Hs_mN-Hs_mNm1)/(4*dy), 2.0));

    // If Hs in left cell is greater
    if (Hs_m1 > Hs_) {
      // Check water depth condition
      if (h_m1 > hmin && (h_m1+h_)/2.0 > hmin && abs(Ssw) > delta) {
        Kwest = pow((h_m1+h_)/2.0, 5.0/3.0)/(0.5*(mann_m1+mann_)*sqrt(Ssw));
      } else {
        Kwest = K0;
      }

    } else {
      // Check water depth condition
      if (h_ > hmin && (h_m1+h_)/2.0 > hmin && abs(Ssw) > delta) {
        Kwest = pow((h_m1+h_)/2.0, 5.0/3.0)/(0.5*(mann_m1+mann_)*sqrt(Ssw));
      } else {
        Kwest = K0;
      }
    }
    K2w[I] = Kwest;

    // Calculate slope Ss north and K north - Eqn (11)
    Ssn = sqrt(pow((Hs_pN - Hs_)/dy, 2.0)
        + pow((Hs_pNp1+Hs_p1-Hs_pNm1-Hs_m1)/(4*dx), 2.0));

    // If Hs in upper (north) cell is greater
    if (Hs_pN > Hs_) {
      // Check water depth condition
      if (h_pN > hmin && (h_pN+h_)/2.0 > hmin && abs(Ssn) > delta) {
        Knorth = pow((h_pN+h_)/2.0, 5.0/3.0)/(0.5*(mann_pN+mann_)*sqrt(Ssn));
      } else {
        Knorth = K0;
      }
    } else {
      // check water depth condition
      if (h_ > hmin && (h_pN+h_)/2.0 > hmin && abs(Ssn) > delta) {
        Knorth = pow((h_pN+h_)/2.0, 5.0/3.0)/(0.5*(mann_pN+mann_)*sqrt(Ssn));
      } else {
        Knorth = K0;
      }
    }

    // Calculate slope Ss south and K south - Eqn (11)
    Sss = sqrt(pow((Hs_ - Hs_mN)/dy, 2.0)
        + pow((Hs_p1+Hs_mNp1-Hs_m1-Hs_mNm1)/(4*dx), 2.0));

    // If Hs in lower (south) cell is greater
    if (Hs_mN > Hs_) {
      // Check water depth condition
      if (h_mN > hmin && (h_mN+h_)/2.0 > hmin && abs(Sss) > delta) {
        Ksouth = pow((h_mN+h_)/2.0, 5.0/3.0)/(0.5*(mann_mN+mann_)*sqrt(Sss));
      } else {
        Ksouth = K0;
      }
    } else {
      // Check water depth condition
      if (h_ > hmin && (h_mN+h_)/2.0 > hmin && abs(Sss) > delta) {
        Ksouth = pow((h_mN+h_)/2.0, 5.0/3.0)/(0.5*(mann_mN+mann_)*sqrt(Sss));
      } else {
        Ksouth = K0;
      }
    }

    dA = dx*dy;

    // Sub-diagonal, main-diagonal, sup-diagonal
    a_x[k] = -0.5 * dt/dA * Kwest;
    b_x[k] = 1 + 0.5 * dt/dA * (Keast+Kwest);
    c_x[k] = -0.5 * dt/dA * Keast;

    // IN and ET must be smaller than water depth h
    PPTX = eRF[I];    
    if (h_ > 0.0) {
      if (h_ > ET[t]) {
        ETX = ET[t];
      } else {
        ETX = h_;
      }
      
      if (h_ > IN[I]) {
        INX = IN[I];
      } else {
        INX = h_;
      }      
    } else {
      ETX = 0.0;
      INX = 0.0;
    }

    // Vector right hand side
    if (i == 1) {
      d_x[k] = 0.5*dt/dA*Ksouth*Hs_mN + (1-0.5*dt/dA*(Ksouth+Knorth))*Hs_
          + 0.5*dt/dA*Knorth*Hs_pN + 0.5*dt/dA*Kwest*Hs_m1 + 0.5*(PPTX-INX-ETX);
    } else if (i == N-2) {
      d_x[k] = 0.5*dt/dA*Ksouth*Hs_mN + (1-0.5*dt/dA*(Ksouth+Knorth))*Hs_
          + 0.5*dt/dA*Knorth*Hs_pN + 0.5*dt/dA*Keast*Hs_p1 + 0.5*(PPTX-INX-ETX);
    } else {
      d_x[k] = 0.5*dt/dA*Ksouth*Hs_mN + (1-0.5*dt/dA*(Ksouth+Knorth))*Hs_
          + 0.5*dt/dA*Knorth*Hs_pN + 0.5*(PPTX-INX-ETX);
    }
  }
}

__device__ void SweYMatrix(double *a_y, double *b_y, double *c_y, double *d_y,
      int M, int N, int i, int t, double *Hs, double *h, double *mann,
      double *eRF, double *IN, double *ET, double *K2n, double *K2s) {
  double Keast, Kwest, Knorth, Ksouth;
  double Sse, Ssw, Ssn, Sss;
  double INY, ETY, PPTY, dA;
  int I, k;
  double h_, h_p1, h_m1, h_pN, h_mN, mann_, mann_p1, mann_m1, mann_pN, mann_mN;
  double Hs_, Hs_p1, Hs_pN, Hs_pNp1, Hs_pNm1, Hs_m1, Hs_mN, Hs_mNp1, Hs_mNm1;

  for (int j = 1; j < M-1; j++) {
    k = j-1;
    I = j*N+i;

    // Transfer data from global to local memory for speed
    // p1: +1, m1: -1, pN: +N, mN: -N
    // water depth h
    h_ = h[I];
    h_p1 = h[I+1];
    h_pN = h[I+N];
    h_m1 = h[I-1];
    h_mN = h[I-N];

    // water elevation Hs=h+z
    Hs_ = Hs[I];
    Hs_p1 = Hs[I+1];
    Hs_pN = Hs[I+N];
    Hs_pNp1 = Hs[I+N+1];
    Hs_pNm1 = Hs[I+N-1];
    Hs_m1 = Hs[I-1];
    Hs_mN = Hs[I-N];
    Hs_mNp1 = Hs[I-N+1];
    Hs_mNm1 = Hs[I-N-1];

    // Manning's coefficients
    mann_ = mann[I];
    mann_p1 = mann[I+1];
    mann_pN = mann[I+N];
    mann_m1 = mann[I-1];
    mann_mN = mann[I-N];

    // Calculate slope Ss east and K east - Eqn (11)
    Sse = sqrt(pow((Hs_p1 - Hs_)/dx, 2.0)
        + pow((Hs_pNp1+Hs_pN-Hs_mNp1-Hs_mN)/(4*dy), 2.0));

    // If Hs in right cell is greater
    if (Hs_p1 > Hs_) {
      // Check water depth condition
      if (h_p1 > hmin && (h_p1+h_)/2.0 > hmin && abs(Sse) > delta) {
        Keast = pow((h_p1+h_)/2.0, 5.0/3.0)/(0.5*(mann_p1+mann_)*sqrt(Sse));
      } else {
        Keast = K0;
      }
    } else {
      // Check water depth condition
      if (h_ > hmin && (h_p1+h_)/2.0 > hmin && abs(Sse) > delta) {
        Keast = pow((h_p1 + h_)/2.0, 5.0/3.0)/(0.5*(mann_p1+mann_)*sqrt(Sse));
      } else {
        Keast = K0;
      }
    }

    // Calculate slope Ss west and K west - Eqn (11)
    Ssw = sqrt(pow((Hs_ - Hs_m1)/dx, 2.0)
        + pow((Hs_pN+Hs_pNm1-Hs_mN-Hs_mNm1)/(4*dy), 2.0));

    // If Hs in left cell is greater
    if (Hs_m1 > Hs_) {
      // Check water depth condition
      if (h_m1 > hmin && (h_m1+h_)/2.0 > hmin && abs(Ssw) > delta) {
        Kwest = pow((h_m1+h_)/2.0, 5.0/3.0)/(0.5*(mann_m1+mann_)*sqrt(Ssw));
      } else {
        Kwest = K0;
      }
    } else {
      // check water depth condition
      if (h_ > hmin && (h_m1+h_)/2.0 > hmin && abs(Ssw) > delta) {
        Kwest = pow((h_m1+h_)/2.0, 5.0/3.0)/(0.5*(mann_m1+mann_)*sqrt(Ssw));
      } else {
        Kwest = K0;
      }
    }

    // Calculate slope Ss north and K north - Eqn (11)
    Ssn = sqrt(pow((Hs_pN - Hs_)/dy, 2.0)
        + pow((Hs_pNp1+Hs_p1-Hs_pNm1-Hs_m1)/(4*dx), 2.0));

    // If Hs in upper (north) cell is greater
    if (Hs_pN > Hs_) {
      // Check water depth condition
      if (h_pN > hmin && (h_pN+h_)/2.0 > hmin && abs(Ssn) > delta) {
        Knorth = pow((h_pN+h_)/2.0, 5.0/3.0)/(0.5*(mann_pN+mann_)*sqrt(Ssn));
      } else {
        Knorth = K0;
      }

    } else {
      // Check water depth condition
      if (h_ > hmin && (h_pN+h_)/2.0 > hmin && abs(Ssn) > delta) {
        Knorth = pow((h_pN+h_)/2.0, 5.0/3.0)/(0.5*(mann_pN+mann_)*sqrt(Ssn));
      } else {
        Knorth = K0;
      }
    }
    K2n[I] = Knorth;

    // Calculate slope Ss south and K south - Eqn (11)
    Sss = sqrt(pow((Hs_ - Hs_mN)/dy, 2.0)
        + pow((Hs_p1+Hs_mNp1-Hs_m1-Hs_mNm1)/(4*dx), 2.0));

    // If Hs in lower (south) cell is greater
    if (Hs_mN > Hs_) {
      // Check water depth condition
      if (h_mN > hmin && (h_mN+h_)/2.0 > hmin && abs(Sss) > delta) {
        Ksouth = pow((h_mN+h_)/2.0, 5.0/3.0)/(0.5*(mann_mN+mann_)*sqrt(Sss));
      } else {
        Ksouth = K0;
      }

    } else {
      // check water depth condition
      if (h_ > hmin && (h_mN+h_)/2.0 > hmin && abs(Sss) > delta) {
        Ksouth = pow((h_mN+h_)/2.0, 5.0/3.0)/(0.5*(mann_mN+mann_)*sqrt(Sss));
      } else {
        Ksouth = K0;
      }
    }
    K2s[I] = Ksouth;

    dA = dx*dy;

    // Sub-diagonal, main-diagonal, sup-diagonal
    a_y[k] = -0.5 * dt/(dx*dy) * Ksouth;
    b_y[k] = 1 + 0.5 * dt/(dx*dy) * (Knorth+Ksouth);
    c_y[k] = -0.5 * dt/(dx*dy) * Knorth;

    // IN, and ET must be smaller than water depth h
    PPTY = eRF[I];
    //INY = IN[I];
    if (h_ > 0.0) {
      if (h_ > ET[t]) {
          ETY = ET[t];
      } else {
          ETY = h_;
      }

      if (h_ > IN[I]) {
          INY = IN[I];
      } else {
          INY = h_;
      }      
    } else {
      ETY = 0.0;
      INY = 0.0;
    }

    if (j == 1) {
      d_y[k] = 0.5*dt/dA*Kwest*Hs_m1 + (1-0.5*dt/dA*(Kwest+Keast))*Hs_
          + 0.5*dt/dA*Keast*Hs_p1 + 0.5*dt/dA*Ksouth*Hs_mN + 0.5*(PPTY-INY-ETY);

    } else if (j == M-2) {
      d_y[k] = 0.5*dt/dA*Kwest*Hs_m1 + (1-0.5*dt/dA*(Kwest+Keast))*Hs_
          + 0.5*dt/(dx*dy)*Keast*Hs_p1 + 0.5*dt/dA*Knorth*Hs_pN + 0.5*(PPTY-INY-ETY);

    } else {
      d_y[k] = 0.5*dt/dA*Kwest*Hs_m1 + (1-0.5*dt/dA*(Kwest+Keast))*Hs_
          + 0.5*dt/dA*Keast*Hs_p1 + 0.5*(PPTY-INY-ETY);
    }
  }
}


// -----------------------------------------------------------------------------
// :::::::::::::::::::::::::::::::: END OF FILE ::::::::::::::::::::::::::::::::
// -----------------------------------------------------------------------------