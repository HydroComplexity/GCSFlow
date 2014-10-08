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


#include "../include/constantextern.h"

// --------------------------------------------------------------------
// vanGenuchten()
//    Estimates theta, K, and C from h using closed form formulation by
//    van Genuchten (1980).
//    Input: Vector pressure head h; Scalar parameters n, m
//    Output: Vector theta, K, C.
// --------------------------------------------------------------------

__device__ void vanGenuchten(double *C, double *K, double *Ksat, double *theta,
    double *h, double n, int size) {
  double Se, h_, theta_, m;
  m = 1.0 - 1.0/n;

  for (int i = 0; i < size; i++) {
    // Convert pressure unit from [m] to [cm]
    h_ = h[i] * 100;

    // Compute the volumetric moisture content [eqn 21]
    if (h_ < 0) {
      // for unsaturated soil conditions
      theta_ = (theta_S - theta_R)/pow(1.0 + pow((-h_*alpha), n), m) + theta_R;
    } else {
      // for saturated soil condition h_ >= 0
      theta_ = theta_S;
    }
    theta[i] = theta_;

    // Compute the effective saturation [eqn 2]
    Se = (theta_ - theta_R)/(theta_S - theta_R);

    // Compute the hydraulic conductivity [eqn 8] - [Convert to unit: m/hr
    K[i] = Ksat[i] * sqrt(Se) * (1.0 - pow(1.0-pow(Se, 1.0/m), m))
           * (1.0 - pow(1.0-pow(Se, 1.0/m), m));

    // Compute the specific moisture storage derivative of eqn (21).
    // So we have to calculate C = d(theta)/dh. Then the unit is converted
    // into [1/m].
    if (h_ < 0) {
      C[i] = 100 * -alpha * n * -1 * (1.0/n-1.0)*pow(alpha*abs(h_), n-1)
            * (theta_R-theta_S) * pow(pow(alpha*abs(h_), n)+1, 1.0/n-2.0);
    } else {
      C[i] = 0.0;
    }
  }
}


// -----------------------------------------------------------------------------
// :::::::::::::::::::::::::::::::: END OF FILE ::::::::::::::::::::::::::::::::
// -----------------------------------------------------------------------------