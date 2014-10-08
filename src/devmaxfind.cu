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
// MaxError()
//    Finds the maximum value of the input array (X & Y directions).
// --------------------------------------------------------------------
__device__ double MaxError(double *Array_in, int SIZE) {
  int i;
  double maxval = abs(Array_in[0]);

  for (i = 1; i < SIZE; ++i) {
    if ( abs(Array_in[i]) > maxval ) {
      maxval = abs(Array_in[i]);
    }
  }
  return maxval;
}


// --------------------------------------------------------------------
// MaxErrorZ()
//    Finds the maximum value of the input array in Z direction.
// --------------------------------------------------------------------
__device__ double MaxErrorZ(double *Array_in, int size) {
  int i;
  double maxval = abs(Array_in[1]);

  for (i = 2; i < size-1; ++i) {
    if ( abs(Array_in[i]) > maxval ) {
      maxval = abs(Array_in[i]);
    }
  }
  return maxval;
}


// -----------------------------------------------------------------------------
// :::::::::::::::::::::::::::::::: END OF FILE ::::::::::::::::::::::::::::::::
// -----------------------------------------------------------------------------