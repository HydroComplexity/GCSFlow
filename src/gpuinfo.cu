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


// In this file, we provide a function to check GPU capability using CUDA APIs.
// The file loads device information and check compute capability. Only devices
// that have compute capability 2.0 or higher are eligible for using in GCSFlow
// model due to the usage of double precision and some other functions.

#include <stdio.h>
#include <cmath>
#include "../include/colorcode.h"

// --------------------------------------------------------------------
// GPUGetInfo()
//    Gets and prints GPU/device information for compute capability.
//    Eligible devices require compute capability 2.0 or higher.
//
//    - num_device: Number of GPU devices 
//    - num_eligible_devices: Number of eligible devices (>=2.0)
// --------------------------------------------------------------------
void GPUGetInfo(int *device) {
  int num_devices;
  int num_eligible_devices;
  cudaDeviceProp prop;

  num_devices = 0;
  num_eligible_devices = 0;

  // Get number of GPU devices
  cudaGetDeviceCount(&num_devices);
  if (num_devices == 0) {
    printf("The system does not have a CUDA capable device\n");
    return;
  }

  // Print out GPU information and eligible devices
  printf("\n");
  printf(L1 " GPU DEVICES INFORMATION \n" RS);
  printf(" ----------------------------------------------------------------\n");
  printf("            Device         Clock rate     Compute                \n");
  printf(" ID          Name             (GHz)      Capability  Eligibility \n");
  printf(" ----------------------------------------------------------------\n");

  for (int i = 0; i < num_devices; i++) {
    cudaGetDeviceProperties(&prop, i);
    printf("%3d", i);
    printf("%20s        ", prop.name);
    printf("%4.2f      ", prop.clockRate* 1e-6f);
    printf("%6d.%d         ", prop.major, prop.minor);

    // Check compute capability
    if (prop.major >= 2) {
      num_eligible_devices += 1;
      printf("x");
    } else {
      printf("");
    }
    printf("\n");
  }
  printf(" ----------------------------------------------------------------\n");
  printf(L2 " Egibility requires compute capability 2.0 or higher\n" RS);
  printf("\n");
  printf(L1 " NUMBER OF ELIGIBLE DEVICES: " RS "%d\n", num_eligible_devices);
  printf("\n");

  // Return device info via pointer
  device[0] = num_devices;
  device[1] = num_eligible_devices;
}

// -----------------------------------------------------------------------------
// :::::::::::::::::::::::::::::::: END OF FILE ::::::::::::::::::::::::::::::::
// -----------------------------------------------------------------------------
