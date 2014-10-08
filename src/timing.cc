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
// GetTime()
//    Gets the time if called. This functions is used for measuring
//    time durations in Windows or Unix
// --------------------------------------------------------------------
#if defined(WIN32) || defined(_WIN32) || defined(__WIN32) || defined(CYGWIN)
#include <windows.h>
double GetTime() {
  SYSTEMTIME time;
  GetSystemTime(&time);
  return static_cast<double>(time.wSecond + (time.wMilliseconds/1000.0));
}

#else
#include <sys/time.h>
#include <cstdlib>
double GetTime() {
  struct timeval time;
  gettimeofday(&time, NULL);
  return static_cast<double>(time.tv_sec+(time.tv_usec/1000000.0));
}
#endif


