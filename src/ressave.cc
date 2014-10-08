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


// This file provide funcitons for Saving Results/Outputs from GCSFlow model.

#include <netcdf.h>     // For using NetCDF format


// --------------------------------------------------------------------
// SaveOutput2D()
//    Saves two 2D double type variables into a NetCDF file.
// --------------------------------------------------------------------
void SaveOutput2D(const char *file, int Nx, int My, const char *data_name1,
    double *data_in1, const char *data_name2, double *data_in2) {
  int ncid, x_dimid, y_dimid, varid1, varid2, dimids[2];

  nc_create(file, NC_CLOBBER, &ncid);
  nc_def_dim(ncid, "x", Nx, &x_dimid);
  nc_def_dim(ncid, "y", My, &y_dimid);
  dimids[0] = y_dimid;
  dimids[1] = x_dimid;
  nc_def_var(ncid, data_name1, NC_DOUBLE, 2, dimids, &varid1);
  nc_def_var(ncid, data_name2, NC_DOUBLE, 2, dimids, &varid2);
  nc_enddef(ncid);
  nc_put_var_double(ncid, varid1, &data_in1[0]);
  nc_put_var_double(ncid, varid2, &data_in2[0]);
  nc_close(ncid);
}

// --------------------------------------------------------------------
// SaveOneIntOutput2D()
//    Saves one/single integer type 2D variable into a NetCDF file.
// --------------------------------------------------------------------
void SaveOneIntOutput2D(const char *file, int Nx, int My, const char *data_name,
    int *data_in) {
  int ncid, x_dimid, y_dimid, varid, dimids[2];

  nc_create(file, NC_CLOBBER, &ncid);
  nc_def_dim(ncid, "x", Nx, &x_dimid);
  nc_def_dim(ncid, "y", My, &y_dimid);
  dimids[0] = y_dimid;
  dimids[1] = x_dimid;
  nc_def_var(ncid, data_name, NC_INT, 2, dimids, &varid);
  nc_enddef(ncid);
  nc_put_var_int(ncid, varid, &data_in[0]);
  nc_close(ncid);
}

// --------------------------------------------------------------------
// SaveOneOutput2D()
//    Saves one/single double type 2D variable into a NetCDF file.
// --------------------------------------------------------------------
void SaveOneOutput2D(const char *file, int Nx, int My, const char *data_name,
    double *data_in) {
  int ncid, x_dimid, y_dimid, varid, dimids[2];

  nc_create(file, NC_CLOBBER, &ncid);
  nc_def_dim(ncid, "x", Nx, &x_dimid);
  nc_def_dim(ncid, "y", My, &y_dimid);
  dimids[0] = y_dimid;
  dimids[1] = x_dimid;
  nc_def_var(ncid, data_name, NC_DOUBLE, 2, dimids, &varid);
  nc_enddef(ncid);
  nc_put_var_double(ncid, varid, &data_in[0]);
  nc_close(ncid);
}

// --------------------------------------------------------------------
// SaveOutput3D()
//    Saves two 3D double type variables into a NetCDF file.
// --------------------------------------------------------------------
void SaveOutput3D(const char *file, int Nx, int My, int Pz, 
    const char *data_name1, double *data_in1, const char *data_name2, 
    double *data_in2) {
  int ncid, x_dimid, y_dimid, z_dimid, varid1, varid2, dimids[3];

  nc_create(file, NC_CLOBBER, &ncid);
  nc_def_dim(ncid, "x", Nx, &x_dimid);
  nc_def_dim(ncid, "y", My, &y_dimid);
  nc_def_dim(ncid, "z", Pz, &z_dimid);
  dimids[0] = z_dimid;
  dimids[1] = y_dimid;
  dimids[2] = x_dimid;
  nc_def_var(ncid, data_name1, NC_DOUBLE, 3, dimids, &varid1);
  nc_def_var(ncid, data_name2, NC_DOUBLE, 3, dimids, &varid2);
  nc_enddef(ncid);
  nc_put_var_double(ncid, varid1, &data_in1[0]);
  nc_put_var_double(ncid, varid2, &data_in2[0]);
  nc_close(ncid);
}


// --------------------------------------------------------------------
// SaveIterationMap()
//    Saves number of iteration (integer type) into a NetCDF file.
// --------------------------------------------------------------------
void SaveIterationMap(const char *file, int Nx, int My, int Pz, 
    const char *nameZ, int *Z_in, const char *nameX, int *X_in, 
    const char *nameY, int *Y_in) {
  int ncid, x_dimid, y_dimid, z_dimid, varidZ, varidX, varidY, dimids[2];

  nc_create(file, NC_CLOBBER, &ncid);
  nc_def_dim(ncid, "x", Nx, &x_dimid);
  nc_def_dim(ncid, "y", My, &y_dimid);
  nc_def_dim(ncid, "z", Pz, &z_dimid);

  dimids[0] = y_dimid;
  dimids[1] = x_dimid;
  nc_def_var(ncid, nameZ, NC_INT, 2, dimids, &varidZ);

  dimids[0] = z_dimid;
  dimids[1] = y_dimid;
  nc_def_var(ncid, nameX, NC_INT, 2, dimids, &varidX);

  dimids[0] = z_dimid;
  dimids[1] = x_dimid;
  nc_def_var(ncid, nameY, NC_INT, 2, dimids, &varidY);
  nc_enddef(ncid);

  nc_put_var_int(ncid, varidZ, &Z_in[0]);
  nc_put_var_int(ncid, varidX, &X_in[0]);
  nc_put_var_int(ncid, varidY, &Y_in[0]);
  nc_close(ncid);
}


// --------------------------------------------------------------------
// SaveVariables1D()
//    Saves 1D double type variables into a NetCDF file.
// --------------------------------------------------------------------
void SaveVariables1D(const char *file, int time_steps, const char *var_name, 
    double *var) {
  int ncid, t_dimid, varid, dimids[1];

  nc_create(file, NC_CLOBBER, &ncid);
  nc_def_dim(ncid, "x", time_steps, &t_dimid);
  dimids[0] = t_dimid;

  nc_def_var(ncid, var_name, NC_DOUBLE, 1, dimids, &varid);

  nc_enddef(ncid);

  nc_put_var_double(ncid, varid, &var[0]);
  nc_close(ncid);
}


// --------------------------------------------------------------------
// SaveK2D()
//    Saves K values in overland flow model into a NetCDF file.
// --------------------------------------------------------------------
void SaveK2D(const char *file, int My, int Nx, const char *nameW, 
    double *data_W, const char *nameE, double *data_E, const char *nameN, 
    double *data_N, const char *nameS, double *data_S) {
  int ncid, x_dimid, y_dimid, varidW, varidE, varidN, varidS, dimids[2];

  nc_create(file, NC_CLOBBER, &ncid);
  nc_def_dim(ncid, "x", Nx, &x_dimid);
  nc_def_dim(ncid, "y", My, &y_dimid);
  dimids[0] = y_dimid;
  dimids[1] = x_dimid;
  nc_def_var(ncid, nameW, NC_DOUBLE, 2, dimids, &varidW);
  nc_def_var(ncid, nameE, NC_DOUBLE, 2, dimids, &varidE);
  nc_def_var(ncid, nameN, NC_DOUBLE, 2, dimids, &varidN);
  nc_def_var(ncid, nameS, NC_DOUBLE, 2, dimids, &varidS);
  nc_enddef(ncid);
  nc_put_var_double(ncid, varidW, &data_W[0]);
  nc_put_var_double(ncid, varidE, &data_E[0]);
  nc_put_var_double(ncid, varidN, &data_N[0]);
  nc_put_var_double(ncid, varidS, &data_S[0]);
  nc_close(ncid);
}

// -----------------------------------------------------------------------------
// :::::::::::::::::::::::::::::::: END OF FILE ::::::::::::::::::::::::::::::::
// -----------------------------------------------------------------------------
