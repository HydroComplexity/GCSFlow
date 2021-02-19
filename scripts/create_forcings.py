import numpy as np
from netCDF4 import Dataset
import scipy.io as sio
import os

nc_fid = Dataset('forcings_7536.nc','r')
PPT_in = nc_fid.variables['Precipitation'][:]
num_steps_in = len(PPT_in)
ET_in = 0.001 * np.ones((num_steps_in))

# Check/create folder and forcing file name
folder = 'output'
if not os.path.exists(folder):
    os.makedirs(folder)

forcing = Dataset(folder + '/forcings_goosebw.nc', 'w', format='NETCDF4')
forcing.description = 'Forcing Data for GCSFlow'

# Set up dimension
forcing.createDimension('time', num_steps_in)	# Time series variable
forcing.createDimension('scalar', 1)			# Scalar variable


"""
Create variables in the netcdf file
	var = netcdf.createVariable('Var_name', 'var_type', ('dimension_type'))
"""
PPT       = forcing.createVariable('PPT_in', 'f8', ('time'))
ET        = forcing.createVariable('Evap_in', 'f8', ('time'))
num_steps = forcing.createVariable('num_steps', 'i4', ('scalar'))


# Assign data to variables in NetCDF file
PPT[:]       = PPT_in
ET[:]        = ET_in
num_steps[:] = num_steps_in

# Close the file
forcing.close()
