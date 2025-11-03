# -*- coding: utf-8 -*-
"""
Created on Tue Apr  9 20:08:47 2024

@author: PanosotG
"""

import h5py
from pathlib import Path
import numpy as np
import numpy.ma as ma
import datetime as dt
# import itertools
from tuflow_to_hecras_mesh import HecRasMesh

# sample_hec_file = r'Sample_HEC_data\Hydraulic_Data\Calibrated\Calibrated.hdf'
# f_hec = h5py.File(sample_hec_file, 'r')

tuflow_res_file = r'D:\hydraulicmodels\KolanRiver_135002A\50_Hydraulics\510_TUFLOW\results\NC_HYD002_HR02_BCLS\KolanRiver_BCLS_HR02_HYD002_HGS_005_NC_CC.nc'  # TODO: user input
tlf_filepath = r'D:\hydraulicmodels\KolanRiver_135002A\50_Hydraulics\510_TUFLOW\runs\log\NC_HYD002_HR02_BCLS\KolanRiver_BCLS_HR02_HYD002_HGS_005_NC.tlf'
out_folder = r'C:\PythonProjects\tuflow_to_lifesim\Output'

tuflow_run_id = Path(tuflow_res_file).stem      # File name (without extension). Requires Python 3.4+

test_num = '13'
output_id = tuflow_run_id + test_num

area_id = 'Domain_001'
    
def tuflow_to_ras_hdf():
    output_file = r'Output\{}.hdf'.format(output_id)                        # TODO: user input directory
    
    with h5py.File(output_file, 'w') as f_out:
        with h5py.File(tuflow_res_file, 'r') as f_tuflow:
            
            # Write file attributes:
            f_out.attrs['File Type'] = f_tuflow.attrs['references']         # e.g. b'TUFLOW NetCDF Raster Output Format (https://wiki.tuflow.com/TUFLOW_NetCDF_Raster_Output_Format)'
            f_out.attrs['File Version'] = f_tuflow.attrs['source']          # e.g. b'TUFLOW Build: 2023-03-AD-iSP-w64'
            f_out.attrs['Projection'] = f_tuflow['crs'].attrs['crs_wkt']    # e.g. b'PROJCS["GDA2020_MGA_Zone_56",GEOGCS["GCS_GDA2020",DATUM["GDA2020",SPHEROID["GRS_1980",6378137.0,298.257222101]],PRIMEM["Greenwich",0.0],UNIT["Degree",0.0174532925199433]],PROJECTION["Transverse_Mercator"],PARAMETER["False_Easting",500000.0],PARAMETER["False_Northing",10000000.0],PARAMETER["Central_Meridian",153.0],PARAMETER["Scale_Factor",0.9996],PARAMETER["Latitude_Of_Origin",0.0],UNIT["Meter",1.0]]'
            f_out.attrs['Units System'] = np_bytes('SI Units')              # TODO: hard-coded for now
    
            # Read data axes
            time_in_hours = np.array(f_tuflow['time'])
            model_start_time = get_start_time(f_tuflow)
            # x = np.array(f_tuflow['x'])         # X axis
            # y = np.array(f_tuflow['y'])         # Y axis
            
            # Read hydraulic time series grids
            nc_water_level = np.array(f_tuflow['water_level'])              # 3 dim array of time series water level grid
            max_water_level = np.array(f_tuflow['maximum_water_level'])     # 3 dim array (time length = 1)
            max_depth = np.array(f_tuflow['maximum_depth'])
            nc_velocity_dir = np.array(f_tuflow['direction_of_velocity'])
            nc_velocity_mag = np.array(f_tuflow['magnitude_of_velocity'])
    
        # Write run information
        f_out.create_group('Plan Data/Plan Information')
        f_out['Plan Data/Plan Information'].attrs['Plan ShortID'] = np_bytes(output_id)
        
        ### Write geometry
        names = f_out.create_dataset('Geometry/2D Flow Areas/Names', shape = (1, ), dtype = 'S16', data = area_id.ljust(16), chunks=True)
        names.attrs['Column'] = np_bytes('Name')
        names.attrs['Row'] = np_bytes('Feature')
        
        mesh_group = f_out.create_group('Geometry/2D Flow Areas/{}/'.format(area_id))        
        mesh_geometry = HecRasMesh(tlf_filepath)
        for key, value in mesh_geometry.data.items():
            mesh_group.create_dataset(key, data = value, chunks=True, compression="gzip", compression_opts=1)
        
        # Write time axis
        time_in_days = time_in_hours / 24.0
        time_axis = f_out.create_group('Results/Unsteady/Output/Output Blocks/Base Output/Unsteady Time Series/')
        time_axis.create_dataset('Time', data = time_in_days, chunks=True, compression="gzip", compression_opts=1)
        time_axis['Time'].attrs['Number of actual Time Steps'] = [len(time_in_days)]
        time_axis['Time'].attrs['Time Stamp'] = [np_bytes(model_start_time.strftime('%d%b%Y %H%M'))]
        time_axis['Time'].attrs['Time'] = ['days']
        
        time_stamp = []
        for timestep in time_in_hours:
            stamp = model_start_time + dt.timedelta(hours = timestep)
            time_stamp.append(stamp.strftime('%d%b%Y %H%M').upper())
        
        time_axis.create_dataset('Time Date Stamp', data = time_stamp, chunks=True, compression="gzip", compression_opts=1)
        
        # Write hydraulic data
        data = []
        for grid in nc_water_level:                                     # loop for each time step
            grid = np.maximum(grid, 150)                                  # replace -999.0 with ground elevation (except dry cells)
            row = grid.flatten()
            data.append(row)
        hdf_water_level = np.vstack(data)                               # row axis = time; column axis = cells
        hdf_water_level = np.pad(hdf_water_level, ((0, 0), (0 , mesh_geometry.number_external_cells)))
        
        hdf_velocity = mesh_geometry.zero_velocity(len(time_in_days))
        
        # Reformat velocity data
        # data = []
        # for grid in nc_velocity_mag:
        #     wet_cells = ma.masked_array(grid, mask = dry_cells)
        #     row = wet_cells.compressed()
        #     data.append(row)
        # hdf_velocity = np.vstack(data)
            
        # Write hydraulic data
        hydraulics = f_out.create_group('Results/Unsteady/Output/Output Blocks/Base Output/Unsteady Time Series/2D Flow Areas/{}/'.format(area_id))
        f_out['Results/Unsteady/'].attrs['Short ID'] = np_bytes(output_id)
        hydraulics.create_dataset('Water Surface', data = hdf_water_level, chunks=True, compression="gzip", compression_opts=1)
        hydraulics['Water Surface'].attrs['WSEL'] = np_bytes('Meters')
        hydraulics.create_dataset('Face Velocity', data = hdf_velocity, chunks=True, compression="gzip", compression_opts=1)
        hydraulics['Face Velocity'].attrs['Velocity'] = np_bytes('m/s')
    return


def superseded_code():
    output_file = r'Output\{}.hdf'.format(output_id)
    with h5py.File(output_file, 'w') as f_out:       
        
        # Process hydraulic data
        dry_cells = max_water_level[0] < 0
        
        # Minimum water level surface as proxy for ground surface
        ground = ma.masked_less_equal(nc_water_level, -999.0).min(axis = 0).data
        
        # Reformat water level data
        data = []
        for grid in nc_water_level:                                     # loop for each time step
            grid = np.maximum(grid, ground)                             # replace -999.0 with ground elevation (except dry cells)
            wet_cells = ma.masked_array(grid, mask = dry_cells)         # mask dry cells
            row = wet_cells.compressed()                                # wet cells to row
            data.append(row)
        hdf_water_level = np.vstack(data)                               # row axis = time; column axis = cells
            
        
        # Reformat velocity data
        data = []
        for grid in nc_velocity_mag:
            wet_cells = ma.masked_array(grid, mask = dry_cells)
            row = wet_cells.compressed()
            data.append(row)
        hdf_velocity = np.vstack(data)
            
        # Write hydraulic data
        hydraulics = f_out.create_group('Results/Unsteady/Output/Output Blocks/Base Output/Unsteady Time Series/2D Flow Areas/A/')
        f_out['Results/Unsteady/'].attrs['Short ID'] = np_bytes(output_id)
        hydraulics.create_dataset('Water Surface', data = hdf_water_level)
        hydraulics['Water Surface'].attrs['WSEL'] = np_bytes('Meters')
        hydraulics.create_dataset('Face Velocity', data = hdf_velocity)
        hydraulics['Face Velocity'].attrs['Velocity'] = np_bytes('m/s')
        
        # # Write geometry
        # yx = np.array(list(itertools.product(y, x)))
        # xy = yx[:, [1, 0]]                  # Swap X-Y columns
        # mask = dry_cells.flatten()
        # xy = xy[~mask]

        # coords = f_out.create_group('Geometry/2D Flow Areas/A/').create_dataset('Cells Center Coordinate', data = xy)
        # coords.attrs['Column'] = ['X', 'Y']
        # coords.attrs['Row'] = np_bytes('Cell')    
        
        return

def np_bytes(string):
    return np.bytes_(bytes(string,'utf-8'))
    
def get_start_time(netcdf_file):
    start_time = netcdf_file['time'].attrs['units']                 # e.g. b'hours since 2000-01-01 00:00'
    start_time = start_time[12:]                                    # Drop first 12 chars ("hours since ")
    start_time = np.datetime64(start_time)                          # Convert to numpy date time
    start_time = start_time.astype(dt.datetime)
    # start_time = start_time.strftime('%d%b%Y %H%M')                 # Convert to string
    # start_time = np_bytes(start_time)                               # Convert to numpy.bytes_
    
    return start_time
    