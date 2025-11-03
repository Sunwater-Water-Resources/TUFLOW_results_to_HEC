# -*- coding: utf-8 -*-
"""
Created on Wed Apr 17 12:43:41 2024

@author: PanosotG
"""

import h5py
from pathlib import Path
import numpy as np
import numpy.ma as ma
import datetime as dt
import itertools

sample_hec_file = r'C:\Users\panosotg\Work\RAS_Projects\2D Unsteady Flow Hydraulics\BaldEagleCrkMulti2D\BaldEagleDamBrk.p03.hdf'
joso_hec_file = r'Sample_HEC_data\Hydraulic_Data\zz_Calibrated\Calibrated.hdf'

tuflow_res_file = r'C:\PythonProjects\tuflow_to_lifesim\TUFLOW_kolan\KolanRiver_BCLS_HR02_HYD002_HGS_005_NC.nc'
# tuflow_run_id = Path(tuflow_res_file).stem      # File name (without extension). Requires Python 3.4+


test_num = '55'
output_id_baldeagle = 'BaldEagleDamBrk_' + test_num

output_id_joso = 'Joso_09'

ls_required_geometry = ['Cells Center Coordinate',
                        'Cells Face and Orientation Info',
                        'Cells Face and Orientation Values',
                        'Cells Minimum Elevation',
                        'FacePoints Coordinate',
                        'FacePoints Face and Orientation Info',
                        'FacePoints Face and Orientation Values',
                        'Faces Cell Indexes',
                        'Faces FacePoint Indexes',
                        'Faces Perimeter Info',
                        'Faces Perimeter Values',
                        'Perimeter']

ls_geo_attributes = ['Extents']

def copy_pieces(hec_ras = 'BaldEagleCr'):
    
    if hec_ras == 'BaldEagleCr':
        f_hec = h5py.File(sample_hec_file, 'r')
        output_id = output_id_baldeagle
        output_file = r'Output\c{}.hdf'.format(output_id)
        area_id = 'BaldEagleCr'
    elif hec_ras == 'Joso': 
        f_hec = h5py.File(joso_hec_file, 'r')
        output_id = output_id_joso
        output_file = r'Output\{}.hdf'.format(output_id_joso)
        area_id = 'Joso Area'
    else:
        raise Exception()
        
    
    f_out = h5py.File(output_file, 'a')
    f_tuflow1 = h5py.File(tuflow_res_file, 'r')
    
    # Write file attributes:
    f_out.attrs['File Type'] = f_tuflow1.attrs['references']         # e.g. b'TUFLOW NetCDF Raster Output Format (https://wiki.tuflow.com/TUFLOW_NetCDF_Raster_Output_Format)'
    f_out.attrs['File Version'] = f_tuflow1.attrs['source']          # e.g. b'TUFLOW Build: 2023-03-AD-iSP-w64'
    f_out.attrs['Projection'] = f_hec.attrs['Projection']    # e.g. b'PROJCS["GDA2020_MGA_Zone_56",GEOGCS["GCS_GDA2020",DATUM["GDA2020",SPHEROID["GRS_1980",6378137.0,298.257222101]],PRIMEM["Greenwich",0.0],UNIT["Degree",0.0174532925199433]],PROJECTION["Transverse_Mercator"],PARAMETER["False_Easting",500000.0],PARAMETER["False_Northing",10000000.0],PARAMETER["Central_Meridian",153.0],PARAMETER["Scale_Factor",0.9996],PARAMETER["Latitude_Of_Origin",0.0],UNIT["Meter",1.0]]'
    f_out.attrs['Units System'] = f_hec.attrs['Units System']
    
    f_tuflow1.close()
    
    # Write run information
    f_out.create_group('Plan Data/Plan Information')
    f_out['Plan Data/Plan Information'].attrs['Plan ShortID'] = np_bytes(output_id)
    
    # Write 
    geo = f_out.create_group('Geometry/2D Flow Areas/{}/'.format(area_id))
    
    names = f_out.create_dataset('Geometry/2D Flow Areas/Names', shape = (1, ), dtype = 'S16', data = area_id.ljust(16), chunks=True)
    names.attrs['Column'] = np_bytes('Name')
    names.attrs['Row'] = np_bytes('Feature')
        
    grp_hec_2d_area = f_hec['Geometry/2D Flow Areas/{}/'.format(area_id)]
    for dataset in ls_required_geometry:
        grp_hec_2d_area.copy(dataset, geo)
    
    geo.attrs['Extents'] = grp_hec_2d_area.attrs['Extents']
    
    # hydraulics
    res = f_out.create_group('Results/Unsteady/Output/Output Blocks/Base Output/Summary Output/2D Flow Areas/{}/'.format(area_id))
    f_hec.copy('Results/Unsteady/Output/Output Blocks/Base Output/Summary Output/2D Flow Areas/{}/Maximum Water Surface'.format(area_id), res)
    f_hec.copy('Results/Unsteady/Output/Output Blocks/Base Output/Summary Output/2D Flow Areas/{}/Maximum Face Velocity'.format(area_id), res)
    res = f_out.create_group('Results/Unsteady/Output/Output Blocks/Base Output/Unsteady Time Series/2D Flow Areas/{}/'.format(area_id))
    f_hec.copy('Results/Unsteady/Output/Output Blocks/Base Output/Unsteady Time Series/2D Flow Areas/{}/Water Surface'.format(area_id), res)
    f_hec.copy('Results/Unsteady/Output/Output Blocks/Base Output/Unsteady Time Series/2D Flow Areas/{}/Face Velocity'.format(area_id), res)
    f_hec.copy('Results/Unsteady/Output/Output Blocks/Base Output/Unsteady Time Series/Time', f_out['Results/Unsteady/Output/Output Blocks/Base Output/Unsteady Time Series/'])
    f_hec.copy('Results/Unsteady/Output/Output Blocks/Base Output/Unsteady Time Series/Time Date Stamp', f_out['Results/Unsteady/Output/Output Blocks/Base Output/Unsteady Time Series/'])
    
    f_out.close()
    f_hec.close()


def copy_strip_hdf():
    f_out = h5py.File(output_file, 'a')
    f_tuflow1 = h5py.File(tuflow_res_file, 'r')
    
    # Write file attributes:
    f_out.attrs['File Type'] = f_tuflow1.attrs['references']         # e.g. b'TUFLOW NetCDF Raster Output Format (https://wiki.tuflow.com/TUFLOW_NetCDF_Raster_Output_Format)'
    f_out.attrs['File Version'] = f_tuflow1.attrs['source']          # e.g. b'TUFLOW Build: 2023-03-AD-iSP-w64'
    f_out.attrs['Projection'] = f_hec.attrs['Projection']    # e.g. b'PROJCS["GDA2020_MGA_Zone_56",GEOGCS["GCS_GDA2020",DATUM["GDA2020",SPHEROID["GRS_1980",6378137.0,298.257222101]],PRIMEM["Greenwich",0.0],UNIT["Degree",0.0174532925199433]],PROJECTION["Transverse_Mercator"],PARAMETER["False_Easting",500000.0],PARAMETER["False_Northing",10000000.0],PARAMETER["Central_Meridian",153.0],PARAMETER["Scale_Factor",0.9996],PARAMETER["Latitude_Of_Origin",0.0],UNIT["Meter",1.0]]'
    f_out.attrs['Units System'] = f_hec.attrs['Units System']
    
    # Write run information
    f_out.create_group('Plan Data/Plan Information')
    f_out['Plan Data/Plan Information'].attrs['Plan ShortID'] = np_bytes(output_id)
    
    # Write time axis
    # time_axis = f_out.create_group('Results/Unsteady/Output/Output Blocks/Base Output/Unsteady Time Series/')    
    # f_hec.copy('Results/Unsteady/Output/Output Blocks/Base Output/Unsteady Time Series/Time', time_axis)
    # time_axis['Time'].attrs['Number of actual Time Steps'] = [433]
    
    # # Write hydraulic data
    # f_out['Results/Unsteady/'].attrs['Short ID'] = np_bytes(output_id)
    # hydraulics = f_out.create_group('Results/Unsteady/Output/Output Blocks/Base Output/Unsteady Time Series/2D Flow Areas/BaldEagleCr/')
    # f_hec.copy('Results/Unsteady/Output/Output Blocks/Base Output/Unsteady Time Series/2D Flow Areas/BaldEagleCr/Water Surface', hydraulics)
    # hydraulics['Water Surface'].attrs['WSEL'] = np_bytes('Meters')
    # f_hec.copy('Results/Unsteady/Output/Output Blocks/Base Output/Unsteady Time Series/2D Flow Areas/BaldEagleCr/Face Velocity', hydraulics)
    # hydraulics['Face Velocity'].attrs['Velocity'] = np_bytes('m/s')
    
    # hydraulics = f_out.create_group('Results/Unsteady/Output/Output Blocks/Base Output/Unsteady Time Series/2D Flow Areas/BaldEagleCr/')
    res = f_out.create_group('Results/Unsteady/')
    f_hec.copy('Results/Unsteady/Output/', res)
    # f_hec.copy('Geometry/', f_out)
    # 
    
    # Write 
    geo = f_out.create_group('Geometry/2D Flow Areas/')
    f_hec.copy('Geometry/2D Flow Areas/BaldEagleCr', geo)
    # geo.create_dataset('Attributes', np.array(['BaldEagleCr']))
    # f_hec.copy('Geometry/2D Flow Areas/Attributes', geo)
    
    # DO NOT uncomment
    # del f_out['Geometry/2D Flow Areas/BaldEagleCr/Cells Minimum Elevation']
    # del f_out['Geometry/2D Flow Areas/BaldEagleCr/Cells Face and Orientation Values']
    # del f_out['Geometry/2D Flow Areas/BaldEagleCr/Cells Face and Orientation Info']
    # del f_out['Geometry/2D Flow Areas/BaldEagleCr/Faces Perimeter Info']
    # del f_out['Geometry/2D Flow Areas/BaldEagleCr/Faces Perimeter Values']
    # del f_out['Geometry/2D Flow Areas/BaldEagleCr/Perimeter']
    # del f_out['Geometry/2D Flow Areas/BaldEagleCr/FacePoints Coordinate']
    # del f_out['Geometry/2D Flow Areas/BaldEagleCr/Faces Cell Indexes']
    # del f_out['Geometry/2D Flow Areas/BaldEagleCr/Faces FacePoint Indexes']
    # del f_out['Geometry/2D Flow Areas/BaldEagleCr/FacePoints Face and Orientation Info']
    # del f_out['Geometry/2D Flow Areas/BaldEagleCr/FacePoints Face and Orientation Values']
    
    del f_out['Geometry/2D Flow Areas/BaldEagleCr/Cells Center Manning\'s n']
    del f_out['Geometry/2D Flow Areas/BaldEagleCr/Cells FacePoint Indexes']
    del f_out['Geometry/2D Flow Areas/BaldEagleCr/Cells Surface Area']
    del f_out['Geometry/2D Flow Areas/BaldEagleCr/Cells Volume Elevation Info']
    del f_out['Geometry/2D Flow Areas/BaldEagleCr/Cells Volume Elevation Values']
    del f_out['Geometry/2D Flow Areas/BaldEagleCr/Infiltration']
    del f_out['Geometry/2D Flow Areas/BaldEagleCr/Percent Impervious']
    
    del f_out['Geometry/2D Flow Areas/BaldEagleCr/Faces Area Elevation Values']
    del f_out['Geometry/2D Flow Areas/BaldEagleCr/Faces Area Elevation Info']
    
    del f_out['Geometry/2D Flow Areas/BaldEagleCr/Faces Low Elevation Centroid']
    del f_out['Geometry/2D Flow Areas/BaldEagleCr/Faces Minimum Elevation']
    del f_out['Geometry/2D Flow Areas/BaldEagleCr/Faces NormalUnitVector and Length']
    
    del f_out['Geometry/2D Flow Areas/BaldEagleCr/FacePoints Cell Index Values']
    del f_out['Geometry/2D Flow Areas/BaldEagleCr/FacePoints Cell Info']
    del f_out['Geometry/2D Flow Areas/BaldEagleCr/FacePoints Is Perimeter']
    
    f_out.close()
    f_hec.close()
    f_tuflow1.close()
    
def asdfasdf():    
    output_file = r'Output\{}.hdf'.format(output_id)
    f_out = h5py.File(output_file, 'a')
    
    # Write file attributes:
    f_out.attrs['File Type'] = f_tuflow1.attrs['references']         # e.g. b'TUFLOW NetCDF Raster Output Format (https://wiki.tuflow.com/TUFLOW_NetCDF_Raster_Output_Format)'
    f_out.attrs['File Version'] = f_tuflow1.attrs['source']          # e.g. b'TUFLOW Build: 2023-03-AD-iSP-w64'
    f_out.attrs['Projection'] = f_tuflow1['crs'].attrs['crs_wkt']    # e.g. b'PROJCS["GDA2020_MGA_Zone_56",GEOGCS["GCS_GDA2020",DATUM["GDA2020",SPHEROID["GRS_1980",6378137.0,298.257222101]],PRIMEM["Greenwich",0.0],UNIT["Degree",0.0174532925199433]],PROJECTION["Transverse_Mercator"],PARAMETER["False_Easting",500000.0],PARAMETER["False_Northing",10000000.0],PARAMETER["Central_Meridian",153.0],PARAMETER["Scale_Factor",0.9996],PARAMETER["Latitude_Of_Origin",0.0],UNIT["Meter",1.0]]'
    f_out.attrs['Units System'] = np_bytes('SI Units')
    
    # Write run information
    f_out.create_group('Plan Data/Plan Information')
    f_out['Plan Data/Plan Information'].attrs['Plan ShortID'] = np_bytes(output_id)
    
    # Read time axis
    time_in_days = np.array(f_hec['Results/Unsteady/Output/Output Blocks/Base Output/Unsteady Time Series/Time'])
    
    # Write time axis
    time_axis = f_out.create_group('Results/Unsteady/Output/Output Blocks/Base Output/Unsteady Time Series/')
    time_axis.create_dataset('Time', data = time_in_days)
    time_axis['Time'].attrs['Number of actual Time Steps'] = [len(time_in_days)]
    # time_axis['Time'].attrs['Time Stamp'] = [get_timestamp(f_tuflow)]
    # time_axis['Time'].attrs['Variables'] = ['Time|Day', '']         # Mimics HEC-RAS attribute. Not sure why there's a blank
    
    hdf_water_level = np.array(f_hec['Results/Unsteady/Output/Output Blocks/Base Output/Unsteady Time Series/2D Flow Areas/BaldEagleCr/Water Surface'])
    
    # Write hydraulic data
    hydraulics = f_out.create_group('Results/Unsteady/Output/Output Blocks/Base Output/Unsteady Time Series/2D Flow Areas/A/')
    f_out['Results/Unsteady/'].attrs['Short ID'] = np_bytes(output_id)
    hydraulics.create_dataset('Water Surface', data = hdf_water_level)
    hydraulics['Water Surface'].attrs['WSEL'] = np_bytes('Meters')
    # hydraulics.create_dataset('Face Velocity', data = hdf_velocity)
    # hydraulics['Face Velocity'].attrs['Velocity'] = np_bytes('m/s')
    
    # Write geometry
    # xy = xy[~mask]
    
    # coords = f_out.create_group('Geometry/2D Flow Areas/A/').create_dataset('Cells Center Coordinate', data = xy)
    # coords.attrs['Column'] = ['X', 'Y']
    # coords.attrs['Row'] = np_bytes('Cell')    
    # Close
    f_out.close()
    f_hec.close()
    f_tuflow1.close()

def np_bytes(string):
    return np.bytes_(bytes(string,'utf-8'))
    
def get_timestamp(netcdf_file):
    start_time = netcdf_file['time'].attrs['units']                 # e.g. b'hours since 2000-01-01 00:00'
    start_time = start_time[12:]                                    # Drop first 12 chars ("hours since ")
    start_time = np.datetime64(start_time)                          # Convert to numpy date time
    start_time = start_time.astype(dt.datetime)
    start_time = start_time.strftime('%d%b%Y %H%M')                 # Convert to string
    start_time = np_bytes(start_time)                               # Convert to numpy.bytes_
    
    return start_time
    
    
