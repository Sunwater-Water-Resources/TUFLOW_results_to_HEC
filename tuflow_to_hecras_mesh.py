# -*- coding: utf-8 -*-
"""
Created on Sun Apr 21 13:36:19 2024

@author: PanosotG
"""

import numpy as np
import math
import pandas as pd
from pathlib import Path

# tlf_filepath = r'J:\BW_WaterResources\hydraulicmodels\callide\TUFLOW\runs\log\Callide_Debbie_211.tlf'
# tlf_filepath = r'D:\hydraulicmodels\KolanRiver_135002A\50_Hydraulics\510_TUFLOW\runs\log\NC_HYD002_HR02_BCLS\KolanRiver_BCLS_HR02_HYD002_HGS_005_NC.tlf'
# out_folder = r'C:\PythonProjects\tuflow_to_lifesim\TUFLOW_kolan\chk_gis'

class HecRasMesh:
    def __init__(self, tlf_filepath, out_folder = None):
        self.tlf_filepath = tlf_filepath
        self.n_star = 0
        self.m_star = 0
        self.data = {}
        self.number_external_cells = 0

    # def tuflow_as_mesh(self, out_folder): #tlf_filepath = tlf_filepath, out_folder = out_folder):
        
        mesh_geometry = {}
        tlf_filepath = self.tlf_filepath
        out_file_prefix = Path(tlf_filepath).stem
        if out_folder is not None:
            Path(out_folder).mkdir(parents=True, exist_ok=True)
        
        # parse .TLF file
        try:
            tlf = open(tlf_filepath, 'r') 
        except FileNotFoundError:
            print(r'Unable to open .TLF file: ' + tlf_filepath)
            return
        else:
            with tlf:
                for line in tlf:
                    line = line.rstrip()                                        # Remove newline char at end of line
                    # if line.startswith(r'Write Check Files == '):
                    #     check_prefix = line.split('== ')[1]
                    #     dom_check = check_prefix + '_dom_check_R.shp'
                    if line.lower().startswith('tgc>> cell size == '):
                        cell_size = line.split('== ')[1]
                        cell_size = float(cell_size)
                    elif line.startswith('Domain_001  : Grid Dimensions (N,M):'):
                        line = line.split(': ')[2].split(',')
                        n = int(line[0])                                        # Number of rows and columns in grid
                        m = int(line[1])
                    elif line.startswith('Domain_001  : Geographic Origin:'): 
                        line = line.split(': ')[2].split()
                        x0 = float(line[0])                                     # Grid origin coordinates
                        y0 = float(line[1])
                    elif line.startswith('Domain_001  : Geographic Rotation Angle: '):
                        line = line.split(': ')[2]
                        grid_rotation_angle = float(line)                       # The angle is in degrees relative to east (e.g. X-axis directly north would be 90).
                        grid_rotation_angle = math.radians(grid_rotation_angle) # Convert to radians
                    elif line.startswith('...Reduced output grid by '):
                        line = line.split('[')
                        lower_left_corner = line[1].split(']')[0].split(',')
                        top_right_corner = line[2].split(']')[0].split(',')
                        i_min = int(lower_left_corner[0]) - 1
                        j_min = int(lower_left_corner[1]) - 1
                        i_max = int(top_right_corner[0]) - 1
                        j_max = int(top_right_corner[1]) - 1
                        
        # # get 2d_domain vertices
        # try:
        #     dom = open(dom_check, 'r') 
        # except FileNotFoundError:
        #     print(r'Unable to open _dom_check file: ' + dom_check)
        # else:
        #     with dom:
        #         x0, y0, x1, y1, x3, y3 = 0  # TODO: code
                
        # Calculate coordinates
        deg90_as_rads = math.radians(90)
        x1 = x0 + n * cell_size * math.cos(grid_rotation_angle + deg90_as_rads)          # top left corner of the 2D grid.
        y1 = y0 + n * cell_size * math.sin(grid_rotation_angle + deg90_as_rads)
        
        x3 = x0 + m * cell_size * math.cos(grid_rotation_angle)                 # bottom right corner of the grid.
        y3 = y0 + m * cell_size * math.sin(grid_rotation_angle)
        
        xi = np.linspace(x0, x1, n + 1)
        xj = np.linspace(x0, x3, m + 1)
        yi = np.linspace(y0, y1, n + 1)
        yj = np.linspace(y0, y3, m + 1)
        
        x_shift = 0.5 * cell_size * (math.cos(grid_rotation_angle) + math.cos(grid_rotation_angle + deg90_as_rads))      # Position of cell centroid relative to cell bottom left corner
        y_shift = 0.5 * cell_size * (math.sin(grid_rotation_angle) + math.sin(grid_rotation_angle + deg90_as_rads))
        
        FacePoints_coords = []
        CellCenter_coords = []
        for i in range(i_min, i_max + 1):
            for j in range(j_min, j_max + 1):
                zh_point_x = -x0 + xi[i] + xj[j]                             # coordinate of ZH point (cell corner vertex)
                zh_point_y = -y0 + yi[i] + yj[j]
                FacePoints_coords.append([zh_point_x, zh_point_y])
                
                zc_point_x = zh_point_x + x_shift                            # coordinate of ZC point (cell centre)
                zc_point_y = zh_point_y + y_shift
                CellCenter_coords.append([zc_point_x, zc_point_y])
                
            zh_point_x = -x0 + xi[i] + xj[j_max + 1]                         # coordinate of ZH point at end of row
            zh_point_y = -y0 + yi[i] + yj[j_max + 1]
            FacePoints_coords.append([zh_point_x, zh_point_y])
            
        for j in range(j_min, j_max + 2):
            zh_point_x = -x0 + xi[i_max + 1] + xj[j]                         # coordinate of ZH point (final row)
            zh_point_y = -y0 + yi[i_max + 1] + yj[j]
            FacePoints_coords.append([zh_point_x, zh_point_y])
        
        FacePoints_coords = np.array(FacePoints_coords)
        if out_folder is not None:
            pd.DataFrame(FacePoints_coords, columns = ['X', 'Y']).to_csv('{0}\{1}_FacePoint.csv'.format(out_folder, out_file_prefix), index_label='Index')
        
        # Add external 'cells' along perimeter
        x_shift = 0.5 * cell_size * math.cos(grid_rotation_angle)           # shift x, y coordinates along x^ direction
        y_shift = 0.5 * cell_size * math.sin(grid_rotation_angle)
        
        for j in range(j_min, j_max + 1):                                   # external 'cells' along row 0 perimeter
            zc_point_x = -x0 + xi[i_min] + xj[j] + x_shift
            zc_point_y = -y0 + yi[i_min] + yj[j] + y_shift
            CellCenter_coords.append([zc_point_x, zc_point_y])
        
        for j in range(j_min, j_max + 1):                                   # external 'cells' along row n* perimeter
            zc_point_x = -x0 + xi[i_max + 1] + xj[j] + x_shift
            zc_point_y = -y0 + yi[i_max + 1] + yj[j] + y_shift
            CellCenter_coords.append([zc_point_x, zc_point_y])
        
        x_shift = 0.5 * cell_size * math.cos(grid_rotation_angle + deg90_as_rads)       # shift x, y coordinates along y^ direction
        y_shift = 0.5 * cell_size * math.sin(grid_rotation_angle + deg90_as_rads)
        
        for i in range(i_min, i_max + 1):                                   # external 'cells' along column 0 perimeter
            zc_point_x = -x0 + xi[i] + xj[j_min] + x_shift
            zc_point_y = -y0 + yi[i] + yj[j_min] + y_shift
            CellCenter_coords.append([zc_point_x, zc_point_y])
        
        for i in range(i_min, i_max + 1):                                   # external 'cells' along column m* perimeter
            zc_point_x = -x0 + xi[i] + xj[j_max + 1] + x_shift
            zc_point_y = -y0 + yi[i] + yj[j_max + 1] + y_shift
            CellCenter_coords.append([zc_point_x, zc_point_y])
        
        CellCenter_coords = np.array(CellCenter_coords)
        if out_folder is not None:
            pd.DataFrame(CellCenter_coords, columns = ['X', 'Y']).to_csv('{0}\{1}_CellCentre.csv'.format(out_folder, out_file_prefix), index_label='Index')
        
        # Define parameters
        m_star = i_max - i_min + 1                                  # Number of rows and columns in reduced grid
        n_star = j_max - j_min + 1
        self.m_star = m_star
        self.n_star = n_star
        
        # number_of_facepoints = (n_star + 1) * (m_star + 1)
        # number_of_cells = (n_star + 2) * (m_star + 2) - 4           # includes external 'cells'
        number_of_faces = (2 * m_star + 1) * n_star + m_star
        
        index_of_external_cells = {'lower': m_star * n_star,
                                   'upper': m_star * (n_star + 1),
                                   'left' : m_star * (n_star + 2),
                                   'right': m_star * (n_star + 2) + n_star}
        
        # Extract perimeter points
        Perimeter = [FacePoints_coords[0, ], 
                     FacePoints_coords[n_star, ], 
                     FacePoints_coords[-1, ], 
                     FacePoints_coords[(n_star + 1) * m_star, ],
                     FacePoints_coords[0, ]]
        
        Perimeter = np.array(Perimeter)
        if out_folder is not None:
            pd.DataFrame(Perimeter, columns = ['X', 'Y']).to_csv('{0}\{1}_Perimeter.csv'.format(out_folder, out_file_prefix), index_label='Index')
        
        mesh_geometry['Cells Center Coordinate'] = CellCenter_coords
        mesh_geometry['FacePoints Coordinate'] = FacePoints_coords
        
        mesh_geometry['Perimeter'] = Perimeter
        mesh_geometry['Faces Perimeter Info'] = np.zeros(shape=(number_of_faces, 2))
        mesh_geometry['Faces Perimeter Values'] = np.empty(shape=(0, 0))
        
        # Face indices    
        Faces_FacePoint_index = []
        for i in range(n_star):
            for j in range (m_star): 
                facepoint_index = i * (m_star + 1) + j
                Faces_FacePoint_index.append([facepoint_index, facepoint_index + m_star + 1])       # U face
                Faces_FacePoint_index.append([facepoint_index, facepoint_index + 1])                # V face
            
            j = m_star                                                                              # last facepoint on each row, U face only
            facepoint_index = i * (m_star + 1) + j
            Faces_FacePoint_index.append([facepoint_index, facepoint_index + m_star + 1])
        
        i = n_star                              # last row, V face only
        for j in range(m_star):
            facepoint_index = i * (m_star + 1) + j
            Faces_FacePoint_index.append([facepoint_index, facepoint_index + 1])
            
        Faces_FacePoint_index = np.array(Faces_FacePoint_index)
        
        
        # Faces cell index
        Faces_Cell_index = []
        
        ## First row
        i = 0; j = 0
        cell_index = 0
        Faces_Cell_index.append([index_of_external_cells['left'], 0])
        Faces_Cell_index.append([index_of_external_cells['lower'], 0])
        
        i = 0
        for j in range(1, m_star):
            cell_index = j
            Faces_Cell_index.append([cell_index - 1, cell_index])
            Faces_Cell_index.append([index_of_external_cells['lower'] + j, j])
        
        i = 0; j = m_star - 1
        Faces_Cell_index.append([j, index_of_external_cells['right']])
        
        ## Subsequent rows
        for i in range(1, n_star):
            j = 0
            cell_index = i * m_star + j
            Faces_Cell_index.append([index_of_external_cells['left'] + i, cell_index])
            Faces_Cell_index.append([cell_index - m_star, cell_index])
            
            for j in range(1, m_star):
                cell_index = i * m_star + j
                Faces_Cell_index.append([cell_index - 1, cell_index])
                Faces_Cell_index.append([cell_index - m_star, cell_index])
            
            j = m_star - 1
            cell_index = i * m_star + j
            Faces_Cell_index.append([cell_index, index_of_external_cells['right'] + i])
        
        ## final row, V face only
        i = n_star - 1
        for j in range(m_star):
            cell_index = i * m_star + j
            Faces_Cell_index.append([cell_index, index_of_external_cells['upper'] + j])
        
        Faces_Cell_index = np.array(Faces_Cell_index)
        
        # for i in range(m_star - 1):
        #     # first cell in each row
            
            
        #     for j in range(n_star - 1):
        #         cell_index = i * n_star + j
        #         top_right_facepoint_index = (i + 1) * (n_star + 1) + j + 1
        #         Faces_Cell_index.append([cell_index, cell_index + 1])
        #         Faces_FacePoint_index.append([top_right_facepoint_index - (n_star + 1), top_right_facepoint_index])
        #         Faces_Cell_index.append([cell_index, cell_index + n_star])
        #         Faces_FacePoint_index.append([top_right_facepoint_index - 1, top_right_facepoint_index])
            
        #     # last cell in each row
        #     j = n_star - 1
        #     cell_index = i * n_star + j
        #     top_right_facepoint_index = (i + 1) * (n_star + 1) + j + 1
        #     Faces_Cell_index.append([cell_index, cell_index + n_star])
        #     Faces_FacePoint_index.append([top_right_facepoint_index - 1, top_right_facepoint_index])
        
        # ## last row of cells
        # i = m_star - 1
        # for j in range(n_star - 1):
        #     cell_index = i * n_star + j
        #     top_right_facepoint_index = (i + 1) * (n_star + 1) + j + 1
        #     Faces_Cell_index.append([cell_index, cell_index + 1])
        #     Faces_FacePoint_index.append([top_right_facepoint_index - (n_star + 1), top_right_facepoint_index])
        
        mesh_geometry['Faces Cell Indexes'] = Faces_Cell_index
        mesh_geometry['Faces FacePoint Indexes'] = Faces_FacePoint_index
        
        # Faces_FacePoint_index = []
        # i = 0
        # for j in range (1, n_star):
        #     facepoint_index = j
        #     Faces_FacePoint_index.append([facepoint_index, facepoint_index + n_star + 1])
        
        # for i in range(1, m_star):
        #     # j = 0
        #     facepoint_index = i * (n_star + 1) # + j
        #     Faces_FacePoint_index.append([facepoint_index, facepoint_index + 1])
            
        #     for j in range (1, n_star):
        #         facepoint_index = i * (n_star + 1) + j
        #         Faces_FacePoint_index.append([facepoint_index, facepoint_index + 1])
        #         Faces_FacePoint_index.append([facepoint_index, facepoint_index + n_star + 1])
                
        
        # Cells face orientation info and values
        number_internal_cells = m_star * n_star
        number_external_cells = 2 * (m_star + n_star)
        self.number_external_cells = number_external_cells
        counter_internal = np.arange(0, number_internal_cells * 4, 4).reshape(m_star * n_star, 1)        # Column [0, 4, 8, ... ]
        faces_internal = np.ones(shape = (number_internal_cells, 1), dtype = np.int32) * 4               # Column [4, 4, 4, ... ] representing 4 faces for each cell
        internal_cells = np.concatenate((counter_internal, faces_internal), axis = 1)
        
        counter_external = np.arange(number_external_cells).reshape(number_external_cells, 1) + (number_internal_cells * 4)   # 
        faces_external = np.ones(shape = (number_external_cells, 1), dtype = np.int32)                   # Column of ones, representing 1 face for each external 'cell'
        external_cells = np.concatenate((counter_external, faces_external), axis = 1)
        
        mesh_geometry['Cells Face and Orientation Info'] = np.concatenate((internal_cells, external_cells), axis = 0)
        
        internal_face_orientation = np.tile([1, 1, -1, -1], number_internal_cells).reshape((number_internal_cells * 4, 1))
        external_face_orientation = np.concatenate((np.repeat([-1, 1], m_star), np.repeat([-1, 1], n_star))).reshape(number_external_cells, 1)
        face_orientation = np.concatenate((internal_face_orientation, external_face_orientation), axis = 0)
        # face_orientation = -1 * face_orientation
        
        face_index = []
        for i in range(n_star - 1):
            for j in range(m_star):
                left_face = i * (2*m_star + 1) + 2*j
                lower_face = left_face + 1
                right_face = left_face + 2
                upper_face = lower_face + (2*m_star + 1)
                
                face_index.append(left_face)
                face_index.append(lower_face)
                face_index.append(right_face)
                face_index.append(upper_face)
        
        # last row
        i = n_star - 1
        for j in range(m_star):
            left_face = i * (2*m_star + 1) + 2*j
            lower_face = left_face + 1
            right_face = left_face + 2
            upper_face = n_star * (2*m_star + 1) + j
            
            face_index.append(left_face)
            face_index.append(lower_face)
            face_index.append(right_face)
            face_index.append(upper_face)
        
        # external "cells"
        for j in range(m_star):
            lower_ext_face = 2*j + 1
            face_index.append(lower_ext_face)
        
        i = n_star
        for j in range(m_star):
            upper_ext_face = i * (2*m_star + 1) + j
            face_index.append(upper_ext_face)
        
        for i in range(n_star):
            left_ext_face = i * (2*m_star + 1)
            face_index.append(left_ext_face)
        
        j = m_star
        for i in range(n_star):
            right_ext_face = i * (2*m_star + 1) + 2*j
            face_index.append(right_ext_face)
            
        face_index = np.array(face_index, dtype = np.int32).reshape((-1 , 1))
        
        # return face_index, face_orientation
        
        mesh_geometry['Cells Face and Orientation Values'] = np.concatenate((face_index, face_orientation), axis = 1)
        
        # Facepoints face orientation info and values
        face_index = []             # redundant – for clarity of code only
        face_orientation = []
        counter = []
        count = []
        
        ## first row
        i = 0; j = 0
        pos_x_face = 1
        face_index = [0, 1]
        face_orientation = [-1, -1]
        counter = [0]
        count = [2]
        
        for j in range(1, m_star):
            neg_x_face = pos_x_face
            pos_y_face = 2*j
            pos_x_face = pos_y_face + 1
            
            face_index.append(neg_x_face)
            face_index.append(pos_y_face)
            face_index.append(pos_x_face)
            
            face_orientation.extend([1, -1, -1])
            counter.append(sum(count))
            count.append(3)
            
            
        j = m_star
        neg_x_face = pos_x_face
        pos_y_face = 2*j
        face_index.append(neg_x_face)
        face_index.append(pos_y_face)
        
        face_orientation.extend([1, -1])
        counter.append(sum(count))
        count.append(2)
        
        ## Subsequent rows
        for i in range(1, n_star):
            j = 0
            neg_y_face = (i-1) * (2*m_star + 1)
            pos_y_face = i * (2*m_star + 1)
            pos_x_face = pos_y_face + 1
            
            face_index.append(neg_y_face)
            face_index.append(pos_y_face)
            face_index.append(pos_x_face)
            
            face_orientation.extend([1, -1, -1])
            counter.append(sum(count))
            count.append(3)
            
            for j in range(1, m_star):
                neg_y_face = (i-1) * (2*m_star + 1) + 2*j
                neg_x_face = pos_x_face
                pos_y_face = i * (2*m_star + 1) + 2*j
                pos_x_face = pos_y_face + 1
                
                face_index.append(neg_y_face)
                face_index.append(neg_x_face)
                face_index.append(pos_y_face)
                face_index.append(pos_x_face)
                
                face_orientation.extend([1, 1, -1, -1])
                counter.append(sum(count))
                count.append(4)
            
            j = m_star
            neg_y_face = (i-1) * (2*m_star + 1) + 2*j
            neg_x_face = pos_x_face
            pos_y_face = i * (2*m_star + 1) + 2*j
            
            face_index.append(neg_y_face)
            face_index.append(neg_x_face)
            face_index.append(pos_y_face)
            
            face_orientation.extend([1, 1, -1])
            counter.append(sum(count))
            count.append(3)
        
        ## last row
        i = n_star
        j = 0
        neg_y_face = (i-1) * (2*m_star + 1)
        pos_x_face = i * (2*m_star + 1)
        
        face_index.append(neg_y_face)
        face_index.append(pos_x_face)
        
        face_orientation.extend([1, -1])
        counter.append(sum(count))
        count.append(2)
        
        for j in range(1, m_star):
            neg_y_face = (i-1) * (2*m_star + 1) + 2*j
            neg_x_face = pos_x_face
            pos_x_face += 1
            
            face_index.append(neg_y_face)
            face_index.append(neg_x_face)
            face_index.append(pos_x_face)
            
            face_orientation.extend([1, 1, -1])
            counter.append(sum(count))
            count.append(3)
        
        j = m_star
        neg_y_face = (i-1) * (2*m_star + 1) + 2*j
        neg_x_face = pos_x_face
        
        face_index.append(neg_y_face)
        face_index.append(neg_x_face)
        
        face_orientation.extend([1, 1])
        counter.append(sum(count))
        count.append(2)
        
        mesh_geometry['FacePoints Face and Orientation Info'] = np.array([counter, count]).transpose()
        mesh_geometry['FacePoints Face and Orientation Values'] = np.array([face_index, face_orientation]).transpose()
        
        self.data = mesh_geometry
        return # mesh_geometry
    
    def zero_velocity(self, number_of_timesteps):
        n_star = self.n_star
        m_star = self.m_star
        number_of_faces = (2*m_star + 1) * n_star + m_star
        
        velocity = np.zeros(shape=(number_of_timesteps, number_of_faces))
        
        return velocity
        
    
    def map_velocity(magnitude_of_velocity, direction_of_velocity, n_star, m_star):
        Face_Velocity = np.zeros((n_star, m_star))
        timesteps = 0
        for t in timesteps:
            
            # First row
            i = 0; j = 0
            mag_velocity = magnitude_of_velocity(t, j, i)                   # velocity at cell i, j at time t
            if mag_velocity > 0:
                dir_velocity = math.radians(direction_of_velocity(t, j, i))
                u = math.cos(dir_velocity) * mag_velocity
                v = math.sin(dir_velocity) * mag_velocity
                
                if u > 0:
                    Face_Velocity[i * (2*n_star -1) + 2*j + 1] += u         # no else condition
                if v > 0:
                    Face_Velocity[i * (2*n_star -1) + 2*j + 2] += v         # no else condition
            
            for j in range(1, m_star - 1):
                mag_velocity = magnitude_of_velocity(t, j, i)                   # velocity at cell i, j at time t
                if mag_velocity > 0:
                    dir_velocity = math.radians(direction_of_velocity(t, j, i))
                    u = math.cos(dir_velocity) * mag_velocity
                    v = math.sin(dir_velocity) * mag_velocity
                    
                    if u > 0:
                        Face_Velocity[i * (2*n_star -1) + 2*j + 1] += u
                    else:
                        Face_Velocity[i * (2*n_star -1) + 2*j - 1] += u
                    if v > 0:
                        Face_Velocity[i * (2*n_star -1) + 2*j + 2] += v         # no else condition
            
            j = m_star - 1
            mag_velocity = magnitude_of_velocity(t, j, i)                   # velocity at cell i, j at time t
            if mag_velocity > 0:
                dir_velocity = math.radians(direction_of_velocity(t, j, i))
                u = math.cos(dir_velocity) * mag_velocity
                v = math.sin(dir_velocity) * mag_velocity
                
                if u < 0:  #LESS THAN
                    Face_Velocity[i * (2*n_star -1) + 2*j - 1] += u
                if v > 0:
                    Face_Velocity[i * (2*n_star -1) + 2*j + 2] += v         # no else condition
            
            # Subsequent rows
            for i in range(1, n_star - 1):
                j = 0
                mag_velocity = magnitude_of_velocity(t, j, i)                   # velocity at cell i, j at time t
                if mag_velocity > 0:
                    dir_velocity = math.radians(direction_of_velocity(t, j, i))
                    u = math.cos(dir_velocity) * mag_velocity
                    v = math.sin(dir_velocity) * mag_velocity
                    
                    if u > 0:
                        Face_Velocity[i * (2*n_star -1) + 2*j + 1] += u         # No else condition
                    if v > 0:
                        Face_Velocity[i * (2*n_star -1) + 2*j + 2] += v
                    else:
                        Face_Velocity[(i-1)* (2*n_star -1) + 2*j + 2] += v
                
                for j in range(1, m_star - 1):
                    mag_velocity = magnitude_of_velocity(t, j, i)                   # velocity at cell i, j at time t
                    if mag_velocity > 0:
                        dir_velocity = math.radians(direction_of_velocity(t, j, i))
                        u = math.cos(dir_velocity) * mag_velocity
                        v = math.sin(dir_velocity) * mag_velocity
                        
                        if u > 0:
                            Face_Velocity[i * (2*n_star -1) + 2*j + 1] += u
                        else:
                            Face_Velocity[i * (2*n_star -1) + 2*j - 1] += u
                        if v > 0:
                            Face_Velocity[i * (2*n_star -1) + 2*j + 2] += v
                        else:
                            Face_Velocity[(i-1)* (2*n_star -1) + 2*j + 2] += v
                
                j = m_star - 1
                mag_velocity = magnitude_of_velocity(t, j, i)                   # velocity at cell i, j at time t
                if mag_velocity > 0:
                    dir_velocity = math.radians(direction_of_velocity(t, j, i))
                    u = math.cos(dir_velocity) * mag_velocity
                    v = math.sin(dir_velocity) * mag_velocity
                    
                    if u < 0:  #LESS THAN
                        Face_Velocity[i * (2*n_star -1) + 2*j - 1] += u
                    if v > 0:
                        Face_Velocity[i * (2*n_star -1) + 2*j + 2] += v
                    else:
                        Face_Velocity[(i-1)* (2*n_star -1) + 2*j + 2] += v
            
            # last row
            i = n_star - 1
            j = 0
            mag_velocity = magnitude_of_velocity(t, j, i)                   # velocity at cell i, j at time t
            if mag_velocity > 0:
                dir_velocity = math.radians(direction_of_velocity(t, j, i))
                u = math.cos(dir_velocity) * mag_velocity
                v = math.sin(dir_velocity) * mag_velocity
                
                if u > 0:
                    Face_Velocity[i * (2*n_star -1) + 2*j + 1] += u         # No else condition
                if v > 0:
                    pass
                else:
                    Face_Velocity[(i-1)* (2*n_star -1) + 2*j + 2] += v
                    
            for j in range(1, m_star - 1):
                mag_velocity = magnitude_of_velocity(t, j, i)                   # velocity at cell i, j at time t
                if mag_velocity > 0:
                    dir_velocity = math.radians(direction_of_velocity(t, j, i))
                    u = math.cos(dir_velocity) * mag_velocity
                    v = math.sin(dir_velocity) * mag_velocity
                    
                    if u > 0:
                        Face_Velocity[i * (2*n_star -1) + 2*j + 1] += u
                    else:
                        Face_Velocity[i * (2*n_star -1) + 2*j - 1] += u
                    if v > 0:
                        pass
                    else:
                        Face_Velocity[(i-1)* (2*n_star -1) + 2*j + 2] += v
                        
            j = m_star - 1
            mag_velocity = magnitude_of_velocity(t, j, i)                   # velocity at cell i, j at time t
            if mag_velocity > 0:
                dir_velocity = math.radians(direction_of_velocity(t, j, i))
                u = math.cos(dir_velocity) * mag_velocity
                v = math.sin(dir_velocity) * mag_velocity
                
                if u > 0:
                    pass
                else:
                    Face_Velocity[i * (2*n_star -1) + 2*j - 1] += u
                if v > 0:
                    pass
                else:
                    Face_Velocity[(i-1)* (2*n_star -1) + 2*j + 2] += v
        
        return Face_Velocity
    