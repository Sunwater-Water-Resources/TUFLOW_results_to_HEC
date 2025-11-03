# -*- coding: utf-8 -*-
"""
Created on Thu May  2 15:00:25 2024

@author: PanosotG
"""

import h5py

file1 = r'C:\PythonProjects\tuflow_to_lifesim\Output\cBaldEagleDamBrk_55.hdf'
file2 = r'C:\PythonProjects\tuflow_to_lifesim\Output\KolanRiver_BCLS_HR02_HYD002_HGS_005_NC_CC11.hdf'

with h5py.File(file2) as hf:
    hf.visit(print)