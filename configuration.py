# !/Users/maximejuventin/anaconda3/bin/python




import os
import sys
import getpass

import xarray as xr
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import datetime
on_disk = True
# work_dir = os.path.dirname(os.path.dirname(os.path.realpath(__file__))) + '/'
work_dir = '/Users/maximejuventin/Desktop/Analysis_miRNA/'
# print(work_dir)
if on_disk :
    scoring_data_dir = '/Volumes/TOSHIBA EXT/scripts_MirUp/scoring/'
    edf_data_dir = '/Volumes/TOSHIBA EXT/scripts_MirUp/MiRNA_edf/'

excel_dir = work_dir + '/excels/'
if getpass.getuser() in ('samuel.garcia', 'maxime.juventin') and  sys.platform.startswith('linux'):
    precompute_dir = '/mnt/data/CAmatPion201801_intra_respi_Maxime/amalo/precompute/'
else :
    precompute_dir = work_dir + '/precompute/'
    # precompute_dir = 'Z:/precompute_edf/'
print(precompute_dir)
RT_PCR_execption = ['B2878','B5850', 'B2879', 'B2877','B5921', 'B3720']

if __name__ == '__main__':
    print(work_dir)
