#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Sep 15 09:24:54 2019

@author: ckychang
"""

"""
Converts csv files downloaded from the Ameriflux archive into a CF-compliant netCDF4 file.
"""
import glob
import numpy as np
from netCDF4 import Dataset

# define the io function
def nc_io_function(filename, var_name, units):
    # list of files
    files = glob.glob(filename)
    # print(files)
    matching = str([s for s in files if "2100" in s])
    matching = matching.replace('[','').replace(']','').replace("'",'')
    idx = files.index(matching)
    for ii in range(0, idx+1):
        dataset = Dataset(files[ii])
        model_time = np.array(dataset.variables['time'])
        model_lon = np.array(dataset.variables['lon'])
        model_lat = np.array(dataset.variables['lat'])
        data = np.array(dataset.variables[var_name])
        data[data<=-999.0] = np.nan
        data[data>=10.0**20] = np.nan
        # aggregate the time series
        if (ii==0):
            time = model_time
            data_agg = data
        else:
            time = np.concatenate((time, model_time))
            data_agg = np.concatenate((data_agg, data))
    # Open a dataset for writing
    filename = filename.replace('*','')
    output = './ensemble_members/'+filename+'.nc'
    dset = Dataset(output,mode="w")
    # Create netCDF dimensions
    dset.createDimension("time",size=time.size)
    dset.createDimension("lat",size=model_lat.size)
    dset.createDimension("lon",size=model_lon.size)
    # Create netCDF variables
    T  = dset.createVariable("time"       ,time.dtype   ,("time"       ))
    X  = dset.createVariable("lat"        ,model_lat.dtype ,("lat"       ))
    Y  = dset.createVariable("lon"        ,model_lon.dtype ,("lon"       ))
    C_out  = dset.createVariable(var_name,data_agg.dtype,("time","lat","lon"))
    # Load data and encode attributes
    T [...]    = time
    T.units    = "days since 1850-01-01"
    T.calendar = "noleap"
    X [...]    = model_lat
    X.units    = "degrees_north"
    Y [...]    = model_lon
    Y.units    = "degrees_east"
    C_out[...]     = data_agg
    C_out.units    = units
    dset.close()
# declare parameters
mon = np.array([31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31])
mon_doy = np.cumsum(mon)
# calculate the start and end DOY for each month
mon_start_doy = np.zeros(12)
mon_end_doy = np.zeros(12)
for m in range(0,12):
    if m==0:
        mon_start_doy[m] = 1
        mon_end_doy[m] = mon_doy[m]
    else:
        mon_start_doy[m] = mon_doy[m-1]+1
        mon_end_doy[m] = mon_doy[m]

# start importing files
nc_io_function("evspsbl_Amon_E3SM-1-1-ECA_ssp585_r1i1p1f1_gr*", 'evspsbl', "kg m-2 s-1")
nc_io_function("evspsbl_Amon_EC-Earth3-CC_ssp585_r1i1p1f1_gr*", 'evspsbl', "kg m-2 s-1")
nc_io_function("evspsbl_Amon_EC-Earth3-Veg_ssp585_r1i1p1f1_gr*", 'evspsbl', "kg m-2 s-1")
nc_io_function("evspsbl_Amon_FGOALS-g3_ssp585_r1i1p1f1_gn*", 'evspsbl', "kg m-2 s-1")
nc_io_function("evspsbl_Amon_GISS-E2-1-G_ssp585_r1i1p3f1_gn*", 'evspsbl', "kg m-2 s-1")
nc_io_function("evspsbl_Amon_GISS-E2-1-H_ssp585_r2i1p1f2_gn*", 'evspsbl', "kg m-2 s-1")
nc_io_function("evspsbl_Amon_GISS-E2-2-G_ssp585_r3i1p3f1_gn*", 'evspsbl', "kg m-2 s-1")
nc_io_function("evspsbl_Amon_MPI-ESM1-2-LR_ssp585_r11i1p1f1_gn*", 'evspsbl', "kg m-2 s-1")
nc_io_function("evspsbl_Amon_UKESM1-0-LL_ssp585_r1i1p1f2_gn*", 'evspsbl', "kg m-2 s-1")
#  gpp
nc_io_function("gpp_Lmon_EC-Earth3-CC_ssp585_r1i1p1f1_gr*", 'gpp', "kg m-2 s-1")
nc_io_function("gpp_Lmon_GISS-E2-1-G_ssp585_r1i1p3f1_gn*", 'gpp', "kg m-2 s-1")
nc_io_function("gpp_Lmon_GISS-E2-1-H_ssp585_r1i1p1f2_gn*", 'gpp', "kg m-2 s-1")
nc_io_function("gpp_Lmon_GISS-E2-2-G_ssp585_r1i1p3f1_gn*", 'gpp', "kg m-2 s-1")
nc_io_function("gpp_Lmon_MPI-ESM1-2-LR_ssp585_r11i1p1f1_gn*", 'gpp', "kg m-2 s-1")
nc_io_function("gpp_Lmon_UKESM1-0-LL_ssp585_r1i1p1f2_gn*", 'gpp', "kg m-2 s-1")
#  hurs
nc_io_function("hurs_Amon_AWI-CM-1-1-MR_ssp585_r1i1p1f1_gn*", 'hurs', "%")
nc_io_function("hurs_Amon_EC-Earth3-CC_ssp585_r1i1p1f1_gr*", 'hurs', "%")
nc_io_function("hurs_Amon_FGOALS-g3_ssp585_r1i1p1f1_gn*", 'hurs', "%")
nc_io_function("hurs_Amon_GISS-E2-1-G_ssp585_r1i1p3f1_gn*", 'hurs', "%")
nc_io_function("hurs_Amon_GISS-E2-1-H_ssp585_r2i1p1f2_gn*", 'hurs', "%")
nc_io_function("hurs_Amon_GISS-E2-2-G_ssp585_r3i1p3f1_gn*", 'hurs', "%")
nc_io_function("hurs_Amon_MPI-ESM1-2-HR_ssp585_r2i1p1f1_gn*", 'hurs', "%")
nc_io_function("hurs_Amon_MPI-ESM1-2-LR_ssp585_r11i1p1f1_gn_*", 'hurs', "%")
nc_io_function("hurs_Amon_UKESM1-0-LL_ssp585_r1i1p1f2_gn*", 'hurs', "%")
#  npp
nc_io_function("npp_Lmon_EC-Earth3-CC_ssp585_r1i1p1f1_gr*", 'npp', "kg m-2 s-1")
nc_io_function("npp_Lmon_GISS-E2-1-G_ssp585_r1i1p3f1_gn*", 'npp', "kg m-2 s-1")
nc_io_function("npp_Lmon_GISS-E2-1-H_ssp585_r3i1p1f2_gn*", 'npp', "kg m-2 s-1")
nc_io_function("npp_Lmon_GISS-E2-2-G_ssp585_r1i1p3f1_gn*", 'npp', "kg m-2 s-1")
nc_io_function("npp_Lmon_MPI-ESM1-2-LR_ssp585_r11i1p1f1_gn*", 'npp', "kg m-2 s-1")
nc_io_function("npp_Lmon_UKESM1-0-LL_ssp585_r1i1p1f2_gn*", 'npp', "kg m-2 s-1")
#  tas
nc_io_function("tas_Amon_E3SM-1-1-ECA_ssp585_r1i1p1f1_gr*", 'tas', "K")
nc_io_function("tas_Amon_EC-Earth3-CC_ssp585_r1i1p1f1_gr*", 'tas', "K")
nc_io_function("tas_Amon_FGOALS-g3_ssp585_r1i1p1f1_gn*", 'tas', "K")
nc_io_function("tas_Amon_GISS-E2-1-G_ssp585_r1i1p3f1_gn*", 'tas', "K")
nc_io_function("tas_Amon_GISS-E2-1-H_ssp585_r2i1p1f2_gn*", 'tas', "K")
nc_io_function("tas_Amon_GISS-E2-2-G_ssp585_r3i1p3f1_gn*", 'tas', "K")
nc_io_function("tas_Amon_MPI-ESM1-2-LR_ssp585_r11i1p1f1_gn*", 'tas', "K")
nc_io_function("tas_Amon_UKESM1-0-LL_ssp585_r1i1p1f2_gn*", 'tas', "K")
# tran
nc_io_function("tran_Lmon_E3SM-1-1-ECA_ssp585_r1i1p1f1_gr*", 'tran', "kg m-2 s-1")
nc_io_function("tran_Lmon_EC-Earth3-CC_ssp585_r1i1p1f1_gr*", 'tran', "kg m-2 s-1")
nc_io_function("tran_Lmon_GISS-E2-1-G_ssp585_r1i1p3f1_gn*", 'tran', "kg m-2 s-1")
nc_io_function("tran_Lmon_GISS-E2-1-H_ssp585_r1i1p1f2_gn*", 'tran', "kg m-2 s-1")
nc_io_function("tran_Lmon_GISS-E2-2-G_ssp585_r1i1p3f1_gn*", 'tran', "kg m-2 s-1")
nc_io_function("tran_Lmon_HadGEM3-GC31-LL_ssp585_r2i1p1f3_gn*", 'tran', "kg m-2 s-1")
nc_io_function("tran_Lmon_MPI-ESM1-2-LR_ssp585_r11i1p1f1_gn*", 'tran', "kg m-2 s-1")
nc_io_function("tran_Lmon_UKESM1-0-LL_ssp585_r1i1p1f2_gn*", 'tran', "kg m-2 s-1")

print("finish encoding")

