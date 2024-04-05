#-------------------------------------------------------------------------------
# Name:        getClimGenFns.py
# Purpose:     additional functions for getClimGenNC.py
# Author:      s03mm5
# Created:     08/02/2018
# Copyright:   (c) s03mm5 2015
# Licence:     <your licence>
#-------------------------------------------------------------------------------
#!/usr/bin/env python

__prog__ = 'getClimGenFns.py'
__author__ = 's03mm5'

import time
import sys
from glob import glob
import netCDF4 as cdf
import numpy as np
import math

NGRANULARITY = 120
DELTA_DEG = 0.001   # ensures creation of enclosing climate grid

def associate_climate(site_rec, climgen, pettmp_hist, pettmp_fut):
    """
    this function associates each soil grid point with the proximate upper right gridded climate data
    at the time of writing (Dec 2015) soil data is on 30 second grid whereas climate data is half or eighth of a degree.
    """
    func_name =  __prog__ + ' associate_climate'

    resol_lat = climgen.resolution_lat
    resol_lon = climgen.resolution_lon

    proximate_keys = {}
    gran_lat, gran_lon, latitude, longitude, dummy, dummy = site_rec

    # determine nearest climate location - this code is repeated in function genLocalGrid, getClimGenNC.py
    # ==================================
    num_lats = math.ceil( abs((latitude - climgen.lat_min)/resol_lat) )
    latMax = abs(num_lats*resol_lat) + climgen.lat_min
    gran_lat_ur = int( round((90.0 - latMax)*NGRANULARITY))

    # if soil grid cell latitude coincides then shift latitude slightly southwards to complete enclosure
    if gran_lat_ur == gran_lat:
        latitude = latitude - DELTA_DEG

    num_lats = int( abs((latitude - climgen.lat_min)/resol_lat) )
    latMin = abs(num_lats*resol_lat) + climgen.lat_min
    gran_lat_ll = int( round((90.0 - latMin)*NGRANULARITY))

    num_lons = math.ceil((longitude - climgen.lon_min)/resol_lon)
    lonMax = num_lons*resol_lon + climgen.lon_min
    gran_lon_ur = int( round((180.0 + lonMax)*NGRANULARITY))

    # if soil grid cell longitude coincides then shift longitude slightly westwards to complete enclosure
    if gran_lon_ur == gran_lon:
        longitude = longitude - DELTA_DEG

    num_lons = int((longitude - climgen.lon_min)/resol_lon)
    lonMin = num_lons*resol_lon + climgen.lon_min
    gran_lon_ll = int( round((180.0 + lonMin)*NGRANULARITY))

    # collect other proximate climate locations
    # =========================================
    pettmp_fut_precip = pettmp_fut['precipitation']
    pettmp_fut_temp   = pettmp_fut['temperature']
    key_ur = '{:0=5d}_{:0=5d}'.format(gran_lat_ur, gran_lon_ur)
    if key_ur in pettmp_fut_precip and key_ur in pettmp_fut_temp:
        if len(pettmp_fut_precip[key_ur]) > 0 and len(pettmp_fut_temp[key_ur]) > 0:
            proximate_keys['key_ur'] = list([key_ur, gran_lat_ur, gran_lon_ur])

    key_ll = '{:0=5d}_{:0=5d}'.format(gran_lat_ll, gran_lon_ll)
    if key_ll in pettmp_fut_precip and key_ll in pettmp_fut_temp:
        if len(pettmp_fut_precip[key_ll]) > 0 and len(pettmp_fut_temp[key_ll]) > 0:
            proximate_keys['key_ll'] = list([key_ll, gran_lat_ll, gran_lon_ll])

    key_ul = '{:0=5d}_{:0=5d}'.format(gran_lat_ur, gran_lon_ll)
    if key_ul in pettmp_fut_precip and key_ul in pettmp_fut_temp:
        if len(pettmp_fut_precip[key_ul]) > 0 and len(pettmp_fut_temp[key_ul]) > 0:
            proximate_keys['key_ul'] = list([key_ul, gran_lat_ur, gran_lon_ll])

    key_lr = '{:0=5d}_{:0=5d}'.format(gran_lat_ll, gran_lon_ur)
    if key_lr in pettmp_fut_precip and key_lr in pettmp_fut_temp:
        if len(pettmp_fut_precip[key_lr]) > 0 and len(pettmp_fut_temp[key_lr]) > 0:
            proximate_keys['key_lr'] = list([key_lr, gran_lat_ll, gran_lon_ur])

    if len(proximate_keys) == 0:
        mess = 'No weather keys assigned for site record with granular coordinates: '
        mess += '{} {}\tand lat/lon: {} {}'.format(gran_lat, gran_lon, round(latitude,4), round(longitude,4))
        climgen.lgr.info(mess)
        # print(mess)
        ret_code = None
    else:
        ret_code = _apply_weather(climgen, proximate_keys, gran_lat, gran_lon, pettmp_hist, pettmp_fut)

    return ret_code

def _apply_weather(climgen, proximate_keys, gran_lat_cell, gran_lon_cell, pettmp_hist, pettmp_fut, use_min_val = True):
    """
    use the inverse distance squared downscaling method to determine weather for specified grid cell
    """

    # first stanza: calculate the squares of the distances in granular units between the grid cell and weathers cells
    # =============
    dist = {}
    metric_list = pettmp_fut.keys()
    for key in proximate_keys:
        lookup_key, gran_lat, gran_lon = proximate_keys[key]

        # situation where grid cell is coincidental with weather cell
        # ===========================================================
        if gran_lat == gran_lat_cell and gran_lon == gran_lon_cell:
            dist[lookup_key] = 0.0
        else:
            dist[lookup_key] = (gran_lat - gran_lat_cell)**2 + (gran_lon - gran_lon_cell)**2

    # find key corresponding to the minimum value using conversions to lists
    # ======================================================================
    minval = sorted(dist.values())[0]
    lookup_key = list(dist.keys())[list(dist.values()).index(minval)]
    pettmp_fut_final = {}
    pettmp_hist_final = {}
    for metric in metric_list:
        pettmp_fut_final[metric]  = pettmp_fut[metric][lookup_key]
        pettmp_hist_final[metric] = pettmp_hist[metric][lookup_key]

    return pettmp_hist_final, pettmp_fut_final, lookup_key

def check_clim_nc_limits(form, wthr_rsrce, bbox_aoi = None):
    """
    this function checks that the specified bounding box lies within extent of the requested weather dataset
    """
    func_name =  __prog__ + ' check_clim_nc_limits'

    limits_ok_flag = True
    if wthr_rsrce == 'NASA':
        return limits_ok_flag

    lon_ll_aoi = float(form.w_ll_lon.text())
    lat_ll_aoi = float(form.w_ll_lat.text())
    lon_ur_aoi = float(form.w_ur_lon.text())
    lat_ur_aoi = float(form.w_ur_lat.text())

    wthr_rsrce = wthr_rsrce + '_Mnth'
    lat_ur_dset = form.weather_sets[wthr_rsrce]['lat_ur']
    lon_ur_dset = form.weather_sets[wthr_rsrce]['lon_ur']
    lat_ll_dset = form.weather_sets[wthr_rsrce]['lat_ll']
    lon_ll_dset = form.weather_sets[wthr_rsrce]['lon_ll']

    # similar functionality in lu_extract_fns.py in LU_extract project
    # ================================================================
    if (lon_ll_dset < lon_ll_aoi and lon_ur_dset > lon_ur_aoi) and \
                    (lat_ll_dset < lat_ll_aoi and lat_ur_dset > lat_ur_aoi):
        print('AOI lies within ' + wthr_rsrce + ' weather datasets')
    else:
        print('AOI lies outwith ' + wthr_rsrce + ' weather datasets - LL long/lat: {} {}\tUR long/lat: {} {}'
              .format(lon_ll_dset, lat_ll_dset, lon_ur_dset, lat_ur_dset))
        limits_ok_flag = False

    return limits_ok_flag

def update_progress_clim_soil(last_time, nsoilres, pt_key, ncsv_lines, skipped = 0, failed = 0):

    """Update progress bar."""
    new_time = time.time()
    if new_time - last_time > 5:

        mess = '\rSize of soil list: {}\tpt_key: {}\tNumber of sites remaining: {}'\
                                            .format(nsoilres, pt_key, ncsv_lines - nsoilres)
        sys.stdout.flush()
        sys.stdout.write(mess)
        last_time = new_time

    return last_time

def update_progress_clim(last_time, ngrid_cells, total_num, no_data):
    """Update progress bar."""
    new_time = time.time()
    if new_time - last_time > 2:
        mess = '\rCompleted checking of: {} climate cells\tNo data: {}\tRemaining: {}'\
                                            .format(ngrid_cells, no_data, 1 + total_num - ngrid_cells)
        sys.stdout.flush()
        sys.stdout.write(mess)
        last_time = new_time

    return last_time
