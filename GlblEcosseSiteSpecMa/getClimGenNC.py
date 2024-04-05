#-------------------------------------------------------------------------------
# Name:        getClimGenNC.py
# Purpose:     read netCDF files comprising ClimGen data
# Author:      s03mm5
# Created:     08/12/2015
# Copyright:   (c) s03mm5 2015
# Licence:     <your licence>
#-------------------------------------------------------------------------------
#!/usr/bin/env python

__prog__ = 'getClimGenNC.py'
__author__ = 's03mm5'

from os import remove, mkdir
from os.path import join, normpath, isdir, isfile, split
from glob import glob
from copy import copy
from csv import writer as csv_writer

from _datetime import datetime
from netCDF4 import Dataset, date2index
import math

from numpy import seterr, ma, arange
from numpy.ma.core import MaskedConstant
from warnings import simplefilter

from thornthwaite import thornthwaite

SPACER_LEN = 12

_MONTHDAYS = [31,28,31,30,31,30,31,31,30,31,30,31]
LEAP_MONTHDAYS = [31,29,31,30,31,30,31,31,30,31,30,31]

numSecsDay = 3600*24
NGRANULARITY = 120
METRICS = list(['precipitation', 'temperature'])
ERROR_STR = '*** Error *** '

def _consistency_check(pettmp, varnams_mapped):
    """
    make sure for a each key if one metric is zero length then all other metrics for that key are also blank
    TODO: this function only works for two metrics and is unpythonic!
    """
    metric_list = list(varnams_mapped.values())
    metric0 = metric_list[0]
    metric1 = metric_list[1]
    for key in pettmp[metric0]:
        len_key0 = len(pettmp[metric0][key])

        if len_key0 == 0:
            pettmp[metric1][key] = []

        len_key1 = len(pettmp[metric1][key])
        if len_key1 == 0:
            pettmp[metric0][key] = []

    return pettmp

def _check_list_for_none(metric_list):
    """
    if a None is found then return an empty list
    """
    for indx, val in enumerate(metric_list):
        if val is None:
            return []

    return metric_list

def _input_txt_line_layout(data, comment):

        spacer_len = max(SPACER_LEN - len(data), 2)
        spacer = ' ' * spacer_len
        return '{}{}# {}\n'.format(data, spacer, comment)

class ClimGenNC(object,):

    def __init__(self, form, start_from_1801):
        """
        typically form.inpnc_dir = r'E:\mark2mike\climgenNC'  (get climgen future climate netCDF4 data from here)
                 form.inp_hist_dir = r'E:\mark2mike\fut_data'  (get CRU historic climate netCDF4 data from here)
        """
        func_name =  __prog__ +  ' ClimGenNC __init__'

        # determine user choices
        # ======================
        wthr_rsrce = form.combo10w.currentText()
        self.wthr_rsrce = wthr_rsrce
        wthr_rsrce_mnth = wthr_rsrce + '_Mnth'
        wthr_rsrce_day  = wthr_rsrce + '_Day'

        if form.w_mnthly.isChecked():
            self.fut_wthr_mnthly_flag = True
        else:
            self.fut_wthr_mnthly_flag = False

        # assign historic and future weather sets
        # =======================================
        hist_start_year = int(form.combo09s.currentText())
        hist_end_year = int(form.combo09e.currentText())
        lat = 'latitude'
        lon = 'longitude'
        hist_weather_set = form.weather_sets[wthr_rsrce_mnth]

        if self.fut_wthr_mnthly_flag:
            fut_weather_set  = form.weather_sets[wthr_rsrce_mnth]
        else:
            fut_weather_set  = form.weather_sets[wthr_rsrce_day]

        fut_clim_scen = '26'      # use IPCC scenarios

        # create weather resource directory if necessary
        # ==============================================
        rgn_wthr_dir = wthr_rsrce
        clim_dir = normpath(join(form.sttngs['sims_dir'], rgn_wthr_dir))
        if not isdir(clim_dir):
            mkdir(clim_dir)
            print('\tcreated: ' + clim_dir)

        self.rgn_wthr_dir = rgn_wthr_dir

        hist_start_year = max(hist_weather_set['start_year'], hist_start_year)
        hist_end_year   = min(hist_weather_set['end_year'], hist_end_year)

        num_hist_years = hist_end_year - hist_start_year + 1
        self.num_hist_years = num_hist_years
        self.hist_start_year = hist_start_year
        self.hist_end_year   = hist_end_year
        self.months = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec']
        self.fut_clim_scen = fut_clim_scen

        # future data
        # ===========
        self.fut_precip_fname = fut_weather_set['ds_precip']
        self.fut_tas_fname = fut_weather_set['ds_tas']
        self.resolution_lon = fut_weather_set['resol_lat']
        self.resolution_lat = fut_weather_set['resol_lon']

        # granularity
        # ===========
        self.lon = lon
        self.lat = lat
        self.lon_min = fut_weather_set['lat_ll']
        self.lat_min = fut_weather_set['lon_ll']
        self.longitudes = fut_weather_set['longitudes']
        self.latitudes =  fut_weather_set['latitudes']
        self.pettmp = {}        # dictionary whose keys will reference the climate grid, pt_grid
        self.lgr = form.lgr
        self.fobjs = form.fobjs
        self.zeros_file = form.zeros_file

        # past (monthly) data
        # ===================
        self.hist_precip_fname = hist_weather_set['ds_precip']
        self.hist_tas_fname    = hist_weather_set['ds_tas']
        self.latitudes_hist    = hist_weather_set['latitudes']

        # New stanza to facilitate option when user selects "use average weather"
        # =======================================================================
        # num_years_str = '{:0=3d}'.format(num_hist_years)
        met_ave_file = 'met' + str(hist_start_year) + '_' + str(hist_end_year) + 'a.txt'
        self.met_ave_file = met_ave_file

        # only the years for which we we have historic data will be taken into account
        self.num_ave_wthr_years = hist_end_year - hist_start_year + 1

        # simulation period
        # =================
        self.start_from_1801 = start_from_1801
        if form.years_from_flag == 'sowing':
            fut_start_year = int(form.crop_type_sowing['frst_year'])
            fut_end_year   = int(form.crop_type_sowing['last_year'])
        else:
            fut_start_year = int(form.combo11s.currentText())
            fut_end_year = int(form.combo11e.currentText())

        self.fut_start_year = fut_start_year
        self.fut_end_year   = fut_end_year
        self.num_fut_years = fut_end_year - fut_start_year + 1
        fut_ave_file = 'met{}_to_{}_ave.txt'.format(fut_start_year, fut_end_year)
        self.fut_ave_file = fut_ave_file

        req_met_fnames = ['met{}s.txt'.format(year) for year in range(fut_start_year, fut_end_year + 1)]
        if start_from_1801:
            add_met_fnames = [met_ave_file for year in range(1801, fut_start_year)]
            self.req_met_fnames = add_met_fnames + req_met_fnames
        else:
            self.req_met_fnames = req_met_fnames

        self.nsim_yrs = len(self.req_met_fnames)

    def genLocalGrid(self, bbox, snglPntFlag = False):
        """
        # return the weather indices for the area which encloses the supplied bounding box
        # this function does not alter the ClimGenNC (self) object
        """
        func_name =  __prog__ +  ' genLocalGrid'
        junk = seterr(all='ignore') # switch off warning messages

        if snglPntFlag:
            bbLonMin, bbLatMin = bbox
            bbLonMax = bbLonMin
            bbLatMax = bbLatMin
        else:
            bbLonMin, bbLatMin, bbLonMax, bbLatMax = bbox

        # determine bounds for climate grid which will enclose the supplied bounding box
        # ==============================================================================
        resol_lat = self.resolution_lat   # negative for future CRU data
        lat_indices = []
        lat_indices_hist = []
        clim_lat_min = self.lat_min
        num_lats = math.ceil( abs((bbLatMax - clim_lat_min)/resol_lat) )
        latMax = round(abs(num_lats*resol_lat) + clim_lat_min, 4)
        lat_indices.append(self.latitudes.index(latMax))
        lat_indices_hist.append(self.latitudes_hist.index(latMax))

        num_lats = int( abs((bbLatMin - clim_lat_min)/resol_lat) )
        latMin = round(abs(num_lats*resol_lat) + clim_lat_min, 4)
        try:
            lat_indices.append(self.latitudes.index(latMin))
        except ValueError as err:
            print(ERROR_STR + 'Latitude: ' + str(err) + ' in ' + func_name)
            return None

        lat_indices_hist.append(self.latitudes_hist.index(latMin))

        # longitudes
        # ==========
        lon_indices = []
        resol_lon = self.resolution_lon
        clim_lon_min = self.lon_min
        num_lons = math.ceil((bbLonMax - clim_lon_min)/resol_lon)
        lonMax = round(num_lons*resol_lon + clim_lon_min, 4)
        try:
            lon_indices.append(self.longitudes.index(lonMax))
        except ValueError as err:
            print(ERROR_STR + 'Longitude: ' + str(err) + ' in ' + func_name)
            return None

        num_lons = int((bbLonMin - clim_lon_min)/resol_lon)
        lonMin = round(num_lons*resol_lon + clim_lon_min, 4)
        lon_indices.append(self.longitudes.index(lonMin))

        # generate ClimGen grid    NB need to add one when slicing!!!
        # =====================    ==================================
        alons = arange(lonMin, lonMax, abs(resol_lon))
        alats = arange(latMin, latMax, abs(resol_lat))
        nlats = len(alats)
        nlons = len(alons)

        granlons = arange(nlons)
        for ic, lon in enumerate(alons):
            granlons[ic] = (180.0 + lon)*NGRANULARITY
        granlons.sort()

        granlats = arange(nlats)
        for ic, lat in enumerate(alats):
            granlats[ic] = (90.0 - lat)*NGRANULARITY
        granlats.sort()

        # must be in correct order
        # ========================
        lat_indices.sort()
        lat_indices_hist.sort()
        lon_indices.sort()

        aoi_indices_fut = lat_indices + lon_indices
        aoi_indices_hist = lat_indices_hist + lon_indices
        return aoi_indices_fut, aoi_indices_hist

    def check_ecss_wthr_data(self, sims_dir, climgen, pettmp_hist, num_band):
        """
        check to see if ECOSSE weather data already exists
        """
        clim_dir = normpath(join(sims_dir, climgen.rgn_wthr_dir))
        metric = METRICS[0]

        # check that number and names of met files match expectation
        # ==========================================================
        ncoords = len(pettmp_hist[metric])
        pettmp = {}
        n_met_data = 0
        for gran_coord in pettmp_hist[metric]:
            pettmp[gran_coord] = [False]
            gran_dir = join(clim_dir, gran_coord)
            if isdir(gran_dir):
                met_fnames = [split(met_fname)[1] for met_fname in glob(join(gran_dir, gran_coord) + '/met*s.txt')]
                if len(met_fnames) == climgen.nsim_yrs:
                    if met_fnames == climgen.req_met_fnames:
                        pettmp[gran_coord] = [True]
                        n_met_data += 1

        if n_met_data == ncoords:
            mess = 'Number of weather directories: '
            mess += '{}\tnumber with met data: {}\tfor band {} '.format(ncoords, n_met_data, num_band)
            print(mess)
            climgen.lgr.info(mess)

        pettmp_dict = {'precipitation': pettmp, 'temperature': pettmp}
        return pettmp_dict

    def fetch_ewembi_NC_data(self, aoi_indices, num_band, future_flag = True):
        """
        get precipitation or temperature data for a given variable and lat/long index for all times
        CRU uses NETCDF4 format
        """
        func_name = __prog__ +  ' fetch_fut_future_NC_data'
        simplefilter('default')

        num_key_masked = 0
        lat_indx_min, lat_indx_max, lon_indx_min, lon_indx_max = aoi_indices
        pettmp = {}

        # process future climate
        # ======================
        wthr_rsrce =  self.weather_resource
        varnams_mapped = {'pr':'precipitation','tas':'temperature'}
        varnams = sorted(varnams_mapped.keys())

        for varname, fname in zip(varnams, list([self.fut_precip_fname, self.fut_tas_fname])):
            varnam_map = varnams_mapped[varname]
            pettmp[varnam_map] = {}
            nc_dset = Dataset(fname, mode='r')

            # collect readings for all time values
            # ====================================
            slice = nc_dset.variables[varname][:, lat_indx_min:lat_indx_max + 1, lon_indx_min:lon_indx_max + 1]

            if ma.is_masked(slice):
                slice_is_masked_flag = True
                self.lgr.info('Future slice is masked in band {}'.format(num_band))
            else:
                slice_is_masked_flag = False

            # reform slice
            # ============
            for ilat, lat_indx in enumerate(range(lat_indx_min, lat_indx_max + 1)):
                gran_lat = round((90.0 - self.latitudes[lat_indx])*NGRANULARITY)

                for ilon, lon_indx in enumerate(range(lon_indx_min, lon_indx_max + 1)):
                    gran_lon = round((180.0 + self.longitudes[lon_indx])*NGRANULARITY)
                    key = '{:0=5d}_{:0=5d}'.format(int(gran_lat), int(gran_lon))

                    # validate values
                    # ===============
                    pettmp[varnam_map][key] = None
                    if slice_is_masked_flag :
                        val = slice[ilat,ilon,0]
                        if val is ma.masked:
                            self.lgr.info('val is ma.masked for key ' + key)
                            pettmp[varnam_map][key] = None
                            num_key_masked += 1

                    # add data for this coordinate
                    # ============================
                    if pettmp[varnam_map][key] is None:
                        if varname == 'pr':
                            pettmp[varnam_map][key] = [round(val*numSecsDay, 3) for val in slice[:,ilat,ilon]]
                        elif varname == 'tas':
                            pettmp[varnam_map][key] = [round(val - 273.15, 3) for val in slice[:,ilat,ilon]]

            # close netCDF file
            nc_dset.close()
            print('# masked weather keys: {}'.format(num_key_masked))

        return pettmp

    def fetch_eobs_NC_data(self, wthr_set_defn, aoi_indices, num_band, future_flag = True):
        """
        get precipitation or temperature data for a given variable and lat/long index for all times
        EObs uses NETCDF format
        historic data is monthly therefore no need to get time indices
        """
        func_name = __prog__ +  ' fetch_eobs_NC_data'
        simplefilter('default')

        lat_indx_min, lat_indx_max, lon_indx_min, lon_indx_max = aoi_indices
        pettmp = {}
        if future_flag == True:
            precip_fname = self.fut_precip_fname
            tas_fname    = self.fut_tas_fname
            strt_date = datetime(int(self.fut_start_year), 1, 1)    # year, month, dat, hour, minute, second
            end_date = datetime(int(self.fut_end_year), 12, 31)
        else:
            precip_fname = self.hist_precip_fname
            tas_fname    = self.hist_tas_fname

        # process future climate
        # ======================
        varnams_mapped = {'rr':'precipitation', 'tg':'temperature'}
        varnams = sorted(varnams_mapped.keys())

        for varname, fname in zip(varnams, list([precip_fname, tas_fname])):
            varnam_map = varnams_mapped[varname]
            pettmp[varnam_map] = {}
            nc_dset = Dataset(fname, mode='r')
            if wthr_set_defn['resol_time'] == 'Daily':
                strt_indx = date2index(strt_date, nc_dset.variables['time'])
                end_indx  = date2index(end_date,  nc_dset.variables['time'])

            # collect readings for all time values
            # ====================================
            try:
                if wthr_set_defn['resol_time'] == 'Daily':
                    slice = nc_dset.variables[varname][strt_indx:end_indx + 1, lat_indx_min:lat_indx_max + 1,
                                                                                    lon_indx_min:lon_indx_max + 1]
                else:
                    slice = nc_dset.variables[varname][:, lat_indx_min:lat_indx_max + 1,
                                                                                    lon_indx_min:lon_indx_max + 1]
            except RuntimeWarning as err:
                print(err)

            # reform slice
            # ============
            for ilat, lat_indx in enumerate(range(lat_indx_min, lat_indx_max + 1)):
                latitude = self.latitudes[lat_indx]
                gran_lat = round((90.0 - latitude)*NGRANULARITY)

                for ilon, lon_indx in enumerate(range(lon_indx_min, lon_indx_max + 1)):
                    longitude = self.longitudes[lon_indx]
                    gran_lon = round((180.0 + longitude)*NGRANULARITY)
                    key = '{:0=5d}_{:0=5d}'.format(int(gran_lat), int(gran_lon))

                    # validate values
                    # ===============
                    pettmp[varnam_map][key] = None
                    val = slice[0,ilat,ilon]
                    if type(val) is MaskedConstant:      # implies all values are masked constants with val.item() = 0.0
                        pettmp[varnam_map][key] = [] # implies all values are masked constants with val.item() = 0.0
                    else:
                        # add data for this coordinate
                        # ============================
                        if pettmp[varnam_map][key] is None:
                            if varname == 'rr':
                                rainfall = slice[:,ilat,ilon]
                                pettmp[varnam_map][key] = _check_list_for_none(rainfall.tolist())
                            elif varname == 'tg':
                                temperature = slice[:,ilat,ilon]
                                pettmp[varnam_map][key] = _check_list_for_none(temperature.tolist())

            nc_dset.close() # close netCDF file

        pettmp = _consistency_check(pettmp, varnams_mapped)
        return pettmp

    def fetch_nasa_NC_data(self, aoi_indices, num_band, future_flag = True):
        """
        get precipitation or temperature data for a given variable and lat/long index for all times
        """
        func_name = __prog__ +  ' fetch_nasa_NC_data'
        simplefilter('default')

        lat_indx_min, lat_indx_max, lon_indx_min, lon_indx_max = aoi_indices
        pettmp = {}
        if future_flag == True:
            precip_fname = self.fut_precip_fname
            tas_fname = self.fut_tas_fname
        else:
            precip_fname = self.hist_precip_fname
            tas_fname    = self.hist_tas_fname

        # process future climate
        # ======================
        varnams_mapped = {'prate':'precipitation','tmax':'temperature'}
        varnams = sorted(varnams_mapped.keys())

        for varname, fname in zip(varnams, list([precip_fname, tas_fname])):
            varnam_map = varnams_mapped[varname]
            pettmp[varnam_map] = {}
            nc_dset = Dataset(fname, mode='r')

            # collect readings for all time values
            # ====================================
            try:
                slice = nc_dset.variables[varname][:, lat_indx_min:lat_indx_max + 1, lon_indx_min:lon_indx_max + 1]
            except RuntimeWarning as err:
                pettmp = None
                print(err)
                break

            # reform slice
            # ============
            if pettmp is not None:
                for ilat, lat_indx in enumerate(range(lat_indx_min, lat_indx_max + 1)):
                    gran_lat = round((90.0 - self.latitudes[lat_indx])*NGRANULARITY)

                    for ilon, lon_indx in enumerate(range(lon_indx_min, lon_indx_max + 1)):
                        gran_lon = round( self.longitudes[lon_indx]*NGRANULARITY )
                        key = '{:0=5d}_{:0=5d}'.format(int(gran_lat), int(gran_lon))

                        # add data for this coordinate
                        # ============================
                        if varname == 'prate':
                            pettmp[varnam_map][key] = [round(val*numSecsDay, 3) for val in slice[:,ilat,ilon]]
                        elif varname == 'tmax':
                            pettmp[varnam_map][key] = [round(val, 3) for val in slice[:,ilat,ilon]]

            nc_dset.close()     # close netCDF file

        return pettmp

    def fetch_climate_NC_data(self, wthr_rsrce, wthr_set_defn, aoi_indices, climgen, num_band, future_flag = True):
        """
        get precipitation or temperature data for a given variable and lat/long index for requested times
        """
        func_name = __prog__ +  ' fetch_climate_NC_data'
        simplefilter('default')

        lat_indx_min, lat_indx_max, lon_indx_min, lon_indx_max = aoi_indices
        pettmp = {}
        if future_flag == True:
            precip_fname = self.fut_precip_fname
            tas_fname = self.fut_tas_fname
        else:
            precip_fname = self.hist_precip_fname
            tas_fname    = self.hist_tas_fname
        precip_var = wthr_set_defn['precip_var']
        tas_var = wthr_set_defn['tas_var']

        # year, month, dat (also hour, minute, second)
        # ===========================================
        strt_date = datetime(int(climgen.fut_start_year), 1, 1)
        end_date  = datetime(int(climgen.fut_end_year), 12, 31)

        # process future climate
        # ======================
        varnams_mapped = {precip_var: 'precipitation', tas_var: 'temperature'}
        varnams = sorted(varnams_mapped.keys())

        for varname, fname in zip(varnams, list([precip_fname, tas_fname])):
            varnam_map = varnams_mapped[varname]
            pettmp[varnam_map] = {}
            nc_dset = Dataset(fname, mode='r')

            if wthr_set_defn['resol_time'] == 'Daily':
                strt_indx = date2index(strt_date, nc_dset.variables['time'])
                end_indx  = date2index(end_date,  nc_dset.variables['time'])

            # collect readings for all time values
            # ====================================
            try:
                if wthr_set_defn['resol_time'] == 'Daily':
                    slice = nc_dset.variables[varname][strt_indx:end_indx+1, lat_indx_min:lat_indx_max + 1,
                                                                                lon_indx_min:lon_indx_max + 1]
                else:
                    slice = nc_dset.variables[varname][:, lat_indx_min:lat_indx_max + 1,
                                                                                lon_indx_min:lon_indx_max + 1]
            except RuntimeWarning as err:
                print(err)

            # reform slice
            # ============
            for ilat, lat_indx in enumerate(range(lat_indx_min, lat_indx_max + 1)):
                lat = self.latitudes[lat_indx]
                gran_lat = round((90.0 - lat)*NGRANULARITY)

                for ilon, lon_indx in enumerate(range(lon_indx_min, lon_indx_max + 1)):
                    lon = self.longitudes[lon_indx]
                    gran_lon = round((180.0 + lon) * NGRANULARITY)

                    # check mask
                    # ==========
                    lsmask = nc_dset.variables['lsmask'][lat_indx, lon_indx]
                    if lsmask == 0:
                        mess = 'Land sea mask is zero at lat: {} lon: {}'.format(lat, lon)
                        climgen.lgr.info(mess)
                        continue

                    key = '{:0=5d}_{:0=5d}'.format(int(gran_lat), int(gran_lon))

                    # identify masked elements
                    # ========================
                    nmasked = 0
                    for indx, val in enumerate(slice[:, ilat, ilon]):
                        if val is ma.masked:
                            nmasked += 1

                    if nmasked >  0:
                        mess = 'No usable weather at lat: {} lon: {}\t'.format(lat, lon)
                        mess += 'granular key: {}\tN masked: {}'.format(key, nmasked)
                        climgen.lgr.info(mess)
                        continue

                    # add data for this coordinate
                    # ============================
                    if varnam_map == 'precipitation':
                        if wthr_rsrce == 'ECLIPS45V2' or wthr_rsrce == 'ECLIPS85V2':
                            pettmp[varnam_map][key] = [float(val) for val in slice[:,ilat,ilon]]
                        else:
                            pettmp[varnam_map][key] = [round(val*numSecsDay, 3) for val in slice[:,ilat,ilon]]

                    if varnam_map == 'temperature':
                        if wthr_rsrce == 'ECLIPS45V2' or wthr_rsrce == 'ECLIPS85V2':
                            pettmp[varnam_map][key] = [float(val) for val in slice[:,ilat,ilon]]
                        else:
                            pettmp[varnam_map][key] = [round(val - 273.15, 3) for val in slice[:,ilat,ilon]]

            nc_dset.close()     # close netCDF file

        return pettmp

    def create_lta_avrg_met_file(self, clim_dir, lat_inp, site, lta_precip, lta_tmean):
        """
        use prexisting metyyyys.txt files to generate a text file of average weather which will subsequently
        be included in the input.txt file
        also create a climate file for each of the simulation years based on average weather from the CRU year range
        """
        func_name =  'create_lta_avrg_met_file'
        over_write_flag = True

        # ==========
        random_year = 2001
        lta_pet = thornthwaite(lta_tmean, lat_inp, random_year)
        site.lta_pet = [round(pet, 2) for pet in lta_pet]
        site.lta_precip = [round(precip, 2) for precip in lta_precip]
        site.lta_tmean = [round(tmean, 2) for tmean in lta_tmean]

        # write long term average climate file
        # ====================================
        met_ave_file = join(normpath(clim_dir), self.met_ave_file)
        if not isfile(met_ave_file) or over_write_flag:
            output = []
            imnth = 1
            for precip, pet, tmean in zip(site.lta_precip, site.lta_pet, site.lta_tmean):
                output.append([imnth, precip, pet, tmean])
                imnth += 1

            with open(met_ave_file, 'w', newline='') as fpout:
                writer = csv_writer(fpout, delimiter='\t')
                writer.writerows(output)

            self.lgr.info('Successfully wrote average weather file {}'.format(met_ave_file))

        return 0