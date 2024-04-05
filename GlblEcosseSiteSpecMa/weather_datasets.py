#-------------------------------------------------------------------------------
# Name:        weather_datasets.py
# Purpose:     script to create weather object and other functions
# Author:      Mike Martin
# Created:     31/07/2015
# Licence:     <your licence>
#-------------------------------------------------------------------------------
#!/usr/bin/env python

__prog__ = 'weather_datasets.py'
__version__ = '0.0.0'

# Version history
# ---------------
# 
from os import remove
from os.path import join, normpath, isfile, lexists
from numpy import zeros
from netCDF4 import Dataset, num2date
from glob import glob
from unidecode import unidecode

daycent_null = -99.99   # null value for Daycent
sleepTime = 5

def _average_slice(slice):

    ndays, nlats, nlons = slice.shape
    new_slice = zeros(ndays)

    for indx_day in range(ndays):
        day_slice = slice[indx_day, :, :]
        new_val = day_slice.mean()
        new_slice[indx_day] = new_val

    return new_slice

def _open_output_files(admin_div, out_dir, start_year, end_year):

    # delete files if exists, then open
    # =================================
    admin_div = unidecode(admin_div.replace(' ','_'))     # remove spaces and awkward characters
    fname = join(out_dir, admin_div + '_' + str(start_year) + '_' + str(end_year) + '.100')
    clim_file = normpath(fname)
    if isfile(clim_file):
        remove(clim_file)
        print('\nDeleted ' + clim_file)

    fhand_clim = open(clim_file, 'w')

    return fhand_clim, clim_file

def _write_slices_to_file(slices, fhand_clim, clim_file):
    '''
    write all records for this slice set to a tab separated ASCII file
    '''

    # step through slices
    #====================
    print('Writing to climate file ' + clim_file)
    for ic, date_set in enumerate(slices['daycent_dates']):
        tmax = round(slices['ds_tmax'][ic], 1)
        tmin = round(slices['ds_tmin'][ic], 1)
        precip = round(slices['ds_precip'][ic]/10.0, 2)     # convert to cms
        rec = date_set + list([tmax, tmin, precip, daycent_null, daycent_null, daycent_null])
        str_rec = [str(val) for val in rec]             # use list comprehension
        fhand_clim.write( '\t'.join(str_rec) + '\n' )

    # close file
    # ===========
    fhand_clim.close()
    print('Wrote {} records to climate file {}'.format(ic, clim_file))

    return

def _fetch_weather_nc_parms(nc_fname, rsrc_name, wthr_rsrce, resol_time, time_var_name = 'time'):
    '''
    create a data record and lat/lon lists from weather datasets
    '''
    nc_fname = normpath(nc_fname)
    nc_dset = Dataset(nc_fname, 'r')
    time_var = nc_dset.variables[time_var_name]
    if 'calendar' in time_var.ncattrs():
        calendar_attr = time_var.calendar
    else:
        calendar_attr = 'standard'

    # standard names
    # ==============
    if rsrc_name == 'NASA' or rsrc_name == 'EObs':
        lat = 'latitude'
        lon = 'longitude'
    else:
        lat = 'lat'
        lon = 'lon'

    lat_var = nc_dset.variables[lat]
    lon_var = nc_dset.variables[lon]

    # create lat/lon lists taking into account particularities of each resource
    # =========================================================================
    if rsrc_name == 'EObs':
        lats = [round(float(lat), 3) for lat in lat_var]
        lons = [round(float(lon), 3) for lon in lon_var]
    else:
        lats = [float(lat) for lat in lat_var]
        lons = [float(lon) for lon in lon_var]

    # bounding box
    # ============
    lon_ll = min(lons)
    lon_ur = max(lons)
    lat_ll = min(lats)
    lat_ur = max(lats)

    # resolutions
    # ===========
    '''
    resol_lon = round((lon_var[-1] - lon_var[0])/(len(lon_var) - 1), 5)
    resol_lat = round((lat_var[-1] - lat_var[0])/(len(lat_var) - 1), 5)
    '''
    resol_lon = (lons[-1] - lons[0])/(len(lons) - 1)
    resol_lat = (lats[-1] - lats[0])/(len(lats) - 1)
    if abs(resol_lat) != abs(resol_lon):
        print('Warning - weather resource {} has different lat/lon resolutions: {} {}'
                                                        .format(wthr_rsrce, resol_lat, resol_lon))

    # Get the start and end date of the time series (as datetime objects):
    # ====================================================================
    if rsrc_name == 'NASA':
        time_var_units = time_var.units.replace('0 01-','0-')   # patch for NASA datasets
    else:
        time_var_units = time_var.units
    start_day = int(time_var[0].item())
    try:
        start_date = num2date(start_day, units = time_var_units, calendar = calendar_attr)
    except (TypeError) as err:
        print('Error deriving start and end year for dataset: ' + nc_fname)
        return None, None, None

    end_day = int(time_var[-1].item())
    end_date = num2date(end_day, units = time_var_units, calendar = calendar_attr)

    # convert time variable to Daycent format, namely day of month, month, year and	Julian day (Day of year)
    # =======================================
    daycent_dates = []
    '''
    if resol_time == 'Monthly':
        last_year = daycent_null
        for time_val in time_var:
            date_val = num2date(time_val, units = time_var_units, calendar = calendar_attr)

            this_year = date_val.year
            if this_year != last_year:
                julian_day = 1
                last_year = this_year
            else:
                julian_day += 1

            daycent_dates.append([date_val.day, date_val.month, date_val.year, julian_day])
    '''
    nc_dset.close()

    data_rec = {'start_year': start_date.year,  'end_year': end_date.year,
            'resol_lat': resol_lat, 'lat_frst': lats[0], 'lat_last': lats[-1], 'lat_ll': lat_ll, 'lat_ur': lat_ur,
            'resol_lon': resol_lon, 'lon_frst': lons[0], 'lon_last': lons[-1], 'lon_ll': lon_ll, 'lon_ur': lon_ur,
                    'resol_time': resol_time,  'daycent_dates': daycent_dates}

    print('{} start and end year: {} {}\tresolution: {} degrees'
            .format(wthr_rsrce, data_rec['start_year'],  data_rec['end_year'], abs(data_rec['resol_lat'])))

    return data_rec, lons, lats

def read_weather_dsets_detail(form):
    '''
    ascertain the year span for historic datasets
    potential weather dataset resources are:
    NASA - Monthly:  from 1980-01-16 to 2010-12-16      Daily:  from 1980-01-01 to 2010-12-31
    EObs - Monthly:  from 1980-01-31 to 2017-12-31      Daily:  from 1950-01-01 to 2017-12-31
    '''
    weather_sets = {}
    weather_dir = form.weather_dir

    form.scenarios = sorted(list(['RCP6.0', 'RCP4.5', 'RCP2.6','RCP8.5']))
    valid_wthr_dset_rsrces = []

    precip_vars = list(['Precipalign', 'Precipalign'])
    tas_vars = list(['Tairalign', 'Tairalign'])
    root_dirs = list(['ECLIPS45V2_45', 'ECLIPS85V2_85'])  # TODO: NASA function requires finishing
    glob_strs = list(['1961_2100', '1961_2100'])

    period_names   = list(['Monthly', 'Daily'])
    period_abbrevs = list(['Mnth',    'Day'])

    for root_dir, glob_str, precip_var, tas_var in zip(root_dirs, glob_strs, precip_vars, tas_vars):

        rsrc_name = root_dir.split('_')[0]

        # check daily and monthly dsets
        # =============================
        resource_valid_flag = True
        for period_name, period_abbrev in zip(period_names, period_abbrevs):

            period_dir = weather_dir + '\\' + root_dir + '\\' + period_name
            if lexists(period_dir):

                wthr_rsrce = rsrc_name + '_' + period_abbrev
                nc_fnames = glob(period_dir + '/*' + glob_str + '*.nc')
                if len(nc_fnames) > 0:
                    data_rec, lons, lats = _fetch_weather_nc_parms(nc_fnames[0], rsrc_name, wthr_rsrce, period_name)
                    if data_rec is None:
                        break
                    else:
                        weather_sets[wthr_rsrce] = data_rec
                        weather_sets[wthr_rsrce]['longitudes'] = lons
                        weather_sets[wthr_rsrce]['latitudes']  = lats
                        weather_sets[wthr_rsrce]['base_dir']   = period_dir
                        weather_sets[wthr_rsrce]['ds_precip']  = nc_fnames[0]
                        weather_sets[wthr_rsrce]['precip_var'] = precip_var
                        weather_sets[wthr_rsrce]['ds_tas']     = nc_fnames[1]
                        weather_sets[wthr_rsrce]['tas_var']    = tas_var
                else:
                    print('No ' + rsrc_name + ' ' + period_name.lower() + ' datasets present in ' + period_dir)
                    resource_valid_flag = False
                    break
            else:
                print(rsrc_name + ' ' + period_name.lower() + ' folder ' + period_dir + ' does not exist')
                if period_abbrev == 'Mnth':
                    resource_valid_flag = False
                    break

        if resource_valid_flag:
            valid_wthr_dset_rsrces.append(rsrc_name)

    form.valid_wthr_dset_rsrces = sorted(valid_wthr_dset_rsrces)
    form.weather_sets = weather_sets

    print('')
    return

def _patch_nasa_time_units(form):
    '''
    patch to overcome bug in NASA data - use only once, not called
    '''
    import netCDF4 as cdf
    from glob import glob

    nasa_dir = 'F:\\GlobalEcosseData\\NASA_NCs'
    nasa_fnames = glob(nasa_dir + '/*AgMERRA*.nc')
    for fname in nasa_fnames:
        nc_dset = Dataset(fname,'r+')

        time_var = nc_dset.variables['time']
        raw_time_units = time_var.units
        time_var.units = raw_time_units.replace('0 01-','0-')

        nc_dset.close()

    return

def report_aoi_size(form, lon_ll, lat_ll, lon_ur, lat_ur):
    '''
    write ASCII climate files
    '''
    func_name =  __prog__ + ' report_aoi_size'

    # this will be initially only NASA
    # ================================
    resource = form.combo10w.currentText()
    weather_set = form.weather_sets[resource]
    resol_lat = weather_set['resol_lat']
    lat0 = weather_set['lat0']
    resol_lon = weather_set['resol_lon']
    lon0 = weather_set['lon0']

    lat_indx_ll = int(round((lat_ll - lat0)/resol_lat))
    lon_indx_ll = int(round((lon_ll - lon0)/resol_lon))

    lat_indx_ur = int(round((lat_ur - lat0)/resol_lat))
    lon_indx_ur = int(round((lon_ur - lon0)/resol_lon))

    lat_indx_min = min(lat_indx_ll, lat_indx_ur)
    lat_indx_max = max(lat_indx_ll, lat_indx_ur)
    nlats = lat_indx_max - lat_indx_min + 1

    lon_indx_min = min(lon_indx_ll, lon_indx_ur)
    lon_indx_max = max(lon_indx_ll, lon_indx_ur)
    nlons = lon_indx_max - lon_indx_min + 1

    # get slice for each dataset metric
    # =================================
    mess = 'will retrieve weather for {} locations - nlats/nlons: {} x {} '.format(nlats*nlons, nlats, nlons)

    print(mess)

    return
