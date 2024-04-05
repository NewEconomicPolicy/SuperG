"""
#-------------------------------------------------------------------------------
# Name:
# Purpose:     consist of high level functions invoked by main GUI
# Author:      Mike Martin
# Created:     11/12/2015
# Licence:     <your licence>
#-------------------------------------------------------------------------------
#
"""

__prog__ = 'crop_soil_fert_funcs.py'
__version__ = '0.0.1'
__author__ = 's03mm5'

import os
from math import ceil
from calendar import isleap
from locale import LC_ALL, setlocale, format_string
from pandas import read_csv, DataFrame

from pedotransfer import boelter, bss
from make_site_specific_files import Crop, FertiliserApplication, ManureApplication, Cultivation

sleepTime = 5
max_lines = 10

def _identify_separator(csv_file):
    """
    reads the first n (e.g. 5) lines and identifes file by:
        a) the delimiter used           - tab or comma
        b) the number of header records - 0 or 1
        c) the number of fields,
        d) for each field their type (float, integer or string)
        e) the range for each float or integer field
    """
    if not os.path.isfile(csv_file):
        print('File {} does not exist'.format(csv_file))
        return None

    # read first n lines to estimate number of lines and number of fields
    # ===================================================================
    fobj = open(csv_file, 'r')
    file_size_read = 0
    line_recs = []
    for nline in range(0, max_lines):
        raw_line = fobj.readline()

        # empty string indicates an end of line
        # =====================================
        if raw_line == '':
            break
        file_size_read += len(raw_line)
        line = raw_line.rstrip('\n')
        line_recs.append(line)
    fobj.close()

    if file_size_read == 0:
        print('File {} is empty'.format(csv_file))
        return -1

    # identify type of separator
    # ==========================
    num_tabs   = line.count('\t')
    num_commas = line.count(',')
    num_spaces = line.count(' ')
    separators = {'tabs': num_tabs, 'commas': num_commas, 'spaces': num_spaces}
    separator_maps = {'tabs': '\t', 'commas': ',', 'spaces': ' '}

    sep_name = max(separators, key=lambda key: separators[key])
    print(csv_file + ' has ' + sep_name + ' separator')

    return separator_maps[sep_name]

def read_sowing_dates_file(form, sowing_dates_fname = None):
    """
    # read the CSV of sowing dates
    # ============================
    """
    setlocale(LC_ALL, '')

    lookup_frame = DataFrame()
    headers = list(['lon','lat', 'year', 'sow_day', 'harv_day'])

    if sowing_dates_fname is None:
        sowing_dates_fname = form.w_lbl15.text()

    separator = _identify_separator(sowing_dates_fname)
    if separator == None:
        return lookup_frame, DataFrame()

    data_frame = read_csv(sowing_dates_fname, sep = separator, names = headers)
    nlines = len(data_frame)
    nlines_str = format_string('%d', nlines, grouping=True)
    print('Read ' + nlines_str + ' lines from sowing dates file')
    if 'lon' not in data_frame.columns or 'lat' not in data_frame.columns:
        print('Sowing dates file ' + sowing_dates_fname + ' must have fields lon and lat')
        return lookup_frame, DataFrame()

    # determine number of years
    # =========================
    frst_year = data_frame['year'].values[0]
    last_year = frst_year
    for ic, year in enumerate(data_frame['year'].values[1:]):
        if year == frst_year:
            break
        last_year = year

    nyears = ic + 1
    print('# years {}\tfirst: {}\tlast: {}'.format(nyears, frst_year, last_year))
    form.crop_type_sowing['frst_year'] = frst_year
    form.crop_type_sowing['last_year'] = last_year

    # create set of points to assist with location
    # ============================================
    lats = []; lons = []; lookup_ids = []
    recid= 0
    for year, lat, lon in zip(data_frame['year'], data_frame['lat'], data_frame['lon']):
        if year == frst_year:
            lats.append(lat)
            lons.append(lon)
            lookup_ids.append(recid)
        recid += 1

    lookup_frame['lat']  = lats
    lookup_frame['lon']  = lons
    lookup_frame['lookup_id']  = lookup_ids
    lookup_frame['point'] = [(y, x) for y, x in zip(lookup_frame['lat'], lookup_frame['lon'])]

    lats = sorted(lookup_frame['lat'].unique())
    lons = sorted(lookup_frame['lon'].unique())
    resols = []
    lat1 = lats[0]
    for lat2 in lats[1:]:
        step = lat2 - lat1
        resols.append(step)
        lat1 = lat2
    resols.sort()
    if len(resols) == 0:
        print('Could not determine latitude resolution of sowing dates file ' + sowing_dates_fname + ' please check file')
        return lookup_frame, DataFrame()

    print('Resolution: {}'.format(resols[0]))

    return lookup_frame, data_frame

def read_fertiliser_file(form):
    """
    # read the CSV of fertilisers
    # ===========================
    """
    fertiliser_fname = form.w_lbl14.text()
    data_frame = read_csv(fertiliser_fname, sep = '\t')
    nlines = len(data_frame)
    nlines_str = format_string('%d', nlines, grouping=True)
    print('Read ' + nlines_str + ' lines from fertiliser file')
    if 'lon' not in data_frame.columns or 'lat' not in data_frame.columns:
        print('Fertiliser file ' + fertiliser_fname + ' must have fields lon and lat')
        return DataFrame()

    # create set of points to assist with location
    # ============================================
    data_frame['point'] = [(y, x) for y, x in zip(data_frame['lat'], data_frame['lon'])]

    lats = sorted(data_frame['lat'].unique())
    lons = sorted(data_frame['lon'].unique())
    resols = []
    lat1 = lats[0]
    for lat2 in lats[1:]:
        step = lat2 - lat1
        resols.append(step)
        lat1 = lat2
    resols.sort()
    print('Resolution: {}'.format(resols[0]))

    return data_frame

def create_site_soil_layers(site, soil_list, depths):
    """
    rework soil list into an intelligible dictionary
    """
    # soc, Bulk density [g/cm3], pH ,  % clay by weight, % silt by weight, % sand by weight
    soil_metric_list = list(['soc', 'Bulk_Density', 'ph', 'clay', 'silt', 'sand'])

    if len(soil_list) == 7:
        num_lyrs = 1
        lyr_depths = depths[0:1]
    else:
        num_lyrs= 2
        lyr_depths = depths[:]

    soil = {}
    for metric in soil_metric_list:
        soil[metric] = []

    for lyr_num in range(num_lyrs):    # typically: 2
        for ic, metric in enumerate(soil_metric_list):
            strt_indx = 6*lyr_num
            soil[metric].append(soil_list[strt_indx + ic])

    site.nlyrs = num_lyrs
    site.lyr_depths = lyr_depths

    site.soil_name = 'xxxx'  # TODO: ideas?
    lu_names = site._luts

    # use first layer
    # ===============
    ilayer = 0
    site.iom_c = 0.0    # TODO: estimate Inert organic matter?

    for ilayer in range(site.nlyrs):
        for lu in range(len(lu_names)):
                site.add_soil_lyr(
                    lu_names[lu],
                    soil['soc'][ilayer],    # TODO: check this value, soil carbon content [kgC/ha]
                    soil['Bulk_Density'][ilayer],
                    soil['ph'][ilayer],
                    soil['clay'][ilayer],
                    soil['silt'][ilayer],
                    soil['sand'][ilayer]
                )

    # Sozanka equation
    # ================
    site.min_n_lvl = round(5 + (soil['clay'][0] * 10.0 / 65.0), 1)

    # Calculate soil water capacities
    # ===============================
    # TODO: check whether top soil or sub-soil required
    bulk_dens = soil['Bulk_Density'][0]
    org_carb = soil['soc'][0]/(bulk_dens * site.lyr_depths[0] * 1000.0)

    if org_carb > 30.0:
        # Peat, so use Boelter equations
        wc_wp, awc_fc, awc_sat = boelter(bulk_dens, site.soil_depth * 10.0)

    else:
        # Non-peat so use British Soil Service equations
        wc_wp, awc_fc, awc_sat = bss(
            soil['sand'][0], soil['silt'][0], soil['clay'][0], org_carb, bulk_dens, site.soil_depth * 10.0, True)


    site.wc_wp = round(wc_wp, 1)
    site.awc_fc = round(awc_fc, 1)
    site.awc_sat = round(awc_sat, 1)

    # only required for output
    # ========================
    site.soc  = soil['soc'][0]
    site.clay = soil['clay'][0]
    site.ph   = soil['ph'][0]
    site.bulk_dens = bulk_dens

    return

def site_data_modify_daily(site, latitude, longitude, crop_code, crop_name, manure_flag, fert_amnts,
                                                                sow_days, harv_days, yields, met_fnames):
    """
    daily
    """
    func_name =  __prog__ +  ' site_data_modify_daily'

    # for use in main loop
    # ===================
    start_from_1801 = site.start_from_1801
    if start_from_1801:
        start_year = 1801
    else:
        start_year = site.start_year

    if len(sow_days) == 0:
        sow_doy = 45
        harv_doy = 228
    else:
        if start_from_1801:
            sow_doy = sow_days[0]
            harv_doy = harv_days[0]
        else:
            sow_doy, harv_doy = None, None

    ntime_steps = 0
    for year in range(start_year, site.end_year):
        if isleap(year):
            ntime_steps += 366
        else:
            ntime_steps += 365

        site.ntime_steps = ntime_steps

    site.toc = 2.7     # total organic matter in top 50cm of soil [kg C/ha] varies between 0.39 to 5.27 (not used)
    site.iom_c = 0.0   # TODO: estimate this? Inert organic matter in top 50cm of soil [kg C/ha]
    site.met_fnames = met_fnames

    # transfer data from grid cell to site
    # ====================================
    site.latitude = round(float(latitude), 3)
    site.longitude = round(float(longitude), 3)

    site.drain_class = 2
    site.depth_imperm_lyr = 3
    site.wtr_tbl_dpth = 300
    site.prev_lu = 1  # Arable
    site.prev_crop_code = 1
    site.yield_prev_crop = yields[0]
    site.prev_crop_harvest_doy = 0      # must not be more than 12 in case of monthly timestep

    cult_type  = 0         # Type of cultivation: 0, 1, 2 or 3� see table B1.1.3 in user manual
    cult_vigor = 0.5       # Vigour of cultivation: 0.0 � 1.0

    # main loop: assume single annual crop
    # ====================================
    iyr = 0
    naccum_tsteps = 0
    for year in range(start_year, site.end_year):
        if sow_doy is None:
            sowing_ts  = sow_days[iyr]  + naccum_tsteps
            harvest_ts = harv_days[iyr] + naccum_tsteps
        else:
            sowing_ts  = sow_doy + naccum_tsteps
            harvest_ts = harv_doy + naccum_tsteps

        cult_ts     = sowing_ts - 7     # Tillage/cultivation one week before sowing
        frst_app_ts = sowing_ts - 1     # First fertiliser application one day before sowing
        scnd_app_ts = sowing_ts + 28    # Second fertiliser application 28 days after sowing

        crop = Crop()
        crop.code = crop_code
        crop.crop_name = crop_name
        crop.sowing_doy  = sowing_ts
        crop.harvest_doy = harvest_ts
        crop.n_uptake = fert_amnts[iyr]
        crop.exp_yield = yields[iyr]
        crop.residues_inc = 1

        fert_amnt = fert_amnts[iyr]
        if manure_flag:
            crop.nfert_apps = 0
            manure = ManureApplication(fert_amnt/4.0, frst_app_ts)
            crop.manure_apps.append(manure)
            manure = ManureApplication(3.0*fert_amnt/4.0, scnd_app_ts)
            crop.manure_apps.append(manure)
            crop.nmanure_apps = len(crop.manure_apps)
        else:
            # args: amount, app_doy, no3_pc, nh4_pc, urea_pc, non_amm_sulphate_salts=0, labelled=0
            # ====================================================================================
            crop.nmanure_apps = 0
            fert = FertiliserApplication(fert_amnt/4.0, frst_app_ts, 50, 50, 0)
            crop.fert_apps.append(fert)
            fert = FertiliserApplication(3.0*fert_amnt/4.0, scnd_app_ts, 50, 50, 0)
            crop.fert_apps.append(fert)
            crop.nfert_apps = len(crop.fert_apps)

        site.crops.append(crop)

        # cultivation
        # ===========
        cultivation = Cultivation(cult_ts, cult_type, cult_vigor)
        cultivation.cult_ts   =  cult_ts

        site.cultivations.append(cultivation)

        if isleap(year):
            ndays = 365     # TODO: make sure ECOSSE works on 365 days only
        else:
            ndays = 365

        naccum_tsteps += ndays
        iyr += 1

    # this stanza creates new list of met files for each year of simulation
    # =====================================================================
    if start_from_1801:
        nyears = len(site.crops)
        nsets = ceil(nyears/len(met_fnames))

        all_met = []
        for iset in range(nsets):
            all_met += met_fnames

        n_chop = len(all_met) - nyears
        site.met_fnames = all_met[n_chop:]
        
        # should be same
        # ==============
        site.nyears = nyears
        site.ncrops = len(site.crops)

    return site

def site_data_modify_mnthly(site, latitude, longitude, crop_code, crop_name, ndepos,
                    manure_flag, manure_code, manure_amnts, fert_flag, fert_amnts, graze_flag, graze_n_amnts, yields):
    """
    monthly
    """
    func_name = __prog__ + ' site_data_modify_mnthly'

    start_year = site.start_year
    site.ntime_steps = site.nyears * 12

    # site.toc = 2.7  # total organic matter in top 50cm of soil [kg C/ha] varies between 0.39 to 5.27 (not used)
    site.iom_c = 0.0  # TODO: estimate this? Inert organic matter in top 50cm of soil [kg C/ha]

    # transfer data from grid cell to site
    # ====================================
    site.latitude = round(float(latitude), 3)
    site.longitude = round(float(longitude), 3)
    # site.atmos_n_dep = round(sum(ndepos)/len(ndepos), 3)
    site.atmos_n_dep = round(ndepos[-1], 3)

    site.drain_class = 2
    site.depth_imperm_lyr = 3
    site.wtr_tbl_dpth = 300
    site.prev_lu = 1  # Arable
    site.prev_crop_code = crop_code
    site.yield_prev_crop = yields[0]
    site.prev_crop_harvest_doy = 0  # must not be more than 12 in case of monthly timestep
    site.crops = []
    site.cultivations = []

    cult_type = 0  # Type of cultivation: 0, 1, 2 or 3� see table B1.1.3 in user manual
    cult_vigor = 0.5  # Vigour of cultivation: 0.0 � 1.0

    # main loop: assume single annual crop
    # ====================================
    naccum_tsteps = 0
    for iyr, year in enumerate(range(start_year, site.end_year + 1)):
        if year >= site.sim_strt_yr:
            mnr_flag = manure_flag
            grz_flag = graze_flag
            frt_flag = fert_flag
        else:
            mnr_flag = True
            grz_flag = True
            frt_flag = True

        sowing_ts = 1 + naccum_tsteps
        harvest_ts = 12 + naccum_tsteps
        cult_ts = sowing_ts
        apply_ts_frst_n = sowing_ts + 3     # April
        apply_ts_scnd_n = sowing_ts + 5     # June
        apply_ts_frst_mnr = sowing_ts + 6   # July
        apply_ts_scnd_mnr = sowing_ts + 7   # August

        crop = Crop()
        crop.code = crop_code
        crop.crop_name = crop_name
        crop.sowing_doy = sowing_ts
        crop.harvest_doy = harvest_ts
        crop.exp_yield = yields[iyr]
        crop.residues_inc = 1

        mnr_n_amnt, grz_n_amnt, fert_amnt = 3*[0]

        if frt_flag and fert_amnts is not None:
            # args: amount, app_doy, no3_pc, nh4_pc, urea_pc, non_amm_sulphate_salts=0, labelled=0
            # ====================================================================================
            fert_amnt = fert_amnts[iyr]
            if fert_amnt > 0:
                fert = FertiliserApplication(fert_amnt / 4.0, apply_ts_frst_n, 50, 50, 0)
                crop.fert_apps.append(fert)
                fert = FertiliserApplication(3.0 * fert_amnt / 4.0, apply_ts_scnd_n, 50, 50, 0)
                crop.fert_apps.append(fert)

        if grz_flag or mnr_flag:
            mnr_n_amnt = 0
            # add grazing N to manure
            # =======================
            if grz_flag:
                mnr_n_amnt += graze_n_amnts[iyr]
            if mnr_flag:
                mnr_n_amnt += manure_amnts[iyr]

            if mnr_n_amnt > 0:
                manure = ManureApplication(mnr_n_amnt / 2.0, apply_ts_frst_mnr, manure_code)
                crop.manure_apps.append(manure)
                manure = ManureApplication(mnr_n_amnt / 2.0, apply_ts_scnd_mnr, manure_code)
                crop.manure_apps.append(manure)

        crop.nfert_apps = len(crop.fert_apps)
        crop.nmanure_apps = len(crop.manure_apps)
        crop.n_uptake = fert_amnt + grz_n_amnt + mnr_n_amnt  # TODO: check with MA and MR

        site.crops.append(crop)

        # cultivation
        # ===========
        cultivation = Cultivation(cult_ts, cult_type, cult_vigor)
        cultivation.cult_ts = cult_ts

        site.cultivations.append(cultivation)

        naccum_tsteps += 12

        iyr += 1

    return site