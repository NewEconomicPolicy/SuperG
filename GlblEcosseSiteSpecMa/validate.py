#-------------------------------------------------------------------------------
# Name:        validate.py
# Purpose:     Class to organise run parameters
# Author:      Mike Martin
# Created:     28/09/2015
# Licence:     <your licence>
#-------------------------------------------------------------------------------
#!/usr/bin/env python

__prog__ = 'validate.py'
__version__ = '0.0.0'

# Version history
# ---------------
# 
from os.path import isfile
from csv import writer as csv_writer

ERROR_STR = '*** Error *** '

def grid_info_hdr():
    """

    """
    hdr  = '{:8s}\t{:8s}\t'.format('Nsyn', 'amnt')
    hdr += '{:8s}\t{:8s}\t{:8s}\t{:8s}\t'.format('Ngraz', 'amnt', 'ar grass', 'ar cell')
    hdr += '{:8s}\t{:8s}\t'.format('manure', 'amnt')
    hdr += '{:8s}\t{:8s}\t{:8s}\n'.format('lat', 'lon', 'gran_coord')

    return hdr

def nitro_info_hdr(climgen):
    """

    """
    hdr  = '{:8s}\t{:8s}\t{:8s}\t'.format('lat', 'lon', 'gran_coord')
    for yr in range(climgen.fut_start_year, climgen.fut_end_year + 1):
        hdr += '\t{:s}'.format(str(yr))
    hdr += '\n'
    
    return hdr

def report_nitro_inpts_lite(fobj, lat, lon, gran_lat, gran_lon, strt_yr_indx, syn_frt_n_amnts, grz_n_amnts, mnr_amnts):
    """
    report data availability for this grid cell
    """
    gran_coord = '{0:0=5g}_{1:0=5g}'.format(gran_lat, gran_lon)
    mess = '{:.2f}\t{:.2f}\t{:s}'.format(lat, lon, gran_coord)

    for val in [(grz_n + frt_n + mnr_n) for grz_n, frt_n, mnr_n in zip(grz_n_amnts[strt_yr_indx:],
                                                    syn_frt_n_amnts[strt_yr_indx:], mnr_amnts[strt_yr_indx:])]:
        mess += '\t{:.3f}'.format(val)

    fobj.write(mess + '\n')

    return

def report_nitro_inpts(fobj, lat, lon, gran_lat, gran_lon, strt_yr_indx, nyears,
                                                                manure_flag, n_mnr_amnts, n_fert_flag, n_fert_amnts):
    """
    report data availability for this grid cell
    """
    gran_coord = '{0:0=5g}_{1:0=5g}'.format(gran_lat, gran_lon)
    mess = '{:.2f}\t{:.2f}\t{:s}'.format(lat, lon, gran_coord)

    if not manure_flag:
        n_mnr_amnts = nyears*[0]

    if not n_fert_flag or n_fert_amnts is None:
        n_fert_amnts = nyears*[0]

    for val in [(mnr + nfrt) for mnr, nfrt in zip(n_mnr_amnts[strt_yr_indx:], n_fert_amnts[strt_yr_indx:])]:
        mess += '\t{:.2f}'.format(val)

    fobj.write(mess + '\n')

    return

def report_data_avail(fobj, counters, syn_n_amnts, graze_n_amnts, lu_areas, manure_amnts,
                      lat, lon, gran_lat, gran_lon):
    """
    report data availability for this grid cell
    """
    gran_coord = '{0:0=5g}_{1:0=5g}'.format(gran_lat, gran_lon)
    counters.all += 1

    syn_n_state, graze_n_state, manure_state = 3*['ok']
    if all(item == 0 for item in syn_n_amnts):
        syn_n_state = 'no data'
        counters.n_syn += 1

    cell_area = lu_areas['area_cell']
    graze_area = lu_areas['area_graz']
    if all(item == 0 for item in graze_n_amnts):
        graze_n_state = 'no data'
        counters.n_grz += 1

    if all(item == 0 for item in manure_amnts):
        manure_state = 'no data'
        counters.manure += 1

    mess = '{:8s}\t{:8.2f}\t'.format(syn_n_state, syn_n_amnts[-1])
    mess += '{:8s}\t{:8.3f}\t{:8.1f}\t{:8.1f}\t'.format(graze_n_state, graze_n_amnts[-1], graze_area, cell_area)
    mess += '{:8s}\t{:8.2f}\t'.format(manure_state, manure_amnts[-1])
    mess += '{:8.2f}\t{:8.2f}\t{:8s}\n'.format(lat, lon, gran_coord)
    fobj.write(mess)

    return

def check_harv_doys_adjust(sow_days, harv_days):
    '''
    add 365 days if sowing date later than harvest
    '''
    if sow_days[0] > harv_days[0]:
        harv_days_adust = [doy + 365 for doy in harv_days]
    else:
        harv_days_adust = harv_days

    return harv_days_adust

def check_phenology_fnames(form):
    '''
    fertiliser and crop input files should look include the crop type in their names
    '''

    sowing_crop_type = 'unknown'
    fert_crop_type = 'unknown';
    fert_field = ''
    form.crop_type_fert = {'type': fert_crop_type, 'field': fert_field}
    form.crop_type_sowing = {'type': sowing_crop_type, 'frst_year': -999, 'last_year': -999}

    fert_fname = form.w_lbl14.text()
    if not isfile(fert_fname):
        return 'Fertiliser input file ' + fert_fname + ' does not exist'
    try:
        with open(fert_fname, 'r') as fobj:
            print('Fertiliser file ' + fert_fname + ' is readable')

    except (OSError, IOError) as err:
        print(err)
        return 'Could not read fertiliser input file'

    fert_crop_type = 'spring wheat'
    # fert_crop_type = 'maize'

    # sowing dates
    # ============
    sowing_fname = form.w_lbl15.text()
    if not isfile(sowing_fname):
        return 'Sowing_dates input file ' + sowing_fname + ' does not exist'
    try:
        with open(sowing_fname, 'r') as fobj:
            print('Sowing dates input file ' + sowing_fname + ' is readable')

    except (OSError, IOError) as err:
        print(err)
        return 'Could not read sowing dates input file'

    sowing_crop_type = 'spring wheat'
    #    sowing_crop_type = 'maize'

    form.w_create_files.setEnabled(False)
    # form.w_create_files.setEnabled(True)

    form.crop_type_fert = {'type': fert_crop_type, 'field': fert_field}
    form.crop_type_sowing = {'type': sowing_crop_type, 'frst_year': -999, 'last_year': -999}

    # disable simulation years
    # ========================
    if form.years_from_flag == 'sowing':
        form.combo11s.setEnabled(False)
        form.combo11e.setEnabled(False)

    return 'tout est bon'

def check_fert_sowing_fnames(form):
    '''
    fertiliser and crop input files should look include the crop type in their names
    '''

    sowing_crop_type = 'unknown'
    fert_crop_type = 'unknown'; fert_field = ''
    form.crop_type_fert   = {'type': fert_crop_type,   'field': fert_field}
    form.crop_type_sowing = {'type': sowing_crop_type, 'frst_year': -999, 'last_year': -999}

    fert_fname = form.w_lbl14.text()
    if not isfile(fert_fname):
        return 'Fertiliser input file ' + fert_fname + ' does not exist'
    try:
        with open(fert_fname, 'r') as fobj:
            print('Fertiliser file ' + fert_fname + ' is readable')

    except (OSError, IOError) as err:
            print(err)
            return 'Could not read fertiliser input file'

    if fert_fname.lower().rfind('swhe') > 0:
        fert_crop_type = 'spring wheat'
        fert_field = 'SWHE'
    elif fert_fname.lower().rfind('maiz') > 0:
        fert_crop_type = 'maize'
        fert_field = 'MAIZ'

    # sowing dates
    # ============
    sowing_fname = form.w_lbl15.text()
    if not isfile(sowing_fname):
        return 'Sowing_dates input file ' + sowing_fname + ' does not exist'
    try:
        with open(sowing_fname, 'r') as fobj:
            print('Sowing dates input file ' + sowing_fname + ' is readable - simulation years will be disabled')

    except (OSError, IOError) as err:
            print(err)
            return 'Could not read sowing dates input file'
    
    if sowing_fname.lower().rfind('swheat') > 0:
        sowing_crop_type = 'spring wheat'
    elif sowing_fname.lower().rfind('maize') > 0:
        sowing_crop_type = 'maize'

    # check file contents
    # ===================
    mess = ''
    if fert_crop_type == 'unknown' or sowing_crop_type == 'unknown':
        mess = ERROR_STR + 'unknown crop type or sowing crop type'
        if fert_crop_type == 'unknown':
            print(ERROR_STR + 'fertiliser crop type is unknown - cannot enable simulation file creation')
        if sowing_crop_type == 'unknown':
            print(ERROR_STR + 'sowing crop type is unknown - cannot enable simulation file creation')
        form.w_create_files.setEnabled(False)
    else:
        mess = 'fertiliser for crop ' + fert_crop_type + ' and sowing dates for crop ' + sowing_crop_type + \
                                                                                            ' input files are valid'
        form.w_create_files.setEnabled(True)

    form.crop_type_fert   = {'type': fert_crop_type,   'field': fert_field}
    form.crop_type_sowing = {'type': sowing_crop_type, 'frst_year': -999, 'last_year': -999}

    # disable simulation years
    # ========================
    if form.years_from_flag == 'sowing':
        form.combo11s.setEnabled(False)
        form.combo11e.setEnabled(False)

    return mess