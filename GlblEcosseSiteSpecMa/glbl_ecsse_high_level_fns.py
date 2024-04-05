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

__prog__ = 'glbl_ecsse_high_level_fns.py'
__version__ = '0.0.1'
__author__ = 's03mm5'

from time import time
from operator import itemgetter
from copy import copy

from getClimGenNC import ClimGenNC
import hwsd_bil_v1

from hwsd_mu_globals_fns import gen_grid_cells_for_band
from prepare_ecosse_files import update_progress, make_ecosse_files
from commonCmpntsGUI import calculate_grid_cell
from getClimGenFns import check_clim_nc_limits, associate_climate
from mngmnt_fns_and_class import check_mask_location, create_proj_data_defns, open_proj_NC_sets, close_proj_NC_sets
from glbl_ecss_cmmn_funcs import create_dump_files
from validate import grid_info_hdr, nitro_info_hdr

MASK_VARNAMES = ['cropland', 'forest', 'grassland', 'other', 'pasture']
ERROR_STR = '*** Error *** '

from collections import namedtuple
counters = namedtuple("MngmtInpts", "all yields n_syn n_grz manure")

def _fetch_future_wthr(wthr_rsrce, climgen, wthr_set_defn, lookup_key, site_rec, num_band):
    """

    """
    climgen.lgr.info('Retrieving future weather for lookup key: {}'.format(lookup_key))

    bbox = [site_rec[3], site_rec[2]]
    aoi_indices, dummy = climgen.genLocalGrid(bbox, snglPntFlag=True)

    if wthr_rsrce == 'NASA':
        pettmp_fut = climgen.fetch_nasa_NC_data(aoi_indices, num_band)
    
    elif wthr_rsrce == 'EObs':
        pettmp_fut = climgen.fetch_eobs_NC_data(wthr_set_defn, aoi_indices, num_band)
    
    elif wthr_rsrce == 'EWEMBI':
        pettmp_fut = climgen.fetch_ewembi_NC_data(aoi_indices, num_band)
    
    else:
        pettmp_fut = climgen.fetch_climate_NC_data(wthr_rsrce, wthr_set_defn, aoi_indices, climgen, num_band)

    pettmp_fut_grid_cell = {}
    for metric in pettmp_fut:
        pettmp_fut_grid_cell[metric] = pettmp_fut[metric][lookup_key]

    return  pettmp_fut_grid_cell

def _generate_ecosse_files(form, climgen, mask_defn, fert_defn, fert_era_defn, yield_defn, graze_defn, ndepos_defn,
                           counters, num_band):
    """
    Main loop for generating ECOSSE outputs
    """
    func_name =  __prog__ + '\t_generate_ecosse_files'

    fut_strt_yr = climgen.fut_start_year
    fut_end_yr = climgen.fut_end_year

    study = form.w_study.text()
    print('Gathering soil and climate data for study ' + study + '...')

    max_cells = int(form.w_max_cells.text())

    # instantiate a soil grid and climate objects
    hwsd = hwsd_bil_v1.HWSD_bil(form)

    # add requested grid resolution attributes to the form object
    calculate_grid_cell(form, hwsd.granularity)
    bbox = form.bbox
    nvals_read = hwsd.read_bbox_hwsd_mu_globals(bbox, form.hwsd_mu_globals, form.req_resol_upscale)

    # retrieve dictionary consisting of mu_globals (keys) and number of occurences (values)
    mu_globals = hwsd.get_mu_globals_dict()
    if len(mu_globals) == 0:
        print('No soil records for this area\n')
        return

    mess = 'Retrieved {} values  of HWSD grid consisting of {} rows and {} columns: ' \
          '\n\tnumber of unique mu_globals: {}'.format(nvals_read, hwsd.nlats, hwsd.nlons, len(mu_globals))
    form.lgr.info(mess)

    # for each grid point in the band return quintuples of integer Lat/Lon 30 seconds coords, Lat/Lons and mu_global
    # ==============================================================================================================
    hwsd.bad_muglobals = form.hwsd_mu_globals.bad_mu_globals
    aoi_res, bbox = gen_grid_cells_for_band(hwsd, form.req_resol_upscale)
    if form.w_use_high_cover.isChecked():
        aoi_res = _simplify_aoi(aoi_res)

    lon_ll_aoi, lat_ll_aoi, lon_ur_aoi, lat_ur_aoi = bbox
    n_meta_cells = len(aoi_res)
    print('Band extent - LL lat/lon: {} {}\tUR lat/lon: {} {}\tN meta cells: {}'
                            .format(lat_ll_aoi, lon_ll_aoi, lat_ur_aoi, lon_ur_aoi, n_meta_cells))
    if n_meta_cells == 0:
        mess = 'No aoi_res recs therefore unable to create simulation files... \n'
        print(mess); form.lgr.info(mess)
        return

    # open NC datasets
    # ================
    open_proj_NC_sets(mask_defn, yield_defn, ndepos_defn, fert_defn, fert_era_defn, graze_defn)

    # user feedback
    # =============
    mess = 'Band {} will generate up to {} grid cells '.format(num_band, n_meta_cells)
    if form.w_use_dom_soil.isChecked():
        mess += 'and simulations'
    else:
        # 4.5 = estimated mean number of dominant soils per cell
        # ======================================================
        mess += 'and an estimated {} simulations'.format(int(n_meta_cells * 4.5))

    form.lgr.info(mess); print(mess)

    # generate weather dataset indices which enclose the AOI for this band
    # ====================================================================
    ret_code = climgen.genLocalGrid(bbox)
    if ret_code is None:
        return
    aoi_indices_fut, aoi_indices_hist = ret_code

    wthr_rsrce = climgen.wthr_rsrce
    if climgen.fut_wthr_mnthly_flag:
        wthr_set_fut_defn = form.weather_sets[wthr_rsrce + '_Mnth']
    else:
        wthr_set_fut_defn = form.weather_sets[wthr_rsrce + '_Day']

    print('Retrieving weather data')
    #      ========================
    wthr_set_hist_defn = form.weather_sets[wthr_rsrce + '_Mnth']
    if wthr_rsrce  == 'NASA':
        pettmp_hist = climgen.fetch_nasa_NC_data(aoi_indices_hist, num_band, future_flag = False)

    elif wthr_rsrce  == 'EObs':
        pettmp_hist = climgen.fetch_eobs_NC_data(wthr_set_hist_defn, aoi_indices_fut, num_band, future_flag = False)

    elif wthr_rsrce  == 'EWEMBI':
        pettmp_hist = climgen.fetch_ewembi_NC_data(aoi_indices_hist, num_band, future_flag = False)

    else:
        pettmp_hist = climgen.fetch_climate_NC_data(wthr_rsrce, wthr_set_hist_defn,
                                                         aoi_indices_hist, climgen, num_band, future_flag = False)

    pettmp_fut = climgen.check_ecss_wthr_data(form.sttngs['sims_dir'], climgen, pettmp_hist, num_band)
  
    print('Creating simulation files for years {} to {} for band {}...'.format(fut_strt_yr, fut_end_yr, num_band))
    #      =========================================
    last_time = time()
    start_time = time()
    ncompleted = 0
    no_climate = 0
    no_pasture = 0
    max_cells_flag = False
    copy_fut_from_hist = True

    # generate sets of Ecosse files for each site where each site has one or more soils
    # each soil can have one or more dominant soils
    # =======================================================================
    for site_rec in aoi_res:
        last_time = update_progress(last_time, start_time, ncompleted, n_meta_cells, no_climate, no_pasture)

        lu_found, lu_areas = check_mask_location(mask_defn, site_rec, MASK_VARNAMES, form.req_resol_deg)
        if not lu_found:
            no_pasture += 1
            continue

        # yield_set = associate_yield(form)
        ret_code = associate_climate(site_rec, climgen, pettmp_hist, pettmp_fut)
        if ret_code is None:
            no_climate += 1
            continue

        pettmp_hist_grid_cell, pettmp_fut_grid_cell, lookup_key = ret_code

        # generate or pull simulation weather as required
        # ===============================================
        if copy_fut_from_hist:
            pettmp_fut_grid_cell = copy(pettmp_hist_grid_cell)
        else:
            pettmp_fut_grid_cell = _fetch_future_wthr(wthr_rsrce, climgen, wthr_set_fut_defn,
                                                                                lookup_key, site_rec, num_band)

        make_ecosse_files(form, climgen, site_rec, study, fert_defn, fert_era_defn, yield_defn, graze_defn,
                                    ndepos_defn, pettmp_hist_grid_cell, pettmp_fut_grid_cell, counters, lu_areas)
        ncompleted += 1
        if ncompleted >= max_cells:
            max_cells_flag = True
            break

    # tidy up and leave
    # =================
    close_proj_NC_sets(mask_defn, yield_defn, ndepos_defn, fert_defn, fert_era_defn, graze_defn)

    mess = '\nSummary for band: {}\tsites generated: {}'.format(num_band, ncompleted)
    mess += '\tsites with no pasture: {}   no climate: {}'.format(no_pasture, no_climate)
    print(mess)

    return max_cells_flag

def generate_banded_sims(form):
    """
    called from GUI
    """
    study = form.w_study.text()

    if form.hwsd_mu_globals == None:
        print('Undetermined HWSD aoi - please select a valid HSWD csv file')
        return -1

    if form.w_use_dom_soil.isChecked():
        use_dom_soil_flag = True
    else:
        use_dom_soil_flag = False

    if form.w_strt_1801.isChecked():
        start_from_1801 = True
    else:
        start_from_1801 = False

    # write header
    # ============
    counters.all, counters.yields, counters.n_syn, counters.n_grz, counters.manure = 5*[0]
    create_dump_files(form)
    header  = '{:8s}\t{:8s}\t'.format('Nsyn', 'amnt')
    header += '{:8s}\t{:8s}\t{:8s}\t{:8s}\t'.format('Ngraz', 'amnt', 'ar grass', 'ar cell')
    header += '{:8s}\t{:8s}\t'.format('manure', 'amnt')
    header += '{:8s}\t{:8s}\t{:8s}\n'.format('lat', 'lon', 'gran_coord')
    form.fobjs['grid'].write(header)

    # verify mask, sowing, yields, fertiliser NC files
    # ================================================
    crop_name = form.combo12.currentText()
    proj_data_defns = create_proj_data_defns(form.sttngs['proj_path'], form.sttngs['mask_fn'], crop_name)
    if proj_data_defns is None:
        print(ERROR_STR + 'verifing NC files for study ' + study)
        return -1

    mask_defn, yield_defn, pheno_defn, fert_defn, fert_era_defn, graze_defn, ndepos_defn = proj_data_defns
    del (proj_data_defns)

    # user defined bounding box
    # =========================
    lon_ll_gui = float(form.w_ll_lon.text())
    lat_ll_gui = float(form.w_ll_lat.text())
    lon_ur_gui = float(form.w_ur_lon.text())
    lat_ur_gui = float(form.w_ur_lat.text())
    form.bbox =  list([lon_ll_gui, lat_ll_gui, lon_ur_gui, lat_ur_gui])

    # lat_ll_aoi is the most southerly, i.e. least latitude, of the HWSD aoi
    # ============================================================
    lat_ll_aoi = form.hwsd_mu_globals.lat_ll_aoi
    lon_ll_aoi = form.hwsd_mu_globals.lon_ll_aoi
    lat_ur_aoi = form.hwsd_mu_globals.lat_ur_aoi
    lon_ur_aoi = form.hwsd_mu_globals.lon_ur_aoi
    bbox_aoi = list([lon_ll_aoi,lat_ll_aoi,lon_ur_aoi,lat_ur_aoi])

    # extract required values from the HWSD database and simplify if requested
    # ========================================================================
    hwsd = hwsd_bil_v1.HWSD_bil(form)
    soil_recs = hwsd.get_soil_recs(form.hwsd_mu_globals.mu_global_list)  # list is already sorted
    form.hwsd_mu_globals.soil_recs = _simplify_soil_recs(hwsd.bad_muglobals, soil_recs, use_dom_soil_flag)

    form.hwsd_mu_globals.bad_mu_globals = [0] +  hwsd.bad_muglobals
    del(hwsd)

    # weather choice - CRU is default
    # ==============
    wthr_rsrce = form.combo10w.currentText()

    # check requested AOI coordinates against extent of the weather resource dataset
    # ==============================================================================
    if check_clim_nc_limits(form, wthr_rsrce, bbox_aoi):
        print('Selected ' + wthr_rsrce)
    else:
        return -1

    climgen = ClimGenNC(form, start_from_1801)
    form.fobjs['grid'].write(grid_info_hdr())
    form.fobjs['nitro'].write(nitro_info_hdr(climgen))

    # main banding loop
    # =================
    start_at_band = form.start_at_band
    stop_after_band = form.stop_after_band
    print('Starting at band {}\tstopping after band: {}'.format(start_at_band, stop_after_band))

    lat_step = 0.5
    nsteps = int((lat_ur_gui - lat_ll_gui)/lat_step) + 1
    lat_ur = copy(lat_ur_gui)
    for isec in range(nsteps):
        lat_ll_new = lat_ur - lat_step
        if lat_ll_new < lat_ll_gui:
            lat_ll_new = lat_ll_gui     # avoid generating grid cells beyond the user specified southerly extent

        num_band = isec + 1
        if num_band > stop_after_band:
            break

        # if the latitude floor of the band has not reached the ceiling of the HWSD aoi then skip this band
        # =================================================================================================
        if lat_ll_new > form.hwsd_mu_globals.lat_ur_aoi or num_band < start_at_band:
            print('Skipping out of area band {} of {} with latitude extent of min: {}\tmax: {}'
                                                    .format(num_band, nsteps, round(lat_ll_new,6), round(lat_ur, 6)))
        else:
            form.bbox = list([lon_ll_gui, lat_ll_new, lon_ur_gui, lat_ur])
            print('\nProcessing band {} of {} with latitude extent of min: {}\tmax: {}'
                                                    .format(num_band, nsteps, round(lat_ll_new,4), round(lat_ur, 4)))
            # does actual work
            # ================
            max_cells_flag = _generate_ecosse_files(form, climgen, mask_defn, fert_defn, fert_era_defn,
                                                    yield_defn, graze_defn, ndepos_defn, counters, num_band)
            if max_cells_flag:
                break

        # check to see if the last band is completed
        # ==========================================
        if lat_ll_aoi > lat_ll_new or num_band == nsteps:
            break

        lat_ur = lat_ll_new

    print('\nFinished processing after {} bands of latitude extents\n'.format(isec + 1))
    #        ===================
    if hasattr(form, 'fstudy'):
        for ichan in range(len(form.fstudy)):
            form.fstudy[ichan].close()

    print('Total # grid cells: {}\tno yield: {}\tno Nsyn: {}\tno Ngraze: {}\tno manure: {}\n'.format(counters.all,
                                                    counters.yields, counters.n_syn, counters.n_grz, counters.manure))
    for key in form.fobjs:
        form.fobjs[key].close()
        print('Wrote: ' + form.fobjs[key].name)

    return 0

def _simplify_soil_recs(bad_muglobals, soil_recs, use_dom_soil_flag):
    """
    compress soil records if duplicates are present
    simplify soil records if requested
    each mu_global points to a group of soils
    a soil group can have up to ten soils
    """
    func_name =  __prog__ + ' _simplify_soil_recs'

    # clean out 7000, 7003
    # ====================
    for mu_global in bad_muglobals:
        if mu_global in soil_recs:
            del soil_recs[mu_global]

    num_raw = 0 # total number of sub-soils
    num_compress = 0 # total number of sub-soils after compressions

    new_soil_recs = {}
    for mu_global in soil_recs:

        # no processing necessary
        # =======================
        num_sub_soils = len(soil_recs[mu_global])
        num_raw += num_sub_soils
        if num_sub_soils == 1:
            num_compress += 1
            new_soil_recs[mu_global] = soil_recs[mu_global]
            continue

        # check each soil for duplicates
        # ==============================
        new_soil_group = []
        soil_group = sorted(soil_recs[mu_global])

        first_soil = soil_group[0]
        metrics1 = first_soil[:-1]
        share1   = first_soil[-1]
        for soil in soil_group[1:]:
            metrics2 = soil[:-1]
            share2 =   soil[-1]
            if metrics1 == metrics2:
                share1 += share2
            else:
                new_soil_group.append(metrics1 + [share1])
                metrics1 = metrics2
                share1 = share2

        new_soil_group.append(metrics1 + [share1])
        num_sub_soils = len(new_soil_group)
        num_compress += num_sub_soils
        if num_sub_soils == 1:
            new_soil_recs[mu_global] = new_soil_group
            continue

        if use_dom_soil_flag:
            # assign 100% to the first entry of sorted list
            # =============================================
            dom_soil = copy(sorted(new_soil_group, reverse = True, key=itemgetter(-1))[0])
            dom_soil[-1] = 100.0
            new_soil_recs[mu_global] = list([dom_soil])

    mess = '\nExiting {}\trecords in: {} out: {}'.format(func_name, len(soil_recs),len(new_soil_recs))
    print(mess + '\tnum raw sub-soils: {}\tafter compression: {}'.format(num_raw, num_compress))
    return new_soil_recs

def _simplify_aoi(aoi_res):
    """
    simplify AOI records
    """
    aoi_res_new = []
    j = 0
    for site_rec in aoi_res:
        content = site_rec[-1]
        npairs = len(content)
        if npairs == 0:
            print('No soil information for AOI cell {} - will skip'.format(site_rec))
        elif npairs == 1:
            aoi_res_new.append(site_rec)
        else:
            site_rec_list = list(site_rec)  # convert tuple to a list so we can edit last element
            new_content = sorted(content.items(), reverse = True, key = itemgetter(1))  # sort content so we can pick up most dominant mu_global
            total_proportion = sum(content.values())    # add up proportions
            site_rec_list[-1] = {new_content[0][0]: total_proportion}       # create a new single mu global with summed proportions

            aoi_res_new.append(tuple(site_rec_list)) # convert list to tuple

    return aoi_res_new
