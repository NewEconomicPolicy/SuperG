#-------------------------------------------------------------------------------
# Name:
# Purpose:
# Author:      s03mm5
# Created:     08/12/2015
# Copyright:   (c) s03mm5 2015
# Licence:     <your licence>
#-------------------------------------------------------------------------------
#

__version__ = '1.0.00'
__prog__ = 'prepare_ecosse_files.py'

import csv
from os import makedirs
from os.path import lexists, join, normpath, basename, exists, split as split_dir
from time import time
import sys
from calendar import isleap
from shutil import copyfile
from copy import copy
from glob import glob

from thornthwaite import thornthwaite
from make_site_specific_files import MakeSiteFiles
from crop_soil_fert_funcs import site_data_modify_daily, site_data_modify_mnthly, create_site_soil_layers
from mngmnt_fns_and_class import (fetch_vals_from_dset, fetch_contrived_phenology, fetch_vals_from_graze_dset,
                                        fetch_manure_code, fetch_yields_from_dset)
from validate import report_data_avail, report_nitro_inpts, report_nitro_inpts_lite

from glbl_ecss_cmmn_funcs import write_kml_file, write_signature_file, write_manifest_file, input_txt_line_layout

NGRANULARITY = 120

_MONTHDAYS = [31,28,31,30,31,30,31,31,30,31,30,31]
_LEAP_MONTHDAYS = [31,29,31,30,31,30,31,31,30,31,30,31]
SPACER_LEN = 12

# syn_fert_flag, manure_flag, graze_flag
# ======================================
from collections import namedtuple
RunCnfg = namedtuple("RunCnfg", "syn_fert_flag manure_flag graze_flag")
RUN_CNFGS = {'all': RunCnfg(True, True, True), 'nsyn_grz': RunCnfg(True, False, True),
            'nsyn_mnr': RunCnfg(True, True, False), 'grz_mnr': RunCnfg(False, True, True),
                                                                        'no_inpts': RunCnfg(False, False, False)}

def _populate_equilib_period(fut_strt_year, strt_yr, syn_fert_n_amnts, graze_n_amnts, manure_amnts):
    """
	During 1700-1960 use the same inputs (synthetic N, slurry, and grazing slurry) as we use for 1961 in the
    all-treatments (baseline) scenario.
    """
    nequil_yrs = fut_strt_year - strt_yr

    if manure_amnts is None:
        res_mnr = None
    else:
        val_use = manure_amnts[nequil_yrs]
        res_mnr = nequil_yrs * [val_use] + manure_amnts[nequil_yrs:]

    if syn_fert_n_amnts is None:
        res_fert = None
    else:
        val_use = syn_fert_n_amnts[nequil_yrs]
        res_fert = nequil_yrs * [val_use] + syn_fert_n_amnts[nequil_yrs:]

    if graze_n_amnts is None:
        res_grz_n = None
    else:
        val_use = graze_n_amnts[nequil_yrs]
        res_grz_n = nequil_yrs * [val_use] + graze_n_amnts[nequil_yrs:]

    return res_fert, res_grz_n, res_mnr

def _make_hist_ave_wthr(climgen, pettmp_hist, wthr_set, lat, lon):
    """
    calculate historic average weather
    """
    dset_strt_year = wthr_set['start_year']

    hist_strt_year = climgen.hist_start_year
    indx_start = 12*(hist_strt_year - dset_strt_year)

    hist_end_year = climgen.hist_end_year
    indx_end   = 12*(hist_end_year - dset_strt_year + 1) # end year includes all 12 months - TODO: check

    # use dict-comprehension to initialise precip. and temperature dictionaries
    # =========================================================================
    hist_precip = {mnth: 0.0 for mnth in climgen.months}
    hist_tmean  = {mnth: 0.0 for mnth in climgen.months}

    for indx in range(indx_start, indx_end, 12):

        for imnth, month in enumerate(climgen.months):
            try:
                hist_precip[month] += pettmp_hist['precipitation'][indx + imnth]
                hist_tmean[month]  += pettmp_hist['temperature'][indx + imnth]
            except (IndexError, TypeError) as err:
                print('\n' + str(err) + 'indx: {}\timnth: {}\tat lat/lon: {} {}'.format(indx, imnth, lat, lon))
                return

    # write stanza for input.txt file consisting of long term average climate
    # =======================================================================
    hist_wthr_recs = []
    num_hist_years = hist_end_year - hist_strt_year + 1
    hist_lta_precip = []
    for month in climgen.months:
        ave_precip = hist_precip[month]/num_hist_years
        data_txt = '{}'.format(round(ave_precip,1))
        comment ='{} long term average monthly precipitation [mm]'.format(month)
        hist_wthr_recs.append(input_txt_line_layout(data_txt, comment))
        hist_lta_precip.append(ave_precip)

    hist_lta_tmean = []
    for month in climgen.months:
        ave_tmean = hist_tmean[month]/num_hist_years
        data_txt = '{}'.format(round(ave_tmean,2))
        comment = '{} long term average monthly temperature [degC]'.format(month)
        hist_wthr_recs.append(input_txt_line_layout(data_txt, comment))
        hist_lta_tmean.append(ave_tmean)

    return hist_wthr_recs, hist_lta_precip, hist_lta_tmean

def make_ecosse_files(form, climgen, site_rec, province, fert_defn, fert_era_defn, yield_defn,
                        graze_defn, ndepos_defn, pettmp_hist, pettmp_fut, counters, lu_areas):
    """
    generate sets of Ecosse files for each site
    where each site has one or more soils and each soil can have one or more dominant soils
    pettmp_grid_cell is climate data for this soil grid point
    """
    func_name = 'make_ecosse_files'

    # simulation start and end year
    # =============================

    end_yr = climgen.fut_end_year
    fut_strt_yr = climgen.fut_start_year
    if climgen.start_from_1801:
        strt_year = 1801
        strt_yr_indx = fut_strt_yr - strt_year
    else:
        strt_year = fut_strt_yr
        strt_yr_indx = 0

    ntot_years = end_yr - strt_year + 1

    # crop and weather
    # ================
    manure_code = fetch_manure_code(form.combo13, form.manure_types)
    study = form.w_study.text()
    crop_name = form.combo12.currentText()
    crop_code = form.crop_codes[crop_name]
    wthr_rsrce = climgen.wthr_rsrce
    if wthr_rsrce not in list(['ECLIPS45V2', 'ECLIPS85V2']):
        print('Weather resource {} not recognised in function {}'.format(wthr_rsrce, func_name))
        return

    if not isinstance(pettmp_hist, dict) or not isinstance(pettmp_hist, dict):
        print('Bad input to ' + func_name)
        return

    # ==============================
    gran_lat, gran_lon, lat, lon, area, mu_globals_props = site_rec
    lat_rnd = round(lat, 3)
    lon_rnd = round(lat, 3)
    sims_dir = form.sttngs['sims_dir']
    fut_clim_scen = climgen.fut_clim_scen

    # N deposition
    # ============
    ndepos = fetch_vals_from_dset(ndepos_defn, lat, lon, strt_year, end_yr)

    # Synthetic fertiliser - use coarser dataset as a fallback
    # ========================================================
    syn_fert_n_amnts = fetch_vals_from_dset(fert_era_defn, lat, lon, strt_year, end_yr, 'Nmineral')
    if all(item == 0 for item in syn_fert_n_amnts):
        syn_fert_n_amnts = fetch_vals_from_dset(fert_defn, lat, lon, strt_year, end_yr, 'Nmineral')

    # Manure
    # ======
    manure_amnts = fetch_vals_from_dset(fert_defn, lat, lon, strt_year, end_yr, 'Nmanure')

    # retrieve N from grazing and aggregate with manure
    # =================================================
    graze_n_amnts = fetch_vals_from_graze_dset(form.lgr, graze_defn, lu_areas, lat, lon,
                                                                            form.req_resol_deg, strt_year, end_yr)

    syn_fert_n_amnts, graze_n_amnts, manure_amnts = _populate_equilib_period(climgen.fut_start_year,
                                                            strt_year, syn_fert_n_amnts, graze_n_amnts, manure_amnts)

    # report data for this grid cell
    # ==============================
    report_data_avail(form.fobjs['grid'], counters, syn_fert_n_amnts, graze_n_amnts, lu_areas,
                                                                            manure_amnts, lat, lon, gran_lat, gran_lon)
    report_nitro_inpts_lite(form.fobjs['nitro'], lat, lon, gran_lat, gran_lon, strt_yr_indx,
                                                                        syn_fert_n_amnts, graze_n_amnts, manure_amnts)
    if form.w_report.isChecked():
        return

    # Yields
    # ======
    yields = fetch_yields_from_dset(yield_defn, lat, lon, strt_year, end_yr, climgen, pettmp_fut)
    if all(item == 0 for item in yields):
        counters.yields += 1
        mess = 'Yields for grid cell at lat: {}\tlong: {} are zero - will skip'.format(lat_rnd, lon_rnd)
        form.lgr.info(mess); print(mess)
        return

    sow_days, harv_days = fetch_contrived_phenology(ntot_years, 'monthly')

    # calculate historic average weather
    # ==================================
    wthr_set = form.weather_sets[wthr_rsrce + '_Mnth']
    hist_wthr_recs, hist_lta_precip, hist_lta_tmean = _make_hist_ave_wthr(climgen, pettmp_hist, wthr_set, lat, lon)

    # write a single set of met files for all simulations for this grid cell
    # ======================================================================
    gran_coord = '{0:0=5g}_{1:0=5g}'.format(gran_lat, gran_lon)
    met_rel_path = '..\\..\\' + climgen.rgn_wthr_dir + '\\' + gran_coord + '\\'
    clim_dir = normpath(join(sims_dir, climgen.rgn_wthr_dir, gran_coord))
    met_fnames = make_met_files(clim_dir, lat, climgen, pettmp_fut)     # future weather
    nyears = ntot_years

    site = MakeSiteFiles(form, climgen, comments = True)
    # TODO: these need sorted
    site.equil_mode = form.w_equimode.text()
    site.soil_code = 1
    site.prev_crop_harvest_doy = 0

    site.nyears = nyears
    site.ncrops = nyears # one crop per year
    site.ncultivations = nyears # one cultivation per year

    # create long-term average met file
    # =================================
    irc = climgen.create_lta_avrg_met_file(clim_dir, lat, site, hist_lta_precip, hist_lta_tmean)
    site.start_year = strt_year
    site.sim_strt_yr = climgen.fut_start_year
    site.end_year = end_yr

    #------------------------------------------------------------------
    # Create a set of simulation input files for each dominant
    # soil-land use type combination
    #------------------------------------------------------------------
    # construct directory name with all dominant soils

    for pair in mu_globals_props.items():
        mu_global, proportion = pair
        area_for_soil = area*proportion
        soil_list = form.hwsd_mu_globals.soil_recs[mu_global]

        for soil_num, soil in enumerate(soil_list):
            identifer = 'lat{0:0=7d}_lon{1:0=7d}_mu{2:0=5d}_s{3:0=2d}'.format(gran_lat, gran_lon, mu_global, soil_num + 1)

            # =======
            for sim_nm in RUN_CNFGS:
                syn_fert_flag, manure_flag, graze_flag = RUN_CNFGS[sim_nm].syn_fert_flag, \
                                                         RUN_CNFGS[sim_nm].manure_flag, RUN_CNFGS[sim_nm].graze_flag

                sim_dir = join(sims_dir, study + '_' + sim_nm, identifer)

                if not lexists(sim_dir):
                    makedirs(sim_dir)

                this_site = copy(site)
                create_site_soil_layers(this_site, soil, form.depths)

                this_site = site_data_modify_mnthly(this_site, lat, lon, crop_code, crop_name, ndepos, manure_flag,
                    manure_code, manure_amnts, syn_fert_flag, syn_fert_n_amnts, graze_flag, graze_n_amnts, yields)

                this_site.write_sim_files(sim_dir, soil, lat, hist_wthr_recs, met_rel_path)

                # write kml file if requested and signature file
                # ==============================================
                if form.sttngs['kml_flag'] and soil_num == 0:
                    write_kml_file(sim_dir, str(mu_global), mu_global, lat, lon)

                write_signature_file(sim_dir, mu_global, soil, lat, lon, province)

                # copy across Model Switches and fnames files
                # ===========================================
                for ecss_key in list(['model_switches', 'fnames', 'crop_sun']):
                    ecss_fn = form.dflt_ecss_fns[ecss_key]
                    copyfile(ecss_fn, join(sim_dir, basename(ecss_fn)))

                # manifest file is essential for subsequent processing
                # ====================================================
                write_manifest_file(study, fut_clim_scen, sim_dir, soil_list, mu_global, lat, lon, area_for_soil)

    # end of Soil loop
    # ================

    return

def make_met_files(clim_dir, lat, climgen, pettmp_fut, skip_fut_wthr_flag = True):
    """
    feed annual temperatures to Thornthwaite equations to estimate Potential Evapotranspiration [mm/month]
    """
    strt_year = climgen.fut_start_year
    end_year = climgen.fut_end_year

    # do not repeat met file creation if already in existence
    # =======================================================
    if lexists(clim_dir):
        if skip_fut_wthr_flag:
            met_fnames = [split_dir(met_fname)[1] for met_fname in glob(clim_dir + '/met*s.txt')]
            nyears = end_year - strt_year + 1
            if len(met_fnames) >= nyears:
                return met_fnames
    else:
        makedirs(clim_dir)

    precip = pettmp_fut['precipitation']
    temp   = pettmp_fut['temperature']

    indx1 = 0
    met_fnames = []
    for year in range(strt_year, end_year + 1):
        fname = 'met{}s.txt'.format(year)
        met_fnames.append(fname)
        met_path = join(clim_dir, fname)

        if climgen.fut_wthr_mnthly_flag:
            ntime_incrs = 12
        else:
            if isleap(year):
                ntime_incrs = 366
            else:
                ntime_incrs = 365

        indx2 = indx1 + ntime_incrs

        # precipitation and temperature
        precipitation = precip[indx1:indx2]
        temp_mean     = temp[indx1:indx2]

        # pet
        # ===
        if climgen.fut_wthr_mnthly_flag:
            pet = thornthwaite(temp_mean, lat, year)
        else:
            pet = _get_thornthwaite(temp_mean, lat, year)

        # TODO: do something about occasional runtime warning...
        pot_evapotrans = [round(p, 2) for p in pet]
        precip_out     = [round(p, 2) for p in precipitation]
        tmean_out      = [round(t, 2) for t in temp_mean]

        # write file
        output = []
        for tstep, mean_temp in enumerate(tmean_out):
            output.append([tstep+1, precip_out[tstep], pot_evapotrans[tstep], mean_temp])

        with open(met_path, 'w', newline='') as fpout:
            writer = csv.writer(fpout, delimiter='\t')
            writer.writerows(output)
            fpout.close()

        indx1 += ntime_incrs

    return  met_fnames

def _get_thornthwaite(temp_mean, lat, year):
    """
    feed daily annual temperatures to Thornthwaite equations to estimate Potential Evapotranspiration [mm/month]
    """
    func_name =  __prog__ + ' _get_thornthwaite'

    ntime_steps = len(temp_mean)

    if ntime_steps == 365:
        month_days = _MONTHDAYS
    else:
        month_days =  _LEAP_MONTHDAYS

    indx1 = 0
    monthly_t = []
    for ndays in month_days:
        indx2 = indx1 + ndays
        accum_temp = 0.0
        for temp in temp_mean[indx1:indx2]:
            accum_temp +=temp
        monthly_t.append(accum_temp/ndays)
        indx1 = indx2

    pet_mnthly = thornthwaite(monthly_t, lat, year)

    # now inflate pet to daily
    # ========================
    pet_daily = []
    for pet_month, ndays in zip(pet_mnthly, month_days):
        pet_day = pet_month/ndays
        for ic in range(ndays):
            pet_daily.append(pet_day)

    # print('{} {} {} {}'.format(year, ntime_steps, len(pet_daily),func_name))
    return pet_daily

def update_progress(last_time, start_time, completed, est_num_sims, no_climate, no_pasture):
    """
    Update progress bar             num_meta_cells, no_climate, no_pasture)
    """
    new_time = time()
    if new_time - last_time > 5:
        # used to be: Estimated remaining
        mess = '\rCompleted: {:<5d}\tRemaining: {:<5d}\tNo climate: {:<4d}\tNot pasture: {:<5d}'\
                                    .format(completed, est_num_sims - completed, no_climate, no_pasture)
        sys.stdout.flush()
        sys.stdout.write(mess)
        last_time = new_time

    return last_time