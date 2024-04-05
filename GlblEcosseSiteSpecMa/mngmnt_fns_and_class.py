#-------------------------------------------------------------------------------
# Name:        mngmnt_fns_and_class.py
# Purpose:     script to create objects describing NC data sets
# Author:      Mike Martin
# Created:     31/05/2020
# Licence:     <your licence>
#-------------------------------------------------------------------------------
#!/usr/bin/env python

__prog__ = 'mngmnt_fns_and_class.py'
__version__ = '0.0.0'

# Version history
# ---------------
# 
from os.path import join, normpath, isfile, splitext
from netCDF4 import Dataset
from glob import glob
from numpy import seterr, ma, arange, isnan, all as np_all, nan as np_nan
from numpy.ma.core import MaskedConstant

from glbl_ecss_bio_fns import MiamiDyceAnnualNpps
from shape_funcs import calculate_area

ERROR_STR = '*** Error *** '
NDEP_FACTOR =  10000 / 1000000    # factor for N deposition - from: mg N m-2 yr-1    to: kg N ha-1

# factor to convert to tonnes C ha-1 for yields
# =============================================
NYLDS_GRASS_FACTOR = 1000               # Grass Biomass Use: kgDM/ha
# NYLDS_NPP_FACTOR = 10000 / 1000000     # NPP Grass: g C m-2
NYLDS_NPP_FACTOR = 0.0073333 # from 10 * 100 / (1000 * 45) * (33 / 100)    NPP Grass: g C m-2

NDAYS_YEAR = 365
KM2_TO_HA = 100

CROP_NAME_DICT = {'Grassland': 'Grassland', 'Winter Wheat': 'Wheat_Winter', 'Spring Wheat': 'Wheat_Spring',
                        'Grain Maize': 'Maize', 'Spring Oats': 'Oats_Spring', 'Spring Barley': 'Barley_Spring'}

def integrate_n_manure(manure_amnts, graze_n_amnts):
    """

    """
    n_manure_amnts = [amnt for amnt in manure_amnts]
    graze_aggr_mnr_amnts = []

    # add cattle, goats and sheep contributions
    # =========================================
    for indx in range(len(manure_amnts)):
        all_anmls_amnt = 0.0
        for amnl_type in graze_n_amnts:
            all_anmls_amnt += graze_n_amnts[amnl_type][indx]

        n_manure_amnts[indx] +=  all_anmls_amnt
        graze_aggr_mnr_amnts.append(all_anmls_amnt)

    return n_manure_amnts, graze_aggr_mnr_amnts


def integrate_n_fert(syn_fert_n_amnts, graze_n_amnts):
    """

    """
    n_fert_amnts = [amnt for amnt in syn_fert_n_amnts]
    graze_aggr_n_amnts = []

    # add cattle, goats and sheep contributions
    # =========================================
    for indx in range(len(n_fert_amnts)):
        all_amnt = 0.0
        for amnl_type in graze_n_amnts:
            all_amnt += graze_n_amnts[amnl_type][indx]

        n_fert_amnts[indx] +=  all_amnt
        graze_aggr_n_amnts.append(all_amnt)

    return n_fert_amnts, graze_aggr_n_amnts

def fetch_manure_code(combo13, manure_types):
    """

    """
    manure_type = combo13.currentText()
    mnr_indx = manure_types['Manure Type'].index(manure_type)
    manure_code = manure_types['Code'][mnr_indx]

    return manure_code

def fetch_contrived_phenology(ntot_years, time_step ='daily'):
    """

    """
    if time_step == 'daily':
        sow_times = ntot_years*[60]
        harv_times = ntot_years * [235]
    else:
        sow_times = ntot_years * [1]
        harv_times = ntot_years * [12]

    return sow_times, harv_times

def check_mask_location(mask_defn, site_rec, land_uses, resol_deg):
    """
    resol_deg is size of cell in degrees
    masked elements indicate zero, i.e. supposing the landuse is  pasture, then if all the elements are zero then the
                                                                                        cell is counted as not pasture
    """
    gran_lat, gran_lon, lat, lon, dummy, dummy = site_rec

    res_d2 = resol_deg/2

    lat_ur = lat + res_d2
    lon_ur = lon + res_d2
    lat_ur_indx, lon_ur_indx, ret_code = mask_defn.get_nc_coords(lat_ur, lon_ur)

    lat_ll = lat - res_d2
    lon_ll = lon - res_d2
    lat_ll_indx, lon_ll_indx, ret_code = mask_defn.get_nc_coords(lat_ll, lon_ll)
    nsub_cells = (lat_ur_indx - lat_ll_indx + 1)*(lon_ur_indx - lon_ll_indx + 1)

    # area in km**2 of cell
    # =====================
    bbox = list([0.0, lat_ll, resol_deg, lat_ur])
    area_cell = calculate_area(bbox)

    lu_areas = {}
    for land_use in land_uses:
        vals = mask_defn.nc_dset.variables[land_use][lat_ll_indx:lat_ur_indx + 1, lon_ll_indx:lon_ur_indx + 1]
        nvals = vals.count()
        if ma.is_masked(vals):
            if nvals > 0:       # Count the non-masked elements of the array along the given axis
                lu_found = True
        else:
            val_mean = vals.mean()
            if val_mean > 0:
                pass

        lu_areas[land_use] = nvals

    ncrop = lu_areas['cropland']
    ngraz = lu_areas['grassland'] + lu_areas['pasture']

    if ngraz == 0:
        lu_found_rc = False
        lu_areas['prop_graz'] = 0
        lu_areas['area_graz'] = 0
    else:
        lu_found_rc = True
        lu_areas['prop_graz'] = ngraz / (ncrop + ngraz)
        lu_areas['area_graz'] = area_cell*(ngraz/nsub_cells)

    lu_areas['area_cell'] = area_cell

    return lu_found_rc, lu_areas

def _calc_area_graz_cell(resol_deg, lat):
    """

    """
    res_d2 = resol_deg/2
    lat_ur = lat + res_d2
    lat_ll = lat - res_d2

    # area in km**2 of cell
    # =====================
    bbox = list([0.0, lat_ll, resol_deg, lat_ur])
    area_cell = calculate_area(bbox)

    return area_cell

def fetch_vals_from_graze_dset(lggr, graze_defn, lu_areas, lat, lon, resol_deg, sim_strt_yr, sim_end_yr):
    """
    read data from graze NC dataset from FAO which is derived from 10km grid cell
    NB  this dataset has no time dimension
    """
    func_name = 'fetch_vals_from_graze_dset'

    area = lu_areas['area_graz']    # km**2

    area_graz = _calc_area_graz_cell(graze_defn.resol_lat, lat)
    ratio_areas = area/area_graz

    res_d2 = resol_deg / 2

    lat_ur = lat + res_d2
    lon_ur = lon + res_d2
    lat_ur_indx, lon_ur_indx, ret_code = graze_defn.get_nc_coords(lat_ur, lon_ur)

    lat_ll = lat - res_d2
    lon_ll = lon - res_d2
    lat_ll_indx, lon_ll_indx, ret_code = graze_defn.get_nc_coords(lat_ll, lon_ll)
    ncells = (lat_ur_indx - lat_ll_indx + 1) * (lon_ur_indx - lon_ll_indx + 1)

    varnames = graze_defn.varnames
    val_tmp = 0
    for varname in varnames:
        vals_raw = graze_defn.nc_dset.variables[varname][lat_ll_indx:lat_ur_indx + 1, lon_ll_indx:lon_ur_indx + 1]
        val_mean = vals_raw.mean()
        val = val_mean * ratio_areas / KM2_TO_HA  # convert to kg per Ha

        if type(val) is MaskedConstant:
            val = 0.0
            lggr.info('MaskedConstant in {} lat: {} lon: {}'.format(func_name, round(lat,5), round(lon,5)))

        val_tmp += float(val)

    vals = (sim_end_yr - sim_strt_yr + 1) * [val_tmp]      # expand data

    return vals

def _check_array(vals_raw):
    """
    if all values are nan then return list of zeros
    """
    vals = None
    for val in vals_raw:
        if not isnan(val):
            vals = [round(val_raw.item(), 5) for val_raw in vals_raw]
            break

    if vals is None:
        vals = len(vals_raw)*[0]

    return vals

def fetch_yields_from_dset(metric_defn, lat, lon, sim_strt_yr, sim_end_yr, climgen, pettmp_fut):
    """
    read data from NC file
    """
    lat_indx, lon_indx, ret_code = metric_defn.get_nc_coords(lat, lon)

    miami_obj = MiamiDyceAnnualNpps(climgen, pettmp_fut)

    if metric_defn.varname is None:
        varname = metric_defn.varnames[0]
    else:
        varname = metric_defn.varname

    vals_raw = metric_defn.nc_dset.variables[varname][:, lat_indx, lon_indx]

    if metric_defn.varnames[0] == 'Grass_Biomass_Use':
        orchidee_npps = [val_raw.item()*NYLDS_GRASS_FACTOR for val_raw in vals_raw]
    else:
        orchidee_npps = [val_raw.item()*NYLDS_NPP_FACTOR for val_raw in vals_raw]    # NPP Grass

    # stretch values to cover simulation period
    # =========================================
    nyrs_add = metric_defn.start_year - sim_strt_yr - 1
    if nyrs_add > 0:
        orchidee_npps = nyrs_add * [orchidee_npps[0]] + orchidee_npps      # prepend first year of retrieved data
    elif nyrs_add < 0:
        nyrs_clip = abs(nyrs_add)
        orchidee_npps = orchidee_npps[nyrs_clip:]      # clip years from start of retrieved data

    # obtain MIAMI NPP in 2012
    # ========================
    miami_indx = metric_defn.end_year - miami_obj.strt_yr
    miami_npp = miami_obj.npps[miami_indx]

    orchidee_miami_ratio = orchidee_npps[-1]/miami_npp

    nyrs_add = miami_obj.end_yr - metric_defn.end_year

    orchidee_npps = orchidee_npps + [miami_npp * orchidee_miami_ratio for miami_npp in miami_obj.npps[miami_indx:]]

    return orchidee_npps

def fetch_vals_from_dset(metric_defn, lat, lon, sim_strt_yr, sim_end_yr, varname = None):
    """
    read data from NC file depending on resource
    NB  cnvrsn_fctr = metric_defn.cnvrsn_fctr
        vals = [round(val_raw.item() * cnvrsn_fctr, 3) for val_raw in vals_raw]
    """
    lat_indx, lon_indx, ret_code = metric_defn.get_nc_coords(lat, lon)
    rsrce = metric_defn.rsrce

    if varname is None:
        if metric_defn.varname is None:
            varname = metric_defn.varnames[0]
        else:
            varname = metric_defn.varname

    if varname is None:
        print(ERROR_STR + 'cannot ascertain NetCDF variable for resource ' + rsrce)
        return None

    vals_raw = metric_defn.nc_dset.variables[varname][:, lat_indx, lon_indx]

    if rsrce == 'fertiliser':
        vals = _check_array(vals_raw)

    elif rsrce == 'Ndeposition':

        # derive annual value from monthly values
        # =======================================
        vals = []
        for indx in range(0, len(vals_raw), 12):
            vals.append(sum(vals_raw[indx:indx + 12])/12.0)

        vals = [val*NDEP_FACTOR for val in vals]   #  change units from mg N m-2 yr-1 to Kg N ha-1 yr-1
        pass
    else:
        vals = [round(val_raw.item(), 5) for val_raw in vals_raw]

    # stretch values to cover simulation period
    # =========================================
    nyrs_add = metric_defn.start_year - sim_strt_yr
    if nyrs_add > 0:
        vals = nyrs_add * [vals[0]] + vals      # prepend first year of retrieved data
    elif nyrs_add < 0:
        nyrs_clip = abs(nyrs_add)
        vals = vals[nyrs_clip:]      # clip years from start of retrieved data

    nyrs_add = sim_end_yr - metric_defn.end_year
    if nyrs_add > 0:
        vals = vals + nyrs_add * [vals[-1]]
    elif nyrs_add < 0:
        vals = vals[:nyrs_add]      # clip years from end of retrieved data

    return vals

def identify_datasets(project_path, mask_fname):
    """
    called at startup - builds description of required NC files
    """
    descriptor = '\tNC files - mask file: '

    if isfile(mask_fname):
        descriptor += 'OK\t'
        print('\tmask file: ' + mask_fname)
    else:
        descriptor += 'none\t'
        print('\tmask file: none')

    # yields  NB yield is a Python reserved word
    # ==========================================
    resource = 'yields'
    descriptor += resource + ': '

    fns_path = join(project_path, 'Yield data')
    yield_fnames = glob(fns_path + '\\*.nc')
    if len(yield_fnames) == 0:
        descriptor += 'none\t'
        print('\t' + resource + ' file: none')
    else:
        descriptor += 'OK\t'
        print('\t' + resource + ' file: ' + yield_fnames[0])

    # fertiliser
    # ==========
    resource = 'fertiliser'
    descriptor += resource + ': '

    fert_path = join(project_path, 'Fertiliser')
    fert_fnames = glob(fert_path + '\\erawatch_eu_grm*.nc')
    if len(fert_fnames) == 0:
        descriptor += 'none\t'
        print('\t' + resource + ' file: none')
    else:
        descriptor += 'OK\t'
        print('\t' + resource + ' file: ' + fert_fnames[0])

    # N from livestock
    # ================
    resource = 'Nlivestock'
    descriptor += resource + ': '

    n_lvstck_path = join(project_path, 'N_from_livestock')
    n_lvstck_fnames = glob(n_lvstck_path + '\\n_available_livestock.nc')
    if len(n_lvstck_fnames) == 0:
        descriptor += 'none\t'
        print('\t' + resource + ' file: none')
    else:
        descriptor += 'OK\t'
        print('\t' + resource + ' file: ' + n_lvstck_fnames[0])

    # sowing_harvest
    # ==============
    resource = 'phenology'
    descriptor += resource + ': '

    pheno_path = join(project_path, 'Phenology')
    pheno_fnames = glob(pheno_path + '\\' + '*.nc')
    nfiles = len(pheno_fnames)
    if nfiles == 0:
        descriptor += 'none\t'
    else:
        descriptor += 'OK\t'
        print('\t' + resource + ' first of {} files: {}'.format(nfiles, pheno_fnames[0]))

    return descriptor

def create_proj_data_defns(project_path, mask_fname, crop_name, npp_yields_flag = True):
    """

    """
    if crop_name not in CROP_NAME_DICT:
        print(ERROR_STR + crop_name + ' not mappable to phenology NC name')
        return None
    crop_mapped = CROP_NAME_DICT[crop_name]

    resource = 'cropmasks'
    mask_defn = ManagementSet(mask_fname, resource)

    resource = 'Ndeposition'
    fns_path = join(project_path, resource)
    ndepos_fname = glob(fns_path + '\\*.nc')
    if len(ndepos_fname) == 0:
        print(ERROR_STR + fns_path + ': no N deposition file found for ' + crop_mapped)
        return None

    ndepos_defn = ManagementSet(ndepos_fname[0], resource, 'noy', units = 'N mg/m2/yr')

    # for now we just use the mean of the yields from 2000 to 2014;  NB yield is a Python reserved word
    # =================================================================================================
    resource = 'yields'
    fns_path = join(project_path, 'Yield data')
    if npp_yields_flag:
        yield_fname = glob(fns_path + '\\NPP_Grassland*.nc')
    else:
        yield_fname = glob(fns_path + '\\Chang*_Grass_Biomass_Use_*.nc')
    if len(yield_fname) == 0:
        print(ERROR_STR + fns_path + ': no yield file found for ' + crop_mapped)
        return None

    yield_defn = ManagementSet(yield_fname[0], resource, rqrd_varname = None, units = 't/ha')

    # fertiliser - 0.5 deg resolution
    # ================================
    resource = 'fertiliser'
    fert_path = join(project_path, 'Fertiliser')
    fert_fname = glob(fert_path + '\\Chang*_Fertilization_Rate_*.nc')
    if len(fert_fname) == 0:
        print(ERROR_STR + fert_path + ': no ' + resource + ' files found')
        return None

    fert_defn = ManagementSet(fert_fname[0], resource, rqrd_varname = None, units = 'N kg/ha')

    if len(fert_fname) == 0:
        print(ERROR_STR + fert_path + ': no ' + resource + ' files found')
        return None
    fert_fname = glob(fert_path + '\\erawatch_eu_grm_input_*.nc')
    fert_era_defn = ManagementSet(fert_fname[0], resource, rqrd_varname=None, units='N kg/ha')

    # N from livestock
    # ================
    resource = 'Nlivestock'
    n_lvstck_path = join(project_path, 'N_from_livestock')
    n_lvstck_fname = glob(n_lvstck_path + '\\n_available_livestock.nc')
    if len(n_lvstck_fname) == 0:
        print(ERROR_STR + fert_path + ': no ' + resource + ' files found')
        return None

    n_lvstck_defn = ManagementSet(n_lvstck_fname[0], resource, rqrd_varname=None, units='N kg/ha')

    # sowing harvest - always January and December
    # ============================================
    resource = 'phenology'
    if crop_name == 'Grassland':
        pheno_defn = None
    else:
        pheno_path = join(project_path, 'Phenology')
        pheno_fname = glob(pheno_path + '\\' + crop_mapped + '*.nc')
        if len(pheno_fname) == 0:
            print(ERROR_STR + pheno_path + ': no phenology file found')
            return None

        pheno_defn = ManagementSet(pheno_fname[0], resource)

    return mask_defn, yield_defn, pheno_defn, fert_defn, fert_era_defn, n_lvstck_defn, ndepos_defn

def open_proj_NC_sets(mask_defn, yield_defn, ndepos_defn, fert_defn, fert_era_defn, graze_defn):
    """
    """
    mask_defn.nc_dset  = Dataset(mask_defn.nc_fname, mode='r')
    yield_defn.nc_dset = Dataset(yield_defn.nc_fname, mode='r')
    ndepos_defn.nc_dset = Dataset(ndepos_defn.nc_fname, mode='r')
    fert_defn.nc_dset = Dataset(fert_defn.nc_fname, mode='r')
    fert_era_defn.nc_dset = Dataset(fert_era_defn.nc_fname, mode='r')
    graze_defn.nc_dset = Dataset(graze_defn.nc_fname, mode='r')

    return

def close_proj_NC_sets(mask_defn, yield_defn, ndepos_defn, fert_defn, fert_era_defn, graze_defn):
    """
    """
    mask_defn.nc_dset.close()
    yield_defn.nc_dset.close()
    ndepos_defn.nc_dset.close()
    fert_defn.nc_dset.close()
    fert_era_defn.nc_dset.close()
    graze_defn.nc_dset.close()

    return

class ManagementSet(object, ):
    """

    """
    def __init__(self, nc_fname, resource, rqrd_varname = None, units = None, cnvrsn_fctr = None):
        """

        """
        nc_fname = normpath(nc_fname)

        nc_dset = Dataset(nc_fname, mode='r')
        if 'lat' in nc_dset.variables:
            lat_var = 'lat'
            lon_var = 'lon'
        else:
            lat_var = 'latitude'
            lon_var = 'longitude'
        lats = nc_dset.variables[lat_var][:]
        lons = nc_dset.variables[lon_var][:]

        # record var names
        # ================
        exclude_vars = list([lat_var, lon_var, 'time'])
        start_year = None
        end_year   = None
        varnames  = []
        for var in nc_dset.variables:
            if var not in exclude_vars:
                varnames.append(var)

            if var == 'time':
                time_var = nc_dset.variables[var]
                if hasattr(time_var, 'units'):
                    time_units = time_var.units
                else:
                    time_units = 'year'

                if resource == 'yields' or resource == 'fertiliser':
                    '''
                    typical file - formerly: E:\\datasets\\Fertiliser\\Fertilization_Rate_Grassland_1860-2012.nc
                    typical file - currently: E:\\datasets\\Fertiliser\\erawatch_eu_grm_input_1901-2010.nc
                    '''
                    year_span = splitext(nc_fname)[0].split('_')[-1]
                    start_year, end_year = year_span.split('-')

                elif resource == 'Ndeposition':          # TODO: more elegant to use datetime i.e. days since 1860-1-1
                    fn_cmpnts = splitext(nc_fname)[0].split('_')
                    start_year = fn_cmpnts[-2]
                    end_year = fn_cmpnts[-1]
                else:
                    '''
                    since_time = time_units.split('since')[1]
                    start_year = int(since_time.split('-')[0]) + time_var[0]    # messy way to get to 1961
                    '''
                    start_year = 1981
                    end_year = start_year + len(time_var) - 1

        nc_dset.close()

        lat_frst = float(lats[0])
        lon_frst = float(lons[0])
        lat_last = float(lats[-1])
        lon_last = float(lons[-1])

        # required for bounding box
        # =========================
        if lat_last > lat_frst:
            lat_ll = lat_frst
            lat_ur = lat_last
        else:
            lat_ll = lat_last
            lat_ur = lat_frst

        if lon_last > lon_frst:
            lon_ll = lon_frst
            lon_ur = lon_last
        else:
            lon_ll = lon_last
            lon_ur = lon_frst

        self.lat_frst = float(lats[0])
        self.lon_frst = float(lons[0])
        self.lat_last = float(lats[-1])
        self.lon_last = float(lons[-1])

        self.lat_var = lat_var
        self.lon_var = lon_var
        self.bbox = lon_ll, lat_ll, lon_ur, lat_ur

        self.nc_fname = nc_fname
        self.varnames = varnames
        self.varname = rqrd_varname
        self.units = units
        self.cnvrsn_fctr = cnvrsn_fctr
        self.nc_dset = None

        # resolutions
        # ===========
        self.resol_lon = (lons[-1] - lons[0])/(len(lons) - 1)
        self.resol_lat = (lats[-1] - lats[0])/(len(lats) - 1)
        self.max_lat_indx = len(lats) - 1
        self.max_lon_indx = len(lons) - 1

        #
        self.lats = list(lats)
        self.lons = list(lons)

        if start_year is None:
            self.start_year = start_year
            self.end_year = end_year
        else:
            self.start_year = int(start_year)
            self.end_year   = int(end_year)

        self.rsrce = resource

    def get_nc_coords(self, latitude, longitude):

        ret_code = 'OK'

        lat_indx = int(round((latitude -  self.lat_frst)/self.resol_lat))
        lon_indx = int(round((longitude - self.lon_frst)/self.resol_lon))

        if lat_indx < 0 or lat_indx > self.max_lat_indx:
            ret_code = '*** Warning *** latitude index {} out of bounds for latitude {}\tmax indx: {}'.format(lat_indx,
                                                                                round(latitude, 4), self.max_lat_indx)
        if lon_indx < 0 or lon_indx > self.max_lon_indx:
            ret_code = '*** Warning *** longitude index {} out of bounds for longitude {}\tmax indx: {}'.format(lon_indx,
                                                                                round(longitude, 4), self.max_lon_indx)
        return lat_indx, lon_indx, ret_code
