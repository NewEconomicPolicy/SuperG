"""
#-------------------------------------------------------------------------------
# Name:        initialise_funcs.py
# Purpose:     script to read read and write the setup and configuration files
# Author:      Mike Martin
# Created:     31/07/2015
# Licence:     <your licence>
#-------------------------------------------------------------------------------
"""

__prog__ = 'initialise_funcs.py'
__version__ = '0.0.0'

# Version history
# ---------------
# 
from os import getcwd, remove, makedirs, mkdir, name as os_name
from os.path import exists, normpath, splitext, isfile, isdir, join, lexists

from json import dump as json_dump, load as json_load
from json import JSONDecodeError

from time import sleep
from sys import exit
from glob import glob

from set_up_logging import set_up_logging
from glbl_ecss_cmmn_funcs import read_crop_pars_codes, read_manure_types, check_sims_dir, check_runsites
from glbl_ecss_cmmn_cmpntsGUI import print_resource_locations
from hwsd_mu_globals_fns import HWSD_mu_globals_csv
from shape_funcs import format_bbox, calculate_area
from weather_datasets import read_weather_dsets_detail
from mngmnt_fns_and_class import identify_datasets

BBOX_DEFAULT = [116.90045, 28.2294, 117.0, 29.0]    # bounding box default - somewhere in SE Europe
sleepTime = 5

ERROR_STR = '*** Error *** '
GLBL_ECSSE_STR = 'global_ecosse_config_site_specific_'
RUNSITES_CONFIG_FN = 'global_ecosse_site_spec_runsites_config.json'
APPLIC_STR = 'glbl_ecss_site_spec_ma'

MIN_GUI_LIST = ['autoRunEcFlag', 'wthrRsrce', 'aveWthrFlag', 'bbox', 'fertiliserFname', 'sowingDatesFname',
                                                                        'dailyMode',  'maxCells', 'hwsdCsvFname']
CMN_GUI_LIST = ['study', 'cruStrtYr', 'cruEndYr', 'climScnr', 'futStrtYr', 'futEndYr',  'cropIndx',
                                                                        'manureIndx', 'gridResol', 'eqilMode']
REQUIRED_CROPS = ['Grassland', 'Winter Wheat', 'Spring Wheat', 'Grain Maize', 'Spring Oats', 'Spring Barley']

def initiation(form):
    """
    this function is called to initiate the programme to process non-GUI settings.
    """
    # retrieve settings
    form.sttngs = _read_setup_file(form)
    form.parms_settings = _read_site_specific_parms()

    # set years_from_flag to gui or 'sowing' to take simulation period from sowing dates file
    # =======================================================================================
    form.years_from_flag = 'gui'
    form.req_resol_deg = None

    # facilitate multiple config file choices
    # =======================================
    glbl_ecsse_str = form.sttngs['config_dir'] + '/' + GLBL_ECSSE_STR
    config_files = glob(glbl_ecsse_str + '*.json')
    studies = []
    for fname in config_files:
        dummy, remainder = fname.split(GLBL_ECSSE_STR)
        study, dummy = splitext(remainder)
        if study != '':
            studies.append(study)
    form.studies = studies

    if len(config_files) > 0:
        form.config_file = config_files[0]
    else:
        form.config_file = glbl_ecsse_str + 'dummy.json'

    set_up_logging(form, APPLIC_STR)

    # Nitpars, fnames.dat and model_switches must be present - crop_pars for limited data mode only
    # =============================================================================================
    form.dflt_ecss_fns = {}
    for ecosse_fname in list(['Model_Switches', 'fnames', 'CROP_SUN', 'Nitpars']):
        dflt_ecss_fn  = join(form.sttngs['ecosse_fpath'], ecosse_fname + '.dat')
        if isfile(dflt_ecss_fn):
            form.dflt_ecss_fns[ecosse_fname.lower()] = dflt_ecss_fn
        else:
            print('{} file must exist'.format(dflt_ecss_fn))
            sleep(sleepTime)
            exit(0)

    form.crop_codes = read_crop_pars_codes(form.dflt_ecss_fns['crop_sun'], REQUIRED_CROPS, nset_lines=10)

    manure_types_fn = join(form.sttngs['ecosse_fpath'], 'manure_types.xlsx')
    form.manure_types = read_manure_types(manure_types_fn)

    if not check_sims_dir(form.lgr, form.sttngs['sims_dir']):
        sleep(sleepTime)
        exit(0)

    return

def change_config_file(form, new_study = None):
    """
    identify and read the new configuration file
    """
    if new_study is None:
        new_study = form.combo00s.currentText()

    new_config = 'global_ecosse_config_site_specific_' + new_study
    config_file = form.sttngs['config_dir'] + '/' + new_config + '.json'

    if isfile(config_file):
        form.config_file = config_file
        read_config_file(form)
        form.study = new_study
        form.w_study.setText(new_study)
    else:
        print('Could not locate ' + config_file)

def _default_parms_settings ():

    _default_parms = {
          'site' : {
            'iom_c': 0.0,
            'toc': 2.7,
            'drain_class': 2,
            'depth_imperm_lyr': 3,
            'wtr_tbl_dpth': 300,
            'prev_lu': 1,
            'prev_crop_code': 1,
            'yield_prev_crop': 8.0,
            'prev_crop_harvest_doy': 0,
            'equil_mode': 2
          },
          'cultivation': {
            'cult_type': 3,
            'cult_vigor': 0.5
          }
        }
    return _default_parms

def _read_site_specific_parms():
    """
    read programme run settings from the parameters file, if it exists
    """
    func_name =  __prog__ +  '  _read_site_specific_parms'

    # look for setup file here...
    parms_setup_file = join(getcwd(), 'additional_setup', 'site_specific_parms.json')

    if exists(parms_setup_file):
        with open(parms_setup_file, 'r') as fsetup:
            try:
                parms_settings = json_load(fsetup)
            except (JSONDecodeError, OSError, IOError) as err:
                print(err)
    else:
        parms_settings = _default_parms_settings ()

    return parms_settings

def _read_setup_file(form):
    """
    read settings used for programme from the setup file, if it exists,
    or create setup file using default values if file does not exist
    """
    func_name =  __prog__ +  ' _read_setup_file'

    # validate setup file
    # ===================
    fname_setup = APPLIC_STR + '_setup.json'
    setup_file = join(getcwd(), fname_setup)

    if exists(setup_file):
        with open(setup_file, 'r') as fsetup:
            try:
                settings = json_load(fsetup)
            except (JSONDecodeError, OSError, IOError) as err:
                print(ERROR_STR + 'reading setup file ' + setup_file + '\n\t' + str(err))
                sleep(sleepTime)
                exit(0)
    else:
        settings = _write_default_setup_file(setup_file)
        print('Read setup file ' + setup_file)

    # validate setup file
    # ===================
    grp = 'setup'
    settings_list = ['config_dir', 'fname_png', 'hwsd_dir', 'log_dir', 'mask_fn', 'proj_path',
                                                'python_exe', 'runsites_py', 'sims_dir', 'weather_dir']
    for key in settings_list:
        if key not in settings[grp]:
            print(ERROR_STR + 'setting {} is required in setup file {} '.format(key, setup_file))
            sleep(sleepTime)
            exit(0)

    # initialise vars
    # ===============
    config_dir      = settings[grp]['config_dir']
    form.hwsd_dir   = settings[grp]['hwsd_dir']
    log_dir         = settings[grp]['log_dir']
    mask_fn         = settings[grp]['mask_fn']
    proj_path  = settings[grp]['proj_path']
    python_exe = settings[grp]['python_exe']
    runsites_py     = settings[grp]['runsites_py']
    sims_dir   = settings[grp]['sims_dir']
    weather_dir     = settings[grp]['weather_dir']

    print('\nProject path: ' + proj_path)

    # Hilda land use
    # ==============
    if isfile(mask_fn):
        print('Will use Hilda land use mask: ' + mask_fn)
    else:
        print(ERROR_STR + 'in setup file {}\n\tHilda land use mask {} must exist'.format(setup_file, mask_fn))
        sleep(sleepTime)
        exit(0)

    # make sure directories exist for configuration and log files
    # ===========================================================
    if not lexists(log_dir):
        makedirs(log_dir)

    if not lexists(config_dir):
        makedirs(config_dir)

    # check ECOSSE can be run
    # =======================
    check_runsites(setup_file, settings[grp], config_dir, RUNSITES_CONFIG_FN, python_exe, runsites_py)

    # weather is crucial
    # ==================
    if lexists(weather_dir):
        form.weather_dir = weather_dir
    else:
        print(ERROR_STR + 'reading {}\tClimate directory {} must exist'.format(setup_file, weather_dir))
        sleep(sleepTime)
        exit(0)

    # check weather data
    # ==================
    read_weather_dsets_detail(form)
    if len(form.weather_sets) == 0:
        print(ERROR_STR + 'no recognised weather sets in ' + weather_dir)
        sleep(sleepTime)
        exit(0)

    # location of Ecosse files e.g. fnames.dat
    # ========================================
    ecosse_fpath = join(proj_path, 'Ecosse_input_files')
    if not lexists(ecosse_fpath):
        print(ERROR_STR + 'Ecosse input files path {} must exist'.format(ecosse_fpath))
        sleep(sleepTime)
        exit(0)
    print('Ecosse input files path: ' + ecosse_fpath)
    settings[grp]['ecosse_fpath'] = ecosse_fpath

    # TODO: most of these are not used
    # ================================
    grp = 'run_settings'
    try:
        check_space_every = settings[grp]['check_space_every']
        completed_max = settings[grp]['completed_max']
        settings['setup']['kml_flag'] = settings[grp]['kml_flag']
        soilTestFlag = settings[grp]['soil_test_flag']
        space_remaining_limit = settings[grp]['space_remaining_limit']
        form.zeros_file   = settings[grp]['zeros_file']
    except KeyError:
        print(ERROR_STR + 'in group: ' + grp)
        sleep(sleepTime)
        exit(0)

    lta_nc_fname = None
    print_resource_locations(setup_file, config_dir, form.hwsd_dir, weather_dir, lta_nc_fname, sims_dir, log_dir)

    form.settings = {'log_dir': log_dir}
    return settings['setup']

def _write_default_setup_file(setup_file):
    """
    #  stanza if setup_file needs to be created
    """
    # Windows only for now
    # =====================
    if os_name != 'nt':
        print('Operating system is ' + os_name + 'should be nt - cannot proceed with writing default setup file')
        sleep(sleepTime)
        exit(0)

    # return list of drives
    # =====================
    import win32api

    drives = win32api.GetLogicalDriveStrings()
    drives = drives.split('\000')[:-1]
    if 'S:\\' in drives:
        root_dir_app  = 'S:\\tools\\'     # Global Ecosse installed here
        root_dir_user = 'H:\\'            # user files reside here
    else:
        root_dir_app  = 'E:\\'
        root_dir_user = 'C:\\AbUniv\\'

    suite_path = root_dir_app + 'GlobalEcosseSuite\\'
    data_path  = root_dir_app + 'GlobalEcosseData\\'
    outputs_path  = root_dir_app + 'GlobalEcosseOutputs\\'
    root_dir_user += 'GlobalEcosseSuite\\'
    runsites_py = ''

    _default_setup = {
        'setup': {
            'root_dir_user' : root_dir_user,
            'fname_png'     : join(suite_path + 'Images', 'Tree_of_life.PNG'),
            'images_dir'    : outputs_path + 'images',
            'python_exe'   : 'E:\\Python36\\python.exe',
            'runsites_py'  : runsites_py,
            'sims_dir'      : outputs_path + 'EcosseSims',
            'log_dir'       : root_dir_user + 'logs',
            'config_dir'    : root_dir_user + 'config',
            'hwsd_csv_fname': '',
            'hwsd_dir'   : data_path + 'HWSD_NEW',
            'weather_dir': 'E:\\GlobalEcosseData'
        },
        'run_settings': {
            'completed_max'         : 5000000000,
            'check_space_every'     : 10,
            'space_remaining_limit' : 1270,
            'kml_flag'      : True,
            'soil_test_flag': False,
            'zeros_file'    : False
        }
    }
    # create setup file
    # =================
    with open(setup_file, 'w') as fsetup:
        json_dump(_default_setup, fsetup, indent=2, sort_keys=True)
        fsetup.close()
        return _default_setup

def _write_default_config_file(config_file):
    """
    #        ll_lon,    ll_lat  ur_lon,ur_lat
    # stanza if config_file needs to be created
    """
    _default_config = {
        'minGUI': {
            'bbox': BBOX_DEFAULT,
            'wthrRsrce': 0,
            'aveWthrFlag': False,
            'hwsdCsvFname': '',
            'fertiliserFname': '',
            'daily_mode': True
        },
        'cmnGUI': {
            'climScnr' : 0,
            'cropIndx' : 0,
            'cruStrtYr': 0,
            'cruEndYr' : 0,
            'eqilMode' : 9.5,
            'futStrtYr': 0,
            'futEndYr' : 0,
            'gridResol': 0,
            'study'    : ''
        }
    }
    # if config file does not exist then create it...
    with open(config_file, 'w') as fconfig:
        json_dump(_default_config, fconfig, indent=2, sort_keys=True)
        fconfig.close()
        return _default_config

def read_config_file(form):
    """
    read widget settings used in the previous programme session from the config file, if it exists,
    or create config file using default settings if config file does not exist
    """
    func_name =  __prog__ +  ' read_config_file'
    config_file = form.config_file
    if exists(config_file):
        try:
            with open(config_file, 'r') as fconfig:
                config = json_load(fconfig)
                print('Read config file ' + config_file)
        except (OSError, IOError) as err:
                print(err)
                return False
    else:
        config = _write_default_config_file(config_file)

        print('Wrote configuration file ' + config_file)

    grp = 'minGUI'
    for key in MIN_GUI_LIST:
        if key not in config[grp]:
            print(ERROR_STR + 'attribute {} required for group {} in config file {}'.format(key, grp, config_file))
            sleep(sleepTime)
            exit(0)

    wthr_rsrce_indx = config[grp]['wthrRsrce']
    auto_run_ec = config[grp]['autoRunEcFlag']
    ave_weather = config[grp]['aveWthrFlag']
    form.bbox = config[grp]['bbox']
    fertiliser_fname = config[grp]['fertiliserFname']
    sowing_dates_fname = config[grp]['sowingDatesFname']
    max_cells = config[grp]['maxCells']
    daily_mode =  config[grp]['dailyMode']
    hwsd_csv_fname = config[grp]['hwsdCsvFname']

    # patch to permit backwards compatibility
    # =======================================
    if not isinstance(wthr_rsrce_indx, int):
        wthr_rsrce_indx = 0   # sets to CRU, the default

    # make sure index is within the permissable range of entries
    # ==========================================================
    nitems = form.combo10w.count()
    if wthr_rsrce_indx >= 0 and wthr_rsrce_indx < nitems:
        form.combo10w.setCurrentIndex(wthr_rsrce_indx)

    form.w_max_cells.setText(max_cells)

    # redundant - fertiliser and sowing dates file feedback
    # ======================================================
    descr = identify_datasets(form.sttngs['proj_path'], form.sttngs['mask_fn'])  # displays file info
    print(descr)

    if auto_run_ec:
        form.w_auto_run_ec.setCheckState(2)
    else:
        form.w_auto_run_ec.setCheckState(0)

    # common area
    # ===========
    grp = 'cmnGUI'
    for key in CMN_GUI_LIST:
        if key not in config[grp]:
            if key == 'manureIndx':
                config[grp]['manureIndx'] = 0
            else:
                print(ERROR_STR + 'attribute {} is required for group {} in config file {}'.format(key, grp, config_file))
                sleep(sleepTime)
                exit(0)

    form.w_study.setText(str(config[grp]['study']))
    form.combo09s.setCurrentIndex(config[grp]['cruStrtYr'])
    form.combo09e.setCurrentIndex(config[grp]['cruEndYr'])
    form.combo10.setCurrentIndex(config[grp]['climScnr'])
    form.combo11s.setCurrentIndex(config[grp]['futStrtYr'])
    form.combo11e.setCurrentIndex(config[grp]['futEndYr'])
    form.combo12.setCurrentIndex(config[grp]['cropIndx'])
    form.combo13.setCurrentIndex(config[grp]['manureIndx'])
    form.combo16.setCurrentIndex(config[grp]['gridResol'])
    form.w_equimode.setText(str(config[grp]['eqilMode']))

    # bounding box set up
    area = calculate_area(form.bbox)
    ll_lon, ll_lat, ur_lon, ur_lat = form.bbox
    form.w_ll_lon.setText(str(ll_lon))
    form.w_ll_lat.setText(str(ll_lat))
    form.w_ur_lon.setText(str(ur_lon))
    form.w_ur_lat.setText(str(ur_lat))
    form.lbl03.setText(format_bbox(form.bbox,area))

    # reset widgets associated with the HWSD file
    # ===========================================
    form.fstudy = ''
    if hwsd_csv_fname != '':
        form.w_lbl06.setText(hwsd_csv_fname)
        if isfile(hwsd_csv_fname):
            # read CSV file using pandas and create obj
            form.hwsd_mu_globals = HWSD_mu_globals_csv(form, hwsd_csv_fname)
            if form.success_flag:
                form.w_lbl07.setText(form.hwsd_mu_globals.aoi_label)
        else:
            print('HWSD csv file ' + hwsd_csv_fname + ' does not exist')
            hwsd_csv_fname = ''

    if hwsd_csv_fname == '':
        form.hwsd_mu_globals = None
        form.w_lbl06.setText('')
        form.w_lbl07.setText('')

    form.hwsd_csv_fname = hwsd_csv_fname

    # set check boxes
    # ===============
    if ave_weather:
        form.w_strt_1801.setCheckState(2)
    else:
        form.w_strt_1801.setCheckState(0)

    if daily_mode:
        form.w_daily.setChecked(True)
    else:
        form.w_mnthly.setChecked(True)

    if form.sttngs['runsites_cnfg_fn'] == None:
        form.w_auto_run_ec.setEnabled(False)
        form.w_process_files.setEnabled(False)

    # avoids errors when exiting
    # ==========================
    form.req_resol_deg = None
    form.req_resol_granul = None

    form.w_use_dom_soil.setChecked(True)
    form.w_use_high_cover.setChecked(True)
    form.w_use_mask.setChecked(True)

    DFLT_MAX_BAND  = 999

    start_at_band = 0
    stop_after_band = DFLT_MAX_BAND
    # stop_after_band = 8
    if start_at_band != 0 or stop_after_band < DFLT_MAX_BAND:
        print('***Warning*** will start simulations at band {} and stop after {}'.format(start_at_band, stop_after_band))
    form.start_at_band = start_at_band
    form.stop_after_band = stop_after_band

    return True

def write_config_file(form):
    """
    write current selections to config file
    """
    study = form.w_study.text()
    if study == '':
        study = form.combo00s.currentText()

    # facilitate multiple config file choices
    config_file = join(form.sttngs['config_dir'], GLBL_ECSSE_STR + study + '.json')

    # prepare the bounding box
    try:
        ll_lon = float(form.w_ll_lon.text())
        ll_lat = float(form.w_ll_lat.text())
        ur_lon = float(form.w_ur_lon.text())
        ur_lat = float(form.w_ur_lat.text())
    except ValueError:
        ll_lon = 0.0
        ll_lat = 0.0
        ur_lon = 0.0
        ur_lat = 0.0
    form.bbox =  list([ll_lon,ll_lat,ur_lon,ur_lat])

    # TODO:
    # print('Weather choice Id: {}'.format(form.w_weather_choice.checkedId()))
    config = {
        'minGUI': {
            'autoRunEcFlag': form.w_auto_run_ec.isChecked(),
            'bbox': form.bbox,
            'wthrRsrce': form.combo10w.currentIndex(),
            'aveWthrFlag': form.w_strt_1801.isChecked(),
            'hwsdCsvFname': normpath(form.w_lbl06.text()),
            'dailyMode': form.w_daily.isChecked(),
            'maxCells': form.w_max_cells.text(),
            'fertiliserFname' : None,
            'sowingDatesFname' : None
        },
        'cmnGUI': {
            'study'    : form.w_study.text(),
            'cruStrtYr': form.combo09s.currentIndex(),
            'cruEndYr' : form.combo09e.currentIndex(),
            'climScnr' : form.combo10.currentIndex(),
            'futStrtYr': form.combo11s.currentIndex(),
            'futEndYr' : form.combo11e.currentIndex(),
            'cropIndx' : form.combo12.currentIndex(),
            'manureIndx': form.combo13.currentIndex(),
            'eqilMode' : form.w_equimode.text(),
            'gridResol': form.combo16.currentIndex()
            }
        }
    if isfile(config_file):
        descriptor = 'Overwrote existing'
    else:
        descriptor = 'Wrote new'
    if study != '':
        with open(config_file, 'w') as fconfig:
            json_dump(config, fconfig, indent=2, sort_keys=True)
            print('\n' + descriptor + ' configuration file ' + config_file)
    return