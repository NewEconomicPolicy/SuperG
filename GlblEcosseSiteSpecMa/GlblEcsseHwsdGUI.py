#-------------------------------------------------------------------------------
# Name:
# Purpose:     Creates a GUI with five adminstrative levels plus country
# Author:      Mike Martin
# Created:     11/12/2015
# Licence:     <your licence>
#-------------------------------------------------------------------------------
#!/usr/bin/env python

__prog__ = 'GlblEcsseHwsdGUI.py'
__version__ = '0.0.1'
__author__ = 's03mm5'

import sys
from os import system

from PyQt5.QtCore import Qt
from PyQt5.QtGui import QPixmap
from PyQt5.QtWidgets import (QLabel, QWidget, QApplication, QHBoxLayout, QVBoxLayout, QGridLayout, QLineEdit,
                                            QComboBox, QRadioButton, QButtonGroup, QPushButton, QCheckBox, QFileDialog)

from shape_funcs import format_bbox, calculate_area
from commonCmpntsGUI import exit_clicked, commonSection, calculate_grid_cell, save_clicked
from glbl_ecss_cmmn_cmpntsGUI import grid_resolutions
from glbl_ecsse_high_level_fns import generate_banded_sims

import hwsd_mu_globals_fns
from initialise_funcs import read_config_file, initiation, change_config_file, write_config_file
from glbl_ecss_cmmn_funcs import build_and_display_studies, write_runsites_cnfg_file
from hilda_fns import generate_masks, write_countries_nc
from filter_hwsd_fns import filter_hwsd_csv, expand_phenology

GLBL_ECSSE_STR = 'global_ecosse_config_site_specific_'
STD_FLD_SIZE_60 = 60
STD_FLD_SIZE_80 = 80
STD_FLD_SIZE_100 = 100
STD_FLD_SIZE_160 = 160
NGRANULARITY = 120

from collections import namedtuple
RunCnfg = namedtuple("RunCnfg", "syn_fert_flag manure_flag graze_flag")
RUN_CNFGS = {'all': RunCnfg(True, True, True), 'nsyn_grz': RunCnfg(True, False, True),
             'nsyn_mnr': RunCnfg(True, True, False), 'grz_mnr': RunCnfg(False, True, True),
                                                                        'no_inpts': RunCnfg(False, False, False)}
class Form(QWidget):

    def __init__(self, parent=None):

        super(Form, self).__init__(parent)

        self.version = 'HWSD_grid'
        initiation(self)
        '''
        define two vertical boxes, in LH vertical box put the painter and in RH put the grid
        define horizon box to put LH and RH vertical boxes in
        '''
        hbox = QHBoxLayout()
        hbox.setSpacing(10)

        # left hand vertical box consists of png image
        # ============================================
        lh_vbox = QVBoxLayout()

        # LH vertical box contains image only
        lbl20 = QLabel()
        pixmap = QPixmap(self.sttngs['fname_png'])
        lbl20.setPixmap(pixmap)

        lh_vbox.addWidget(lbl20)

        # add LH vertical box to horizontal box
        hbox.addLayout(lh_vbox)

        # right hand box consists of combo boxes, labels and buttons
        # ==========================================================
        rh_vbox = QVBoxLayout()

        # The layout is done with the QGridLayout
        grid = QGridLayout()
        grid.setSpacing(10)	# set spacing between widgets

        # line 0
        # ======
        irow = 0
        lbl00 = QLabel('Study:')
        lbl00.setAlignment(Qt.AlignRight)
        grid.addWidget(lbl00, irow, 0)

        w_study = QLineEdit()
        helpText = 'Study name should not have spaces'        
        w_study.setToolTip(helpText)
        w_study.setFixedWidth(STD_FLD_SIZE_160)
        grid.addWidget(w_study, irow, 1)
        self.w_study = w_study

        lbl00s = QLabel('studies:')
        lbl00s.setAlignment(Qt.AlignRight)
        helpText = 'list of studies'
        lbl00s.setToolTip(helpText)
        grid.addWidget(lbl00s, irow, 2)

        combo00s = QComboBox()
        for study in self.studies:
            combo00s.addItem(str(study))
        grid.addWidget(combo00s, irow, 3)
        combo00s.currentIndexChanged[str].connect(self.changeConfigFile)
        self.combo00s = combo00s

        # lon/lats - first line
        # =====================
        irow += 1
        lbl02a = QLabel('Upper right longitude:')
        lbl02a.setAlignment(Qt.AlignRight)
        grid.addWidget(lbl02a, irow, 0)

        w_ur_lon = QLineEdit()
        w_ur_lon.setFixedWidth(STD_FLD_SIZE_60)
        w_ur_lon.textChanged[str].connect(self.bboxTextChanged)
        grid.addWidget(w_ur_lon, irow, 1)
        self.w_ur_lon = w_ur_lon

        lbl02b = QLabel('latitude:')
        lbl02b.setAlignment(Qt.AlignRight)
        grid.addWidget(lbl02b, irow, 2)

        w_ur_lat = QLineEdit()
        w_ur_lat.setFixedWidth(STD_FLD_SIZE_60)
        w_ur_lat.textChanged[str].connect(self.bboxTextChanged)
        grid.addWidget(w_ur_lat, irow, 3)
        self.w_ur_lat = w_ur_lat

        # second line
        # ===========
        irow += 1
        lbl01a = QLabel('Lower left longitude:')
        lbl01a.setAlignment(Qt.AlignRight)
        grid.addWidget(lbl01a, irow, 0)

        w_ll_lon = QLineEdit()
        w_ll_lon.setFixedWidth(STD_FLD_SIZE_60)
        w_ll_lon.textChanged[str].connect(self.bboxTextChanged)
        grid.addWidget(w_ll_lon, irow, 1)
        self.w_ll_lon = w_ll_lon

        lbl01b = QLabel('latitude:')
        lbl01b.setAlignment(Qt.AlignRight)
        grid.addWidget(lbl01b, irow, 2)

        w_ll_lat = QLineEdit()
        w_ll_lat.setFixedWidth(STD_FLD_SIZE_60)
        w_ll_lat.textChanged[str].connect(self.bboxTextChanged)
        grid.addWidget(w_ll_lat, irow, 3)
        self.w_ll_lat = w_ll_lat

        # report on AOI and spin-up mode
        # ==============================
        irow += 1
        lbl03a = QLabel('Study bounding box:')
        lbl03a.setAlignment(Qt.AlignRight)
        grid.addWidget(lbl03a, irow, 0)

        self.lbl03 = QLabel()
        grid.addWidget(self.lbl03, irow, 1, 1, 4)

        lbl03b = QLabel('Equilibrium mode:')
        lbl03b.setAlignment(Qt.AlignRight)
        helpText = 'mode of equilibrium run, generally OK with 6'
        lbl03b.setToolTip(helpText)
        grid.addWidget(lbl03b, irow, 5)

        w_equimode = QLineEdit()
        w_equimode.setFixedWidth(STD_FLD_SIZE_60)
        w_equimode.setText('')
        self.w_equimode = w_equimode
        grid.addWidget(w_equimode, irow, 6)

        # link RadioButtons with ButtonGroup
        # ==================================
        irow += 1
        w_lbl04 = QLabel('Timestep:')
        w_lbl04.setAlignment(Qt.AlignRight)
        grid.addWidget(w_lbl04, irow, 0)

        w_daily = QRadioButton("Daily")
        helpText = 'If this option is selected, then Daily timestep is used'
        w_daily.setToolTip(helpText)
        grid.addWidget(w_daily, irow, 1)
        self.w_daily  = w_daily

        w_mnthly = QRadioButton("Monthly")
        helpText = 'If this option is selected, then Monthly timestep is used'
        w_mnthly.setToolTip(helpText)
        w_mnthly.setChecked(True)
        grid.addWidget(w_mnthly, irow, 2)
        self.w_mnthly = w_mnthly

        # assign check values to radio buttons
        # ====================================
        w_timestep_choice = QButtonGroup()
        w_timestep_choice.addButton(w_daily)
        w_timestep_choice.setId(w_daily, 1)
        w_timestep_choice.addButton(w_mnthly)
        w_timestep_choice.setId(w_mnthly, 2)
        self.w_timestep_choice = w_timestep_choice

        # simplification options
        # ======================
        w_use_dom_soil = QCheckBox('Use most dominant soil')
        helpText = 'Each HWSD grid cell can have up to 10 soils. Select this option to use most dominant soil and\n' \
                ' discard all others. The the most dominant soil is defined as having the highest percentage coverage '\
                ' of all the soils for that grid cell'
        w_use_dom_soil.setToolTip(helpText)
        grid.addWidget(w_use_dom_soil, irow, 3)
        w_use_dom_soil.setEnabled(False)
        self.w_use_dom_soil = w_use_dom_soil

        w_use_high_cover = QCheckBox('Use highest coverage soil')
        helpText = 'Each meta-cell has one or more HWSD mu global keys with each key associated with a coverage expressed \n'\
                ' as a proportion of the area of the meta cell. Select this option to use the mu global with the highest coverage,\n' \
                ' discard the others and aggregate their coverages to the selected mu global'
        w_use_high_cover.setToolTip(helpText)
        grid.addWidget(w_use_high_cover, irow, 4, 1, 2)
        w_use_high_cover.setEnabled(False)
        self.w_use_high_cover = w_use_high_cover

        # mask and report only options
        # ============================
        irow += 1
        w_report = QCheckBox('Report data only')
        helpText = 'Only report grid cell inputs, namely Nsyn, Ngraz and manure amounts and lat/lon, to a CSV file'
        helpText += ' to the logs directory - no other processing is performed'
        w_report.setToolTip(helpText)
        grid.addWidget(w_report, irow, 0)
        self.w_report = w_report

        w_use_mask = QCheckBox('Use HILDA+ pasture mask')
        helpText = 'Use mask derived from HILDA+ Land use/land cover data vGLOB-1.0'
        w_use_mask.setToolTip(helpText)
        grid.addWidget(w_use_mask, irow, 3)
        w_use_mask.setEnabled(True)
        self.w_use_mask = w_use_mask

        # spacer if required
        # grid.addWidget(QLabel(''), 5, 0)

        # enable file of HWSD cells
        # =========================
        irow += 1
        w_use_csv_file = QPushButton("HWSD CSV file")
        helpText = 'Option to enable user to select a CSV file comprising latitude, longitiude and HWSD mu_global.'
        w_use_csv_file.setToolTip(helpText)
        grid.addWidget(w_use_csv_file, irow, 0)
        w_use_csv_file.clicked.connect(self.fetchCsvFile)

        w_lbl06 = QLabel('')  # label for HWSD csv file name
        grid.addWidget(w_lbl06, irow, 1, 1, 5)
        self.w_lbl06 = w_lbl06

        irow += 1
        w_lbl07a = QLabel('HWSD bounding box:')
        grid.addWidget(w_lbl07a, irow, 0)
        self.w_lbl07a = w_lbl07a

        self.w_lbl07 = QLabel('')     # label for HWSD AOI bounding box detail
        grid.addWidget(self.w_lbl07, irow, 1, 1, 5)

        # create weather rows
        # ===================
        irow = commonSection(self, grid, irow)

        # ===================
        irow += 1
        lbl05 = QLabel('Maximum cells:')
        lbl05.setAlignment(Qt.AlignRight)
        helpText = 'Maximum number of simulation cells to generate for each region'
        lbl05.setToolTip(helpText)
        grid.addWidget(lbl05, irow, 0)

        w_max_cells = QLineEdit()
        w_max_cells.setText('')
        w_max_cells.setFixedWidth(STD_FLD_SIZE_80)
        self.w_max_cells = w_max_cells
        grid.addWidget(w_max_cells, irow, 1)

        irow = grid_resolutions(self, grid, irow)    # grid resolution

        # =======
        w_nc_inp_data = QCheckBox("Create NC file of input data")
        helpText = 'If this option is selected, then an NC file of the input data will be created - weather and soil'
        w_nc_inp_data.setChecked(False)
        w_nc_inp_data.setEnabled(False)
        w_nc_inp_data.setToolTip(helpText)
        # grid.addWidget(w_nc_inp_data, 17, 4, 1, 2)

        # spacer then action buttons
        # ==========================
        irow += 1
        grid.addWidget(QLabel(), irow, 0)

        irow += 1
        w_create_files = QPushButton("Create sim files")
        helpText = 'Generate ECOSSE simulation file sets corresponding to ordered HWSD global mapping unit set in CSV file'
        w_create_files.setToolTip(helpText)
        # w_create_files.setEnabled(False)
        grid.addWidget(w_create_files, irow, 0)
        w_create_files.clicked.connect(self.createSimsClicked)
        self.w_create_files = w_create_files

        w_auto_run_ec  = QCheckBox('Auto run Ecosse')
        helpText = 'Option to automatically run the ECOSSE program after the completion of input file creation phase'
        w_auto_run_ec.setToolTip(helpText)
        grid.addWidget(w_auto_run_ec, irow, 1)
        self.w_auto_run_ec = w_auto_run_ec

        w_process_files = QPushButton("Run ECOSSE")
        helpText = 'run the ECOSSE program using the input files created for this study'
        w_process_files.setToolTip(helpText)
        w_process_files.setFixedWidth(STD_FLD_SIZE_100)
        grid.addWidget(w_process_files, irow, 2)
        w_process_files.clicked.connect(self.runEcosse)
        self.w_process_files = w_process_files

        icol = 4
        w_save = QPushButton("Save")
        helpText = 'Save configuration and study definition files'
        w_save.setToolTip(helpText)
        w_save.setFixedWidth(STD_FLD_SIZE_100)
        grid.addWidget(w_save, irow, icol)
        w_save.clicked.connect(self.saveClicked)

        icol += 1
        w_spec = QPushButton("Cancel")
        helpText = 'Leaves GUI without saving configuration and study definition files'
        w_spec.setToolTip(helpText)
        grid.addWidget(w_spec, irow, icol)
        w_spec.clicked.connect(self.cancelClicked)

        icol += 1
        w_exit = QPushButton("Exit", self)
        grid.addWidget(w_exit, irow, icol)
        w_exit.clicked.connect(self.exitClicked)

        # secondary actions
        # =================
        irow += 1
        icol = 3
        w_phnlgy = QPushButton("Expand phenology")
        helpText = 'read CSV file and add this phenology to existing NC file'
        w_phnlgy.setToolTip(helpText)
        w_phnlgy.setFixedWidth(STD_FLD_SIZE_100)
        w_phnlgy.setEnabled(False)
        grid.addWidget(w_phnlgy, irow, icol)
        w_phnlgy.clicked.connect(self.expndPhnlgyClicked)

        icol += 1
        w_mk_mask = QPushButton("Make masks")
        helpText = 'Make land use masks for pasture and croplands from HILDA Global Land Use Change datasets\n' + \
                   'NB a pre-existing land use mask file will be deleted'
        w_mk_mask.setToolTip(helpText)
        grid.addWidget(w_mk_mask, irow, icol)
        w_mk_mask.clicked.connect(self.makeMasksClicked)

        icol += 1
        w_csv_wthr = QPushButton("CSV weather")
        helpText = 'Generate monthly and daily CSV weather'
        w_csv_wthr.setToolTip(helpText)
        grid.addWidget(w_csv_wthr, irow, icol)
        w_csv_wthr.clicked.connect(self.genCsvWthrClicked)
        self.w_csv_wthr = w_csv_wthr

        icol += 1
        w_mk_cntries = QPushButton("Make countries")
        helpText = 'Make countries NC file from CSV file supplied by Astley'
        w_mk_cntries.setToolTip(helpText)
        grid.addWidget(w_mk_cntries, irow, icol)
        w_mk_cntries.clicked.connect(self.makeCntriesNcClicked)

        icol += 1
        w_fltr_hwsd = QPushButton("Filter hwsd")
        helpText = 'Step through user HWSD file and check lat/lon from each row against HILDA mask file'
        w_fltr_hwsd.setToolTip(helpText)
        grid.addWidget(w_fltr_hwsd, irow, icol)
        w_fltr_hwsd.clicked.connect(self.filterHwsdClicked)

        # Layout main window
        # ==================
        rh_vbox.addLayout(grid)     # add grid to RH vertical box
        hbox.addLayout(rh_vbox)     # vertical box goes into horizontal box
        self.setLayout(hbox)        # the horizontal box fits inside the window

        self.setGeometry(300, 300, 690, 250)    # posx, posy, width, height
        self.setWindowTitle('Global Ecosse Site Specific - generate sets of ECOSSE input files based on HWSD grid')

        # read and set values from last run
        # =================================
        read_config_file(self)

    def expndPhnlgyClicked(self):

        expand_phenology(self)

    def filterHwsdClicked(self):

        filter_hwsd_csv(self)

    def makeCntriesNcClicked(self):

        write_countries_nc(self)

    def makeMasksClicked(self):

        calculate_grid_cell(self)
        generate_masks(self)

    def genCsvWthrClicked(self):

        # generate_weather_only(self)
        pass

    def saveClicked(self):
        """

        """
        func_name =  __prog__ + ' saveClicked'

        # check for spaces
        # ================
        study = self.w_study.text()
        if study == '':
            print('study cannot be blank')
        else:
            if study.find(' ') >= 0:
                print('*** study name must not have spaces ***')
            else:
                save_clicked(self)
        return

    def fetchCsvFile(self):
        """
        QFileDialog returns a tuple for Python 3.5, 3.6
        """
        fname = self.w_lbl06.text()
        fname, dummy = QFileDialog.getOpenFileName(self, 'Open file', fname, 'CSV files (*.csv)')
        if fname != '':
            self.hwsd_mu_globals = hwsd_mu_globals_fns.HWSD_mu_globals_csv(self, fname)
            if self.success_flag:
                self.w_lbl06.setText(fname)
                self.mu_global_list = self.hwsd_mu_globals.mu_global_list
                self.w_lbl07.setText(self.hwsd_mu_globals.aoi_label)

        return

    def resolutionChanged(self):

        calculate_grid_cell(self, NGRANULARITY)

    def bboxTextChanged(self):

        try:
            bbox = list([float(self.w_ll_lon.text()), float(self.w_ll_lat.text()),
                float(self.w_ur_lon.text()), float(self.w_ur_lat.text())])
            area = calculate_area(bbox)
            self.lbl03.setText(format_bbox(bbox, area))
            self.bbox = bbox
        except ValueError as err:
            pass

    def createSimsClicked(self):

        func_name =  __prog__ + ' createSimsClicked'

        hist_start_year = int(self.combo09s.currentText())
        hist_end_year   = int(self.combo09e.currentText())
        if hist_start_year > hist_end_year:
            print('Historic end year must be greater or equal to the start year')
            return

        fut_start_year = int(self.combo11s.currentText())
        fut_end_year   = int(self.combo11e.currentText())
        if fut_start_year > fut_end_year:
            print('Simulation end year must be greater or equal to the start year')
            return

        # overides ECOSSE file creation
        # =============================
        study = self.w_study.text()
        if study == '':
            print('study cannot be blank')
        else:
            self.study = study
            ret_code = generate_banded_sims(self)

            # run ECOSSE program
            # ==================
            if ret_code == 0:
                if self.w_auto_run_ec.isChecked():
                    self.runEcosse()

    def runEcosse(self):
        '''
        run the make simulations script
        '''
        func_name =  __prog__ + ' runEcosse'

        #  make sure config settings are saved
        # ====================================
        write_config_file(self)
        for sim_nm in RUN_CNFGS:
            if not write_runsites_cnfg_file(self, sim_nm):
                return

            cmd_str = self.sttngs['python_exe'] + ' ' + self.sttngs['runsites_py'] + ' ' + self.sttngs['runsites_cnfg_fn']
            system(cmd_str)

        print('\n*** Completed processing {} simulations ***'.format(len(RUN_CNFGS)))

    def cancelClicked(self):
        '''
        C
        '''
        exit_clicked(self, write_config_flag = False)

    def exitClicked(self):
        '''
        exit cleanly
        '''
        exit_clicked(self)

    def changeConfigFile(self):
        '''
        permits change of configuration file
        '''
        change_config_file(self)

def main():
    """

    """
    app = QApplication(sys.argv)  # create QApplication object
    form = Form() # instantiate form
    # display the GUI and start the event loop if we're not running batch mode
    form.show()             # paint form
    sys.exit(app.exec_())   # start event loop

if __name__ == '__main__':
    main()
