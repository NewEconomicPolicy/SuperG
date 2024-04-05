#-------------------------------------------------------------------------------
# Name:        make_site_specific_files.py
# Purpose:     Generate ECOSSE site mode inpt files
# Author:      Mark Richards, University of Aberdeen
# Created:     30/05/2012
# Licence:     <your licence>
#-------------------------------------------------------------------------------
#!/usr/bin/env python


__version__ = '1.0.00'

from os.path import join, normpath
from numpy.ma.core import MaskedConstant

import validate as vldt

CALC_N_UPTK_INTERNAL = True

class MakeSiteFiles(object):
    _luts = ['ara', 'gra', 'for', 'nat', 'mis', 'src']
    months = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec']

    def __init__(self, form,  climgen, comments=True, spacer_len=12, nan=-999):
        self.comments = comments      # True = write comments, False = leave them out
        self.spacer_len = spacer_len  # Number of spaces between data and comment
        self.nan = nan

        #-----------------------------------------------------------------------
        # Input data
        #-----------------------------------------------------------------------
        # Model settings
        # Timestep (0=30 min, 1=daily, 2=weekly, 3=monthly)
        if climgen.fut_wthr_mnthly_flag:
            self.timestep = 3
        else:
            self.timestep = 1
        self.crop_model = 0             # 0=SUNDIAL, 1=MAGEC
        self.nyears = nan               # Number of years in simulation
        self.start_from_1801 = climgen.start_from_1801
        self.start_year = nan           # Year the simulation starts
        self.start_doy = 0              # Day of the year that the simulation starts on
        self.ntime_steps = nan           # Number of timesteps in simulation
        self.fixed_sim_end = 0          # Fixed end of simulation? (0=no, 1=yes)
        self.equil_mode = nan          # Mode of equilibrium run
        self.met_fnames = climgen.req_met_fnames    # Names of the met files

        # Site
        self.latitude = nan             # Latitude [decimal degrees]

        # Soil
        self.nlyrs = nan
        self.lyr_depths = []
        self.soil_lyrs = {}
        for lut in self._luts:
            self.soil_lyrs[lut] = []
        self.soil_name = 'xxxx'         # Name of soil
        self.soil_code = nan            # Soil code number
        self.drain_class = nan          # Drainage class (1=low, 2=moderate, 3=high) Not currently implemented.
        self.depth_imperm_lyr = nan     # Depth to impermeable layer (1=50cm, 2=100cm, 3=150cm
        self.date_reaches_fc = 1        # Date soil reaches field capacity (1=01/01; 2=01/06)
        self.wtr_tbl_dpth = nan         # Water Table depth [cm],  if > 150 cm than there is no effect
        self.wt_max_stand = 0.0         # Max standing water [cm]
        self.clay = nan                 # clay content [proportion]
        self.bulk_density = nan         # Bulk desnity [g/cm3]
        self.swc_wp = nan               # Soil water content at wilting point [mm/25cm]
        self.awc_fc = nan               # Available water content at field capacity [mm/25cm]
        self.awc_sat = nan              # Available water content at saturation [mm/25cm]
        self.biohum_stable_n2c = 0.118    # Stable N:C ratio of biomass and humus
        self.biohum_frac = 0.29          # Fraction of BIO + HUM formed: (biomass + humus) / total decomposition
        self.biohum_from_bio = 0.85      # Fraction of BIO + HUM produced from biomass decomposition
        self.biohum_from_hum = 0.85      # Fraction of BIO + HUM produced from humus decomposition
        self.bio_frac = 0.028             # Fraction of biomass in total organic C
        self.min_n_lvl = 0.0            # Minimum level of nitrate in soil [kg N/ha/50cm layer]
        self.k_bio = 0.635                # Rate constant for biomass decomposition [/year]
        self.k_hum = 0.02                # Rate constant for humus decomposition [/year]
        self.toc = nan                   # Total organic matter in top 50cm of soil [kg C/ha]
        self.iom_c = nan                  # Inert organic matter in top 50cm of soil [kg C/ha]
        self.soil_depth = 25.0           # Depth of soil used in initialisation [cm]
        self.ph = nan                   # pH of soil
        self.ph_zero_decomp = 2.5       # pH at which decomposition has declined to zero
        self.ph_decline_decomp = 5.5    # pH at which decomposition starts to decline

        # Climate
        self.lta_precip = None
        self.lta_tmean = None
        self.lta_pet = None

        # Land use
        self.prev_lu = nan                   # Land use before equilibrium run (1=arable, 2=grass, 3=forestry, 4=natural/scrub
        self.future_lu = []             # List of future land use codes (1 per year)

        # Crops & plant inputs
        self.prev_crop_code = nan       # Previous crop code
        self.yield_prev_crop = nan      # Yield of previous crop [t/ha]
        self.ncrops = 1                 # Number of crops
        self.crops = []                 # List of Crop objects
        self.cultivations = []          # List of Cultivation objects
        self.future_pi = []             # List of future annual plant inputs (1 element per year), 0 = estimate using MIAMI
        self.plant_inputs = {}
        for lut in self._luts:
            self.plant_inputs[lut] = 0
        self.prev_crop_harvest_doy = None

        # Inputs
        self.atmos_n_dep = 0.0

        # Not yet implemented
        self.c_accum_b4_change = 0.0
        self.ch4_b4_change = 0.0
        self.co2_b4_change = 0.0
        self.doc_loss_b4_change = 0.0

        self.del_soil_lyrs()

    def add_soil_lyr(self, lut_name, c_content, bulk_density, ph, clay_pc, silt_pc, sand_pc):
        self.soil_lyrs[lut_name.lower()].append(
                                SoilLyr(c_content, bulk_density, ph, clay_pc, silt_pc, sand_pc, self.nan))

    def del_soil_lyrs(self):
        for lut in self._luts:
            self.soil_lyrs[lut] = []

    def _line(self, data, comment):
        spacer_len = max(self.spacer_len - len(data), 2)
        spacer = ' ' * spacer_len
        return '{0}{1}# {2}\n'.format(data, spacer, comment)

    def validate(self):
        # Model settings
        assert(self.timestep in range(0,4))
        assert(self.crop_model in [0,1])
        assert(self.nyears > 0)
        assert(0 <= self.start_doy <= 366)
        assert(self.ntime_steps > 0)
        assert(self.fixed_sim_end in [0,1])
        assert(self.met_fnames != [])
        assert(0 < self.equil_mode < 10)

        # Site
        vldt.latitude([self.latitude])

        # Soil
        assert(self.soil_name != '')
        assert(self.soil_code >= 0)
        assert(self.drain_class in [1,2,3])
        assert(self.depth_imperm_lyr in [1,2,3])
        assert(self.date_reaches_fc in [1,2])
        assert(0 <= self.wtr_tbl_dpth <= 300)
        assert(0.0 <= self.clay <= 1.0)
        vldt.bulk_density([self.bulk_density])
        assert(self.wc_wp < self.awc_fc)
        assert(self.awc_fc < self.awc_sat)
        assert(10 <= self.wc_wp <= 100)  # mm/25cm
        assert(20 <= self.awc_fc <= 210)  # mm/25cm
        assert(30 <= self.awc_sat <= 230)  # mm/25cm
        assert(self.biohum_stable_n2c > 0.0)  # Stable N:C ratio of biomass and humus
        assert(self.biohum_frac < 1.0)        # Fraction of BIO + HUM formed: (biomass + humus) / total decomposition
        assert(self.biohum_from_bio < 1.0)    # Fraction of BIO + HUM produced from biomass decomposition
        assert(self.biohum_from_hum < 1.0)    # Fraction of BIO + HUM produced from humus decomposition
        assert(self.bio_frac < 1.0)           # Fraction of biomass in total organic C
        assert(self.min_n_lvl >= 0.0)         # Minimum level of nitrate in soil [kg N/ha/50cm layer]
        assert(self.k_bio > 0.0)              # Rate constant for biomass decomposition [/year]
        assert(self.k_hum > 0.0)              # Rate constant for humus decomposition [/year]
        assert(0.0 < self.toc < 500000)        # Total organic matter in top 50cm of soil [kg C/ha]
        assert(0.0 <= self.iom_c < 25000)       # Inert organic matter in top 50cm of soil [kg C/ha]
        assert(5 <= self.soil_depth <= 300)   # Depth of soil used in initialisation [cm]
        vldt.soil_ph([self.ph])                  # pH of soil
        vldt.soil_ph([self.ph_zero_decomp])       # pH at which decomposition has declined to zero
        assert(self.ph_decline_decomp > self.ph_zero_decomp)
        vldt.soil_ph([self.ph_decline_decomp])    # pH at which decomposition starts to decline
        # Soil layers
        assert(0 < self.nlyrs < 11)
        for i in range(len(self.lyr_depths)):
            assert(0 < self.lyr_depths[i] < 300)
            assert(self.lyr_depths[i] % 5 == 0)  # Depth should be a multiple of 5
            if i > 0:   # If not the top layer check that its depth is greater than the layer above
               assert(self.lyr_depths[i] > self.lyr_depths[i-1])
        for lut in self._luts:
            for lyr_num, lyr in enumerate(self.soil_lyrs[lut]):
                lyr.validate()
                assert(lyr_num < self.num_lyrs)

        # Land use
        assert(self.lu in range(1, len(self._luts) + 1))
        for lu in self.future_lu:
            assert(1 <= lu <= len(self._luts))
        assert(len(self.future_lu) == self.nyears)
        assert(len(self.future_pi) == self.nyears)

        # Climate'
        vldt.rainfall(self.lta_precip, 'daily')
        vldt.air_temp(self.lta_tmean)

        # Crop & plant inputs
        assert(self.prev_crop_code >= 0)       # Previous crop code
        assert(0.0 <= self.yield_prev_crop <= 15)      # Yield of previous crop [t/ha]
        assert(self.ncrops >= 0)               # Number of crops
        for crop in self.crops:
            crop.validate()
        if self.equil_mode in [1,3]:
            for key in self._luts:
                assert(0 <= self.plant_inputs[key] < 20000) # Upper limit?
        vldt.plant_c_inputs(self.future_pi, 'annual')

        # Inputs
        assert(0.0 <= self.atmos_n_dep <= 100.0)

    def write_sim_files(self, sim_dir, soil, latitude, hist_weather_recs, met_rel_path, mngmnt_content = None):
        """

        """
        if mngmnt_content is None:
            self._write_management_file(sim_dir, met_rel_path)
        else:
            with open(join(sim_dir, 'management.txt'), 'w') as fobj:
                fobj.writelines(mngmnt_content)

        self._write_soil_file(sim_dir)
        self._write_site_file(sim_dir)
        self.write_avemet_file(sim_dir, 'avemet.txt')
        self.write_avemet_file(sim_dir, 'AVEMET.DAT')
        write_crop_pars(sim_dir)
        write_model_switches(sim_dir)
        write_nitpars(sim_dir)

    def write_avemet_file(self, sim_dir, short_fname):
        """

        """
        with open(join(sim_dir, short_fname), 'w') as fobj:
            for imnth, (precip, pet, temp) in enumerate(
                    zip(self.lta_precip, self.lta_pet, self.lta_tmean)):
                fobj.write('{} {} {} {}\n'.format(imnth + 1, precip, pet, temp))

    def _write_fnames_file(self, sim_dir):

        dmodel = 1    # Denitrification model chosen: Bradbury, 1 or NEMIS, 2
        icmodel = 2   # Initialisation of SOM pools: fixed initialisation = 1, ICFIXED; RothC equilibrium run = 2, ICROTHCEQ
        inmodel = 1   # Initialisation of N: Bradbury's assumption of stable C:N ratio = 1, INSTABLECN; passed C:N ratio of DPM and RP = 2, INPASSCN
        docmodel = 2  # DOC model: DOC model on =1, DOC_ON; DOC model off = 2, DOC_OFF
        cnmodel = 1   # CN model: C:N ratio obtained by method of MAGEC = 1, CNMAGEC; C:N ratio obtained by method of Foereid = 2, CNFOEREID
        sparmodel = 2 # Soil parameter model: Soil parameters read in from soil.txt file = 1, SPARFILE; Soil parameters calculated from TOC = 2, SPARCALC

        # Type of equilibrium run:	EQNPP, EQTOC, EQNPPTOC, EQHILLIER, EQJONES /1,2,3,4,5/
        # =======================
        equil_mode = float(self.equil_mode)
        if equil_mode >= 9:
            eqmodel = 5
        elif equil_mode == 6:
            eqmodel = 4
        else:
            eqmodel = int(equil_mode)

        imfunc = 0    # Calculation of moisture rate modifiers: IMFUNC_ROTHC = 0, IMFUNC_HADLEY = 1
        itfunc = 0    # Calculation of temp.rate modifiers: ITFUNC_ROTHC = 0, ITFUNC_HADLEY = 1
        ch4model = 0  # CH4 model: CH4 model off = 0, CH4_OFF; Richards CH4 model on = 1, CH4_RICHARDS; Aitkenhead CH4 model on = 2, CH4_AITKENHEAD
        ec_eqrun = 0  # ECOSSE equilibrium run, (0 = off, 1 = on)

        with open(join(sim_dir, 'fnames.dat'), 'w') as fobj:
            fobj.write("'management.txt'  'soil.txt'  'site.txt'")
            fobj.write('\n{}		 DMODEL'.format(dmodel))
            fobj.write('\n{}		 ICMODEL'.format(icmodel))
            fobj.write('\n{}		 INMODEL'.format(inmodel))
            fobj.write('\n{}		 DOCMODEL'.format(docmodel))
            fobj.write('\n{}		 CNMODEL'.format(cnmodel))
            fobj.write('\n{}		 SPARMODEL'.format(sparmodel))
            fobj.write('\n{}		 EQMODEL'.format(eqmodel))
            fobj.write('\n{}		 IMFUNC'.format(imfunc))
            fobj.write('\n{}		 ITFUNC'.format(itfunc))
            fobj.write('\n{}		 CH4MODEL'.format(ch4model))
            fobj.write('\n{}		 EC_EQRUN'.format(ec_eqrun))
            fobj.flush()

    def _write_management_file(self, sim_dir, met_rel_path):
        """

        """
        _line = self._line  # Optimisation
        timestep = self.timestep
        lines = []
        lines.append(self._line('{}'.format(self.soil_code),
                                                    'Soil code number'))
        lines.append(self._line('{}'.format(self.drain_class),
                                                    'Drainage class (1=low, 2=moderate, 3=high) '
                                                    'Not currently implemented.'))
        lines.append(self._line('{}'.format(self.depth_imperm_lyr),
                                                    'Depth to impermeable layer (1=50cm, 2=100cm, 3=150cm)'))
        lines.append(self._line('{}'.format(self.prev_crop_code),
                                                    'Previous crop code'))
        lines.append(self._line('{}'.format(round(self.yield_prev_crop,3)),
                                                    'Yield of previous crop [t/ha]'))
        lines.append(self._line('{}'.format(round(self.atmos_n_dep, 3)),
                                                    'Atmospheric N deposition [kg N/ha]'))
        lines.append(self._line('{}'.format(self.date_reaches_fc),
                     'Date field reaches field capacity (1=01/01; 2=01/06)'))
        lines.append(self._line('{}'.format(self.timestep),
                     'Timestep (0=30 min, 1=daily, 2=weekly, 3=monthly)'))
        lines.append(self._line('{}'.format(self.crop_model),
                     'Crop model type (0=SUNDIAL, 1=MAGEC)'))
        lines.append(self._line('{}'.format(self.nyears),
                     'Number of years in simulation'))
        lines.append(self._line('{}'.format(int(self.prev_crop_harvest_doy - self.start_doy)),
                     'Timesteps from 01/01 to harvest of previous crop'))
        lines.append(self._line('{}'.format(self.start_year),
                     'First year of simulation'))
        lines.append(self._line('{}'.format(self.ntime_steps),
                     'End of simulation [number of timesteps]'))
        lines.append(self._line('{}'.format(self.fixed_sim_end),
                     'Fixed end of simulation? (0=no, 1=yes)'))
        lines.append(self._line('{}'.format(self.latitude),
                     'Latitude [decimal degrees]'))
        lines.append(self._line('{}'.format(self.wtr_tbl_dpth),
                     'Water table depth [cm], if > 150 cm there is no '
                     'effect'))
        lines.append(self._line('{}'.format(self.wt_max_stand),'Max standing water [cm]'))

        for iyr, fname in enumerate(self.met_fnames):
            lines.append(self._line("'{}{}'".format(met_rel_path, fname), 'Met file for year {}'.format(iyr+1)))
        lines.append(self._line('{}'.format(self.ncrops), 'Number of crops'))

        # stanza to permit timing of crops and fertiliser application
        # ===========================================================
        for cropnum, crop in enumerate(self.crops):
            lines.append(self._line('{}'.format(crop.code),
                         'CROP {}\t\tsequence {}'.format(crop.crop_name, cropnum + 1)))
            lines.append(self._line('{}'.format(int(crop.sowing_doy)),
                                                                        'Timesteps to sowing date from 01/01/01'))
            # trap error
            # ==========
            if CALC_N_UPTK_INTERNAL:
                n_uptake = 0
            else:
                n_uptake = crop.n_uptake
                if type(n_uptake) is MaskedConstant:
                    n_uptake = n_uptake.item()

            lines.append(self._line('{}'.format(round(n_uptake,2)),
                                                    'Crop N uptake at harvest (0=calculate internally) [kg N/ha]'))
            lines.append(self._line('{}'.format(int(crop.harvest_doy)),
                                                    'Timesteps to harvest date from 01/01/01'))
            lines.append(self._line('{}'.format(round(crop.exp_yield, 3)),
                                                    'Expected yield [t/ha]'))
            lines.append(self._line('{}'.format(crop.residues_inc),
                         'Crop residues incorporated (0=No, 1=Yes'))
            lines.append(self._line('{}'.format(crop.nfert_apps),
                         'Number of fertiliser applications'))
            lines.append(self._line('{}'.format(crop.nmanure_apps),
                         'Number of organic manure applications'))

            # stanza for fertiliser applications
            # ==================================
            if crop.nfert_apps >= 1:
                for fert in crop.fert_apps:
                    lines.append(_line('{}'.format(round(fert.amount,2)),
                        'Amount of fertiliser applied [kg N/ha]'))
                    lines.append(_line('{}'.format(int(fert.app_doy)),
                        'Timesteps to fertiliser application'))
                    lines.append(_line('{}'.format(fert.no3_pc),
                        'Percentage NO3'))
                    lines.append(_line('{}'.format(fert.nh4_pc),
                        'Percentage NH4'))
                    lines.append(_line('{}'.format(fert.urea_pc),
                        'Percentage urea'))
                    lines.append(_line('{}'.format(fert.non_amm_sulphate_salts),
                        'Does fert.contain ammonium salts other than ammonium '
                        'sulphate (0=No, 1=Yes)'))
                    lines.append(_line('{}'.format(fert.labelled),
                        'Has fertiliser been labelled (0=No, 1=Yes)'))

            # stanza for manure applications
            # ==============================
            if crop.nmanure_apps >= 1:
                for manure in crop.manure_apps:
                    lines.append(_line('{}'.format(round(manure.amount,2)),
                        'Amount of manure applied [kg N/ha]'))
                    lines.append(_line('{}'.format(int(manure.app_doy)),
                        'Timesteps to manure application'))
                    lines.append(_line('{}'.format(manure.type),
                        'Type of manure'))
                    lines.append(_line('{}'.format(manure.labelled),
                        'Has manure been labelled (0=No, 1=Yes)'))

        # stanza for cultivations
        # =======================
        lines.append(self._line('{}'.format(self.ncultivations), 'Number of cultivations'))
        for icult, cult in enumerate(self.cultivations):
            lines.append(_line('{}'.format(int(cult.cult_doy)),
                        'Timesteps from 01/01 when cultivation occurred'))
            lines.append(_line('{}'.format(cult.cult_type),
                        'Type of cultivation'))
            lines.append(_line('{}'.format(cult.cult_vigor),
                        'Vigour of cultivation'))

        with open(join(sim_dir, 'management.txt'), 'w') as fhand:
            fhand.writelines(lines)

    def _write_site_file(self, sim_dir, site_filename='site.txt'):
        _line = self._line
        output = []
        # Soil parameters
        output.append(_line('{0}'.format(self.equil_mode),
                      'Mode of equilibrium run'))
        output.append(_line('{0}'.format(self.nlyrs),
                      'Number of soil layers (max 10)'))
        for lyr_num, lyr_depth in enumerate(self.lyr_depths):
            output.append(_line('{0}'.format(lyr_depth), 'Depth of bottom of '
                          'SOM layer {} [cm]'.format(lyr_num+1)))
        for key in self._luts:
            for lyr_num in range(self.nlyrs):
                output.append(_line('{}'.format(self.soil_lyrs[key][lyr_num].soc),
                              'C content [kgC/ha] for this soil under {0} in '
                              'SOM layer {1}'.format(key, lyr_num + 1)))

                bulk_dens = round(float(self.soil_lyrs[key][lyr_num].bulk_dens), 3)
                output.append(_line('{}'.format(bulk_dens),
                      'Bulk density [g/cm3] for this soil under {} in SOM '
                      'layer {}'.format(key, lyr_num+1)))

                ph = round(float(self.soil_lyrs[key][lyr_num].ph), 1)
                output.append(_line('{}'.format(ph),
                      'pH for this soil under {} in SOM layer {}'
                      .format(key, lyr_num+1)))

                clay_pc = round(float(self.soil_lyrs[key][lyr_num].clay_pc), 1)
                output.append(_line('{}'.format(clay_pc),
                      '% clay by weight for this soil under {} in SOM layer {}'
                      .format(key, lyr_num+1)))

                silt_pc = round(float(self.soil_lyrs[key][lyr_num].silt_pc), 1)
                output.append(_line('{}'.format(silt_pc),
                      '% silt by weight for this soil under {} in SOM layer {}'
                      .format(key, lyr_num+1)))

                sand_pc = round(float(self.soil_lyrs[key][lyr_num].sand_pc), 1)
                output.append(_line('{}'.format(sand_pc),
                      '% sand by weight for this soil under {} in SOM layer {}'
                      .format(key, lyr_num+1)))
        for key in self._luts:
            output.append(_line('{0}'.format(self.plant_inputs[key]),
                          '{} long term average plant C input [kgC/ha/yr] '
                          '(used in modes 1 & 3 only)'.format(key)))

        # Long term average climate
        for precip, month in zip(self.lta_precip, self.months):
            output.append(_line('{0}'.format(precip),
                          '{} long term average monthly precipitation [mm]'
                          .format(month)))
        for tmean, month in zip(self.lta_tmean, self.months):
            output.append(_line('{0}'.format(tmean),
                          '{} long term average monthly temperature [mm]'
                          .format(month)))

        # Other bits and bobs
        output.append(_line('{0}'.format(self.latitude),
                      'Latitude [decimal deg]'))
        output.append(_line('{0}'.format(self.wtr_tbl_dpth),
                      'Water table depth at start [cm]'))
        output.append(_line('{0}'.format(self.drain_class), 'Drainage class'))
        output.append(_line('{0}'.format(self.c_accum_b4_change),
                     'C accumulated before change [kgC/ha/yr] (only for mode '
                     '4 - if not use a dummy a value)'))
        output.append(_line('{0}'.format(self.ch4_b4_change),
                      'CH4 emission before change [kgC/ha/yr] (not used yet)'))
        output.append(_line('{0}'.format(self.co2_b4_change),
                      'CO2 emission before change [kgC/ha/yr] (not used yet)'))
        output.append(_line('{0}'.format(self.doc_loss_b4_change),
                      'DOC loss before change [kgC/ha/yr] (not used yet)'))
        output.append(_line('{0}'.format(self.nyears),
                      'Number of growing seasons to simulate'))

        # Future land use and plant inputs
        yr_num = 1
        # for lu, pi in zip(self.future_lu, self.future_pi):
        lu = 1
        pi = 1000.0*self.yield_prev_crop
        pi = 0.0
        for iyear in range(self.nyears):
            output.append(_line('{}, {}'.format(lu, pi),
                         'Year {} land use and plant C input [kgC/ha/yr] (Not '
                         'used in mode 7. If plant input set to zero it is '
                         'obtained from RothC instead)'.format(yr_num)))
            yr_num += 1

        # Climate file names
        for year_num, fname in enumerate(self.met_fnames):
            output.append(_line('{}'.format(fname), 'Year {} climate file'.format(year_num + 1)))

        path = join(normpath(sim_dir), site_filename)
        with open(path, 'w', newline='') as fobj:
            fobj.writelines(output)


    def _write_soil_file(self, sim_dir):

        lines = []
        lines.append(self._line('{}'.format(self.soil_name), 'Soil name'))
        lines.append(self._line('{}'.format(self.soil_code), 'Soil code (used '
                                'in management file)'))
        lines.append(self._line('{}'.format(self.awc_fc), 'Available water '
                                'content at field capacity [mm/0-25cm]'))
        lines.append(self._line('{}'.format(self.biohum_stable_n2c), 'Stable '
                                'N:C of biomass and humus pools'))
        lines.append(self._line('{}'.format(self.biohum_frac), 'Fraction of '
                                'BIO + HUM formed: (BIO + HUM) / Total '
                                'decomposition'))
        lines.append(self._line('{}'.format(self.biohum_from_bio), 'Biomass/'
                                'Humus produced from biomass decomposition'))
        lines.append(self._line('{}'.format(self.biohum_from_hum), 'Biomass/'
                                'Humus produced from humus decomposition'))
        lines.append(self._line('{}'.format(self.bio_frac), 'Fraction of '
                                'biomass in total organic C'))
        lines.append(self._line('{}'.format(self.min_n_lvl), 'Minimum level of '
                                'nitrate in soil [kgN/ha/50cm layer]'))
        lines.append(self._line('{}'.format(self.k_bio), 'Rate constant for '
                                'biomass decomposition [/year]'))
        lines.append(self._line('{}'.format(self.k_hum), 'Rate constant for '
                                'humus decomposition [/year]'))
        lines.append(self._line('{}'.format(self.soc), 'Total organic C in '
                                'top 50 cm of soil [kgC/ha]'))
        lines.append(self._line('{}'.format(self.iom_c), 'Inert organic matter'
                                'in top 50 cm of soil [kgC/ha]'))
        lines.append(self._line('{}'.format(self.prev_lu), 'Land use before '
                                'equilibrium run (1=arable, 2=grass, '
                                '3=forestry, 4=natural'))
        lines.append(self._line('{}'.format(round(float(self.clay),1)), 'Clay content '
                                '[proportion]'))
        lines.append(self._line('{}'.format(self.soil_depth), 'Depth of soil '
                                'used in initialisation [cm]'))
        lines.append(self._line('{}'.format(round(float(self.ph),1)), 'pH of soil '))
        lines.append(self._line('{}'.format(self.ph_zero_decomp), 'pH at '
                                'which decomposition has declined to zero'))
        lines.append(self._line('{}'.format(self.ph_decline_decomp), 'pH at '
                                'which decomposition starts to decline'))
        lines.append(self._line('{}'.format(round(float(self.bulk_dens),1)), 'bulk denisty '
                                '[g/cm3]'))
        lines.append(self._line('{}'.format(self.awc_sat), 'Available water '
                                'content at saturation [mm/0-25cm]'))
        lines.append(self._line('{}'.format(self.wc_wp), 'Water content at '
                                'wilting point [mm/0-25cm]'))

        with open(join(sim_dir, 'soil.txt'), 'w') as fobj:
            fobj.writelines(lines)

class Crop(object):
    def __init__(self, nan=-999):
        self.code = nan
        self.sowing_doy = nan               # Day of the year when crop is sown
        self.harvest_doy = nan              # Day of the year when crop is harvested
        self.n_uptake = 0                   # Crop N uptake at harvest (0=calculate internally) [kg N/ha]
        self.exp_yield = nan                # Expected yield [t/ha]
        self.residues_inc = nan             # Crop residues incorporated (0=No, 1=Yes)
        self.nfert_apps = nan               # Number of fertiliser applications
        self.fert_apps = []
        self.nmanure_apps = nan             # Number of organic manure applications
        self.manure_apps = []
        # self.ngrz_n_apps = nan  # Number of N grazing applications
        # self.grz_n_apps = []

    def validate(self):
        assert(self.code >= 0)
        assert(0 < self.sowing_doy <= 366)   # Day of the year when crop is sown
        assert(0 < self.harvest_doy <= 366)  # Day of the year when crop is harvested
        assert(0 <= self.n_uptake <= 1000)   # Crop N uptake at harvest (0=calculate internally) [kg N/ha]
        assert(0.0 < self.exp_yield <= 25.0) # Expected yield [t/ha]
        assert(self.residues_inc in [0,1])   # Crop residues incorporated (0=No, 1=Yes)
        assert(0 <= self.nfert_apps < 10)    # Number of fertiliser applications
        assert(0 <= self.nmanure_apps < 10)  # Number of organic manure applications

class ManureApplication(object):
    def __init__(self, amount, app_doy, type = 1, labelled = 0):

        self.amount = amount
        self.app_doy = app_doy
        self.type = type
        self.labelled = labelled

class FertiliserApplication(object):
    def __init__(self, amount, app_doy, no3_pc, nh4_pc, urea_pc, non_amm_sulphate_salts=0, labelled=0):

        self.amount = amount
        self.app_doy = app_doy
        self.no3_pc = no3_pc
        self.nh4_pc = nh4_pc
        self.urea_pc = urea_pc
        self.non_amm_sulphate_salts = non_amm_sulphate_salts
        self.labelled = labelled

    def validate(self):
        assert(self.amount >= 0.0 and self.amount < 1000.0)
        assert(self.app_doy > 0 and self.app_doy < 367)
        assert(self.no3_pc >= 0.0 and self.no3_pc <= 100.0)
        assert(self.nh4_pc >= 0.0 and self.nh4_pc <= 100.0)
        assert(self.urea_pc >= 0.0 and self.urea_pc <= 100.0)
        assert(self.no3_pc + self.nh4_pc + self.urea_pc <= 100.0)
        assert(self.non_amm_sulphate_salts in [0, 1])
        assert(self.labelled in [0, 1])

class Cultivation(object):
    def __init__(self, cult_doy, cult_type, cult_vigor):
        self.cult_doy   = cult_doy             # Day of the year when cultivation occurred
        self.cult_type  = cult_type            # Type of cultivation � see table B1.1.3 in user manual
        self.cult_vigor = cult_vigor           # Vigour of cultivation (0.0 � 1.0)

    def validate(self):
        assert(0 < self.cult_doy <= 366)      # Day of the year when cultivation occurred
        assert(0 < self.cult_type <= 3)       # Type of cultivation � see table B1.1.3 in user manual
        assert(0.0 <= self.cult_vigor <= 1.0) # Vigour of cultivation (0.0 � 1.0)


class SoilLyr(object):
    def __init__(self, soc, bulk_dens, ph, clay_pc, silt_pc, sand_pc, no_data=-999):
        self.bulk_dens = bulk_dens
        self.ph = ph
        self.soc = soc     # C content
        self.clay_pc = clay_pc
        self.silt_pc = silt_pc
        self.sand_pc = sand_pc
        self.no_data = no_data

    def validate(self):
        if self.c != self.no_data: vldt.total_soil_carbon([self.c], depth='1 m')
        if self.bulk_dens != self.no_data: vldt.bulk_density([self.bulk_dens])
        if self.ph != self.no_data: vldt.soil_ph([self.ph])
        if self.clay_pc != self.no_data: vldt.percent([self.clay_pc])
        if self.silt_pc != self.no_data: vldt.percent([self.silt_pc])
        if self.sand_pc != self.no_data: vldt.percent([self.sand_pc])
        total = 0
        for val in [self.clay_pc, self.silt_pc, self.sand_pc]:
            if val != self.no_data:
                total += val
        vldt.percent([total])

def main():
    pass

if __name__ == '__main__':
    main()

def write_crop_pars(sim_dir):
    strng = """25,100,188,100,94,188,80,50,188                    CtoN       CROP_GIS_LIM   Ara,Gra,For,Semi_Nat,Misc,SRC,OSR,Beet,SRF
2,1,1,1,1,1,2,2,1                                  Sowmonth   CROP_GIS_LIM   Ara,Gra,For,Semi_Nat,Misc,SRC,OSR,Beet,SRF
7,12,12,12,12,12,9,9,12                            Harvmonth  CROP_GIS_LIM   Ara,Gra,For,Semi_Nat,Misc,SRC,OSR,Beet,SRF
4,4,0,0,5,0,4,4,4                                  Fertmonth  CROP_GIS_LIM   Ara,Gra,For,Semi_Nat,Misc,SRC,OSR,Beet,SRF
7,11,11,11,11,11,11,11,11                          Anthmonth  CROP_GIS_LIM   Ara,Gra,For,Semi_Nat,Misc,SRC,OSR,Beet,SRF

1,1,20,17,6,5,0,0,0                               Years to equilibrium      RUN1_GIS_CROP   Ara,Gra,For,Semi_Nat,Misc,SRC,OSR,Beet,SRF
0.0,0.1,0.4,0.2,0.0,0.5,0,0,0                     Fraction change to PI after land use   RUN1_GIS_CROP   Ara,Gra,For,Semi_Nat,Misc,SRC,OSR,Beet,SRF
0.29,1.00,1.00,1.00,0.67,1.00,0,0,0               Prop.of total N req.that goes below ground before harvest  RUN1_GIS_CROP   Ara,Gra,For,Semi_Nat,Misc,SRC,OSR,Beet,SRF
1.00,1.00,1.00,1.00,1.00,1.00,0,0,0               Ratio of C that goes below ground as debris before harvest to total C in above group crop at harvest   Ara,Gra,For,Semi_Nat,Misc,SRC,OSR,Beet,SRF
0.1,0.1,0.1,0.1,0.1,0.1,0,0,0                     Rate factor for N in non-cartable deb  RUN1_GIS_CROP   Ara,Gra,For,Semi_Nat,Misc,SRC,OSR,Beet,SRF
0.15,0.15,0.15,0.15,0.15,0.15,0,0,0               Rate factor for C in non-cartable deb  RUN1_GIS_CROP   Ara,Gra,For,Semi_Nat,Misc,SRC,OSR,Beet,SRF
0.085,0.0,0.0,0.0,0.0,0.0,0,0,0		              Proportion stubble in total N requirement  RUN1_GIS_CROP   Ara,Gra,For,Semi_Nat,Misc,SRC,OSR,Beet,SRF
0.0,0.0,0.0,0.0,0.0,0.0,0,0,1                     Ratio of C that goes below ground as stubble to total C in above group crop at harvest  RUN1_GIS_CROP   Ara,Gra,For,Semi_Nat,Misc,SRC,OSR,Beet,SRF

0.003,0.003,0.003,0.003,0.003,0.003,0.003,0.003,0.003    Crop uptake parameter RUN2_GIS_CROP   Ara,Gra,For,Semi_Nat,Misc,SRC,OSR,Beet,SRF
150,50,300,50,50,300,150,50,300                          Maximum rooting depth RUN2_GIS_CROP   Ara,Gra,For,Semi_Nat,Misc,SRC,OSR,Beet,SRF
5,20,15,15,20,15,5,5,5                                   Weeks from anthesis to end of growing season RUN2_GIS_CROP   Ara,Gra,For,Semi_Nat,Misc,SRC,OSR,Beet,SRF
5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0                      Rate of root growth (cm/week) RUN2_GIS_CROP   Ara,Gra,For,Semi_Nat,Misc,SRC,OSR,Beet,SRF
25.0,25.0,25.0,25.0,25.0,25.0,25.0,25.0,25.0             Root exploration factor RUN2_GIS_CROP   Ara,Gra,For,Semi_Nat,Misc,SRC,OSR,Beet,SRF
4,0,0,0,0,0,4,4,0                                        Seed N (arable only- others LUs perenniel) RUN2_GIS_CROP   Ara,Gra,For,Semi_Nat,Misc,SRC,OSR,Beet,SRF
2.5,2.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5                      Crop uptake parameter RUN2_GIS_CROP   Ara,Gra,For,Semi_Nat,Misc,SRC,OSR,Beet,SRF
"""
    with open(join(sim_dir, 'crop_pars.dat'), 'w', newline='') as fobj:
        fobj.write(strng)


def write_model_switches(sim_dir):
    strng = """2        Denitrification model chosen (Bradbury, 1 or NEMIS, 2)
0    	 Crop model type: 0=SUNDIAL, 1=MAGEC
1        Soil parameter model (from file, 1 or calc, 2)
0		 Metahne Model 0 = off 1 = on
2 	     DOC model 1 = off, 2 = on
0        N limitation spin-up used? 0 = No, 1 =Yes
0        Modify plant inputs according to full run using long term average weather data
0        Modify decomposition according to full run using long term average weather data
0        Choice of moisture rate modifier (0=ROTHC or 1=HADLEY)
0        Choice of temperature rate modifier (0=ROTHC or 1=HADLEY)
0        Use full equilibrium run of ECOSSE to initialise or not (0 = off, 1 = on)
2    	 Initisalisation of N by assuming steady state (after Bradbury = 1) or initialisation of N by passing the C:N ratio of DPM and RPM (=2)
1      	 0 = pH set to neutral or passed from parameter file, 1 = pH is read in from input file and does not change,  2 = pH is calculated using VSD (a very simple dynamic version of the MAGIC model by Ed Rowe & Chris Evans, CEH, Bangor)
0        Output full 0 = false, 1 = true
1        Output summary 0 = false, 1 = true
1        Plant Input from input or estimated from Total Organic Carbon (PI input = 1, PI est from TOC = 2)
"""
    with open(join(sim_dir, 'Model_Switches.dat'), 'w', newline='') as fobj:
        fobj.write(strng)


def write_nitpars(sim_dir):
    strng = """50            satnit        ! Max Conc of NH4 to be processed in a time step
0.02          PE1           ! Prop.gas.loss by nitrificn
0.4           PE3           ! Prop.gas.loss by nitrificn that is NO
0.02          FCGAS			! Prop.N2O loss by partial nitrificn at F
0.6           KNITRIF		! Rate constant for nitrification (/week)
16.5          CD1			! Amount of nitrate at which denitrification is 50% of its full potential (kg N / ha)
0.62          CD2           ! Fitted constant describing water modifier for denitrification
0.38          CD3			! Fitted constant describing water modifier for denitrification
1.74          CD4			! Fitted constant describing water modifier for denitrification
0.55          PN2FC			! Proportion of N2 produced at field capacity
10            PN2ZERO		! Proportion of N2 produced at wilting point
200           PART50        ! Soil N at which 50:50 N2:N2O is released
"""
    with open(join(sim_dir, 'Nitpars.dat'), 'w', newline='') as fobj:
        fobj.write(strng)
