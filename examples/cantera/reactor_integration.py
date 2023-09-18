#!/usr/bin/env python3

################################################################################
# Raffaele Cheula, LCCP, Politecnico di Milano, cheula.raffaele@gmail.com
################################################################################

import os
import csv
import timeit
import cantera as ct
import numpy as np
import copy as cp
import matplotlib.pyplot as plt
from pymatgen.core.lattice import Lattice
from shape.thermochemistry.wulff_with_support import WulffShape
from shape.units import *
from shape.ase_utils import convert_miller_index
from shape.thermochemistry.ab_initio_thermodynamics import (
    wulff_shape_from_diameter,
)
from shape.cantera_utils import (
    get_H0_0K_dict                  ,
    get_deltamu_dict                ,
    get_std_gibbs_dict              ,
    get_enthalpies_dict             ,
    get_std_entropies_dict          ,
    update_kinetics                 ,
    update_kinetics_facet           ,
    advance_to_steady_state         ,
    reaction_path_analysis          ,
    calculate_steady_state_coverages,
    get_reactions_data              ,
    get_energy_path                 ,
    plot_energy_paths               ,
)

################################################################################
# MEASURE TIME START
################################################################################

measure_time = True

if measure_time is True:
    start = timeit.default_timer()

################################################################################
# CALCULATIONS OPTIONS
################################################################################

reaction = 'revWGS'

results_dir = 'results'
csv_file    = 'results.csv'
RPA_file    = 'RPA.txt'
DRC_file    = 'DRC.txt'

"""
Discretization parameters.
"""

n_cstr   = 100+1
n_frames = 10

cstr_length_fix = 1e-08

"""
Composition parameters.
"""

fixed_comp    = False
conversion    = 0.20
eta_reaction  = None # 1.e+02
concentration = None # (1.00*atm*0.10)/(Rgas*273.15*Kelvin)
z_analysis    = 0.022 * centimeter

"""
Reaction path analysis (RPA).
"""

RPA = True

if reaction == 'WGS':
    react_select = [2, 6, 12, 13]
    spec_select_RPA = ['CO2(*)']

elif reaction == 'revWGS':
    react_select = [3, 7, 12, 14]
    spec_select_RPA = ['H2O(*)']

elif reaction == 'SR':
    react_select = []
    spec_select_RPA = ['C(*)']

elif reaction == 'DR':
    react_select = []
    spec_select_RPA = ['C(*)']

elif reaction == 'MET':
    react_select = []
    spec_select_RPA = ['CO(*)']

"""
Degree of rate control (DRC)
"""

DRC       = True
DRC_pts   = 1
DRC_thold = 0.01

""""
Apparent activation energy (Eapp).
"""
Eapp      = False
deltaT    = 5. # [K]

"""
Wulff construction.
"""

n_wulff    = 1
plot_wulff = False

"""
Print outlet composition and DRCs.
"""

var_type = None
var_num  = 1

print_outputs = False
print_csv     = False

if reaction in ('WGS', 'revWGS'):
    gas_selected = ['CO', 'H2O', 'CO2', 'H2']
else:
    gas_selected = ['CO', 'H2O', 'CO2', 'H2', 'CH4']

ads_selected = ['Rh', 'CO', 'H', 'H2O', 'CO2', 'OH', 'tCOOH']

step_selected = [6, 7, 8, 12, 14]

"""
Plots parameters.
"""

plot_profiles = True
all_plots     = False

plot_rates    = False
plot_paths    = False

y_max_plot = 0.02
width      = 0.8
tick_size  = 14
label_size = 16

"""
Parameters of the model.
"""

lateral_inter = True
update_kin    = True
fixed_shape   = False

"""
File with thermodynamic and kinetics parameters.
"""

cti_file = 'WGS_on_Rh_HRHT.cti'

################################################################################
# OPERATIVE CONDITIONS
################################################################################

temperature_celsius = 450.

temperature = Celsius_to_Kelvin(temperature_celsius)
pressure    =  1.00 * atm

gas_molar_fracs = {}
molar_fracs     = 0.01

if reaction == 'WGS':

    main_reactant = 'H2O'
    main_product  = 'H2'

    gas_molar_fracs['CO']  = 1*molar_fracs
    gas_molar_fracs['H2O'] = 1*molar_fracs
    gas_molar_fracs['CO2'] = 0.0
    gas_molar_fracs['H2']  = 0.0
    gas_molar_fracs['O2']  = 0.0
    gas_molar_fracs['CH4'] = 0.0

elif reaction == 'revWGS':

    main_reactant = 'CO2'
    main_product  = 'CO'

    gas_molar_fracs['CO']  = 0.0
    gas_molar_fracs['H2O'] = 0.0
    gas_molar_fracs['CO2'] = 1*molar_fracs
    gas_molar_fracs['H2']  = 1*molar_fracs
    gas_molar_fracs['O2']  = 0.0
    gas_molar_fracs['CH4'] = 0.0

elif reaction == 'SR':

    main_reactant = 'CH4'
    main_product  = 'H2'
    
    gas_molar_fracs['H2O'] = 1*molar_fracs
    gas_molar_fracs['CH4'] = 1*molar_fracs
    gas_molar_fracs['CO']  = 0.0
    gas_molar_fracs['CO2'] = 0.0
    gas_molar_fracs['H2']  = 0.0
    gas_molar_fracs['O2']  = 0.0

elif reaction == 'DR':

    main_reactant = 'CH4'
    main_product  = 'H2'

    gas_molar_fracs['CH4'] = 1*molar_fracs
    gas_molar_fracs['CO2'] = 1*molar_fracs
    gas_molar_fracs['CO']  = 0.0
    gas_molar_fracs['H2O'] = 0.0
    gas_molar_fracs['H2']  = 0.0
    gas_molar_fracs['O2']  = 0.0


elif reaction == 'MET':

    main_reactant = 'CO'
    main_product  = 'CH4'

    gas_molar_fracs['CO']  = 1*molar_fracs
    gas_molar_fracs['H2']  = 3*molar_fracs
    gas_molar_fracs['H2O'] = 0.0
    gas_molar_fracs['CO2'] = 0.0
    gas_molar_fracs['CH4'] = 0.0
    gas_molar_fracs['O2']  = 0.0

else:
    raise NameError('Wrong reaction name')

gas_molar_fracs['N2'] = 1.-sum([gas_molar_fracs[s]
                                for s in gas_molar_fracs if s != 'N2'])

################################################################################
# REACTOR PARAMETERS
################################################################################

"""
Set the reactor length and the gas flow velocity.
"""
reactor_length = 2.20 * centimeter
gas_velocity   = 4.40 * centimeter/second

"""
Set internal and external diameter of the annular reactor.
"""
diameter_ext = 5.00 * millimeter
diameter_int = 4.00 * millimeter

"""
Calculate the cross section and the total volume of the reactor.
"""
cross_section  = np.pi*(diameter_ext**2-diameter_int**2)/4.
reactor_volume = cross_section*reactor_length

################################################################################
# CATALYST PARAMETERS
################################################################################

"""
Set the catalyst mass and the fraction of active phase.
"""
catalyst_mass  = 10.00 * milligram
cat_percentage =  0.04

"""
Calculate gas velocity from volumetric flow rate.
"""

volumetric_flow_rate = gas_velocity*cross_section

molar_flow_rate = volumetric_flow_rate*pressure/Rgas/temperature # [kmol/s]
molar_flow_rate *= (kilomole/mole) * (minute/second) # [mol/min]
molar_flow_rate *= 1/(catalyst_mass/gram) # [mol/min/gram_cat]

"""
Set the catalyst nanoparticles diameter.
"""
wulff_diameter = 4. * nanometer

"""
Set the catalyst parameters.
"""
lattice_constant = 0.383 * nanometer
cat_MW           = 102.9 * gram/mole
cat_site_density = 2.49e-09 * mole/centimeter**2
cat_dispersion   = 0.05

"""
Calculate the catalyst active area per unit of reactor volume (alpha_cat),
employed if we do not know the catalyst shape.
"""
cat_moles = catalyst_mass*cat_percentage/cat_MW
cat_area  = cat_moles*cat_dispersion/cat_site_density
alpha_cat = cat_area/reactor_volume

"""
Calculate the catalyst volume per unit of reactor volume (beta_cat),
employed if we have a model for the catalyst shape.
"""
lattice_volume = lattice_constant**3/4.
cat_volume     = lattice_volume*cat_moles*Navo
beta_cat       = cat_volume/reactor_volume

################################################################################
# PHASES
################################################################################

"""
Create a solution representing the gas phase.
"""
gas = ct.Solution(cti_file, 'gas')
gas.TPX = temperature, pressure, gas_molar_fracs

"""
Create interfaces representing the catalyst facets.
"""
active_phases = ['Rh100', 'Rh110', 'Rh111', 'Rh311', 'Rh331']
site_names    = [  '100',   '110',   '111',   '311',   '331']
miller_list   = [  '100',   '110',   '111',   '311',   '331']

free_site = 'Rh'

cat_list = []
for phase in active_phases:
    cat = ct.Interface(cti_file, phase, [gas])
    cat.TP = temperature, pressure
    cat_list += [cat]

n_ads_species = max([len(cat.species_names) for cat in cat_list])

ads_names = [s.split('(')[0] for s in cat_list[0].species_names]

colors = {'100': 'r'     ,
          '110': 'orange',
          '111': 'b'     ,
          '311': 'm'     ,
          '331': 'g'     }

colors_dict = {}
custom_colors = {}
for facet in miller_list:
    custom_colors[convert_miller_index(facet)] = colors[facet]
    colors_dict[free_site+facet] = colors[facet]

################################################################################
# INITIAL CATALYST MORPHOLOGY
################################################################################

"""
Import the functions which represent how the specific surface free energy
of the catalyst facets change with the change of the gas phase chemical
potentials.
"""

from e_surf_funs_CO_on_Rh import e_surf_fun_dict, miller_supp, e_surf_supp

H0_0K_dict = get_H0_0K_dict(phase = gas, units = eV/molecule)

if len(miller_list) == 1:
    miller_supp = miller_list[0]
    e_surf_supp = None

for i in range(len(miller_list)):
    miller_list[i] = convert_miller_index(miller_list[i])

"""
Given the nanoparticle diameter, calculate its volume. The gas chemical
potentials have to be the ones of the gas phase during the measurement of
the nanoparticle diameter.
"""

deltamu_CO = -2.4 # [eV]

e_surf_list = np.zeros(len(miller_list))

for i in range(len(miller_list)):

    facet = miller_list[i]
    e_surf_list[i] = e_surf_fun_dict[facet](deltamu_CO)

lattice = Lattice(np.identity(3))

wulff = wulff_shape_from_diameter(lattice        = lattice       ,
                                  miller_list    = miller_list   ,
                                  e_surf_list    = e_surf_list   ,
                                  miller_supp    = miller_supp   ,
                                  e_surf_supp    = e_surf_supp   ,
                                  wulff_diameter = wulff_diameter)

wulff_volume = wulff.volume

################################################################################
# NUMERICAL PARAMETERS
################################################################################

""" 
The pfr is solved using a series of cstr, which represent the control
volumes that the gas encounter proceding along the reactor. The inlet
composition of each cstr is the composition of the previous cstr (or the
reactor inlet composition for the first cstr). After each integration,
the control volume representing each cstr is reinitialized, then it is
riutilized fort the next integration.
At each timestep, the coverages of adsorbates in equilibrium with the gas
phase are calculated and the kinetic parameters are updated because of
coverage effects. Then, the cstr is solved (changing also the gas phase
composition) with the kinetics obtained at the beginning of the integration
step. 
""" 

cstr_length = reactor_length/(n_cstr-1)
cstr_volume = cross_section*cstr_length

cstr_cat_area   = alpha_cat*cstr_volume
cstr_cat_volume = beta_cat*cstr_volume

mass_flow_rate = gas_velocity*gas.density*cross_section

"""
Create a reactor representing the control volumes to be integrated.
We assign cstr_cat_area as area of the facets but then we substitute it.
"""

cstr = ct.IdealGasReactor(gas, energy = 'off')
cstr.volume = cstr_volume

surf_list = []
for i in range(len(cat_list)):
    surf = ct.ReactorSurface(cat_list[i], cstr, A = cstr_cat_area)
    surf_list += [surf]

upstream = ct.Reservoir(gas, name = 'upstream')
master = ct.MassFlowController(upstream, cstr, mdot = mass_flow_rate)

downstream = ct.Reservoir(gas, name = 'downstream')
pcontrol = ct.PressureController(cstr, downstream, master = master, K = 1e-6)

water = ct.Water()
water.TP = temperature, pressure

#heater = ct.Reservoir(water, name = 'heater')
#tcontrol = ct.Wall(cstr, heater, U = 1e+6)

"""
Define the parameters of the simulation.
"""
sim = ct.ReactorNet([cstr])
sim.rtol = 1.0e-14
sim.atol = 1.0e-22
sim.max_err_test_fails = 1e9

if fixed_comp is False:
    sim.rtol = 1.0e-09
    sim.atol = 1.0e-18

################################################################################
# ETA REACTION
################################################################################

if eta_reaction is not None:

    x_tot = 1.-gas_molar_fracs['N2']

    gas.equilibrate('TP')

    if reaction == 'WGS':
    
        Keq = (gas['CO'].X[0]*gas['H2O'].X[0])/(gas['CO2'].X[0]*gas['H2'].X[0])
        
        denominator = 2.+2./np.sqrt(Keq*eta_reaction)
    
        gas_molar_fracs['CO']  = x_tot/denominator
        gas_molar_fracs['H2O'] = x_tot/denominator
        gas_molar_fracs['CO2'] = x_tot/np.sqrt(Keq*eta_reaction)/denominator
        gas_molar_fracs['H2']  = x_tot/np.sqrt(Keq*eta_reaction)/denominator
    
    elif reaction == 'revWGS':
    
        Keq = (gas['CO2'].X[0]*gas['H2'].X[0])/(gas['CO'].X[0]*gas['H2O'].X[0])
        
        denominator = 2.+2./np.sqrt(Keq*eta_reaction)
        
        gas_molar_fracs['CO']  = x_tot/np.sqrt(Keq*eta_reaction)/denominator
        gas_molar_fracs['H2O'] = x_tot/np.sqrt(Keq*eta_reaction)/denominator
        gas_molar_fracs['CO2'] = x_tot/denominator
        gas_molar_fracs['H2']  = x_tot/denominator

    gas.TPX = temperature, pressure, gas_molar_fracs

################################################################################
# MODIFY GAS COMPOSITION
################################################################################

if var_type in gas_selected:

    gas_molar_fracs[var_type] *= (1.+var_num)

    gas_molar_fracs['N2'] = 1.-sum([gas_molar_fracs[s]
                                    for s in gas_molar_fracs if s != 'N2'])

    gas.TPX = temperature, pressure, gas_molar_fracs

    ads_selected = ['Rh', 'CO', 'H']

################################################################################
# FIXED COMPOSITION
################################################################################

if fixed_comp is True:

    plot_profiles = False

    z_analysis = 0.000 * centimeter

    n_cstr = 1

    cstr_length = cstr_length_fix

    cstr_volume = cross_section*cstr_length
    cstr.volume = cstr_volume

    cstr_cat_volume = beta_cat*cstr_volume

################################################################################
# LATERAL INTERACTIONS
################################################################################

"""
Update kinetic parameters of the surface reactions, calculated from
differences in Gibbs free energies between transition states and
reactants.
"""

from lateral_interactions import lat_dict

if lateral_inter is False:

    for facet in lat_dict:
    
        for lat in lat_dict[facet]:
    
            lat.coeffs = [0., 0.]
            lat.x0 = 0.
            lat.m  = 0.

"""
Create transition states objects and store the NASA coefficients
of the adsorbates before accounting for the lateral interactions.
"""

TS_list = []
coeffs_dict = {}

for spec in gas.species_names:

    index = gas.species_names.index(spec)
    species = gas.species(index)

    coeffs = species.thermo.coeffs
    coeffs_dict[spec] = coeffs

for i in range(len(active_phases)):

    cat  = cat_list[i]
    site = site_names[i]

    TS = ct.Solution(cti_file, active_phases[i]+'-TS', [gas])
    TS.TP = temperature, pressure
    TS_list += [TS]

    for spec in cat.species_names:

        index = cat.species_names.index(spec)
        species = cat.species(index)

        coeffs = species.thermo.coeffs
        coeffs_dict[spec] = coeffs

    for spec in TS.species_names:

        index = TS.species_names.index(spec)
        species = TS.species(index)

        coeffs = species.thermo.coeffs
        coeffs_dict[spec] = coeffs

################################################################################
# CALCULATE WULFF SHAPE
################################################################################

"""
Calculate facets areas from the nanoparticle Wulff shape in reacting
conditions. The scale_factor parameter is used to mantain the same
volume (number of atoms) when changin the catalyst shape.
"""

deltamu = get_deltamu_dict(phase      = gas        ,
                           H0_0K_dict = H0_0K_dict ,
                           units      = eV/molecule)

#e_surf_list = surface_energies(e_surf_funs = e_surf_funs  ,
#                               deltamu_A   = deltamu['CO'],
#                               deltamu_B   = deltamu['H2'])

e_surf_list = np.zeros(len(miller_list))

for i in range(len(miller_list)):
    e_surf_list[i] = e_surf_fun_dict[miller_list[i]](deltamu['CO'])

wulff = WulffShape(lattice     = lattice    ,
                   miller_list = miller_list,
                   e_surf_list = e_surf_list,
                   miller_supp = miller_supp,
                   e_surf_supp = e_surf_supp)

scale_factor = (wulff_volume/wulff.volume)**(1./3.)

for i in range(len(miller_list)):
    facet = miller_list[i]
    spec_area = wulff.miller_area_dict[facet]*scale_factor**2/wulff_volume
    area = cstr_cat_volume*spec_area
    surf_list[i].area = area

################################################################################
# CALCULATE INITIAL KINETICS
################################################################################

"""
Define the parameters for the calculation of the steady state coverages.
"""

x_min_CO = 1e-4

mu_gas_dict = {}

mu_gas_dict['CO'] = lambda gas: gas['CO'].chemical_potentials[0]
mu_gas_dict['H']  = lambda gas: gas['H2'].chemical_potentials[0]/2.

method = 'df-sane'
options = {'maxfev': 10000}

"""
Update the kinetics at the inlet of the reactor according to gas composition.
"""

if update_kin is True:

    TDY = gas.TDY
    
    mari_list = []
    
    if gas['CO'].X[0] > x_min_CO:
        mari_list += ['CO']
    
    calculate_steady_state_coverages(gas         = gas        ,
                                     cat_list    = cat_list   ,
                                     TS_list     = TS_list    ,
                                     site_names  = site_names ,
                                     lat_dict    = lat_dict   ,
                                     coeffs_dict = coeffs_dict,
                                     free_site   = free_site  ,
                                     mari_list   = mari_list  ,
                                     mu_gas_dict = mu_gas_dict,
                                     options     = options    ,
                                     method      = method     )
    
    update_kinetics(gas           = gas          ,
                    cat_list      = cat_list     ,
                    TS_list       = TS_list      ,
                    site_names    = site_names   ,
                    active_phases = active_phases,
                    lat_dict      = lat_dict     ,
                    coeffs_dict   = coeffs_dict  )
    
    gas.TDY = TDY

################################################################################
# PRINT INLET COMPOSITION
################################################################################

if print_outputs is True:

    os.makedirs(results_dir, exist_ok = True)

    filename = results_dir+'/INLET_{:+06.2f}'.format(var_num)

    f = open(filename, 'w+')

    for spec in gas_selected:
        print('{0:10s} {1:+11.9f}'.format(spec, *gas[spec].X), file = f)

    f.close()

################################################################################
# PRE PROCESSING
################################################################################

print('\n OPERATIVE CONDITIONS \n')
print('temperature     = {:9.2f} K'.format(temperature))
print('pressure        = {:9.2e} atm'.format(pressure/atm))
print('molar flow rate = {:9.4f} mol/min/gcat'.format(molar_flow_rate))

if fixed_comp is True:

    print('\ncstr_length  = {:9.2e}'.format(cstr_length))

    if eta_reaction is not None:
        print('eta_reaction = {:9.2e}'.format(eta_reaction))

    elif conversion is not None:
        print('conversion   = {:9.2e}'.format(conversion))

"""
Prepare solution vectors and plots.
"""

if print_csv is True:
    fileobj = open(csv_file, 'w')
    writer = csv.writer(fileobj)
    writer.writerow(['Distance [m]', 'T [K]', 'P [Pa]'] + gas.species_names +
                    [spec for cat in cat_list for spec in cat.species_names])

print('\n AXIAL PROFILES \n')

string = 'distance [m]'.rjust(14)
for spec in gas.species_names:
    string += ('x_'+spec+' [-]').rjust(12)
print(string)

z_vector = np.array(range(n_cstr))*reactor_length/n_cstr
gas_molar_fracs_vectors = np.zeros([n_cstr, len(gas.species_names)])
ads_coverages_matrices = np.zeros([n_cstr, len(cat_list), n_ads_species])
cat_areas_vectors = np.zeros([n_cstr, len(cat_list)])

fig_num = 0

if plot_wulff is True:
    fig_wulff = plt.figure(fig_num)
    fig_wulff.set_size_inches(8, 8)
    fig_wulff.show()

################################################################################
# INTEGRATION
################################################################################

wulff_frames = []

extra = ['z'] #, 'cat_areas', 'ads_coverages']
solution = ct.SolutionArray(gas, extra = extra)

for n in range(n_cstr):

    z = n*cstr_length

    """
    The gas composition of the inlet of the cstr is set as equal to the
    composition of the previous cstr.
    """

    if n > 0:
        gas.TDY = cstr.thermo.TDY
        upstream.syncState()

    """
    Calculate facets areas from the nanoparticle Wulff shape in reacting
    conditions. The scale_factor parameter is used to mantain the same
    volume (number of atoms) when changin the catalyst shape.
    """

    deltamu = get_deltamu_dict(phase      = gas        ,
                               H0_0K_dict = H0_0K_dict ,
                               units      = eV/molecule)
    
    #e_surf_list = surface_energies(e_surf_funs = e_surf_funs  ,
    #                               deltamu_A   = deltamu['CO'],
    #                               deltamu_B   = deltamu['H2'])
    
    if fixed_shape is True:
        deltamu['CO'] = -2.4
    
    if deltamu['CO'] < -2.4:
        deltamu['CO'] = -2.4
    elif deltamu['CO'] > 0.0:
        deltamu['CO'] = 0.0

    for i in range(len(miller_list)):
        e_surf_list[i] = e_surf_fun_dict[miller_list[i]](deltamu['CO'])
    
    wulff = WulffShape(lattice     = lattice    ,
                       miller_list = miller_list,
                       e_surf_list = e_surf_list,
                       miller_supp = miller_supp,
                       e_surf_supp = e_surf_supp)
    
    scale_factor = (wulff_volume/wulff.volume)**(1./3.)
    
    for i in range(len(miller_list)):
        facet = miller_list[i]
        spec_area = wulff.miller_area_dict[facet]*scale_factor**2/wulff_volume
        area = cstr_cat_volume*spec_area
        surf_list[i].area = area
        cat_areas_vectors[n,i] = area

    """
    Store information of the catalyst morphology and composition of gas
    and catalyst facets.
    """

    gas_molar_fracs_vectors[n] = gas.X

    for j in range(len(cat_list)):
        coverage = cat_list[j].coverages
        coverage.resize(n_ads_species)
        ads_coverages_matrices[n,j,:] = coverage

    if fixed_comp is True or not n % int((n_cstr-1)/n_frames):
        
        """
        Print gaseous species axial profiles.
        """

        string_old = cp.deepcopy(string)
        string = '  {0:12f}'.format(z)
        for spec in gas.species_names:
            string += '  {0:10f}'.format(*gas[spec].X)
        print(string)
    
    if plot_wulff is True:

        if n_wulff != 1 and not n % int((n_cstr-1)/n_wulff):

            """
            Plot Wulff construction shape as animation.
            """

            wulff_frames.append(wulff.get_plot(azim          = -135.        ,
                                               elev          = -100.        ,
                                               custom_colors = custom_colors,
                                               fig           = fig_wulff    ))

            plt.pause(0.001)

        elif n_wulff == 1 and abs(z_analysis-z) < cstr_length-1e-7:
    
            """
            Plot Wulff construction shape.
            """

            wulff.get_plot(azim          = -135.        ,
                           elev          = -100.        ,
                           custom_colors = custom_colors,
                           fig           = fig_wulff    )
    
            plt.pause(0.001)

    """
    Store solution in array and csv file.
    """
    
    solution.append(gas.state,
                    z             = z                        ,
                    cat_areas     = cat_areas_vectors[n]     ,
                    ads_coverages = ads_coverages_matrices[n])
    
    if print_csv is True:
        writer.writerow([z, gas.T, gas.P] + list(gas.X) +
                        [x for s in surf_list for x in s.coverages])

    """
    Integrate the cstr.
    """

    if update_kin is True:

        mari_list = []
        
        if gas['CO'].X[0] > x_min_CO:
            mari_list += ['CO']
        
        calculate_steady_state_coverages(gas         = gas        ,
                                         cat_list    = cat_list   ,
                                         TS_list     = TS_list    ,
                                         site_names  = site_names ,
                                         lat_dict    = lat_dict   ,
                                         coeffs_dict = coeffs_dict,
                                         free_site   = free_site  ,
                                         mari_list   = mari_list  ,
                                         mu_gas_dict = mu_gas_dict,
                                         options     = options    ,
                                         method      = method     )
        
        update_kinetics(gas           = gas          ,
                        cat_list      = cat_list     ,
                        TS_list       = TS_list      ,
                        site_names    = site_names   ,
                        active_phases = active_phases,
                        lat_dict      = lat_dict     ,
                        coeffs_dict   = coeffs_dict  )
    
    advance_to_steady_state(sim = sim)

if print_csv is True:
    fileobj.close()
    print('\nPrinted profiles into file: {}'.format(csv_file))

################################################################################
# PRINT OUTLET COMPOSITION
################################################################################

gas.equilibrate('TP')

if print_outputs is True:

    filename = results_dir+'/OUTLET_{:+06.2f}'.format(var_num)

    f = open(filename, 'w+')

    for spec in gas_selected:
        print('{0:10s} {1:+11.9f}'.format(spec, *gas[spec].X), file = f)

    for i in range(len(cat_list)):

        cat = cat_list[i]

        for spec in ads_selected:

            site_name = site_names[i]
            ads_name = spec+'('+site_name+')'

            print('{0:10s} {1:+7.5e}'.format(ads_name, 
                  *cat[ads_name].coverages), file = f)

    f.close()

################################################################################
# POST PROCESSING
################################################################################

if plot_profiles is True and all_plots is False:

    fig_num += 1
    
    fig = plt.figure(fig_num)
    plt.style.context('bmh')

    """
    Plot gas species molar fractions profiles.
    """

    plt.axis([0., max(z_vector), 0., y_max_plot])
    plt.title('gas molar fractions')
    for i in range(len(gas.species_names)):
        plt.plot(z_vector, gas_molar_fracs_vectors[:,i],
                 label = gas.species_names[i])
    plt.legend(bbox_to_anchor = (-0.15, 1), loc = 1, borderaxespad = 0.)

elif plot_profiles is True and all_plots is True:

    fig_num += 1

    fig = plt.figure(fig_num)
    fig.set_size_inches(16, 16)
    plt.style.context('bmh')

    """
    Plot gas species molar fractions profiles.
    """

    plt.subplot(331)
    #plt.subplot(321)
    plt.axis([0., max(z_vector), 0., y_max_plot])
    plt.title('gas molar fractions')
    for i in range(len(gas.species_names)):
        plt.plot(z_vector, gas_molar_fracs_vectors[:,i],
                 label = gas.species_names[i])
    plt.legend(bbox_to_anchor = (-0.15, 1), loc = 1, borderaxespad = 0.)

    """
    Plot catalyst shape. [TODO]
    """

    #plt.subplot(332)

    """
    Plot catalyst facets areas.
    """

    plt.subplot(333)
    #plt.subplot(322)
    plt.axis([0., max(z_vector), 0., 1.1*np.max(cat_areas_vectors)])
    plt.title('catalyst facets areas')
    for i in range(len(cat_list)):
        plt.plot(z_vector[1:], cat_areas_vectors[1:,i],
                 label = active_phases[i])
    plt.legend(bbox_to_anchor = (1.05, 1), loc = 2, borderaxespad = 0.)

    """
    Plot adsorbates coverages profiles.
    """

    for j in range(len(cat_list)):
        plt.subplot(334+j)
        #plt.subplot(323+j)
        cat = cat_list[j]
        plt.axis([0., max(z_vector), 0., 1.])
        plt.title('{} adsorbates coverages'.format(cat.name))
        for i in range(len(cat.species_names)):
            plt.plot(z_vector[1:], ads_coverages_matrices[1:,j,i],
                     label = cat.species_names[i])
        if j == 0:
            plt.legend(bbox_to_anchor = (-0.15, 1), loc = 1, borderaxespad = 0.)

################################################################################
# THERMODYNAMIC COMPOSITION
################################################################################

"""
Calculate and print the thermodynamic composition of the gas phase.
"""

print('\n TD EQ COMPOSITION \n')

gas.equilibrate('TP')

string = 'distance [m]'.rjust(14)
for spec in gas.species_names:
    string += ('x_'+spec+' [-]').rjust(12)
print(string)

string = 'infinity'.rjust(14)
for spec in gas.species_names:
    string += '  {0:10f}'.format(*gas[spec].X)
print(string)

################################################################################
# POST PROCESSING
################################################################################

print('\n\n POST PROCESSING \n')

cov_min_print = 1e-4

sim.rtol = 1.0e-14
sim.atol = 1.0e-22

cstr.volume = cross_section*cstr_length_fix

index = np.argmin(np.abs(solution.z-z_analysis))

residence_time = solution.z[index]/gas_velocity

print('Distance       = {:.3e} m'.format(solution.z[index]))
print('Residence time = {:.3e} s'.format(residence_time))

gas.TDY = solution[index].TDY
upstream.syncState()

area_spec_dict = {}

for i in range(len(surf_list)):

    #area = solution.cat_areas[index][i]*cstr.volume/cstr_volume
    area = cat_areas_vectors[index][i]*cstr.volume/cstr_volume

    phase = active_phases[i]

    surf_list[i].area = area

    #cat.coverages = solution.ads_coverages[index][i,:cat.n_species]
    cat.coverages = ads_coverages_matrices[index][i,:cat.n_species]
    
    area_spec_dict[phase] = area/cstr.volume

if update_kin is True:

    mari_list = []
    
    if gas['CO'].X[0] > x_min_CO:
        mari_list += ['CO']
    
    calculate_steady_state_coverages(gas         = gas        ,
                                     cat_list    = cat_list   ,
                                     TS_list     = TS_list    ,
                                     site_names  = site_names ,
                                     lat_dict    = lat_dict   ,
                                     coeffs_dict = coeffs_dict,
                                     free_site   = free_site  ,
                                     mari_list   = mari_list  ,
                                     mu_gas_dict = mu_gas_dict,
                                     options     = options    ,
                                     method      = method     )
    
    update_kinetics(gas           = gas          ,
                    cat_list      = cat_list     ,
                    TS_list       = TS_list      ,
                    site_names    = site_names   ,
                    active_phases = active_phases,
                    lat_dict      = lat_dict     ,
                    coeffs_dict   = coeffs_dict  )

advance_to_steady_state(sim = sim)

print('\n GAS MOLAR FRACTIONS \n')

for spec in gas_selected:
    print('{0:4s} = {1:13.10f}'.format(spec, gas[spec].X[0]))

print('\n CATALYST FACETS AREAS \n')

print(' '*13+'[m^2/m^3]'+' '*5+'[%]')

area_spec_sum = sum([area_spec_dict[phase] for phase in area_spec_dict])

for i in range(len(cat_list)):

    phase = active_phases[i]

    print('{:4s} ->'.format(phase), end = '  ')
    print('{:12.4e}'.format(area_spec_dict[phase]), end = ' ')
    print('{:7.3f}'.format(area_spec_dict[phase]/area_spec_sum*100))

print(' tot  ->', end = '  ')
print('{:12.4e}'.format(area_spec_sum), end = ' ')
print('{:7.3f}'.format(100.))

print('\n MARIS COVERAGES \n')

cove_dict = {}
cove_max_dict = {}
cove_ave_dict = {}

for spec in ads_names:

    cove_dict[spec] = {}
    cove_max_dict[spec] = 0.
    cove_ave_dict[spec] = 0.

for i in range(len(cat_list)):

    cat = cat_list[i]
    phase = active_phases[i]

    coverages = cat.coverages
    coverages.resize(n_ads_species)

    for j in range(len(ads_names)):

        spec = ads_names[j]
        cove_dict[spec][phase] = coverages[j]

        if coverages[j] > cove_max_dict[spec]:
            cove_max_dict[spec] = coverages[j]

print(' '*9, end = ' ')

for j in range(len(ads_names)):

    spec = ads_names[j]

    if cove_max_dict[spec] > cov_min_print:

        print('{:7s}'.format(spec), end = ' ')

for i in range(len(cat_list)):

    phase = active_phases[i]

    print('\n{:4s} ->'.format(phase), end = ' ')

    for j in range(len(ads_names)):
    
        spec = ads_names[j]
    
        cove_ave_dict[spec] += (cove_dict[spec][phase] * 
                                area_spec_dict[phase]/area_spec_sum)
    
        if cove_max_dict[spec] > cov_min_print:
    
            print('{:7.4f}'.format(cove_dict[spec][phase]), end = ' ')

print('\n ave  ->', end = ' ')

for j in range(len(ads_names)):
    
    spec = ads_names[j]

    if cove_max_dict[spec] > cov_min_print:

        print('{:7.4f}'.format(cove_ave_dict[spec]), end = ' ')

"""

print('\n\n PUNCTUAL REACTION RATES')

reac_dict = {}
TOFs_dict = {}
prod_dict = {}

for s in range(len(gas.species_names)):

    spec = gas.species_names[s]

    reac_dict[spec] = {}
    TOFs_dict[spec] = {}
    prod_dict[spec] = {}

    for i in range(len(cat_list)):

        cat = cat_list[i]
        phase = active_phases[i]
        
        area_spec = area_spec_dict[phase]

        reacs = cat.get_net_production_rates(gas)
        TOFs  = cat.get_net_production_rates(gas)/cat.site_density
        prods = cat.get_net_production_rates(gas)*area_spec

        reac_dict[spec][phase] = reacs[s]
        TOFs_dict[spec][phase] = TOFs[s]
        prod_dict[spec][phase] = prods[s]

for spec in gas_selected:

    print('\nProduction of {} \n'.format(spec))
    print(' '*11+'[kmol/m^2/s]'+' '*9+'[1/s]'+' '*1,
          '[kmol/m^3/s]'+' '*5+'[%]')

    reac_sum = sum([reac_dict[spec][phase] for phase in reac_dict[spec]])
    TOFs_sum = sum([TOFs_dict[spec][phase] for phase in TOFs_dict[spec]])
    prod_sum = sum([prod_dict[spec][phase] for phase in prod_dict[spec]])

    for i in range(len(cat_list)):

        phase = active_phases[i]

        print('{:4s} ->'.format(phase), end = '  ')
        print('{:+13.4e}'.format(reac_dict[spec][phase]), end = ' ')
        print('{:+13.4e}'.format(TOFs_dict[spec][phase]), end = ' ')
        print('{:+13.4e}'.format(prod_dict[spec][phase]), end = ' ')
        print('{:7.3f}'.format(prod_dict[spec][phase]/prod_sum*100))

    print(' tot  ->', end = '  ')
    print('{:+13.4e}'.format(reac_sum), end = ' ')
    print('{:+13.4e}'.format(TOFs_sum), end = ' ')
    print('{:+13.4e}'.format(prod_sum), end = ' ')
    print('{:7.3f}'.format(100.))

    if spec == main_product:
        reac_sum_main = reac_sum
        TOFs_sum_main = TOFs_sum
        prod_sum_main = prod_sum

print('\nTotal production of {} \n'.format(main_product))

print('Rtot = {:+13.4e} kmol/m^3/s'.format(prod_sum_main))

"""

print('\n ADSORPTION EQUILIBRIUM CONSTANTS')

for i in range(len(cat_list)):

    cat = cat_list[i]
    site = site_names[i]

    teta_free = cat[free_site+'('+site+')'].coverages[0]
    teta_CO = cat['CO('+site+')'].coverages[0]
    teta_H = cat['H('+site+')'].coverages[0]

    P_CO = gas['CO'].X[0]*pressure/atm
    P_H2 = gas['H2'].X[0]*pressure/atm

    Keq_CO = teta_CO/teta_free/P_CO
    Keq_H2 = teta_H**2/teta_free**2/P_H2

    print('\nKeq CO Rh({0}) = {1:.4e} 1/atm'.format(site, Keq_CO))
    print('Keq H2 Rh({0}) = {1:.4e} 1/atm'.format(site, Keq_H2))

"""

if plot_rates is True:

    fig_num += 1

    fig = plt.figure(fig_num)
    fig.set_size_inches(20, 8)

    alpha = 0.9

    plt.subplot(131)
    
    ylabel = 'Area [%]'
    
    for i in range(len(cat_list)):
        phase = active_phases[i]
        x     = phase
        y     = area_spec_dict[phase]/area_spec_sum*100
        color = colors_dict[phase]
        plt.bar(x, y, width, alpha = alpha, color = color,
                edgecolor = 'k', linewidth = 1.)
    
    plt.yticks(fontsize = tick_size)
    plt.xticks(fontsize = tick_size)

    plt.axis([-0.5, len(cat_list)-0.5, 0., 100.])

    plt.ylabel(ylabel, fontsize = label_size)
    
    plt.subplot(132)
    
    ylabel = 'TOF [%]'
    
    for i in range(len(cat_list)):
        phase = active_phases[i]
        x     = phase
        y     = TOFs_dict[main_product][phase]/TOFs_sum_main*100
        color = colors_dict[phase]
        plt.bar(x, y, width, alpha = alpha, color = color,
                edgecolor = 'k', linewidth = 1.)
    
    plt.yticks(fontsize = tick_size)
    plt.xticks(fontsize = tick_size)
    
    plt.axis([-0.5, len(cat_list)-0.5, 0., 100.])
    
    plt.ylabel(ylabel, fontsize = label_size)

    plt.subplot(133)
    
    ylabel = 'Prod [%]'
    
    for i in range(len(cat_list)):
        phase = active_phases[i]
        x     = phase
        y     = prod_dict[main_product][phase]/prod_sum_main*100
        color = colors_dict[phase]
        plt.bar(x, y, width, alpha = alpha, color = color,
                edgecolor = 'k', linewidth = 1.)
    
    plt.yticks(fontsize = tick_size)
    plt.xticks(fontsize = tick_size)
    
    plt.axis([-0.5, len(cat_list)-0.5, 0., 100.])
    
    plt.ylabel(ylabel, fontsize = label_size)

if print_outputs is True:

    filename = results_dir+'/PRODUCTIVITIES_{:+06.2f}'.format(var_num)

    f = open(filename, 'w+')

    for i in range(len(cat_list)):

        phase = active_phases[i]
        area_rel = area_spec_dict[phase]/area_spec_sum

        print('Area_rel_{0:6s} {1:9.7f}'.format(phase, area_rel), file = f)

    for i in range(len(cat_list)):

        phase = active_phases[i]
        TOFs_rel = TOFs_dict[main_product][phase]/TOFs_sum_main

        print('TOFs_rel_{0:6s} {1:9.7f}'.format(phase, TOFs_rel), file = f)

    for i in range(len(cat_list)):

        phase = active_phases[i]
        prod_rel = prod_dict[main_product][phase]/prod_sum_main

        print('Prod_rel_{0:6s} {1:9.7f}'.format(phase, prod_rel), file = f)

    for i in range(len(cat_list)):

        phase = active_phases[i]
        area_abs = area_spec_dict[phase]

        print('Area_abs_{0:6s} {1:9.3e}'.format(phase, area_abs), file = f)

    print('Area_abs_total  {0:9.3e}'.format(area_spec_sum), file = f)

    for i in range(len(cat_list)):

        phase = active_phases[i]
        TOFs_abs = TOFs_dict[main_product][phase]

        print('TOFs_abs_{0:6s} {1:8.3e}'.format(phase, TOFs_abs), file = f)

    print('TOFs_abs_total  {0:9.3e}'.format(TOFs_sum_main), file = f)

    for i in range(len(cat_list)):

        phase = active_phases[i]
        prod_abs = prod_dict[main_product][phase]

        print('Prod_abs_{0:6s} {1:9.3e}'.format(phase, prod_abs), file = f)

    print('Prod_abs_total  {0:9.3e}'.format(prod_sum_main), file = f)

    f.close()

"""

################################################################################
# REACTION PATH ANALYSIS
################################################################################

print('\n REACTION PATH ANALYSIS \n')

if RPA is True:

    if print_outputs is True:
        filename_select = results_dir+'/PATHS_{:+06.2f}'.format(var_num)
    else:
        filename_select = None

    for i in range(len(spec_select_RPA)):
        spec = spec_select_RPA[i]
        if '*' in spec:
            del spec_select_RPA[i]
            for site_name in site_names:
                spec_select_RPA += [spec.replace('*', site_name)]

    reaction_path_analysis(gas             = gas            ,
                           cat_list        = cat_list       ,
                           surf_list       = surf_list      ,
                           filename        = RPA_file       ,
                           spec_select     = spec_select_RPA,
                           react_select    = react_select   ,
                           filename_select = filename_select)

################################################################################
# CALCUlATE GO
################################################################################

units = eV/molecule

G0_dict = get_std_gibbs_dict(gas, units = units)

G0_ref_dict = {}

if reaction == 'WGS':

    G0_ref_dict['O'] = G0_dict['H2O']-G0_dict['H2']
    G0_ref_dict['C'] = G0_dict['CO']-G0_ref_dict['O']
    G0_ref_dict['H'] = (G0_dict['H2O']-G0_ref_dict['O'])/2.

elif reaction == 'revWGS':

    G0_ref_dict['O'] = G0_dict['CO2']-G0_dict['CO']
    G0_ref_dict['C'] = G0_dict['CO2']-2*G0_ref_dict['O']
    G0_ref_dict['H'] = G0_dict['H2']/2.

G0_rel_TS_dict = {}
G0_rel_TS_selected = {}

for i in range(len(active_phases)):
    
    site_name = site_names[i]
    
    cat = cat_list[i]
    cat.TP = gas.TP
    
    TS = TS_list[i]
    TS.TP = gas.TP

    reactions_data = get_reactions_data(gas         = gas        ,
                                        cat         = cat        ,
                                        TS          = TS         ,
                                        units       = eV/molecule,
                                        G0_ref_dict = G0_ref_dict)

    for n in range(cat.n_reactions):

        G0_rel_TS = reactions_data['G0 rel'][n]

        step = cat.reaction_equation(n)

        G0_rel_TS_dict[step] = G0_rel_TS

        if n+1 in step_selected:
            G0_rel_TS_selected[step] = G0_rel_TS

if print_outputs is True:

    filename = results_dir+'/G0_TS_{:+06.2f}'.format(var_num)

    f = open(filename, 'w+')
        
    for step in G0_rel_TS_selected:

        name = step.replace(' ', '')

        print('{0:40s} {1:+9.6f}'.format(name, G0_rel_TS_dict[step]), file = f)

    f.close()

################################################################################
# PRINT REACTIONS DETAILS
################################################################################

if print_outputs is True:

    filename = results_dir+'/H0_G0_{:+06.2f}'.format(var_num)

    f = open(filename, 'w+')

    head_spec = '{0:40s} {1:12s} {2:12s} {3:12s}'
    head_reax = head_spec+' {4:12s} {5:12s} {6:12s} {7:12s}'

    print(head_spec.format('gas species',
                           'H0 [eV]'    ,
                           'S0 [eV/K]'  ,
                           'G0 [eV]'    ), file = f)

    H0_dict = get_enthalpies_dict(gas, units = units)
    S0_dict = get_std_entropies_dict(gas, units = units)
    G0_dict = get_std_gibbs_dict(gas, units = units)

    string_spec = '{0:40s} {1:+12.8f} {2:+12.4e} {3:+12.8f}'
    string_reax = string_spec+' {4:+12.8f} {5:+12.8f} {6:+12.8f} {7:+12.8f}'

    for s in range(gas.n_species):
        
        spec = gas.species_names[s]
        
        print(string_spec.format(spec         ,
                                 H0_dict[spec],
                                 S0_dict[spec],
                                 G0_dict[spec]), file = f)

    for i in range(len(active_phases)):
    
        site_name = site_names[i]
    
        cat = cat_list[i]
        cat.TP = gas.TP
    
        TS = TS_list[i]
        TS.TP = gas.TP

        H0_dict.update(get_enthalpies_dict(cat, units = units))
        H0_dict.update(get_enthalpies_dict(TS, units = units))
    
        S0_dict.update(get_std_entropies_dict(cat, units = units))
        S0_dict.update(get_std_entropies_dict(TS, units = units))

        G0_dict.update(get_std_gibbs_dict(cat, units = units))
        G0_dict.update(get_std_gibbs_dict(TS, units = units))

        reactions_data = get_reactions_data(gas         = gas        ,
                                            cat         = cat        ,
                                            TS          = TS         ,
                                            units       = eV/molecule,
                                            G0_ref_dict = G0_ref_dict)

        H0_act_vect  = reactions_data['H0 act'] 
        deltaH0_vect = reactions_data['deltaH0']
        
        G0_act_vect  = reactions_data['G0 act'] 
        deltaG0_vect = reactions_data['deltaG0']

        print('', file = f)

        print(head_spec.format('{} adsorbates'.format(site_name),
                              'H0 [eV]'  ,
                              'S0 [eV/K]',
                              'G0 [eV]'  ), file = f)

        for s in range(cat.n_species):
            
            spec = cat.species_names[s]
            
            print(string_spec.format(spec         ,
                                     H0_dict[spec],
                                     S0_dict[spec],
                                     G0_dict[spec]), file = f)

        print('', file = f)

        print(head_reax.format('{} reactions'.format(site_name),
                              'H0 [eV]'     ,
                              'S0 [eV/K]'   ,
                              'G0 [eV]'     ,
                              'H0 act [eV]' ,
                              'deltaH0 [eV]',
                              'G0 act [eV]' ,
                              'deltaG0 [eV]'), file = f)

        for s in range(TS.n_species):
            
            step = cat.reaction_equation(s)
            name = step.replace(' ', '')
            
            spec = TS.species_names[s]
            
            print(string_reax.format(name           ,
                                     H0_dict[spec]  ,
                                     S0_dict[spec]  ,
                                     G0_dict[spec]  ,
                                     H0_act_vect[s] ,
                                     deltaH0_vect[s],
                                     G0_act_vect[s] ,
                                     deltaG0_vect[s]), file = f)

    f.close()

################################################################################
# CALCUlATE REACTION PATHS
################################################################################

energy_type = 'G0'

plot_custom_path = False

plot_CO_O_path = True
plot_COOH_path = True
plot_OHOH_path = True

plot_C_O_path = False
plot_COH_path = False

if plot_paths is True:

    # CUSTOM
    
    if reaction == 'WGS':
        reaction_sequence = [3, 7, 8]
    elif reaction == 'revWGS':
        reaction_sequence = [2, 6]
    
    spec_sequences = {}
    H0_sequences = {}
    G0_sequences = {}

    for i in range(len(active_phases)):
    
        site_name = site_names[i]
    
        cat = cat_list[i]
        cat.TP = gas.TP
    
        TS = TS_list[i]
        TS.TP = gas.TP
    
        species_dict = {}
    
        species_dict['CO']  = 1. if reaction == 'WGS' else 0.
        species_dict['H2O'] = 1. if reaction == 'WGS' else 0.
        species_dict['CO2'] = 1. if reaction == 'revWGS' else 0.
        species_dict['H2']  = 1. if reaction == 'revWGS' else 0.
        species_dict['Rh('+site_name+')'] = 2.
    
        react_path_data = get_energy_path(gas               = gas              ,
                                          cat               = cat              ,
                                          TS                = TS               ,
                                          species_dict      = species_dict     ,
                                          reaction_sequence = reaction_sequence,
                                          units             = units            ,
                                          print_data        = False            )
    
        spec_sequences[site_name] = react_path_data['spec sequence']
        H0_sequences[site_name]   = react_path_data['H0 sequence']
        G0_sequences[site_name]   = react_path_data['G0 sequence']

    if energy_type == 'H0':
        E_sequences = H0_sequences
        y_min, y_max = -3., 2.
    elif energy_type == 'G0':
        E_sequences = G0_sequences
        y_min, y_max = -3., 2.

    if plot_custom_path is True:

        fig_num += 1
    
        plot_energy_paths(site_names     = site_names    ,
                          E_sequences    = E_sequences   ,
                          spec_sequences = spec_sequences,
                          colors         = colors        ,
                          y_min          = y_min         ,
                          y_max          = y_max         ,
                          show_plot      = False         ,
                          fig_num        = fig_num       )

    # CO* + O* <=> CO2**
    
    if reaction == 'WGS':
        reaction_sequence = [1, 3, 7, 8, -6, -2, -4]
    elif reaction == 'revWGS':
        reaction_sequence = [4, 2, 6, -8, -7, -3, -1]
    
    spec_sequences = {}
    H0_sequences = {}
    G0_sequences = {}

    for i in range(len(active_phases)):
    
        site_name = site_names[i]
    
        cat = cat_list[i]
        cat.TP = gas.TP
    
        TS = TS_list[i]
        TS.TP = gas.TP
    
        species_dict = {}
    
        species_dict['CO']  = 1. if reaction == 'WGS' else 0.
        species_dict['H2O'] = 1. if reaction == 'WGS' else 0.
        species_dict['CO2'] = 1. if reaction == 'revWGS' else 0.
        species_dict['H2']  = 1. if reaction == 'revWGS' else 0.
        species_dict['Rh('+site_name+')'] = 4.
    
        react_path_data = get_energy_path(gas               = gas              ,
                                          cat               = cat              ,
                                          TS                = TS               ,
                                          species_dict      = species_dict     ,
                                          reaction_sequence = reaction_sequence,
                                          units             = units            ,
                                          print_data        = False            )
    
        spec_sequences[site_name] = react_path_data['spec sequence']
        H0_sequences[site_name]   = react_path_data['H0 sequence']
        G0_sequences[site_name]   = react_path_data['G0 sequence']

    if energy_type == 'H0':
        E_sequences = H0_sequences
        y_min, y_max = -4., 1.
    elif energy_type == 'G0':
        E_sequences = G0_sequences
        y_min, y_max = -2., 3.

    if plot_CO_O_path is True:

        fig_num += 1
    
        plot_energy_paths(site_names     = site_names    ,
                          E_sequences    = E_sequences   ,
                          spec_sequences = spec_sequences,
                          colors         = colors        ,
                          y_min          = y_min         ,
                          y_max          = y_max         ,
                          show_plot      = False         ,
                          fig_num        = fig_num       )

    # CO* + OH* <=> COOH**
    
    if reaction == 'WGS':
        reaction_sequence = [1, 3, 7, 7, 11, 10, 12, -2, -4]
    elif reaction == 'revWGS':
        reaction_sequence = [4, 2, -12, -10, -11, -7, -7, -3, -1]
    
    spec_sequences = {}
    H0_sequences = {}
    G0_sequences = {}
    
    for i in range(len(active_phases)):
    
        site_name = site_names[i]
    
        cat = cat_list[i]
        cat.TP = gas.TP
    
        TS = TS_list[i]
        TS.TP = gas.TP
    
        species_dict = {}
    
        species_dict['CO']  = 1. if reaction == 'WGS' else 0.
        species_dict['H2O'] = 1. if reaction == 'WGS' else 0.
        species_dict['CO2'] = 1. if reaction == 'revWGS' else 0.
        species_dict['H2']  = 1. if reaction == 'revWGS' else 0.
        species_dict['Rh('+site_name+')'] = 4.
        species_dict['H2O('+site_name+')'] = 1.
    
        react_path_data = get_energy_path(gas               = gas              ,
                                          cat               = cat              ,
                                          TS                = TS               ,
                                          species_dict      = species_dict     ,
                                          reaction_sequence = reaction_sequence,
                                          units             = units            ,
                                          print_data        = False            )
    
        spec_sequences[site_name] = react_path_data['spec sequence']
        H0_sequences[site_name]   = react_path_data['H0 sequence']
        G0_sequences[site_name]   = react_path_data['G0 sequence']
    
    if energy_type == 'H0':
        E_sequences = H0_sequences
        y_min, y_max = -4., 1.
    elif energy_type == 'G0':
        E_sequences = G0_sequences
        y_min, y_max = -2., 3.

    if plot_COOH_path is True:

        fig_num += 1
    
        plot_energy_paths(site_names     = site_names    ,
                          E_sequences    = E_sequences   ,
                          spec_sequences = spec_sequences,
                          colors         = colors        ,
                          y_min          = y_min         ,
                          y_max          = y_max         ,
                          show_plot      = False         ,
                          fig_num        = fig_num       )
    
    # OH* + OH* <=> H2O* + O*
    
    if reaction == 'WGS':
        reaction_sequence = [1, 3, -6, -14, 8, 8, -2, -4]
    elif reaction == 'revWGS':
        reaction_sequence = [4, 2, -8, -8, 14, 6, -1, -3]
    
    spec_sequences = {}
    H0_sequences = {}
    G0_sequences = {}
    
    for i in range(len(active_phases)):
    
        site_name = site_names[i]
    
        cat = cat_list[i]
        cat.TP = gas.TP
    
        TS = TS_list[i]
        TS.TP = gas.TP
    
        species_dict = {}
    
        species_dict['CO']  = 1. if reaction == 'WGS' else 0.
        species_dict['H2O'] = 1. if reaction == 'WGS' else 0.
        species_dict['CO2'] = 1. if reaction == 'revWGS' else 0.
        species_dict['H2']  = 1. if reaction == 'revWGS' else 0.
        species_dict['Rh('+site_name+')'] = 4.
        species_dict['O('+site_name+')'] = 2.
    
        react_path_data = get_energy_path(gas               = gas              ,
                                          cat               = cat              ,
                                          TS                = TS               ,
                                          species_dict      = species_dict     ,
                                          reaction_sequence = reaction_sequence,
                                          units             = units            ,
                                          print_data        = False            )
    
        spec_sequences[site_name] = react_path_data['spec sequence']
        H0_sequences[site_name]   = react_path_data['H0 sequence']
        G0_sequences[site_name]   = react_path_data['G0 sequence']
    
    if energy_type == 'H0':
        E_sequences = H0_sequences
        y_min, y_max = -4., 1.
    elif energy_type == 'G0':
        E_sequences = G0_sequences
        y_min, y_max = -2., 3.

    if plot_OHOH_path is True:

        fig_num += 1
    
        plot_energy_paths(site_names     = site_names    ,
                          E_sequences    = E_sequences   ,
                          spec_sequences = spec_sequences,
                          colors         = colors        ,
                          y_min          = y_min         ,
                          y_max          = y_max         ,
                          show_plot      = False         ,
                          fig_num        = fig_num       )
    
if plot_paths is True and 'CH4' in gas_selected:
    
    # CO* <=> CO* + O*
    
    reaction_sequence = [16, -17, -18, -19, -20]
    
    spec_sequences = {}
    H0_sequences = {}
    G0_sequences = {}
    
    for i in range(len(active_phases)):
    
        site_name = site_names[i]
    
        cat = cat_list[i]
        cat.TP = gas.TP
    
        TS = TS_list[i]
        TS.TP = gas.TP
    
        species_dict = {}
    
        species_dict['CO('+site_name+')']  = 1.
        species_dict['H('+site_name+')'] = 6.
        species_dict['Rh('+site_name+')'] = 0.
    
        react_path_data = get_energy_path(gas               = gas              ,
                                          cat               = cat              ,
                                          TS                = TS               ,
                                          species_dict      = species_dict     ,
                                          reaction_sequence = reaction_sequence,
                                          units             = units            ,
                                          print_data        = False            )
    
        spec_sequences[site_name] = react_path_data['spec sequence']
        H0_sequences[site_name]   = react_path_data['H0 sequence']
        G0_sequences[site_name]   = react_path_data['G0 sequence']
    
    if energy_type == 'H0':
        E_sequences = H0_sequences
        y_min, y_max = -3., 2.
    elif energy_type == 'G0':
        E_sequences = G0_sequences
        y_min, y_max = -1., 4.
    
    if plot_C_O_path is True:

        fig_num += 1
    
        plot_energy_paths(site_names     = site_names    ,
                          E_sequences    = E_sequences   ,
                          spec_sequences = spec_sequences,
                          colors         = colors        ,
                          y_min          = y_min         ,
                          y_max          = y_max         ,
                          show_plot      = False         ,
                          fig_num        = fig_num       )

    # CO* + H* <=> COH* + *

    # TBD

################################################################################
# DEGREE OF RATE CONTROL
################################################################################

if DRC is True:

    if DRC_pts == 1:
        z_DRC_vect = [z_analysis]
    else:
        l = reactor_length
        z_DRC_vect = np.arange(0., l*(1+1e-4), l/DRC_pts)

    print('\n\n DEGREE OF RATE CONTROL ')

    DRC_multiplier = 1.05

    cstr.volume = cross_section*cstr_length_fix

    DRC_vect_dict = {}
    DRC_reduced   = {}
    for cat in cat_list:
        for j in range(cat.n_reactions):
            reaction_name = cat.reaction_equation(j)
            DRC_vect_dict[reaction_name] = []
    
    main_index = gas.species_names.index(main_product)
    
    for z_DRC in z_DRC_vect:

        index = np.argmin(np.abs(solution.z-z_DRC))
        
        residence_time = solution.z[index]/gas_velocity

        print('\nDistance       = {:.3e} m'.format(solution.z[index]))
        print('Residence time = {:.3e} s \n'.format(residence_time))
        
        gas.TDY = solution[index].TDY
        upstream.syncState()

        for i in range(len(surf_list)):

            #area = solution.cat_areas[index][i]*cstr.volume/cstr_volume
            area = cat_areas_vectors[index][i]*cstr.volume/cstr_volume

            surf_list[i].area = area

            cat = cat_list[i]
            TS = TS_list[i]
            site = site_names[i]

            #cat.coverages = solution.ads_coverages[index][i,:cat.n_species]
            cat.coverages = ads_coverages_matrices[index][i,:cat.n_species]

            if update_kin is True:

                update_kinetics_facet(gas         = gas        ,
                                      cat         = cat        ,
                                      TS          = TS         ,
                                      site        = site       ,
                                      lat_dict    = lat_dict   ,
                                      coeffs_dict = coeffs_dict)

        n_dot_zero = mass_flow_rate/sum(np.multiply(gas.molecular_weights,
                                                    gas.X))

        ni_dot_zero = np.multiply(n_dot_zero, gas.X)

        advance_to_steady_state(sim = sim)
    
        n_dot = mass_flow_rate/sum(np.multiply(gas.molecular_weights, gas.X))

        ni_dot = np.multiply(n_dot, gas.X)
    
        dn_original = ni_dot-ni_dot_zero
        dn_original = np.array([x if abs(x) > 0. else 1e-50
                                for x in dn_original]) 
    
        DRC_dict = {}
        step_selected_names = []

        for i in range(len(cat_list)):

            cat = cat_list[i]

            for j in range(cat.n_reactions):

                reaction_name = cat.reaction_equation(j)

                cat.set_multiplier(value = DRC_multiplier, i_reaction = j)
    
                gas.TDY = solution[index].TDY
                upstream.syncState()

                advance_to_steady_state(sim = sim)
    
                n_dot = (mass_flow_rate/sum(np.multiply(gas.molecular_weights,
                                                        gas.X)))

                ni_dot = np.multiply(n_dot, gas.X)
    
                dn_modified = ni_dot-ni_dot_zero
                dn_modified = np.array([x if abs(x) > 0. else 1e-50 
                                        for x in dn_modified])
    
                DRC_vector = np.multiply(dn_modified-dn_original,
                                         1/dn_original/(DRC_multiplier-1.))

                DRC_dict[reaction_name] = DRC_vector[main_index]
    
                cat.set_multiplier(value = 1.0, i_reaction = j)

                DRC_vect_dict[reaction_name] += [DRC_vector[main_index]]

                if j+1 in step_selected:
                    step_selected_names += [reaction_name]

        DRC_tot = 0.
        for step in DRC_dict:
            DRC_tot += DRC_dict[step]
            if abs(DRC_dict[step]) > DRC_thold:
                print('{0:50s} DRC({1}) = '.format(step, main_product),
                      '{:+.3f}'.format(DRC_dict[step]))
    
        print('\nDRC tot = {:.4f}'.format(DRC_tot))

    if print_outputs is True:
    
        filename = results_dir+'/DRCS_{:+06.2f}'.format(var_num)
        
        f = open(filename, 'w+')
    
        for step in [step for step in DRC_dict if step in step_selected_names]:

            name = step.replace(' ', '')

            print('{0:40s} {1:+9.6f}'.format(name, DRC_dict[step]), file = f)
    
        f.close()

    if DRC_pts > 1:

        """
        Plot the DRC along the reactor.
        """
        fig_num += 1
        
        fig = plt.figure(fig_num)
        
        ylabel = 'DRC [-]'
        
        ax = plt.subplot(111)
        
        plt.axis([0., max(z_vector), 0., 1.])
        plt.title('degree of rate control')
        for step in DRC_vect_dict:
            if max(DRC_vect_dict[step]) > DRC_thold:
                ax.plot(z_DRC_vect, DRC_vect_dict[step], label = step)
                DRC_reduced[step] = DRC_vect_dict[step]
        
        box = ax.get_position()
        ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
        
        plt.legend(loc = 'center left', bbox_to_anchor = (1, 0.5))
        
        plt.ylabel(ylabel, fontsize = label_size)
        
        plt.yticks(fontsize = tick_size)
        plt.xticks(fontsize = tick_size)
        
        if print_csv is True:

            """
            Write DRC profiles to file.
            """
        
            fileobj = open(DRC_file, 'w+')
            
            writer = csv.writer(fileobj)
            writer.writerow(['Distance [m]'] + [step for step in DRC_reduced])
            
            for i in range(len(z_DRC_vect)):

                writer.writerow([z_DRC_vect[i]] + [DRC_reduced[step][i] 
                                 for step in DRC_reduced])
        
            print('\nDRC data printed on file: {}'.format(DRC_file))
    
        fileobj.close()

################################################################################
# APPARENT ACTIVATION ENERGY
################################################################################

if Eapp is True:

    cstr.volume = cross_section*cstr_length_fix

    print('\n\n APPARENT ACTIVATION ENERGY \n')

    index = np.argmin(np.abs(solution.z-z_analysis))

    residence_time = solution.z[index]/gas_velocity

    print('Distance       = {:.3e} m'.format(solution.z[index]))
    print('Residence time = {:.3e} s \n'.format(residence_time))

    gas.TDY = solution[index].TDY
    upstream.syncState()

    area_spec_dict = {}

    for i in range(len(surf_list)):

        #area = solution.cat_areas[index][i]*cstr.volume/cstr_volume
        area = cat_areas_vectors[index][i]*cstr.volume/cstr_volume

        phase = active_phases[i]

        surf_list[i].area = area

        #cat.coverages = solution.ads_coverages[index][i,:cat.n_species]
        cat.coverages = ads_coverages_matrices[index][i,:cat.n_species]
        
        area_spec_dict[phase] = area/cstr.volume
    
    if update_kin is True:
    
        mari_list = []
        
        if gas['CO'].X[0] > x_min_CO:
            mari_list += ['CO']
        
        calculate_steady_state_coverages(gas         = gas        ,
                                         cat_list    = cat_list   ,
                                         TS_list     = TS_list    ,
                                         site_names  = site_names ,
                                         lat_dict    = lat_dict   ,
                                         coeffs_dict = coeffs_dict,
                                         free_site   = free_site  ,
                                         mari_list   = mari_list  ,
                                         mu_gas_dict = mu_gas_dict,
                                         options     = options    ,
                                         method      = method     )
        
        update_kinetics(gas           = gas          ,
                        cat_list      = cat_list     ,
                        TS_list       = TS_list      ,
                        site_names    = site_names   ,
                        active_phases = active_phases,
                        lat_dict      = lat_dict     ,
                        coeffs_dict   = coeffs_dict  )
    
    advance_to_steady_state(sim = sim)
    
    main_index = gas.species_names.index(main_product)
    
    prod_main_dict_old = {}
    prod_main_sum_old = 0.
    
    for i in range(len(cat_list)):
    
        cat = cat_list[i]
        phase = active_phases[i]
        
        area_spec = area_spec_dict[phase]
    
        prods = cat.get_net_production_rates(gas)*area_spec
    
        prod_main_dict_old[phase] = prods[main_index]

        prod_main_sum_old += prod_main_dict_old[phase]
    
    T, D, Y = solution[index].TDY

    T_new = T+deltaT
    
    water.TP = T_new, pressure
    
    heater_new = ct.Reservoir(water, name = 'heater_new')
    
    tcontrol_new = ct.Wall(cstr, heater_new, U = 1e+2)
    
    cstr.energy_enabled = True
    
    while abs(cstr.thermo.T-T_new) > 0.1:
        advance_to_steady_state(sim = sim)

    tcontrol_new.heat_transfer_coeff = 1e-6

    T, D, Y = solution[index].TDY

    gas.TDY = [T_new, D, Y]
    upstream.syncState()
    
    for cat in cat_list:
        cat.TP = T_new, pressure

    for TS in TS_list:
        TS.TP = T_new, pressure

    cstr.energy_enabled = False

    if update_kin is True:

        calculate_steady_state_coverages(gas         = gas        ,
                                         cat_list    = cat_list   ,
                                         TS_list     = TS_list    ,
                                         site_names  = site_names ,
                                         lat_dict    = lat_dict   ,
                                         coeffs_dict = coeffs_dict,
                                         free_site   = free_site  ,
                                         mari_list   = mari_list  ,
                                         mu_gas_dict = mu_gas_dict,
                                         options     = options    ,
                                         method      = method     )
    
        update_kinetics(gas           = gas          ,
                        cat_list      = cat_list     ,
                        TS_list       = TS_list      ,
                        site_names    = site_names   ,
                        active_phases = active_phases,
                        lat_dict      = lat_dict     ,
                        coeffs_dict   = coeffs_dict  )

    advance_to_steady_state(sim = sim)

    prod_main_sum_new = 0.
    prod_main_dict_new = {}
    
    for i in range(len(cat_list)):
    
        cat = cat_list[i]
        phase = active_phases[i]
        
        area_spec = area_spec_dict[phase]
    
        prods = cat.get_net_production_rates(gas)*area_spec
    
        prod_main_dict_new[phase] = prods[main_index]
    
        prod_main_sum_new += prod_main_dict_new[phase]

    Eapp = kB/eV*T**2*(np.log(prod_main_sum_new/prod_main_sum_old))/deltaT

    print('Eapp calculated from the SDMKM = {:+7.4f} eV'.format(Eapp))

    Eapp = 0.

    for i in range(len(cat_list)):
    
        cat = cat_list[i]
        cat.TP = gas.TP
        
        for n in range(cat.n_reactions):
    
            step = cat.reaction_equation(n)
    
            Eapp += G0_rel_TS_dict[step]*DRC_dict[step]

    print('Eapp calculated from the DRCS  = {:+7.4f} eV'.format(Eapp))

################################################################################
# MEASURE TIME END
################################################################################

plt.show()

if measure_time is True:
    stop = timeit.default_timer()
    print('\nExecution time = {0:6.3} s\n'.format(stop-start))

################################################################################
# END
################################################################################
