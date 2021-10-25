################################################################################
# Raffaele Cheula, LCCP, Politecnico di Milano, raffaele.cheula@polimi.it
################################################################################

from __future__ import absolute_import, division, print_function
import warnings
import numpy as np
import copy as cp
import cantera as ct
from scipy import optimize
from collections import OrderedDict
from shape.units import J, kmol, Rgas, atm, kB, hP, eV, molecule

################################################################################
# GET DELTAMU
################################################################################

def get_H0_0K_dict(phase, units = J/kmol):

    T, P = phase.TP

    phase.TP = 1e-10, P
    H0_0K_dict = get_enthalpies_dict(phase, units = units)

    phase.TP = T, P

    return H0_0K_dict

################################################################################
# GET DELTAMU
################################################################################

def get_deltamu_dict(phase, H0_0K_dict, units = J/kmol):

    deltamu = {}

    for spec in phase.species_names:
        deltamu[spec] = (phase[spec].chemical_potentials[0]*J/kmol/units -
                         H0_0K_dict[spec])

    return deltamu

################################################################################
# GET ENTHALPIES DICT
################################################################################

def get_enthalpies_dict(phase, units = J/kmol):

    H_dict = {}

    for s in range(phase.n_species):
        spec = phase.species_names[s]
        H_dict[spec] = (phase.standard_enthalpies_RT[s]*ct.gas_constant *
                        phase.T*J/kmol/units)

    return H_dict

################################################################################
# GET STD ENTROPIES DICT
################################################################################

def get_std_entropies_dict(phase, units = J/kmol):

    S0_dict = {}

    for s in range(phase.n_species):
        spec = phase.species_names[s]
        S0_dict[spec] = (phase.standard_entropies_R[s]*ct.gas_constant * 
                         J/kmol/units)

    return S0_dict

################################################################################
# GET STD GIBBS DICT
################################################################################

def get_std_gibbs_dict(phase, units = J/kmol):

    G0_dict = {}

    for s in range(phase.n_species):
        spec = phase.species_names[s]
        G0_dict[spec] = (phase.standard_gibbs_RT[s]*ct.gas_constant *
                         phase.T*J/kmol/units)

    return G0_dict

################################################################################
# GET ENTROPIES DICT
################################################################################

def get_entropies_dict(phase, units = J/kmol, xmin = 1e-100):

    S_dict = {}

    for s in range(phase.n_species):
        spec = phase.species_names[s]
        xs = phase[spec].X[0]+xmin
        S_dict[spec] = ((phase.standard_entropies_R[s]-np.log(xs)) *
                         ct.gas_constant*J/kmol/units)

    return S_dict

################################################################################
# GET GIBBS DICT
################################################################################

def get_gibbs_dict(phase, units = J/kmol, xmin = 1e-100):

    G_dict = {}

    for s in range(phase.n_species):
        spec = phase.species_names[s]
        xs = phase[spec].X[0]+xmin
        G_dict[spec] = ((phase.standard_gibbs_RT[s]+np.log(xs)) *
                         ct.gas_constant*phase.T*J/kmol/units)

    return G_dict

################################################################################
# MODIFY ADSORPTION RXN
################################################################################

def modify_adsorption_rxn(gas, cat, num, spec, stick = 1., beta = 0.5,
                          Eact = 0.):
    
    MW = gas[spec].molecular_weights[0]
    k0 = stick*np.sqrt(Rgas/(2*np.pi*MW))/cat.site_density
    
    modify_rxn(cat = cat, num = num, k0 = k0, beta = beta, Eact = Eact)

################################################################################
# MODIFY DESORPTION RXN
################################################################################

def modify_desorption_rxn(gas, cat, num, spec, stick = 1., beta = -0.5,
                          Eact = 0.):

    MW = gas[spec].molecular_weights[0]
    k0 = stick/np.sqrt(2*np.pi*MW*Rgas)*cat.P/cat.site_density
    
    stoich = (cat.product_stoich_coeffs()-cat.reactant_stoich_coeffs())[:,num]
    k0 += np.dot(stoich, np.append(gas.standard_entropies_R,
                                   cat.standard_entropies_R))
    Eact = np.dot(stoich, np.append(gas.standard_enthalpies_RT,
                                    cat.standard_enthalpies_RT))
    
    modify_rxn(cat = cat, num = num, k0 = k0, beta = beta, Eact = Eact)

################################################################################
# MODIFY RXN
################################################################################

def modify_rxn(cat, num, k0, beta, Eact):

    rxn = ct.Reaction.fromCti("""surface_reaction('{0}',
    [{1:.6E}, {2:.1f}, {3:.1f}])""".format(cat.reaction(num), k0, beta, Eact))
    cat.modify_reaction(num, rxn)

################################################################################
# PRINT RPA DETAILS
################################################################################

def print_RPA_details(data, fileobj = None):

    lines = data.splitlines()
    
    reactants = OrderedDict()
    r_cons = OrderedDict()
    
    for line in lines[2:]:
        entries = line.split(' ')
        if float(entries[2]) > 0.:
            r_for, r_rev = float(entries[2]), -float(entries[3])
            if r_for > r_rev:
                reac, prod = entries[0], entries[1]
            else:
                reac, prod = entries[1], entries[0]
                r_for, r_rev = r_rev, r_for
            try: reactants[reac] += [[prod, r_for, r_rev]]
            except KeyError: reactants[reac] = [[prod, r_for, r_rev]]
            try: r_cons[reac] += r_for-r_rev
            except KeyError: r_cons[reac] = r_for-r_rev

    for reac in reactants:
        first = True
        for step in reactants[reac]:
            prod, r_for, r_rev = step
            phi = r_for/(r_for+r_rev)
            try: r_perc = (r_for-r_rev)/r_cons[reac]*100
            except ZeroDivisionError: r_perc = 0.
            if first is True:
                string = '{:12s} -> '.format(reac)
            else:
                string = '{:12s} -> '.format('')
            format_str = '{0:12s} {1:6.2f} %   r = {2:4.2E}   phi = {3:4.2f}'
            string += format_str.format(prod, r_perc, r_for-r_rev, phi)
            print(string, file = fileobj)
            first = False

################################################################################
# UPDATE KINETICS
################################################################################

def update_kinetics(gas, cat_list, TS_list, site_names, active_phases,
                    lat_dict, coeffs_dict, T_low = 200., T_high = 3000.,
                    P_ref = 1*atm):

    for i in range(len(active_phases)):

        phase = active_phases[i]

        cat = cat_list[i]
        cat.TP = gas.TP

        TS = TS_list[i]
        TS.TP = gas.TP

        site = site_names[i]

        update_kinetics_facet(gas, cat, TS, site, lat_dict, coeffs_dict)

################################################################################
# UPDATE KINETICS FACET
################################################################################

def update_kinetics_facet(gas, cat, TS, site, lat_dict, coeffs_dict, 
                          T_low = 200., T_high = 3000., P_ref = 1*atm):

    lat_list = lat_dict[site]

    corr_dict = {}

    lat_coverages = {}

    for spec in cat.species_names:

        corr = 0.

        for lat in [lat for lat in lat_list if lat.ads_name == spec]:

            try:
                coverage = lat_coverages[lat.cov_dep]
            except:
                coverage = cat[lat.cov_dep].coverages[0]
                lat_coverages[lat.cov_dep] = coverage

            corr += lat.get_coeff_correction(coverage = coverage)

        corr_dict[spec] = corr

        index = cat.species_names.index(spec)
        species = cat.species(index)

        coeffs = cp.deepcopy(coeffs_dict[spec])
    
        coeffs[6]  += corr
        coeffs[13] += corr

        species.thermo = ct.NasaPoly2(T_low  = T_low ,
                                      T_high = T_high,
                                      P_ref  = P_ref ,
                                      coeffs = coeffs)

        cat.modify_species(index, species)

    for n in [n for n in range(cat.n_reactions) 
              if not cat.reaction(n).is_sticking_coefficient]:

        corr = 0.

        spec = TS.species_names[n]
    
        #for react in [react for react in cat.reaction(n).reactants
        #              if react in cat.species_names]:
        #    corr += corr_dict[react]
    
        for react in [react for react in cat.reaction(n).reactants
                      if react in cat.species_names]:
            corr += corr_dict[react]/2.
    
        for prod in [prod for prod in cat.reaction(n).products
                     if prod in cat.species_names]:
            corr += corr_dict[prod]/2.
    
        index = TS.species_names.index(spec)
        species = TS.species(index)
    
        coeffs = cp.deepcopy(coeffs_dict[spec])
    
        coeffs[6]  += corr
        coeffs[13] += corr
    
        species.thermo = ct.NasaPoly2(T_low  = T_low ,
                                      T_high = T_high,
                                      P_ref  = P_ref ,
                                      coeffs = coeffs)
    
        TS.modify_species(index, species)

    G0_dict = get_std_gibbs_dict(gas, units = J/kmol)
    G0_dict.update(get_std_gibbs_dict(cat, units = J/kmol))
    G0_dict.update(get_std_gibbs_dict(TS, units = J/kmol))

    for n in [n for n in range(cat.n_reactions) 
              if not cat.reaction(n).is_sticking_coefficient]:

        spec = TS.species_names[n]
        G0_act = G0_dict[spec]

        reactants = cat.reaction(n).reactants
        products  = cat.reaction(n).products

        reactants_gas = []
        reactants_ads = []
        
        for spec in reactants:
            if spec in gas.species_names:
                reactants_gas += [spec]*int(reactants[spec])
            elif spec in cat.species_names:
                reactants_ads += [spec]*int(reactants[spec])

        for react in reactants:
            G0_act -= G0_dict[react]*reactants[react]

        k = kB*cat.T/hP*np.exp(-G0_act/Rgas/cat.T)

        k *= (Rgas*cat.T/cat.P)**(len(reactants_gas))
        k *= cat.site_density**(1-len(reactants_ads))

        rxn = ct.InterfaceReaction(reactants = reactants,
                                   products  = products )

        rxn.rate = ct.Arrhenius(A = k, b = 0., E = 0.)

        cat.modify_reaction(n, rxn)

################################################################################
# ADVANCE AND UPDATE KINETICS
################################################################################

def advance_and_update_kinetics(sim, gas, cat_list, TS_list, site_names,
                                active_phases, lat_dict, coeffs_dict,
                                timestep_zero = 1e-4, timestep_min = 1e-8, 
                                tot_err_thold = 1e-8, derr_dt_thold = 1e-2,
                                afterfail_fac = 0.1, print_err = False):

    sim.reinitialize()

    results_dict = {}
    for spec in gas.species_names:
        results_dict[spec] = gas[spec].X[0]
    for cat in cat_list:
        for ads in cat.species_names:
            results_dict[ads] = cat[ads].coverages[0]

    timestep = timestep_zero

    err = tot_err_thold*10
    derr_dt = err/timestep

    time = 0.
    count = 0
    while err > tot_err_thold or derr_dt > derr_dt_thold:

        try:

            time += timestep

            update_kinetics(gas           = gas          ,
                            cat_list      = cat_list     ,
                            TS_list       = TS_list      ,
                            site_names    = site_names   ,
                            active_phases = active_phases,
                            lat_dict      = lat_dict     ,
                            coeffs_dict   = coeffs_dict  )

            sim.advance(time)

            string = ''

            for cat in cat_list:
                string += '{:.6f} '.format(cat.coverages[0])
    
            err = 0.
            for spec in gas.species_names:
                err += abs(gas[spec].X[0]-results_dict[spec])
                results_dict[spec] = gas[spec].X[0]
            for cat in cat_list:
                for ads in cat.species_names:
                    err += abs(cat[ads].coverages[0]-results_dict[ads])
                    results_dict[ads] = cat[ads].coverages[0]
    
            derr_dt = err/timestep
    
            if print_err is True and not count % 100:
                string +='   tot_err = {:.4e}'.format(err)
                string +='   derr_dt = {:.4e}'.format(derr_dt)
                print(string)
    
            count += 1

        except: 

            if timestep > timestep_min:
                timestep *= afterfail_fac
            else:
                timestep += timestep_min

            if print_err is True:
                print('\n timestep = {}\n'.format(timestep))

            count = 0

    string = '\n steady state free sites (CALCULATED): \n'
    for cat in cat_list:
        string += '{:.6f} '.format(cat.coverages[0])
    string +='   tot_err = {:.4e}'.format(err)
    string +='   derr_dt = {:.4e}'.format(derr_dt)
    print(string)

    advance_to_steady_state(sim = sim, n_max_err = 10000)

    err = 0.
    for spec in gas.species_names:
        err += abs(gas[spec].X[0]-results_dict[spec])
        results_dict[spec] = gas[spec].X[0]
    for cat in cat_list:
        for ads in cat.species_names:
            err += abs(cat[ads].coverages[0]-results_dict[ads])
            results_dict[ads] = cat[ads].coverages[0]

    string = '\n steady state free sites (CHECK): \n'
    for cat in cat_list:
        string += '{:.6f} '.format(cat.coverages[0])
    string +='   error   = {:.4e}'.format(err)
    print(string)

################################################################################
# ADVANCE TO STEADY STATE
################################################################################

def advance_to_steady_state(sim, n_max_err = 10000):

    finish = False
    count = 0
    while finish is False:
        try:
            sim.reinitialize()
            sim.advance_to_steady_state()
            finish = True
        except:
            count += 1
            pass
        if count > n_max_err:
            print('\n REACHED MAXIMUM AMOUNT OF ERRORS \n')
            break

################################################################################
# REACTOR ODE
################################################################################

class CstrLines:

    def __init__(self, gas, cat_list, surf_list, TS_list, site_names,
                 active_phases, lat_dict, coeffs_dict, mass_flow_rate,
                 cstr_volume):

        self.gas            = gas
        self.cat_list       = cat_list
        self.surf_list      = surf_list
        self.TS_list        = TS_list
        self.site_names     = site_names
        self.active_phases  = active_phases
        self.lat_dict       = lat_dict
        self.coeffs_dict    = coeffs_dict
        self.mass_flow_rate = mass_flow_rate
        self.cstr_volume    = cstr_volume

        self.Y_inlet = gas.Y

    def __call__(self, t, A):

        string = '  {0:12f}'.format(t)
        for spec in self.gas.species_names:
            string += '  {0:10f}'.format(*self.gas[spec].X)
        print(string)

        n = self.gas.n_species
        m = self.gas.n_species

        self.gas.set_unnormalized_mass_fractions(A[:n])

        dAdt = np.zeros(n)
        dYdt = np.zeros(n)

        for i in range(len(self.cat_list)):

            cat = self.cat_list[i]
            area = self.surf_list[i].area

            n += cat.n_species

            cat.set_unnormalized_coverages(A[m:n])

            dtetadt = cat.get_net_production_rates(cat)/cat.site_density*1e5

            dAdt = np.hstack((dAdt, dtetadt))

            dYdt += ((self.mass_flow_rate*(self.Y_inlet-self.gas.Y)
                     +np.multiply(cat.get_net_production_rates(self.gas),
                                  self.gas.molecular_weights)*area)
                     /self.cstr_volume/self.gas.density)

            m += cat.n_species
            
        n = self.gas.n_species

        dAdt[:n] = dYdt

        update_kinetics(gas           = self.gas          ,
                        cat_list      = self.cat_list     ,
                        TS_list       = self.TS_list      ,
                        site_names    = self.site_names   ,
                        active_phases = self.active_phases,
                        lat_dict      = self.lat_dict     ,
                        coeffs_dict   = self.coeffs_dict  )

        return dAdt

################################################################################
# REACTION PATH ANALYSIS
################################################################################

def reaction_path_analysis(gas, cat_list, surf_list, filename = None,
                           perc_thold = 0., filename_select = None,
                           react_select = [], spec_select = [],
                           units = eV/molecule):

    if filename is not None:
        fileobj = open(filename, 'w+')
    else:
        fileobj = None

    phi = OrderedDict()
    deltaG = OrderedDict()
    react_rates = OrderedDict()
    react_names = OrderedDict()
    react_rates_for = OrderedDict()
    react_rates_rev = OrderedDict()

    for i in range(gas.n_reactions):
        react = gas.reaction(i)
        r_net = gas.net_rates_of_progress[i]
        r_for = gas.forward_rates_of_progress[i]
        r_rev = gas.reverse_rates_of_progress[i]
        if r_for < r_rev:
            r_for, r_rev = r_rev, r_for
        phi[react] = r_for/(r_for+r_rev+1e-80)
        react_rates[react] = r_net
        react_rates_for[react] = r_for
        react_rates_rev[react] = r_rev
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            deltaG[react] = -Rgas*gas.T*np.log(r_for/r_rev)/units

    for c in range(len(cat_list)):
        cat = cat_list[c]
        area = surf_list[c].area
        for i in range(cat.n_reactions):
            react = cat.reaction(i)
            r_net = cat.net_rates_of_progress[i]*area
            r_for = cat.forward_rates_of_progress[i]*area
            r_rev = cat.reverse_rates_of_progress[i]*area
            if r_for < r_rev:
                r_for, r_rev = r_rev, r_for
            phi[react] = r_for/(r_for+r_rev+1e-50)
            react_rates[react] = r_net
            react_names[react] = cat.reaction_equation(i)
            react_rates_for[react] = r_for
            react_rates_rev[react] = r_rev
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                deltaG[react] = -Rgas*cat.T*np.log(r_for/r_rev)/units

    rates_dict = OrderedDict()

    for spec in gas.species_names:
        rates_dict[spec] = OrderedDict()

    for react in react_rates:

        for spec in react.reactants:
            nu = -react.reactants[spec]
            try: 
                rates_dict[spec][react] = nu*react_rates[react]
            except:
                rates_dict[spec] = OrderedDict()
                rates_dict[spec][react] = nu*react_rates[react]

        for spec in react.products:
            nu = react.products[spec]
            try: 
                rates_dict[spec][react] = nu*react_rates[react]
            except:
                rates_dict[spec] = OrderedDict()
                rates_dict[spec][react] = nu*react_rates[react]

    for spec in rates_dict:

        prod_sum = sum([rates_dict[spec][react] for react in rates_dict[spec]
                        if rates_dict[spec][react] > 0.])

        cons_sum = sum([rates_dict[spec][react] for react in rates_dict[spec]
                        if rates_dict[spec][react] < 0.])

        string = '{0} {1} {2}\n\n'.format('#'*5, spec, '#'*(83-len(spec)))

        for react in sorted(rates_dict[spec], 
                            key = lambda react: rates_dict[spec][react]):

            if rates_dict[spec][react] >= 0.:
                prod_perc = rates_dict[spec][react]/(prod_sum+1e-50)*100
            else:
                prod_perc = -rates_dict[spec][react]/(cons_sum+1e-50)*100

            if abs(prod_perc) > perc_thold:

                string += '{:50s}'.format(react_names[react])
                string += 'r = {:+10.4e}'.format(rates_dict[spec][react])
                string += '{:+9.2f} %'.format(prod_perc)
                string += '  phi = {:6.4f}'.format(phi[react])
                string += '  r_for = {:+10.4e}'.format(react_rates_for[react])
                string += '  r_rev = {:+10.4e}'.format(react_rates_rev[react])
                string += '  deltaG = {:+6.4f} eV'.format(deltaG[react])
                string += '\n'

        print(string, file = fileobj)

    if filename is not None:
        fileobj.close()
        print('RPA data printed on file: {}'.format(filename))

    string_select = ''

    for spec in spec_select:
        for cat in cat_list:
            for i in react_select:
                react_eq = cat.reaction_equation(i-1)
                for react in rates_dict[spec]:
                    if react_eq == react.equation:
                        react_rate = rates_dict[spec][react]
                        react_name = react_names[react].replace(' ', '')
                        string_select += '{:50s}'.format(react_name)
                        string_select += '{:+10.4e}'.format(react_rate)
                        string_select += '\n'

    if filename_select is not None:

        fileobj_select = open(filename_select, 'w+')

        print(string_select, file = fileobj_select)

        fileobj.close()

################################################################################
# STEADY STATE COVERAGES
################################################################################

def steady_state_coverages(coverages_list, gas, cat_list, TS_list, site_names,
                           active_phases, lat_dict, coeffs_dict, free_site, 
                           mari_list, mu_gas_dict, epsi = 1.e-20,
                           scale_factor = 1.e10):

    mari_coverages = [{}]*len(active_phases)
    coverage_free = [1.]*len(active_phases)

    err = [0.]*len(active_phases)*len(mari_list)

    for i in range(len(active_phases)):

        cat = cat_list[i]
        site = site_names[i]

        coverage_tot = 0.
        mari_coverages_i = {}

        for j in range(len(mari_list)):

            spec = mari_list[j]

            coverage = coverages_list[i*len(mari_list)+j]

            if coverage >= 1.-epsi:
                err[i*len(mari_list)+j] += scale_factor*(coverage-1.+epsi)
                coverage = 1.-epsi
            
            elif coverage <= 0.+epsi:
                err[i*len(mari_list)+j] += scale_factor*(coverage-0.-epsi)
                coverage = 0.+epsi

            mari_coverages_i[spec] = coverage
            coverage_tot += coverage
            
        if coverage_tot >= 1.-epsi:
            for j in range(len(mari_list)):
                err[i*len(mari_list)+j] += scale_factor*(coverage_tot-1.+epsi)
            coverage_tot = 1.-epsi

        coverage_free[i] = 1.-coverage_tot

        mari_coverages[i] = mari_coverages_i

        coverages = np.zeros(cat.n_species)

        for j in range(cat.n_species):

            if cat.species_names[j] == free_site+'('+site+')':
                coverages[j] = coverage_free[i]

            else:
                for spec in mari_list:
                    if cat.species_names[j] == spec+'('+site+')':
                        coverages[j] += mari_coverages[i][spec]

        cat.coverages = coverages

    update_kinetics(gas           = gas          ,
                    cat_list      = cat_list     ,
                    TS_list       = TS_list      ,
                    site_names    = site_names   ,
                    active_phases = active_phases,
                    lat_dict      = lat_dict     ,
                    coeffs_dict   = coeffs_dict  )

    for i in range(len(active_phases)):

        cat = cat_list[i]
        site = site_names[i]

        for j in range(len(mari_list)):

            spec = mari_list[j]

            coverage = mari_coverages[i][spec]

            mu0_ads = (cat[spec+'('+site+')'].standard_gibbs_RT[0] *
                       ct.gas_constant*cat.T)
            mu_gas = mu_gas_dict[spec](gas)

            if err[i*len(mari_list)+j] == 0.:
                err[i*len(mari_list)+j] = (mu0_ads-mu_gas+
                    Rgas*cat.T*np.log(coverage/coverage_free[i]))

    return err

################################################################################
# STEADY STATE COVERAGE
################################################################################

def steady_state_coverage(coverages_list, gas, cat, TS, site,
                          lat_dict, coeffs_dict, free_site, mari_list, 
                          mu_gas_dict, epsi = 1.e-20, scale_factor = 1.e12):

    err = [0.]*len(mari_list)

    coverage_tot = 0.
    mari_coverages = {}

    for j in range(len(mari_list)):

        spec = mari_list[j]

        coverage = coverages_list[j]

        if coverage >= 1.-epsi:
            err[j] += scale_factor*(coverage-1.+epsi)
            coverage = 1.-epsi
        
        elif coverage <= 0.+epsi:
            err[j] += scale_factor*(coverage-0.-epsi)
            coverage = 0.+epsi

        mari_coverages[spec] = coverage
        coverage_tot += coverage
        
    if coverage_tot >= 1.-epsi:
        for j in range(len(mari_list)):
            err[j] += scale_factor*(coverage_tot-1.+epsi)/len(mari_list)
        coverage_tot = 1.-epsi

    coverage_free = 1.-coverage_tot

    coverages = np.zeros(cat.n_species)

    for j in range(cat.n_species):

        if cat.species_names[j] == free_site+'('+site+')':
            coverages[j] = coverage_free

        else:
            for spec in mari_list:
                if cat.species_names[j] == spec+'('+site+')':
                    coverages[j] += mari_coverages[spec]

    cat.coverages = coverages

    update_kinetics_facet(gas         = gas        ,
                          cat         = cat        ,
                          TS          = TS         ,
                          site        = site       ,
                          lat_dict    = lat_dict   ,
                          coeffs_dict = coeffs_dict)

    for j in range(len(mari_list)):

        spec = mari_list[j]

        coverage = mari_coverages[spec]

        mu0_ads = (cat[spec+'('+site+')'].standard_gibbs_RT[0] * 
                   ct.gas_constant*cat.T)
        mu_gas = mu_gas_dict[spec](gas)

        if err[j] == 0.:
            err[j] = (mu0_ads-mu_gas+Rgas*cat.T*np.log(coverage/coverage_free))

    return err

################################################################################
# CALCULATE STEADY STATE COVERAGE
################################################################################

def calculate_steady_state_coverages(gas, cat_list, TS_list, site_names, 
                                     lat_dict, coeffs_dict, free_site,
                                     mari_list, mu_gas_dict, options, method):

    for i in range(len(cat_list)):
    
        coverages_list = [0.2]*len(mari_list)
    
        cat = cat_list[i]
        cat.TP = gas.TP
    
        TS = TS_list[i]
        TS.TP = gas.TP
    
        site = site_names[i]
    
        sol = optimize.root(steady_state_coverage, coverages_list,
                            options = options, method = method, args = (gas,
                            cat, TS, site, lat_dict, coeffs_dict,
                            free_site, mari_list, mu_gas_dict))

        if sol.success is False:

            print(sol)

################################################################################
# MODIFY E ACT
################################################################################

def modify_activation_energies(cat, TS, corr, coeffs_dict, T_low = 200.,
                               T_high = 3000., P_ref = 1*atm):

    for n in [n for n in range(cat.n_reactions) 
              if not cat.reaction(n).is_sticking_coefficient]:

        spec = TS.species_names[n]
    
        index = TS.species_names.index(spec)
        species = TS.species(index)
    
        coeffs = coeffs_dict[spec]
    
        coeffs[6]  += corr
        coeffs[13] += corr
    
        species.thermo = ct.NasaPoly2(T_low  = T_low ,
                                      T_high = T_high,
                                      P_ref  = P_ref ,
                                      coeffs = coeffs)
    
        TS.modify_species(index, species)

        coeffs_dict[spec] = coeffs

    return coeffs_dict

################################################################################
# GET REACTIONS DATA
################################################################################

def get_reactions_data(gas, cat, TS, units = eV/molecule, G0_ref_dict = {},
                       get_G = False):

    H0_dict = get_enthalpies_dict(gas, units = units)
    H0_dict.update(get_enthalpies_dict(cat, units = units))
    H0_dict.update(get_enthalpies_dict(TS, units = units))

    G0_dict = get_std_gibbs_dict(gas, units = units)
    G0_dict.update(get_std_gibbs_dict(cat, units = units))
    G0_dict.update(get_std_gibbs_dict(TS, units = units))

    H0_rel_vect = []
    H0_act_vect = []
    deltaH0_vect = []

    G0_rel_vect = []
    G0_act_vect = []
    deltaG0_vect = []

    if get_G is True:
        G_rel_vect = []
        G_act_vect = []
        deltaG_vect = []

    for n in range(cat.n_reactions):

        reactants = cat.reaction(n).reactants

        reactants_gas = []
        reactants_ads = []

        for spec in reactants:
            if spec in gas.species_names:
                reactants_gas += [spec]*int(reactants[spec])
            elif spec in cat.species_names:
                reactants_ads += [spec]*int(reactants[spec])

        products = cat.reaction(n).products

        products_gas = []
        products_ads = []

        for spec in products:
            if spec in gas.species_names:
                products_gas += [spec]*int(products[spec])
            elif spec in cat.species_names:
                products_ads += [spec]*int(products[spec])

        ts = TS.species_names[n]
        
        k_for = cat.forward_rate_constants[n]
        k_for /= (Rgas*cat.T/cat.P)**(len(reactants_gas))
        k_for /= cat.site_density**(1-len(reactants_ads))

        k_rev = cat.reverse_rate_constants[n]
        k_rev /= (Rgas*cat.T/cat.P)**(len(products_gas))
        k_rev /= cat.site_density**(1-len(products_ads))
        
        if cat.reaction(n).is_sticking_coefficient is True:

            H0_dict[ts] = cat.reaction(n).rate.activation_energy/units

            for spec in reactants:
                H0_dict[ts] += H0_dict[spec]*reactants[spec]

            G0_dict[ts] = -Rgas*cat.T*np.log(k_for/(kB*cat.T/hP))/units
        
            for spec in reactants:
                G0_dict[ts] += G0_dict[spec]*reactants[spec]
        
        H0_act = H0_dict[ts]
        deltaH0 = 0.

        for spec in reactants:
            H0_act -= H0_dict[spec]*reactants[spec]
            deltaH0 -= H0_dict[spec]*reactants[spec]
        for spec in products:
            deltaH0 += H0_dict[spec]*products[spec]
            
        H0_act_vect += [H0_act]
        deltaH0_vect += [deltaH0]

        G0_act = G0_dict[ts]
        deltaG0 = 0.

        for spec in reactants:
            G0_act -= G0_dict[spec]*reactants[spec]
            deltaG0 -= G0_dict[spec]*reactants[spec]
        for spec in products:
            deltaG0 += G0_dict[spec]*products[spec]

        G0_act_vect += [G0_act]
        deltaG0_vect += [deltaG0]

        G0_ref = G0_dict[ts]

        for element in G0_ref_dict:
            for spec in reactants_gas:
                G0_ref -= gas.n_atoms(spec, element)*G0_ref_dict[element]
            for spec in reactants_ads:
                G0_ref -= cat.n_atoms(spec, element)*G0_ref_dict[element]

        G0_rel_vect += [G0_ref]

        if get_G is True:

            r_for = cat.forward_rates_of_progress[n]
            r_rev = cat.reverse_rates_of_progress[n]
            
            G_act = -Rgas*cat.T*np.log(r_for/(kB*cat.T/hP))/units
            deltaG = -Rgas*cat.T*np.log(r_for/r_rev)/units
            
            G_act_vect += [G_act]
            deltaG_vect += [deltaG]

    reactions_data = {}

    reactions_data['H0 act']  = H0_act_vect
    reactions_data['deltaH0'] = deltaH0_vect
 
    reactions_data['G0 act']  = G0_act_vect
    reactions_data['deltaG0'] = deltaG0_vect
    reactions_data['G0 rel']  = G0_rel_vect
 
    if get_G is True:
        reactions_data['G act']  = G_act_vect
        reactions_data['deltaG'] = deltaG_vect

    return reactions_data

################################################################################
# PRINT ENERGY PATH
################################################################################

def get_energy_path(gas, cat, TS, species_dict, reaction_sequence,
                    print_data = True, reactions_data = None,
                    units = eV/molecule):

    H0_tot = 0.
    G0_tot = 0.
    #G_tot = 0.

    H0_sequence = [H0_tot, H0_tot]
    G0_sequence = [G0_tot, G0_tot]
    #G_sequence  = [G_tot, G_tot]
    spec_sequence = [cp.deepcopy(species_dict)]

    if reactions_data is None:
        reactions_data = get_reactions_data(gas, cat, TS, units = eV/molecule)

    H0_act_vect  = reactions_data['H0 act']
    deltaH0_vect = reactions_data['deltaH0']
    G0_act_vect  = reactions_data['G0 act']
    deltaG0_vect = reactions_data['deltaG0']
    #G_act_vect   = reactions_data['G act']
    #deltaG_vect  = reactions_data['deltaG']

    if print_data is True:
        print('\nspecies'.ljust(60)+'H0'.ljust(11)+'G0'.ljust(11) + 
            'H0 act'.ljust(11)+'G0 act')

    string = ''
    for spec in [spec for spec in species_dict if species_dict[spec] > 0]:
        string += '{0:2.0f} {1:8s} '.format(species_dict[spec], spec)
    string = string.ljust(58)
    string += '{0:+7.3f} eV {1:+7.3f} eV '.format(H0_tot, G0_tot)

    if print_data is True:
        print(string)

    for m in reaction_sequence:

        if m > 0:
            n = m-1
            reactants = cat.reaction(n).reactants
            products  = cat.reaction(n).products
            H0_act    = H0_act_vect[n]
            deltaH0   = deltaH0_vect[n]
            G0_act    = G0_act_vect[n]
            deltaG0   = deltaG0_vect[n]
            #G_act     = G_act_vect[n]
            #deltaG    = deltaG_vect[n]
        else:
            n = -m-1
            products  = cat.reaction(n).reactants
            reactants = cat.reaction(n).products
            H0_act    = H0_act_vect[n]-deltaH0_vect[n]
            deltaH0   = -deltaH0_vect[n]
            G0_act    = G0_act_vect[n]-deltaG0_vect[n]
            deltaG0   = -deltaG0_vect[n]
            #G_act     = G_act_vect[n]-deltaG_vect[n]
            #deltaG    = -deltaG_vect[n]

        ts = TS.species_names[n]

        if H0_act < deltaH0:
            H0_act = deltaH0
        if G0_act < deltaG0:
            G0_act = deltaG0
        
        if H0_act < 0.:
            H0_act = 0.
        if G0_act < 0.:
            G0_act = 0.

        H0_tot_TS = H0_tot+H0_act
        G0_tot_TS = G0_tot+G0_act
        #G_tot_TS = G_tot+G_act

        H0_tot += deltaH0
        G0_tot += deltaG0
        #G_tot += deltaG

        for spec in reactants:
            species_dict[spec] -= reactants[spec]

        string = '  '+ts+': '+cat.reaction_equation(n)
        string = string.ljust(58)
        string += '{0:+7.3f} eV {1:+7.3f} eV '.format(H0_tot_TS, G0_tot_TS)
        string += '{0:+7.3f} eV {1:+7.3f} eV'.format(H0_act, G0_act)
        
        if print_data is True:
            print(string)

        for spec in products:
            try: species_dict[spec] += products[spec]
            except: species_dict[spec] = products[spec]

        string = ''
        for spec in [spec for spec in species_dict if species_dict[spec] > 0]:
            string += '{0:2.0f} {1:8s} '.format(species_dict[spec], spec)
        string = string.ljust(58)
        string += '{0:+7.3f} eV {1:+7.3f} eV'.format(H0_tot, G0_tot)

        if print_data is True:
            print(string)

        H0_sequence += [H0_tot_TS]
        H0_sequence += [H0_tot, H0_tot]
        G0_sequence += [G0_tot_TS]
        G0_sequence += [G0_tot, G0_tot]
        #G_sequence += [G_tot_TS]
        #G_sequence += [G_tot, G_tot]
        spec_sequence += [cp.deepcopy(species_dict)]

    react_path_data = {}
    
    react_path_data['spec sequence'] = spec_sequence
    react_path_data['H0 sequence']   = H0_sequence
    react_path_data['G0 sequence']   = G0_sequence
    #react_path_data['G sequence']    = G_sequence

    return react_path_data

################################################################################
# PLOT ENERGY PATH
################################################################################

def plot_energy_paths(site_names, E_sequences, colors, spec_sequences = None,
                      y_min = -3., y_max = 3., show_plot = True, fig_num = 1):

    import matplotlib.pyplot as plt
    from scipy.interpolate import make_interp_spline

    plt.figure(fig_num)

    for i in range(len(spec_sequences[site_names[0]])):

        spec_sequence = spec_sequences[site_names[0]]

        x_text = i*3+0.18

        j = 0

        species_dict = spec_sequence[i]

        for spec in [spec for spec in species_dict 
                     if species_dict[spec] > 0]:

            y_text = y_min+0.18*j+0.1

            if species_dict[spec] > 1:
                string = '{:.0f} '.format(species_dict[spec])
            else:
                string = ''
            string += spec.replace('('+site_names[0]+')', '*')

            plt.text(x_text, y_text, string, fontsize = 12)
            
            j += 1

    y_max_pts = cp.deepcopy(E_sequences[site_names[0]])
    y_min_pts = cp.deepcopy(E_sequences[site_names[0]])

    for site in site_names:
    
        y_plot = E_sequences[site]
        x_vect = range(len(y_plot))
    
        for i in x_vect:
            if y_plot[i] > y_max_pts[i]:
                y_max_pts[i] = y_plot[i]
            if y_plot[i] < y_min_pts[i]:
                y_min_pts[i] = y_plot[i]

    j = 0

    for site in site_names:

        color = colors[site]
    
        y_plot = E_sequences[site]
    
        #plt.figure(1)
    
        x_vect = range(len(y_plot))
    
        j += 1
    
        x_new = []
        y_new = []
    
        delta = 0.01
    
        for i in x_vect:
    
            if (i+1) % 3 == 0:
    
                x_spline = np.linspace(i-1, i+delta, 10)
    
                x = [i-1, i-delta, i+delta]
                y = [y_plot[i-1], y_plot[i], y_plot[i]]
    
                spl = make_interp_spline(x, y, k = 2)
                y_spline = spl(x_spline)
    
                x_new += list(x_spline)
                y_new += list(y_spline)
    
                x_spline = np.linspace(i-delta, i+1, 10)
    
                x = [i-delta, i+delta, i+1]
                y = [y_plot[i], y_plot[i], y_plot[i+1]]
    
                spl = make_interp_spline(x, y, k = 2)
                y_spline = spl(x_spline)
    
                x_new += list(x_spline)
                y_new += list(y_spline)
    
                text = '{:.2f}'.format(y_plot[i])
                x_text = i-0.26
                #y_text = y_min_pts[i]-j*0.1-0.1
                y_ave = (4.*y_min_pts[i]+y_min_pts[i-1]+y_min_pts[i+1])/6.
                y_text = y_ave-j*0.1-0.1
                #y_text = y_min+1.0+j*0.1
    
                plt.text(x_text, y_text, text, color = color, fontsize = 10)
    
                text = '{:.2f}'.format(y_plot[i]-y_plot[i-1])
                x_text = i-0.26
                y_text = y_max_pts[i]-j*0.1+0.6
    
                plt.text(x_text, y_text, text, color = color, fontsize = 10)
    
            else:
                x_new += [i]
                y_new += [y_plot[i]]

                if i % 3 == 0:
    
                    text = '{:.2f}'.format(y_plot[i])
                    x_text = i+0.24
                    y_text = y_min_pts[i]-j*0.1-0.1
                    #y_text = y_min+1.0+j*0.1
        
                    plt.text(x_text, y_text, text, color = color, fontsize = 10)
    
        plt.plot(x_new, y_new, color = color)
        
    plt.axis([0, len(x_vect)-1, y_min, y_max])

    if show_plot is True:
        plt.show()

################################################################################
# END
################################################################################
