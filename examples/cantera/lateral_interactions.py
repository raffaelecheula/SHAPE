################################################################################
# Raffaele Cheula, LCCP, Politecnico di Milano, cheula.raffaele@gmail.com
################################################################################

from shape.units import *

################################################################################
# LATERAL INTERACTIONS CLASS
################################################################################

class LateralInteraction():

    def __init__(self, ads_name, cov_dep, coeffs, units):

        self.ads_name = ads_name
        self.cov_dep  = cov_dep
        self.coeffs   = coeffs
        self.units    = units

        self.m = coeffs[0]

        if coeffs[0] != 0.:
            self.x0 = -coeffs[1]/coeffs[0]
        else:
            self.x0 = 0. 

    def get_coeff_correction(self, coverage):

        if coverage > self.x0:
            corr = self.m/Rgas*(coverage-self.x0) * self.units
        else:
            corr = 0.

        return corr

################################################################################
# LATERAL INTERACTIONS
################################################################################

lat_dict = {}

"""
Lateral interactions for the adsorbates on Rh(100)
"""

lat_dict['100'] = []

lat_dict['100'] += [LateralInteraction(ads_name = 'CO2(100)',
                                       cov_dep  = 'CO(100)',
                                       coeffs   = [1.96, -0.87],
                                       units    = eV/molecule)]

#lat_dict['100'] += [LateralInteraction(ads_name = 'H(100)',
#                                       cov_dep  = 'CO(100)',
#                                       coeffs   = [0.62, -0.28],
#                                       units    = eV/molecule)]

lat_dict['100'] += [LateralInteraction(ads_name = 'H(100)',
                                       cov_dep  = 'CO(100)',
                                       coeffs   = [1.72, -0.86],
                                       units    = eV/molecule)]

lat_dict['100'] += [LateralInteraction(ads_name = 'CO(100)',
                                       cov_dep  = 'CO(100)',
                                       coeffs   = [0.83, -0.36],
                                       units    = eV/molecule)]

lat_dict['100'] += [LateralInteraction(ads_name = 'H2O(100)',
                                       cov_dep  = 'CO(100)',
                                       coeffs   = [0.47, -0.25],
                                       units    = eV/molecule)]

lat_dict['100'] += [LateralInteraction(ads_name = 'OH(100)',
                                       cov_dep  = 'CO(100)',
                                       coeffs   = [0.70, -0.31],
                                       units    = eV/molecule)]

lat_dict['100'] += [LateralInteraction(ads_name = 'O(100)',
                                       cov_dep  = 'CO(100)',
                                       coeffs   = [1.46, -0.48],
                                       units    = eV/molecule)]

lat_dict['100'] += [LateralInteraction(ads_name = 'cCOOH(100)',
                                       cov_dep  = 'CO(100)',
                                       coeffs   = [2.02, -0.85],
                                       units    = eV/molecule)]

lat_dict['100'] += [LateralInteraction(ads_name = 'tCOOH(100)',
                                       cov_dep  = 'CO(100)',
                                       coeffs   = [2.02, -0.85],
                                       units    = eV/molecule)]

lat_dict['100'] += [LateralInteraction(ads_name = 'HCOO(100)',
                                       cov_dep  = 'CO(100)',
                                       coeffs   = [1.14, -0.57],
                                       units    = eV/molecule)]

lat_dict['100'] += [LateralInteraction(ads_name = 'C(100)',
                                       cov_dep  = 'CO(100)',
                                       coeffs   = [2.42, -1.21],
                                       units    = eV/molecule)]

lat_dict['100'] += [LateralInteraction(ads_name = 'CH(100)',
                                       cov_dep  = 'CO(100)',
                                       coeffs   = [1.02, -0.45],
                                       units    = eV/molecule)]

lat_dict['100'] += [LateralInteraction(ads_name = 'CH2(100)',
                                       cov_dep  = 'CO(100)',
                                       coeffs   = [0.62, -0.28],
                                       units    = eV/molecule)]

lat_dict['100'] += [LateralInteraction(ads_name = 'CH3(100)',
                                       cov_dep  = 'CO(100)',
                                       coeffs   = [1.66, -0.73],
                                       units    = eV/molecule)]

lat_dict['100'] += [LateralInteraction(ads_name = 'OCH(100)',
                                       cov_dep  = 'CO(100)',
                                       coeffs   = [1.09, -0.54],
                                       units    = eV/molecule)]

"""
Lateral interactions for the adsorbated on Rh(110)
"""

lat_dict['110'] = []

lat_dict['110'] += [LateralInteraction(ads_name = 'CO2(110)',
                                       cov_dep  = 'CO(110)',
                                       coeffs   = [1.82, -0.87],
                                       units    = eV/molecule)]

lat_dict['110'] += [LateralInteraction(ads_name = 'H(110)',
                                       cov_dep  = 'CO(110)',
                                       coeffs   = [0.07, -0.03],
                                       units    = eV/molecule)]

lat_dict['110'] += [LateralInteraction(ads_name = 'CO(110)',
                                       cov_dep  = 'CO(110)',
                                       coeffs   = [0.56, -0.25],
                                       units    = eV/molecule)]

lat_dict['110'] += [LateralInteraction(ads_name = 'H2O(110)',
                                       cov_dep  = 'CO(110)',
                                       coeffs   = [0.44, -0.22],
                                       units    = eV/molecule)]

lat_dict['110'] += [LateralInteraction(ads_name = 'OH(110)',
                                       cov_dep  = 'CO(110)',
                                       coeffs   = [0.62, -0.31],
                                       units    = eV/molecule)]

lat_dict['110'] += [LateralInteraction(ads_name = 'O(110)',
                                       cov_dep  = 'CO(110)',
                                       coeffs   = [1.03, -0.48],
                                       units    = eV/molecule)]

lat_dict['110'] += [LateralInteraction(ads_name = 'cCOOH(110)',
                                       cov_dep  = 'CO(110)',
                                       coeffs   = [0.67, -0.33],
                                       units    = eV/molecule)]

lat_dict['110'] += [LateralInteraction(ads_name = 'tCOOH(110)',
                                       cov_dep  = 'CO(110)',
                                       coeffs   = [0.67, -0.33],
                                       units    = eV/molecule)]

lat_dict['110'] += [LateralInteraction(ads_name = 'HCOO(110)',
                                       cov_dep  = 'CO(110)',
                                       coeffs   = [0.82, -0.41],
                                       units    = eV/molecule)]

lat_dict['110'] += [LateralInteraction(ads_name = 'C(110)',
                                       cov_dep  = 'CO(110)',
                                       coeffs   = [0.24, -0.10],
                                       units    = eV/molecule)]

lat_dict['110'] += [LateralInteraction(ads_name = 'CH(110)',
                                       cov_dep  = 'CO(110)',
                                       coeffs   = [0.18, -0.08],
                                       units    = eV/molecule)]

lat_dict['110'] += [LateralInteraction(ads_name = 'CH2(110)',
                                       cov_dep  = 'CO(110)',
                                       coeffs   = [0.55, -0.28],
                                       units    = eV/molecule)]

lat_dict['110'] += [LateralInteraction(ads_name = 'CH3(110)',
                                       cov_dep  = 'CO(110)',
                                       coeffs   = [1.46, -0.73],
                                       units    = eV/molecule)]

lat_dict['110'] += [LateralInteraction(ads_name = 'OCH(110)',
                                       cov_dep  = 'CO(110)',
                                       coeffs   = [0.40, -0.20],
                                       units    = eV/molecule)]

"""
Lateral interactions for the adsorbated on Rh(111)
"""

lat_dict['111'] = []

lat_dict['111'] += [LateralInteraction(ads_name = 'CO2(111)',
                                       cov_dep  = 'CO(111)',
                                       coeffs   = [2.38, -0.87],
                                       units    = eV/molecule)]

lat_dict['111'] += [LateralInteraction(ads_name = 'H(111)',
                                       cov_dep  = 'CO(111)',
                                       coeffs   = [0.30, -0.13],
                                       units    = eV/molecule)]

lat_dict['111'] += [LateralInteraction(ads_name = 'CO(111)',
                                       cov_dep  = 'CO(111)',
                                       coeffs   = [1.26, -0.36],
                                       units    = eV/molecule)]

lat_dict['111'] += [LateralInteraction(ads_name = 'H2O(111)',
                                       cov_dep  = 'CO(111)',
                                       coeffs   = [0.65, -0.25],
                                       units    = eV/molecule)]

lat_dict['111'] += [LateralInteraction(ads_name = 'OH(111)',
                                       cov_dep  = 'CO(111)',
                                       coeffs   = [0.81, -0.31],
                                       units    = eV/molecule)]

lat_dict['111'] += [LateralInteraction(ads_name = 'O(111)',
                                       cov_dep  = 'CO(111)',
                                       coeffs   = [1.34, -0.48],
                                       units    = eV/molecule)]

lat_dict['111'] += [LateralInteraction(ads_name = 'cCOOH(111)',
                                       cov_dep  = 'CO(111)',
                                       coeffs   = [2.34, -0.85],
                                       units    = eV/molecule)]

lat_dict['111'] += [LateralInteraction(ads_name = 'tCOOH(111)',
                                       cov_dep  = 'CO(111)',
                                       coeffs   = [2.34, -0.85],
                                       units    = eV/molecule)]

lat_dict['111'] += [LateralInteraction(ads_name = 'HCOO(111)',
                                       cov_dep  = 'CO(111)',
                                       coeffs   = [2.28, -0.84],
                                       units    = eV/molecule)]

lat_dict['111'] += [LateralInteraction(ads_name = 'C(111)',
                                       cov_dep  = 'CO(111)',
                                       coeffs   = [0.96, -0.38],
                                       units    = eV/molecule)]

lat_dict['111'] += [LateralInteraction(ads_name = 'CH(111)',
                                       cov_dep  = 'CO(111)',
                                       coeffs   = [1.17, -0.45],
                                       units    = eV/molecule)]

lat_dict['111'] += [LateralInteraction(ads_name = 'CH2(111)',
                                       cov_dep  = 'CO(111)',
                                       coeffs   = [0.72, -0.28],
                                       units    = eV/molecule)]

lat_dict['111'] += [LateralInteraction(ads_name = 'CH3(111)',
                                       cov_dep  = 'CO(111)',
                                       coeffs   = [1.20, -0.27],
                                       units    = eV/molecule)]

lat_dict['111'] += [LateralInteraction(ads_name = 'OCH(111)',
                                       cov_dep  = 'CO(111)',
                                       coeffs   = [2.83, -1.10],
                                       units    = eV/molecule)]

"""
Lateral interactions for the adsorbated on Rh(211)
"""

lat_dict['211'] = []

lat_dict['211'] += [LateralInteraction(ads_name = 'CO2(211)',
                                       cov_dep  = 'CO(211)',
                                       coeffs   = [1.13, -0.57],
                                       units    = eV/molecule)]

lat_dict['211'] += [LateralInteraction(ads_name = 'H(211)',
                                       cov_dep  = 'CO(211)',
                                       coeffs   = [0.01, -0.00],
                                       units    = eV/molecule)]

lat_dict['211'] += [LateralInteraction(ads_name = 'CO(211)',
                                       cov_dep  = 'CO(211)',
                                       coeffs   = [1.07, -0.51],
                                       units    = eV/molecule)]

lat_dict['211'] += [LateralInteraction(ads_name = 'H2O(211)',
                                       cov_dep  = 'CO(211)',
                                       coeffs   = [0.16, -0.08],
                                       units    = eV/molecule)]

lat_dict['211'] += [LateralInteraction(ads_name = 'OH(211)',
                                       cov_dep  = 'CO(211)',
                                       coeffs   = [0.33, -0.16],
                                       units    = eV/molecule)]

lat_dict['211'] += [LateralInteraction(ads_name = 'O(211)',
                                       cov_dep  = 'CO(211)',
                                       coeffs   = [0.87, -0.44],
                                       units    = eV/molecule)]

lat_dict['211'] += [LateralInteraction(ads_name = 'cCOOH(211)',
                                       cov_dep  = 'CO(211)',
                                       coeffs   = [0.44, -0.22],
                                       units    = eV/molecule)]

lat_dict['211'] += [LateralInteraction(ads_name = 'tCOOH(211)',
                                       cov_dep  = 'CO(211)',
                                       coeffs   = [0.44, -0.22],
                                       units    = eV/molecule)]

lat_dict['211'] += [LateralInteraction(ads_name = 'HCOO(211)',
                                       cov_dep  = 'CO(211)',
                                       coeffs   = [0.87, -0.43],
                                       units    = eV/molecule)]

lat_dict['211'] += [LateralInteraction(ads_name = 'C(211)',
                                       cov_dep  = 'CO(211)',
                                       coeffs   = [0.10, -0.05],
                                       units    = eV/molecule)]

lat_dict['211'] += [LateralInteraction(ads_name = 'CH(211)',
                                       cov_dep  = 'CO(211)',
                                       coeffs   = [0.25, -0.13],
                                       units    = eV/molecule)]

lat_dict['211'] += [LateralInteraction(ads_name = 'CH2(211)',
                                       cov_dep  = 'CO(211)',
                                       coeffs   = [0.05, -0.02],
                                       units    = eV/molecule)]

lat_dict['211'] += [LateralInteraction(ads_name = 'CH3(211)',
                                       cov_dep  = 'CO(211)',
                                       coeffs   = [0.62, -0.31],
                                       units    = eV/molecule)]

lat_dict['211'] += [LateralInteraction(ads_name = 'OCH(211)',
                                       cov_dep  = 'CO(211)',
                                       coeffs   = [0.05, -0.03],
                                       units    = eV/molecule)]

"""
Lateral interactions for the adsorbated on Rh(311)
"""

lat_dict['311'] = []

lat_dict['311'] += [LateralInteraction(ads_name = 'CO2(311)',
                                       cov_dep  = 'CO(311)',
                                       coeffs   = [1.13, -0.57],
                                       units    = eV/molecule)]

lat_dict['311'] += [LateralInteraction(ads_name = 'H(311)',
                                       cov_dep  = 'CO(311)',
                                       coeffs   = [0.01, -0.00],
                                       units    = eV/molecule)]

lat_dict['311'] += [LateralInteraction(ads_name = 'CO(311)',
                                       cov_dep  = 'CO(311)',
                                       coeffs   = [1.07, -0.51],
                                       units    = eV/molecule)]

lat_dict['311'] += [LateralInteraction(ads_name = 'H2O(311)',
                                       cov_dep  = 'CO(311)',
                                       coeffs   = [0.16, -0.08],
                                       units    = eV/molecule)]

lat_dict['311'] += [LateralInteraction(ads_name = 'OH(311)',
                                       cov_dep  = 'CO(311)',
                                       coeffs   = [0.33, -0.16],
                                       units    = eV/molecule)]

lat_dict['311'] += [LateralInteraction(ads_name = 'O(311)',
                                       cov_dep  = 'CO(311)',
                                       coeffs   = [0.67, -0.44],
                                       units    = eV/molecule)]

lat_dict['311'] += [LateralInteraction(ads_name = 'cCOOH(311)',
                                       cov_dep  = 'CO(311)',
                                       coeffs   = [0.24, -0.22],
                                       units    = eV/molecule)]

lat_dict['311'] += [LateralInteraction(ads_name = 'tCOOH(311)',
                                       cov_dep  = 'CO(311)',
                                       coeffs   = [0.24, -0.22],
                                       units    = eV/molecule)]

lat_dict['311'] += [LateralInteraction(ads_name = 'HCOO(311)',
                                       cov_dep  = 'CO(311)',
                                       coeffs   = [0.87, -0.43],
                                       units    = eV/molecule)]

lat_dict['311'] += [LateralInteraction(ads_name = 'C(311)',
                                       cov_dep  = 'CO(311)',
                                       coeffs   = [0.10, -0.05],
                                       units    = eV/molecule)]

lat_dict['311'] += [LateralInteraction(ads_name = 'CH(311)',
                                       cov_dep  = 'CO(311)',
                                       coeffs   = [0.25, -0.13],
                                       units    = eV/molecule)]

lat_dict['311'] += [LateralInteraction(ads_name = 'CH2(311)',
                                       cov_dep  = 'CO(311)',
                                       coeffs   = [0.05, -0.02],
                                       units    = eV/molecule)]

lat_dict['311'] += [LateralInteraction(ads_name = 'CH3(311)',
                                       cov_dep  = 'CO(311)',
                                       coeffs   = [0.62, -0.31],
                                       units    = eV/molecule)]

lat_dict['311'] += [LateralInteraction(ads_name = 'OCH(311)',
                                       cov_dep  = 'CO(311)',
                                       coeffs   = [0.05, -0.03],
                                       units    = eV/molecule)]

"""
Lateral interactions for the adsorbated on Rh(331)
"""

lat_dict['331'] = []

lat_dict['331'] += [LateralInteraction(ads_name = 'CO2(331)',
                                       cov_dep  = 'CO(331)',
                                       coeffs   = [1.13, -0.57],
                                       units    = eV/molecule)]

lat_dict['331'] += [LateralInteraction(ads_name = 'H(331)',
                                       cov_dep  = 'CO(331)',
                                       coeffs   = [0.01, -0.00],
                                       units    = eV/molecule)]

lat_dict['331'] += [LateralInteraction(ads_name = 'CO(331)',
                                       cov_dep  = 'CO(331)',
                                       coeffs   = [1.07, -0.51],
                                       units    = eV/molecule)]

lat_dict['331'] += [LateralInteraction(ads_name = 'H2O(331)',
                                       cov_dep  = 'CO(331)',
                                       coeffs   = [0.16, -0.08],
                                       units    = eV/molecule)]

lat_dict['331'] += [LateralInteraction(ads_name = 'OH(331)',
                                       cov_dep  = 'CO(331)',
                                       coeffs   = [0.33, -0.16],
                                       units    = eV/molecule)]

lat_dict['331'] += [LateralInteraction(ads_name = 'O(331)',
                                       cov_dep  = 'CO(331)',
                                       coeffs   = [0.87, -0.44],
                                       units    = eV/molecule)]

lat_dict['331'] += [LateralInteraction(ads_name = 'cCOOH(331)',
                                       cov_dep  = 'CO(331)',
                                       coeffs   = [0.44, -0.22],
                                       units    = eV/molecule)]

lat_dict['331'] += [LateralInteraction(ads_name = 'tCOOH(331)',
                                       cov_dep  = 'CO(331)',
                                       coeffs   = [0.44, -0.22],
                                       units    = eV/molecule)]

lat_dict['331'] += [LateralInteraction(ads_name = 'HCOO(331)',
                                       cov_dep  = 'CO(331)',
                                       coeffs   = [0.87, -0.43],
                                       units    = eV/molecule)]

lat_dict['331'] += [LateralInteraction(ads_name = 'C(331)',
                                       cov_dep  = 'CO(331)',
                                       coeffs   = [0.10, -0.05],
                                       units    = eV/molecule)]

lat_dict['331'] += [LateralInteraction(ads_name = 'CH(331)',
                                       cov_dep  = 'CO(331)',
                                       coeffs   = [0.25, -0.13],
                                       units    = eV/molecule)]

lat_dict['331'] += [LateralInteraction(ads_name = 'CH2(331)',
                                       cov_dep  = 'CO(331)',
                                       coeffs   = [0.05, -0.02],
                                       units    = eV/molecule)]

lat_dict['331'] += [LateralInteraction(ads_name = 'CH3(331)',
                                       cov_dep  = 'CO(331)',
                                       coeffs   = [0.62, -0.31],
                                       units    = eV/molecule)]

lat_dict['331'] += [LateralInteraction(ads_name = 'OCH(331)',
                                       cov_dep  = 'CO(331)',
                                       coeffs   = [0.05, -0.03],
                                       units    = eV/molecule)]

"""
Lateral interactions for the adsorbated on Rh(int)
"""

lat_dict['int'] = []

lat_dict['int'] += [LateralInteraction(ads_name = 'CO2(int)',
                                       cov_dep  = 'CO(int)',
                                       coeffs   = [2.38, -0.87],
                                       units    = eV/molecule)]

lat_dict['int'] += [LateralInteraction(ads_name = 'H(int)',
                                       cov_dep  = 'CO(int)',
                                       coeffs   = [0.30, -0.13],
                                       units    = eV/molecule)]

lat_dict['int'] += [LateralInteraction(ads_name = 'CO(int)',
                                       cov_dep  = 'CO(int)',
                                       coeffs   = [1.26, -0.36],
                                       units    = eV/molecule)]

lat_dict['int'] += [LateralInteraction(ads_name = 'H2O(int)',
                                       cov_dep  = 'CO(int)',
                                       coeffs   = [0.65, -0.25],
                                       units    = eV/molecule)]

lat_dict['int'] += [LateralInteraction(ads_name = 'OH(int)',
                                       cov_dep  = 'CO(int)',
                                       coeffs   = [0.81, -0.31],
                                       units    = eV/molecule)]

lat_dict['int'] += [LateralInteraction(ads_name = 'O(int)',
                                       cov_dep  = 'CO(int)',
                                       coeffs   = [1.34, -0.48],
                                       units    = eV/molecule)]

lat_dict['int'] += [LateralInteraction(ads_name = 'cCOOH(int)',
                                       cov_dep  = 'CO(int)',
                                       coeffs   = [2.34, -0.85],
                                       units    = eV/molecule)]

lat_dict['int'] += [LateralInteraction(ads_name = 'tCOOH(int)',
                                       cov_dep  = 'CO(int)',
                                       coeffs   = [2.34, -0.85],
                                       units    = eV/molecule)]

lat_dict['int'] += [LateralInteraction(ads_name = 'HCOO(int)',
                                       cov_dep  = 'CO(int)',
                                       coeffs   = [2.28, -0.84],
                                       units    = eV/molecule)]

lat_dict['int'] += [LateralInteraction(ads_name = 'C(int)',
                                       cov_dep  = 'CO(int)',
                                       coeffs   = [0.96, -0.38],
                                       units    = eV/molecule)]

lat_dict['int'] += [LateralInteraction(ads_name = 'CH(int)',
                                       cov_dep  = 'CO(int)',
                                       coeffs   = [1.17, -0.45],
                                       units    = eV/molecule)]

lat_dict['int'] += [LateralInteraction(ads_name = 'CH2(int)',
                                       cov_dep  = 'CO(int)',
                                       coeffs   = [0.72, -0.28],
                                       units    = eV/molecule)]

lat_dict['int'] += [LateralInteraction(ads_name = 'CH3(int)',
                                       cov_dep  = 'CO(int)',
                                       coeffs   = [1.20, -0.27],
                                       units    = eV/molecule)]

lat_dict['int'] += [LateralInteraction(ads_name = 'OCH(int)',
                                       cov_dep  = 'CO(int)',
                                       coeffs   = [2.83, -1.10],
                                       units    = eV/molecule)]

################################################################################
# END
################################################################################
