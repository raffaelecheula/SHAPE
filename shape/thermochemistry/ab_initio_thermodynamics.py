################################################################################
# Raffaele Cheula, cheula.raffaele@gmail.com
################################################################################

import numpy as np
import ase
from scipy.optimize import least_squares
from scipy.interpolate import RectBivariateSpline
from shape.units import eV, molecule, Rgas
from shape.thermochemistry.wulff_with_support import WulffShape

################################################################################
# STANDARD CONCENTRATION ENTROPY
################################################################################

def standard_concentration_entropy(temperature, sitedensity):

    area = 1./sitedensity/100.0**2

    N_over_A = np.exp(1./3)*(10.0**5/(ase.units._k*temperature))**(2./3)
    S_c = (1 - np.log(N_over_A) - np.log(area))*ase.units.kB

    return S_c

################################################################################
# ENTHALPY AND ENTROPY FUNCTIONS
################################################################################

def enthalpy_and_entropy_over_Rgas(thermo, n_points, temper_range,
                                   coeff_clean = None):
    
    T  = np.zeros(n_points)
    H0 = np.zeros(n_points)
    S0 = np.zeros(n_points)
    
    H0_ref = np.zeros(n_points)
    S0_ref = np.zeros(n_points)
    
    for t in range(n_points):
        T[t] = temper_range[0]+t*(temper_range[1]-temper_range[0])/n_points
        H0[t] = thermo.get_internal_energy(temperature = T[t], verbose = False)
        H0[t] *= eV/molecule/Rgas
        S0[t] = thermo.get_entropy(temperature = T[t], verbose = False)
        S0[t] *= eV/molecule/Rgas
        if coeff_clean is not None:
            args = [T[t]] + coeff_clean
            H0[t] -= H0_fitting(*args)
            S0[t] -= S0_fitting(*args)

    return T, H0, S0

################################################################################
# H0 FITTING
################################################################################

def H0_fitting(t, a1, a2, a3, a4, a5, a6, a7):

    H0_fit = a1*t + a2/2*t**2 + a3/3*t**3 + a4/4*t**4 + a5/5*t**5 + a6
    
    return H0_fit

################################################################################
# S0 FITTING
################################################################################

def S0_fitting(t, a1, a2, a3, a4, a5, a6, a7):

    S0_fit = a1*np.log(t) + a2*t + a3/2*t**2 + a4/3*t**3 + a5/4*t**4 + a7

    return S0_fit

################################################################################
# NASA COEFF FITTING
################################################################################

def NASA_coefficients_fitting(x, T, H0, S0):

    a1, a2, a3, a4, a5, a6, a7 = x

    H0_fit = H0_fitting(T, a1, a2, a3, a4, a5, a6, a7)
    S0_fit = S0_fitting(T, a1, a2, a3, a4, a5, a6, a7)

    err = H0-H0_fit
    err = np.append(err, S0-S0_fit)
    
    return err

################################################################################
# NASA COEFFICIENTS
################################################################################

def NASA_coefficients(T, H0, S0, coeff0 = None, print_coeff = False,
                      plot_fit = False):

    coeff = least_squares(NASA_coefficients_fitting, np.ones(7),
                          args = (T, H0, S0)).x

    if coeff0 is not None:

        coeff -= coeff0

    if print_coeff is True:

        print('[{:15.8E}, {:15.8E}, {:15.8E},'.format(*coeff[0:3]))
        print(' {:15.8E}, {:15.8E}, {:15.8E},'.format(*coeff[3:6]))
        print(' {:15.8E}]'.format(coeff[6]))
        print('')

    if plot_fit is True:

        import matplotlib.pyplot as plt

        plt.figure(1)

        plt.subplot(211)
        plt.plot(T, H0, label = 'H0')
        plt.plot(T, H0_fitting(T, *coeff), label = 'H0 fit')
        plt.xlabel('temperature [K]')
        plt.ylabel('enthalpy')
        plt.legend()

        plt.subplot(212)
        plt.plot(T, S0, label = 'S0')
        plt.plot(T, S0_fitting(T, *coeff), label = 'S0 fit')
        plt.xlabel('temperature [K]')
        plt.ylabel('entropy')
        plt.legend()

        plt.show()

    return coeff

################################################################################
# PRINT NASA BLOCK
################################################################################

def NASA_block(name, atom_dict, temper_list, coeff_list, print_block = True):

    block = """species(name   = "{0}",
        atoms  = {1},""".format(name, atom_dict)

    for i in range(len(temper_list)):
        if i == 0:
            block += """
        thermo = (NASA([{0:.2f}, {1:.2f}],""".format(*temper_list[i])
        else:
            block += """,
                  NASA([{0:.2f}, {1:.2f}],""".format(*temper_list[i])
            
        block += """
                       [{0:15.8E}, {1:15.8E}, {2:15.8E},
                        {3:15.8E}, {4:15.8E}, {5:15.8E},
                        {6:15.8E}])""".format(*coeff_list[i])
    block += '))'

    if print_block is True:
        print(block+'\n')

    return block

################################################################################
# FACET SURFACE ENERGY
################################################################################

def surface_energy_fun(planes, facet, rangemu_A, rangemu_B, n_points = 101):

    x_grid = np.arange(n_points)
    y_grid = np.arange(n_points)
    x_vect = rangemu_A[0] + (rangemu_A[1]-rangemu_A[0])/(n_points-1)*x_grid
    y_vect = rangemu_B[0] + (rangemu_B[1]-rangemu_B[0])/(n_points-1)*y_grid
    
    z_vect = np.ones((n_points, n_points))*100
    
    for pln in planes[facet]:
        for x, y in [ (x, y) for x in x_grid for y in y_grid ]:
            z = pln[0]+(pln[1]-pln[0])*x/n_points+(pln[2]-pln[0])*y/n_points
            if z < z_vect[x, y]:
                z_vect[x, y] = z
    
    fun = RectBivariateSpline(x_vect, y_vect, z_vect)

    return fun

################################################################################
# SURFACE ENERGIES FUNS
################################################################################

def surface_energies_funs(planes, miller_list, rangemu_A, rangemu_B,
                          n_points = 101):

    e_surf_funs = []

    for i in range(len(miller_list)):
        e_surf_funs.append(surface_energy_fun(planes    = planes        ,
                                              facet     = miller_list[i],
                                              rangemu_A = rangemu_A     ,
                                              rangemu_B = rangemu_B     ,
                                              n_points  = n_points      ))

    return e_surf_funs

################################################################################
# SURFACE ENERGIES
################################################################################

def surface_energies(e_surf_funs, deltamu_A, deltamu_B):

    e_surf_list = np.zeros(len(e_surf_funs))
    
    for i in range(len(e_surf_funs)):
        e_surf_list[i] = float(e_surf_funs[i](deltamu_A, deltamu_B)[0])

    return e_surf_list

################################################################################
# SCALE FACTOR FROM DIAMETER
################################################################################

def scale_factor_from_diameter(lattice, miller_list, e_surf_list, miller_supp,
                               e_surf_supp, wulff_diameter):

    if e_surf_supp is None or e_surf_supp > 0.:
        e_surf_supp_zero = 1e-10
    else:
        e_surf_supp_zero = e_surf_supp
    
    wulffcut = WulffShape(lattice     = lattice         ,
                          miller_list = miller_list     ,
                          e_surf_list = e_surf_list     ,
                          miller_supp = miller_supp     ,
                          e_surf_supp = e_surf_supp_zero)
    
    wulffcut_diameter = (4./np.pi*wulffcut.miller_area_dict['(support)'])**0.5

    scale_factor = wulff_diameter/wulffcut_diameter

    return scale_factor

################################################################################
# WULFF SHAPE FROM DIAMETER
################################################################################

def wulff_shape_from_diameter(lattice, miller_list, e_surf_list, miller_supp,
                              e_surf_supp, wulff_diameter):

    scale_factor = scale_factor_from_diameter(lattice        = lattice       ,
                                              miller_list    = miller_list   ,
                                              e_surf_list    = e_surf_list   ,
                                              miller_supp    = miller_supp   ,
                                              e_surf_supp    = e_surf_supp   ,
                                              wulff_diameter = wulff_diameter)
    
    e_surf_list = e_surf_list*scale_factor

    if e_surf_supp is not None:
        e_surf_supp = e_surf_supp*scale_factor
    
    wulff = WulffShape(lattice     = lattice    ,
                       miller_list = miller_list,
                       e_surf_list = e_surf_list,
                       miller_supp = miller_supp,
                       e_surf_supp = e_surf_supp)

    return wulff

################################################################################
# SCALE FACTOR FROM VOLUME
################################################################################

def scale_factor_from_volume(lattice, miller_list, e_surf_list,
                             miller_supp, e_surf_supp, wulff_volume):

    wulff = WulffShape(lattice     = lattice    ,
                       miller_list = miller_list,
                       e_surf_list = e_surf_list,
                       miller_supp = miller_supp,
                       e_surf_supp = e_surf_supp)
    
    scale_factor = (wulff_volume/wulff.volume)**(1./3.)

    return scale_factor

################################################################################
# WULFF SHAPE AT GIVEN VOLUME
################################################################################

def wulff_shape_from_volume(lattice, miller_list, e_surf_list,
                            miller_supp, e_surf_supp, wulff_volume):

    scale_factor = scale_factor_from_volume(lattice      = lattice     ,
                                            miller_list  = miller_list ,
                                            e_surf_list  = e_surf_list ,
                                            miller_supp  = miller_supp ,
                                            e_surf_supp  = e_surf_supp ,
                                            wulff_volume = wulff_volume)
    
    e_surf_list = e_surf_list * scale_factor
    e_surf_supp = e_surf_supp * scale_factor
    
    wulff = WulffShape(lattice     = lattice    ,
                       miller_list = miller_list,
                       e_surf_list = e_surf_list,
                       miller_supp = miller_supp,
                       e_surf_supp = e_surf_supp)

    return wulff

################################################################################
# WULFF RELATIVE AREAS
################################################################################

def wulff_relative_areas(lattice, miller_list, miller_supp, e_surf_funs,
                         e_surf_supp, deltamu_A, deltamu_B):

    e_surf_list = surface_energies(e_surf_funs = e_surf_funs,
                                   deltamu_A   = deltamu_A,
                                   deltamu_B   = deltamu_B)

    wulff = WulffShape(lattice     = lattice    ,
                       miller_list = miller_list,
                       e_surf_list = e_surf_list,
                       miller_supp = miller_supp,
                       e_surf_supp = e_surf_supp)

    relative_areas = []
    for facet in miller_list:
        relative_areas.append(wulff.miller_area_dict[facet]/wulff.surface_area)

    return relative_areas

################################################################################
# WULFF SPECIFIC AREAS
################################################################################

def wulff_specific_areas(lattice, miller_list, miller_supp, e_surf_funs,
                         e_surf_supp, deltamu_A, deltamu_B, wulff_volume):

    e_surf_list = surface_energies(e_surf_funs = e_surf_funs,
                                   deltamu_A   = deltamu_A  ,
                                   deltamu_B   = deltamu_B  )

    wulff = WulffShape(lattice     = lattice    ,
                       miller_list = miller_list,
                       e_surf_list = e_surf_list,
                       miller_supp = miller_supp,
                       e_surf_supp = e_surf_supp)
    
    scale_factor = (wulff_volume/wulff.volume)**(1./3.)

    specific_areas = []
    for facet in miller_list:
        spec_area = wulff.miller_area_dict[facet]*scale_factor**2/wulff_volume
        specific_areas.append(spec_area)

    return specific_areas

################################################################################
# PLOT DOUBLE CHEMICAL POTENTIAL
################################################################################

def plot_double_chemical_potential(x_vect, y_vect, z_vect):

    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D

    fig = plt.figure()
    ax = fig.gca(projection = '3d')
    
    x_plot, y_plot = np.meshgrid(x_vect, y_vect)
    
    ax.plot_surface(x_plot, y_plot, z_vect)
    
    plt.show()

################################################################################
# END
################################################################################
