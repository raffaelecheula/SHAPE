#!/usr/bin/env python3

################################################################################
# Raffaele Cheula, LCCP, Politecnico di Milano, raffaele.cheula@polimi.it
################################################################################

import matplotlib.pyplot as plt
from shape.ase_utils import get_symbols_list
from shape.qe_utils import (
    ReadQeOut,
    read_projwfc,
    get_atoms_details,
    print_atoms_details,
    plot_band_levels,
    scale_band_energies,
    get_pdos,
)

################################################################################
# READ QE OUT
################################################################################

kpoint = 1

# TODO: EXTEND TO ALL K POINTS

qe_out = ReadQeOut(filename = 'pw.out')

atoms = qe_out.get_atoms()

qe_out.read_bands(scale_band_energies = True)

e_fermi = qe_out.e_fermi

bands_energies = qe_out.e_bands_dict[kpoint]

symbols_list = get_symbols_list(atoms)

################################################################################
# PLOT BAND LEVELS
################################################################################

plot_bands = True

atom_num_list = [i+1 for i in range(len(atoms))]

delta_e = 0.

e_min = -10.
e_max =   5.

num_min_print = 0.

color_list = ['limegreen', 'darkorange', 'royalblue', 'crimson']

if plot_bands is True:

    filename = 'projwfc/projwfc.out'

    states_list, bands_list = read_projwfc(filename = filename,
                                           kpoint   = kpoint  )

    bands_list = scale_band_energies(bands_list, e_fermi)

    color_dict = {}
    for i in range(len(symbols_list)):
        color_dict[symbols_list[i]] = color_list[i]

    atoms_pp_list = get_atoms_details(states_list   = states_list  ,
                                      bands_list    = bands_list   ,
                                      atom_num_list = atom_num_list,
                                      delta_e       = delta_e      ,
                                      color_dict    = color_dict   )

    print_atoms_details(atoms_pp_list = atoms_pp_list     ,
                        filename      = 'atom_details.out')
    
    fig = plot_band_levels(atoms_pp_list  = atoms_pp_list ,
                           num_min_print  = num_min_print ,
                           bands_energies = bands_energies,
                           e_min          = e_min         ,
                           e_max          = e_max         )

    plt.plot([0., 1.], [0.]*2, color = 'red')

    #plt.text(0.9, 0.1, 'E Fermi', size = 'xx-small', color = 'red')

    plt.ylabel('energy [eV]')

    plt.savefig('band_levels.png', dpi = 300)

################################################################################
# PLOT PDOS
################################################################################

plot_pdos = True

if plot_pdos is True:

    fig = plt.figure(1)
    
    filename = 'projwfc/pdos.pdos_tot'
    
    energy, dos = get_pdos(filename = filename,
                           e_fermi  = e_fermi )
    
    x_max = max(dos)+1.
    
    plt.plot(dos, energy)
    
    filename_str = 'projwfc/pdos.pdos_atm#{0}({1})_wfc#{2}({3})'
    
    filename_list = []
    
    alpha_dict = {}
    alpha_dict['s'] = 0.50
    alpha_dict['p'] = 0.75
    alpha_dict['d'] = 1.00
    
    for atom in atoms_pp_list:

        for state in atom.states:
    
            orbital = state.orbital[:1]
    
            filename = filename_str.format(atom.atom_num    ,
                                           atom.element     ,
                                           state.shell_num  ,
                                           orbital          )
    
            if filename not in filename_list:
    
                filename_list += [filename]
    
                energy, dos = get_pdos(filename = filename,
                                       e_fermi  = e_fermi )
    
                alpha = alpha_dict[orbital]
    
                plt.plot(dos, energy, color = atom.color, alpha = alpha)
    
    plt.plot([0., x_max], [0.]*2, color = 'red')
    
    #plt.text(0.9, 0.1, 'E Fermi', size = 'xx-small', color = 'red')
    
    plt.xlim([0., x_max])
    plt.ylim([e_min, e_max])
    
    plt.xlabel('dos [states/eV]')
    plt.ylabel('energy [eV]')
    
    plt.savefig('pdos.png', dpi = 300)

################################################################################
# END
################################################################################
