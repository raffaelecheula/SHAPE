#!/usr/bin/env python3

################################################################################
# Raffaele Cheula, LCCP, Politecnico di Milano, raffaele.cheula@polimi.it
################################################################################

from __future__ import absolute_import, division, print_function
import os, argparse
import numpy as np
from distutils.util import strtobool
from collections import OrderedDict
from qe_utils import write_projwfc_input

################################################################################
# PARSE CMD LINE ARGUMENTS
################################################################################

parser = argparse.ArgumentParser()

fbool = lambda x: bool(strtobool(x))

parser.add_argument('--change_dir', '-c',
                    type     = fbool,
                    nargs    = 1    ,
                    required = False,
                    default  = [True])

parser.add_argument('--write_input', '-w',
                    type     = fbool,
                    nargs    = 1    ,
                    required = False,
                    default  = [True])

parser.add_argument('--run_qe_bin', '-r',
                    type     = fbool,
                    nargs    = 1    ,
                    required = False,
                    default  = [True])

parser.add_argument('--postprocess', '-p',
                    type     = fbool,
                    nargs    = 1    ,
                    required = False,
                    default  = [True])

parsed_args = parser.parse_args()

################################################################################
# CHANGE DIRECTORY
################################################################################

change_dir = parsed_args.change_dir[0]

if change_dir is True:

    try: os.mkdir('projwfc')
    except: pass

    os.chdir('projwfc')

    pw_out_dir = '../'

else:

    pw_out_dir = './'

################################################################################
# PRINT PROJWFC INP
################################################################################

prefix = 'calc'
outdir = 'tmp'

write_input = parsed_args.write_input[0]

if write_input is True:

    proj_data = OrderedDict()

    proj_data['prefix']  = prefix
    proj_data['outdir']  = pw_out_dir+outdir
    proj_data['ngauss']  = -99
    proj_data['Emin']    = -30.
    proj_data['Emax']    = +20.
    proj_data['DeltaE']  = 0.01
    proj_data['lsym']    = True
    proj_data['filpdos'] = 'pdos'

    write_projwfc_input(proj_data = proj_data, filename = 'projwfc.inp')

################################################################################
# RUN PROJWFC
################################################################################

run_qe_bin = parsed_args.run_qe_bin[0]

if run_qe_bin is True:
    
    os.system('projwfc.x < projwfc.inp > projwfc.out')

################################################################################
# POSTPROCESS
################################################################################

postprocess = parsed_args.postprocess[0]

if postprocess is True:

    ############################################################################
    # READ QE OUT
    ############################################################################

    kpoint = 1
    
    # TODO: EXTEND TO ALL K POINTS
    
    qe_out = ReadQeOut(filename = 'pw.out')
    
    atoms = qe_out.get_atoms()
    
    qe_out.read_bands(scale_band_energies = True)
    
    e_fermi = qe_out.e_fermi
    
    bands_energies = qe_out.e_bands_dict[kpoint]
    
    symbols_list = get_symbols_list(atoms)

    atom_num_list = [i+1 for i in range(len(atoms))]

    ############################################################################
    # PLOT BAND LEVELS
    ############################################################################

    delta_e = 0.
    
    e_min = -10.
    e_max =   5.
    
    num_min_print = 0.
    
    color_list = ['limegreen', 'darkorange', 'royalblue', 'crimson']

    filename = 'projwfc.out'

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

    ############################################################################
    # PLOT PDOS
    ############################################################################

    fig = plt.figure(1)
    
    filename = 'pdos.pdos_tot'
    
    energy, dos = get_pdos(filename = filename,
                           e_fermi  = e_fermi )
    
    x_max = max(dos)+1.
    
    plt.plot(dos, energy)
    
    filename_str = 'pdos.pdos_atm#{0}({1})_wfc#{2}({3})'
    
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
# CHANGE DIRECTORY
################################################################################

if change_dir is True:

    os.chdir('..')

################################################################################
# END
################################################################################
