################################################################################
# Raffaele Cheula, cheula.raffaele@gmail.com
################################################################################

import ast
import re
import os
import pickle
import numpy as np
import copy as cp
from ase import Atoms
from ase.io.espresso import SSSP_VALENCE
from ase.data import atomic_numbers
from ase.units import create_units
from ase.constraints import FixAtoms, FixCartesian
from ase.calculators.espresso import Espresso
from shape.ase_utils import get_symbols_list, get_symbols_dict

################################################################################
# READ QUANTUM ESPRESSO OUT
################################################################################

def read_qe_out(filename):

    units = create_units('2006')

    atoms = Atoms(pbc = True)
    cell = np.zeros((3, 3))
    energy = None
    spin_pol_inp = False
    spin_pol_out = False

    with open(filename, 'rU') as fileobj:
        lines = fileobj.readlines()

    n_mag_inp = []

    for n, line in enumerate(lines):
        if 'positions (alat units)' in line:
            atomic_pos_units = 'alat'
            n_pos = n+1
        elif 'ATOMIC_POSITIONS' in line and 'crystal' in line:
            atomic_pos_units = 'crystal'
            n_pos = n+1
        elif 'ATOMIC_POSITIONS' in line and 'angstrom' in line:
            atomic_pos_units = 'angstrom'
            n_pos = n+1
        elif 'celldm(1)' in line:
            celldm = float(line.split()[1]) * units['Bohr']
        elif 'crystal axes: (cart. coord. in units of alat)' in line:
            cell_units = 'alat'
            n_cell = n+1
        elif 'CELL_PARAMETERS' in line and 'angstrom' in line:
            cell_units = 'angstrom'
            n_cell = n+1
        elif '!' in line:
            n_nrg = n
        elif 'Final energy' in line:
            n_fin = n
        elif 'Starting magnetic structure' in line and spin_pol_out is False:
            spin_pol_out = True
            n_mag_out = n+2
        elif 'starting_magnetization' in line:
            spin_pol_inp = True
            n_mag_inp += [n]
        elif 'ATOMIC_SPECIES' in line:
            n_atom_list = n+1

    for i in range(3):
        line = lines[n_cell+i]
        if cell_units == 'alat':
            cell[i] = [float(c)*celldm for c in line.split()[3:6]]
        elif cell_units == 'angstrom':
            cell[i] = [float(c) for c in line.split()[:3]]

    atoms.set_cell(cell)

    try: energy = float(lines[n_nrg].split()[4])*units['Ry']
    except: energy = None

    try: energy = float(lines[n_fin].split()[3])*units['Ry']
    except: pass

    index = 0
    indices = []
    constraints = []
    translate_constraints = {0: True, 1: False}

    magmoms_dict = {}

    atoms_list = []

    if spin_pol_inp is True:

        for line in lines[n_atom_list:]:
    
            if len(line.split()) == 0:
                break

            atoms_list += [line.split()[0]]

        for n in n_mag_inp:
    
            num = ''
            read = False

            for i in lines[n]:
                if i == ')':
                    read = False
                if read is True:
                    num += i
                if i == '(':
                    read = True

            magmoms_dict[atoms_list[int(num)-1]] = float(lines[n].split()[2])

    if spin_pol_out is True:

        for line in lines[n_mag_out:]:
    
            if len(line.split()) == 0 or line.split()[0] == 'End':
                break
    
            magmoms_dict[line.split()[0]] = float(line.split()[1])

    for line in lines[n_pos:]:

        if len(line.split()) == 0 or line.split()[0] == 'End':
            break

        if atomic_pos_units == 'alat':
            name = line.split()[1]
            positions = [[float(i)*celldm for i in line.split()[6:9]]]
            fix = [False, False, False]
        else:
            name = line.split()[0]
            positions = [[float(i) for i in line.split()[1:4]]]
            fix = [translate_constraints[int(i)] for i in line.split()[4:]]

        symbol = ''
        magmom_tag = ''

        for i in range(len(name)):
            if name[i].isdigit():
                magmom_tag += name[i]
            else:
                symbol += name[i]
        
        if spin_pol_inp is True or spin_pol_out is True:
            magmom = magmoms_dict[name]
            magmom *= SSSP_VALENCE[atomic_numbers[symbol]]
            magmoms = [magmom]

        else:
            magmoms = [0.]*len(positions)

        if atomic_pos_units == 'crystal':
            atoms += Atoms(symbols          = symbol   ,
                           scaled_positions = positions,
                           magmoms          = magmoms  )
        else:
            atoms += Atoms(symbols   = symbol   ,
                           positions = positions,
                           magmoms   = magmoms  )

        if fix == [True, True, True]:
            indices.append(index)
        elif True in fix:
            constraints.append(FixCartesian([index], fix))
        index += 1

    constraints.append(FixAtoms(indices = indices))
    atoms.set_constraint(constraints)
    results = {'energy' : energy}
    calc = Espresso()
    calc.results = results
    atoms.set_calculator(calc)
    atoms.potential_energy = energy

    return atoms

################################################################################
# READ QUANTUM ESPRESSO OUT
################################################################################

class ReadQeOut:

    def __init__(self, filename):

        self.filename         = filename
        self.atoms            = None
        self.potential_energy = None

    def get_atoms(self, cell = None):

        atoms = read_qe_out(self.filename)
        if cell is not None:
            atoms.set_cell(cell)
        self.atoms = atoms

        return atoms

    def get_potential_energy(self):

        if self.atoms is None:
            atoms = self.get_atoms()
        else:
            atoms = self.atoms

        try: potential_energy = atoms.potential_energy
        except: potential_energy = atoms.calc.results['energy']

        self.potential_energy = potential_energy

        return potential_energy

    def read_bands(self, scale_band_energies = True):

        e_bands_dict, e_fermi = read_pw_bands(self.filename      ,
                                              scale_band_energies)

        self.e_bands_dict = e_bands_dict
        self.e_fermi      = e_fermi

        return e_bands_dict, e_fermi

################################################################################
# READ QUANTUM ESPRESSO INP
################################################################################

def read_qe_inp(filename):

    with open(filename, 'rU') as fileobj:
        lines = fileobj.readlines()

    n_as = 0
    n_kp = 0
    gamma = False
    for n, line in enumerate(lines):
        if 'ATOMIC_SPECIES' in line:
            n_as = n
        elif 'K_POINTS' in line:
            if 'gamma' in line:
                gamma = True
            n_kp = n

    input_data = {}
    for n, line in enumerate(lines):
        if ('ATOMIC_SPECIES' in line 
            or 'ATOMIC_POSITIONS' in line
            or 'K_POINTS' in line 
            or 'CELL_PARAMETERS' in line):
            break
        if len(line.strip()) == 0 or line == '\n':
            pass
        elif line[0] in ('&', '/'):
            pass
        else:
            keyword, argument = line.split('=')
            keyword = re.sub(re.compile(r'\s+'), '', keyword)
            argument = re.sub(re.compile(r'\s+'), '', argument)
            if '.true.' in argument:
                argument = True
            elif '.false.' in argument:
                argument = False
            else:
                argument = ast.literal_eval(argument)
            if type(argument) is tuple: argument = argument[0]
            input_data[keyword] = argument

    pseudos = {}
    for n, line in enumerate(lines[n_as+1:]):
        if len(line.strip()) == 0 or line == '\n':
            break
        element, MW, pseudo = line.split()
        pseudos[element] = pseudo

    if gamma:
        kpts = (1,1,1)
        koffset = (0,0,0)
    else:
        kpts = [ int(i) for i in lines[n_kp+1].split()[:3] ]
        koffset = [ int(i) for i in lines[n_kp+1].split()[3:] ]

    return input_data, pseudos, kpts, koffset

################################################################################
# READ QUANTUM ESPRESSO INP
################################################################################

class ReadQeInp:

    def __init__(self, filename):

        self.filename = filename
        self.atoms    = None
        self.label    = filename.split('.')[0]

    def get_data_pseudos_kpts(self):

        input_data, pseudos, kpts, koffset = read_qe_inp(self.filename)

        self.input_data = input_data
        self.pseudos    = pseudos
        self.kpts       = kpts
        self.koffset    = koffset

        return input_data, pseudos, kpts, koffset

    def get_calculator(self):
        
        input_data, pseudos, kpts, koffset = self.get_data_pseudos_kpts()
        
        calc = Espresso(input_data       = input_data,
                        pseudopotentials = pseudos   ,
                        kpts             = kpts      ,
                        koffset          = koffset   ,
                        label            = self.label)
        
        return calc

    def get_atoms(self):

        atoms = read_qe_out(self.filename)

        self.atoms = atoms

        return atoms

    def get_input_data_dicts(self, remove_keywords = True):

        hubbard_U_dict    = {}
        hubbard_J0_dict   = {}
        init_charges_dict = {}

        if self.atoms is None:
            self.get_atoms()

        symbols_list = get_symbols_list(self.atoms)

        del_keywords = []

        for keyword in self.input_data:

            if 'Hubbard_U' in keyword:

                n = int(keyword.split('(', ')')[1])
                symbol = symbols_list[n-1]
                hubbard_U_dict[symbol] = self.input_data[keyword]
                del_keywords += [keyword]

            elif 'Hubbard_J0' in keyword:

                n = int(keyword.split('(', ')')[1])
                symbol = symbols_list[n-1]
                hubbard_J0_dict[symbol] = self.input_data[keyword]
                del_keywords += [keyword]
                
            elif 'starting_charge' in keyword:

                n = int(re.split(r'\(|\)', keyword)[1])
                symbol = symbols_list[n-1]
                init_charges_dict[symbol] = self.input_data[keyword]
                del_keywords += [keyword]

                if ('tot_charge' in self.input_data
                    and 'tot_charge' not in del_keywords):
                    del_keywords += ['tot_charge']

        if remove_keywords is True:
            for keyword in del_keywords:
                del self.input_data[keyword]

        self.hubbard_U_dict    = hubbard_U_dict
        self.hubbard_J0_dict   = hubbard_J0_dict
        self.init_charges_dict = init_charges_dict

################################################################################
# WRITE NEB DAT
################################################################################

def write_neb_dat(neb_data, filename = 'neb.dat', mode = 'w+'):

    neb_dict = {
        'string_method'        : None,
        'restart_mode'         : None,
        'nstep_path'           : None,
        'num_of_images'        : None,
        'opt_scheme'           : None,
        'CI_scheme'            : None,
        'first_last_opt'       : None,
        'minimum_image'        : None,
        'temp_req'             : None,
        'ds'                   : None,
        'k_max'                : None,
        'k_min'                : None,
        'path_thr'             : None,
        'use_masses'           : None,
        'use_freezing'         : None,
        'lfcpopt'              : None,
        'fcp_mu'               : None,
        'fcp_tot_charge_first' : None,
        'fcp_tot_charge_last'  : None,
    }

    for arg in neb_data:
        neb_dict[arg] = neb_data[arg]

    with open(filename, mode) as f:

        f.write('&PATH\n')
        for arg in [arg for arg in neb_dict if neb_dict[arg] is not None]:
            if isinstance(neb_dict[arg], str):
                neb_dict[arg] = '\''+neb_dict[arg]+'\''
            elif neb_dict[arg] is True:
                neb_dict[arg] = '.true.'
            elif neb_dict[arg] is False:
                neb_dict[arg] = '.false.'
            f.write('   {0} = {1}\n'.format(str(arg).ljust(16), neb_dict[arg]))
        f.write('/')

################################################################################
# WRITE NEB INP
################################################################################

def write_neb_inp(neb_data, images, calc, filename = 'neb.pwi'):

    calc = cp.deepcopy(calc)
    calc.label = 'tmp'

    with open(filename, 'w+') as f:
        f.write('BEGIN\n')
        f.write('BEGIN_PATH_INPUT\n')

    write_neb_dat(neb_data, filename, mode = 'a+')

    with open(filename, 'a+') as f:

        f.write('\nEND_PATH_INPUT\n')
        f.write('BEGIN_ENGINE_INPUT\n')

        for i in range(len(images)):

            calc.write_input(images[i])

            with open('tmp.pwi', 'rU') as g:
                lines = g.readlines()

            for n, line in enumerate(lines):
                if 'ATOMIC_POSITIONS' in line:
                    break

            if i == 0:
                for line in lines[:n]:
                    f.write(line)
                f.write('BEGIN_POSITIONS\n')
                f.write('FIRST_IMAGE\n')
            elif i == len(images)-1:
                f.write('LAST_IMAGE\n')
            else:
                f.write('INTERMEDIATE_IMAGE\n')

            for line in lines[n:]:
                f.write(line)

        os.remove('tmp.pwi')

        f.write('END_POSITIONS\n')
        f.write('END_ENGINE_INPUT\n')
        f.write('END\n')

################################################################################
# READ NEB CRD
################################################################################

def read_neb_crd(images, filename = 'pwscf.crd'):

    with open(filename, 'rU') as fileobj:
        lines = fileobj.readlines()
    
    n_atoms = len(images[0])

    num = 2
    for image in images:
        positions = []
        for line in lines[num:num+n_atoms]:
            positions.append(line.split()[1:4])
        image.set_positions(positions)
        num += n_atoms+2

    return images

################################################################################
# UPDATE PSEUDOPOTENTIALS
################################################################################

def update_pseudos(pseudos, filename):

    _, pseudos_new, _, _ = read_qe_inp(filename)
    pseudos.copy()
    pseudos.update(pseudos_new)
    
    return pseudos

################################################################################
# READ PW BANDS
################################################################################

def read_pw_bands(filename = 'pw.pwo', scale_band_energies = True):

    with open(filename, 'rU') as fileobj:
        lines = fileobj.readlines()

    kpt = 0
    e_bands_dict = {}
    n_kpts = 0
    kpts_list = []
    n_spin = 1
    read_bands = False

    for line in lines:
        if 'End of self-consistent calculation' in line:
            kpt = 0
            e_bands_dict = {}
            kpts_list = []
        if 'number of k points' in line:
            n_kpts = int(line.split()[4])
        if 'SPIN UP' in line:
            n_spin = 2
        if ' k =' in line:
            read_bands = True
            count = 0
            kpt += 1
            e_bands_dict[kpt] = []
            kpts_list += [kpt]
        if read_bands is True:
            if count == 1:
                for i in range(8):
                    if len(line) > 9*(i+1)+2:
                        e_bands_dict[kpt] += [float(line[9*i+2:9*(i+1)+2])]
            if len(line.strip()) == 0 or line == '\n':
                count += 1
        if 'the Fermi energy is' in line:
            e_fermi = float(line.split()[4])

    n_kpts *= n_spin

    #assert n_kpts == kpts_list[-1]

    if scale_band_energies is True:
        for kpt in e_bands_dict:
            for i in range(len(e_bands_dict[kpt])):
                e_bands_dict[kpt][i] -= e_fermi

    return e_bands_dict, e_fermi

################################################################################
# CREATE PP_INP
################################################################################

def create_pp_inp(filename   = 'pp.pwi',
                  outname    = 'pw.pwo',
                  band_num   = 'homo',
                  datname    = 'charge',
                  pp_data    = {},
                  plot_data  = {},
                  delta_e    = 0.001,
                  kpts       = None,
                  print_summ = False,
                  ):

    e_bands_dict, e_fermi = read_pw_bands(filename = outname)
    kpts_list = [kpt for kpt in e_bands_dict]
    n_bands = max([len(e_bands_dict[kpt]) for kpt in e_bands_dict])

    if kpts is not None:
        kpts_list = kpts

    n_kpts = len(kpts_list)
    n_min = n_bands
    n_max = 0

    for kpt in e_bands_dict:

        e_bands = e_bands_dict[kpt]

        if band_num == 'homo':
            for i in range(len(e_bands)):
                if e_bands[i] < e_fermi and e_fermi-e_bands[i] < delta_e:
                    if i > n_max:
                        n_max = i
                    if i < n_min:
                        n_min = i

        elif band_num == 'lumo':
            for i in range(len(e_bands)):
                if e_bands[i] > e_fermi and e_bands[i]-e_fermi < delta_e:
                    if i > n_max:
                        n_max = i
                    if i < n_min:
                        n_min = i

        else:
            n_min, n_max = band_num

    pp_data['filplot']   = datname
    pp_data['kpoint(1)'] = 1

    if n_kpts > 1:
        pp_data['kpoint(2)'] = n_kpts
    if n_max == n_min:
        pp_data['kband(1)'] = n_max
    else:
        pp_data['kband(1)'] = n_min
        pp_data['kband(2)'] = n_max

    plot_data['nfile']     = 1
    plot_data['filepp(1)'] = datname

    write_pp_input(pp_data   = pp_data  ,
                   plot_data = plot_data,
                   filename  = filename )

    if print_summ is True:
        print('number of bands = {}\n'.format(n_bands))

################################################################################
# MERGE CHARGE FILES
################################################################################

def merge_charge_files(files_in, file_out):

    for num, filename in enumerate(files_in):

        with open(filename, 'rU') as fileobj:
            lines = fileobj.readlines()

        if num == 0:
            new_lines = cp.deepcopy(lines)
            for n, line in enumerate(lines):
                if 'BEGIN_BLOCK_DATAGRID_3D' in line:
                    n_den = n
                if 'END_DATAGRID_3D' in line:
                    n_end = n
            n_grid = [int(n) for n in lines[n_den+3].split()]
            n_points = n_grid[0]*n_grid[1]*n_grid[2]
            density = np.zeros(n_points)

        i = 0
        for line in lines[n_den+8:n_end]:
            for l in line.split():
                density[i] += float(l)
                i += 1

    with open(file_out, 'w+') as f:

        for line in new_lines[:n_den+8]:
            f.write(line)

        for i in range(n_points):
            if i != 0 and i % 6 == 0:
                f.write('\n')
            f.write('{:14.6E}'.format(density[i]))

        f.write('\n')
        for line in new_lines[n_end:]:
            f.write(line)

################################################################################
# CLASS STATE
################################################################################

class State:

    def __init__(self, state_num, atom_num, element, shell_num, l, m):

        self.state_num = state_num
        self.element   = element
        self.atom_num  = atom_num
        self.shell_num = shell_num
        self.l         = l
        self.m         = m

        p_dict = {1: 'p z', 2: 'p x', 3: 'p y'}
        d_dict = {1: 'd z^2', 2: 'd zx', 3: 'd zy', 4: 'd x^2-y^2', 5: 'd xy'}

        if l == 0:
            self.orbital = 's'

        elif l == 1:
            self.orbital = p_dict[self.m]

        elif l == 2:
            self.orbital = d_dict[self.m]

################################################################################
# CLASS BAND
################################################################################

class Band:

    def __init__(self, band_num, energy, state_nums, weights):

        self.band_num   = band_num
        self.energy     = energy
        self.state_nums = state_nums
        self.weights    = weights

################################################################################
# CLASS ATOM PP
################################################################################

class AtomPP:

    def __init__(self,
                 atom_num, 
                 element,
                 states  = [],
                 bands   = [],
                 weights = [],
                 color   = None,
                 ):

        self.atom_num = atom_num
        self.element  = element
        self.states   = states
        self.bands    = bands
        self.weights  = weights
        self.color    = color

################################################################################
# READ PROJWFC
################################################################################

def read_projwfc(filename, kpoint, print_summary = False):

    with open(filename, 'rU') as fileobj:
        lines = fileobj.readlines()

    states_list = []
    bands_list = []

    read = False
    kpt = 0
    band_num = 0

    for line in lines:

        if 'state #' in line:
            state_num = int(line[12:16])
            atom_num = int(line[22:26])
            element = line[28:30].strip()
            shell_num = int(line[38:40])
            l = int(line[44:45])
            m = int(line[48:50])
            state = State(state_num = state_num,
                          atom_num  = atom_num ,
                          element   = element  ,
                          shell_num = shell_num,
                          l         = l        ,
                          m         = m        )
            states_list += [state]

        if ' k = ' in line:
            kpt += 1

        if kpoint == kpt:
            if '    |psi|^2' in line:
                read = False
                band = Band(band_num   = band_num  ,
                            energy     = energy    ,
                            state_nums = state_nums,
                            weights    = weights   )
                bands_list += [band]

            if read is True:
                for i in range(5):
                    weight = line[11+14*i:16+14*i]
                    state_num = line[19+14*i:23+14*i]
                    try:
                        weights += [float(weight)]
                        state_nums += [int(state_num)]
                    except: pass

            if '==== e(' in line:
                read = True
                band_num = int(line[7:11])
                energy = float(line[14:26])
                weights = []
                state_nums = []

            if '     e =' in line:
                read = True
                band_num += 1
                energy = float(line[8:20])
                weights = []
                state_nums = []

    if print_summary is True:
        print(f'n states = {len(states_list)} - n bands = {len(bands_list)}')

    return states_list, bands_list

################################################################################
# SCALE BAND ENERGIES
################################################################################

def scale_band_energies(band_list, e_fermi):

    for band in band_list:
        band.energy -= e_fermi

    return band_list

################################################################################
# GET ATOMS DETAILS
################################################################################

def get_atoms_details(states_list,
                      bands_list,
                      atom_num_list,
                      delta_e    = 0.,
                      color_dict = None,
                      ):

    atoms_pp_list = []
    
    for atom_num in atom_num_list:

        states = [s for s in states_list if s.atom_num == atom_num]
        element = states[0].element
        color = color_dict[element] if color_dict is not None else None

        atom = AtomPP(atom_num = atom_num,
                      element  = element ,
                      states   = []      ,
                      bands    = []      ,
                      weights  = []      ,
                      color    = color   )
        atoms_pp_list += [atom]

        for band in bands_list:
            for state in states:
                if state.state_num in band.state_nums:
                    index = band.state_nums.index(state.state_num)
                    weight = band.weights[index]
                    if weight > delta_e:
                        atom.states += [state]
                        atom.bands += [band]
                        atom.weights += [weight]

    return atoms_pp_list

################################################################################
# PRINT ATOMS DETAILS
################################################################################

def print_atoms_details(atoms_pp_list, filename = 'atom_details.out'):

    fileobj = open(filename, 'w+')

    count_bands = {}

    for atom in atoms_pp_list:

        atom_num = atom.atom_num

        print('\n atom {0} {1}'.format(atom.element, str(atom_num)),
              file = fileobj)
        print('| state num | orbital type | band num | weight |  energy  |',
              file = fileobj)

        for i in range(len(atom.states)):

            state = atom.states[i]
            band = atom.bands[i]
            weight = atom.weights[i]

            print(' {0:10d}   {1:9s} {2:13d} {3:8.4f}   {4:+8.3f}'.format(
                  state.state_num, state.orbital, band.band_num, weight,
                  band.energy), file = fileobj)

            try:
                if state.atom_num not in count_bands[band.band_num]:
                    count_bands[band.band_num] += [state.atom_num]
            except:
                count_bands[band.band_num] = [state.atom_num]

    print('\n SHARED BANDS \n', file = fileobj)
    
    for num in sorted([n for n in count_bands if len(count_bands[n]) > 1]):
    
        all_bands = sum([a.bands for a in atoms_pp_list], [])
    
        band = [b for b in all_bands if b.band_num == num][0]
    
        string = 'band {0:4d} ({1:+8.3f}) shared by atoms: '.format(
                  num, band.energy)

        for atom_num in sorted(count_bands[num]):

            atom = [a for a in atoms_pp_list if a.atom_num == atom_num][0]
            string += ' {0:2s}{1:3d}  - '.format(atom.element, atom_num)

        print(string[:-2], file = fileobj)

    print('', file = fileobj)

    fileobj.close()

################################################################################
# PLOT ENERGY LEVELS
################################################################################

def plot_band_levels(atoms_pp_list,
                     num_min_print,
                     bands_energies,
                     e_min,
                     e_max,
                     ):

    import matplotlib.pyplot as plt

    fig = plt.figure(0)

    if bands_energies is not None:

        for energy in bands_energies:

            plt.plot([0., 1.], [energy]*2, color = 'whitesmoke')

    weight_cum = {}

    for atom in atoms_pp_list:
    
        color = atom.color
    
        for i in range(len(atom.bands)):
    
            band = atom.bands[i]
            weight = atom.weights[i]
    
            try:
                w_min = weight_cum[band]
                weight_cum[band] += weight
            except:
                w_min = 0.
                weight_cum[band] = weight
    
            w_max = weight_cum[band]
    
            plt.plot([w_min, w_max], [band.energy]*2, color = color)

    for band in weight_cum:

        if (weight_cum[band] > num_min_print and e_min < band.energy < e_max):

            x_text = 1.01
            y_text = band.energy-0.05

            plt.text(x_text, y_text, band.band_num, size = 'xx-small')

    plt.xlim([0., 1.])
    plt.ylim([e_min, e_max])

    return fig

################################################################################
# GET DOS
################################################################################

def get_dos(filename):

    with open(filename, 'rU') as fileobj:
        lines = fileobj.readlines()

    if 'dosup(E)' in lines[0]:
        nspin = 2
    else:
        nspin = 1

    energy = np.zeros(len(lines)-1)
    if nspin == 2:
        dos = np.zeros((len(lines)-1, 2))
    else:
        dos = np.zeros(len(lines)-1)

    e_fermi = float(lines[0].split()[-2])
    for i, line in enumerate(lines[1:]):
        energy[i] = float(line.split()[0])-e_fermi
        if nspin == 2:
            dos[i, 0] = float(line.split()[1])
            dos[i, 1] = float(line.split()[2])
        else:
            dos[i] = float(line.split()[1])

    return energy, dos

################################################################################
# GET PDOS
################################################################################

def get_pdos(filename, e_fermi):

    with open(filename, 'rU') as fileobj:
        lines = fileobj.readlines()

    if 'dosup(E)' in lines[0] or 'ldosup(E)' in lines[0]:
        nspin = 2
    else:
        nspin = 1

    energy = np.zeros(len(lines)-1)
    if nspin == 2:
        pdos = np.zeros((len(lines)-1, 2))
    else:
        pdos = np.zeros(len(lines)-1)

    for i, line in enumerate(lines[1:]):
        energy[i] = float(line.split()[0])-e_fermi
        if nspin == 2:
            pdos[i, 0] = float(line.split()[1])
            pdos[i, 1] = float(line.split()[2])
        else:
            pdos[i] = float(line.split()[1])

    return energy, pdos

################################################################################
# GET PDOS VECT
################################################################################

def get_pdos_vect(atoms, e_fermi, filename = 'projwfc.pwo'):

    states_list, _ = read_projwfc(filename = filename,
                                  kpoint   = None    )

    pdos_vect = np.array([None]*len(atoms), dtype = object)
    names_list = []
    for state in states_list:
        atom_num = state.atom_num
        orbital_type = state.orbital[:1]
        name = 'pdos.pdos_atm#{0}({1})_wfc#{2}({3})'.format(atom_num       ,
                                                            state.element  ,
                                                            state.shell_num,
                                                            orbital_type   )
        if name not in names_list:
            names_list += [name]
            energy, pdos = get_pdos(filename = name   ,
                                    e_fermi  = e_fermi)
            if pdos_vect[atom_num-1] is None:
                pdos_vect[atom_num-1] = {}
            if orbital_type in pdos_vect[atom_num-1]:
                pdos_vect[atom_num-1][orbital_type] += pdos
            else:
                pdos_vect[atom_num-1][orbital_type] = pdos

    return energy, pdos_vect

################################################################################
# GET FEATURES BANDS
################################################################################

def get_features_bands(atoms, 
                       energy,
                       pdos_vect,
                       delta_e     = 0.1, 
                       save_pickle = True,
                       ):

    i_zero = np.argmin(np.abs(energy))
    i_minus = np.argmin(np.abs(energy+delta_e))
    i_plus = np.argmin(np.abs(energy-delta_e))

    features = np.zeros((len(atoms),8))

    for i, _ in enumerate(atoms):
        pdos_dict = pdos_vect[i]
        for orbital in pdos_dict:
            if len(pdos_dict[orbital].shape) > 1:
                pdos_dict[orbital] = np.sum(pdos_dict[orbital], axis = 1)
        
        pdos_sp = pdos_dict['s']
        pdos_sp += pdos_dict['p']
        pdos_sp = pdos_dict['p']
        sp_filling = np.trapz(y = pdos_sp[:i_zero], x = energy[:i_zero])
        sp_density = (np.sum(pdos_sp[i_minus:i_plus]) /
                      len(pdos_sp[i_minus:i_plus]))
        
        if 'd' in pdos_dict:
            pdos_d = pdos_dict['d']
            d_filling = np.trapz(y = pdos_d[:i_zero], x = energy[:i_zero])
            d_density = (np.sum(pdos_d[i_minus:i_plus]) / 
                         len(pdos_d[i_minus:i_plus]))
            d_centre = np.trapz(pdos_d*energy, energy)/np.trapz(pdos_d, energy)
            d_mom_2 = (np.trapz(pdos_d*np.power(energy-d_centre,2), energy) / 
                       np.trapz(pdos_d, energy))
            d_width = np.sqrt(d_mom_2)
            d_mom_3 = (np.trapz(pdos_d*np.power(energy-d_centre,3), energy) / 
                       np.trapz(pdos_d, energy))
            d_skewness = d_mom_3/np.power(d_width,3)
            d_mom_4 = (np.trapz(pdos_d*np.power(energy-d_centre,4), energy) / 
                       np.trapz(pdos_d, energy))
            d_kurtosis = d_mom_4/np.power(d_width,4)
        else:
            d_filling  = np.nan
            d_density  = np.nan
            d_centre   = np.nan
            d_width    = np.nan
            d_skewness = np.nan
            d_kurtosis = np.nan
        
        features[i,0] = d_filling
        features[i,1] = d_centre
        features[i,2] = d_width
        features[i,3] = d_skewness
        features[i,4] = d_kurtosis
        features[i,5] = sp_filling
        features[i,6] = d_density
        features[i,7] = sp_density
        
        if save_pickle is True:
            with open('features_bands.pickle', 'wb') as fileobj:
                pickle.dump(features, fileobj)

    return features

################################################################################
# WRITE FEATURES OUT
################################################################################

def write_features_out(atoms, features_names, features, filename):

    with open(filename, 'w+') as fileobj:
        print("Calculated Features", file=fileobj)
        assert len(features_names) == features.shape[1]
        print(f'{"symbol":7s}', end='', file=fileobj)
        for i in range(features.shape[1]):
            print(f'  {features_names[i]:11s}', end='', file=fileobj)
        print('', file=fileobj)
        for i in range(features.shape[0]):
            print(f'{atoms[i].symbol:7s}', end='', file=fileobj)
            for feature in features[i,:]:
                print(f'{feature:+13.4e}', end='', file=fileobj)
            print('', file=fileobj)

################################################################################
# ASSIGN HUBBARD U
################################################################################

def assign_hubbard_U(atoms, pw_data, hubbard_U_dict, hubbard_J0_dict = None):

    symbols_list = get_symbols_list(atoms)
    for i in range(len(symbols_list)):
        symbol = symbols_list[i]
        pw_data['Hubbard_U({})'.format(i+1)] = hubbard_U_dict[symbol]
        if hubbard_J0_dict is not None:
            pw_data['Hubbard_J0({})'.format(i+1)] = hubbard_J0_dict[symbol]
    pw_data['lda_plus_u'] = True

    return pw_data

################################################################################
# ASSIGN INIT CHARGES
################################################################################

def assign_init_charges(atoms, pw_data, init_charges_dict):

    symbols_list = get_symbols_list(atoms)
    symbols_dict = get_symbols_dict(atoms)

    i = 0
    charge_dict = {}
    tot_charge = 0.
    for symbol in symbols_list:
        charge = init_charges_dict[symbol]
        if charge != 0.:
            charge_dict['starting_charge({})'.format(i+1)] = charge
            tot_charge += symbols_dict[symbol]*charge
        i += 1

    pw_data['tot_charge'] = tot_charge
    if 'system' in pw_data:
        pw_data['system'].update(charge_dict)
    else:
        pw_data['system'] = charge_dict

    return pw_data

################################################################################
# CREATE EOS INPUTS
################################################################################

def create_eos_inputs(atoms, delta_x, npoints, run_cmd = None):

    cell_zero = atoms.get_cell()
    x = 1.-delta_x
    for i in range(npoints):
    
        dirname = 'volumes_{:02d}'.format(i)
        try: os.mkdir(dirname)
        except: pass
    
        os.chdir(dirname)
        v = x**(1./3.)
        atoms.set_cell(v*cell_zero, scale_atoms = True)
        calc = atoms.get_calculator()
        calc.write_input(atoms)
        x += (2.*delta_x)/(npoints-1)
        if run_cmd is not None:
            os.system(run_cmd)
        os.chdir('..')

################################################################################
# READ EOS OUTPUTS
################################################################################

def read_eos_outputs(npoints, filename = 'pw.pwo', get_cells = False):

    volumes  = []
    energies = []
    cells    = []
    
    for i in range(npoints):
    
        dirname = 'volumes_{:02d}'.format(i)
    
        os.chdir(dirname)
        atoms = read_qe_out(filename)
        cells.append(atoms.cell)
        volumes.append(atoms.get_volume())
        energies.append(atoms.potential_energy)
        os.chdir('..')

    with open('volumes_energies.txt', 'w+') as f:
        f.write('volumes     energies\n')
        for i in range(len(volumes)):
            f.write('{0:10.7f}  {1:14.7f}\n'.format(volumes[i], energies[i]))

    if get_cells is True:
        return volumes, energies, cells
    else:
        return volumes, energies

################################################################################
# REORDER NEB IMAGES
################################################################################

def reorder_neb_images(first, last):

    coupled = cp.deepcopy(last)
    spared = [a for a in last]
    
    del coupled [range(len(coupled))]
    
    for a in first:
        distances = [10.]*len(last)
        for b in [ b for b in last if b.symbol == a.symbol ]:
            distances[b.index] = np.linalg.norm(b.position-a.position)
        if np.min(distances) > 0.5:
            first += first.pop(i = a.index)
    
    for a in first:
        distances = [10.]*len(last)
        for b in [ b for b in last if b.symbol == a.symbol ]:
            distances[b.index] = np.linalg.norm(b.position-a.position)
        if np.min(distances) < 0.5:
            index = np.argmin(distances)
            coupled += last[index]
            spared[index] = None
    
    spared = Atoms([a for a in spared if a is not None])
    last = coupled+spared

    return first, last

################################################################################
# GET NUMBER OF ELECTRONS
################################################################################

def get_number_of_electrons(atoms):

    n_electrons = 0
    for a in atoms:
        n_electrons += SSSP_VALENCE[atomic_numbers[a.symbol]]

    return n_electrons

################################################################################
# WRITE QE INPUT BLOCK
################################################################################

def write_qe_input_block(fileobj, block_name, block_data, col = 23):

    print('&'+block_name, file = fileobj)

    for arg in [arg for arg in block_data if block_data[arg] is not None]:

        if type(block_data[arg]) == str:
            string = '\''+block_data[arg]+'\''
        elif block_data[arg] is True:
            string = '.true.'
        elif block_data[arg] is False:
            string = '.false.'
        else:
            string = block_data[arg]

        print('   {0} = {1}'.format(arg.ljust(col), string), file = fileobj)

    print('/', file = fileobj)

################################################################################
# WRITE DOS INPUT
################################################################################

def write_dos_input(dos_data, filename = 'dos.pwi', col = 23):

    with open(filename, 'w+') as fileobj:
        write_qe_input_block(fileobj    = fileobj   ,
                             block_name = 'DOS'     ,
                             block_data = dos_data  ,
                             col        = col       )

################################################################################
# WRITE PP INPUT
################################################################################

def write_pp_input(pp_data, plot_data, filename = 'pp.pwi', col = 23):

    with open(filename, 'w+') as fileobj:
    
        write_qe_input_block(fileobj    = fileobj   ,
                             block_name = 'INPUTPP' ,
                             block_data = pp_data   ,
                             col        = col       )

        write_qe_input_block(fileobj    = fileobj   ,
                             block_name = 'PLOT'    ,
                             block_data = plot_data ,
                             col        = col       )

################################################################################
# WRITE PROJWFC INPUT
################################################################################

def write_projwfc_input(proj_data, filename = 'projwfc.pwi', col = 23):

    with open(filename, 'w+') as fileobj:
        write_qe_input_block(fileobj    = fileobj   ,
                             block_name = 'PROJWFC' ,
                             block_data = proj_data ,
                             col        = col       )

################################################################################
# WRITE DOS INPUT
################################################################################

def write_environ_input(env_dict,
                        bon_dict,
                        ele_dict,
                        filename = 'environ.in',
                        reg_list = [],
                        col      = 23,
                        ):

    fileobj = open(filename, 'w+')

    block_name = 'ENVIRON'

    write_qe_input_block(fileobj    = fileobj   ,
                         block_name = block_name,
                         block_data = env_dict  ,
                         col        = col       )

    block_name = 'BOUNDARY'

    write_qe_input_block(fileobj    = fileobj   ,
                         block_name = block_name,
                         block_data = bon_dict  ,
                         col        = col       )

    block_name = 'ELECTROSTATIC'

    write_qe_input_block(fileobj    = fileobj   ,
                         block_name = block_name,
                         block_data = ele_dict  ,
                         col        = col       )

    if len(reg_list) > 0:

        print('\nDIELECTRIC_REGIONS angstrom', file = fileobj)

        for reg in reg_list:
            
            print('{:.4f}'.format(reg.eps_stat), end = ' ', file = fileobj)
            print('{:.4f}'.format(reg.eps_opt), end = ' ', file = fileobj)
            print('{:.14f}'.format(reg.position[0]), end = ' ', file = fileobj)
            print('{:.14f}'.format(reg.position[1]), end = ' ', file = fileobj)
            print('{:.14f}'.format(reg.position[2]), end = ' ', file = fileobj)
            print('{:.14f}'.format(reg.width), end = ' ', file = fileobj)
            print('{:.4f}'.format(reg.spread), end = ' ', file = fileobj)
            print(reg.dim, end = ' ', file = fileobj)
            print(reg.axis, file = fileobj)
            
    fileobj.close()

################################################################################
# DIELECTRIC REGION
################################################################################

class DielectricRegion:

    def __init__(self,
                 eps_stat = None,
                 eps_opt  = None,
                 position = None,
                 width    = None,
                 spread   = None,
                 dim      = None,
                 axis     = None):

        self.eps_stat = eps_stat
        self.eps_opt  = eps_opt
        self.position = position
        self.width    = width
        self.spread   = spread
        self.dim      = dim
        self.axis     = axis

################################################################################
# END
################################################################################
