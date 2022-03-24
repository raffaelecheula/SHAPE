################################################################################
# Raffaele Cheula, cheula.raffaele@gmail.com
################################################################################

import numpy as np
from math import sin, pi, sqrt
from ase import Atom, Atoms
from ase.units import kB
from ase.io.espresso import SSSP_VALENCE
from ase.parallel import paropen

################################################################################
# GET SYMBOLS LIST
################################################################################

def get_symbols_list(atoms, check_magmoms = True):

    symbols = atoms.get_chemical_symbols()
    magmoms = [a.magmom for a in atoms]
    
    if check_magmoms is True:
        magmoms = [a.magmom for a in atoms]
    else:
        magmoms = [0. for a in atoms]
    
    if len(symbols) > 1:
        for i in range(len(symbols)-1, 0, -1):
            for j in range(i):
                if symbols[j] == symbols[i] and magmoms[i] == magmoms[j]:
                    del symbols[i]
                    break

    return symbols

################################################################################
# GET ATOM LIST
################################################################################

def get_atom_list(atoms):

    return get_symbols_list(atoms)

################################################################################
# GET SYMBOLS DICT
################################################################################

def get_symbols_dict(atoms):

    symbols_dict = {}
    symbols = atoms.get_chemical_symbols()
    for symbol in symbols:
        if symbol in symbols_dict:
            symbols_dict[symbol] += 1
        else:
            symbols_dict[symbol] = 1

    return symbols_dict

################################################################################
# GET FORMULA REPETITIONS
################################################################################

def get_formula_repetitions(atoms):

    symbols_dict = get_symbols_dict(atoms)
    formula_repetitions = min([symbols_dict[i] for i in symbols_dict])

    return formula_repetitions

################################################################################
# ARRAY FROM DICT
################################################################################

def array_from_dict(symbol, array_dict):

    array = []*len(symbol)
    for i in range(len(symbol)):
        if symbol[i] in array_dict:
            array[i] = array_dict[symbol[i]]

    return array

################################################################################
# ATOMS FIXED
################################################################################

def atoms_fixed(atoms):

    fixed = np.concatenate([a.__dict__['index'] for a in atoms.constraints 
                            if a.__class__.__name__ == 'FixAtoms'])

    return fixed

################################################################################
# ATOMS NOT FIXED
################################################################################

def atoms_not_fixed(atoms):

    fixed = atoms_fixed(atoms)
    not_fixed = [i for i in range(len(atoms)) if i not in fixed]
    
    return not_fixed

################################################################################
# GET VALENCE ELECTRONS
################################################################################

def get_valence_electrons(atoms):

    n_electrons = 0
    for a in atoms:
        n_electrons += SSSP_VALENCE[a.number]
    
    return n_electrons

################################################################################
# PRINT AXSF
################################################################################

def print_axsf(filename, animation, variable_cell = False, parallel = False):

    if parallel is True:
        f = paropen(filename, 'w+')
    else:
        f = open(filename, 'w+')

    print(' ANIMSTEP', len(animation), file = f)
    print(' CRYSTAL', file = f)

    if variable_cell is False:

        cell = animation[0].cell

        print(' PRIMVEC', file = f)
        print("{0:14.8f} {1:14.8f} {2:14.8f}".format(*cell[0]), file = f)
        print("{0:14.8f} {1:14.8f} {2:14.8f}".format(*cell[1]), file = f)
        print("{0:14.8f} {1:14.8f} {2:14.8f}".format(*cell[2]), file = f)

    for i, atoms in enumerate(animation):

        if variable_cell is True:

            cell = atoms.cell
    
            print(' PRIMVEC', i+1, file = f)
            print("{0:14.8f} {1:14.8f} {2:14.8f}".format(*cell[0]), file = f)
            print("{0:14.8f} {1:14.8f} {2:14.8f}".format(*cell[1]), file = f)
            print("{0:14.8f} {1:14.8f} {2:14.8f}".format(*cell[2]), file = f)

        print(' PRIMCOORD', i+1, file = f)
        print(len(atoms), len(get_atom_list(atoms)), file = f)
        for a in atoms:
            print("{0:3s} {1:14.8f} {2:14.8f} {3:14.8f}".format(a.symbol, 
                  a.position[0], a.position[1], a.position[2]), file = f)

    f.close()

################################################################################
# READ AXSF
################################################################################

def read_axsf(filename):

    with open(filename, 'rU') as fileobj:
        lines = fileobj.readlines()

    for line in lines:
        if 'PRIMCOORD' in line:
            key = 'PRIMCOORD'
            break
        elif 'ATOMS' in line:
            key = 'ATOMS'
            break

    if key == 'PRIMCOORD':
        for n, line in enumerate(lines):
            if 'PRIMVEC' in line:
                break
        cell_vectors = np.zeros((3, 3))
        for i, line in enumerate(lines[n+1:n+4]):
            entries = line.split()
            cell_vectors[i][0] = float(entries[0])
            cell_vectors[i][1] = float(entries[1])
            cell_vectors[i][2] = float(entries[2])
        atoms_zero = Atoms(cell = cell_vectors, pbc = (True, True, True))
        increment = 2

    elif key == 'ATOMS':
        atoms_zero = Atoms(pbc = (False, False, False))
        increment = 1

    key = 'PRIMCOORD'
    animation = []
    for n, line in enumerate(lines):
        if key in line:
            atoms = Atoms(cell = cell_vectors, pbc = (True, True, True))
            for line in lines[n+increment:]:
                entr = line.split()
                if entr[0] == key:
                    break
                symbol = entr[0]
                position = (float(entr[1]), float(entr[2]), float(entr[3]))
                atoms += Atom(symbol, position = position)
            animation += [atoms]

    return animation

################################################################################
# READ VIB ENERGIES
################################################################################

def read_vib_energies(filename = 'log.txt', imaginary = False):

    vib_energies = []
    fileobj = open(filename, 'rU')
    lines = fileobj.readlines()
    fileobj.close()

    for i in range(3, len(lines)):
        
        if lines[i][0] == '-':
            break
        
        string = lines[i].split()[1]
        if string[-1] == 'i':
            if imaginary is True:
                vib_energies.append(complex(0., float(string[:-1])*1e-3))
        else:
            vib_energies.append(complex(float(string)*1e-3))

    return vib_energies

################################################################################
# READ MODES AXSF
################################################################################

def write_modes_axsf(vib, kT = kB*300, nimages = 30):

    for index, energy in enumerate(vib.get_energies()):

        if abs(energy) > 1e-5:
        
            animation = []
            mode = vib.get_mode(index) * sqrt(kT / abs(vib.hnu[index]))
            p = vib.atoms.positions.copy()
            index %= 3 * len(vib.indices)
            for x in np.linspace(0, 2 * pi, nimages, endpoint = False):
                vib.atoms.set_positions(p + sin(x) * mode)
                animation += [vib.atoms.copy()]
            vib.atoms.set_positions(p)

            print_axsf('mode_{}.axsf'.format(str(index)), animation)

################################################################################
# GET MOMENT OF INERTIA XYZ
################################################################################

def get_moments_of_inertia_xyz(atoms, center = None):

    if center is None:
        center = atoms.get_center_of_mass()

    positions = atoms.get_positions()-center
    masses = atoms.get_masses()

    I = np.zeros(3)

    for i in range(len(atoms)):

        x, y, z = positions[i]
        m = masses[i]

        I[0] += m*(y**2+z**2)
        I[1] += m*(x**2+z**2)
        I[2] += m*(x**2+y**2)

    return I

################################################################################
# END
################################################################################
