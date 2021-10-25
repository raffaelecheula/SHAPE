#!/usr/bin/env python3

################################################################################
# Raffaele Cheula, LCCP, Politecnico di Milano, raffaele.cheula@polimi.it
################################################################################

from __future__ import division, print_function
import signal
import numpy as np
from ase.io import read
from ase.units import kB
from ase.vibrations import Vibrations
from ase.calculators.espresso import Espresso
from qe_utils import ReadQeInp, ReadQeOut, write_modes_axsf
from supercell_utils import atoms_not_fixed

################################################################################
# MAX SECONDS
################################################################################

max_seconds = 55050
signal.alarm(max_seconds)

################################################################################
# VIBRATIONS
################################################################################

atoms_to_vib = 'all' # 'all' | 'only_relaxed' | n last atoms | list of indices

qe_inp = ReadQeInp('pw.inp')

input_data, pseudos, kpts, koffset = qe_inp.get_data_pseudos_kpts()

constraints = qe_inp.get_atoms().constraints

if input_data['calculation'] == 'scf':
    qe_out = ReadQeOut('pw.inp')
else:
    qe_out = ReadQeOut('pw.out')

atoms = qe_out.get_atoms()

atoms.constraints = constraints

conv_thr = 1e-07

calc = Espresso(input_data       = input_data    ,
                pseudopotentials = pseudos       ,
                kpts             = kpts          ,
                koffset          = koffset       ,
                tprnfor          = True          ,
                tstress          = True          ,
                calculation      = 'scf'         ,
                restart_mode     = 'from_scratch',
                max_seconds      = max_seconds   ,
                conv_thr         = conv_thr      )

atoms.set_calculator(calc)

if atoms_to_vib == 'all':
    indices = [a.index for a in atoms]

elif atoms_to_vib == 'only_relaxed':
    indices = atoms_not_fixed(atoms)

elif type(atoms_to_vib) == int:
    indices = [a.index for a in atoms][len(atoms)-atoms_to_vib:]

else:
    indices = atoms_to_vib

atoms.set_constraint()

print('indices of atoms to vibrate = {}'.format(indices))

vib = Vibrations(atoms, indices = indices, delta = 0.01, nfree = 2)

################################################################################
# RUN VIB
################################################################################

run_vib = True

if run_vib is True:

    vib.clean(empty_files = True)
    vib.run()
    open('log.txt', 'w').close()
    vib.summary(method = 'standard', log = 'log.txt')

    write_modes_axsf(vib, kT = kB * 300, nimages = 30)

################################################################################
# WRITE HESSIAN
################################################################################

write_hessian = False

if write_hessian is True:

    fileobj = open('Hessian.txt', 'w+')

    for line in vib.H:
        for num in line:
            print('{:7.2f}'.format(num), end = '', file = fileobj)
        print('', file = fileobj)

    fileobj.close()

################################################################################
# END
################################################################################
