#!/usr/bin/env python3

################################################################################
# Raffaele Cheula, LCCP, Politecnico di Milano, raffaele.cheula@polimi.it
################################################################################

import os, argparse
from ase.calculators.espresso import Espresso
from shape.qe_utils import ReadQeInp, ReadQeOut

################################################################################
# PARSE CMD LINE ARGUMENTS
################################################################################

parser = argparse.ArgumentParser()

parser.add_argument('--iteration', '-i',
                    type     = int ,
                    nargs    = 1   ,
                    required = True,
                    help     = 'iteration of simulation')

parsed_args = parser.parse_args()

################################################################################
# ITERATIONS
################################################################################

iteration = parsed_args.iteration[0]

if iteration == 0:

    ecutwfc       = 40.0
    ecutrho       = 320.0
    etot_conv_thr = 0.001
    forc_conv_thr = 0.01
    conv_thr      = 1e-05
    gamma         = True
    restart       = False

if iteration == 1:

    ecutwfc       = None
    ecutrho       = None
    etot_conv_thr = None
    forc_conv_thr = None
    conv_thr      = None
    gamma         = False
    restart       = False

################################################################################
# RELAX
################################################################################

if os.path.exists('./pw.inp.bckp') is False:
    os.rename('pw.inp', 'pw.inp.bckp')

qe_inp = ReadQeInp('pw.inp.bckp')

pw_data, pseudos, kpts, koffset = qe_inp.get_data_pseudos_kpts()

if ecutwfc is not None:
    pw_data['ecutwfc'] = ecutwfc

if ecutrho is not None:
    pw_data['ecutrho'] = ecutrho

if etot_conv_thr is not None:
    pw_data['etot_conv_thr'] = etot_conv_thr

if forc_conv_thr is not None:
    pw_data['forc_conv_thr'] = forc_conv_thr

if conv_thr is not None:
    pw_data['conv_thr'] = conv_thr

if gamma is True:
    kpts    = (1, 1, 1)
    koffset = (0, 0, 0)

if restart is True:
    pw_data['restart_mode'] = 'restart'
else:
    pw_data['restart_mode'] = 'from_scratch'

constraints = qe_inp.get_atoms().constraints

if iteration == 0:
    atoms = qe_inp.get_atoms()
else:
    qe_out = ReadQeOut('pw.out')
    atoms = qe_out.get_atoms()

atoms.constraints = constraints

calc = Espresso(input_data       = pw_data,
                pseudopotentials = pseudos,
                kpts             = kpts   ,
                koffset          = koffset)

atoms.set_calculator(calc)

calc.write_input(atoms)
os.rename('espresso.pwi', 'pw.inp')

################################################################################
# END
################################################################################
