#!/usr/bin/env python3

################################################################################
# Raffaele Cheula, LCCP, Politecnico di Milano, raffaele.cheula@polimi.it
################################################################################

from __future__ import absolute_import, division, print_function
import os, argparse
from math import ceil
import numpy as np
from ase.neb import NEB
from ase.calculators.espresso import Espresso
from shape.qe_utils import (write_neb_inp, print_axsf, read_qe_inp, read_qe_out,
                            read_neb_crd, reorder_neb_images)

################################################################################
# MAX SECONDS
################################################################################

max_seconds = 14000

iteration = 0

################################################################################
# PARSE CMD LINE ARGUMENTS
################################################################################

parser = argparse.ArgumentParser()

parser.add_argument('--iteration', '-i',
                    type     = int,
                    nargs    = '+',
                    required = False,
                    help     = 'iteration of simulation')

args = parser.parse_args()

################################################################################
# ITERATIONS
################################################################################

try: iteration = args.iteration[0]
except: pass

if iteration == 0:
    nstep_path       = 20
    opt_scheme       = 'quick-min'
    ci_scheme        = 'no-CI'
    k_max            = 0.4
    k_min            = 0.2
    path_thr         = 0.5
    use_masses       = True
    ecutwfc          = 25.0
    ecutrho          = 200.0
    conv_thr         = 1e-05
    nstep            = 300
    gamma            = True
    restart_from_crd = False

elif iteration == 1:
    nstep_path       = 20
    opt_scheme       = 'quick-min'
    ci_scheme        = 'no-CI'
    k_max            = 0.4
    k_min            = 0.2
    path_thr         = 0.3
    use_masses       = True
    ecutwfc          = None
    ecutrho          = None
    conv_thr         = None
    nstep            = 300
    gamma            = False
    restart_from_crd = True

elif iteration == 2:
    nstep_path       = 20
    opt_scheme       = 'broyden'
    ci_scheme        = 'auto'
    k_max            = 0.8
    k_min            = 0.4
    path_thr         = 0.05
    use_masses       = False
    ecutwfc          = None
    ecutrho          = None
    conv_thr         = None
    nstep            = 300
    gamma            = False
    restart_from_crd = True

elif iteration == 3:
    nstep_path       = 30
    opt_scheme       = 'broyden'
    ci_scheme        = 'auto'
    k_max            = 1.0
    k_min            = 0.6
    path_thr         = 0.05
    use_masses       = False
    ecutwfc          = None
    ecutrho          = None
    conv_thr         = None
    nstep            = 3000
    gamma            = False
    restart_from_crd = True

################################################################################
# NUDGED ELASTIC BAND
################################################################################

n_images = 10

input_data, pseudos, kpts, koffset = read_qe_inp('first/pw.inp')

del input_data['calculation']
del input_data['restart_mode']

input_data['max_seconds'] = max_seconds
input_data['nstep']       = nstep
#input_data['wf_collect']  = False

if ecutwfc is not None:
    input_data['ecutwfc'] = ecutwfc

if ecutrho is not None:
    input_data['ecutrho'] = ecutrho

if conv_thr is not None:
    input_data['conv_thr'] = conv_thr

if gamma is True:
    kpts    = (1, 1, 1)
    koffset = (0, 0, 0)

calc = Espresso(input_data       = input_data,
                pseudopotentials = pseudos   ,
                kpts             = kpts      ,
                koffset          = koffset   )

first = read_qe_out('first/pw.out')
last = read_qe_out('last/pw.out')
images = [first]
images += [ first.copy() for i in range(n_images-2) ]
images += [last]

if restart_from_crd is True:
    images = read_neb_crd(images, 'pwscf.crd')
else:
    neb = NEB(images)
    neb.interpolate('idpp')
    print_axsf('pwscf.axsf', images)

neb_data = {}

neb_data['string_method'] = 'neb'
neb_data['restart_mode']  = 'from_scratch'
neb_data['nstep_path']    = nstep_path
neb_data['num_of_images'] = n_images
neb_data['opt_scheme']    = opt_scheme
neb_data['CI_scheme']     = ci_scheme
neb_data['ds']            = 1
neb_data['k_max']         = k_max
neb_data['k_min']         = k_min
neb_data['path_thr']      = path_thr
neb_data['use_masses']    = use_masses

write_neb_inp(neb_data, images, calc)

################################################################################
# RUN QUANTUM ESPRESSO
################################################################################

run_qe = False

if run_qe is True:
    os.system("run neb -n=4 -t=4")

################################################################################
# END
################################################################################
