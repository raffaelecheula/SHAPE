#!/usr/bin/env python3

################################################################################
# Raffaele Cheula, LCCP, Politecnico di Milano, raffaele.cheula@polimi.it
################################################################################

import os, argparse
import matplotlib.pyplot as plt
from distutils.util import strtobool
from collections import OrderedDict
from shape.qe_utils import write_dos_input, get_dos

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

    try: os.mkdir('dos')
    except: pass

    os.chdir('dos')

    pw_out_dir = '../'

else:

    pw_out_dir = './'

################################################################################
# PRINT DOS INPUT
################################################################################

prefix = 'calc'
outdir = 'tmp'
Emin   = -30.
Emax   = +20.
DeltaE = 0.01

write_input = parsed_args.write_input[0]

if write_input is True:

    dos_data = OrderedDict()

    dos_data['prefix']  = prefix
    dos_data['outdir']  = pw_out_dir+outdir
    dos_data['ngauss']  = 0
    dos_data['degauss'] = 0.001
    dos_data['Emin']    = Emin
    dos_data['Emax']    = Emax
    dos_data['DeltaE']  = DeltaE

    write_dos_input(dos_data = dos_data, filename = 'dos.inp')

################################################################################
# RUN DOS
################################################################################

run_qe_bin = parsed_args.run_qe_bin[0]

if run_qe_bin is True:

    os.system('dos.x < dos.inp > dos.out')

################################################################################
# PLOT DOS
################################################################################

postprocess = parsed_args.postprocess[0]

if postprocess is True:

    energy, dos = get_dos(filename = '{}.dos'.format(prefix))
    
    fig = plt.figure(1)
    
    x_max = max(dos)+1.
    
    plt.plot(dos, energy)
    plt.plot([0., x_max], [0.]*2, color = 'red')
    
    plt.xlim([0., x_max])
    plt.ylim([Emin, Emax])
    
    plt.xlabel('dos [states/eV]')
    plt.ylabel('energy [eV]')
    
    plt.savefig('dos.png', dpi = 300)

################################################################################
# CHANGE DIRECTORY
################################################################################

if change_dir is True:

    os.chdir('..')

################################################################################
# END
################################################################################
