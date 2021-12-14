#!/usr/bin/env python3

################################################################################
# Raffaele Cheula, LCCP, Politecnico di Milano, raffaele.cheula@polimi.it
################################################################################

import os, argparse
from distutils.util import strtobool
from shape.qe_utils import write_pp_input

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

    try: os.mkdir('bader')
    except: pass

    os.chdir('bader')

    pw_out_dir = '../'

else:

    pw_out_dir = './'

################################################################################
# PRINT BADER INPUT
################################################################################

prefix = 'calc'
outdir = 'tmp'

write_input = parsed_args.write_input[0]

if write_input is True:

    pp_data   = {}
    plot_data = {}
    
    pp_data['prefix']          = prefix
    pp_data['outdir']          = pw_out_dir+outdir
    
    plot_data['nfile']         = 1
    plot_data['iflag']         = 3
    plot_data['output_format'] = 6
    
    pp_data['filplot']         = 'filplot_val'
    pp_data['plot_num']        = 0
    plot_data['fileout']       = 'pp_val.cube'

    write_pp_input(pp_data, plot_data, filename = 'pp_val.inp')

    pp_data['filplot']         = 'filplot_all'
    pp_data['plot_num']        = 21
    plot_data['fileout']       = 'pp_all.cube'

    write_pp_input(pp_data, plot_data, filename = 'pp_all.inp')

################################################################################
# RUN PP AND BADER
################################################################################

run_qe_bin = parsed_args.run_qe_bin[0]

if run_qe_bin is True:
    
    os.system('pp.x < pp_val.inp > pp_val.out')
    os.system('pp.x < pp_all.inp > pp_all.out')

    bader_bin = '/galileo/home/userexternal/rcheula0/bader/bader'

    os.system(bader_bin+' pp_val.cube -ref pp_all.cube -vac auto > bader.out')

################################################################################
# CHANGE DIRECTORY
################################################################################

if change_dir is True:

    os.chdir('..')

################################################################################
# END
################################################################################
