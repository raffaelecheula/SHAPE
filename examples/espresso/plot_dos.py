#!/usr/bin/env python3

################################################################################
# Raffaele Cheula, LCCP, Politecnico di Milano, raffaele.cheula@polimi.it
################################################################################

from __future__ import absolute_import, division, print_function
import os
import numpy as np
import matplotlib.pyplot as plt
from collections import OrderedDict
from shape.qe_utils import ReadQeOut, write_dos_input, get_dos

################################################################################
# PRINT DOS INPUT
################################################################################

print_dos_inp = True

dos_data = OrderedDict()

dos_data['prefix']  = 'calc'
dos_data['outdir']  = 'tmp'
dos_data['ngauss']  = 0
dos_data['degauss'] = 0.001
dos_data['Emin']    = -30.
dos_data['Emax']    = +20.
dos_data['DeltaE']  = 0.01

if print_dos_inp is True:

    write_dos_input(dos_data = dos_data, filename = 'dos.inp')

################################################################################
# RUN DOS
################################################################################

run_dos = True

if run_dos is True:

    os.system('dos.x < dos.inp > dos.out')

################################################################################
# READ DOS
################################################################################

plot_dos = True

if plot_dos is True:

    filename = '{0}.dos'.format(dos_data['prefix'])
    
    energy, dos = get_dos(filename = filename)
    
    fig = plt.figure(1)
    
    x_max = max(dos)+1.
    
    plt.plot(dos, energy)
    plt.plot([0., x_max], [0.]*2, color = 'red')
    
    plt.xlim([0., x_max])
    plt.ylim([dos_data['Emin'], dos_data['Emax']])
    
    plt.xlabel('dos [states/eV]')
    plt.ylabel('energy [eV]')
    
    plt.savefig('dos.png', dpi = 300)

################################################################################
# END
################################################################################
