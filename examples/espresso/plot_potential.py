#!/usr/bin/env python3

################################################################################
# Raffaele Cheula, LCCP, Politecnico di Milano, raffaele.cheula@polimi.it
################################################################################

import os
import matplotlib.pyplot as plt
from collections import OrderedDict
from shape.qe_utils import write_pp_input

################################################################################
# PRINT PP INPUT
################################################################################

print_pp_inp = True

pp_data   = OrderedDict()
plot_data = OrderedDict()

pp_data['prefix']   = 'calc'
pp_data['outdir']   = 'tmp'
pp_data['filplot']  = 'filplot'
pp_data['plot_num'] = 11

if print_pp_inp is True:

    write_pp_input(pp_data, plot_data, filename = 'pp.inp')

################################################################################
# PRINT AVE INPUT
################################################################################

print_ave_inp = True

n_files  = 1
filplots = [pp_data['filplot']]
weigths  = [1.]
n_points = 1000
plane    = 3
window   = 1.

if print_ave_inp is True:

    fileobj = open('ave.inp', 'w+')
    
    print(n_files, file = fileobj)
    for i in range(n_files):
        print(filplots[i], file = fileobj)
        print(weigths[i], file = fileobj)
    print(n_points, file = fileobj)
    print(plane, file = fileobj)
    print(window, file = fileobj)
    
    fileobj.close()

################################################################################
# RUN PP
################################################################################

run_pp = True

if run_pp is True:

    os.system('pp.x < pp.inp > pp.out')
    
    os.system('average.x < ave.inp > ave.out')

################################################################################
# PLOT POTENTIAL
################################################################################

plot_potential = True

if plot_potential is True:

    fileobj = open('avg.dat', 'r')
    
    lines = fileobj.readlines()
    
    fileobj.close()
    
    x_vec = []
    y_tot = []
    y_ave = []
    
    for line in lines:
        
        line_split = line.split()
        
        x_vec += [float(line_split[0])]
        y_tot += [float(line_split[1])]
        y_ave += [float(line_split[2])]
    
    plt.xlim([min(x_vec), max(x_vec)])
    
    plt.xlabel('z [A]')
    plt.ylabel('V [eV]')
    
    plt.plot(x_vec, y_tot)
    
    plt.savefig('potential.png', dpi = 300)

################################################################################
# END
################################################################################
