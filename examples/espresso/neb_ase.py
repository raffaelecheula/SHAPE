#!/usr/bin/env python3

################################################################################
# Raffaele Cheula, LCCP, Politecnico di Milano, raffaele.cheula@polimi.it
################################################################################

import signal
#from mpi4py import MPI
from ase.io import read
from ase.optimize import BFGS, FIRE, MDMin
from ase.neb import NEB, NEBTools
from ase.io.trajectory import TrajectoryWriter
from ase.calculators.espresso import Espresso
from ase.parallel import world, parprint, paropen
from shape.qe_utils import read_qe_inp, read_qe_out, print_axsf

################################################################################
# MAX SECONDS
################################################################################

iteration = 0

max_seconds = 14000
signal.alarm(max_seconds)

#rank = MPI.COMM_WORLD.Get_rank()

################################################################################
# PARAMETERS
################################################################################

run_neb_ase  = True
run_neb_ml   = False
run_neb_auto = False

first = read_qe_out('first/pw.out')
last  = read_qe_out('last/pw.out')

################################################################################
# ITERATIONS
################################################################################

if iteration == 0:
    climb    = False
    fmax     = 0.100 # [eV/A]
    k_spring = 0.100 # [eV/A]
    ecutwfc  =  25.0
    ecutrho  = 200.0
    conv_thr = 1e-05
    gamma    = True
    restart  = False

elif iteration == 1:
    climb    = False
    fmax     = 0.100 # [eV/A]
    k_spring = 0.100 # [eV/A]
    ecutwfc  = None
    ecutrho  = None
    conv_thr = 1e-05
    gamma    = True
    restart  = True

elif iteration == 2:
    climb    = False
    fmax     = 0.050 # [eV/A]
    k_spring = 0.100 # [eV/A]
    ecutwfc  = None
    ecutrho  = None
    conv_thr = 1e-06
    gamma    = False
    restart  = True

elif iteration == 3:
    climb    = True
    fmax     = 0.050 # [eV/A]
    k_spring = 0.100 # [eV/A]
    ecutwfc  = None
    ecutrho  = None
    conv_thr = 1e-06
    gamma    = False
    restart  = True

################################################################################
# QUANTUM ESPRESSO
################################################################################

pw_data, pseudos, kpts, koffset = read_qe_inp('first/pw.inp')

pw_data['calculation']  = 'scf'
pw_data['restart_mode'] = 'from_scratch'
pw_data['max_seconds']  = max_seconds
pw_data['tprnfor']      = True
pw_data['tstress']      = True

if ecutwfc is not None:
    pw_data['ecutwfc'] = ecutwfc

if ecutrho is not None:
    pw_data['ecutrho'] = ecutrho

if conv_thr is not None:
    pw_data['conv_thr'] = conv_thr

if gamma is True:
    kpts    = None
    koffset = None

################################################################################
# NEB ASE
################################################################################

if run_neb_ase is True:

    parprint('\n NEB ASE \n')
    
    n_images  = 10
    
    method    = 'eb'   # aseneb | improvedtangent | eb
    optimizer = 'BFGS' # BFGS | MDMin | FIRE
    
    trajname  = 'neb_ase.traj'
    
    only_initialization = True
    
    if restart is False:
        images = [first]
        images += [first.copy() for i in range(n_images-2)]
        images += [last]
        neb = NEB(images)
        neb.interpolate('idpp')
    else:
        images = read(trajname)
    
    neb = NEB(images   = images  ,
              k        = k_spring,
              climb    = climb   ,
              parallel = True    ,
              world    = world   ,
              method   = method  )
    
    if optimizer == 'BFGS':
        opt = BFGS(neb, trajectory = trajname)
    
    elif optimizer == 'MDMin':
        opt = MDMin(neb, trajectory = trajname)
    
    elif optimizer == 'FIRE':
        opt = FIRE(neb, trajectory = trajname)
    
    for image in images:
        image.calc = Espresso(input_data       = pw_data,
                              pseudopotentials = pseudos,
                              kpts             = kpts   ,
                              koffset          = koffset)
    
    if only_initialization is False:
    
        images[0].get_potential_energy()
        images[-1].get_potential_energy()
        
        opt.run(fmax = fmax)
        
        nebtools = NEBTools(images)
        
        parprint(nebtools.get_barrier())
    
    print_axsf('pwscf.axsf', images, parallel = True)

################################################################################
# NEB ML
################################################################################

if run_neb_ml is True:

    parprint('\n NEB ML \n')
    
    from catlearn.optimize.mlneb import MLNEB
    from catlearn.optimize.tools import plotneb
    
    first.write('first.traj')
    last.write('last.traj')
    
    n_images = 10
    
    trajname = 'neb_ml.traj'
    
    ase_calc = Espresso(input_data       = pw_data,
                        pseudopotentials = pseudos,
                        kpts             = kpts   ,
                        koffset          = koffset)
    
    neb_catlearn = MLNEB(start         = 'first.traj',
                         end           = 'last.traj' ,
                         ase_calc      = ase_calc    ,
                         n_images      = n_images    ,
                         interpolation = 'idpp'      ,
                         restart       = restart     )
    
    neb_catlearn.run(fmax = fmax, trajectory = trajname)

################################################################################
# NEB AUTO
################################################################################

if run_neb_auto is True:

    parprint('\n NEB AUTO \n')
    
    from ase.autoneb import AutoNEB
    
    iter_folder   = 'AutoNEB'
    prefix        = 'image'
    n_simul       =  1
    n_images      = 10
    n_images_init =  4
    smooth_curve  = False
    
    trajname      = 'neb_auto.traj'
    
    if restart is False:
        images = [first]
        images += [first.copy() for i in range(n_images_init-2)]
        images += [last]
        if n_images_init > 2:
            neb = NEB(images)
            neb.interpolate('idpp')
    else:
        images = read(trajname)
    
    for i in range(len(images)):
        
        images[i].calc = Espresso(input_data       = pw_data,
                                  pseudopotentials = pseudos,
                                  kpts             = kpts   ,
                                  koffset          = koffset)
        
        images[i].get_potential_energy()
        images[i].write('{0}{1:03d}.traj'.format(prefix, i))
    
    def attach_calculators(images):
        for image in images:
            image.calc = Espresso(input_data       = pw_data,
                                  pseudopotentials = pseudos,
                                  kpts             = kpts   ,
                                  koffset          = koffset)
    
    neb_auto = AutoNEB(attach_calculators = attach_calculators,
                       prefix             = prefix            ,
                       n_simul            = n_simul           ,
                       n_max              = n_images          ,
                       iter_folder        = iter_folder       ,
                       fmax               = fmax              ,
                       maxsteps           = 10000             ,
                       k                  = 0.1               ,
                       method             = 'eb'              ,
                       optimizer          = 'BFGS'            ,
                       climb              = False             ,
                       space_energy_ratio = 0.5               ,
                       world              = world             ,
                       parallel           = True              ,
                       smooth_curve       = smooth_curve      ,
                       interpolate_method ='idpp'             )
    
    neb_auto.run()
    
    images = neb_auto.all_images
    
    for i in range(len(images)):
        images[i].write('{0}{1:03d}.traj'.format(prefix, i))
    
    with TrajectoryWriter(trajname, mode = 'w') as traj:
        for image in images:
            traj.write(image)

    nebtools = NEBTools(images)

    parprint(nebtools.get_barrier())

    print_axsf('pwscf.axsf', images, parallel = True)

################################################################################
# END
################################################################################
