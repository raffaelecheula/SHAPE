#!/bin/bash

################################################################################
# USER VARIABLES
################################################################################

PROJECTNAME='IscrC_STEFANIE' # IscrC_SHAPING | IscrC_FINENESS
PRODUCTION_PARTITION='m100_usr_prod'
SERIAL_PARTITION='m100_all_serial'
DEBUG_QOS='m100_qos_dbg'

n_proc_per_node=32

run_num_max_pw=5
run_num_max_neb=10
run_num_max_vib=8

bash_scripts_dir=$HOME/BashScripts
python_scripts_dir=$HOME/PythonScripts

status_file='STATUS.txt'

env_file='environ.in'

################################################################################
# READ ARGUMENTS
################################################################################

CONSECUTIVES=0

JOBNAME='automatic'
QOS='automatic'

N_POOL='automatic'
N_POOL_AUTO=4

restart_options=false
pool_options=false
band_options=false
environ_options=false
diag_options=false
copy_py_file=false
autorestart_options=false

for i in "$@"
do
case $i in
    -r|--restart)
    RESTART=true
    shift
    ;;
    -fs|--from_scratch)
    FROM_SCRATCH=true
    shift
    ;;
    -env=*|--environ=*)
    ENVIRON="${i#*=}"
    shift
    ;;
    -n=*|--nodes=*)
    NODES="${i#*=}"
    shift
    ;;
    -np=*|--n_pool=*)
    N_POOL="${i#*=}"
    shift
    ;;
    -nb=*|--n_band=*)
    N_BAND="${i#*=}"
    shift
    ;;
    -c=*|--consecutives=*)
    CONSECUTIVES="${i#*=}"
    shift
    ;;
    -aj=*|--after_job=*)
    AFTER_JOB="${i#*=}"
    shift
    ;;
    -t=*|--time=*)
    TIME="${i#*=}"
    shift
    ;;
    -p=*|--projectname=*)
    PROJECTNAME="${i#*=}"
    shift
    ;;
    -q=*|--queue=*)
    QOS="${i#*=}"
    shift
    ;;
    -j=*|--jobname=*)
    JOBNAME="${i#*=}"
    shift
    ;;
    -i=*|--iteration=*)
    ITERATION="${i#*=}"
    shift
    ;;
    -d=*|--diagonalization=*)
    DIAGONALIZATION="${i#*=}"
    shift
    ;;
    -l|--local)
    LOCAL=true
    shift
    ;;
    --check)
    CHECK=true
    shift
    ;;
    -*)
    WRONG="${i#*}"
    shift
    ;;
    *)
    CALCULATION="${i#*}"
    shift
esac
done

################################################################################
# WRONG ARGUMENTS
################################################################################

if [ $WRONG ]; then
    echo "Wrong argument:" $WRONG
    exit 1
fi
if [ ! $CALCULATION ]; then
    echo "Specify calculation type!"
    exit 1
fi

################################################################################
# CALCULATION
################################################################################

if [ "$CALCULATION" == "pw" ]; then
    cmd_file='pw_cmd'
    inp_file='pw.inp'
    py_file='pw.py'
    restart_options=true
    autorestart_options=true
    run_num_max=$run_num_max_pw
elif [ "$CALCULATION" == "cp" ]; then
    cmd_file='cp_cmd'
    inp_file='cp.inp'
    restart_options=true
elif [ "$CALCULATION" == "ph" ]; then
    cmd_file='ph_cmd'
    inp_file='ph.inp'
    restart_options=true
elif [ "$CALCULATION" == "neb" ]; then
    cmd_file='neb_cmd'
    inp_file='neb.inp'
    py_file='neb.py'
    restart_options=true
    autorestart_options=true
    run_num_max=$run_num_max_neb
    copy_py_file=true
elif [ "$CALCULATION" == "eos" ]; then
    cmd_file='eos_cmd'
    inp_file='eos.py'
elif [ "$CALCULATION" == "vib" ]; then
    cmd_file='vib_cmd'
    inp_file='vib.py'
    py_file='vib.py'
    autorestart_options=true
    run_num_max=$run_num_max_vib
    copy_py_file=true
elif [ "$CALCULATION" == "phon" ]; then
    cmd_file='phon_cmd'
    inp_file='phon.py'
elif [ "$CALCULATION" == "hop" ]; then
    cmd_file='hop_cmd'
    inp_file='hop.py'
elif [ "$CALCULATION" == "md" ]; then
    cmd_file='md_cmd'
    inp_file='md.py'
    restart_options=true
elif [ "$CALCULATION" == "cij" ]; then
    cmd_file='cij_cmd'
    inp_file='cij.py'
elif [ "$CALCULATION" == "projwfc" ]; then
    cmd_file='projwfc_cmd'
    inp_file='projwfc.py'
    py_file='projwfc.py'
    copy_py_file=true
elif [ "$CALCULATION" == "bader" ]; then
    cmd_file='bader_cmd'
    py_file='bader.py'
    copy_py_file=true
elif [ "$CALCULATION" == "dos" ]; then
    cmd_file='dos_cmd'
    py_file='dos.py'
    copy_py_file=true
elif [ "$CALCULATION" == "potential" ]; then
    cmd_file='potential_cmd'
    py_file='potential.py'
    copy_py_file=true
elif [ "$CALCULATION" == "bonds" ]; then
    cmd_file='bonds_cmd'
    py_file='bonds.py'
    copy_py_file=true
else
    echo "Wrong calculation keyword:" $CALCULATION
    exit 1
fi

################################################################################
# RUN LOCAL AND OPTIONS
################################################################################

if [ $LOCAL ]; then
    RUN_MPI=false
    RUN_LOC=true
    cmd_options=false
else
    RUN_MPI=true
    RUN_LOC=false
    cmd_options=true
fi

if [ "$CALCULATION" == "pw" ] || [ "$CALCULATION" == "neb" ]; then
    pool_options=true
    band_options=true
    environ_options=true
    diag_options=true
fi

################################################################################
# CHECK ONLY
################################################################################

if [ $CHECK ]; then
    RUN_MPI=false
    RUN_LOC=false
fi

################################################################################
# CHECK FILES
################################################################################

if [ $copy_py_file = true ]; then
    if [ ! -f $py_file ]; then
        cp $python_scripts_dir/$py_file .
    fi
fi

if [ $cmd_options = true ]; then
    if [ ! -f $inp_file ]; then
        echo "No input file"
        exit 1
    fi
    if [ ! -f $cmd_file ]; then
        cp $bash_scripts_dir/$cmd_file .
    fi
fi

if [ -f $env_file ]; then
    ENVIRON=true
fi

################################################################################
# AUTOMATIC RESTART
################################################################################

if [ $autorestart_options = true ]; then
    touch $status_file
    run_str=$(grep "RUN" $status_file)
    if [[ -n $run_str ]]; then
        run_num=${run_str:4}
        run_num_new=$(echo "$run_num" | awk '{printf "%02.0f\n", $1+1}')
        if [[ $run_num_new -gt $run_num_max ]]; then
            re_str=$(grep 'autorestart=' $cmd_file | head -1)
            re_str_new="autorestart=false"
            sed -i -e "s/$re_str/$re_str_new/" $cmd_file
        fi
    else
        echo 'RUN 00' >> $status_file
    fi
fi

################################################################################
# FROM_SCRATCH
################################################################################

if [ $FROM_SCRATCH ] && [ $restart_options = true ]; then
    if [ "$CALCULATION" == "pw" ] || [ "$CALCULATION" == "neb" ]; then
        fs_str=$(grep 'restart_mode' $inp_file)
        fs_newstr="   restart_mode     = 'from_scratch'"
    elif [ "$CALCULATION" == "ph" ]; then
        fs_str=$(grep 'recover' $inp_file)
        fs_newstr="   recover          = .false."
    elif [ "$CALCULATION" == "md" ]; then
        fs_str=$(grep 'restart = ' $inp_file)
        fs_newstr="restart = False"
    fi
    sed -i -e "s/$fs_str/$fs_newstr/" $inp_file
    echo 'RUN 00' > $status_file
    echo 'from_scratch' >> $status_file
    if [ -f $env_file ]; then
        env_str=$(grep 'environ_restart' $env_file)
        env_newstr="   environ_restart         = .false."
        sed -i -e "s/$env_str/$env_newstr/" $env_file
    fi
fi

################################################################################
# RESTART
################################################################################

if [ $RESTART ] && [ $restart_options = true ]; then
    if [ "$CALCULATION" == "pw" ] || [ "$CALCULATION" == "neb" ]; then
        r_str=$(grep 'restart_mode' $inp_file)
        r_newstr="   restart_mode     = 'restart'"
    elif [ "$CALCULATION" == "ph" ]; then
        r_str=$(grep 'recover' $inp_file)
        r_newstr="   recover          = .true."
    elif [ "$CALCULATION" == "md" ]; then
        r_str=$(grep 'restart = ' $inp_file)
        r_newstr='restart = True'
    fi
    sed -i -e "s/$r_str/$r_newstr/" $inp_file
    echo 'restart' >> $status_file
    if [ -f $env_file ]; then
        env_str=$(grep 'environ_restart' $env_file)
        env_newstr="   environ_restart         = .true."
        sed -i -e "s/$env_str/$env_newstr/" $env_file
    fi
fi

################################################################################
# ENVIRON
################################################################################

if [ $ENVIRON ] && [ $environ_options = true ]; then
    e_str=$(grep 'environ=' $cmd_file)
    e_newstr="environ=true"
    sed -i -e "s/$e_str/$e_newstr/" $cmd_file
fi

################################################################################
# ITERATION
################################################################################

if [ $ITERATION ] && [ "$CALCULATION" == "pw" ]; then
    if [[ $ITERATION == 0 ]]; then
        touch ITERATION_0
        echo 'ITERATION 0' >> $status_file
        if [ ! -f $py_file ]; then
            cp $python_scripts_dir/$py_file .
        fi
        python3 $py_file -i 0
    elif [[ $ITERATION == 1 ]]; then
        rm -f ITERATION_0
        echo 'ITERATION 1' >> $status_file
        python3 $py_file -i 1
    else
        echo 'Wrong iteration keyword:' $ITERATION
        exit 1
    fi
fi

################################################################################
# N_POOL
################################################################################

if [ $N_POOL ] && [ $pool_options = true ]; then
    if [ $N_POOL == 'automatic' ]; then
        gamma=$(grep 'K_POINTS gamma' $inp_file)
        if [[ -n $gamma ]]; then
            N_POOL=1
        else
            N_POOL=$N_POOL_AUTO
        fi
    fi
    np_str=$(grep 'N_POOL=' $cmd_file)
    np_newstr='N_POOL='$N_POOL
    sed -i -e "s/$np_str/$np_newstr/" $cmd_file
fi

################################################################################
# N_BAND
################################################################################

if [ $N_BAND ] && [ $band_options = true ]; then
    nb_str=$(grep 'N_BAND=' $cmd_file)
    nb_newstr='N_BAND='$N_BAND
    sed -i -e "s/$nb_str/$nb_newstr/" $cmd_file
fi

################################################################################
# NODES
################################################################################

if [ $NODES ] && [ $cmd_options = true ]; then
    n_proc=$(echo "$NODES" "$n_proc_per_node" | awk '{printf "%.0f\n", $1*$2}')
    NODES_new=$(echo "$NODES" | awk '{printf "%.0f\n", $1+0.49}')
    n_str=$(grep '#SBATCH --nodes=' $cmd_file)
    n_newstr="#SBATCH --nodes="$NODES_new
    sed -i -e "s/$n_str/$n_newstr/" $cmd_file
fi

################################################################################
# TIME
################################################################################

if [ $TIME ] && [ $cmd_options = true ]; then
    time=$(echo "$TIME" | awk '{printf "%.0f\n", $1*3400-200}')
    if [ "$CALCULATION" == "pw" ]; then
        t_str=$(grep 'max_seconds' $inp_file)
        t_newstr="   max_seconds      = $time"
        sed -i -e "s/$t_str/$t_newstr/" $inp_file
    elif [ "$CALCULATION" == "neb" ]; then
        t_str=$(grep 'max_seconds' $inp_file)
        t_newstr="   max_seconds      = $time"
        sed -i -e "s/$t_str/$t_newstr/" $inp_file
        #if [ -f $py_file ]; then
        #    t_str=$(grep 'max_seconds = ' $py_file)
        #    t_newstr='max_seconds = $time'
        #    sed -i -e "s/$t_str/$t_newstr/" $py_file
        #fi
    elif [ "$CALCULATION" == "vib" ] || [ "$CALCULATION" == "md" ] ||
         [ "$CALCULATION" == "phon" ]; then
        t_str=$(grep 'max_seconds = ' $inp_file)
        t_newstr="max_seconds = $time"
        sed -i -e "s/$t_str/$t_newstr/" $inp_file
    fi
    hours=$(echo "$TIME" | cut -d. -f1 | awk '{printf "%02d\n", $1}')
    minutes=$(echo "$TIME" "$hours" | awk '{printf "%02d\n", ($1-$2)*60}')
    wt_str=$(grep '#SBATCH --time=' $cmd_file)
    wt_newstr='#SBATCH --time='$hours':'$minutes':00'
    sed -i -e "s/$wt_str/$wt_newstr/" $cmd_file
fi

################################################################################
# PROJECTNAME
################################################################################

if [ $PROJECTNAME ] && [ $cmd_options = true ]; then
    p_str=$(grep '#SBATCH --account=' $cmd_file)
    p_newstr='#SBATCH --account='$PROJECTNAME
    sed -i -e "s/$p_str/$p_newstr/" $cmd_file
fi

################################################################################
# JOBNAME
################################################################################

if [ $JOBNAME ] && [ $cmd_options = true ]; then
    j_str=$(grep '#SBATCH --job-name=' $cmd_file)
    if [ $JOBNAME == 'automatic' ]; then
        JOBNAME=$(echo $(basename "$PWD"))
    fi
    j_newstr='#SBATCH --job-name='$JOBNAME
    sed -i -e "s/$j_str/$j_newstr/" $cmd_file
fi

################################################################################
# PARTITION
################################################################################

if [ $QOS ] && [ $cmd_options = true ]; then
    if [ $QOS == 'automatic' ] && [ $NODES ] && [ $TIME ]; then
        pa_str=$(grep '#SBATCH --qos=' $cmd_file)
        if (( $(echo "$NODES <= 4" | bc -l) )) && 
            (( $(echo "$TIME <= 2." | bc -l) )); then
            QOS=$DEBUG_QOS
            pa_newstr='#SBATCH --qos='$QOS
        else
            pa_newstr='##SBATCH --qos='
        fi
        sed -i -e "s/$pa_str/$pa_newstr/" $cmd_file
    fi
fi

################################################################################
# DIAGONALIZATION
################################################################################

if [ $DIAGONALIZATION ] && [ diag_options = true ]; then
    if [[ $DIAGONALIZATION == david ]]; then
        d_addstr="   diago_david_ndim = 2"
    elif [[ $DIAGONALIZATION == cg ]]; then
        d_addstr="   diago_cg_maxiter = 100"
    else
        echo 'Wrong diagonalization'
        exit 1
    fi
    d_str=$(grep 'diagonalization' $inp_file)
    d_newstr="  diagonalization = $DIAGONALIZATION"
    sed -i -e "s/$d_str/$d_newstr/" $inp_file
    d_david=$(grep 'diago_david_ndim' $inp_file)
    d_cg=$(grep 'diago_cg_maxiter' $inp_file)
    if [[ -n $d_david ]]; then
        sed -i -e "s/$d_david/$d_addstr/" $inp_file
    elif [[ -n $d_cg ]]; then
        sed -i -e "s/$d_cg/$d_addstr/" $inp_file
    fi
fi

################################################################################
# RUN JOB
################################################################################

if [ $RUN_MPI = true ]; then

    if [ $autorestart_options = true ]; then
        run_str=$(grep "RUN" $status_file)
        run_num=${run_str:4}
        run_num_new=$(echo "$run_num" | awk '{printf "%02.0f\n", $1+1}')
        run_newstr='RUN '$run_num_new
        sed -i -e "s/$run_str/$run_newstr/" $status_file
    fi

    if [ $AFTER_JOB ]; then
        sub=$(sbatch --dependency=afterany:$AFTER_JOB $cmd_file)
        echo $sub
    else
        sub=$(sbatch $cmd_file)
        echo $sub
    fi
    
    i=0
    echo 'running job:' ${sub:20} $CALCULATION >> $status_file
    
    while [ $i -lt $CONSECUTIVES ]
    do
        sub=$(sbatch --dependency=afterany:${sub:20} $cmd_file)
        echo $sub
        echo 'running job:' ${sub:20} $CALCULATION >> $status_file
        i=$[i+1]
    done

elif [ $RUN_LOC = true ]; then

    module load autoload profile/global qe

    python3 $py_file

fi

exit 0

################################################################################
# END
################################################################################

