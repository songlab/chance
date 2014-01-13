#!/bin/bash

#PBS -e chance_test.err
#PBS -N chance
#PBS -m ae
#PBS -l nodes=1:ppn=12:bigmem


/songlab/aaron/research/fun_genom/active_work/chance/chance_com/chance/distrib/run_chance.sh /songlab/aaron/MATLAB/MATLAB_Compiler_Runtime/v80/ batch -p /songlab/aaron/research/fun_genom/active_work/chance/chance_com/example_parameter_files/param_batch.txt -o /songlab/aaron/research/fun_genom/active_work/chance/chance_com/chance/distrib/chance_out.txt