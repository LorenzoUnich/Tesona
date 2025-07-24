#!/bin/bash

DIR_out=/sps/km3net/users/widrissi/JRunAnalyzer_plots/output_job

for run in {10000..10100}
do 

sbatch -L sps --time=0-24:00 --mem-per-cpu 4G -n 8 -o ${DIR_out}/o_${DATETIME}_${run}.txt -e ${DIR_out}/e_${DATETIME}_${run}.txt jra_ARCA21.sh ${run}

done
