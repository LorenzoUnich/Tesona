#!/bin/sh

export CLING_DEBUG=1

DIR_out=/sps/km3net/users/lunich/binned/arca-ps-aart_update_bands/scripts/pickle_files


 sbatch -L sps --time=0-24:00 --mem-per-cpu 8G -n 10 -o ${DIR_out}/o_${DATETIME}_100000_PE_with_multiprocessing_aysn_different_chunksize_four_less_GIGA_200000PE_SECOND_2PI.txt -e ${DIR_out}/e_${DATETIME}.txt  data_sample_133.py
#quindi quello che fai e che gli chiedi di dividere il lavor oin tanti nodiiiiiiii gnammete



