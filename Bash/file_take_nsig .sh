#!/bin/bash

for i in $(seq 0 31); do
  pval=$(grep "TOTAL background:" o__${i}_100000_PE_with_multiprocessing_aysn_different_chunksize_four_less_GIGA_200000PE_SECOND_2PI.txt | awk -F':' '{print $2}')
  echo "$i $pval"
done > nback_total.txt
