#!/bin/bash

xCLING_DEBUG=1

folder=$1
foldername=$2
bands=$3
python arca_updated_binned.py -c ${folder} -p ${foldername} -q True -f ARCA_21_all_period_bdt095 -A ANTARES_2022 -W Expr -D False -L True -B ${bands}  # -p e mettici path outuput lo s>

#questo script lancia effettivamente arca binned
#python arca_binned.py -c ${folder} -q Tue -W Expr -D True -M 6 -m 7 
#L è lo stcking (lascia sta), Q è il primo periodo. per le prove di ora va bene. -F e -A lascia sta.  -W sarebbe la law di potenza. -C è la folder, e non è altro che la sorgente i-ma.>
#

