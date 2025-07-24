#!/bin/sh 

#for run in {14500..14600}
#do 
run=$1

echo ${run}

JRunAnalyzer -f root://ccxroot:1999//hpss/in2p3.fr/group/km3net/data/raw/sea/KM3NeT_00000049/10/KM3NeT_00000049_000${run}.root -a /sps/km3net/users/widrissi/JRunAnalyzer_plots/KM3NeT_00000049_20220116.detx -o /sps/km3net/users/widrissi/JRunAnalyzer_plots/plots/JRA_KM3NeT_0000049_000${run}.root -@ "p0=1"

#done 
