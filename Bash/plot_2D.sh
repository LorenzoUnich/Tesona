#!/bin/bash 

for run in {14..15000}
do
    for du in {1..32}
    do	
        for dom in {1..18}
	do  
#	    JPlot2D -B -Z -f /sps/km3net/users/widrissi/JRunAnalyzer_plots/plots/JRA_KM3NeT_00000133_000${run}.root:Detector/DU${du}/F${dom}/h_pmt_rate_distributions_Timeslice -O "colz" -o /sps/km3net/users/widrissi/JRunAnalyzer_plots/plots_2d_rate/L0_RATE_run${run}_du${du}_dom${dom}.png
#	    dom_id = ${dom}
#	    du_id = ${du}
	    module=$(JPrintDetector -a /sps/km3net/users/widrissi/JRunAnalyzer_plots/KM3NeT_00000133_20220917.detx -O modules | awk -F' ' '{if ($3=="'$du'" && $4=="'$dom'")print $2}')
	    echo ${module}
	    if [ -z "${module//[0-9]}" ]
	    then
		  JPlot2D -B -Z -f /sps/km3net/users/widrissi/JRunAnalyzer_plots/plots/JRA_KM3NeT_00000133_000${run}.root:Detector/DU${du}/F${dom}/${module}_Timeslice_2SToT -O "colz" -o /sps/km3net/users/widrissi/JRunAnalyzer_plots/plot_2d_tot/L0_TOT_run${run}_du${du}_dom${dom}_module${module}.png
            else 
		  continue 		
            fi
	
	done
    done

done
