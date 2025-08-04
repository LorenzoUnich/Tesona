  GNU nano 5.6.1                                                                         run_cands.sh                                                                         Modified  
#!/bin/bash

#BE CAREFUL WIHTH FOLDER AND THE PATH

export CLING_DEBUG=1
band="$1"
echo "hope you've checked the path. It's not my fault you forget things..."
echo "if you  type - scancel -u $USER - it deletes all your queing jobs"
echo band
DIR_out=/sps/km3net/users/lunich/binned/arca-ps-aart_update_bands/scripts/new_catalouge_analysis/$band
#for source_num in {0..8}   ; do
    #for folder_num in 5; do

#  be careful with the setting of the minimum and maximum index of the catalogue
#by this i mean, remember to do so and remember that here indices start from 0. 
if [[ "$band" == "all" ]]; then
  min=0
  max=297
elif [[ "$band" == "1006" ]]; then
  min=0
  max=2
elif [[ "$band" == "0602" ]]; then
  min=3
  max=4
elif [[ "$band" == "0202" ]]; then
  min=5
  max=5 
elif [[ "$band" == "0206" ]]; then
  min=6
  max=7
elif [[ "$band" == "4" ]]; then
 min=4
 max=4
 band="0602"
elif [[ "$band" == "auto" ]]; then
 min=0
 max=297
 echo "I HOPE YOU REMEMBER TO PUT THE CORRECT SOURCE NUMBER"
else
  echo "Valore band non riconosciuto: $band"
  exit 1
fi

echo $band
#here be careful with the  parameters for parallel computing
for (( source_num=min; source_num<=max; source_num++ )); do
    sbatch --partition=htc_highmem  --time=0-24:00 --mem-per-cpu 16G -n 8 -o ${DIR_out}/o_${DATETIME}_${source_num}_100000_PE_with_multiprocessing_aysn_different_chunksize_four_less_G>
done


#sbatch --partition=htc_highmem  --time=0-24:00 --mem-per-cpu 8G -n 10 -o ${DIR_out}/o_${DATETIME}_${source_num}_100000_PE_with_multiprocessing_aysn_different_chunksize_four_less_GIGA>
#quindi quello che fai e che gli chiedi di dividere il lavor oin tanti nodiiiiiiii gnammete
#done


