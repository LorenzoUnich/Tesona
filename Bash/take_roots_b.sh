#script to take the root files of the analysis without any band 

remote_base="ulorenzo@cca.in2p3.fr:/sps/km3net/users/lunich/binned/arca-ps-aart_update_bands/scripts/new_fit_analysis"
folders=(
  "1006"
  "0602"
  "0202"
  "0206"
  )
files=(
  "Source_cand_0_first_declination_histograms_without_2pi.root"
  "Source_cand_1_second_declination_histograms_without_2pi.root"
  "Source_cand_2_third_declination_histograms_without_2pi.root"
  "Source_cand_3_four_declination_histograms_without_2pi.root"
  "Source_cand_4_fifth_declination_histograms_without_2pi.root"
  "Source_cand_5_six_declination_histograms_without_2pi.root"
  "Source_cand_6_seven_declination_histograms_without_2pi.root"
  "Source_cand_7_eight_declination_histograms_without_2pi.root"
)

for file in "${files[@]}"; do

  for fol in "${folders[@]}"; do
  	scp "$remote_base/$fol/$file" .
  done
done
