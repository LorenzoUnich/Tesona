#!/bin/bash

# Lista di cartelle dove cercare
rm -f extracted_summary_bands.txt  # elimina il file se esiste
# Intestazione della tabella
echo -e "i\tss90\tdiscovery\ttotal_signal\ttotal_noise\tfolder" > extracted_summary_bands.txt

# Loop su ciascuna cartella
for file in o__*_100000_PE_with_multiprocessing_aysn_different_chunksize_four_less_GIGA_200000PE_SECOND_2PI.txt; do
    # Controlla che il file esista (evita errori se non matcha nulla)
    [[ -e "$file" ]] || continue

    filename=$(basename "$file")
    i=$(echo "$filename" | grep -oP '(?<=o__)\d+')

    ss90=$(grep -m1 "ss90:" "$file" | awk '{print $2}')
    discovery=$(grep -m1 "discovery:" "$file" | awk '{print $2}')
    total_signal=$(grep -m1 "TOTAL signal:" "$file" | awk '{print $3}')
    total_noise=$(grep -m1 "TOTAL background:" "$file" | awk '{print $3}')
    echo -e "$i\t$ss90\t$discovery\t$total_signal\t$total_noise\t$dir"
done >> extracted_summary_bands.txt
echo " extracted_summary_bands.txt file created"
