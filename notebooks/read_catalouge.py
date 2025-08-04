import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import re




# Header, li mantengo qua per memoria futura
IR_header = ['Index', 'CommonName', 'IRASName', 'alpha', 'delta', 'S12',
       'S25', 'S60', 'S100', 'DC", "log', 'goals']
header_radio = ["Index", "Name", "S", "S_P", "S_I", "phi_M", "phi_m", "PA", 
          "alpha", "delta", "theta_M", "theta_m", "BPA", "RUN"]




def decl(datum):
    piece = str(datum)[1:]
    if piece == "an":
        return np.nan
    else:
        decl = float(piece[0:2])+float(piece[3:5])/60+   float(piece[6:10])/(3600)

        return  decl

def asc(datum):
    piece = str(datum)
    if piece == "nan":
        return np.nan
    else:
        asc = float(piece[0:2])+float(piece[3:5])/60+   float(piece[6:10])/(3600)
        return asc
def sum_smart(l):
    result = 0
    for i in l:
        if not np.isnan(i):
            result += i
    return result

def extract_fields(line, positions):
        fields = []
        for i in range(len(positions) - 1):
            start = positions[i]
            end = positions[i + 1]
            field = line[start:end].strip()
            if field == '':
                fields.append(np.nan)
            else:
                # prova a convertire in float se possibile
                try:
                    fields.append(float(field))
                except ValueError:
                    fields.append(field)
        return fields
def make_weights(array, filename):
    l = list(array)
    s = sum_smarts(l)
    with open(filename, "w+") as file:
        for w_i in l:
            file.write(str(w_i/s)+"\n")
    


#la funzione restituisce un dataframe, dopo che gli hai dato il nome del file
#e mette come header le stringhe che sono nella lista
def make_df_from_file(filename, header_list, start_idx = 30, convert_coordinates= True, weights = False):

    # Trova la prima riga utile con tutti i campi (es. riga 30)
    with open(filename, "r") as f:
        raw_lines = f.readlines()
    reference_line = raw_lines[start_idx].rstrip("\n")

    # Trova posizioni di inizio dei campi usando regex (più di 1 spazio)
    matches = list(re.finditer(r' {2,}', reference_line))
    start_positions = [0] + [m.end() for m in matches]
    # Aggiungi anche la fine della riga per poter fare gli slice finali
    start_positions.append(len(reference_line) + 1)

    # Funzione per estrarre i campi da una riga data
    

    # Applica a tutte le righe a partire da start_idx
    parsed_lines = []
    for line in raw_lines[start_idx:]:
        line = line.rstrip("\n")
        fields = extract_fields(line, start_positions)
        parsed_lines.append(fields)

    # Crea il DataFrame finale
    dframe = pd.DataFrame(parsed_lines, columns=header_list)
    if convert_coordinates:
        dframe["alpha"] = dframe["alpha"].apply(asc)
        dframe["delta"] = dframe["delta"].apply(decl)
    if weights:
        #print("It is possible to generate a txt file with weights")
        #print("These weights can be calculated directly trough one of the columns")
        keys = dframe.keys() #If it ain't you baby :)
        if weights in keys:
            col = dframe[weights]
            make_weights(col, "weights_on+"+ str(col)+".txt")
    print("Dataframe created with head:")
    print(dframe.head())
    return dframe

#df_radio = make_df_from_file("table_radio.txt", header_radio)
#df_IR = make_df_from_file("table_IR.txt", IR_header)


