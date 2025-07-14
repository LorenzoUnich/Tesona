##code to help reformatting the catalogue

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt



datastr = "catalogo.txt"

def count_commas(string):
    """
    Conta le virgole in ciascuna stringa della lista.

    Args:
        strings (list of str): Lista di stringhe da analizzare.

    Returns:
        list of int: Lista con il numero di virgole per ogni stringa.
    """
    return string.count(',') 
file = open(datastr)
lines = file.readlines()

csvlines = []

for index, line in enumerate(lines):
    csvline = ""
    char_n = 0
    l = len(line)
    chunk = ""

    while char_n < l:
        char = line[char_n]

        if char == " ":
            try:
                next_char = line[char_n + 1]

                if next_char == " ":
                    # doppio spazio → fine campo
                    csvline += chunk.strip() + ","
                    chunk = ""
                    # salta tutti gli spazi successivi
                    while char_n < l and line[char_n] == " ":
                        char_n += 1
                else:
                    chunk += char
                    char_n += 1
            except IndexError:
                break
        else:
            chunk += char
            char_n += 1

    # Aggiungi l'ultimo campo, se c'è
    if chunk:
        csvline += chunk.strip()
    count = count_commas(csvline)
    if  count== 10:
        csvlines.append(csvline)
    else:
        print(count)
        csvline+=", "
        csvlines.append(csvline)
file.close()


with open(datastr) as file:
    lines = file.readlines()

csvlines = []

for index, line in enumerate(lines):

    csvline = ""
    char_n = 0
    l = len(line)
    chunk = ""

    while char_n < l:
        char = line[char_n]

        if char == " ":
            try:
                next_char = line[char_n + 1]

                if next_char == " ":
                    csvline += chunk.strip() + ","
                    chunk = ""
                    while char_n < l and line[char_n] == " ":
                        char_n += 1
                else:
                    chunk += char
                    char_n += 1
            except IndexError:
                break
        else:
            chunk += char
            char_n += 1

    if chunk:
        csvline += chunk.strip()


    if csvline.strip().strip(","):
        csvlines.append(csvline)


with open("catalogue_gpt.csv", "w", newline='') as f:
    writer = csv.writer(f)
    for line in csvlines:
        writer.writerow(line.split(","))
no_empty_lines = []
for line in csvlines:
    if  not line =="\n":
        no_empty_lines.append(line)
with open("catalogue.csv", "w", encoding="utf-8") as f:
    for line in no_empty_lines:
        f.write(line.strip() + "\n")  # Rimuove eventuali spazi o newline extra
f.close()
