#!/usr/bin/env python
import numpy as np
import pandas as pd
from statistics import mean
import os
import shutil
import glob
from pathlib import Path
import argparse
import re

def get_args():
    """Function to pass in arguments"""
    getfiles = argparse.ArgumentParser(description="This script takes in files produced by BetaMining.py and sorts them into FASTA formats and provides annotation.")
    getfiles.add_argument("-f", "--files", help = "Prefixes or filenames of the output datatables to parse and sort. Wildcards accepted.", required = True, nargs = "*")
    getfiles.add_argument("-a", "--annotation", help = "File with annotations for Uniprot IDs. First column should be the UniProtIDs.", required = True)
    getfiles.add_argument("-o", "--output", help = "Output path for the data table generated. Current format is comma-separated values.", required = False, default = "_annotation.csv")
    return getfiles.parse_args()
args = get_args()

#test commands: ./BetaMiningAnnotate.py -f betaHits_YEAST.csv -a resources/UP000002311_annotations.csv -o betaHits_YEAST_annotation_12072021.csv

#./BetaMiningAnnotate.py -f betaHits_HUMAN0*_12072021.csv -a resources/UP000005640_annotations.csv -o betaHits_HUMAN_annotation_12072021.csv

# columns of datatables to parse are: "proteome", "taxonomy", "organism", "proteinID", "fragment","version", "res","res_n", "phi", "psi", "omega","sec_str","X", "Y", "Z", "pLDDT", "contacts","beta_neighbors"

#first add KEGG, GO, GeneNames to datatables.

#establish uniprot annotation dataframes
df_annotations = pd.read_csv(args.annotation)
df_annotations.fillna("", inplace=True)
df_annotations.rename(columns = {df_annotations.columns[0]: "proteinID"}, inplace = True)
annotateList = list(df_annotations.columns) #will be used to grab column info for the wide format dataframe
#import the csv outputs from BetaMining.py


df_n = 1
for item in args.files:
    if df_n == 1:
        df_output = pd.read_csv(item)
        df_n = df_n + 1
    else:
        df_temp = pd.read_csv(item)
        df_output = pd.concat([df_output, df_temp])




df_output.fillna("", inplace = True)
df_output.sort_values(by = ["proteinID", "fragment","res_n"], ascending = [True, True, True], inplace = True)

df_outputAnnotated = pd.merge(df_output, df_annotations, on = "proteinID", how = "left")
df_outputAnnotated.fillna("", inplace = True)

df_outputAnnotated.to_csv("COMPLETE_" + args.output, index=False)

letterDict = {"ARG": "R", "HIS": "H", "LYS": "K", "ASP": "D", "GLU": "E", "SER": "S", "THR": "T", "GLN": "Q", "CYS": "C", "GLY": "G", "PRO": "P", "ALA": "A", "VAL": "V", "ILE": "I", "LEU": "L", "MET": "M", "PHE": "F", "TYR": "Y", "TRP": "W", "ASN": "N"}

letterDF = pd.DataFrame.from_dict(letterDict, orient = "index", columns=["AA"])
letterDF.reset_index(inplace=True)
letterDF = letterDF.rename(columns = {"index": "res"})

#add index identifiers for each sequence to use to group by
df_outputAnnotated["SeqID"] = ""
SeqCount = 0
ResNumber = -2
proteome, taxonomy, organism, proteinID, fragment, version = "", "", "", "", "", ""
idxList = []

for index, row in df_outputAnnotated.iterrows():
    if (proteome, taxonomy, organism, proteinID, fragment, version) != (row["proteomeID"], row["taxonomy"], row["organism"],row["proteinID"],row["fragment"], row["version"]) or row["res_n"] - ResNumber != 1:
        SeqCount = SeqCount + 1
        idxList.append(SeqCount)
        df_outputAnnotated.loc[index, "SeqID"] = SeqCount
        ResNumber = row["res_n"]
        proteome, taxonomy, organism, proteinID, fragment, version = row["proteomeID"], row["taxonomy"], row["organism"],row["proteinID"], row["fragment"], row["version"]
    else:
        df_outputAnnotated.loc[index, "SeqID"] = SeqCount
        ResNumber = row["res_n"]


# twist interval -5 to 30 on at least 3 sequential beta sheet residues
pd.options.mode.chained_assignment = None
for i in idxList:
    tmp_df = df_outputAnnotated[(df_outputAnnotated["SeqID"] == i)]
    tmp_df
    #tmp_df["mask"] = 1
    for index, row in tmp_df.iterrows():
        if (-5 <= row["twist"] <= 30) and row["SecondaryStructure"] == "β-sheet":
            tmp_df.loc[index,"mask"] = 0
        else:
            tmp_df.loc[index, "mask"] = 1

    tmp_df["CumSum"] = tmp_df["mask"].cumsum()
    tmp_df_group = tmp_df[tmp_df["mask"] == 0].groupby("CumSum")
    size = tmp_df_group.size()
    if len(size[size > 3]) == 0:
        idxList.remove(i)
df_outputAnnotated = pd.merge(df_outputAnnotated, pd.DataFrame(idxList, columns = ["SeqID"]), on = "SeqID", how = "right")
df_outputAnnotated.to_csv("TWISTFILTER_" + args.output, index=False)

#convert to wide and fasta format
wideCols = ["SeqID","proteomeID", "taxonomy", "organism", "proteinID", "fragment", "version", "pLDDTmean", "offtargetpLDDTmean", "targetpLDDTmean"] + annotateList
print(wideCols)
df_outputAnnotated = pd.merge(df_outputAnnotated, letterDF, on = "res", how = "left")
df_outputWide = pd.merge(df_outputAnnotated[wideCols].drop_duplicates(), df_outputAnnotated.groupby("SeqID")["res_n"].apply(list).apply(min).reset_index().rename(columns = {"res_n": "Start"}), on = "SeqID", how = "left")
#df_resMin = df_resMin.rename(columns = {"res_n": "Start"})

df_outputWide = pd.merge(df_outputWide, df_outputAnnotated.groupby("SeqID")["res_n"].apply(list).apply(max).reset_index().rename(columns = {"res_n": "Stop"}), on = "SeqID", how = "left")

df_outputWide = pd.merge(df_outputWide, df_outputAnnotated.groupby("SeqID")["AA"].apply(list).apply("".join).reset_index(), on = "SeqID", how = "left").drop(columns = ["SeqID"])

df_outputWide.to_csv("TWISTFILTER_WIDE_" + args.output, index=False)
#df_outputWIDE = pd.merge(df_)
#print(df_outputAnnotated)
