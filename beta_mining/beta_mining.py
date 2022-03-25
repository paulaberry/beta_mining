#!/usr/bin/python3

import prody
import filter_proteins
import numpy as np
import pandas as pd
import os
import tarfile
import gzip
import shutil
import glob
from pathlib import Path
import argparse
import re

def get_args():
    """Function to pass in AlphaFold2 .tar locations."""
    getfiles = argparse.ArgumentParser(description="This script will take .tar files downloaded from AlphaFold2 database and produce a datatable of possible beta sheet targets.")
    getfiles.add_argument("-t", "--tarball", help = "Filepath(s) to .tar files downloaded from AlphaFold2.", required = True, nargs = "*")
    getfiles.add_argument("-o", "--outputDir", help = "Output path for the data table generated. Current format is comma-separated values.", required = False, default = "")
    return getfiles.parse_args()
args = get_args()


#test command: beta_mining -t /home/pberry/Halfmann/AlphaFold2/BM0000_0000_TEST.tar
# required arguments
tarList = args.tarball

if args.outputDir == "":
    outputPath = "betaMiningHits.csv"
else:
    outputPath = str(args.outputDir)

# set the parameters
contactsAngstroms = 6
contactsExclusion = 2
betasAngstroms = 14
betasExclusion = 4
gap = 6
flank = 5
complete_CSV = True
betaRegEx = re.compile("β-sheetβ-sheetβ-sheet.+β-sheetβ-sheetβ-sheet")# regex for betasheet patterns

# create output file and column headers
outputTable = open(outputPath, "w")
outputTable.write("proteomeID,taxonomy,organism,proteinID,fragment,version,res,res_n,phi,psi,omega,twist,SecondaryStructure,X,Y,Z,pLDDT,contacts,beta_neighbors\n")
outputTable.close()

# make temporary outdirectory to hold extracted .gz files
if os.path.exists("temp_BetaMining") == False:
    os.mkdir("temp_BetaMining")

# extract individual .gz files to the temp directory for processing into a temporary directory
for archive in tarList:
    pathMembers = archive.split("/")
    organism = pathMembers[-1].split("_")[-1][0:-4]
    proteomeID = pathMembers[-1].split("_")[0]
    taxonomy = pathMembers[-1].split("_")[1]
    t = tarfile.open(archive, "r")
    # extract individual .gz files to the temp directory for processing
    filecount = 1
    filenames = list(filter(lambda element: "pdb.gz" in element, t.getnames()))
    totalcount = len(filenames)
    for member in filenames:
        protein = member.split("-")[1]
        fragment = member.split("-")[2]
        modelN = member.split("-")[3][0:-7]
        print("Processing model " + str(filecount) + " of " + str(totalcount) + ", in " + proteomeID + " (" + organism + "): " + protein + ", " + fragment + ", " + modelN + ".")
        t.extract(member, "temp_BetaMining")
        with gzip.open("temp_BetaMining/" + member, "rb") as file_in:
            with open("temp_BetaMining/" + member[0:-3], "wb") as file_out:
                shutil.copyfileobj(file_in, file_out)
        model = prody.parsePDB("temp_BetaMining/" + member[0:-3])
        protein_dataframe, protein_dictionary = filter_proteins.protein_df(proteomeID, taxonomy, organism, protein, fragment, modelN, model, contactsAngstroms, contactsExclusion, betasAngstroms, betasExclusion)
        if complete_CSV == True:
            protein_dataframe.to_csv(protein + "_complete.csv", index = False)
        betasheet_dataframe = protein_dataframe[protein_dataframe["res_n"].isin(filter_proteins.filter_df(protein_dataframe, protein_dictionary, betaRegEx, gap, flank))]
        betasheet_dataframe.to_csv(outputPath, mode = "a", index=False, header=False)

        os.remove("temp_BetaMining/" + member)
        os.remove("temp_BetaMining/" + member[0:-3])
            #time.stop()
        filecount = filecount + 1
