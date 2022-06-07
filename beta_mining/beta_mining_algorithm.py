#!/usr/bin/env python3
#
# This script provides an interface that combines the functions in the
# "BetaMining" package to analyze .pdb files of interest.
import prody
import numpy as np
import pandas as pd

from Bio.PDB.PDBParser import PDBParser# as p
p = PDBParser(PERMISSIVE=1)

from sklearn.metrics import pairwise_distances

from biopandas.pdb import PandasPdb
import datetime
import os
import sys
import math
import re
import json
import os
import tarfile
import gzip
import shutil
import glob
from pathlib import Path
import importlib.resources

import beta_mining
from beta_mining import beta_mining_algorithm
from beta_mining import beta_mining_functions


def analyze_structure(filename, config, json, output_dictionary):
    """Main beta_mining algorithm.

    Keyword arguments:
    filename -- the .pdb or .pdb.gz structure file being analyzed
    config -- the configuration YAML file
    json -- the json used to guide the analysis
    output_dictionary -- dictionary of output files to generate
    """
    # create ProDy model, meta information dictionary, BioPandas Dataframe and AA sequence list
    model, meta_dictionary, af_object, af_sequence = beta_mining_functions.create_model_metainfo(filename)

    # create the dataframe of dihedral angles, pandas series of secondary structures, and a dictionary with secondary structures as keys and a list of residue numbers as values.
    dihedral_df, secondary_structure_series, secondary_structure_reference = beta_mining_functions.calculation_df(model, meta_dictionary["fragment_offset"], json["secondary_structures"], config["units"])

    # create dataframe of per-residue contacts specified in the features json
    # create a dictionary of boolean Series, keys are "include" or "exclude"
    for target in json["target_region_features"]:
        if target["name"] in config["target_names"]:
            if "contacts" in list(target["include"].keys()) +  list(target["exclude"].keys()):
                contacts_dataframe, mask_dict = beta_mining_functions.contacts_df(af_object, json["target_region_features"], secondary_structure_reference, config["target_names"])
                # replace features from include with the specified symbol first
                # replace features from exclude with the specified symbol last
                for condition in config["conditions"]:
                    for series_symbol in mask_dict[condition]:
                        secondary_structure_series.mask(series_symbol[0], series_symbol[1], inplace = True)
                # left join the calculated dihedral and twist dataframe and the contacts dataframe
                # create a secondary structure symbols string to search with regex
                dihedral_df = pd.merge(dihedral_df, contacts_dataframe, on = "residue_number", how = "left")

    structure_symbols_string = "".join(secondary_structure_series)
    structure_length = len(structure_symbols_string)

    # generate dataframe of entire protein
    residue_df = polymer_df(meta_dictionary, af_object)
    residue_df = pd.merge(residue_dataframe, dihedral_df, on = "residue_number", how = "left")
    residue_df["structure_symbol"] = secondary_structure_series
    if output_dictionary["proteome_aa"] != False:
        output_filename = config["output_filepath"] + config["results_prefix"] + meta_dictionary["accession"] + "_f" + meta_dictionary["fragment"] + "_" + meta_dictionary["depo_date"].lower() + config["output_filename"]
        residue_df.to_csv(output_filename, index = False)


    # look for regex matches in the secondary structure string for each type of target
    for target in json["target_region_features"]:
        if target["name"] in config["target_names"]:
            compiled_regex = re.compile(target["regex"])
            targets_found = re.finditer(compiled_regex, structure_symbols_string)
            fasta_fields = [
                            target["name"],
                            meta_dictionary["id_code"],
                            meta_dictionary["accession"],
                            meta_dictionary["organism_scientific"],
                            meta_dictionary["organism_taxid"],
                            meta_dictionary["full_title"],
                            meta_dictionary["depo_date"],
                            "fragment " + str(meta_dictionary["fragment"])
                        ]
            regex_flank = target["regex_flank"]
            for match_obj in targets_found:
                if attribute_filter(target, match_obj, dihedral_df) == False:
                    print("attribute_filter was False!")
                    continue
                else:
                    # these values are the 0 index of residues for slicing
                    res_idx_start = max(0, match_obj.span()[0] - regex_flank)
                    res_idx_end = min(structure_length - 1, match_obj.span()[1] + regex_flank)

                    # these values are the 1 index of residues for output
                    res_start = res_idx_start + 1
                    res_end = res_idx_start + 1

                    # sequences of amino acids and structures
                    aa_sequence = "".join(af_sequence[res_idx_start:res_idx_end])
                    ss_sequence = structure_symbols_string[res_idx_start:res_idx_end]

                    if output_dictionary["hits_fasta"] != False:
                        fasta_header = ">" + "|".join(fasta_fields) + "residues " + str(res_start) + "-" + str(res_end)
                        output_dictionary["hits_fasta"].writelines([fasta_header, aa_sequence])
                    if output_dictionary["hits_aa"] != False:
                        mean_plddt = 
                        mean_twist =
                        hit_line = fasta_fields +
                            [str(res_start),
                            str(res_end),
                            aa_sequence,
                            ss_sequence,

                            ]



def main(config_settings):
    """Filehandling algorithm.

    Keyword arguments:
    config_dictionary -- the dictionary built from the config YAML file given to the algorithm.
    """
    ### set global run variables from config YAML ###
    input_path = config_settings["input_filepath"]
    output_path = config_settings["output_filepath"]
    results_prefix = config_settings["results_prefix"]
    conditions = config_settings["conditions"]

    ## import parameter json ##
    if config_settings["custom_json"] == None:
        with importlib.resources.open_text("beta_mining", "structure_dictionary.json") as param_json:
            target_parameters = json.load(param_json)
    else:
        with open(config_settings["custom_json"], "r") as param_json:
            target_parameters = json.load(param_json)

    ### prepare output files ###
    #hits_fasta
    if "hits_fasta" in config_settings["save_files"]:
        hits_fasta = open(output_path + results_prefix + ".fasta", "w")
    else:
        hits_fasta = False
    #hits_aa
    if "hits_aa" in config_settings["save_files"]:
        hits_aa = open(output_path + results_prefix + "_aa.csv", "w")
    else:
        hits_aa = False
    #proteome_aa files are generated in analyze_structure function
    if "proteome_aa" in config_settings["save_files"]:
        proteome_aa = True
    else:
        proteome_aa = False

    output_files = {"hits_fasta": hits_fasta, "hits_aa": hits_aa, "proteome_aa": proteome_aa}

    ### read in files and process different formats ###
    if os.path.isdir(output_path) == False:
        os.mkdir(output_path)
    file_list = glob.glob(input_path + "**.pdb*", recursive = True)
    file_number = 1
    file_total = len(file_list)
    for file in file_list:
        print("Mining file " + file + ": "+ str(file_number) + " of " + str(file_total) + "...")
        analyze_structure(file, config_settings, target_parameters, output_files)
        file_number = file_number + 1

    #close the file objects
    hits_fasta.close()
    hits_aa.close()




if __name__ == "__main__":
    main(config_settings)
