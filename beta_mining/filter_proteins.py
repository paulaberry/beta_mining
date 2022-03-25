#!/usr/bin/env python3
#
# This script provides an interface that combines the functions in the
# "BetaMining" package to analyze .pdb files of interest.
import os
import re
import shutil
import prody
import numpy as np
import pandas as pd


def assign_secondary_structure(phi, psi, omega):
    """returns a secondary structure ("α-helix", "polyproline II", "β-sheet",
    "left-handed", "γ", "cis") or an empty string

    Keyword arguments:
    phi -- residue phi angle in degrees
    psi -- residue psi angle in degrees
    omega -- residue omega angle in degrees
    """
    if -200 <= phi <= 0 and -120 <= psi <= 40:
        return "α-helix"
    elif -90 <= phi <= 0 and 40 <= psi <= 240:
        return "polyproline II"
    elif -200 <= phi <= -90 and 40 <= psi <= 240:
        return "β-sheet"
    elif 0 <= phi <= 160 and -90 <= psi <= 110:
        return "left-handed"
    elif 0 <= phi <= 160 and 110 <= psi <= 270:
        return "γ"
    elif -90 <= omega <= 90:
        return "cis"
    else:
        return ""


def find_contacts(model, distance=6, exclude=2):
    """Return a two-column pandas dataframe of residue numbers ("res_n") and
    lists of residues ("contacts") within some distance in angstroms in the same
    chain, excluding a number of flanking residues.

    Keyword arguments:
    model -- an AtomGroup class from ProDy
    distance -- the maximum distance in angstroms between alpha carbons that
    will be considered to be in contact with each other; default value is 6
    angstroms
    exclude -- the number of residues on each side to exclude from the list of
    contacts; default value is 2 residues on each side
    """
    chain_id = model.getChids()[0]
    chain = model.ca[chain_id]
    data = []
    for aa in prody.findNeighbors(chain, distance, chain):
        data.append(
            (
                aa[0].getResnum(),
                aa[1].getResnum(),
                abs(aa[0].getResnum() - aa[1].getResnum()),
            )
        )
    df = pd.DataFrame(data, columns=["res_n", "contacts", "Difference"])
    df = df[df["Difference"] > exclude]  # only keep residue numbers outside of
    # the excluded flank
    df = df.groupby("res_n")["contacts"].apply(list).apply(sorted).reset_index()
    # group dataframe by residue and convert the contacts for the same residue
    # to a list, remove the "Difference" column
    return df


def find_betas(model, dictionary, distance=14, exclude=4):
    """Return a two-column pandas dataframe of residue numbers ("res_n") and a
    list of residues in contact with that residue that have been assigned a
    beta-sheet structure. ("beta_neighbors")

    Keyword arguments:
    model -- an AtomGroup class from ProDy
    dictionary -- a dictionary consisting of residue numbers as keys and their secondary structure assignment
    distance -- the maximum distance in angstroms between alpha carbons that will be considered to be neighbors; default value is 14 angstroms
    exclude -- the number of residues on each side to exclude from the list of contacts; default value is 4 residues on either side
    """
    df = find_contacts(model, distance, exclude)
    for index, row in df.iterrows():
        for r in row["contacts"]:
            if dictionary[r] != "β-sheet":
                row.contacts.remove(r)
        if row.contacts == []:
            df.drop(index, inplace=True)
    df = df.rename(columns={"contacts": "beta_neighbors"})
    return df


def protein_df(
    proteome,
    taxonomy,
    organism,
    protein_id,
    fragment,
    model_version,
    model,
    contact_distance,
    contact_exclude,
    beta_distance,
    beta_exclude,
):
    """Return a pandas dataframe of a protein, consisting of organism identifying columns ("proteomeID", "taxonomy", "organism"), protein and AlphaFold2 model identifying columns ("proteinID", "fragment", "version"), columns with conformation information about individual residues ("res", "res_n", "phi", "psi", "omega", "twist", "sec_str", "X", "Y", "Z") and the AlphaFold2 confidence score. ("pLDDT")

    Keyword arguments:
    proteome -- the UniProt Proteome identifier prefix, taken from the filename of the .tar downloaded from the AlphaFold Protein Structure Database
    taxonomy -- the four digit UniProt Taxon ID, taken from the filename of the .tar downloaded from the AlphaFold Protein Structure Database
    organism -- the common name of the organism, taken from the filename of the .tar downloaded from the AlphaFold Protein Structure Database
    protein_id -- the UniProt code for a given protein
    fragment -- the protein fragment a residue comes from, taken from the filename of the individual protein structure .pdb file
    model_version -- the prediction version present in the downloaded AF2 file
    contact_distance -- the maximum distance in angstroms between two residues' alpha carbons for them to be considered in contact
    contact_exclude -- the number of residues flanking another to exclude from contact calculation results
    beta_distance -- the maximum distance in angstroms between two beta-sheet assigned residues for them to be considered neighbors
    beta_exclude -- the number of beta-assigned residues flanking another beta-assigned residue to be excluded from the neighboring beta-sheet residue list
    """
    data = []
    dictionary = {}
    for residue in model.iterResidues():
        try:
            phi = prody.calcPhi(residue, radian=False, dist=None)
            print(phi)
            psi = prody.calcPsi(residue, radian=False, dist=None)
            print(psi)
            omega = prody.calcOmega(residue, radian=False, dist=None)
            print(omega)
        except:
            dictionary[residue.getResnum()] = ""
            print("No angles!")
            continue
        else:
            data.append(
                (
                    proteome,
                    taxonomy,
                    organism,
                    protein_id,
                    fragment,
                    model_version,
                    residue.getResname(),
                    residue.getResnum(),
                    phi,
                    psi,
                    omega,
                    psi - abs(phi),
                    assign_secondary_structure(phi, psi, omega),
                    prody.getCoords(residue.ca)[0][0],
                    prody.getCoords(residue.ca)[0][1],
                    prody.getCoords(residue.ca)[0][2],
                    residue.ca.getBetas()[0],
                )
            )
            dictionary[int(residue.getResnum())] = assign_secondary_structure(
                phi, psi, omega
            )
    df = pd.DataFrame(
        data,
        columns=[
            "proteomeID",
            "taxonomy",
            "organism",
            "proteinID",
            "fragment",
            "version",
            "res",
            "res_n",
            "phi",
            "psi",
            "omega",
            "twist",
            "sec_str",
            "X",
            "Y",
            "Z",
            "pLDDT",
        ],
    )
    df = pd.merge(
        df,
        find_contacts(model, contact_distance, contact_exclude),
        on="res_n",
        how="left",
    ).fillna("")
    df = pd.merge(
        df,
        find_betas(model, dictionary, beta_distance, beta_exclude),
        on="res_n",
        how="left",
    ).fillna("")
    return df, dictionary


def identify_beta_runs(beta_residue_list, reg_ex, flank):
    """Return a list of sequential residue numbers that fits the allowed beta-sheet run pattern. Return None if there are no residues that fit the allowed
    beta-sheet run pattern.

    Keyword arguments:
    beta_residue_list -- a list of residue numbers that have been assigned beta-
    sheet secondary structure
    reg_ex -- a compiled regular expression object of the desired beta-sheet run
    pattern
    flank -- the number of residues to add on each side of identified beta-sheet
    runs
    """
    beta_residue_range = sorted(set(beta_residue_list))
    secondary_structure_string = ""
    range_position_ceiling = len(beta_residue_list) - 1
    for number in range(range_position_ceiling):
        if beta_residue_range[number + 1] - beta_residue_range[number] == 1:
            secondary_structure_string = secondary_structure_string + "β-sheet"
        else:
            secondary_structure_string = secondary_structure_string + "β-sheetSPACE"
    if reg_ex.search(secondary_structure_string) != None:
        lowest_residue, highest_residue = (
            beta_residue_range[0] - flank,
            beta_residue_range[-1] + flank,
        )
        residue_list = [*range(lowest_residue, highest_residue + 1)]
        return residue_list
    else:
        return None


def filter_df(dataframe, dictionary, reg_ex, allowed_gap, flank):
    """Return a list of residues in a protein that match the desired beta-sheet
    run pattern.

    Keyword arguments:
    dataframe -- a pandas dataframe of protein residues and secondary structure
    information
    dictionary -- a dictionary with residue numbers as keys and the assigned
    secondary structure as values
    reg_ex -- a compiled regular expression object of the desired beta-sheet run
    pattern
    allowed_gap -- an integer value for the largest allowable gap between beta-
    sheet assigned residues
    This function takes in a dataframe, a dictionary of secondary structure by residue number, a compiled regular expression, and an integer for allowable gap in beta-sheet runs.
    """
    df = dataframe[
        (dataframe["contacts"] == "")
        & (dataframe["beta_neighbors"] != "")
        & (dataframe["sec_str"] == "β-sheet")
    ]
    residue_count = len(dictionary)
    beta_residues = list(df["res_n"])
    range_upper_limit = len(beta_residues) - 1
    master_beta_list = []  # heh
    beta_range = []
    for number in range(range_upper_limit):
        if beta_residues[number + 1] - beta_residues[number] > allowed_gap:
            beta_range.append(beta_residues[number])
            pattern_residues = identify_beta_runs(beta_range, reg_ex, flank)
            if pattern_residues != None:
                master_beta_list.extend(pattern_residues)  # HEH
            beta_range = []
        else:
            beta_range.append(beta_residues[number])
            if beta_residues[number + 1] == beta_residues[-1]:
                beta_range.append(beta_residues[number + 1])
                pattern_residues = identify_beta_runs(beta_range, reg_ex, flank)
                if pattern_residues != None:
                    master_beta_list.extend(pattern_residues)
    return master_beta_list
