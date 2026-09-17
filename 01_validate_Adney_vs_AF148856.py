# -*- coding: utf-8 -*-
"""
Created on Thu Aug 27 12:36:05 2026

@author: kdalw
"""

# ============================================================
# LINE-1 FUNCTIONAL ANNOTATION PROJECT
#
# STEP 1:
# Validate Adney et al. trialanine mutants against
# the L1-RP wild-type reference sequence.
#
# WT reference:
# L1-RP
# GenBank accession: AF148856
#
# ORF1p protein: AAD39214.1
# ORF2p protein: AAD39215.1
# ============================================================


# ============================================================
# IMPORT PACKAGES
# ============================================================

import pandas as pd
from pathlib import Path
from Bio import Entrez, SeqIO


# ============================================================
#              EDIT THESE SETTINGS
# ============================================================

# NCBI asks users of Entrez to provide an email address.
# Replace this with your actual email.
Entrez.email = "manaal.dalwadi@cuanschutz.edu"


# Name of the Excel workbook containing Supplemental Table 2.
# If the Excel file is in the SAME folder as this Python script,
# you do not need to change anything here.
INPUT_FILE_NAME = "Adney_Supplemental_Table2_LINE1_Master.xlsx"


# Sheet containing the complete Table 2 dataset.
INPUT_SHEET_NAME = "Master_Table2"


# Name of the output file that will be created.
OUTPUT_FILE_NAME = "Adney_Table2_AF148856_Validated.xlsx"


# ============================================================
#           DO NOT NEED TO EDIT BELOW HERE
# ============================================================


# ============================================================
# DEFINE L1-RP REFERENCE ACCESSIONS
# ============================================================

L1_REFERENCE = "AF148856.1"

REFERENCE_PROTEINS = {
    "ORF1": "AAD39214.1",
    "ORF2": "AAD39215.1"
}


# ============================================================
# SET FILE LOCATIONS
# ============================================================

# Find the folder containing this Python script.
SCRIPT_FOLDER = Path(__file__).resolve().parent

# Input Excel file
INPUT_FILE = SCRIPT_FOLDER / INPUT_FILE_NAME

# Output Excel file
OUTPUT_FILE = SCRIPT_FOLDER / OUTPUT_FILE_NAME


# ============================================================
# CHECK THAT INPUT FILE EXISTS
# ============================================================

if not INPUT_FILE.exists():

    raise FileNotFoundError(
        f"\nCould not find the input Excel file:\n"
        f"{INPUT_FILE}\n\n"
        f"Make sure '{INPUT_FILE_NAME}' is in the same "
        f"folder as this Python script."
    )


print("\n====================================================")
print("LINE-1 / AF148856 VALIDATION")
print("====================================================")

print("\nInput file:")
print(INPUT_FILE)

print("\nReference:")
print("L1-RP =", L1_REFERENCE)


# ============================================================
# FUNCTION:
# DOWNLOAD PROTEIN SEQUENCE FROM NCBI
# ============================================================

def fetch_protein_sequence(accession):

    """
    Download a protein sequence from NCBI.

    Parameters
    ----------
    accession : str
        NCBI protein accession.

    Returns
    -------
    str
        Amino-acid sequence.
    """

    print(f"\nDownloading {accession} from NCBI...")

    with Entrez.efetch(
        db="protein",
        id=accession,
        rettype="fasta",
        retmode="text"
    ) as handle:

        record = SeqIO.read(handle, "fasta")

    sequence = str(record.seq).upper()

    print(
        f"Downloaded {accession}: "
        f"{len(sequence)} amino acids"
    )

    return sequence


# ============================================================
# DOWNLOAD WT L1-RP PROTEINS
# ============================================================

orf1_sequence = fetch_protein_sequence(
    REFERENCE_PROTEINS["ORF1"]
)

orf2_sequence = fetch_protein_sequence(
    REFERENCE_PROTEINS["ORF2"]
)


# Store them in a dictionary so we can select
# the correct sequence based on ORF.
REFERENCE_SEQUENCES = {
    "ORF1": orf1_sequence,
    "ORF2": orf2_sequence
}


# ============================================================
# CHECK EXPECTED PROTEIN LENGTHS
# ============================================================

print("\nChecking reference protein lengths...")


if len(orf1_sequence) != 338:

    print(
        "WARNING:",
        f"ORF1p is {len(orf1_sequence)} aa.",
        "Expected 338 aa."
    )

else:

    print("ORF1p length = 338 aa: PASS")


if len(orf2_sequence) != 1275:

    print(
        "WARNING:",
        f"ORF2p is {len(orf2_sequence)} aa.",
        "Expected 1275 aa."
    )

else:

    print("ORF2p length = 1275 aa: PASS")


# ============================================================
# LOAD SUPPLEMENTAL TABLE 2
# ============================================================

print("\nLoading Supplemental Table 2...")


df = pd.read_excel(
    INPUT_FILE,
    sheet_name=INPUT_SHEET_NAME
)


print(f"Rows loaded: {len(df)}")


# ============================================================
# CHECK REQUIRED COLUMNS
# ============================================================

required_columns = [

    "ORF",
    "Mutant_ID",
    "Start_AA",
    "End_AA",
    "WT_Residues",
    "RetroT_Average"

]


missing_columns = [

    column
    for column in required_columns
    if column not in df.columns

]


if missing_columns:

    raise ValueError(
        "\nThe following required columns are missing "
        "from the Excel file:\n"
        + ", ".join(missing_columns)
    )


# ============================================================
# FUNCTION:
# GET WT RESIDUES FROM AF148856
# ============================================================

def get_reference_residues(row):

    """

    Look up the amino acids from the AF148856-derived
    protein sequence at the positions reported by
    Adney et al.

    """

    protein = row["ORF"]

    start = int(row["Start_AA"])

    end = int(row["End_AA"])


    if protein not in REFERENCE_SEQUENCES:

        return None


    sequence = REFERENCE_SEQUENCES[protein]


    # Protein numbering starts at 1.
    # Python string indexing starts at 0.
    #
    # Example:
    #
    # Residues 2-4
    #
    # Python:
    # sequence[1:4]

    reference_residues = sequence[
        start - 1:end
    ]


    return reference_residues


# ============================================================
# EXTRACT AF148856 WT RESIDUES
# ============================================================

print("\nExtracting WT residues from AF148856...")


df["AF148856_WT_Residues"] = df.apply(
    get_reference_residues,
    axis=1
)


# ============================================================
# COMPARE ADNEY WT RESIDUES TO AF148856
# ============================================================

df["AF148856_Match"] = (

    df["WT_Residues"]
    .astype(str)
    .str.upper()

    ==

    df["AF148856_WT_Residues"]
    .astype(str)
    .str.upper()

)


# ============================================================
# ADD REFERENCE INFORMATION
# ============================================================

df["Reference_L1"] = L1_REFERENCE


df["Reference_Protein"] = df["ORF"].map(
    REFERENCE_PROTEINS
)


# ============================================================
# ADD RETROT CONTEXT
#
# IMPORTANT:
# These are descriptive project labels.
# We are NOT deleting any mutants.
# ============================================================

def assign_retrot_context(value):

    """

    Add a descriptive context to each RetroT value.

    <100%
        Reduced relative to WT

    100-125%
        WT-adjacent

    >125%
        Elevated relative to WT

    """

    if pd.isna(value):

        return "Missing RetroT value"


    value = float(value)


    if value < 100:

        return "Reduced relative to WT"


    elif value <= 125:

        return "WT-adjacent"


    else:

        return "Elevated relative to WT"


df["RetroT_Context"] = (

    df["RetroT_Average"]
    .apply(assign_retrot_context)

)


# ============================================================
# CREATE VALIDATION STATUS COLUMN
# ============================================================

df["Validation_Status"] = df["AF148856_Match"].map(

    {
        True: "PASS",
        False: "MISMATCH"
    }

)


# ============================================================
# GENERATE VALIDATION SUMMARY
# ============================================================

total_mutants = len(df)

matches = int(
    df["AF148856_Match"].sum()
)

mismatches = (
    total_mutants - matches
)


orf1_count = (
    df["ORF"] == "ORF1"
).sum()


orf2_count = (
    df["ORF"] == "ORF2"
).sum()


print("\n====================================================")
print("VALIDATION RESULTS")
print("====================================================")


print(
    f"\nTotal mutants tested: {total_mutants}"
)

print(
    f"ORF1 mutants: {orf1_count}"
)

print(
    f"ORF2 mutants: {orf2_count}"
)

print(
    f"\nAF148856 matches: {matches}"
)

print(
    f"AF148856 mismatches: {mismatches}"
)


# ============================================================
# SHOW ANY MISMATCHES
# ============================================================

mismatch_df = df[
    df["AF148856_Match"] == False
].copy()


if len(mismatch_df) == 0:

    print(
        "\nPASS:"
        "\nAll Adney WT residue windows match "
        "the L1-RP reference sequence."
    )


else:

    print(
        "\nWARNING:"
        "\nSome mutant residue windows do NOT "
        "match AF148856."
    )


    print(
        mismatch_df[
            [
                "ORF",
                "Mutant_ID",
                "Start_AA",
                "End_AA",
                "WT_Residues",
                "AF148856_WT_Residues"
            ]
        ].to_string(index=False)
    )


# ============================================================
# RETROT CATEGORY SUMMARY
# ============================================================

print("\n====================================================")
print("RETROT CONTEXT")
print("====================================================")


context_counts = (

    df["RetroT_Context"]
    .value_counts()

)


print(context_counts)


# ============================================================
# CREATE SUMMARY TABLE FOR EXCEL
# ============================================================

summary_data = {

    "Metric": [

        "WT nucleotide reference",

        "ORF1p protein accession",

        "ORF2p protein accession",

        "Total trialanine mutants",

        "ORF1 mutants",

        "ORF2 mutants",

        "AF148856 matches",

        "AF148856 mismatches"

    ],


    "Value": [

        L1_REFERENCE,

        REFERENCE_PROTEINS["ORF1"],

        REFERENCE_PROTEINS["ORF2"],

        total_mutants,

        orf1_count,

        orf2_count,

        matches,

        mismatches

    ]

}


summary_df = pd.DataFrame(
    summary_data
)


# ============================================================
# SAVE RESULTS TO EXCEL
# ============================================================

print("\nSaving validated workbook...")


with pd.ExcelWriter(
    OUTPUT_FILE,
    engine="openpyxl"
) as writer:


    # Complete dataset
    df.to_excel(

        writer,

        sheet_name="Validated_Table2",

        index=False

    )


    # Validation summary
    summary_df.to_excel(

        writer,

        sheet_name="Validation_Summary",

        index=False

    )


    # Any mismatches
    mismatch_df.to_excel(

        writer,

        sheet_name="Mismatches",

        index=False

    )


print("\n====================================================")
print("DONE")
print("====================================================")


print(
    "\nValidated workbook saved as:"
)

print(
    OUTPUT_FILE
)


print(
    "\nYou can now use this validated dataset "
    "for natural LINE-1 mutation mapping."
)