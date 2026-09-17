# ============================================================
# LINE-1 FUNCTIONAL ANNOTATION PROJECT
#
# SCRIPT 06:
# Map Script 05 amino-acid differences to Adney et al.
# tri-alanine experimental windows
#
# PURPOSE
# ------------------------------------------------------------
# Script 05 tells us:
#
#     "At ORF2 residue 68, L1-RP has T and this locus has Y."
#
# Script 06 asks:
#
#     "Which Adney experimental window contains residue 68,
#      and what RetroT value was measured when THAT WINDOW was
#      replaced with alanines?"
#
# IMPORTANT INTERPRETATION
# ------------------------------------------------------------
# The Adney RetroT value belongs to the tri-alanine (or edge
# alanine) experimental construct. It is NOT the measured
# activity of the naturally observed genomic substitution.
#
# Example:
#
#     Natural difference: T68Y
#     Adney window: 67-69
#     Adney RetroT: 6% WT
#
# Correct interpretation:
#     T68Y occurs inside an experimentally sensitive window
#     whose alanine-substitution construct retained 6% WT
#     retrotransposition.
#
# Incorrect interpretation:
#     "T68Y has 6% activity."
#
# This script therefore produces WINDOW-LEVEL experimental
# context only. Chemical comparison of the natural amino acid
# with alanine is reserved for Script 07.
#
# ALSO IMPORTANT
# ------------------------------------------------------------
# The original Adney workbook contains older project columns
# such as Keep_For_Analysis / Working_Label. Script 06 ignores
# those columns and uses ALL 538 experimental constructs.
# ============================================================


# ============================================================
# IMPORTS
# ============================================================

import random
from pathlib import Path

import numpy as np
import pandas as pd
from Bio import SeqIO


# ============================================================
#                 EDIT THESE SETTINGS
# ============================================================

# Start with TEST.
# After we inspect the output, change to:
#
#     RUN_MODE = "FULL"
#
RUN_MODE = "FULL"

# Number of unique loci to include in TEST mode.
TEST_N_LOCI = 25

# Fixed seed = reproducible random test set.
RANDOM_SEED = 20260914

# These loci are forced into TEST mode when present.
TEST_LOCUS_IDS = [
    "hg38.RE.chr1.100199603-100206089.+",
    "hg38.RE.chr1.104770248-104776279.-",
    "hg38.RE.chr1.104843834-104849865.-",
]

# Project-specific neutral RetroT descriptors.
PROJECT_WT_ADJACENT_LOW = 100.0
PROJECT_WT_ADJACENT_HIGH = 125.0

# Write Excel summary workbook.
WRITE_EXCEL_SUMMARY = True


# ============================================================
# PROJECT FILENAMES / FOLDERS
# ============================================================

SCRIPT_FOLDER = Path(__file__).resolve().parent

SCRIPT05_FOLDER = (
    SCRIPT_FOLDER
    / "05_ORF_variant_alignment"
)

SCRIPT04_FOLDER = (
    SCRIPT_FOLDER
    / "04_ORF_integrity"
)

SCRIPT03_FOLDER = (
    SCRIPT_FOLDER
    / "03_LINE1_structure"
)

OUTPUT_FOLDER = (
    SCRIPT_FOLDER
    / "06_Adney_mapping"
)

OUTPUT_FOLDER.mkdir(
    parents=True,
    exist_ok=True,
)


# ============================================================
# INPUT FILE RESOLUTION
# ============================================================

def first_existing(paths, description):
    for path in paths:
        if path.exists():
            return path

    attempted = "\n".join(
        f"  - {path}"
        for path in paths
    )

    raise FileNotFoundError(
        f"\nCould not find {description}.\n"
        f"Checked:\n{attempted}\n"
    )


AA_DIFFERENCES_FILE = first_existing(
    [
        SCRIPT05_FOLDER
        / "FULL_LINE1_AA_differences_vs_L1RP.csv.gz",

        SCRIPT_FOLDER
        / "FULL_LINE1_AA_differences_vs_L1RP.csv.gz",
    ],
    "Script 05 FULL AA-difference table",
)


ORF_INTEGRITY_FILE = first_existing(
    [
        SCRIPT04_FOLDER
        / "FULL_LINE1_ORF_integrity_results.csv.gz",

        SCRIPT04_FOLDER
        / "FULL_LINE1_ORF_integrity_results.csv",

        SCRIPT_FOLDER
        / "FULL_LINE1_ORF_integrity_results.csv.gz",

        SCRIPT_FOLDER
        / "FULL_LINE1_ORF_integrity_results.csv",
    ],
    "Script 04 FULL ORF-integrity table",
)


# Prefer the previously validated Adney workbook if present.
# Fall back to the master extraction because Script 06 performs
# its own independent WT-sequence validation against AF148856.1.
ADNEY_WORKBOOK = first_existing(
    [
        SCRIPT_FOLDER
        / "Adney_Table2_AF148856_Validated.xlsx",

        SCRIPT_FOLDER
        / "Adney_Supplemental_Table2_LINE1_Master.xlsx",
    ],
    "Adney Supplemental Table 2 workbook",
)


REFERENCE_GB = first_existing(
    [
        SCRIPT03_FOLDER
        / "AF148856.1.gb",

        SCRIPT_FOLDER
        / "AF148856.1.gb",
    ],
    "AF148856.1 GenBank reference",
)


# ============================================================
# OUTPUT FILES
# ============================================================

MODE = RUN_MODE.upper().strip()

if MODE not in {
    "TEST",
    "FULL",
}:
    raise ValueError(
        "RUN_MODE must be either 'TEST' or 'FULL'."
    )


PREFIX = MODE


MAPPED_AA_CSV_GZ = (
    OUTPUT_FOLDER
    / f"{PREFIX}_LINE1_AA_differences_Adney_mapped.csv.gz"
)

PRIORITY_CSV_GZ = (
    OUTPUT_FOLDER
    / f"{PREFIX}_LINE1_Adney_priority_clean_ORFs.csv.gz"
)

LOCUS_ORF_SUMMARY_CSV_GZ = (
    OUTPUT_FOLDER
    / f"{PREFIX}_LINE1_Adney_locus_ORF_summary.csv.gz"
)

LOCUS_SUMMARY_CSV_GZ = (
    OUTPUT_FOLDER
    / f"{PREFIX}_LINE1_Adney_locus_summary.csv.gz"
)

WINDOW_SUMMARY_CSV = (
    OUTPUT_FOLDER
    / f"{PREFIX}_LINE1_Adney_window_summary.csv"
)

SUMMARY_XLSX = (
    OUTPUT_FOLDER
    / f"{PREFIX}_LINE1_Adney_mapping_summary.xlsx"
)


if MODE == "TEST":

    MAPPED_AA_CSV = (
        OUTPUT_FOLDER
        / "TEST_LINE1_AA_differences_Adney_mapped.csv"
    )

    PRIORITY_CSV = (
        OUTPUT_FOLDER
        / "TEST_LINE1_Adney_priority_clean_ORFs.csv"
    )


# ============================================================
# HELPER FUNCTIONS
# ============================================================

def load_adney_table(workbook_path):
    """
    Find the worksheet containing the Adney Table 2 extraction.

    This makes the script compatible with either:
      - Adney_Table2_AF148856_Validated.xlsx
      - Adney_Supplemental_Table2_LINE1_Master.xlsx
    """

    required_columns = {
        "ORF",
        "Mutant_ID",
        "Start_AA",
        "End_AA",
        "WT_Residues",
        "RetroT_Average",
        "Standard_Deviation",
        "N_Measurements",
    }


    excel_file = pd.ExcelFile(
        workbook_path
    )


    preferred_sheet_names = [
        "Validated_Table2",
        "Master_Table2",
        "Table2",
        "Adney_Table2",
    ]


    ordered_sheets = []

    for sheet in preferred_sheet_names:
        if sheet in excel_file.sheet_names:
            ordered_sheets.append(
                sheet
            )


    for sheet in excel_file.sheet_names:
        if sheet not in ordered_sheets:
            ordered_sheets.append(
                sheet
            )


    for sheet in ordered_sheets:

        candidate = pd.read_excel(
            workbook_path,
            sheet_name=sheet,
        )


        if required_columns.issubset(
            candidate.columns
        ):

            return (
                candidate.copy(),
                sheet,
            )


    raise ValueError(
        "\nCould not find an Adney worksheet containing "
        "the required Table 2 columns."
    )


def extract_reference_proteins(
    genbank_path
):
    """
    Read AF148856.1 and obtain the annotated L1-RP ORF1p and
    ORF2p protein sequences.
    """

    record = SeqIO.read(
        genbank_path,
        "genbank",
    )


    proteins = {}


    for feature in record.features:

        if feature.type != "CDS":
            continue


        gene = (
            feature.qualifiers
            .get(
                "gene",
                [""],
            )[0]
            .upper()
        )


        if gene not in {
            "ORF1",
            "ORF2",
        }:
            continue


        translation = (
            feature.qualifiers
            .get(
                "translation",
                [""],
            )[0]
        )


        if not translation:
            raise ValueError(
                f"{gene} CDS has no translation "
                "in the GenBank record."
            )


        proteins[
            gene
        ] = translation


    if set(
        proteins.keys()
    ) != {
        "ORF1",
        "ORF2",
    }:

        raise ValueError(
            "Could not recover both ORF1 and ORF2 "
            "protein sequences from AF148856.1."
        )


    return proteins


def adney_original_activity_class(
    retro_t
):
    """
    Adney Supplemental Table 6 categories:
      Poor     <25% WT
      Reduced  25-80% WT
      High     >80% WT
    """

    value = float(
        retro_t
    )


    if value < 25:
        return "Poor (<25% WT)"

    if value <= 80:
        return "Reduced (25-80% WT)"

    return "High (>80% WT)"


def project_retro_t_context(
    retro_t
):
    """
    Project-specific descriptive labels.

    These are NOT Adney's original bins.
    """

    value = float(
        retro_t
    )


    if value < PROJECT_WT_ADJACENT_LOW:

        return (
            "Reduced relative to WT"
        )


    if value <= PROJECT_WT_ADJACENT_HIGH:

        return "WT-adjacent"


    return (
        "Elevated relative to WT"
    )


def position_region(
    orf,
    aa_position
):
    """
    Adney Table 8 ORF2 region framework.

    ORF1 is intentionally not subdivided in Script 06 because
    this stage is specifically connecting sequence differences
    to the Adney experimental windows.

    ORF2:
      EN        1-239
      DESERT 1  240-379
      Z         380-480
      unassigned 481-497
      RT        498-773
      DESERT 2  774-1275
    """

    position = int(
        aa_position
    )


    if orf == "ORF1":
        return "ORF1"


    if orf != "ORF2":
        return "UNKNOWN_ORF"


    if 1 <= position <= 239:
        return "EN"

    if 240 <= position <= 379:
        return "DESERT 1"

    if 380 <= position <= 480:
        return "Z"

    if 481 <= position <= 497:
        return "Unassigned 481-497"

    if 498 <= position <= 773:
        return "RT"

    if 774 <= position <= 1275:
        return "DESERT 2"


    return "Outside ORF2 protein"


def window_region_span(
    orf,
    start_aa,
    end_aa
):
    """
    Identify all Adney regions touched by an experimental window.
    This matters for windows that cross a region boundary.
    """

    regions = []

    for position in range(
        int(start_aa),
        int(end_aa) + 1,
    ):

        region = position_region(
            orf,
            position,
        )

        if region not in regions:
            regions.append(
                region
            )


    return ";".join(
        regions
    )


def validate_and_expand_adney_table(
    adney_df,
    reference_proteins,
):
    """
    Validate every Adney WT window against the AF148856.1
    protein sequence and expand the experimental windows into
    a one-row-per-amino-acid lookup table.

    This gives Script 06 a direct mapping:

        (ORF, AA_Position)
                ->
        Adney experimental window
    """

    # Use ALL Adney constructs.
    required = [
        "ORF",
        "Mutant_ID",
        "Start_AA",
        "End_AA",
        "WT_Residues",
        "RetroT_Average",
        "Standard_Deviation",
        "N_Measurements",
    ]


    adney = (
        adney_df[
            required
        ]
        .copy()
    )


    adney["ORF"] = (
        adney["ORF"]
        .astype(str)
        .str.upper()
    )


    adney["Start_AA"] = (
        pd.to_numeric(
            adney["Start_AA"],
            errors="raise",
        )
        .astype(int)
    )

    adney["End_AA"] = (
        pd.to_numeric(
            adney["End_AA"],
            errors="raise",
        )
        .astype(int)
    )


    # Expected library size from Supplemental Table 2.
    if len(adney) != 538:

        raise ValueError(
            "\nExpected 538 Adney constructs, "
            f"but found {len(adney):,}."
        )


    expected_counts = {
        "ORF1": 113,
        "ORF2": 425,
    }


    actual_counts = (
        adney["ORF"]
        .value_counts()
        .to_dict()
    )


    if actual_counts != expected_counts:

        raise ValueError(
            "\nUnexpected Adney ORF construct counts.\n"
            f"Expected: {expected_counts}\n"
            f"Observed: {actual_counts}"
        )


    expanded_rows = []

    validation_errors = []


    for _, row in adney.iterrows():

        orf = row[
            "ORF"
        ]

        start = int(
            row[
                "Start_AA"
            ]
        )

        end = int(
            row[
                "End_AA"
            ]
        )

        wt_residues = str(
            row[
                "WT_Residues"
            ]
        ).strip()


        if orf not in reference_proteins:

            validation_errors.append(
                (
                    row["Mutant_ID"],
                    "Unknown ORF",
                )
            )

            continue


        protein = (
            reference_proteins[
                orf
            ]
        )


        if (
            start < 1
            or end < start
            or end > len(protein)
        ):

            validation_errors.append(
                (
                    row["Mutant_ID"],
                    "Invalid AA coordinates",
                )
            )

            continue


        reference_window = (
            protein[
                start - 1:end
            ]
        )


        if (
            reference_window
            != wt_residues
        ):

            validation_errors.append(
                (
                    row["Mutant_ID"],
                    (
                        f"WT mismatch: "
                        f"Adney={wt_residues}, "
                        f"L1-RP={reference_window}"
                    ),
                )
            )

            continue


        window_size = (
            end
            - start
            + 1
        )


        if (
            len(wt_residues)
            != window_size
        ):

            validation_errors.append(
                (
                    row["Mutant_ID"],
                    "WT_Residues length mismatch",
                )
            )

            continue


        alanine_replacement = (
            "A"
            * window_size
        )


        original_class = (
            adney_original_activity_class(
                row[
                    "RetroT_Average"
                ]
            )
        )


        project_context = (
            project_retro_t_context(
                row[
                    "RetroT_Average"
                ]
            )
        )


        region_span = (
            window_region_span(
                orf,
                start,
                end,
            )
        )


        for position in range(
            start,
            end + 1,
        ):

            index_in_window = (
                position
                - start
            )


            expanded_rows.append(
                {
                    "ORF":
                        orf,

                    "AA_Position":
                        position,

                    "Adney_Mutant_ID":
                        row[
                            "Mutant_ID"
                        ],

                    "Adney_Window_Start_AA":
                        start,

                    "Adney_Window_End_AA":
                        end,

                    "Adney_Window_Size_AA":
                        window_size,

                    "Adney_Window_WT_Residues":
                        wt_residues,

                    "Adney_Alanine_Replacement":
                        alanine_replacement,

                    "Adney_Position_Index_In_Window":
                        index_in_window + 1,

                    "Adney_Position_WT_AA":
                        wt_residues[
                            index_in_window
                        ],

                    "Adney_RetroT_pct_WT":
                        float(
                            row[
                                "RetroT_Average"
                            ]
                        ),

                    "Adney_RetroT_SD":
                        float(
                            row[
                                "Standard_Deviation"
                            ]
                        ),

                    "Adney_N_Measurements":
                        int(
                            row[
                                "N_Measurements"
                            ]
                        ),

                    "Adney_Original_Activity_Class":
                        original_class,

                    "Project_RetroT_Context":
                        project_context,

                    "Adney_Position_Region":
                        position_region(
                            orf,
                            position,
                        ),

                    "Adney_Window_Region_Span":
                        region_span,

                    "Adney_Evidence_Scope":
                        (
                            "WINDOW_LEVEL_ALANINE_"
                            "SUBSTITUTION_EVIDENCE"
                        ),

                    "Natural_Variant_Activity_Measured":
                        False,
                }
            )


    if validation_errors:

        preview = "\n".join(
            f"  {mutant}: {message}"
            for mutant, message
            in validation_errors[:20]
        )

        raise ValueError(
            "\nAdney vs AF148856.1 validation failed.\n"
            f"Errors: {len(validation_errors):,}\n"
            f"{preview}"
        )


    position_map = pd.DataFrame(
        expanded_rows
    )


    # Every tested amino-acid position must map to exactly one
    # experimental window.
    duplicate_position_count = int(
        position_map.duplicated(
            subset=[
                "ORF",
                "AA_Position",
            ],
            keep=False,
        ).sum()
    )


    if duplicate_position_count > 0:

        raise ValueError(
            "\nAdney experimental windows overlap. "
            "A residue maps to more than one construct."
        )


    # The library covers every residue except the initiating
    # methionine (AA1) in ORF1 and ORF2.
    for orf, protein in (
        reference_proteins.items()
    ):

        covered = set(
            position_map.loc[
                position_map[
                    "ORF"
                ]
                == orf,
                "AA_Position",
            ]
            .astype(int)
        )


        expected = set(
            range(
                2,
                len(protein) + 1,
            )
        )


        if covered != expected:

            missing = sorted(
                expected
                - covered
            )

            extra = sorted(
                covered
                - expected
            )

            raise ValueError(
                f"\nUnexpected Adney positional coverage "
                f"for {orf}.\n"
                f"Missing: {missing[:20]}\n"
                f"Extra: {extra[:20]}"
            )


    return (
        adney,
        position_map,
    )


def mapping_status(
    row,
    protein_lengths,
):
    """
    Explain why an AA difference did or did not map to an Adney
    experimental window.
    """

    if pd.notna(
        row[
            "Adney_Mutant_ID"
        ]
    ):

        return "MAPPED"


    orf = row[
        "ORF"
    ]

    position = int(
        row[
            "AA_Position"
        ]
    )


    if (
        position == 1
    ):

        return (
            "NOT_TESTED_INITIATOR_AA1"
        )


    if (
        row[
            "Difference_Type"
        ]
        == "STOP_LOST"
    ):

        return (
            "TERMINAL_STOP_OUTSIDE_ADNEY_LIBRARY"
        )


    protein_length = (
        protein_lengths[
            orf
        ]
    )


    if (
        position
        > protein_length
    ):

        return (
            "OUTSIDE_REFERENCE_PROTEIN"
        )


    return (
        "NOT_COVERED_BY_ADNEY_LIBRARY"
    )


def count_equals(
    series,
    value
):
    return int(
        (
            series
            == value
        ).sum()
    )


def aggregate_locus_orf(
    mapped_df
):
    """
    One row per locus + ORF.
    """

    rows = []


    for (
        locus_id,
        orf,
    ), group in mapped_df.groupby(
        [
            "Locus_ID",
            "ORF",
        ],
        sort=False,
    ):

        mapped = group.loc[
            group[
                "Mapping_Status"
            ]
            == "MAPPED"
        ]


        retro_t = (
            mapped[
                "Adney_RetroT_pct_WT"
            ]
            .dropna()
        )


        rows.append(
            {
                "Locus_ID":
                    locus_id,

                "Primary_Subfamily":
                    group[
                        "Primary_Subfamily"
                    ].iloc[0],

                "ORF":
                    orf,

                "Script04_ORF_Status":
                    group[
                        "Script04_ORF_Status"
                    ].iloc[0],

                "Total_AA_Differences":
                    len(group),

                "Adney_Mapped_Differences":
                    len(mapped),

                "Adney_Unmapped_Differences":
                    len(group)
                    - len(mapped),

                "Unique_Adney_Windows_Impacted":
                    mapped[
                        "Adney_Mutant_ID"
                    ].nunique(),

                "Adney_Poor_Differences":
                    count_equals(
                        mapped[
                            "Adney_Original_Activity_Class"
                        ],
                        "Poor (<25% WT)",
                    ),

                "Adney_Reduced_25_80_Differences":
                    count_equals(
                        mapped[
                            "Adney_Original_Activity_Class"
                        ],
                        "Reduced (25-80% WT)",
                    ),

                "Adney_High_gt80_Differences":
                    count_equals(
                        mapped[
                            "Adney_Original_Activity_Class"
                        ],
                        "High (>80% WT)",
                    ),

                "Project_Reduced_Differences":
                    count_equals(
                        mapped[
                            "Project_RetroT_Context"
                        ],
                        "Reduced relative to WT",
                    ),

                "Project_WT_adjacent_Differences":
                    count_equals(
                        mapped[
                            "Project_RetroT_Context"
                        ],
                        "WT-adjacent",
                    ),

                "Project_Elevated_Differences":
                    count_equals(
                        mapped[
                            "Project_RetroT_Context"
                        ],
                        "Elevated relative to WT",
                    ),

                "Minimum_Adney_RetroT_pct_WT":
                    (
                        float(
                            retro_t.min()
                        )
                        if len(
                            retro_t
                        )
                        else np.nan
                    ),

                "Median_Adney_RetroT_pct_WT":
                    (
                        float(
                            retro_t.median()
                        )
                        if len(
                            retro_t
                        )
                        else np.nan
                    ),

                "Maximum_Adney_RetroT_pct_WT":
                    (
                        float(
                            retro_t.max()
                        )
                        if len(
                            retro_t
                        )
                        else np.nan
                    ),
            }
        )


    return pd.DataFrame(
        rows
    )


def aggregate_locus(
    mapped_df,
    integrity_df,
):
    """
    One row per genomic LINE-1 locus.

    Script 04 provides the base so loci with zero Script 05
    nonsynonymous AA differences are still retained.
    """

    base_columns = [
        "Locus_ID",
        "Primary_Subfamily",
        "Alignment_QC",
        "ORF1_Sequence_Status",
        "ORF2_Sequence_Status",
    ]


    base = (
        integrity_df[
            base_columns
        ]
        .drop_duplicates(
            subset=[
                "Locus_ID",
            ]
        )
        .copy()
    )


    if MODE == "TEST":

        selected = set(
            mapped_df[
                "Locus_ID"
            ]
        )

        base = base.loc[
            base[
                "Locus_ID"
            ]
            .isin(
                selected
            )
        ].copy()


    rows = []


    for locus_id, group in mapped_df.groupby(
        "Locus_ID",
        sort=False,
    ):

        mapped = group.loc[
            group[
                "Mapping_Status"
            ]
            == "MAPPED"
        ]


        retro_t = (
            mapped[
                "Adney_RetroT_pct_WT"
            ]
            .dropna()
        )


        rows.append(
            {
                "Locus_ID":
                    locus_id,

                "Total_AA_Differences":
                    len(
                        group
                    ),

                "Adney_Mapped_Differences":
                    len(
                        mapped
                    ),

                "Adney_Unmapped_Differences":
                    len(
                        group
                    )
                    - len(
                        mapped
                    ),

                "Unique_Adney_Windows_Impacted":
                    mapped[
                        "Adney_Mutant_ID"
                    ].nunique(),

                "Adney_Poor_Differences":
                    count_equals(
                        mapped[
                            "Adney_Original_Activity_Class"
                        ],
                        "Poor (<25% WT)",
                    ),

                "Adney_Reduced_25_80_Differences":
                    count_equals(
                        mapped[
                            "Adney_Original_Activity_Class"
                        ],
                        "Reduced (25-80% WT)",
                    ),

                "Adney_High_gt80_Differences":
                    count_equals(
                        mapped[
                            "Adney_Original_Activity_Class"
                        ],
                        "High (>80% WT)",
                    ),

                "Project_Reduced_Differences":
                    count_equals(
                        mapped[
                            "Project_RetroT_Context"
                        ],
                        "Reduced relative to WT",
                    ),

                "Project_WT_adjacent_Differences":
                    count_equals(
                        mapped[
                            "Project_RetroT_Context"
                        ],
                        "WT-adjacent",
                    ),

                "Project_Elevated_Differences":
                    count_equals(
                        mapped[
                            "Project_RetroT_Context"
                        ],
                        "Elevated relative to WT",
                    ),

                "Minimum_Adney_RetroT_pct_WT":
                    (
                        float(
                            retro_t.min()
                        )
                        if len(
                            retro_t
                        )
                        else np.nan
                    ),

                "Median_Adney_RetroT_pct_WT":
                    (
                        float(
                            retro_t.median()
                        )
                        if len(
                            retro_t
                        )
                        else np.nan
                    ),
            }
        )


    aggregate = pd.DataFrame(
        rows
    )


    result = base.merge(
        aggregate,
        on="Locus_ID",
        how="left",
        validate="one_to_one",
    )


    count_columns = [
        "Total_AA_Differences",
        "Adney_Mapped_Differences",
        "Adney_Unmapped_Differences",
        "Unique_Adney_Windows_Impacted",
        "Adney_Poor_Differences",
        "Adney_Reduced_25_80_Differences",
        "Adney_High_gt80_Differences",
        "Project_Reduced_Differences",
        "Project_WT_adjacent_Differences",
        "Project_Elevated_Differences",
    ]


    for column in count_columns:

        result[
            column
        ] = (
            result[
                column
            ]
            .fillna(
                0
            )
            .astype(
                int
            )
        )


    return result


# ============================================================
# STEP 1 — LOAD INPUT DATA
# ============================================================

print("\n============================================================")
print("SCRIPT 06 — MAP AA DIFFERENCES TO ADNEY WINDOWS")
print("============================================================")

print("\nRun mode:", MODE)

print("\nScript 05 AA differences:")
print(AA_DIFFERENCES_FILE)

print("\nScript 04 ORF integrity:")
print(ORF_INTEGRITY_FILE)

print("\nAdney workbook:")
print(ADNEY_WORKBOOK)

print("\nL1-RP reference:")
print(REFERENCE_GB)


aa_df = pd.read_csv(
    AA_DIFFERENCES_FILE
)


integrity_df = pd.read_csv(
    ORF_INTEGRITY_FILE
)


print(
    f"\nScript 05 AA-difference rows loaded: "
    f"{len(aa_df):,}"
)

print(
    f"Script 04 loci loaded: "
    f"{len(integrity_df):,}"
)


# ============================================================
# STEP 2 — LOAD + VALIDATE ADNEY TABLE
# ============================================================

print("\n============================================================")
print("STEP 2 — Validating Adney Table 2 against AF148856.1")
print("============================================================")


adney_raw, adney_sheet = (
    load_adney_table(
        ADNEY_WORKBOOK
    )
)


reference_proteins = (
    extract_reference_proteins(
        REFERENCE_GB
    )
)


adney_windows, position_map = (
    validate_and_expand_adney_table(
        adney_raw,
        reference_proteins,
    )
)


protein_lengths = {
    orf: len(
        sequence
    )
    for orf, sequence
    in reference_proteins.items()
}


print(
    f"Adney worksheet used: "
    f"{adney_sheet}"
)

print(
    f"Adney constructs validated: "
    f"{len(adney_windows):,}"
)

print(
    "  ORF1:",
    int(
        (
            adney_windows[
                "ORF"
            ]
            == "ORF1"
        ).sum()
    )
)

print(
    "  ORF2:",
    int(
        (
            adney_windows[
                "ORF"
            ]
            == "ORF2"
        ).sum()
    )
)

print(
    f"Expanded tested AA positions: "
    f"{len(position_map):,}"
)

print(
    "\nAA1 in ORF1 and ORF2 is intentionally not "
    "covered by the Adney library."
)


# ============================================================
# STEP 3 — SELECT TEST / FULL AA DIFFERENCES
# ============================================================

print("\n============================================================")
print("STEP 3 — Selecting AA differences")
print("============================================================")


if MODE == "FULL":

    work_df = aa_df.copy()


else:

    available_loci = list(
        aa_df[
            "Locus_ID"
        ]
        .drop_duplicates()
    )


    available_set = set(
        available_loci
    )


    selected_loci = []


    for locus_id in TEST_LOCUS_IDS:

        if (
            locus_id
            in available_set
        ):

            selected_loci.append(
                locus_id
            )


    remaining = [
        locus_id
        for locus_id
        in available_loci
        if locus_id
        not in selected_loci
    ]


    rng = random.Random(
        RANDOM_SEED
    )


    additional_needed = max(
        0,
        TEST_N_LOCI
        - len(
            selected_loci
        )
    )


    if (
        additional_needed
        > len(
            remaining
        )
    ):

        additional_needed = len(
            remaining
        )


    selected_loci.extend(
        rng.sample(
            remaining,
            additional_needed,
        )
    )


    work_df = aa_df.loc[
        aa_df[
            "Locus_ID"
        ]
        .isin(
            selected_loci
        )
    ].copy()


    print(
        f"TEST loci selected: "
        f"{len(selected_loci):,}"
    )

    print(
        f"TEST AA-difference rows: "
        f"{len(work_df):,}"
    )


# ============================================================
# STEP 4 — MAP AA POSITIONS TO ADNEY WINDOWS
# ============================================================

print("\n============================================================")
print("STEP 4 — Mapping AA differences to Adney windows")
print("============================================================")


mapped_df = work_df.merge(
    position_map,
    on=[
        "ORF",
        "AA_Position",
    ],
    how="left",
    validate="many_to_one",
)


mapped_df[
    "Mapping_Status"
] = mapped_df.apply(
    lambda row:
        mapping_status(
            row,
            protein_lengths,
        ),
    axis=1,
)


mapped_df[
    "Adney_Mapped"
] = (
    mapped_df[
        "Mapping_Status"
    ]
    == "MAPPED"
)


# ============================================================
# STEP 5 — CROSS-CHECK WT AA NUMBERING
# ============================================================

print("\n============================================================")
print("STEP 5 — Cross-checking Script 05 AA numbering")
print("============================================================")


mapped_rows = mapped_df.loc[
    mapped_df[
        "Adney_Mapped"
    ]
].copy()


wt_mismatch_mask = (
    mapped_rows[
        "L1RP_AA"
    ]
    != mapped_rows[
        "Adney_Position_WT_AA"
    ]
)


wt_mismatch_count = int(
    wt_mismatch_mask.sum()
)


if wt_mismatch_count > 0:

    mismatch_preview = (
        mapped_rows.loc[
            wt_mismatch_mask,
            [
                "Locus_ID",
                "ORF",
                "AA_Position",
                "L1RP_AA",
                "Adney_Position_WT_AA",
                "Adney_Mutant_ID",
            ],
        ]
        .head(
            20
        )
        .to_string(
            index=False
        )
    )


    raise ValueError(
        "\nScript 05 / Adney AA numbering mismatch detected.\n"
        f"Rows affected: {wt_mismatch_count:,}\n\n"
        f"{mismatch_preview}"
    )


print(
    f"Mapped rows with L1-RP WT AA agreement: "
    f"{len(mapped_rows):,} / {len(mapped_rows):,}"
)


# ============================================================
# STEP 6 — CREATE PRIORITY CLEAN-ORF TABLE
# ============================================================

print("\n============================================================")
print("STEP 6 — Creating clean-ORF Adney priority subset")
print("============================================================")


priority_df = mapped_df.loc[
    (
        mapped_df[
            "Adney_Mapped"
        ]
    )
    &
    (
        mapped_df[
            "Script04_ORF_Status"
        ]
        ==
        "FULL_SPAN_NO_FRAME_STOP_DISRUPTION"
    )
].copy()


print(
    f"Priority mapped differences: "
    f"{len(priority_df):,}"
)


# ============================================================
# STEP 7 — BUILD SUMMARY TABLES
# ============================================================

print("\n============================================================")
print("STEP 7 — Building summaries")
print("============================================================")


locus_orf_summary = (
    aggregate_locus_orf(
        mapped_df
    )
)


locus_summary = (
    aggregate_locus(
        mapped_df,
        integrity_df,
    )
)


# One row per Adney experimental construct.
window_summary = (
    mapped_rows
    .groupby(
        [
            "ORF",
            "Adney_Mutant_ID",
            "Adney_Window_Start_AA",
            "Adney_Window_End_AA",
            "Adney_Window_WT_Residues",
            "Adney_Alanine_Replacement",
            "Adney_RetroT_pct_WT",
            "Adney_RetroT_SD",
            "Adney_N_Measurements",
            "Adney_Original_Activity_Class",
            "Project_RetroT_Context",
            "Adney_Window_Region_Span",
        ],
        dropna=False,
    )
    .agg(
        Total_Natural_AA_Difference_Rows=(
            "Locus_ID",
            "count",
        ),

        Unique_Loci_With_Difference=(
            "Locus_ID",
            "nunique",
        ),

        Unique_AA_Positions_Altered=(
            "AA_Position",
            "nunique",
        ),

        Unique_Observed_Locus_AAs=(
            "Locus_AA",
            "nunique",
        ),

        Clean_ORF_Difference_Rows=(
            "Script04_ORF_Status",
            lambda x:
                int(
                    (
                        x
                        ==
                        "FULL_SPAN_NO_FRAME_STOP_DISRUPTION"
                    ).sum()
                ),
        ),
    )
    .reset_index()
)


# Add Adney constructs that had zero observed differences in the
# selected genomic data.
window_base = (
    adney_windows[
        [
            "ORF",
            "Mutant_ID",
            "Start_AA",
            "End_AA",
            "WT_Residues",
            "RetroT_Average",
            "Standard_Deviation",
            "N_Measurements",
        ]
    ]
    .rename(
        columns={
            "Mutant_ID":
                "Adney_Mutant_ID",

            "Start_AA":
                "Adney_Window_Start_AA",

            "End_AA":
                "Adney_Window_End_AA",

            "WT_Residues":
                "Adney_Window_WT_Residues",

            "RetroT_Average":
                "Adney_RetroT_pct_WT",

            "Standard_Deviation":
                "Adney_RetroT_SD",

            "N_Measurements":
                "Adney_N_Measurements",
        }
    )
    .copy()
)


window_base[
    "Adney_Alanine_Replacement"
] = window_base.apply(
    lambda row:
        "A"
        * (
            int(
                row[
                    "Adney_Window_End_AA"
                ]
            )
            - int(
                row[
                    "Adney_Window_Start_AA"
                ]
            )
            + 1
        ),
    axis=1,
)


window_base[
    "Adney_Original_Activity_Class"
] = window_base[
    "Adney_RetroT_pct_WT"
].apply(
    adney_original_activity_class
)


window_base[
    "Project_RetroT_Context"
] = window_base[
    "Adney_RetroT_pct_WT"
].apply(
    project_retro_t_context
)


window_base[
    "Adney_Window_Region_Span"
] = window_base.apply(
    lambda row:
        window_region_span(
            row[
                "ORF"
            ],
            row[
                "Adney_Window_Start_AA"
            ],
            row[
                "Adney_Window_End_AA"
            ],
        ),
    axis=1,
)


window_merge_keys = [
    "ORF",
    "Adney_Mutant_ID",
    "Adney_Window_Start_AA",
    "Adney_Window_End_AA",
    "Adney_Window_WT_Residues",
    "Adney_Alanine_Replacement",
    "Adney_RetroT_pct_WT",
    "Adney_RetroT_SD",
    "Adney_N_Measurements",
    "Adney_Original_Activity_Class",
    "Project_RetroT_Context",
    "Adney_Window_Region_Span",
]


window_summary = (
    window_base
    .merge(
        window_summary,
        on=window_merge_keys,
        how="left",
        validate="one_to_one",
    )
)


window_count_columns = [
    "Total_Natural_AA_Difference_Rows",
    "Unique_Loci_With_Difference",
    "Unique_AA_Positions_Altered",
    "Unique_Observed_Locus_AAs",
    "Clean_ORF_Difference_Rows",
]


for column in window_count_columns:

    window_summary[
        column
    ] = (
        window_summary[
            column
        ]
        .fillna(
            0
        )
        .astype(
            int
        )
    )


# ============================================================
# STEP 8 — SAVE FULL TABLES
# ============================================================

print("\n============================================================")
print("STEP 8 — Saving outputs")
print("============================================================")


mapped_df.to_csv(
    MAPPED_AA_CSV_GZ,
    index=False,
    compression="gzip",
)


priority_df.to_csv(
    PRIORITY_CSV_GZ,
    index=False,
    compression="gzip",
)


locus_orf_summary.to_csv(
    LOCUS_ORF_SUMMARY_CSV_GZ,
    index=False,
    compression="gzip",
)


locus_summary.to_csv(
    LOCUS_SUMMARY_CSV_GZ,
    index=False,
    compression="gzip",
)


window_summary.to_csv(
    WINDOW_SUMMARY_CSV,
    index=False,
)


if MODE == "TEST":

    mapped_df.to_csv(
        MAPPED_AA_CSV,
        index=False,
    )

    priority_df.to_csv(
        PRIORITY_CSV,
        index=False,
    )


# ============================================================
# STEP 9 — SUMMARY METRICS
# ============================================================

mapping_status_counts = (
    mapped_df[
        "Mapping_Status"
    ]
    .value_counts(
        dropna=False
    )
    .rename_axis(
        "Mapping_Status"
    )
    .reset_index(
        name="AA_Difference_Rows"
    )
)


adney_class_counts = (
    mapped_rows[
        "Adney_Original_Activity_Class"
    ]
    .value_counts()
    .rename_axis(
        "Adney_Original_Activity_Class"
    )
    .reset_index(
        name="Mapped_AA_Difference_Rows"
    )
)


project_context_counts = (
    mapped_rows[
        "Project_RetroT_Context"
    ]
    .value_counts()
    .rename_axis(
        "Project_RetroT_Context"
    )
    .reset_index(
        name="Mapped_AA_Difference_Rows"
    )
)


priority_adney_class_counts = (
    priority_df[
        "Adney_Original_Activity_Class"
    ]
    .value_counts()
    .rename_axis(
        "Adney_Original_Activity_Class"
    )
    .reset_index(
        name="Priority_AA_Difference_Rows"
    )
)


priority_project_context_counts = (
    priority_df[
        "Project_RetroT_Context"
    ]
    .value_counts()
    .rename_axis(
        "Project_RetroT_Context"
    )
    .reset_index(
        name="Priority_AA_Difference_Rows"
    )
)


summary_rows = [
    [
        "Run mode",
        MODE,
    ],
    [
        "Script 05 AA-difference rows analyzed",
        len(
            mapped_df
        ),
    ],
    [
        "AA differences mapped to an Adney window",
        int(
            mapped_df[
                "Adney_Mapped"
            ].sum()
        ),
    ],
    [
        "AA differences not mapped",
        int(
            (
                ~mapped_df[
                    "Adney_Mapped"
                ]
            ).sum()
        ),
    ],
    [
        "Adney constructs used",
        len(
            adney_windows
        ),
    ],
    [
        "Adney ORF1 constructs",
        int(
            (
                adney_windows[
                    "ORF"
                ]
                == "ORF1"
            ).sum()
        ),
    ],
    [
        "Adney ORF2 constructs",
        int(
            (
                adney_windows[
                    "ORF"
                ]
                == "ORF2"
            ).sum()
        ),
    ],
    [
        "Mapped rows with WT AA numbering mismatch",
        wt_mismatch_count,
    ],
    [
        "Priority clean-ORF mapped differences",
        len(
            priority_df
        ),
    ],
    [
        "Priority subset missense rows",
        int(
            (
                priority_df[
                    "Difference_Type"
                ]
                == "MISSENSE"
            ).sum()
        ),
    ],
    [
        "Adney evidence scope",
        (
            "Window-level alanine-substitution evidence; "
            "not measured natural-variant activity"
        ),
    ],
    [
        "Project WT-adjacent range",
        (
            f"{PROJECT_WT_ADJACENT_LOW:g}-"
            f"{PROJECT_WT_ADJACENT_HIGH:g}% WT inclusive"
        ),
    ],
]


summary_df = pd.DataFrame(
    summary_rows,
    columns=[
        "Metric",
        "Value",
    ]
)


# ============================================================
# STEP 10 — EXCEL SUMMARY
# ============================================================

if WRITE_EXCEL_SUMMARY:

    with pd.ExcelWriter(
        SUMMARY_XLSX,
        engine="openpyxl",
    ) as writer:

        summary_df.to_excel(
            writer,
            sheet_name="Summary",
            index=False,
        )

        mapping_status_counts.to_excel(
            writer,
            sheet_name="Mapping_Status",
            index=False,
        )

        adney_class_counts.to_excel(
            writer,
            sheet_name="Adney_Original_Bins",
            index=False,
        )

        project_context_counts.to_excel(
            writer,
            sheet_name="Project_RetroT_Context",
            index=False,
        )

        priority_adney_class_counts.to_excel(
            writer,
            sheet_name="Priority_Adney_Bins",
            index=False,
        )

        priority_project_context_counts.to_excel(
            writer,
            sheet_name="Priority_Project_Context",
            index=False,
        )

        adney_windows.to_excel(
            writer,
            sheet_name="Adney_538_Windows",
            index=False,
        )

        position_map.to_excel(
            writer,
            sheet_name="Adney_Position_Map",
            index=False,
        )

        window_summary.to_excel(
            writer,
            sheet_name="Window_Summary",
            index=False,
        )

        locus_orf_summary.to_excel(
            writer,
            sheet_name="Locus_ORF_Summary",
            index=False,
        )

        locus_summary.to_excel(
            writer,
            sheet_name="Locus_Summary",
            index=False,
        )

        if MODE == "TEST":

            mapped_df.head(
                50000
            ).to_excel(
                writer,
                sheet_name="TEST_Mapped_Details",
                index=False,
            )

            priority_df.head(
                50000
            ).to_excel(
                writer,
                sheet_name="TEST_Priority",
                index=False,
            )


# ============================================================
# FINAL REPORT
# ============================================================

print("\n============================================================")
print("SCRIPT 06 COMPLETE")
print("============================================================")

print(
    f"\nAA-difference rows analyzed: "
    f"{len(mapped_df):,}"
)

print(
    f"Mapped to Adney windows: "
    f"{int(mapped_df['Adney_Mapped'].sum()):,}"
)

print(
    f"Not mapped: "
    f"{int((~mapped_df['Adney_Mapped']).sum()):,}"
)

print(
    f"Priority clean-ORF mapped differences: "
    f"{len(priority_df):,}"
)

print(
    "\nMapping status:"
)

print(
    mapping_status_counts.to_string(
        index=False
    )
)

print(
    "\nIMPORTANT:"
)

print(
    "Adney RetroT values describe the alanine-substitution "
    "experimental WINDOW."
)

print(
    "They are not the measured activity of the natural "
    "genomic amino-acid difference."
)

print(
    "\nOutputs:"
)

print(
    "  Detailed mapped AA differences:",
    MAPPED_AA_CSV_GZ,
)

print(
    "  Priority clean-ORF subset:",
    PRIORITY_CSV_GZ,
)

print(
    "  Locus + ORF summary:",
    LOCUS_ORF_SUMMARY_CSV_GZ,
)

print(
    "  Locus summary:",
    LOCUS_SUMMARY_CSV_GZ,
)

print(
    "  Adney-window summary:",
    WINDOW_SUMMARY_CSV,
)

if WRITE_EXCEL_SUMMARY:

    print(
        "  Summary workbook:",
        SUMMARY_XLSX,
    )


if MODE == "TEST":

    print(
        "\nThis was a TEST run."
    )

    print(
        "After QC, change:"
    )

    print(
        '    RUN_MODE = "TEST"'
    )

    print(
        "to:"
    )

    print(
        '    RUN_MODE = "FULL"'
    )
