# ============================================================
# LINE-1 FUNCTIONAL ANNOTATION PROJECT
#
# SCRIPT 08 V2:
# Final locus -> observed sequence event -> predicted functional
# consequence summary
#
# PURPOSE
# ------------------------------------------------------------
# This script combines the major evidence layers produced by
# Scripts 03-07 into two final products:
#
#   1. Mutation/Event consequence table
#      One row per observed amino-acid difference or ORF indel.
#
#   2. Locus summary table
#      One row per genomic LINE-1 locus.
#
# The goal is to directly answer:
#
#   "Which genomic LINE-1 locus has which observed sequence
#    change, and what is the predicted functional consequence?"
#
# IMPORTANT SCIENTIFIC LANGUAGE
# ------------------------------------------------------------
# Older LINE-1 subfamilies can differ substantially from L1-RP.
# Therefore an observed difference vs L1-RP is NOT automatically
# a newly acquired locus-specific mutation.
#
# The output uses:
#
#   Observed_Event_vs_L1RP
#
# rather than assuming mutational origin.
#
# "Predicted_Functional_Consequence" is evidence-based and
# conservative. It describes likely coding consequences and
# experimental context; it does NOT claim measured activity of
# the genomic locus.
#
# Adney RetroT values remain window-level experimental evidence
# from multi-residue alanine constructs, not measured activity
# of individual natural genomic substitutions.
# ============================================================


# ============================================================
# IMPORTS
# ============================================================

from pathlib import Path

import numpy as np
import pandas as pd


# ============================================================
# SETTINGS
# ============================================================

# Start with TEST. After QC, switch to FULL.
RUN_MODE = "FULL"

WRITE_EXCEL_WORKBOOK = True

# In FULL mode the event table can approach ~400,000 rows.
# Excel can hold this, but writing the sheet can take time and
# create a large workbook. Leave True if you want one large
# human-readable workbook in addition to the compressed CSVs.
WRITE_EVENT_TABLE_TO_EXCEL = True


# ============================================================
# PROJECT PATHS
# ============================================================

SCRIPT_FOLDER = Path(__file__).resolve().parent

SCRIPT03_FOLDER = (
    SCRIPT_FOLDER
    / "03_LINE1_structure"
)

SCRIPT04_FOLDER = (
    SCRIPT_FOLDER
    / "04_ORF_integrity"
)

SCRIPT05_FOLDER = (
    SCRIPT_FOLDER
    / "05_ORF_variant_alignment"
)

SCRIPT07_FOLDER = (
    SCRIPT_FOLDER
    / "07_substitution_context"
)

OUTPUT_FOLDER = (
    SCRIPT_FOLDER
    / "08_final_annotation"
)

OUTPUT_FOLDER.mkdir(
    parents=True,
    exist_ok=True,
)


# ============================================================
# HELPERS
# ============================================================

def first_existing(
    paths,
    description,
):

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


MODE = RUN_MODE.upper().strip()

if MODE not in {
    "TEST",
    "FULL",
}:

    raise ValueError(
        "RUN_MODE must be either 'TEST' or 'FULL'."
    )


# ============================================================
# INPUT FILES
# ============================================================

STRUCTURE_FILE = first_existing(
    [
        SCRIPT03_FOLDER
        / "FULL_LINE1_structure_alignment_results.csv.gz",

        SCRIPT03_FOLDER
        / "FULL_LINE1_structure_alignment_results.csv",

        SCRIPT_FOLDER
        / "FULL_LINE1_structure_alignment_results.csv.gz",

        SCRIPT_FOLDER
        / "FULL_LINE1_structure_alignment_results.csv",
    ],
    "Script 03 FULL structure table",
)


INTEGRITY_FILE = first_existing(
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


INDEL_FILE = first_existing(
    [
        SCRIPT05_FOLDER
        / "FULL_LINE1_ORF_indel_events_vs_L1RP.csv.gz",

        SCRIPT_FOLDER
        / "FULL_LINE1_ORF_indel_events_vs_L1RP.csv.gz",
    ],
    "Script 05 FULL ORF-indel table",
)


if MODE == "FULL":

    SUBSTITUTION_FILE = first_existing(
        [
            SCRIPT07_FOLDER
            / "FULL_LINE1_substitution_context.csv.gz",

            SCRIPT_FOLDER
            / "FULL_LINE1_substitution_context.csv.gz",
        ],
        "Script 07 FULL substitution-context table",
    )

else:

    SUBSTITUTION_FILE = first_existing(
        [
            SCRIPT07_FOLDER
            / "TEST_LINE1_substitution_context.csv.gz",

            SCRIPT07_FOLDER
            / "TEST_LINE1_substitution_context.csv",

            SCRIPT_FOLDER
            / "TEST_LINE1_substitution_context.csv.gz",

            SCRIPT_FOLDER
            / "TEST_LINE1_substitution_context.csv",
        ],
        "Script 07 TEST substitution-context table",
    )


# ============================================================
# OUTPUT FILES
# ============================================================

EVENT_OUTPUT = (
    OUTPUT_FOLDER
    / f"{MODE}_LINE1_final_event_consequence_table_V2.csv.gz"
)

LOCUS_OUTPUT = (
    OUTPUT_FOLDER
    / f"{MODE}_LINE1_final_locus_summary_V2.csv.gz"
)

WORKBOOK_OUTPUT = (
    OUTPUT_FOLDER
    / f"{MODE}_LINE1_final_annotation_summary_V2.xlsx"
)


if MODE == "TEST":

    EVENT_OUTPUT_PLAIN = (
        OUTPUT_FOLDER
        / "TEST_LINE1_final_event_consequence_table_V2.csv"
    )

    LOCUS_OUTPUT_PLAIN = (
        OUTPUT_FOLDER
        / "TEST_LINE1_final_locus_summary_V2.csv"
    )


# ============================================================
# POSITION / REGION HELPERS
# ============================================================

def approximate_region(
    orf,
    aa_position,
):

    if pd.isna(
        aa_position
    ):
        return np.nan


    try:
        position = int(
            float(
                aa_position
            )
        )
    except Exception:
        return np.nan


    if orf == "ORF1":

        if 1 <= position <= 51:
            return "ORF1_NTR"

        if 52 <= position <= 153:
            return "ORF1_CC"

        if 157 <= position <= 252:
            return "ORF1_RRM"

        if 264 <= position <= 323:
            return "ORF1_CTD"

        return "ORF1_unassigned_or_interdomain"


    if orf == "ORF2":

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

        return "ORF2_outside_reference_protein"


    return np.nan


# ============================================================
# PREDICTED CONSEQUENCE LOGIC
# ============================================================

def aa_predicted_consequence(
    row
):
    """
    Conservative predicted consequence for an amino-acid event.

    This intentionally avoids statements such as:
        "functional"
        "nonfunctional"
        "X% activity"

    because those are not directly established by these data.
    """

    difference_type = str(
        row[
            "Difference_Type"
        ]
    )


    if difference_type == "START_LOST":

        return (
            "Likely severe coding disruption: "
            "initiator/start codon is lost."
        )


    if difference_type == "NONSENSE":

        return (
            "Likely severe coding disruption: "
            "premature termination codon introduced."
        )


    if difference_type == "STOP_LOST":

        return (
            "Altered termination predicted: "
            "reference terminal stop codon is lost."
        )


    if difference_type != "MISSENSE":

        return (
            "Sequence consequence identified, but no "
            "missense-specific functional interpretation "
            "is assigned."
        )


    mapping_status = str(
        row[
            "Mapping_Status"
        ]
    )


    if mapping_status != "MAPPED":

        return (
            "Missense substitution without direct Adney "
            "position-level experimental context; "
            "functional effect remains uncertain."
        )


    comparability = str(
        row[
            "Position_Level_Adney_Comparability"
        ]
    )


    adney_class = str(
        row[
            "Adney_Original_Activity_Class"
        ]
    )


    retro_t = row[
        "Adney_RetroT_pct_WT"
    ]


    if pd.isna(
        retro_t
    ):
        retro_text = "unknown RetroT"

    else:
        retro_text = (
            f"{float(retro_t):.1f}% WT RetroT"
        )


    if (
        comparability
        ==
        "ADNEY_DID_NOT_CHANGE_WT_ALANINE_POSITION"
    ):

        return (
            "Missense at a WT-alanine position. "
            "Adney did not alter this specific residue, "
            f"so the {retro_text} window phenotype cannot "
            "be attributed to this position."
        )


    if (
        comparability
        ==
        "NATURAL_EQUALS_ADNEY_ALANINE_AT_THIS_POSITION"
    ):

        return (
            "Natural missense matches the alanine residue "
            "used by Adney at this position within a "
            f"{adney_class} window ({retro_text}). "
            "This is the closest residue-level correspondence "
            "to the tested substitution, but Adney altered "
            "the full multi-residue window."
        )


    direction = str(
        row[
            "Metric_Directional_Agreement"
        ]
    )


    if (
        direction
        ==
        "CONCORDANT_NATURAL_MORE_SIMILAR_TO_WT_THAN_ALANINE"
    ):

        metric_text = (
            "Grantham and BLOSUM62 both place the natural "
            "substitution closer to WT than alanine."
        )

    elif (
        direction
        ==
        "CONCORDANT_NATURAL_LESS_SIMILAR_TO_WT_THAN_ALANINE"
    ):

        metric_text = (
            "Grantham and BLOSUM62 both place the natural "
            "substitution farther from WT than alanine."
        )

    elif (
        direction
        ==
        "BOTH_METRICS_TIED_WITH_ALANINE"
    ):

        metric_text = (
            "Grantham and BLOSUM62 score the natural "
            "substitution equivalently to alanine."
        )

    else:

        metric_text = (
            "Grantham and BLOSUM62 provide mixed or "
            "partially tied context relative to alanine."
        )


    return (
        "Missense substitution in an Adney "
        f"{adney_class} window ({retro_text}). "
        f"{metric_text} "
        "The functional effect of the natural substitution "
        "remains uncertain because Adney measured the full "
        "alanine-substituted window rather than this single "
        "natural change."
    )


def aa_consequence_class(
    row
):

    difference_type = str(
        row[
            "Difference_Type"
        ]
    )


    if difference_type == "START_LOST":
        return "LIKELY_SEVERE_CODING_DISRUPTION_START_LOST"

    if difference_type == "NONSENSE":
        return "LIKELY_SEVERE_CODING_DISRUPTION_NONSENSE"

    if difference_type == "STOP_LOST":
        return "ALTERED_TERMINATION_STOP_LOST"

    if difference_type != "MISSENSE":
        return "OTHER_AA_SEQUENCE_EVENT"


    if (
        str(
            row[
                "Position_Level_Adney_Comparability"
            ]
        )
        ==
        "ADNEY_DID_NOT_CHANGE_WT_ALANINE_POSITION"
    ):
        return "MISSENSE_ADNEY_POSITION_NOT_DIRECTLY_PERTURBED"


    if (
        str(
            row[
                "Position_Level_Adney_Comparability"
            ]
        )
        ==
        "NATURAL_EQUALS_ADNEY_ALANINE_AT_THIS_POSITION"
    ):
        return "MISSENSE_MATCHES_ADNEY_ALANINE_AT_POSITION"


    if (
        str(
            row[
                "Mapping_Status"
            ]
        )
        == "MAPPED"
    ):
        return "MISSENSE_WITH_ADNEY_WINDOW_CONTEXT"


    return "MISSENSE_WITHOUT_ADNEY_POSITION_CONTEXT"


def indel_predicted_consequence(
    row
):

    frame_effect = str(
        row[
            "Frame_Effect"
        ]
    )

    indel_type = str(
        row[
            "Indel_Type"
        ]
    ).lower()

    length_bp = int(
        row[
            "Length_bp"
        ]
    )


    if frame_effect == "FRAMESHIFT":

        return (
            f"Likely severe coding disruption: "
            f"{length_bp}-bp {indel_type} causes a "
            "frameshift and is expected to alter downstream "
            "protein coding sequence."
        )


    if frame_effect == "IN_FRAME":

        return (
            f"In-frame protein alteration: "
            f"{length_bp}-bp {indel_type} changes the "
            "protein sequence without shifting the downstream "
            "reading frame. Functional effect remains uncertain."
        )


    return (
        f"ORF indel event ({length_bp}-bp {indel_type}); "
        "frame consequence was not classified."
    )


def indel_consequence_class(
    row
):

    frame_effect = str(
        row[
            "Frame_Effect"
        ]
    )


    if frame_effect == "FRAMESHIFT":
        return "LIKELY_SEVERE_CODING_DISRUPTION_FRAMESHIFT"

    if frame_effect == "IN_FRAME":
        return "IN_FRAME_PROTEIN_ALTERATION"

    return "OTHER_INDEL_EVENT"


# ============================================================
# LOAD INPUT DATA
# ============================================================

print("\n============================================================")
print("SCRIPT 08 V2 — FINAL LINE-1 ANNOTATION SUMMARY")
print("============================================================")

print("\nRun mode:", MODE)

print("\nLoading Script 03...")
structure_df = pd.read_csv(
    STRUCTURE_FILE
)

print("Loading Script 04...")
integrity_df = pd.read_csv(
    INTEGRITY_FILE
)

print("Loading Script 05 indels...")
indel_df = pd.read_csv(
    INDEL_FILE
)

print("Loading Script 07 substitution context...")
aa_df = pd.read_csv(
    SUBSTITUTION_FILE
)


print(
    f"\nStructure loci: "
    f"{len(structure_df):,}"
)

print(
    f"Integrity loci: "
    f"{len(integrity_df):,}"
)

print(
    f"Indel events: "
    f"{len(indel_df):,}"
)

print(
    f"AA-difference rows: "
    f"{len(aa_df):,}"
)


# ============================================================
# TEST MODE FILTERING
# ============================================================

if MODE == "TEST":

    test_loci = set(
        aa_df[
            "Locus_ID"
        ]
    )


    structure_df = structure_df.loc[
        structure_df[
            "Locus_ID"
        ]
        .isin(
            test_loci
        )
    ].copy()


    integrity_df = integrity_df.loc[
        integrity_df[
            "Locus_ID"
        ]
        .isin(
            test_loci
        )
    ].copy()


    indel_df = indel_df.loc[
        indel_df[
            "Locus_ID"
        ]
        .isin(
            test_loci
        )
    ].copy()


    print(
        f"\nTEST loci carried forward: "
        f"{len(test_loci):,}"
    )


# ============================================================
# LOCUS METADATA / STRUCTURAL BACKBONE
# ============================================================

structure_keep = [
    "Locus_ID",
    "Chromosome",
    "Start_1based",
    "End_1based",
    "Strand",
    "Coordinate_Length",
    "Primary_Subfamily",
    "Contained_Within_Other_L1",
    "Contains_Other_L1",
    "Query_Alignment_Span_pct",
    "Alignment_Identity_pct",
    "ORF1_Coverage_pct",
    "ORF2_Coverage_pct",
    "ORF1_NTR_Coverage_pct",
    "ORF1_CC_Coverage_pct",
    "ORF1_RRM_Coverage_pct",
    "ORF1_CTD_Coverage_pct",
    "Adney_ORF2_EN_Coverage_pct",
    "Adney_ORF2_DESERT1_Coverage_pct",
    "Adney_ORF2_Z_Coverage_pct",
    "Adney_ORF2_Unassigned_481_497_Coverage_pct",
    "Adney_ORF2_RT_Coverage_pct",
    "Adney_ORF2_DESERT2_Coverage_pct",
]


structure_base = (
    structure_df[
        structure_keep
    ]
    .drop_duplicates(
        subset=[
            "Locus_ID",
        ]
    )
)


integrity_keep = [
    "Locus_ID",
    "Alignment_QC",
    "Alignment_QC_Reasons",
    "ORF1_Sequence_Status",
    "ORF1_Disruption_Reasons",
    "ORF1_Frameshifting_Indel_Events",
    "ORF1_Premature_Stop_Count",
    "ORF1_First_Premature_Stop_AA",
    "ORF2_Sequence_Status",
    "ORF2_Disruption_Reasons",
    "ORF2_Frameshifting_Indel_Events",
    "ORF2_Premature_Stop_Count",
    "ORF2_First_Premature_Stop_AA",
]


integrity_base = (
    integrity_df[
        integrity_keep
    ]
    .drop_duplicates(
        subset=[
            "Locus_ID",
        ]
    )
)


locus_base = structure_base.merge(
    integrity_base,
    on="Locus_ID",
    how="left",
    validate="one_to_one",
)


# ============================================================
# BUILD AA EVENT TABLE
# ============================================================

print("\nBuilding amino-acid event rows...")


aa_events = pd.DataFrame()


aa_events[
    "Locus_ID"
] = aa_df[
    "Locus_ID"
]

aa_events[
    "Chromosome"
] = aa_df[
    "Chromosome"
]

aa_events[
    "Start_1based"
] = aa_df[
    "Start_1based"
]

aa_events[
    "End_1based"
] = aa_df[
    "End_1based"
]

aa_events[
    "Strand"
] = aa_df[
    "Strand"
]

aa_events[
    "Primary_Subfamily"
] = aa_df[
    "Primary_Subfamily"
]

aa_events[
    "ORF"
] = aa_df[
    "ORF"
]

aa_events[
    "Event_Type"
] = "AA_DIFFERENCE"

aa_events[
    "Observed_Event_vs_L1RP"
] = aa_df[
    "AA_Difference_vs_L1RP"
]

aa_events[
    "Sequence_Consequence"
] = aa_df[
    "Difference_Type"
]

aa_events[
    "AA_Position"
] = aa_df[
    "AA_Position"
]

aa_events[
    "Reference_AA"
] = aa_df[
    "L1RP_AA"
]

aa_events[
    "Observed_AA"
] = aa_df[
    "Locus_AA"
]

aa_events[
    "Reference_Codon"
] = aa_df[
    "L1RP_Codon"
]

aa_events[
    "Observed_Codon"
] = aa_df[
    "Locus_Codon"
]

aa_events[
    "Nucleotide_Differences"
] = aa_df[
    "Nucleotide_Differences_In_Codon"
]

aa_events[
    "Indel_Type"
] = np.nan

aa_events[
    "Indel_Length_bp"
] = np.nan

aa_events[
    "Frame_Effect"
] = np.nan

aa_events[
    "Indel_Sequence"
] = np.nan

aa_events[
    "Indel_Reference_Location"
] = np.nan

aa_events[
    "Script04_ORF_Status"
] = aa_df[
    "Script04_ORF_Status"
]

aa_events[
    "Script03_Alignment_QC"
] = aa_df[
    "Script03_Alignment_QC"
]

aa_events[
    "ORF_Specific_Identity_pct"
] = aa_df[
    "ORF_Specific_Identity_pct"
]

aa_events[
    "Region"
] = aa_df.apply(
    lambda row:
        (
            approximate_region(
                row["ORF"],
                row["AA_Position"],
            )
            if row["ORF"] == "ORF1"
            else (
                row["Adney_Position_Region"]
                if pd.notna(row["Adney_Position_Region"])
                else approximate_region(
                    row["ORF"],
                    row["AA_Position"],
                )
            )
        ),
    axis=1,
)

aa_events[
    "Adney_Mapped"
] = aa_df[
    "Adney_Mapped"
]

aa_events[
    "Adney_Mutant_ID"
] = aa_df[
    "Adney_Mutant_ID"
]

aa_events[
    "Adney_Window_Start_AA"
] = aa_df[
    "Adney_Window_Start_AA"
]

aa_events[
    "Adney_Window_End_AA"
] = aa_df[
    "Adney_Window_End_AA"
]

aa_events[
    "Adney_Window_WT_Residues"
] = aa_df[
    "Adney_Window_WT_Residues"
]

aa_events[
    "Adney_RetroT_pct_WT"
] = aa_df[
    "Adney_RetroT_pct_WT"
]

aa_events[
    "Adney_Original_Activity_Class"
] = aa_df[
    "Adney_Original_Activity_Class"
]

aa_events[
    "Project_RetroT_Context"
] = aa_df[
    "Project_RetroT_Context"
]

aa_events[
    "Mapping_Status"
] = aa_df[
    "Mapping_Status"
]

aa_events[
    "Grantham_WT_to_Natural"
] = aa_df[
    "Grantham_WT_to_Natural"
]

aa_events[
    "Grantham_WT_to_Alanine"
] = aa_df[
    "Grantham_WT_to_Alanine"
]

aa_events[
    "BLOSUM62_WT_to_Natural"
] = aa_df[
    "BLOSUM62_WT_to_Natural"
]

aa_events[
    "BLOSUM62_WT_to_Alanine"
] = aa_df[
    "BLOSUM62_WT_to_Alanine"
]

aa_events[
    "Metric_Directional_Agreement"
] = aa_df[
    "Metric_Directional_Agreement"
]

aa_events[
    "Position_Level_Adney_Comparability"
] = aa_df[
    "Position_Level_Adney_Comparability"
]

aa_events[
    "Consequence_Class"
] = aa_df.apply(
    aa_consequence_class,
    axis=1,
)

aa_events[
    "Predicted_Functional_Consequence"
] = aa_df.apply(
    aa_predicted_consequence,
    axis=1,
)

aa_events[
    "Evidence_Basis"
] = (
    "Sequence consequence"
    " + Adney window experimental context"
    " + Grantham/BLOSUM62 substitution context"
)

aa_events[
    "Interpretation_Caveat"
] = (
    "Observed difference is relative to L1-RP and is not "
    "automatically a locus-specific mutation. Adney RetroT "
    "is window-level evidence and is not measured activity "
    "of this natural genomic sequence event."
)

aa_events[
    "Event_ID"
] = (
    aa_events["Locus_ID"].astype(str)
    + "|"
    + aa_events["ORF"].astype(str)
    + "|AA|"
    + aa_events["AA_Position"].astype(str)
    + "|"
    + aa_events["Observed_Event_vs_L1RP"].astype(str)
)


# ============================================================
# BUILD INDEL EVENT TABLE
# ============================================================

print("Building ORF-indel event rows...")


metadata_for_indels = (
    locus_base[
        [
            "Locus_ID",
            "Chromosome",
            "Start_1based",
            "End_1based",
            "Strand",
        ]
    ]
)


indel_work = indel_df.merge(
    metadata_for_indels,
    on="Locus_ID",
    how="left",
    validate="many_to_one",
)


indel_events = pd.DataFrame()


indel_events[
    "Locus_ID"
] = indel_work[
    "Locus_ID"
]

indel_events[
    "Chromosome"
] = indel_work[
    "Chromosome"
]

indel_events[
    "Start_1based"
] = indel_work[
    "Start_1based"
]

indel_events[
    "End_1based"
] = indel_work[
    "End_1based"
]

indel_events[
    "Strand"
] = indel_work[
    "Strand"
]

indel_events[
    "Primary_Subfamily"
] = indel_work[
    "Primary_Subfamily"
]

indel_events[
    "ORF"
] = indel_work[
    "ORF"
]

indel_events[
    "Event_Type"
] = "INDEL"

indel_events[
    "Observed_Event_vs_L1RP"
] = (
    indel_work[
        "Length_bp"
    ]
    .astype(
        int
    )
    .astype(
        str
    )
    + "bp_"
    + indel_work[
        "Indel_Type"
    ]
    .astype(
        str
    )
    .str.upper()
    + "_at_"
    + indel_work[
        "Reference_Location"
    ]
    .astype(
        str
    )
    + "_near_AA"
    + indel_work[
        "Approx_AA_Position"
    ]
    .astype(
        str
    )
)

indel_events[
    "Sequence_Consequence"
] = indel_work[
    "Frame_Effect"
].map(
    {
        "FRAMESHIFT":
            "FRAMESHIFT_INDEL",

        "IN_FRAME":
            "IN_FRAME_INDEL",
    }
).fillna(
    "INDEL"
)

indel_events[
    "AA_Position"
] = indel_work[
    "Approx_AA_Position"
]

indel_events[
    "Reference_AA"
] = np.nan

indel_events[
    "Observed_AA"
] = np.nan

indel_events[
    "Reference_Codon"
] = np.nan

indel_events[
    "Observed_Codon"
] = np.nan

indel_events[
    "Nucleotide_Differences"
] = np.nan

indel_events[
    "Indel_Type"
] = indel_work[
    "Indel_Type"
]

indel_events[
    "Indel_Length_bp"
] = indel_work[
    "Length_bp"
]

indel_events[
    "Frame_Effect"
] = indel_work[
    "Frame_Effect"
]

indel_events[
    "Indel_Sequence"
] = indel_work[
    "Sequence"
]

indel_events[
    "Indel_Reference_Location"
] = indel_work[
    "Reference_Location"
]

indel_events[
    "Script04_ORF_Status"
] = indel_work[
    "Script04_ORF_Status"
]

indel_events[
    "Script03_Alignment_QC"
] = indel_work[
    "Script03_Alignment_QC"
]

indel_events[
    "ORF_Specific_Identity_pct"
] = indel_work[
    "ORF_Specific_Identity_pct"
]

indel_events[
    "Region"
] = indel_work.apply(
    lambda row:
        approximate_region(
            row[
                "ORF"
            ],
            row[
                "Approx_AA_Position"
            ],
        ),
    axis=1,
)

indel_events[
    "Adney_Mapped"
] = False

indel_events[
    "Adney_Mutant_ID"
] = np.nan

indel_events[
    "Adney_Window_Start_AA"
] = np.nan

indel_events[
    "Adney_Window_End_AA"
] = np.nan

indel_events[
    "Adney_Window_WT_Residues"
] = np.nan

indel_events[
    "Adney_RetroT_pct_WT"
] = np.nan

indel_events[
    "Adney_Original_Activity_Class"
] = np.nan

indel_events[
    "Project_RetroT_Context"
] = np.nan

indel_events[
    "Mapping_Status"
] = "NOT_APPLICABLE_INDEL"

indel_events[
    "Grantham_WT_to_Natural"
] = np.nan

indel_events[
    "Grantham_WT_to_Alanine"
] = np.nan

indel_events[
    "BLOSUM62_WT_to_Natural"
] = np.nan

indel_events[
    "BLOSUM62_WT_to_Alanine"
] = np.nan

indel_events[
    "Metric_Directional_Agreement"
] = np.nan

indel_events[
    "Position_Level_Adney_Comparability"
] = np.nan

indel_events[
    "Consequence_Class"
] = indel_work.apply(
    indel_consequence_class,
    axis=1,
)

indel_events[
    "Predicted_Functional_Consequence"
] = indel_work.apply(
    indel_predicted_consequence,
    axis=1,
)

indel_events[
    "Evidence_Basis"
] = (
    "ORF-specific indel size and reading-frame consequence"
)

indel_events[
    "Interpretation_Caveat"
] = (
    "Observed indel is relative to L1-RP. Sequence-level "
    "coding consequence does not by itself establish "
    "retrotransposition activity of the genomic locus."
)

indel_events[
    "Event_ID"
] = (
    indel_events["Locus_ID"].astype(str)
    + "|"
    + indel_events["ORF"].astype(str)
    + "|INDEL|"
    + indel_events["Indel_Reference_Location"].astype(str)
    + "|"
    + indel_events["Indel_Type"].astype(str)
    + "|"
    + indel_events["Indel_Length_bp"].astype(str)
    + "bp|"
    + indel_events["Indel_Sequence"].astype(str)
)


# ============================================================
# COMBINE ALL OBSERVED EVENTS
# ============================================================

event_columns = list(
    aa_events.columns
)


indel_events = indel_events[
    event_columns
]


event_df = pd.concat(
    [
        aa_events,
        indel_events,
    ],
    ignore_index=True,
)

event_df.insert(
    0,
    "Event_ID",
    event_df.pop("Event_ID"),
)


event_df = event_df.sort_values(
    by=[
        "Chromosome",
        "Start_1based",
        "Locus_ID",
        "ORF",
        "AA_Position",
        "Event_Type",
    ],
    kind="stable",
).reset_index(
    drop=True
)


# ============================================================
# BUILD LOCUS-LEVEL COUNTS
# ============================================================

print("Building locus-level summary...")


def aggregate_locus_events(
    group
):

    aa = group.loc[
        group[
            "Event_Type"
        ]
        == "AA_DIFFERENCE"
    ]

    indels = group.loc[
        group[
            "Event_Type"
        ]
        == "INDEL"
    ]

    missense = aa.loc[
        aa[
            "Sequence_Consequence"
        ]
        == "MISSENSE"
    ]


    return pd.Series(
        {
            "Total_Observed_Events":
                len(
                    group
                ),

            "AA_Difference_Events":
                len(
                    aa
                ),

            "Missense_Events":
                int(
                    (
                        aa[
                            "Sequence_Consequence"
                        ]
                        == "MISSENSE"
                    ).sum()
                ),

            "Nonsense_Events":
                int(
                    (
                        aa[
                            "Sequence_Consequence"
                        ]
                        == "NONSENSE"
                    ).sum()
                ),

            "Start_Lost_Events":
                int(
                    (
                        aa[
                            "Sequence_Consequence"
                        ]
                        == "START_LOST"
                    ).sum()
                ),

            "Stop_Lost_Events":
                int(
                    (
                        aa[
                            "Sequence_Consequence"
                        ]
                        == "STOP_LOST"
                    ).sum()
                ),

            "Indel_Events":
                len(
                    indels
                ),

            "Frameshift_Indel_Events":
                int(
                    (
                        indels[
                            "Frame_Effect"
                        ]
                        == "FRAMESHIFT"
                    ).sum()
                ),

            "In_Frame_Indel_Events":
                int(
                    (
                        indels[
                            "Frame_Effect"
                        ]
                        == "IN_FRAME"
                    ).sum()
                ),

            "Adney_Poor_Window_Missense":
                int(
                    (
                        missense[
                            "Adney_Original_Activity_Class"
                        ]
                        == "Poor (<25% WT)"
                    ).sum()
                ),

            "Adney_Reduced_25_80_Window_Missense":
                int(
                    (
                        missense[
                            "Adney_Original_Activity_Class"
                        ]
                        == "Reduced (25-80% WT)"
                    ).sum()
                ),

            "Adney_High_gt80_Window_Missense":
                int(
                    (
                        missense[
                            "Adney_Original_Activity_Class"
                        ]
                        == "High (>80% WT)"
                    ).sum()
                ),

            "Natural_Equals_Adney_Alanine":
                int(
                    (
                        missense[
                            "Position_Level_Adney_Comparability"
                        ]
                        ==
                        "NATURAL_EQUALS_ADNEY_ALANINE_AT_THIS_POSITION"
                    ).sum()
                ),

            "WT_Alanine_Not_Directly_Perturbed":
                int(
                    (
                        missense[
                            "Position_Level_Adney_Comparability"
                        ]
                        ==
                        "ADNEY_DID_NOT_CHANGE_WT_ALANINE_POSITION"
                    ).sum()
                ),

            "Natural_More_Similar_to_WT_than_Alanine":
                int(
                    (
                        missense[
                            "Metric_Directional_Agreement"
                        ]
                        ==
                        "CONCORDANT_NATURAL_MORE_SIMILAR_TO_WT_THAN_ALANINE"
                    ).sum()
                ),

            "Natural_Less_Similar_to_WT_than_Alanine":
                int(
                    (
                        missense[
                            "Metric_Directional_Agreement"
                        ]
                        ==
                        "CONCORDANT_NATURAL_LESS_SIMILAR_TO_WT_THAN_ALANINE"
                    ).sum()
                ),

            "Mixed_or_Tied_Substitution_Metrics":
                int(
                    (
                        missense[
                            "Metric_Directional_Agreement"
                        ]
                        ==
                        "MIXED_OR_PARTIALLY_TIED_METRICS"
                    ).sum()
                ),
        }
    )


event_summary = (
    event_df
    .groupby(
        "Locus_ID",
        sort=False,
    )
    .apply(
        aggregate_locus_events
    )
    .reset_index()
)


locus_summary = locus_base.merge(
    event_summary,
    on="Locus_ID",
    how="left",
    validate="one_to_one",
)


count_columns = [
    "Total_Observed_Events",
    "AA_Difference_Events",
    "Missense_Events",
    "Nonsense_Events",
    "Start_Lost_Events",
    "Stop_Lost_Events",
    "Indel_Events",
    "Frameshift_Indel_Events",
    "In_Frame_Indel_Events",
    "Adney_Poor_Window_Missense",
    "Adney_Reduced_25_80_Window_Missense",
    "Adney_High_gt80_Window_Missense",
    "Natural_Equals_Adney_Alanine",
    "WT_Alanine_Not_Directly_Perturbed",
    "Natural_More_Similar_to_WT_than_Alanine",
    "Natural_Less_Similar_to_WT_than_Alanine",
    "Mixed_or_Tied_Substitution_Metrics",
]


for column in count_columns:

    locus_summary[
        column
    ] = (
        locus_summary[
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
# LOCUS-LEVEL INTERPRETATION
# ============================================================

def locus_interpretation(
    row
):

    orf1 = str(
        row[
            "ORF1_Sequence_Status"
        ]
    )

    orf2 = str(
        row[
            "ORF2_Sequence_Status"
        ]
    )


    severe_count = (
        int(
            row[
                "Nonsense_Events"
            ]
        )
        +
        int(
            row[
                "Start_Lost_Events"
            ]
        )
        +
        int(
            row[
                "Frameshift_Indel_Events"
            ]
        )
    )


    if severe_count > 0:

        return (
            "Sequence-level coding disruption detected in at "
            "least one ORF (start loss, premature stop and/or "
            "frameshift). This supports reduced coding potential "
            "but does not directly measure locus activity."
        )


    if (
        orf1
        ==
        "FULL_SPAN_NO_FRAME_STOP_DISRUPTION"
        and orf2
        ==
        "FULL_SPAN_NO_FRAME_STOP_DISRUPTION"
    ):

        if (
            int(
                row[
                    "Adney_Poor_Window_Missense"
                ]
            )
            > 0
        ):

            return (
                "Both ORFs retain full-span frame/stop integrity. "
                "Missense differences occur in one or more "
                "experimentally sensitive Adney windows; "
                "natural-substitution effects remain uncertain."
            )


        return (
            "Both ORFs retain full-span frame/stop integrity "
            "with no obvious start/frameshift/premature-stop "
            "disruption detected. This indicates retained "
            "coding structure, not demonstrated activity."
        )


    if (
        "PARTIAL_SPAN"
        in orf1
        or "PARTIAL_SPAN"
        in orf2
    ):

        return (
            "At least one ORF is only partially represented, "
            "so complete coding integrity cannot be assessed."
        )


    return (
        "At least one ORF contains sequence-level frame/stop "
        "disruption according to Script 04. This indicates "
        "reduced coding integrity but does not directly measure "
        "retrotransposition activity."
    )


locus_summary[
    "Both_ORFs_No_Frame_Stop_Disruption"
] = (
    (
        locus_summary[
            "ORF1_Sequence_Status"
        ]
        ==
        "FULL_SPAN_NO_FRAME_STOP_DISRUPTION"
    )
    &
    (
        locus_summary[
            "ORF2_Sequence_Status"
        ]
        ==
        "FULL_SPAN_NO_FRAME_STOP_DISRUPTION"
    )
)


locus_summary[
    "Locus_Level_Functional_Interpretation"
] = locus_summary.apply(
    locus_interpretation,
    axis=1,
)


locus_summary[
    "Activity_Claim"
] = (
    "No direct activity claim. Annotation is based on "
    "sequence structure, coding consequences, Adney window "
    "experiments, and substitution-context metrics."
)


# ============================================================
# QC
# ============================================================

print("\nRunning final QC...")


if event_df[
    "Locus_ID"
].isna().any():

    raise ValueError(
        "Final event table contains missing Locus_ID values."
    )


if event_df[
    "Event_ID"
].isna().any():

    raise ValueError(
        "Final event table contains missing Event_ID values."
    )


duplicate_event_ids = int(
    event_df["Event_ID"].duplicated().sum()
)

if duplicate_event_ids:

    raise ValueError(
        f"Final event table contains "
        f"{duplicate_event_ids:,} duplicate Event_ID values."
    )


if len(
    locus_summary
) != len(
    locus_base
):

    raise ValueError(
        "Locus summary row count changed unexpectedly."
    )


expected_aa = len(
    aa_df
)

observed_aa = int(
    (
        event_df[
            "Event_Type"
        ]
        == "AA_DIFFERENCE"
    ).sum()
)


if observed_aa != expected_aa:

    raise ValueError(
        "AA event count mismatch: "
        f"expected {expected_aa:,}, "
        f"observed {observed_aa:,}."
    )


expected_indels = len(
    indel_df
)

observed_indels = int(
    (
        event_df[
            "Event_Type"
        ]
        == "INDEL"
    ).sum()
)


if observed_indels != expected_indels:

    raise ValueError(
        "Indel event count mismatch: "
        f"expected {expected_indels:,}, "
        f"observed {observed_indels:,}."
    )


# ============================================================
# SAVE CSV OUTPUTS
# ============================================================

print("\nSaving final outputs...")


event_df.to_csv(
    EVENT_OUTPUT,
    index=False,
    compression="gzip",
)


locus_summary.to_csv(
    LOCUS_OUTPUT,
    index=False,
    compression="gzip",
)


if MODE == "TEST":

    event_df.to_csv(
        EVENT_OUTPUT_PLAIN,
        index=False,
    )

    locus_summary.to_csv(
        LOCUS_OUTPUT_PLAIN,
        index=False,
    )


# ============================================================
# EXCEL WORKBOOK
# ============================================================

consequence_legend = pd.DataFrame(
    [
        [
            "LIKELY_SEVERE_CODING_DISRUPTION_START_LOST",
            "Start codon/initiator lost; strong sequence-level "
            "evidence for ORF disruption.",
        ],
        [
            "LIKELY_SEVERE_CODING_DISRUPTION_NONSENSE",
            "Premature stop introduced; strong sequence-level "
            "evidence for truncated protein.",
        ],
        [
            "LIKELY_SEVERE_CODING_DISRUPTION_FRAMESHIFT",
            "Indel shifts reading frame; strong sequence-level "
            "evidence for downstream coding disruption.",
        ],
        [
            "ALTERED_TERMINATION_STOP_LOST",
            "Reference terminal stop lost; termination is altered.",
        ],
        [
            "IN_FRAME_PROTEIN_ALTERATION",
            "Indel preserves downstream reading frame but changes "
            "local protein sequence.",
        ],
        [
            "MISSENSE_WITH_ADNEY_WINDOW_CONTEXT",
            "Missense lies in an Adney-tested window; RetroT is "
            "window-level experimental context only.",
        ],
        [
            "MISSENSE_MATCHES_ADNEY_ALANINE_AT_POSITION",
            "Natural residue equals the alanine introduced by "
            "Adney at that specific position, but the whole "
            "window was mutated experimentally.",
        ],
        [
            "MISSENSE_ADNEY_POSITION_NOT_DIRECTLY_PERTURBED",
            "L1-RP WT residue is already alanine; Adney did not "
            "change that specific position.",
        ],
    ],
    columns=[
        "Consequence_Class",
        "Meaning",
    ],
)


qc_summary = pd.DataFrame(
    [
        [
            "Run mode",
            MODE,
        ],
        [
            "Loci in final summary",
            len(
                locus_summary
            ),
        ],
        [
            "AA-difference events",
            observed_aa,
        ],
        [
            "ORF indel events",
            observed_indels,
        ],
        [
            "Total final event rows",
            len(
                event_df
            ),
        ],
        [
            "Loci with both ORFs no frame/stop disruption",
            int(
                locus_summary[
                    "Both_ORFs_No_Frame_Stop_Disruption"
                ].sum()
            ),
        ],
        [
            "Scientific interpretation",
            (
                "Sequence-based predicted consequence and "
                "experimental context; not measured genomic "
                "locus activity."
            ),
        ],
    ],
    columns=[
        "Metric",
        "Value",
    ],
)


# Curated event columns for the human-readable Excel sheet.
excel_event_columns = [
    "Event_ID",
    "Locus_ID",
    "Chromosome",
    "Start_1based",
    "End_1based",
    "Strand",
    "Primary_Subfamily",
    "ORF",
    "Event_Type",
    "Observed_Event_vs_L1RP",
    "Sequence_Consequence",
    "AA_Position",
    "Region",
    "Indel_Reference_Location",
    "Indel_Type",
    "Indel_Length_bp",
    "Frame_Effect",
    "Script04_ORF_Status",
    "ORF_Specific_Identity_pct",
    "Adney_Mutant_ID",
    "Adney_RetroT_pct_WT",
    "Adney_Original_Activity_Class",
    "Position_Level_Adney_Comparability",
    "Grantham_WT_to_Natural",
    "Grantham_WT_to_Alanine",
    "BLOSUM62_WT_to_Natural",
    "BLOSUM62_WT_to_Alanine",
    "Metric_Directional_Agreement",
    "Consequence_Class",
    "Predicted_Functional_Consequence",
    "Interpretation_Caveat",
]


if WRITE_EXCEL_WORKBOOK:

    with pd.ExcelWriter(
        WORKBOOK_OUTPUT,
        engine="openpyxl",
    ) as writer:

        qc_summary.to_excel(
            writer,
            sheet_name="QC_Summary",
            index=False,
        )

        consequence_legend.to_excel(
            writer,
            sheet_name="Consequence_Legend",
            index=False,
        )

        locus_summary.to_excel(
            writer,
            sheet_name="Locus_Summary",
            index=False,
        )

        if WRITE_EVENT_TABLE_TO_EXCEL:

            event_df[
                excel_event_columns
            ].to_excel(
                writer,
                sheet_name="Mutation_Consequence_Table",
                index=False,
            )


# ============================================================
# FINAL REPORT
# ============================================================

print("\n============================================================")
print("SCRIPT 08 V2 COMPLETE")
print("============================================================")

print(
    f"\nFinal loci: "
    f"{len(locus_summary):,}"
)

print(
    f"AA-difference events: "
    f"{observed_aa:,}"
)

print(
    f"ORF-indel events: "
    f"{observed_indels:,}"
)

print(
    f"Total event rows: "
    f"{len(event_df):,}"
)

print(
    "Duplicate Event_ID values: "
    f"{int(event_df['Event_ID'].duplicated().sum()):,}"
)

print(
    "\nOutputs:"
)

print(
    "  Event consequence table:",
    EVENT_OUTPUT,
)

print(
    "  Locus summary:",
    LOCUS_OUTPUT,
)

if WRITE_EXCEL_WORKBOOK:

    print(
        "  Final Excel workbook:",
        WORKBOOK_OUTPUT,
    )

print(
    "\nIMPORTANT:"
)

print(
    "Predicted consequences are sequence/evidence-based "
    "annotations, not direct measurements of genomic LINE-1 "
    "retrotransposition activity."
)
