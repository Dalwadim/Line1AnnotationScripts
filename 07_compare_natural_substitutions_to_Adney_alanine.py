# ============================================================
# LINE-1 FUNCTIONAL ANNOTATION PROJECT
# SCRIPT 07
#
# Compare natural amino-acid substitutions with the alanine
# substitution used in the Adney et al. alanine-scan library.
#
# Grantham:
#   larger distance = greater physicochemical difference.
#
# BLOSUM62:
#   larger score = substitution is more favored/observed in the
#   BLOSUM62 sequence-evolution framework.
#
# IMPORTANT:
#   These metrics are CONTEXT, not predictions of LINE-1
#   retrotransposition activity.
#
#   Adney RetroT remains WINDOW-LEVEL experimental evidence from
#   a multi-residue alanine construct. Do not assign that RetroT
#   value to an individual natural substitution.
#
# Special case:
#   If the L1-RP WT residue is already alanine (A), Adney did
#   not change that specific residue when making the alanine
#   window. Script 07 flags those rows separately.
#
# References:
#   Grantham R. Science. 1974;185:862-864.
#   DOI: 10.1126/science.185.4154.862
#
#   Henikoff S, Henikoff JG. PNAS. 1992;89:10915-10919.
#   DOI: 10.1073/pnas.89.22.10915
# ============================================================

import random
from pathlib import Path

import numpy as np
import pandas as pd
from Bio.Align import substitution_matrices


# ============================================================
# SETTINGS
# ============================================================

RUN_MODE = "FULL"          # change to "FULL" after QC
TEST_N_LOCI = 25
RANDOM_SEED = 20260922
WRITE_EXCEL_SUMMARY = True

TEST_LOCUS_IDS = [
    "hg38.RE.chr1.100199603-100206089.+",
    "hg38.RE.chr1.104770248-104776279.-",
    "hg38.RE.chr1.104843834-104849865.-",
]


# ============================================================
# PATHS
# ============================================================

SCRIPT_FOLDER = Path(__file__).resolve().parent
SCRIPT06_FOLDER = SCRIPT_FOLDER / "06_Adney_mapping"
OUTPUT_FOLDER = SCRIPT_FOLDER / "07_substitution_context"
OUTPUT_FOLDER.mkdir(parents=True, exist_ok=True)


def first_existing(paths, description):
    for path in paths:
        if path.exists():
            return path

    attempted = "\n".join(f"  - {p}" for p in paths)

    raise FileNotFoundError(
        f"\nCould not find {description}.\n"
        f"Checked:\n{attempted}\n"
    )


MAPPED_INPUT = first_existing(
    [
        SCRIPT06_FOLDER
        / "FULL_LINE1_AA_differences_Adney_mapped.csv.gz",

        SCRIPT_FOLDER
        / "FULL_LINE1_AA_differences_Adney_mapped.csv.gz",
    ],
    "Script 06 FULL mapped AA-difference table",
)


MODE = RUN_MODE.upper().strip()

if MODE not in {"TEST", "FULL"}:
    raise ValueError(
        "RUN_MODE must be either 'TEST' or 'FULL'."
    )


DETAILED_OUTPUT = (
    OUTPUT_FOLDER
    / f"{MODE}_LINE1_substitution_context.csv.gz"
)

PRIORITY_OUTPUT = (
    OUTPUT_FOLDER
    / f"{MODE}_LINE1_priority_clean_ORFs_substitution_context.csv.gz"
)

LOCUS_ORF_SUMMARY_OUTPUT = (
    OUTPUT_FOLDER
    / f"{MODE}_LINE1_substitution_locus_ORF_summary.csv.gz"
)

SUMMARY_XLSX = (
    OUTPUT_FOLDER
    / f"{MODE}_LINE1_substitution_context_summary.xlsx"
)

if MODE == "TEST":
    DETAILED_OUTPUT_PLAIN = (
        OUTPUT_FOLDER
        / "TEST_LINE1_substitution_context.csv"
    )

    PRIORITY_OUTPUT_PLAIN = (
        OUTPUT_FOLDER
        / "TEST_LINE1_priority_clean_ORFs_substitution_context.csv"
    )


# ============================================================
# AMINO-ACID METRICS
# ============================================================

AA_ORDER = list("ARNDCQEGHILKMFPSTWYV")
STANDARD_AA = set(AA_ORDER)

# Standard Grantham distance matrix in AA_ORDER:
# A R N D C Q E G H I L K M F P S T W Y V
GRANTHAM_MATRIX_TEXT = """
0 112 111 126 195 91 107 60 86 94 96 106 84 113 27 99 58 148 112 64
112 0 86 96 180 43 54 125 29 97 102 26 91 97 103 110 71 101 77 96
111 86 0 23 139 46 42 80 68 149 153 94 142 158 91 46 65 174 143 133
126 96 23 0 154 61 45 94 81 168 172 101 160 177 108 65 85 181 160 152
195 180 139 154 0 154 170 159 174 198 198 202 196 205 169 112 149 215 194 192
91 43 46 61 154 0 29 87 24 109 113 53 101 116 76 68 42 130 99 96
107 54 42 45 170 29 0 98 40 134 138 56 126 140 93 80 65 152 122 121
60 125 80 94 159 87 98 0 98 135 138 127 127 153 42 56 59 184 147 109
86 29 68 81 174 24 40 98 0 94 99 32 87 100 77 89 47 115 83 84
94 97 149 168 198 109 134 135 94 0 5 102 10 21 95 142 89 61 33 29
96 102 153 172 198 113 138 138 99 5 0 107 15 22 98 145 92 61 36 32
106 26 94 101 202 53 56 127 32 102 107 0 95 102 103 121 78 110 85 97
84 91 142 160 196 101 126 127 87 10 15 95 0 28 87 135 81 67 36 21
113 97 158 177 205 116 140 153 100 21 22 102 28 0 114 155 103 40 22 50
27 103 91 108 169 76 93 42 77 95 98 103 87 114 0 74 38 147 110 68
99 110 46 65 112 68 80 56 89 142 145 121 135 155 74 0 58 177 144 124
58 71 65 85 149 42 65 59 47 89 92 78 81 103 38 58 0 128 92 69
148 101 174 181 215 130 152 184 115 61 61 110 67 40 147 177 128 0 37 88
112 77 143 160 194 99 122 147 83 33 36 85 36 22 110 144 92 37 0 55
64 96 133 152 192 96 121 109 84 29 32 97 21 50 68 124 69 88 55 0
""".strip()

GRANTHAM_ARRAY = np.array(
    [
        [int(x) for x in line.split()]
        for line in GRANTHAM_MATRIX_TEXT.splitlines()
    ],
    dtype=int,
)

if GRANTHAM_ARRAY.shape != (20, 20):
    raise ValueError(
        f"Unexpected Grantham matrix shape: "
        f"{GRANTHAM_ARRAY.shape}"
    )

AA_INDEX = {
    aa: index
    for index, aa in enumerate(AA_ORDER)
}

BLOSUM62 = substitution_matrices.load("BLOSUM62")


def grantham_distance(aa1, aa2):
    aa1 = str(aa1).upper()
    aa2 = str(aa2).upper()

    if aa1 not in STANDARD_AA or aa2 not in STANDARD_AA:
        return np.nan

    return int(
        GRANTHAM_ARRAY[
            AA_INDEX[aa1],
            AA_INDEX[aa2],
        ]
    )


def blosum62_score(aa1, aa2):
    aa1 = str(aa1).upper()
    aa2 = str(aa2).upper()

    if aa1 not in STANDARD_AA or aa2 not in STANDARD_AA:
        return np.nan

    return float(BLOSUM62[aa1, aa2])


def validate_metric_tables():
    # Symmetry / diagonal checks.
    if not np.array_equal(
        GRANTHAM_ARRAY,
        GRANTHAM_ARRAY.T,
    ):
        raise ValueError("Grantham matrix is not symmetric.")

    if not np.all(
        np.diag(GRANTHAM_ARRAY) == 0
    ):
        raise ValueError("Grantham diagonal is not zero.")

    # Known standard values.
    checks = {
        ("I", "L"): 5,
        ("A", "P"): 27,
        ("D", "E"): 45,
        ("A", "R"): 112,
        ("C", "W"): 215,
    }

    for (aa1, aa2), expected in checks.items():
        observed = grantham_distance(aa1, aa2)

        if observed != expected:
            raise ValueError(
                f"Grantham self-check failed for "
                f"{aa1}-{aa2}: expected {expected}, "
                f"observed {observed}"
            )

    if blosum62_score("A", "A") != 4:
        raise ValueError(
            "Unexpected BLOSUM62 A-A score."
        )

    for aa1 in AA_ORDER:
        for aa2 in AA_ORDER:
            if (
                blosum62_score(aa1, aa2)
                != blosum62_score(aa2, aa1)
            ):
                raise ValueError(
                    f"BLOSUM62 symmetry failed: "
                    f"{aa1}-{aa2}"
                )


# ============================================================
# INPUT VALIDATION
# ============================================================

REQUIRED_COLUMNS = {
    "Locus_ID",
    "Primary_Subfamily",
    "ORF",
    "AA_Position",
    "L1RP_AA",
    "Locus_AA",
    "AA_Difference_vs_L1RP",
    "Difference_Type",
    "Script04_ORF_Status",
    "Adney_Mutant_ID",
    "Adney_Window_Start_AA",
    "Adney_Window_End_AA",
    "Adney_Window_WT_Residues",
    "Adney_Alanine_Replacement",
    "Adney_Position_WT_AA",
    "Adney_RetroT_pct_WT",
    "Adney_Original_Activity_Class",
    "Project_RetroT_Context",
    "Mapping_Status",
    "Adney_Mapped",
}


def validate_input_columns(df):
    missing = sorted(
        REQUIRED_COLUMNS - set(df.columns)
    )

    if missing:
        raise ValueError(
            "Script 06 input is missing required columns:\n"
            + "\n".join(
                f"  - {column}"
                for column in missing
            )
        )


# ============================================================
# TEST SELECTION
# ============================================================

def select_test_loci(df):
    """
    Reproducible 25-locus test set.

    Deliberately includes examples of:
      1. WT already alanine.
      2. Natural residue = alanine.
      3. Clean ORF.
      4. Earlier known QC loci when available.
    """

    available = list(
        df["Locus_ID"].drop_duplicates()
    )
    available_set = set(available)
    selected = []

    def add(locus_id):
        if (
            locus_id in available_set
            and locus_id not in selected
        ):
            selected.append(locus_id)

    for locus_id in TEST_LOCUS_IDS:
        add(locus_id)

    wt_a = df.loc[
        (df["Mapping_Status"] == "MAPPED")
        & (df["Difference_Type"] == "MISSENSE")
        & (df["L1RP_AA"] == "A"),
        "Locus_ID",
    ]

    if len(wt_a):
        add(wt_a.iloc[0])

    natural_a = df.loc[
        (df["Mapping_Status"] == "MAPPED")
        & (df["Difference_Type"] == "MISSENSE")
        & (df["L1RP_AA"] != "A")
        & (df["Locus_AA"] == "A"),
        "Locus_ID",
    ]

    if len(natural_a):
        add(natural_a.iloc[0])

    clean = df.loc[
        (df["Mapping_Status"] == "MAPPED")
        & (df["Difference_Type"] == "MISSENSE")
        & (
            df["Script04_ORF_Status"]
            == "FULL_SPAN_NO_FRAME_STOP_DISRUPTION"
        ),
        "Locus_ID",
    ]

    if len(clean):
        add(clean.iloc[0])

    remaining = [
        locus_id
        for locus_id in available
        if locus_id not in selected
    ]

    rng = random.Random(RANDOM_SEED)

    needed = min(
        max(0, TEST_N_LOCI - len(selected)),
        len(remaining),
    )

    selected.extend(
        rng.sample(remaining, needed)
    )

    return selected


# ============================================================
# ROW SCORING
# ============================================================

def score_row(row):
    """
    Add Grantham/BLOSUM62 context to one Script 06 row.
    """

    out = {
        "Substitution_Scoring_Status": None,
        "Adney_Position_Substitution": None,
        "Adney_Changed_This_Position": False,
        "Natural_Equals_Adney_Alanine": False,
        "Position_Level_Adney_Comparability": None,

        "Grantham_WT_to_Natural": np.nan,
        "Grantham_WT_to_Alanine": np.nan,
        "Grantham_Natural_to_Alanine": np.nan,
        "Grantham_Delta_Natural_minus_Alanine": np.nan,
        "Grantham_Comparison_vs_Alanine": None,

        "BLOSUM62_WT_to_Natural": np.nan,
        "BLOSUM62_WT_to_Alanine": np.nan,
        "BLOSUM62_Natural_to_Alanine": np.nan,
        "BLOSUM62_Delta_Natural_minus_Alanine": np.nan,
        "BLOSUM62_Comparison_vs_Alanine": None,

        "Metric_Directional_Agreement": None,
        "Script07_Interpretation": None,
    }

    mapping_status = str(row["Mapping_Status"])
    difference_type = str(row["Difference_Type"])

    if mapping_status != "MAPPED":
        out["Substitution_Scoring_Status"] = (
            "NOT_APPLICABLE_NO_ADNEY_POSITION"
        )
        out["Position_Level_Adney_Comparability"] = (
            "NO_ADNEY_POSITION_MAPPING"
        )
        out["Script07_Interpretation"] = (
            "No Adney residue-level comparison was made "
            "because this row is outside the Adney-tested "
            "position map."
        )
        return out

    if difference_type != "MISSENSE":
        out["Substitution_Scoring_Status"] = (
            f"NOT_APPLICABLE_{difference_type}"
        )
        out["Position_Level_Adney_Comparability"] = (
            "NONMISSENSE_NOT_SCOREABLE"
        )
        out["Script07_Interpretation"] = (
            "Grantham and BLOSUM62 are not applied because "
            "this is not a standard amino-acid missense "
            "substitution."
        )
        return out

    wt = str(row["L1RP_AA"]).upper()
    natural = str(row["Locus_AA"]).upper()
    adney_wt = str(
        row["Adney_Position_WT_AA"]
    ).upper()

    if (
        wt not in STANDARD_AA
        or natural not in STANDARD_AA
    ):
        out["Substitution_Scoring_Status"] = (
            "NOT_APPLICABLE_NONSTANDARD_AMINO_ACID"
        )
        out["Position_Level_Adney_Comparability"] = (
            "NONSTANDARD_AMINO_ACID"
        )
        return out

    if wt != adney_wt:
        raise ValueError(
            "\nL1-RP / Adney WT residue mismatch:\n"
            f"Locus: {row['Locus_ID']}\n"
            f"ORF: {row['ORF']}\n"
            f"AA: {row['AA_Position']}\n"
            f"L1-RP: {wt}\n"
            f"Adney: {adney_wt}"
        )

    alanine = "A"

    out["Substitution_Scoring_Status"] = "SCORED"
    out["Adney_Position_Substitution"] = (
        f"{wt}{int(row['AA_Position'])}A"
    )
    out["Adney_Changed_This_Position"] = (
        wt != alanine
    )
    out["Natural_Equals_Adney_Alanine"] = (
        natural == alanine
        and wt != alanine
    )

    # Natural substitution metrics.
    g_nat = grantham_distance(wt, natural)
    b_nat = blosum62_score(wt, natural)

    out["Grantham_WT_to_Natural"] = g_nat
    out["BLOSUM62_WT_to_Natural"] = b_nat

    # Natural residue compared directly with alanine.
    out["Grantham_Natural_to_Alanine"] = (
        grantham_distance(natural, alanine)
    )
    out["BLOSUM62_Natural_to_Alanine"] = (
        blosum62_score(natural, alanine)
    )

    # --------------------------------------------------------
    # WT already alanine.
    # --------------------------------------------------------
    if wt == alanine:
        out["Grantham_WT_to_Alanine"] = 0
        out["BLOSUM62_WT_to_Alanine"] = (
            blosum62_score("A", "A")
        )

        out["Position_Level_Adney_Comparability"] = (
            "ADNEY_DID_NOT_CHANGE_WT_ALANINE_POSITION"
        )
        out["Grantham_Comparison_vs_Alanine"] = (
            "NOT_APPLICABLE_WT_ALREADY_ALANINE"
        )
        out["BLOSUM62_Comparison_vs_Alanine"] = (
            "NOT_APPLICABLE_WT_ALREADY_ALANINE"
        )
        out["Metric_Directional_Agreement"] = (
            "NOT_APPLICABLE_WT_ALREADY_ALANINE"
        )
        out["Script07_Interpretation"] = (
            "L1-RP is already alanine at this residue. "
            "The Adney alanine-window construct therefore "
            "did not alter this specific position; the window "
            "RetroT phenotype reflects changes elsewhere in "
            "the multi-residue window."
        )
        return out

    # --------------------------------------------------------
    # WT is not alanine: compare WT->natural with WT->A.
    # --------------------------------------------------------
    g_ala = grantham_distance(wt, alanine)
    b_ala = blosum62_score(wt, alanine)

    g_delta = g_nat - g_ala
    b_delta = b_nat - b_ala

    out["Grantham_WT_to_Alanine"] = g_ala
    out["BLOSUM62_WT_to_Alanine"] = b_ala

    out[
        "Grantham_Delta_Natural_minus_Alanine"
    ] = g_delta

    out[
        "BLOSUM62_Delta_Natural_minus_Alanine"
    ] = b_delta

    if g_delta < 0:
        g_comp = (
            "NATURAL_LOWER_DISTANCE_THAN_ALANINE"
        )
    elif g_delta > 0:
        g_comp = (
            "NATURAL_HIGHER_DISTANCE_THAN_ALANINE"
        )
    else:
        g_comp = (
            "NATURAL_SAME_DISTANCE_AS_ALANINE"
        )

    if b_delta > 0:
        b_comp = (
            "NATURAL_HIGHER_SCORE_THAN_ALANINE"
        )
    elif b_delta < 0:
        b_comp = (
            "NATURAL_LOWER_SCORE_THAN_ALANINE"
        )
    else:
        b_comp = (
            "NATURAL_SAME_SCORE_AS_ALANINE"
        )

    out["Grantham_Comparison_vs_Alanine"] = (
        g_comp
    )
    out["BLOSUM62_Comparison_vs_Alanine"] = (
        b_comp
    )

    if natural == alanine:
        out["Position_Level_Adney_Comparability"] = (
            "NATURAL_EQUALS_ADNEY_ALANINE_AT_THIS_POSITION"
        )
        out["Metric_Directional_Agreement"] = (
            "NATURAL_EQUALS_ALANINE"
        )
        out["Script07_Interpretation"] = (
            "The natural residue is alanine, matching the "
            "alanine residue used by the Adney window at this "
            "position. The Adney RetroT value still belongs to "
            "the entire multi-residue window, not this single "
            "substitution."
        )

    elif g_delta < 0 and b_delta > 0:
        out["Position_Level_Adney_Comparability"] = (
            "DIRECT_RESIDUE_COMPARISON"
        )
        out["Metric_Directional_Agreement"] = (
            "CONCORDANT_NATURAL_MORE_SIMILAR_TO_WT_THAN_ALANINE"
        )
        out["Script07_Interpretation"] = (
            "Both metrics place the natural residue closer to "
            "the WT residue than alanine in their respective "
            "scoring frameworks. This is not a functional "
            "activity prediction."
        )

    elif g_delta > 0 and b_delta < 0:
        out["Position_Level_Adney_Comparability"] = (
            "DIRECT_RESIDUE_COMPARISON"
        )
        out["Metric_Directional_Agreement"] = (
            "CONCORDANT_NATURAL_LESS_SIMILAR_TO_WT_THAN_ALANINE"
        )
        out["Script07_Interpretation"] = (
            "Both metrics place the natural residue farther "
            "from the WT residue than alanine in their "
            "respective scoring frameworks. This is not a "
            "functional activity prediction."
        )

    elif g_delta == 0 and b_delta == 0:
        out["Position_Level_Adney_Comparability"] = (
            "DIRECT_RESIDUE_COMPARISON"
        )
        out["Metric_Directional_Agreement"] = (
            "BOTH_METRICS_TIED_WITH_ALANINE"
        )
        out["Script07_Interpretation"] = (
            "Both metrics score the natural substitution the "
            "same as the WT-to-alanine substitution."
        )

    else:
        out["Position_Level_Adney_Comparability"] = (
            "DIRECT_RESIDUE_COMPARISON"
        )
        out["Metric_Directional_Agreement"] = (
            "MIXED_OR_PARTIALLY_TIED_METRICS"
        )
        out["Script07_Interpretation"] = (
            "Grantham and BLOSUM62 do not provide fully "
            "concordant directionality relative to alanine. "
            "Keep the raw metric values rather than forcing a "
            "single interpretation."
        )

    return out


# ============================================================
# SUMMARY
# ============================================================

def count_value(series, value):
    return int((series == value).sum())


def build_locus_orf_summary(df):
    rows = []

    for (locus_id, orf), group in df.groupby(
        ["Locus_ID", "ORF"],
        sort=False,
    ):
        scored = group.loc[
            group[
                "Substitution_Scoring_Status"
            ]
            == "SCORED"
        ]

        rows.append(
            {
                "Locus_ID": locus_id,
                "Primary_Subfamily":
                    group["Primary_Subfamily"].iloc[0],
                "ORF": orf,
                "Script04_ORF_Status":
                    group["Script04_ORF_Status"].iloc[0],
                "Total_Script06_AA_Difference_Rows":
                    len(group),
                "Scored_Missense_Rows":
                    len(scored),
                "WT_Already_Alanine_Rows":
                    count_value(
                        scored[
                            "Position_Level_Adney_Comparability"
                        ],
                        "ADNEY_DID_NOT_CHANGE_WT_ALANINE_POSITION",
                    ),
                "Natural_Equals_Adney_Alanine_Rows":
                    count_value(
                        scored[
                            "Position_Level_Adney_Comparability"
                        ],
                        "NATURAL_EQUALS_ADNEY_ALANINE_AT_THIS_POSITION",
                    ),
                "Concordant_Natural_More_Similar_Rows":
                    count_value(
                        scored[
                            "Metric_Directional_Agreement"
                        ],
                        "CONCORDANT_NATURAL_MORE_SIMILAR_TO_WT_THAN_ALANINE",
                    ),
                "Concordant_Natural_Less_Similar_Rows":
                    count_value(
                        scored[
                            "Metric_Directional_Agreement"
                        ],
                        "CONCORDANT_NATURAL_LESS_SIMILAR_TO_WT_THAN_ALANINE",
                    ),
                "Both_Metrics_Tied_Rows":
                    count_value(
                        scored[
                            "Metric_Directional_Agreement"
                        ],
                        "BOTH_METRICS_TIED_WITH_ALANINE",
                    ),
                "Mixed_or_Partially_Tied_Rows":
                    count_value(
                        scored[
                            "Metric_Directional_Agreement"
                        ],
                        "MIXED_OR_PARTIALLY_TIED_METRICS",
                    ),
            }
        )

    return pd.DataFrame(rows)


# ============================================================
# LOAD DATA
# ============================================================

print("\n============================================================")
print("SCRIPT 07 — NATURAL AA VS ADNEY ALANINE CONTEXT")
print("============================================================")

print("\nRun mode:", MODE)
print("\nInput:")
print(MAPPED_INPUT)

validate_metric_tables()

print(
    "\nGrantham and BLOSUM62 tables passed "
    "internal validation."
)

input_df = pd.read_csv(MAPPED_INPUT)
validate_input_columns(input_df)

print(
    f"\nScript 06 rows loaded: "
    f"{len(input_df):,}"
)


# ============================================================
# TEST / FULL SELECTION
# ============================================================

if MODE == "FULL":
    work_df = input_df.copy()

else:
    selected_loci = select_test_loci(input_df)

    work_df = input_df.loc[
        input_df["Locus_ID"].isin(
            selected_loci
        )
    ].copy()

    print(
        f"\nTEST loci selected: "
        f"{len(selected_loci):,}"
    )

    print(
        f"TEST rows selected: "
        f"{len(work_df):,}"
    )

    print("\nTEST loci:")

    for locus_id in selected_loci:
        print("  ", locus_id)


# ============================================================
# SCORE
# ============================================================

print("\nScoring substitutions...")

score_df = pd.DataFrame(
    [
        score_row(row)
        for _, row in work_df.iterrows()
    ],
    index=work_df.index,
)

result_df = pd.concat(
    [work_df, score_df],
    axis=1,
)

scored = result_df.loc[
    result_df[
        "Substitution_Scoring_Status"
    ]
    == "SCORED"
].copy()


# ============================================================
# QC
# ============================================================

invalid_scored = scored.loc[
    (scored["Mapping_Status"] != "MAPPED")
    | (scored["Difference_Type"] != "MISSENSE")
]

if len(invalid_scored):
    raise ValueError(
        "A scored row is not a mapped missense."
    )

natural_a = scored.loc[
    (scored["Locus_AA"] == "A")
    & (scored["L1RP_AA"] != "A")
]

if len(natural_a):
    if not (
        natural_a[
            "Grantham_WT_to_Natural"
        ]
        ==
        natural_a[
            "Grantham_WT_to_Alanine"
        ]
    ).all():
        raise ValueError(
            "Natural=A Grantham equality QC failed."
        )

    if not (
        natural_a[
            "BLOSUM62_WT_to_Natural"
        ]
        ==
        natural_a[
            "BLOSUM62_WT_to_Alanine"
        ]
    ).all():
        raise ValueError(
            "Natural=A BLOSUM equality QC failed."
        )

wt_a = scored.loc[
    scored["L1RP_AA"] == "A"
]

if len(wt_a):
    if not (
        wt_a[
            "Position_Level_Adney_Comparability"
        ]
        ==
        "ADNEY_DID_NOT_CHANGE_WT_ALANINE_POSITION"
    ).all():
        raise ValueError(
            "WT=A special-case QC failed."
        )


# ============================================================
# CLEAN-ORF PRIORITY SUBSET
# ============================================================

priority_df = result_df.loc[
    (result_df["Mapping_Status"] == "MAPPED")
    & (result_df["Difference_Type"] == "MISSENSE")
    & (
        result_df["Script04_ORF_Status"]
        ==
        "FULL_SPAN_NO_FRAME_STOP_DISRUPTION"
    )
].copy()


# ============================================================
# SUMMARY TABLES
# ============================================================

locus_orf_summary = build_locus_orf_summary(
    result_df
)

scoring_status_counts = (
    result_df[
        "Substitution_Scoring_Status"
    ]
    .value_counts(dropna=False)
    .rename_axis(
        "Substitution_Scoring_Status"
    )
    .reset_index(name="Rows")
)

comparability_counts = (
    scored[
        "Position_Level_Adney_Comparability"
    ]
    .value_counts(dropna=False)
    .rename_axis(
        "Position_Level_Adney_Comparability"
    )
    .reset_index(name="Scored_Rows")
)

direction_counts = (
    scored[
        "Metric_Directional_Agreement"
    ]
    .value_counts(dropna=False)
    .rename_axis(
        "Metric_Directional_Agreement"
    )
    .reset_index(name="Scored_Rows")
)

priority_direction_counts = (
    priority_df[
        "Metric_Directional_Agreement"
    ]
    .value_counts(dropna=False)
    .rename_axis(
        "Metric_Directional_Agreement"
    )
    .reset_index(name="Priority_Rows")
)

priority_subfamily_counts = (
    priority_df[
        "Primary_Subfamily"
    ]
    .value_counts()
    .rename_axis("Primary_Subfamily")
    .reset_index(name="Priority_Rows")
)

natural_equals_alanine_count = int(
    (
        (scored["Locus_AA"] == "A")
        & (scored["L1RP_AA"] != "A")
    ).sum()
)

summary_df = pd.DataFrame(
    [
        ["Run mode", MODE],
        [
            "Script 06 rows analyzed",
            len(result_df),
        ],
        [
            "Rows scored with Grantham and BLOSUM62",
            len(scored),
        ],
        [
            "WT already alanine scored rows",
            int(
                (scored["L1RP_AA"] == "A").sum()
            ),
        ],
        [
            "Natural residue equals Adney alanine rows",
            natural_equals_alanine_count,
        ],
        [
            "Priority clean-ORF mapped missense rows",
            len(priority_df),
        ],
        [
            "Interpretation rule",
            (
                "Grantham/BLOSUM62 are contextual "
                "substitution metrics, not natural-variant "
                "RetroT predictions"
            ),
        ],
        [
            "WT=A rule",
            (
                "If L1-RP WT is alanine, Adney did not "
                "alter that specific residue"
            ),
        ],
    ],
    columns=["Metric", "Value"],
)


# ============================================================
# SAVE
# ============================================================

print("\nSaving outputs...")

result_df.to_csv(
    DETAILED_OUTPUT,
    index=False,
    compression="gzip",
)

priority_df.to_csv(
    PRIORITY_OUTPUT,
    index=False,
    compression="gzip",
)

locus_orf_summary.to_csv(
    LOCUS_ORF_SUMMARY_OUTPUT,
    index=False,
    compression="gzip",
)

if MODE == "TEST":
    result_df.to_csv(
        DETAILED_OUTPUT_PLAIN,
        index=False,
    )

    priority_df.to_csv(
        PRIORITY_OUTPUT_PLAIN,
        index=False,
    )

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

        scoring_status_counts.to_excel(
            writer,
            sheet_name="Scoring_Status",
            index=False,
        )

        comparability_counts.to_excel(
            writer,
            sheet_name="Position_Comparability",
            index=False,
        )

        direction_counts.to_excel(
            writer,
            sheet_name="Metric_Agreement",
            index=False,
        )

        priority_direction_counts.to_excel(
            writer,
            sheet_name="Priority_Metric_Agreement",
            index=False,
        )

        priority_subfamily_counts.to_excel(
            writer,
            sheet_name="Priority_Subfamilies",
            index=False,
        )

        locus_orf_summary.to_excel(
            writer,
            sheet_name="Locus_ORF_Summary",
            index=False,
        )

        if MODE == "TEST":
            result_df.to_excel(
                writer,
                sheet_name="TEST_Details",
                index=False,
            )

            priority_df.to_excel(
                writer,
                sheet_name="TEST_Priority",
                index=False,
            )


# ============================================================
# FINAL REPORT
# ============================================================

print("\n============================================================")
print("SCRIPT 07 COMPLETE")
print("============================================================")

print(
    f"\nRows analyzed: "
    f"{len(result_df):,}"
)

print(
    f"Rows scored: "
    f"{len(scored):,}"
)

print(
    f"WT already alanine rows: "
    f"{int((scored['L1RP_AA'] == 'A').sum()):,}"
)

print(
    "Natural residue equals Adney alanine rows: "
    f"{natural_equals_alanine_count:,}"
)

print(
    f"Priority clean-ORF rows: "
    f"{len(priority_df):,}"
)

print("\nMetric directional agreement:")

print(
    direction_counts.to_string(
        index=False
    )
)

print("\nIMPORTANT:")
print(
    "Grantham and BLOSUM62 are contextual metrics only."
)
print(
    "They do not predict the measured RetroT of a natural "
    "LINE-1 substitution."
)

print("\nOutputs:")
print("  Detailed:", DETAILED_OUTPUT)
print("  Priority:", PRIORITY_OUTPUT)
print(
    "  Locus+ORF summary:",
    LOCUS_ORF_SUMMARY_OUTPUT,
)

if WRITE_EXCEL_SUMMARY:
    print("  Workbook:", SUMMARY_XLSX)

if MODE == "TEST":
    print("\nThis was a TEST run.")
    print('After QC, change RUN_MODE = "TEST"')
    print('to RUN_MODE = "FULL".')
