# ============================================================
# LINE-1 FUNCTIONAL ANNOTATION PROJECT
#
# SCRIPT 04:
# Characterize ORF1 / ORF2 coding-sequence integrity
#
# Input:
#   Script 03 FULL structural-alignment results
#   hg38 LINE-1 FASTA ZIP
#   AF148856.1 GenBank reference
#
# Main purpose:
#   Move from "is this L1-RP region structurally represented?"
#   to "does the represented ORF contain obvious reading-frame
#   or stop-codon disruptions?"
#
# IMPORTANT:
#   - This script does NOT call a LINE-1 active or inactive.
#   - It does NOT yet interpret missense substitutions.
#   - It does NOT treat all differences from L1-RP as damaging.
#   - "NO_FRAME_STOP_DISRUPTION" means only that this script did
#     not detect an obvious start/stop/frameshift problem in the
#     reconstructed ORF sequence.
# ============================================================


# ============================================================
# IMPORTS
# ============================================================

import csv
import io
import re
import time
import zipfile
from collections import Counter
from pathlib import Path

import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq


# ============================================================
#                 EDIT THESE SETTINGS
# ============================================================

# ------------------------------------------------------------
# FIRST RUN:
#
# Keep RUN_MODE = "TEST" for the first run.
# After reviewing the TEST output, change to:
#
#     RUN_MODE = "FULL"
#
# ------------------------------------------------------------

RUN_MODE = "FULL"

# Number of Script 03 loci to process in TEST mode.
TEST_N = 25


# ------------------------------------------------------------
# ALIGNMENT-COMPLEXITY FLAGS
#
# These are QC flags only. They do NOT mean the locus is
# biologically nonfunctional.
# ------------------------------------------------------------

COMPLEX_QUERY_SPAN_PCT = 95.0
LARGE_INDEL_BP = 1000


# ------------------------------------------------------------
# OUTPUT OPTIONS
# ------------------------------------------------------------

# The full 10,178-row metrics table is small enough for Excel
# because this script does NOT place full ORF sequences in the
# spreadsheet.
WRITE_EXCEL_RESULTS = True


# ------------------------------------------------------------
# PROJECT FILENAMES
# ------------------------------------------------------------

FASTA_ZIP_NAME = "hg38_RE_L1_seqs.zip"

STRUCTURE_FOLDER_NAME = "03_LINE1_structure"
STRUCTURE_RESULTS_GZ_NAME = (
    "FULL_LINE1_structure_alignment_results.csv.gz"
)
STRUCTURE_RESULTS_CSV_NAME = (
    "FULL_LINE1_structure_alignment_results.csv"
)

REFERENCE_GB_NAME = "AF148856.1.gb"

OUTPUT_FOLDER_NAME = "04_ORF_integrity"


# ============================================================
#             DO NOT NEED TO EDIT BELOW HERE
# ============================================================

SCRIPT_FOLDER = Path(__file__).resolve().parent
OUTPUT_FOLDER = SCRIPT_FOLDER / OUTPUT_FOLDER_NAME
OUTPUT_FOLDER.mkdir(parents=True, exist_ok=True)

MODE = RUN_MODE.upper().strip()

if MODE not in {"TEST", "FULL"}:
    raise ValueError(
        "\nRUN_MODE must be either 'TEST' or 'FULL'."
    )

PREFIX = MODE


# ============================================================
# FLEXIBLE INPUT FILE RESOLUTION
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


STRUCTURE_RESULTS = first_existing(
    [
        SCRIPT_FOLDER
        / STRUCTURE_FOLDER_NAME
        / STRUCTURE_RESULTS_GZ_NAME,

        SCRIPT_FOLDER
        / STRUCTURE_FOLDER_NAME
        / STRUCTURE_RESULTS_CSV_NAME,

        SCRIPT_FOLDER
        / STRUCTURE_RESULTS_GZ_NAME,

        SCRIPT_FOLDER
        / STRUCTURE_RESULTS_CSV_NAME,
    ],
    "the Script 03 FULL structural results",
)


REFERENCE_GB = first_existing(
    [
        SCRIPT_FOLDER
        / STRUCTURE_FOLDER_NAME
        / REFERENCE_GB_NAME,

        SCRIPT_FOLDER
        / REFERENCE_GB_NAME,
    ],
    "AF148856.1 GenBank reference",
)


FASTA_ZIP = first_existing(
    [
        SCRIPT_FOLDER / FASTA_ZIP_NAME,
    ],
    "the hg38 LINE-1 FASTA ZIP",
)


# ============================================================
# OUTPUT FILES
# ============================================================

RESULTS_CSV = (
    OUTPUT_FOLDER
    / f"{PREFIX}_LINE1_ORF_integrity_results.csv"
)

RESULTS_GZ = (
    OUTPUT_FOLDER
    / f"{PREFIX}_LINE1_ORF_integrity_results.csv.gz"
)

RESULTS_XLSX = (
    OUTPUT_FOLDER
    / f"{PREFIX}_LINE1_ORF_integrity_results.xlsx"
)

SUMMARY_XLSX = (
    OUTPUT_FOLDER
    / f"{PREFIX}_LINE1_ORF_integrity_summary.xlsx"
)


# ============================================================
# HELPER FUNCTIONS
# ============================================================

def reverse_complement(sequence):
    return str(
        Seq(sequence).reverse_complement()
    )


def parse_cigar(cigar_string):
    """
    Convert a SAM-style CIGAR string such as:

        100=1X50=3I20=2D30=

    into:
        [(100, '='), (1, 'X'), ...]
    """

    if pd.isna(cigar_string):
        return []

    return [
        (int(length), op)
        for length, op
        in re.findall(
            r"(\d+)([MIDNSHP=X])",
            str(cigar_string)
        )
    ]


def interval_overlap_length(
    start_a,
    end_a,
    start_b,
    end_b
):
    """
    Length of overlap between two 0-based half-open intervals.
    """
    start = max(start_a, start_b)
    end = min(end_a, end_b)

    return max(
        0,
        end - start
    )


def all_reference_positions_covered(
    covered_intervals,
    positions0
):
    """
    Return True only when every requested 0-based reference
    position is represented by an aligned query base.
    """

    for position in positions0:

        found = False

        for start0, end0 in covered_intervals:
            if start0 <= position < end0:
                found = True
                break

        if not found:
            return False

    return True


def translate_reconstructed_cds(sequence):
    """
    Translate the reconstructed query CDS in the L1-RP reading
    frame starting at its first base.

    Any trailing 1-2 nt are retained in the length metrics but
    excluded from translation.
    """

    sequence = sequence.upper()

    usable_length = (
        len(sequence) // 3
    ) * 3

    if usable_length == 0:
        return ""

    return str(
        Seq(
            sequence[:usable_length]
        ).translate(
            to_stop=False
        )
    )


def analyze_orf_from_cigar(
    oriented_query,
    cigar_string,
    q_start0,
    ref_start0,
    orf_name,
    orf_start_1based,
    orf_end_1based,
):
    """
    Reconstruct the query sequence corresponding to an L1-RP
    ORF using the Script 03 CIGAR alignment.

    Coordinate conventions:
      - q_start0 and ref_start0 are 0-based alignment starts.
      - ORF coordinates are supplied as 1-based inclusive.
      - Internal calculations use 0-based half-open intervals.
    """

    orf_start0 = orf_start_1based - 1
    orf_end0 = orf_end_1based

    reference_cds_length = (
        orf_end0 - orf_start0
    )

    operations = parse_cigar(
        cigar_string
    )

    qpos = int(q_start0)
    rpos = int(ref_start0)

    reconstructed_parts = []
    covered_reference_intervals = []

    mismatch_bp = 0

    insertion_bp = 0
    insertion_events = 0

    deletion_bp = 0
    deletion_events = 0

    frameshift_events = []
    indel_events = []

    ambiguous_base_count = 0


    for length, op in operations:

        # ----------------------------------------------------
        # M / = / X consume both query and reference.
        # ----------------------------------------------------
        if op in ("M", "=", "X"):

            segment_ref_start = rpos
            segment_ref_end = rpos + length

            overlap_start = max(
                segment_ref_start,
                orf_start0
            )

            overlap_end = min(
                segment_ref_end,
                orf_end0
            )

            if overlap_end > overlap_start:

                offset_start = (
                    overlap_start
                    - segment_ref_start
                )

                offset_end = (
                    overlap_end
                    - segment_ref_start
                )

                query_piece = oriented_query[
                    qpos + offset_start:
                    qpos + offset_end
                ]

                reconstructed_parts.append(
                    query_piece
                )

                covered_reference_intervals.append(
                    (
                        overlap_start,
                        overlap_end
                    )
                )

                ambiguous_base_count += sum(
                    base not in {"A", "C", "G", "T"}
                    for base in query_piece.upper()
                )

                if op == "X":
                    mismatch_bp += (
                        overlap_end
                        - overlap_start
                    )

                elif op == "M":
                    # Parasail generally returns = / X in our
                    # Script 03 output. If M appears, it is an
                    # aligned-but-not-explicitly-classified block.
                    pass

            qpos += length
            rpos += length


        # ----------------------------------------------------
        # I consumes query only.
        #
        # The insertion is located between reference positions.
        # We count it as internal to the ORF only when its
        # reference anchor lies strictly inside the ORF span.
        # ----------------------------------------------------
        elif op == "I":

            if (
                orf_start0
                < rpos
                < orf_end0
            ):

                query_piece = oriented_query[
                    qpos:qpos + length
                ]

                reconstructed_parts.append(
                    query_piece
                )

                insertion_bp += length
                insertion_events += 1

                ambiguous_base_count += sum(
                    base not in {"A", "C", "G", "T"}
                    for base in query_piece.upper()
                )

                approx_aa = (
                    (rpos - orf_start0) // 3
                ) + 1

                event = {
                    "type": "I",
                    "length": length,
                    "ref_position0": rpos,
                    "approx_aa": approx_aa,
                    "frameshifting":
                        (length % 3 != 0),
                }

                indel_events.append(
                    event
                )

                if event["frameshifting"]:
                    frameshift_events.append(
                        event
                    )

            qpos += length


        # ----------------------------------------------------
        # D / N consume reference only.
        # ----------------------------------------------------
        elif op in ("D", "N"):

            segment_ref_start = rpos
            segment_ref_end = rpos + length

            overlap_length = interval_overlap_length(
                segment_ref_start,
                segment_ref_end,
                orf_start0,
                orf_end0
            )

            if overlap_length > 0:

                deletion_bp += overlap_length
                deletion_events += 1

                overlap_start = max(
                    segment_ref_start,
                    orf_start0
                )

                approx_aa = (
                    (overlap_start - orf_start0)
                    // 3
                ) + 1

                event = {
                    "type": "D",
                    "length": overlap_length,
                    "ref_position0":
                        overlap_start,
                    "approx_aa": approx_aa,
                    "frameshifting":
                        (overlap_length % 3 != 0),
                }

                indel_events.append(
                    event
                )

                if event["frameshifting"]:
                    frameshift_events.append(
                        event
                    )

            rpos += length


        # ----------------------------------------------------
        # Soft clipping consumes query only.
        # It occurs outside the local alignment and therefore
        # is not inserted into an L1-RP ORF reconstruction.
        # ----------------------------------------------------
        elif op == "S":
            qpos += length


        # H and P consume neither sequence.
        elif op in ("H", "P"):
            pass


    reconstructed_cds = "".join(
        reconstructed_parts
    ).upper()


    # --------------------------------------------------------
    # ORF boundary coverage
    # --------------------------------------------------------

    start_codon_positions0 = list(
        range(
            orf_start0,
            orf_start0 + 3
        )
    )

    terminal_codon_positions0 = list(
        range(
            orf_end0 - 3,
            orf_end0
        )
    )

    start_codon_covered = (
        all_reference_positions_covered(
            covered_reference_intervals,
            start_codon_positions0
        )
    )

    terminal_codon_covered = (
        all_reference_positions_covered(
            covered_reference_intervals,
            terminal_codon_positions0
        )
    )

    full_orf_bounds_covered = (
        start_codon_covered
        and terminal_codon_covered
    )


    # --------------------------------------------------------
    # Translation-level metrics
    # --------------------------------------------------------

    translation = ""

    start_codon = ""
    start_is_atg = False

    terminal_stop_present = False

    premature_stop_count = None
    first_premature_stop_aa = None
    first_stop_aa = None

    predicted_protein_length_aa = None

    cds_length_delta_bp = None
    cds_length_mod3 = None


    if full_orf_bounds_covered:

        cds_length_delta_bp = (
            len(reconstructed_cds)
            - reference_cds_length
        )

        cds_length_mod3 = (
            len(reconstructed_cds)
            % 3
        )

        if len(reconstructed_cds) >= 3:

            start_codon = (
                reconstructed_cds[:3]
            )

            start_is_atg = (
                start_codon == "ATG"
            )


        translation = (
            translate_reconstructed_cds(
                reconstructed_cds
            )
        )


        if translation:

            stop_positions = [
                index + 1
                for index, aa
                in enumerate(translation)
                if aa == "*"
            ]

            if stop_positions:
                first_stop_aa = (
                    stop_positions[0]
                )

            # A proper terminal stop must be the final translated
            # codon AND the CDS length must be divisible by 3.
            terminal_stop_present = (
                cds_length_mod3 == 0
                and translation.endswith("*")
            )

            if terminal_stop_present:

                premature_positions = [
                    position
                    for position in stop_positions
                    if position < len(translation)
                ]

            else:

                premature_positions = (
                    stop_positions
                )


            premature_stop_count = len(
                premature_positions
            )

            if premature_positions:

                first_premature_stop_aa = (
                    premature_positions[0]
                )


            if stop_positions:

                predicted_protein_length_aa = (
                    stop_positions[0] - 1
                )

            else:

                predicted_protein_length_aa = (
                    len(translation)
                )


    # --------------------------------------------------------
    # Frame / stop evidence
    # --------------------------------------------------------

    disruption_reasons = []


    if full_orf_bounds_covered:

        if not start_is_atg:
            disruption_reasons.append(
                "START_CODON_NOT_ATG"
            )

        if cds_length_mod3 != 0:
            disruption_reasons.append(
                "CDS_LENGTH_NOT_DIVISIBLE_BY_3"
            )

        if len(frameshift_events) > 0:
            disruption_reasons.append(
                "FRAME_SHIFTING_INDEL"
            )

        if (
            premature_stop_count is not None
            and premature_stop_count > 0
        ):
            disruption_reasons.append(
                "PREMATURE_STOP"
            )

        if not terminal_stop_present:
            disruption_reasons.append(
                "TERMINAL_STOP_NOT_DETECTED"
            )

        if ambiguous_base_count > 0:
            disruption_reasons.append(
                "AMBIGUOUS_BASES_PRESENT"
            )


        if disruption_reasons:

            sequence_status = (
                "FULL_SPAN_FRAME_STOP_DISRUPTION"
            )

        else:

            sequence_status = (
                "FULL_SPAN_NO_FRAME_STOP_DISRUPTION"
            )


    else:

        sequence_status = (
            "PARTIAL_SPAN_NOT_FULLY_ASSESSED"
        )

        if not start_codon_covered:
            disruption_reasons.append(
                "ORF_START_NOT_FULLY_COVERED"
            )

        if not terminal_codon_covered:
            disruption_reasons.append(
                "ORF_END_NOT_FULLY_COVERED"
            )


    # --------------------------------------------------------
    # Additional indel summaries
    # --------------------------------------------------------

    net_indel_bp = (
        insertion_bp
        - deletion_bp
    )

    net_indel_mod3 = (
        net_indel_bp % 3
    )

    first_frameshift_aa = None

    if frameshift_events:
        first_frameshift_aa = min(
            event["approx_aa"]
            for event in frameshift_events
        )


    indel_event_string = ";".join(
        (
            f"{event['type']}"
            f"{event['length']}"
            f"@aa~{event['approx_aa']}"
        )
        for event in indel_events
    )


    return {
        f"{orf_name}_Sequence_Status":
            sequence_status,

        f"{orf_name}_Start_Codon_Covered":
            start_codon_covered,

        f"{orf_name}_Terminal_Codon_Covered":
            terminal_codon_covered,

        f"{orf_name}_Full_Bounds_Covered":
            full_orf_bounds_covered,

        f"{orf_name}_Reference_CDS_Length_bp":
            reference_cds_length,

        f"{orf_name}_Reconstructed_CDS_Length_bp":
            (
                len(reconstructed_cds)
                if reconstructed_cds
                else 0
            ),

        f"{orf_name}_CDS_Length_Delta_bp":
            cds_length_delta_bp,

        f"{orf_name}_CDS_Length_Mod3":
            cds_length_mod3,

        f"{orf_name}_Start_Codon":
            start_codon,

        f"{orf_name}_Start_Is_ATG":
            (
                start_is_atg
                if full_orf_bounds_covered
                else None
            ),

        f"{orf_name}_Terminal_Stop_Present":
            (
                terminal_stop_present
                if full_orf_bounds_covered
                else None
            ),

        f"{orf_name}_Mismatch_bp_vs_L1RP":
            mismatch_bp,

        f"{orf_name}_Insertion_Events":
            insertion_events,

        f"{orf_name}_Insertion_bp":
            insertion_bp,

        f"{orf_name}_Deletion_Events":
            deletion_events,

        f"{orf_name}_Deletion_bp":
            deletion_bp,

        f"{orf_name}_Frameshifting_Indel_Events":
            len(frameshift_events),

        f"{orf_name}_Net_Indel_bp":
            net_indel_bp,

        f"{orf_name}_Net_Indel_Mod3":
            net_indel_mod3,

        f"{orf_name}_Approx_First_Frameshift_AA":
            first_frameshift_aa,

        f"{orf_name}_Premature_Stop_Count":
            premature_stop_count,

        f"{orf_name}_First_Premature_Stop_AA":
            first_premature_stop_aa,

        f"{orf_name}_First_Stop_AA":
            first_stop_aa,

        f"{orf_name}_Predicted_Protein_Length_aa":
            predicted_protein_length_aa,

        f"{orf_name}_Ambiguous_Base_Count":
            ambiguous_base_count,

        f"{orf_name}_Disruption_Reasons":
            ";".join(
                disruption_reasons
            ),

        f"{orf_name}_Indel_Events":
            indel_event_string,
    }


def classify_alignment_qc(row):
    """
    Carry forward a simple Script 03 alignment-complexity flag.

    This flag is intentionally based on alignment architecture,
    not evolutionary identity, because lower identity can be
    expected for older LINE-1 subfamilies.
    """

    reasons = []

    query_span = float(
        row["Query_Alignment_Span_pct"]
    )

    insertion_bp = int(
        row["Insertion_bp_vs_L1RP"]
    )

    deletion_bp = int(
        row["Deletion_bp_vs_L1RP"]
    )


    if query_span < COMPLEX_QUERY_SPAN_PCT:
        reasons.append(
            f"QUERY_SPAN_LT_{COMPLEX_QUERY_SPAN_PCT:g}PCT"
        )

    if insertion_bp >= LARGE_INDEL_BP:
        reasons.append(
            f"INSERTION_GE_{LARGE_INDEL_BP}BP"
        )

    if deletion_bp >= LARGE_INDEL_BP:
        reasons.append(
            f"DELETION_GE_{LARGE_INDEL_BP}BP"
        )


    if reasons:
        return (
            "COMPLEX",
            ";".join(reasons)
        )

    return (
        "STANDARD",
        ""
    )


def split_reason_counts(series):
    """
    Count semicolon-separated reason strings.
    """

    counter = Counter()

    for value in series.fillna(""):

        for reason in str(value).split(";"):

            reason = reason.strip()

            if reason:
                counter[reason] += 1


    return pd.DataFrame(
        sorted(
            counter.items(),
            key=lambda x: (
                -x[1],
                x[0]
            )
        ),
        columns=[
            "Reason",
            "Locus_Count",
        ]
    )


# ============================================================
# STEP 1 — LOAD SCRIPT 03 RESULTS
# ============================================================

print("\n============================================================")
print("SCRIPT 04 — CHARACTERIZE ORF INTEGRITY")
print("============================================================")

print("\nRun mode:", MODE)

print("\nScript 03 results:")
print(STRUCTURE_RESULTS)

print("\nReference GenBank:")
print(REFERENCE_GB)

print("\nFASTA ZIP:")
print(FASTA_ZIP)


print("\n============================================================")
print("STEP 1 — Loading Script 03 structural results")
print("============================================================")


structure_df = pd.read_csv(
    STRUCTURE_RESULTS
)


required_columns = [
    "Locus_ID",
    "Chromosome",
    "Start_1based",
    "End_1based",
    "Strand",
    "Coordinate_Length",
    "Primary_Subfamily",
    "FASTA_to_L1RP_Orientation",
    "Query_Alignment_Start_1based",
    "Query_Alignment_End_1based",
    "Query_Alignment_Span_pct",
    "L1RP_Alignment_Start_1based",
    "L1RP_Alignment_End_1based",
    "Alignment_Identity_pct",
    "Insertion_bp_vs_L1RP",
    "Deletion_bp_vs_L1RP",
    "CIGAR",
    "ORF1_Coverage_pct",
    "ORF2_Coverage_pct",
]


missing_columns = [
    column
    for column in required_columns
    if column not in structure_df.columns
]


if missing_columns:
    raise ValueError(
        "\nMissing required Script 03 columns:\n"
        + "\n".join(
            missing_columns
        )
    )


print(
    f"Script 03 loci loaded: "
    f"{len(structure_df):,}"
)


if MODE == "TEST":

    work_df = structure_df.head(
        TEST_N
    ).copy()

    print(
        f"TEST mode: processing "
        f"{len(work_df):,} loci."
    )

else:

    work_df = structure_df.copy()

    print(
        f"FULL mode: processing "
        f"{len(work_df):,} loci."
    )


# ============================================================
# STEP 2 — LOAD AF148856.1 AND PARSE ORF COORDINATES
# ============================================================

print("\n============================================================")
print("STEP 2 — Loading L1-RP ORF coordinates")
print("============================================================")


reference_record = SeqIO.read(
    REFERENCE_GB,
    "genbank"
)


reference_sequence = str(
    reference_record.seq
).upper()


orf_coordinates = {}


for feature in reference_record.features:

    if feature.type != "CDS":
        continue

    gene = (
        feature.qualifiers
        .get("gene", [""])[0]
        .upper()
    )

    if gene not in {
        "ORF1",
        "ORF2",
    }:
        continue

    start_1based = (
        int(feature.location.start)
        + 1
    )

    end_1based = int(
        feature.location.end
    )

    translation = (
        feature.qualifiers
        .get("translation", [""])[0]
    )


    orf_coordinates[gene] = {
        "start": start_1based,
        "end": end_1based,
        "protein_length_aa":
            len(translation),
    }


if not {
    "ORF1",
    "ORF2",
}.issubset(
    orf_coordinates
):
    raise ValueError(
        "\nCould not identify both ORF1 and ORF2 "
        "from AF148856.1."
    )


for orf_name in [
    "ORF1",
    "ORF2",
]:

    info = orf_coordinates[
        orf_name
    ]

    cds_length = (
        info["end"]
        - info["start"]
        + 1
    )

    print(
        f"{orf_name}: "
        f"{info['start']:,}-"
        f"{info['end']:,} "
        f"({cds_length:,} bp; "
        f"{info['protein_length_aa']} aa)"
    )


# ============================================================
# STEP 3 — ADD ALIGNMENT QC FLAGS
# ============================================================

print("\n============================================================")
print("STEP 3 — Carrying forward alignment-complexity QC")
print("============================================================")


alignment_qc = work_df.apply(
    classify_alignment_qc,
    axis=1
)


work_df["Alignment_QC"] = (
    alignment_qc.str[0]
)

work_df["Alignment_QC_Reasons"] = (
    alignment_qc.str[1]
)


print(
    "STANDARD alignments:",
    int(
        (
            work_df["Alignment_QC"]
            == "STANDARD"
        ).sum()
    )
)

print(
    "COMPLEX alignments:",
    int(
        (
            work_df["Alignment_QC"]
            == "COMPLEX"
        ).sum()
    )
)


# ============================================================
# STEP 4 — PREPARE LOCUS LOOKUP
# ============================================================

metadata_columns = [
    "Locus_ID",
    "Chromosome",
    "Start_1based",
    "End_1based",
    "Strand",
    "Coordinate_Length",
    "Primary_Subfamily",

    "FASTA_to_L1RP_Orientation",

    "Query_Alignment_Start_1based",
    "Query_Alignment_End_1based",
    "Query_Alignment_Span_pct",

    "L1RP_Alignment_Start_1based",
    "L1RP_Alignment_End_1based",

    "Alignment_Identity_pct",

    "Insertion_bp_vs_L1RP",
    "Deletion_bp_vs_L1RP",

    "ORF1_Coverage_pct",
    "ORF2_Coverage_pct",

    "Alignment_QC",
    "Alignment_QC_Reasons",

    "CIGAR",
]


work_lookup = (
    work_df[
        metadata_columns
    ]
    .set_index(
        "Locus_ID"
    )
    .to_dict(
        orient="index"
    )
)


target_ids = set(
    work_lookup.keys()
)


# ============================================================
# STEP 5 — STREAM FASTA AND ANALYZE ORFs
# ============================================================

print("\n============================================================")
print("STEP 5 — Reconstructing ORF1 and ORF2 from CIGAR")
print("============================================================")

print(
    "\nNo new sequence alignment is being performed."
)

print(
    "Script 04 reuses the Script 03 CIGAR alignment "
    "and should therefore be much faster than Script 03."
)


start_time = time.time()

results = []

processed = 0


def process_fasta_record(
    fasta_id,
    raw_sequence
):
    global processed

    if fasta_id not in target_ids:
        return

    metadata = work_lookup[
        fasta_id
    ]

    sequence = raw_sequence.upper()

    orientation = metadata[
        "FASTA_to_L1RP_Orientation"
    ]


    if orientation == "+":

        oriented_query = sequence

        q_start0 = (
            int(
                metadata[
                    "Query_Alignment_Start_1based"
                ]
            )
            - 1
        )


    elif orientation == "-":

        oriented_query = (
            reverse_complement(
                sequence
            )
        )

        # Script 03 reported query coordinates back in the
        # ORIGINAL FASTA orientation.
        #
        # For a reverse-complement alignment:
        #
        #   q_begin0 =
        #       query_length
        #       - original_query_end_1based
        #
        q_start0 = (
            len(sequence)
            - int(
                metadata[
                    "Query_Alignment_End_1based"
                ]
            )
        )


    else:

        raise ValueError(
            f"Unexpected alignment orientation "
            f"for {fasta_id}: {orientation}"
        )


    ref_start0 = (
        int(
            metadata[
                "L1RP_Alignment_Start_1based"
            ]
        )
        - 1
    )


    output_row = {
        "Locus_ID":
            fasta_id,

        **metadata,
    }


    for orf_name in [
        "ORF1",
        "ORF2",
    ]:

        info = orf_coordinates[
            orf_name
        ]

        metrics = analyze_orf_from_cigar(
            oriented_query=
                oriented_query,

            cigar_string=
                metadata["CIGAR"],

            q_start0=
                q_start0,

            ref_start0=
                ref_start0,

            orf_name=
                orf_name,

            orf_start_1based=
                info["start"],

            orf_end_1based=
                info["end"],
        )

        output_row.update(
            metrics
        )


    results.append(
        output_row
    )

    processed += 1


    if (
        MODE == "TEST"
        or processed % 500 == 0
    ):

        elapsed = (
            time.time()
            - start_time
        )

        print(
            f"  Processed "
            f"{processed:,} / "
            f"{len(target_ids):,} "
            f"loci "
            f"({elapsed/60:.1f} min)"
        )


with zipfile.ZipFile(
    FASTA_ZIP,
    "r"
) as zf:

    fasta_members = [
        name
        for name in zf.namelist()
        if name.lower().endswith(
            (
                ".fa",
                ".fasta",
                ".fna",
                ".fas",
            )
        )
    ]


    if len(fasta_members) != 1:

        raise ValueError(
            "\nExpected exactly one FASTA file "
            "inside the ZIP.\n"
            f"Found: {fasta_members}"
        )


    fasta_member = fasta_members[0]

    print(
        f"\nFASTA inside ZIP: "
        f"{fasta_member}"
    )


    with zf.open(
        fasta_member,
        "r"
    ) as raw_handle:

        handle = io.TextIOWrapper(
            raw_handle,
            encoding="utf-8",
            errors="replace",
        )


        current_id = None
        sequence_parts = []


        for raw_line in handle:

            line = raw_line.strip()

            if not line:
                continue


            if line.startswith(">"):

                if (
                    current_id is not None
                ):

                    process_fasta_record(
                        current_id,
                        "".join(
                            sequence_parts
                        ),
                    )


                current_id = (
                    line[1:]
                    .split()[0]
                )

                sequence_parts = []


            else:

                if current_id in target_ids:
                    sequence_parts.append(
                        line
                    )


        if current_id is not None:

            process_fasta_record(
                current_id,
                "".join(
                    sequence_parts
                ),
            )


# ============================================================
# STEP 6 — VERIFY FASTA PROCESSING
# ============================================================

print("\n============================================================")
print("STEP 6 — Checking FASTA processing")
print("============================================================")


print(
    f"Requested loci: "
    f"{len(target_ids):,}"
)

print(
    f"Processed loci: "
    f"{processed:,}"
)


if processed != len(target_ids):

    found_ids = {
        row["Locus_ID"]
        for row in results
    }

    missing_ids = sorted(
        target_ids
        - found_ids
    )

    raise ValueError(
        "\nNot every Script 03 locus was found in the FASTA.\n"
        f"Missing count: {len(missing_ids):,}\n"
        f"First missing IDs:\n"
        + "\n".join(
            missing_ids[:20]
        )
    )


# ============================================================
# STEP 7 — BUILD RESULTS DATAFRAME
# ============================================================

print("\n============================================================")
print("STEP 7 — Building ORF-integrity results")
print("============================================================")


results_df = pd.DataFrame(
    results
)


# Remove the enormous raw CIGAR from the main user-facing
# Script 04 table. It remains preserved in Script 03.
if "CIGAR" in results_df.columns:
    results_df = results_df.drop(
        columns=["CIGAR"]
    )


# ============================================================
# STEP 8 — SAVE CSV / CSV.GZ
# ============================================================

print("\n============================================================")
print("STEP 8 — Saving ORF-integrity results")
print("============================================================")


results_df.to_csv(
    RESULTS_CSV,
    index=False
)


results_df.to_csv(
    RESULTS_GZ,
    index=False,
    compression="gzip"
)


print(
    f"\nSaved:\n{RESULTS_CSV}"
)

print(
    f"\nSaved compressed copy:\n{RESULTS_GZ}"
)


# ============================================================
# STEP 9 — SUMMARY TABLES
# ============================================================

print("\n============================================================")
print("STEP 9 — Building summary tables")
print("============================================================")


summary_rows = [
    [
        "Run mode",
        MODE,
    ],
    [
        "Loci requested",
        len(target_ids),
    ],
    [
        "Loci processed",
        processed,
    ],
    [
        "Reference",
        "L1-RP / AF148856.1",
    ],
    [
        "Alignment QC query-span threshold (%)",
        COMPLEX_QUERY_SPAN_PCT,
    ],
    [
        "Alignment QC large-indel threshold (bp)",
        LARGE_INDEL_BP,
    ],
    [
        "STANDARD alignments",
        int(
            (
                results_df[
                    "Alignment_QC"
                ]
                == "STANDARD"
            ).sum()
        ),
    ],
    [
        "COMPLEX alignments",
        int(
            (
                results_df[
                    "Alignment_QC"
                ]
                == "COMPLEX"
            ).sum()
        ),
    ],
]


for orf_name in [
    "ORF1",
    "ORF2",
]:

    status_col = (
        f"{orf_name}_Sequence_Status"
    )

    for status, count in (
        results_df[
            status_col
        ]
        .value_counts(
            dropna=False
        )
        .items()
    ):

        summary_rows.append(
            [
                f"{orf_name}: {status}",
                int(count),
            ]
        )


summary_df = pd.DataFrame(
    summary_rows,
    columns=[
        "Metric",
        "Value",
    ]
)


alignment_qc_counts = (
    results_df[
        "Alignment_QC"
    ]
    .value_counts(
        dropna=False
    )
    .rename_axis(
        "Alignment_QC"
    )
    .reset_index(
        name="Locus_Count"
    )
)


orf1_status_counts = (
    results_df[
        "ORF1_Sequence_Status"
    ]
    .value_counts(
        dropna=False
    )
    .rename_axis(
        "ORF1_Sequence_Status"
    )
    .reset_index(
        name="Locus_Count"
    )
)


orf2_status_counts = (
    results_df[
        "ORF2_Sequence_Status"
    ]
    .value_counts(
        dropna=False
    )
    .rename_axis(
        "ORF2_Sequence_Status"
    )
    .reset_index(
        name="Locus_Count"
    )
)


orf1_reason_counts = split_reason_counts(
    results_df[
        "ORF1_Disruption_Reasons"
    ]
)


orf2_reason_counts = split_reason_counts(
    results_df[
        "ORF2_Disruption_Reasons"
    ]
)


reference_orf_rows = []

for orf_name in [
    "ORF1",
    "ORF2",
]:

    info = orf_coordinates[
        orf_name
    ]

    reference_orf_rows.append(
        [
            orf_name,
            info["start"],
            info["end"],
            (
                info["end"]
                - info["start"]
                + 1
            ),
            info["protein_length_aa"],
        ]
    )


reference_orfs_df = pd.DataFrame(
    reference_orf_rows,
    columns=[
        "ORF",
        "L1RP_Start_1based",
        "L1RP_End_1based",
        "Reference_CDS_Length_bp",
        "Reference_Protein_Length_aa",
    ]
)


# Subfamily summary.
subfamily_summary = (
    results_df
    .groupby(
        "Primary_Subfamily",
        dropna=False
    )
    .agg(
        Locus_Count=(
            "Locus_ID",
            "count"
        ),

        Median_Identity_pct=(
            "Alignment_Identity_pct",
            "median"
        ),

        ORF1_FullSpan_NoFrameStopDisruption=(
            "ORF1_Sequence_Status",
            lambda x:
                (
                    x
                    == "FULL_SPAN_NO_FRAME_STOP_DISRUPTION"
                ).sum()
        ),

        ORF1_FullSpan_Disruption=(
            "ORF1_Sequence_Status",
            lambda x:
                (
                    x
                    == "FULL_SPAN_FRAME_STOP_DISRUPTION"
                ).sum()
        ),

        ORF1_Partial=(
            "ORF1_Sequence_Status",
            lambda x:
                (
                    x
                    == "PARTIAL_SPAN_NOT_FULLY_ASSESSED"
                ).sum()
        ),

        ORF2_FullSpan_NoFrameStopDisruption=(
            "ORF2_Sequence_Status",
            lambda x:
                (
                    x
                    == "FULL_SPAN_NO_FRAME_STOP_DISRUPTION"
                ).sum()
        ),

        ORF2_FullSpan_Disruption=(
            "ORF2_Sequence_Status",
            lambda x:
                (
                    x
                    == "FULL_SPAN_FRAME_STOP_DISRUPTION"
                ).sum()
        ),

        ORF2_Partial=(
            "ORF2_Sequence_Status",
            lambda x:
                (
                    x
                    == "PARTIAL_SPAN_NOT_FULLY_ASSESSED"
                ).sum()
        ),
    )
    .reset_index()
)


complex_df = results_df.loc[
    results_df[
        "Alignment_QC"
    ]
    == "COMPLEX"
].copy()


disruption_examples = results_df.loc[
    (
        results_df[
            "ORF1_Sequence_Status"
        ]
        == "FULL_SPAN_FRAME_STOP_DISRUPTION"
    )
    |
    (
        results_df[
            "ORF2_Sequence_Status"
        ]
        == "FULL_SPAN_FRAME_STOP_DISRUPTION"
    )
].head(
    250
).copy()


# ============================================================
# STEP 10 — SAVE SUMMARY WORKBOOK
# ============================================================

print("\n============================================================")
print("STEP 10 — Saving Excel outputs")
print("============================================================")


with pd.ExcelWriter(
    SUMMARY_XLSX,
    engine="openpyxl"
) as writer:

    summary_df.to_excel(
        writer,
        sheet_name="Summary",
        index=False
    )

    reference_orfs_df.to_excel(
        writer,
        sheet_name="Reference_ORFs",
        index=False
    )

    alignment_qc_counts.to_excel(
        writer,
        sheet_name="Alignment_QC",
        index=False
    )

    orf1_status_counts.to_excel(
        writer,
        sheet_name="ORF1_Status",
        index=False
    )

    orf2_status_counts.to_excel(
        writer,
        sheet_name="ORF2_Status",
        index=False
    )

    orf1_reason_counts.to_excel(
        writer,
        sheet_name="ORF1_Reasons",
        index=False
    )

    orf2_reason_counts.to_excel(
        writer,
        sheet_name="ORF2_Reasons",
        index=False
    )

    subfamily_summary.to_excel(
        writer,
        sheet_name="Subfamily_Summary",
        index=False
    )

    complex_df.to_excel(
        writer,
        sheet_name="Complex_Alignments",
        index=False
    )

    disruption_examples.to_excel(
        writer,
        sheet_name="Disruption_Examples",
        index=False
    )


print(
    f"\nSaved summary workbook:\n"
    f"{SUMMARY_XLSX}"
)


if WRITE_EXCEL_RESULTS:

    results_df.to_excel(
        RESULTS_XLSX,
        index=False,
        engine="openpyxl"
    )

    print(
        f"\nSaved complete Excel results:\n"
        f"{RESULTS_XLSX}"
    )


# ============================================================
# FINAL REPORT
# ============================================================

elapsed = (
    time.time()
    - start_time
)


print("\n============================================================")
print("SCRIPT 04 COMPLETE")
print("============================================================")


print(
    f"\nElapsed time: "
    f"{elapsed/60:.2f} minutes"
)


print(
    f"\nLoci processed: "
    f"{processed:,}"
)


for orf_name in [
    "ORF1",
    "ORF2",
]:

    print(
        f"\n{orf_name} sequence-status counts:"
    )

    print(
        results_df[
            f"{orf_name}_Sequence_Status"
        ]
        .value_counts()
        .to_string()
    )


print(
    "\nIMPORTANT INTERPRETATION:"
)

print(
    "\n'FULL_SPAN_NO_FRAME_STOP_DISRUPTION' does NOT mean "
    "that the LINE-1 is active."
)

print(
    "It means that, within this L1-RP-aligned ORF span, "
    "Script 04 did not detect an obvious start-codon, "
    "premature-stop, terminal-stop, or frameshifting-indel "
    "problem."
)

print(
    "\nMissense substitutions and their experimental "
    "relevance are intentionally left for the next pipeline "
    "steps."
)


if MODE == "TEST":

    print(
        "\nThis was a TEST run."
    )

    print(
        "Review the TEST summary/results before changing:"
    )

    print(
        '\n    RUN_MODE = "FULL"\n'
    )
