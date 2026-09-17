# ============================================================
# LINE-1 FUNCTIONAL ANNOTATION PROJECT
#
# SCRIPT 05 — VERSION 3
# ORF-specific semi-global realignment + exact AA differences
# relative to L1-RP / AF148856.1
#
# WHY V3 EXISTS
# ------------------------------------------------------------
# Two earlier Script 05 test versions failed because Parasail's
# semi-global traceback/CIGAR start-coordinate handling was not
# reliable enough for our validation logic.
#
# V3 removes Parasail from Script 05 entirely.
#
# It uses Biopython PairwiseAligner with:
#   - target  = L1-RP ORF1 or ORF2
#   - query   = genomic LINE-1 ORF-containing sequence
#   - global DP
#   - FREE query overhangs at both ends
#
# This gives us a true semi-global alignment:
#
#       extra genomic flank  ORF-like sequence  extra flank
#       ------------------[===================]------------
#                           | entire L1-RP ORF |
#
# The entire L1-RP ORF is retained as the target, while extra
# genomic bases on either side are not penalized.
#
# Exact CIGAR strings are then built directly from Biopython's
# alignment coordinate path. We do NOT infer traceback start
# coordinates from another representation.
#
# IMPORTANT INTERPRETATION
# ------------------------------------------------------------
# "AA difference vs L1-RP" does NOT automatically mean a
# locus-specific mutation. For older LINE-1 subfamilies, many
# differences are expected evolutionary/subfamily differences.
#
# Script 05 does NOT map differences to Adney yet.
# ============================================================


# ============================================================
# IMPORTS
# ============================================================

import csv
import gzip
import io
import re
import time
import zipfile
from pathlib import Path

import numpy as np
import pandas as pd
from Bio import Align, SeqIO
from Bio.Seq import Seq


# ============================================================
#                 EDIT THESE SETTINGS
# ============================================================

RUN_MODE = "FULL"
TEST_N = 25

# These are always included in TEST mode when present.
TEST_LOCUS_IDS = [
    "hg38.RE.chr1.100199603-100206089.+",
    "hg38.RE.chr1.104770248-104776279.-",
    "hg38.RE.chr1.104843834-104849865.-",
]

# Script 03 is used only to estimate where each ORF lies in the
# genomic locus. We then add flanking sequence and REALIGN.
ORF_QUERY_FLANK_BP = 300

# Only ORFs for which Script 04 said both ORF boundaries were
# represented are sent to exact ORF-specific realignment.
REQUIRE_SCRIPT04_FULL_BOUNDS = True

# Detailed full-run tables are compressed because they may be
# large. TEST mode also writes plain CSV copies.
WRITE_TEST_PLAIN_CSV = True


# ============================================================
# PROJECT FILENAMES / FOLDERS
# ============================================================

FASTA_ZIP_NAME = "hg38_RE_L1_seqs.zip"

STRUCTURE_FOLDER_NAME = "03_LINE1_structure"
STRUCTURE_RESULTS_GZ_NAME = (
    "FULL_LINE1_structure_alignment_results.csv.gz"
)
STRUCTURE_RESULTS_CSV_NAME = (
    "FULL_LINE1_structure_alignment_results.csv"
)

INTEGRITY_FOLDER_NAME = "04_ORF_integrity"
INTEGRITY_RESULTS_GZ_NAME = (
    "FULL_LINE1_ORF_integrity_results.csv.gz"
)
INTEGRITY_RESULTS_CSV_NAME = (
    "FULL_LINE1_ORF_integrity_results.csv"
)

REFERENCE_GB_NAME = "AF148856.1.gb"

OUTPUT_FOLDER_NAME = "05_ORF_variant_alignment"


# ============================================================
# SETUP
# ============================================================

SCRIPT_FOLDER = Path(__file__).resolve().parent
OUTPUT_FOLDER = SCRIPT_FOLDER / OUTPUT_FOLDER_NAME
OUTPUT_FOLDER.mkdir(parents=True, exist_ok=True)

MODE = RUN_MODE.upper().strip()

if MODE not in {"TEST", "FULL"}:
    raise ValueError("RUN_MODE must be 'TEST' or 'FULL'.")

PREFIX = MODE


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
    "Script 03 FULL structural results",
)


INTEGRITY_RESULTS = first_existing(
    [
        SCRIPT_FOLDER
        / INTEGRITY_FOLDER_NAME
        / INTEGRITY_RESULTS_GZ_NAME,

        SCRIPT_FOLDER
        / INTEGRITY_FOLDER_NAME
        / INTEGRITY_RESULTS_CSV_NAME,

        SCRIPT_FOLDER
        / INTEGRITY_RESULTS_GZ_NAME,

        SCRIPT_FOLDER
        / INTEGRITY_RESULTS_CSV_NAME,
    ],
    "Script 04 FULL ORF-integrity results",
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
# OUTPUT PATHS
# ============================================================

ALIGNMENT_RESULTS_GZ = (
    OUTPUT_FOLDER
    / f"{PREFIX}_LINE1_ORF_specific_alignment_results.csv.gz"
)

AA_DIFFERENCES_GZ = (
    OUTPUT_FOLDER
    / f"{PREFIX}_LINE1_AA_differences_vs_L1RP.csv.gz"
)

INDEL_EVENTS_GZ = (
    OUTPUT_FOLDER
    / f"{PREFIX}_LINE1_ORF_indel_events_vs_L1RP.csv.gz"
)

SUMMARY_XLSX = (
    OUTPUT_FOLDER
    / f"{PREFIX}_LINE1_ORF_variant_alignment_summary.xlsx"
)

if MODE == "TEST":
    ALIGNMENT_RESULTS_CSV = (
        OUTPUT_FOLDER
        / "TEST_LINE1_ORF_specific_alignment_results.csv"
    )

    AA_DIFFERENCES_CSV = (
        OUTPUT_FOLDER
        / "TEST_LINE1_AA_differences_vs_L1RP.csv"
    )

    INDEL_EVENTS_CSV = (
        OUTPUT_FOLDER
        / "TEST_LINE1_ORF_indel_events_vs_L1RP.csv"
    )


# ============================================================
# HELPERS
# ============================================================

def clean_dna(sequence):
    sequence = str(sequence).upper()

    return "".join(
        base if base in {"A", "C", "G", "T"}
        else "N"
        for base in sequence
    )


def reverse_complement(sequence):
    return str(
        Seq(sequence).reverse_complement()
    )


def as_bool(value):
    if isinstance(value, bool):
        return value

    if pd.isna(value):
        return False

    return str(value).strip().lower() in {
        "true",
        "1",
        "yes",
        "y",
    }


def parse_cigar(cigar_string):
    if pd.isna(cigar_string):
        return []

    return [
        (int(length), op)
        for length, op
        in re.findall(
            r"(\d+)([MIDNSHP=X])",
            str(cigar_string),
        )
    ]


def add_op(ops, length, op):
    length = int(length)

    if length <= 0:
        return

    if ops and ops[-1][1] == op:
        ops[-1] = (
            ops[-1][0] + length,
            op,
        )

    else:
        ops.append(
            (
                length,
                op,
            )
        )


def ops_to_string(ops):
    return "".join(
        f"{length}{op}"
        for length, op
        in ops
    )


def extended_to_samstyle(ops):
    merged = []

    for length, op in ops:
        sam_op = (
            "M"
            if op in {"=", "X", "M"}
            else op
        )

        add_op(
            merged,
            length,
            sam_op,
        )

    return ops_to_string(
        merged
    )


def oriented_query_start0(row, sequence_length):
    orientation = row[
        "FASTA_to_L1RP_Orientation"
    ]

    if orientation == "+":
        return (
            int(
                row[
                    "Query_Alignment_Start_1based"
                ]
            )
            - 1
        )

    if orientation == "-":
        return (
            sequence_length
            - int(
                row[
                    "Query_Alignment_End_1based"
                ]
            )
        )

    raise ValueError(
        f"Unexpected alignment orientation: "
        f"{orientation}"
    )


def locate_orf_in_query_from_script03(
    cigar_ops,
    whole_query_start0,
    whole_ref_start0,
    orf_ref_start0,
    orf_ref_end0,
):
    """
    Use Script 03 ONLY to locate an approximate ORF-containing
    interval in the oriented genomic LINE-1 sequence.

    Script 03 CIGAR semantics in this pipeline:
        I = query/locus consumes bases
        D = reference/L1-RP consumes bases

    Returns:
        (query_start0, query_end0_inclusive)
    """

    qpos = whole_query_start0
    rpos = whole_ref_start0

    query_positions = []


    for length, op in cigar_ops:

        if op in {"M", "=", "X"}:

            ref_block_start = rpos
            ref_block_end = (
                rpos + length
            )

            overlap_start = max(
                ref_block_start,
                orf_ref_start0,
            )

            overlap_end = min(
                ref_block_end,
                orf_ref_end0,
            )


            if overlap_end > overlap_start:

                q_offset_start = (
                    overlap_start
                    - ref_block_start
                )

                q_offset_end = (
                    overlap_end
                    - ref_block_start
                )

                query_positions.extend(
                    [
                        qpos + q_offset_start,
                        qpos + q_offset_end - 1,
                    ]
                )


            qpos += length
            rpos += length


        elif op == "I":

            # Query/locus-only insertion.
            if (
                orf_ref_start0
                <= rpos
                <= orf_ref_end0
            ):

                query_positions.extend(
                    [
                        qpos,
                        qpos + length - 1,
                    ]
                )

            qpos += length


        elif op in {"D", "N"}:

            rpos += length


        elif op == "S":

            qpos += length


        elif op in {"H", "P"}:

            pass


    if not query_positions:
        return None


    return (
        min(query_positions),
        max(query_positions),
    )


def alignment_to_core_ops(
    alignment,
    reference_cds,
    candidate_sequence,
):
    """
    Convert Biopython PairwiseAligner coordinates directly into
    an extended CIGAR.

    target = L1-RP reference ORF
    query  = candidate genomic sequence

    Coordinate movement:
      target + query move -> aligned bases -> = or X
      target stays, query moves -> I
      target moves, query stays -> D

    Terminal query-only overhangs are FREE in our semi-global
    scoring. They are returned as clip counts and are NOT kept
    as biological ORF insertions.
    """

    coords = np.asarray(
        alignment.coordinates,
        dtype=int,
    )

    raw_ops = []

    segments = []


    for i in range(
        coords.shape[1] - 1
    ):

        t0 = int(
            coords[0, i]
        )

        t1 = int(
            coords[0, i + 1]
        )

        q0 = int(
            coords[1, i]
        )

        q1 = int(
            coords[1, i + 1]
        )

        dt = (
            t1 - t0
        )

        dq = (
            q1 - q0
        )


        segments.append(
            (
                t0,
                t1,
                q0,
                q1,
            )
        )


        if (
            dt > 0
            and dq > 0
        ):

            if dt != dq:
                raise ValueError(
                    "Aligned coordinate block had unequal "
                    "reference/query movement."
                )


            for offset in range(
                dt
            ):

                ref_base = (
                    reference_cds[
                        t0 + offset
                    ]
                )

                query_base = (
                    candidate_sequence[
                        q0 + offset
                    ]
                )

                op = (
                    "="
                    if ref_base == query_base
                    else "X"
                )

                add_op(
                    raw_ops,
                    1,
                    op,
                )


        elif (
            dt == 0
            and dq > 0
        ):

            add_op(
                raw_ops,
                dq,
                "I",
            )


        elif (
            dt > 0
            and dq == 0
        ):

            add_op(
                raw_ops,
                dt,
                "D",
            )


        else:

            raise ValueError(
                "Unexpected zero-length alignment segment."
            )


    # --------------------------------------------------------
    # Remove FREE terminal query overhangs from the core CIGAR.
    # --------------------------------------------------------

    left_free_query_bp = 0
    right_free_query_bp = 0

    core_ops = list(
        raw_ops
    )


    if (
        core_ops
        and core_ops[0][1] == "I"
        and coords[0, 0] == 0
    ):

        left_free_query_bp = (
            core_ops[0][0]
        )

        core_ops = (
            core_ops[1:]
        )


    if (
        core_ops
        and core_ops[-1][1] == "I"
        and coords[0, -1]
        == len(reference_cds)
    ):

        right_free_query_bp = (
            core_ops[-1][0]
        )

        core_ops = (
            core_ops[:-1]
        )


    # Actual aligned query portion inside candidate after free
    # overhangs are removed.
    query_core_start0 = (
        int(coords[1, 0])
        + left_free_query_bp
    )

    query_core_end0_exclusive = (
        int(coords[1, -1])
        - right_free_query_bp
    )


    return {
        "core_ops":
            core_ops,

        "extended_cigar":
            ops_to_string(
                core_ops
            ),

        "samstyle_cigar":
            extended_to_samstyle(
                core_ops
            ),

        "left_free_query_bp":
            left_free_query_bp,

        "right_free_query_bp":
            right_free_query_bp,

        "query_core_start0":
            query_core_start0,

        "query_core_end0_exclusive":
            query_core_end0_exclusive,

        "reference_start0":
            int(
                coords[0, 0]
            ),

        "reference_end0_exclusive":
            int(
                coords[0, -1]
            ),
    }


def build_alignment_maps(
    reference_cds,
    candidate_sequence,
    core_ops,
    query_core_start0,
):
    """
    Build exact reference-position -> genomic-query-base mapping
    from the core ORF-specific CIGAR.
    """

    ref_to_query_base = [
        None
    ] * len(reference_cds)

    ref_to_query_pos = [
        None
    ] * len(reference_cds)

    insertions = []
    deletions = []

    qpos = int(
        query_core_start0
    )

    rpos = 0

    exact_matches = 0
    mismatches = 0
    aligned_pair_bp = 0

    insertion_bp = 0
    deletion_bp = 0

    largest_insertion_bp = 0
    largest_deletion_bp = 0


    for length, op in core_ops:

        if op in {
            "=",
            "X",
            "M",
        }:

            q_piece = (
                candidate_sequence[
                    qpos:qpos + length
                ]
            )

            r_piece = (
                reference_cds[
                    rpos:rpos + length
                ]
            )


            for offset in range(
                length
            ):

                if (
                    rpos + offset
                    < len(reference_cds)
                ):

                    ref_to_query_base[
                        rpos + offset
                    ] = (
                        q_piece[
                            offset
                        ]
                    )

                    ref_to_query_pos[
                        rpos + offset
                    ] = (
                        qpos
                        + offset
                    )


            if op == "=":

                exact_matches += (
                    length
                )

            elif op == "X":

                mismatches += (
                    length
                )

            else:

                for qbase, rbase in zip(
                    q_piece,
                    r_piece,
                ):

                    if qbase == rbase:
                        exact_matches += 1
                    else:
                        mismatches += 1


            aligned_pair_bp += (
                length
            )

            qpos += length
            rpos += length


        elif op == "I":

            inserted_sequence = (
                candidate_sequence[
                    qpos:qpos + length
                ]
            )

            insertions.append(
                {
                    "ref_anchor0":
                        rpos,

                    "length_bp":
                        length,

                    "sequence":
                        inserted_sequence,
                }
            )

            insertion_bp += (
                length
            )

            largest_insertion_bp = max(
                largest_insertion_bp,
                length,
            )

            qpos += length


        elif op == "D":

            deleted_sequence = (
                reference_cds[
                    rpos:rpos + length
                ]
            )

            deletions.append(
                {
                    "ref_start0":
                        rpos,

                    "length_bp":
                        length,

                    "sequence":
                        deleted_sequence,
                }
            )

            deletion_bp += (
                length
            )

            largest_deletion_bp = max(
                largest_deletion_bp,
                length,
            )

            rpos += length


        else:

            raise ValueError(
                f"Unexpected core CIGAR op: {op}"
            )


    if rpos != len(
        reference_cds
    ):

        raise ValueError(
            "Core ORF-specific CIGAR did not consume "
            "the full L1-RP ORF reference."
        )


    identity_pct = (
        None
    )

    if aligned_pair_bp > 0:

        identity_pct = round(
            100.0
            * exact_matches
            / aligned_pair_bp,
            3,
        )


    return {
        "ref_to_query_base":
            ref_to_query_base,

        "ref_to_query_pos":
            ref_to_query_pos,

        "insertions":
            insertions,

        "deletions":
            deletions,

        "exact_matches":
            exact_matches,

        "mismatches":
            mismatches,

        "aligned_pair_bp":
            aligned_pair_bp,

        "identity_pct":
            identity_pct,

        "insertion_bp":
            insertion_bp,

        "deletion_bp":
            deletion_bp,

        "largest_insertion_bp":
            largest_insertion_bp,

        "largest_deletion_bp":
            largest_deletion_bp,
    }


def translate_codon(codon):
    codon = str(
        codon
    ).upper()

    if (
        len(codon) != 3
        or any(
            base not in {
                "A",
                "C",
                "G",
                "T",
            }
            for base in codon
        )
    ):

        return None


    return str(
        Seq(codon).translate(
            to_stop=False
        )
    )


def insertion_anchor_label(
    ref_anchor0,
    ref_length,
):
    if ref_anchor0 <= 0:
        return "before_nt1"

    if ref_anchor0 >= ref_length:
        return (
            f"after_nt{ref_length}"
        )

    return (
        f"between_nt{ref_anchor0}"
        f"_and_nt{ref_anchor0 + 1}"
    )


def analyze_aa_differences(
    locus_id,
    row,
    orf_name,
    reference_cds,
    reference_protein_length,
    maps,
    identity_pct,
):
    """
    Identify simple AA differences vs L1-RP.

    We stop exact AA calling after the first frameshifting indel
    or after a newly introduced nonsense codon.
    """

    ref_to_query_base = (
        maps[
            "ref_to_query_base"
        ]
    )

    ref_to_query_pos = (
        maps[
            "ref_to_query_pos"
        ]
    )

    insertions = (
        maps[
            "insertions"
        ]
    )

    deletions = (
        maps[
            "deletions"
        ]
    )


    # --------------------------------------------------------
    # First frameshift position.
    # --------------------------------------------------------

    frameshift_ref_positions = []


    for event in insertions:

        if (
            event[
                "length_bp"
            ]
            % 3
            != 0
        ):

            frameshift_ref_positions.append(
                event[
                    "ref_anchor0"
                ]
            )


    for event in deletions:

        if (
            event[
                "length_bp"
            ]
            % 3
            != 0
        ):

            frameshift_ref_positions.append(
                event[
                    "ref_start0"
                ]
            )


    first_frameshift_aa = (
        None
    )

    if frameshift_ref_positions:

        first_frameshift_ref0 = min(
            frameshift_ref_positions
        )

        first_frameshift_aa = (
            first_frameshift_ref0
            // 3
        ) + 1


    insertion_anchors = {
        event[
            "ref_anchor0"
        ]
        for event in insertions
    }


    aa_rows = []

    synonymous_count = 0
    nonsynonymous_count = 0
    missense_count = 0
    nonsense_count = 0
    start_lost_count = 0

    uncallable_indel_codons = 0
    uncallable_ambiguous_codons = 0

    first_nonsense_aa = (
        None
    )

    callable_through = (
        reference_protein_length
    )

    call_stop_reason = (
        "REFERENCE_PROTEIN_END"
    )


    if (
        first_frameshift_aa
        is not None
    ):

        callable_through = min(
            callable_through,
            first_frameshift_aa - 1,
        )

        call_stop_reason = (
            "FIRST_FRAMESHIFT"
        )


    # Protein residues only; terminal stop handled separately.
    for aa_position in range(
        1,
        reference_protein_length + 1,
    ):

        if (
            aa_position
            > callable_through
        ):

            break


        codon_start0 = (
            (aa_position - 1)
            * 3
        )

        codon_end0 = (
            codon_start0
            + 3
        )

        ref_codon = (
            reference_cds[
                codon_start0:
                codon_end0
            ]
        )

        ref_aa = (
            translate_codon(
                ref_codon
            )
        )


        query_bases = (
            ref_to_query_base[
                codon_start0:
                codon_end0
            ]
        )

        query_positions = (
            ref_to_query_pos[
                codon_start0:
                codon_end0
            ]
        )


        insertion_inside_codon = any(
            anchor in {
                codon_start0 + 1,
                codon_start0 + 2,
            }
            for anchor
            in insertion_anchors
        )


        if (
            len(query_bases) != 3
            or any(
                base is None
                for base
                in query_bases
            )
            or insertion_inside_codon
        ):

            uncallable_indel_codons += 1
            continue


        if any(
            position is None
            for position
            in query_positions
        ):

            uncallable_indel_codons += 1
            continue


        if not (
            query_positions[1]
            == query_positions[0] + 1
            and query_positions[2]
            == query_positions[1] + 1
        ):

            uncallable_indel_codons += 1
            continue


        query_codon = "".join(
            query_bases
        )

        query_aa = (
            translate_codon(
                query_codon
            )
        )


        if (
            ref_aa is None
            or query_aa is None
        ):

            uncallable_ambiguous_codons += 1
            continue


        if query_codon == ref_codon:
            continue


        nt_changes = []

        for offset, (
            ref_base,
            query_base,
        ) in enumerate(
            zip(
                ref_codon,
                query_codon,
            )
        ):

            if ref_base != query_base:

                nt_position = (
                    codon_start0
                    + offset
                    + 1
                )

                nt_changes.append(
                    f"nt{nt_position}:"
                    f"{ref_base}>{query_base}"
                )


        if query_aa == ref_aa:

            synonymous_count += 1
            continue


        nonsynonymous_count += 1


        if (
            aa_position == 1
            and ref_aa == "M"
            and query_aa != "M"
            and query_aa != "*"
        ):

            difference_type = (
                "START_LOST"
            )

            start_lost_count += 1


        elif query_aa == "*":

            difference_type = (
                "NONSENSE"
            )

            nonsense_count += 1


        else:

            difference_type = (
                "MISSENSE"
            )

            missense_count += 1


        aa_rows.append(
            {
                "Locus_ID":
                    locus_id,

                "Chromosome":
                    row[
                        "Chromosome"
                    ],

                "Start_1based":
                    row[
                        "Start_1based"
                    ],

                "End_1based":
                    row[
                        "End_1based"
                    ],

                "Strand":
                    row[
                        "Strand"
                    ],

                "Primary_Subfamily":
                    row[
                        "Primary_Subfamily"
                    ],

                "ORF":
                    orf_name,

                "AA_Position":
                    aa_position,

                "L1RP_AA":
                    ref_aa,

                "Locus_AA":
                    query_aa,

                "AA_Difference_vs_L1RP":
                    (
                        f"{ref_aa}"
                        f"{aa_position}"
                        f"{query_aa}"
                    ),

                "Difference_Type":
                    difference_type,

                "L1RP_Codon":
                    ref_codon,

                "Locus_Codon":
                    query_codon,

                "Nucleotide_Differences_In_Codon":
                    ";".join(
                        nt_changes
                    ),

                "Script04_ORF_Status":
                    row[
                        f"{orf_name}"
                        "_Sequence_Status"
                    ],

                "Script03_Alignment_QC":
                    row[
                        "Alignment_QC"
                    ],

                "ORF_Specific_Identity_pct":
                    identity_pct,

                "Interpretation":
                    (
                        "DIFFERENCE_VS_L1RP_"
                        "NOT_AUTOMATICALLY_A_"
                        "LOCUS_SPECIFIC_MUTATION"
                    ),
            }
        )


        if (
            difference_type
            == "NONSENSE"
        ):

            first_nonsense_aa = (
                aa_position
            )

            callable_through = (
                aa_position
            )

            call_stop_reason = (
                "FIRST_NONSENSE"
            )

            break


    # --------------------------------------------------------
    # Terminal stop.
    # --------------------------------------------------------

    terminal_stop_status = (
        "NOT_ASSESSED"
    )

    stop_start0 = (
        reference_protein_length
        * 3
    )


    if (
        stop_start0 + 3
        <= len(reference_cds)
        and first_frameshift_aa
        is None
        and first_nonsense_aa
        is None
    ):

        stop_bases = (
            ref_to_query_base[
                stop_start0:
                stop_start0 + 3
            ]
        )

        stop_positions = (
            ref_to_query_pos[
                stop_start0:
                stop_start0 + 3
            ]
        )

        insertion_inside_stop = any(
            anchor in {
                stop_start0 + 1,
                stop_start0 + 2,
            }
            for anchor
            in insertion_anchors
        )


        if (
            len(stop_bases) == 3
            and all(
                base is not None
                for base
                in stop_bases
            )
            and all(
                position is not None
                for position
                in stop_positions
            )
            and not insertion_inside_stop
            and stop_positions[1]
                == stop_positions[0] + 1
            and stop_positions[2]
                == stop_positions[1] + 1
        ):

            query_stop_codon = (
                "".join(
                    stop_bases
                )
            )

            query_stop_aa = (
                translate_codon(
                    query_stop_codon
                )
            )


            if query_stop_aa == "*":

                terminal_stop_status = (
                    "STOP_RETAINED"
                )


            elif query_stop_aa is None:

                terminal_stop_status = (
                    "AMBIGUOUS"
                )


            else:

                terminal_stop_status = (
                    "STOP_LOST"
                )


                ref_stop_codon = (
                    reference_cds[
                        stop_start0:
                        stop_start0 + 3
                    ]
                )


                aa_rows.append(
                    {
                        "Locus_ID":
                            locus_id,

                        "Chromosome":
                            row[
                                "Chromosome"
                            ],

                        "Start_1based":
                            row[
                                "Start_1based"
                            ],

                        "End_1based":
                            row[
                                "End_1based"
                            ],

                        "Strand":
                            row[
                                "Strand"
                            ],

                        "Primary_Subfamily":
                            row[
                                "Primary_Subfamily"
                            ],

                        "ORF":
                            orf_name,

                        "AA_Position":
                            (
                                reference_protein_length
                                + 1
                            ),

                        "L1RP_AA":
                            "*",

                        "Locus_AA":
                            query_stop_aa,

                        "AA_Difference_vs_L1RP":
                            (
                                "*"
                                f"{reference_protein_length + 1}"
                                f"{query_stop_aa}"
                            ),

                        "Difference_Type":
                            "STOP_LOST",

                        "L1RP_Codon":
                            ref_stop_codon,

                        "Locus_Codon":
                            query_stop_codon,

                        "Nucleotide_Differences_In_Codon":
                            "",

                        "Script04_ORF_Status":
                            row[
                                f"{orf_name}"
                                "_Sequence_Status"
                            ],

                        "Script03_Alignment_QC":
                            row[
                                "Alignment_QC"
                            ],

                        "ORF_Specific_Identity_pct":
                            identity_pct,

                        "Interpretation":
                            (
                                "TERMINAL_STOP_"
                                "DIFFERENCE_VS_L1RP"
                            ),
                    }
                )


        else:

            terminal_stop_status = (
                "NOT_SIMPLE_CODON_MAPPING"
            )


    return {
        "aa_rows":
            aa_rows,

        "Synonymous_Codon_Differences":
            synonymous_count,

        "Nonsynonymous_AA_Differences":
            nonsynonymous_count,

        "Missense_AA_Differences":
            missense_count,

        "Nonsense_AA_Differences":
            nonsense_count,

        "Start_Lost_Differences":
            start_lost_count,

        "First_Frameshift_AA":
            first_frameshift_aa,

        "First_Nonsense_AA":
            first_nonsense_aa,

        "AA_Callable_Through_AA":
            max(
                0,
                callable_through,
            ),

        "AA_Call_Stop_Reason":
            call_stop_reason,

        "Terminal_Stop_Status":
            terminal_stop_status,

        "Uncallable_Indel_Affected_Codons":
            uncallable_indel_codons,

        "Uncallable_Ambiguous_Codons":
            uncallable_ambiguous_codons,
    }


def build_indel_rows(
    locus_id,
    row,
    orf_name,
    reference_cds,
    maps,
    identity_pct,
):
    rows = []


    for event in (
        maps[
            "insertions"
        ]
    ):

        ref_anchor0 = (
            event[
                "ref_anchor0"
            ]
        )

        length = (
            event[
                "length_bp"
            ]
        )


        rows.append(
            {
                "Locus_ID":
                    locus_id,

                "Primary_Subfamily":
                    row[
                        "Primary_Subfamily"
                    ],

                "ORF":
                    orf_name,

                "Indel_Type":
                    (
                        "INSERTION_IN_LOCUS_"
                        "VS_L1RP"
                    ),

                "Length_bp":
                    length,

                "Frame_Effect":
                    (
                        "IN_FRAME"
                        if length % 3 == 0
                        else "FRAMESHIFT"
                    ),

                "Reference_Location":
                    insertion_anchor_label(
                        ref_anchor0,
                        len(reference_cds),
                    ),

                "Approx_AA_Position":
                    (
                        ref_anchor0
                        // 3
                    ) + 1,

                "Sequence":
                    event[
                        "sequence"
                    ],

                "Script04_ORF_Status":
                    row[
                        f"{orf_name}"
                        "_Sequence_Status"
                    ],

                "Script03_Alignment_QC":
                    row[
                        "Alignment_QC"
                    ],

                "ORF_Specific_Identity_pct":
                    identity_pct,

                "Interpretation":
                    (
                        "INDEL_VS_L1RP_"
                        "NOT_AUTOMATICALLY_A_"
                        "LOCUS_SPECIFIC_MUTATION"
                    ),
            }
        )


    for event in (
        maps[
            "deletions"
        ]
    ):

        ref_start0 = (
            event[
                "ref_start0"
            ]
        )

        length = (
            event[
                "length_bp"
            ]
        )


        rows.append(
            {
                "Locus_ID":
                    locus_id,

                "Primary_Subfamily":
                    row[
                        "Primary_Subfamily"
                    ],

                "ORF":
                    orf_name,

                "Indel_Type":
                    (
                        "DELETION_IN_LOCUS_"
                        "VS_L1RP"
                    ),

                "Length_bp":
                    length,

                "Frame_Effect":
                    (
                        "IN_FRAME"
                        if length % 3 == 0
                        else "FRAMESHIFT"
                    ),

                "Reference_Location":
                    (
                        f"nt{ref_start0 + 1}"
                        f"-"
                        f"nt{ref_start0 + length}"
                    ),

                "Approx_AA_Position":
                    (
                        ref_start0
                        // 3
                    ) + 1,

                "Sequence":
                    event[
                        "sequence"
                    ],

                "Script04_ORF_Status":
                    row[
                        f"{orf_name}"
                        "_Sequence_Status"
                    ],

                "Script03_Alignment_QC":
                    row[
                        "Alignment_QC"
                    ],

                "ORF_Specific_Identity_pct":
                    identity_pct,

                "Interpretation":
                    (
                        "INDEL_VS_L1RP_"
                        "NOT_AUTOMATICALLY_A_"
                        "LOCUS_SPECIFIC_MUTATION"
                    ),
            }
        )


    return rows


def oriented_to_original_bounds(
    oriented_start0,
    oriented_end0_exclusive,
    sequence_length,
    orientation,
):
    """
    Convert oriented-query interval to 1-based inclusive
    coordinates in the original FASTA record.
    """

    if orientation == "+":

        return (
            oriented_start0 + 1,
            oriented_end0_exclusive,
        )


    if orientation == "-":

        original_start_1based = (
            sequence_length
            - oriented_end0_exclusive
            + 1
        )

        original_end_1based = (
            sequence_length
            - oriented_start0
        )

        return (
            original_start_1based,
            original_end_1based,
        )


    raise ValueError(
        f"Unexpected orientation: {orientation}"
    )


# ============================================================
# OUTPUT FIELDS
# ============================================================

AA_FIELDS = [
    "Locus_ID",
    "Chromosome",
    "Start_1based",
    "End_1based",
    "Strand",
    "Primary_Subfamily",
    "ORF",
    "AA_Position",
    "L1RP_AA",
    "Locus_AA",
    "AA_Difference_vs_L1RP",
    "Difference_Type",
    "L1RP_Codon",
    "Locus_Codon",
    "Nucleotide_Differences_In_Codon",
    "Script04_ORF_Status",
    "Script03_Alignment_QC",
    "ORF_Specific_Identity_pct",
    "Interpretation",
]


INDEL_FIELDS = [
    "Locus_ID",
    "Primary_Subfamily",
    "ORF",
    "Indel_Type",
    "Length_bp",
    "Frame_Effect",
    "Reference_Location",
    "Approx_AA_Position",
    "Sequence",
    "Script04_ORF_Status",
    "Script03_Alignment_QC",
    "ORF_Specific_Identity_pct",
    "Interpretation",
]


# ============================================================
# STEP 1 — LOAD SCRIPT 03 + SCRIPT 04
# ============================================================

print("\n============================================================")
print("SCRIPT 05 V3 — ORF-SPECIFIC SEMI-GLOBAL REALIGNMENT")
print("============================================================")

print("\nRun mode:", MODE)

print("\nScript 03:")
print(STRUCTURE_RESULTS)

print("\nScript 04:")
print(INTEGRITY_RESULTS)

print("\nReference:")
print(REFERENCE_GB)

print("\nFASTA ZIP:")
print(FASTA_ZIP)


integrity_df = pd.read_csv(
    INTEGRITY_RESULTS
)


structure_df = pd.read_csv(
    STRUCTURE_RESULTS,
    usecols=[
        "Locus_ID",
        "CIGAR",
    ],
).rename(
    columns={
        "CIGAR":
            "Script03_Whole_Locus_CIGAR"
    }
)


merged_df = (
    integrity_df
    .merge(
        structure_df,
        on="Locus_ID",
        how="left",
        validate="one_to_one",
    )
)


if (
    merged_df[
        "Script03_Whole_Locus_CIGAR"
    ]
    .isna()
    .any()
):

    raise ValueError(
        "One or more Script 04 loci are missing "
        "their Script 03 whole-locus CIGAR."
    )


print(
    f"\nLoaded loci: "
    f"{len(merged_df):,}"
)


# ============================================================
# STEP 2 — SELECT TEST / FULL LOCI
# ============================================================

if MODE == "FULL":

    work_df = (
        merged_df.copy()
    )


else:

    selected_ids = []

    available_ids = set(
        merged_df[
            "Locus_ID"
        ]
    )


    for locus_id in (
        TEST_LOCUS_IDS
    ):

        if locus_id in available_ids:

            selected_ids.append(
                locus_id
            )


    for locus_id in (
        merged_df[
            "Locus_ID"
        ]
    ):

        if (
            len(selected_ids)
            >= TEST_N
        ):

            break

        if (
            locus_id
            not in selected_ids
        ):

            selected_ids.append(
                locus_id
            )


    work_df = (
        merged_df
        .set_index(
            "Locus_ID"
        )
        .loc[
            selected_ids
        ]
        .reset_index()
    )


print(
    f"Selected loci: "
    f"{len(work_df):,}"
)


# ============================================================
# STEP 3 — LOAD L1-RP ORFS
# ============================================================

reference_record = SeqIO.read(
    REFERENCE_GB,
    "genbank",
)


reference_orfs = {}


for feature in (
    reference_record.features
):

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


    start0 = int(
        feature.location.start
    )

    end0 = int(
        feature.location.end
    )

    cds = clean_dna(
        feature.extract(
            reference_record.seq
        )
    )

    protein = (
        feature.qualifiers
        .get(
            "translation",
            [""],
        )[0]
    )


    reference_orfs[
        gene
    ] = {
        "start0":
            start0,

        "end0":
            end0,

        "start_1based":
            start0 + 1,

        "end_1based":
            end0,

        "cds":
            cds,

        "protein":
            protein,

        "protein_length":
            len(protein),
    }


for orf_name in [
    "ORF1",
    "ORF2",
]:

    info = (
        reference_orfs[
            orf_name
        ]
    )

    translated = str(
        Seq(
            info["cds"]
        ).translate(
            to_stop=False
        )
    )


    if (
        translated
        != info["protein"] + "*"
    ):

        raise ValueError(
            f"Reference translation QC failed "
            f"for {orf_name}."
        )


    print(
        f"{orf_name}: "
        f"{info['start_1based']:,}-"
        f"{info['end_1based']:,} "
        f"({len(info['cds']):,} bp; "
        f"{info['protein_length']} aa "
        f"+ stop)"
    )


# ============================================================
# STEP 4 — BUILD BIOPYTHON SEMI-GLOBAL ALIGNER
# ============================================================

# Start with Biopython's BLASTN scoring preset.
#
# PairwiseAligner convention:
#   target = L1-RP ORF
#   query  = candidate genomic sequence
#
# A target gap with query bases is an insertion relative to the
# target/reference. We make ONLY left/right insertions free,
# which permits extra genomic query flanks without forcing them
# into the ORF alignment.
#
# The reference/L1-RP ORF itself remains fully penalized at both
# ends and is therefore retained end-to-end.

aligner = Align.PairwiseAligner(
    scoring="blastn"
)

aligner.mode = "global"

aligner.open_left_insertion_score = 0
aligner.extend_left_insertion_score = 0

aligner.open_right_insertion_score = 0
aligner.extend_right_insertion_score = 0


print("\nAlignment engine:")
print(
    "Biopython PairwiseAligner "
    "(BLASTN scoring; semi-global query overhangs)"
)

print(
    "Algorithm selected by Biopython:",
    aligner.algorithm,
)


# ============================================================
# STEP 5 — LOOKUPS + OUTPUT STREAMS
# ============================================================

work_lookup = (
    work_df
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


alignment_rows = []
priority_aa_rows = []
priority_indel_rows = []
error_rows = []

processed_loci = 0
completed_orfs = 0
skipped_partial_orfs = 0


aa_gz_handle = gzip.open(
    AA_DIFFERENCES_GZ,
    "wt",
    newline="",
    encoding="utf-8",
)

aa_writer = csv.DictWriter(
    aa_gz_handle,
    fieldnames=AA_FIELDS,
)

aa_writer.writeheader()


indel_gz_handle = gzip.open(
    INDEL_EVENTS_GZ,
    "wt",
    newline="",
    encoding="utf-8",
)

indel_writer = csv.DictWriter(
    indel_gz_handle,
    fieldnames=INDEL_FIELDS,
)

indel_writer.writeheader()


if (
    MODE == "TEST"
    and WRITE_TEST_PLAIN_CSV
):

    aa_plain_handle = open(
        AA_DIFFERENCES_CSV,
        "w",
        newline="",
        encoding="utf-8",
    )

    aa_plain_writer = csv.DictWriter(
        aa_plain_handle,
        fieldnames=AA_FIELDS,
    )

    aa_plain_writer.writeheader()


    indel_plain_handle = open(
        INDEL_EVENTS_CSV,
        "w",
        newline="",
        encoding="utf-8",
    )

    indel_plain_writer = csv.DictWriter(
        indel_plain_handle,
        fieldnames=INDEL_FIELDS,
    )

    indel_plain_writer.writeheader()


else:

    aa_plain_handle = None
    aa_plain_writer = None

    indel_plain_handle = None
    indel_plain_writer = None


start_time = time.time()


# ============================================================
# CORE LOCUS PROCESSOR
# ============================================================

def process_locus(
    locus_id,
    raw_sequence,
):
    global processed_loci
    global completed_orfs
    global skipped_partial_orfs


    row = (
        work_lookup[
            locus_id
        ]
    )


    sequence = clean_dna(
        raw_sequence
    )

    sequence_length = len(
        sequence
    )

    orientation = (
        row[
            "FASTA_to_L1RP_Orientation"
        ]
    )


    if orientation == "+":

        oriented_query = (
            sequence
        )


    elif orientation == "-":

        oriented_query = (
            reverse_complement(
                sequence
            )
        )


    else:

        raise ValueError(
            f"Unexpected orientation for "
            f"{locus_id}: {orientation}"
        )


    whole_ops = parse_cigar(
        row[
            "Script03_Whole_Locus_CIGAR"
        ]
    )

    whole_query_start0 = (
        oriented_query_start0(
            row,
            sequence_length,
        )
    )

    whole_ref_start0 = (
        int(
            row[
                "L1RP_Alignment_Start_1based"
            ]
        )
        - 1
    )


    for orf_name in [
        "ORF1",
        "ORF2",
    ]:

        info = (
            reference_orfs[
                orf_name
            ]
        )


        if (
            REQUIRE_SCRIPT04_FULL_BOUNDS
            and not as_bool(
                row[
                    f"{orf_name}"
                    "_Full_Bounds_Covered"
                ]
            )
        ):

            skipped_partial_orfs += 1


            alignment_rows.append(
                {
                    "Locus_ID":
                        locus_id,

                    "Chromosome":
                        row[
                            "Chromosome"
                        ],

                    "Start_1based":
                        row[
                            "Start_1based"
                        ],

                    "End_1based":
                        row[
                            "End_1based"
                        ],

                    "Strand":
                        row[
                            "Strand"
                        ],

                    "Primary_Subfamily":
                        row[
                            "Primary_Subfamily"
                        ],

                    "ORF":
                        orf_name,

                    "Script04_ORF_Status":
                        row[
                            f"{orf_name}"
                            "_Sequence_Status"
                        ],

                    "Script03_Alignment_QC":
                        row[
                            "Alignment_QC"
                        ],

                    "Alignment_Status":
                        (
                            "SKIPPED_PARTIAL_"
                            "SCRIPT04_ORF_BOUNDS"
                        ),
                }
            )

            continue


        rough_bounds = (
            locate_orf_in_query_from_script03(
                cigar_ops=
                    whole_ops,

                whole_query_start0=
                    whole_query_start0,

                whole_ref_start0=
                    whole_ref_start0,

                orf_ref_start0=
                    info["start0"],

                orf_ref_end0=
                    info["end0"],
            )
        )


        if rough_bounds is None:

            error_rows.append(
                {
                    "Locus_ID":
                        locus_id,

                    "ORF":
                        orf_name,

                    "Error":
                        (
                            "Could not locate approximate "
                            "ORF query interval from "
                            "Script 03."
                        ),
                }
            )

            alignment_rows.append(
                {
                    "Locus_ID":
                        locus_id,

                    "Primary_Subfamily":
                        row[
                            "Primary_Subfamily"
                        ],

                    "ORF":
                        orf_name,

                    "Alignment_Status":
                        (
                            "ERROR_ORF_REGION_"
                            "NOT_LOCATED"
                        ),
                }
            )

            continue


        rough_start0, rough_end0 = (
            rough_bounds
        )


        candidate_start0 = max(
            0,
            rough_start0
            - ORF_QUERY_FLANK_BP,
        )

        candidate_end0 = min(
            sequence_length,
            rough_end0
            + 1
            + ORF_QUERY_FLANK_BP,
        )


        candidate = (
            oriented_query[
                candidate_start0:
                candidate_end0
            ]
        )


        try:

            alignments = (
                aligner.align(
                    info["cds"],
                    candidate,
                )
            )

            alignment = (
                alignments[0]
            )


            cigar_info = (
                alignment_to_core_ops(
                    alignment=
                        alignment,

                    reference_cds=
                        info["cds"],

                    candidate_sequence=
                        candidate,
                )
            )


            if (
                cigar_info[
                    "reference_start0"
                ]
                != 0
                or
                cigar_info[
                    "reference_end0_exclusive"
                ]
                != len(
                    info["cds"]
                )
            ):

                raise ValueError(
                    "Semi-global alignment did not "
                    "retain the complete L1-RP ORF."
                )


            maps = build_alignment_maps(
                reference_cds=
                    info["cds"],

                candidate_sequence=
                    candidate,

                core_ops=
                    cigar_info[
                        "core_ops"
                    ],

                query_core_start0=
                    cigar_info[
                        "query_core_start0"
                    ],
            )


            # -----------------------------------------------
            # Convert the core alignment from candidate coords
            # to the full oriented genomic LINE-1.
            # -----------------------------------------------

            full_oriented_core_start0 = (
                candidate_start0
                + cigar_info[
                    "query_core_start0"
                ]
            )

            full_oriented_core_end0_exclusive = (
                candidate_start0
                + cigar_info[
                    "query_core_end0_exclusive"
                ]
            )


            (
                original_core_start_1based,
                original_core_end_1based,
            ) = oriented_to_original_bounds(
                oriented_start0=
                    full_oriented_core_start0,

                oriented_end0_exclusive=
                    full_oriented_core_end0_exclusive,

                sequence_length=
                    sequence_length,

                orientation=
                    orientation,
            )


            genomic_core_start = (
                int(
                    row[
                        "Start_1based"
                    ]
                )
                + original_core_start_1based
                - 1
            )

            genomic_core_end = (
                int(
                    row[
                        "Start_1based"
                    ]
                )
                + original_core_end_1based
                - 1
            )


            genomic_core_low = min(
                genomic_core_start,
                genomic_core_end,
            )

            genomic_core_high = max(
                genomic_core_start,
                genomic_core_end,
            )


            # -----------------------------------------------
            # Oriented whole-query hard clipping.
            #
            # This is useful if the user BLASTs the full locus
            # against ONLY this WT ORF. It is not expected to
            # equal a BLAST HSP against the entire 6-kb WT L1.
            # -----------------------------------------------

            oriented_5prime_hardclip = (
                full_oriented_core_start0
            )

            oriented_3prime_hardclip = (
                sequence_length
                - full_oriented_core_end0_exclusive
            )


            sam_with_hardclips = ""

            if oriented_5prime_hardclip > 0:

                sam_with_hardclips += (
                    f"{oriented_5prime_hardclip}H"
                )


            sam_with_hardclips += (
                cigar_info[
                    "samstyle_cigar"
                ]
            )


            if oriented_3prime_hardclip > 0:

                sam_with_hardclips += (
                    f"{oriented_3prime_hardclip}H"
                )


            aa_analysis = (
                analyze_aa_differences(
                    locus_id=
                        locus_id,

                    row=
                        row,

                    orf_name=
                        orf_name,

                    reference_cds=
                        info["cds"],

                    reference_protein_length=
                        info[
                            "protein_length"
                        ],

                    maps=
                        maps,

                    identity_pct=
                        maps[
                            "identity_pct"
                        ],
                )
            )


            for aa_row in (
                aa_analysis[
                    "aa_rows"
                ]
            ):

                aa_writer.writerow(
                    aa_row
                )

                if aa_plain_writer:

                    aa_plain_writer.writerow(
                        aa_row
                    )


                if (
                    row[
                        f"{orf_name}"
                        "_Sequence_Status"
                    ]
                    ==
                    "FULL_SPAN_NO_FRAME_STOP_DISRUPTION"
                ):

                    priority_aa_rows.append(
                        aa_row
                    )


            indel_rows = (
                build_indel_rows(
                    locus_id=
                        locus_id,

                    row=
                        row,

                    orf_name=
                        orf_name,

                    reference_cds=
                        info["cds"],

                    maps=
                        maps,

                    identity_pct=
                        maps[
                            "identity_pct"
                        ],
                )
            )


            for indel_row in (
                indel_rows
            ):

                indel_writer.writerow(
                    indel_row
                )

                if indel_plain_writer:

                    indel_plain_writer.writerow(
                        indel_row
                    )


                if (
                    row[
                        f"{orf_name}"
                        "_Sequence_Status"
                    ]
                    ==
                    "FULL_SPAN_NO_FRAME_STOP_DISRUPTION"
                ):

                    priority_indel_rows.append(
                        indel_row
                    )


            alignment_rows.append(
                {
                    "Locus_ID":
                        locus_id,

                    "Chromosome":
                        row[
                            "Chromosome"
                        ],

                    "Start_1based":
                        row[
                            "Start_1based"
                        ],

                    "End_1based":
                        row[
                            "End_1based"
                        ],

                    "Strand":
                        row[
                            "Strand"
                        ],

                    "Primary_Subfamily":
                        row[
                            "Primary_Subfamily"
                        ],

                    "ORF":
                        orf_name,

                    "Script04_ORF_Status":
                        row[
                            f"{orf_name}"
                            "_Sequence_Status"
                        ],

                    "Script03_Alignment_QC":
                        row[
                            "Alignment_QC"
                        ],

                    "Alignment_Status":
                        "ALIGNMENT_COMPLETED",

                    "Alignment_Engine":
                        (
                            "Biopython_PairwiseAligner_"
                            "BLASTN_scoring_semiglobal"
                        ),

                    "ORF_Query_Flank_bp":
                        ORF_QUERY_FLANK_BP,

                    "Rough_ORF_Query_Start_Oriented_1based":
                        rough_start0 + 1,

                    "Rough_ORF_Query_End_Oriented_1based":
                        rough_end0 + 1,

                    "Candidate_Window_Start_Oriented_1based":
                        candidate_start0 + 1,

                    "Candidate_Window_End_Oriented_1based":
                        candidate_end0,

                    "Candidate_Window_Length_bp":
                        len(
                            candidate
                        ),

                    "Free_Query_5prime_bp_In_Candidate":
                        cigar_info[
                            "left_free_query_bp"
                        ],

                    "Free_Query_3prime_bp_In_Candidate":
                        cigar_info[
                            "right_free_query_bp"
                        ],

                    "Core_Query_Start_Original_FASTA_1based":
                        original_core_start_1based,

                    "Core_Query_End_Original_FASTA_1based":
                        original_core_end_1based,

                    "Core_Genomic_Start_1based":
                        genomic_core_low,

                    "Core_Genomic_End_1based":
                        genomic_core_high,

                    "Oriented_Query_5prime_Hardclip_bp":
                        oriented_5prime_hardclip,

                    "Oriented_Query_3prime_Hardclip_bp":
                        oriented_3prime_hardclip,

                    "L1RP_ORF_Reference_Length_bp":
                        len(
                            info["cds"]
                        ),

                    "Reference_ORF_EndToEnd_Span":
                        True,

                    "ORF_Specific_Alignment_Score":
                        float(
                            alignment.score
                        ),

                    "ORF_Specific_Identity_pct":
                        maps[
                            "identity_pct"
                        ],

                    "Exact_Match_bp":
                        maps[
                            "exact_matches"
                        ],

                    "Mismatch_bp_vs_L1RP":
                        maps[
                            "mismatches"
                        ],

                    "Insertion_bp_vs_L1RP":
                        maps[
                            "insertion_bp"
                        ],

                    "Deletion_bp_vs_L1RP":
                        maps[
                            "deletion_bp"
                        ],

                    "Largest_Insertion_bp":
                        maps[
                            "largest_insertion_bp"
                        ],

                    "Largest_Deletion_bp":
                        maps[
                            "largest_deletion_bp"
                        ],

                    "ORF_Specific_Extended_CIGAR":
                        cigar_info[
                            "extended_cigar"
                        ],

                    "ORF_Specific_SAMstyle_CIGAR":
                        cigar_info[
                            "samstyle_cigar"
                        ],

                    "ORF_Specific_SAMstyle_CIGAR_with_Oriented_Hardclips":
                        sam_with_hardclips,

                    "Synonymous_Codon_Differences":
                        aa_analysis[
                            "Synonymous_Codon_Differences"
                        ],

                    "Nonsynonymous_AA_Differences":
                        aa_analysis[
                            "Nonsynonymous_AA_Differences"
                        ],

                    "Missense_AA_Differences":
                        aa_analysis[
                            "Missense_AA_Differences"
                        ],

                    "Nonsense_AA_Differences":
                        aa_analysis[
                            "Nonsense_AA_Differences"
                        ],

                    "Start_Lost_Differences":
                        aa_analysis[
                            "Start_Lost_Differences"
                        ],

                    "First_Frameshift_AA":
                        aa_analysis[
                            "First_Frameshift_AA"
                        ],

                    "First_Nonsense_AA":
                        aa_analysis[
                            "First_Nonsense_AA"
                        ],

                    "AA_Callable_Through_AA":
                        aa_analysis[
                            "AA_Callable_Through_AA"
                        ],

                    "AA_Call_Stop_Reason":
                        aa_analysis[
                            "AA_Call_Stop_Reason"
                        ],

                    "Terminal_Stop_Status":
                        aa_analysis[
                            "Terminal_Stop_Status"
                        ],

                    "Uncallable_Indel_Affected_Codons":
                        aa_analysis[
                            "Uncallable_Indel_Affected_Codons"
                        ],

                    "Uncallable_Ambiguous_Codons":
                        aa_analysis[
                            "Uncallable_Ambiguous_Codons"
                        ],

                    "Difference_Interpretation":
                        (
                            "DIFFERENCES_VS_L1RP_"
                            "NOT_AUTOMATICALLY_"
                            "LOCUS_SPECIFIC_MUTATIONS"
                        ),
                }
            )


            completed_orfs += 1


        except Exception as exc:

            error_rows.append(
                {
                    "Locus_ID":
                        locus_id,

                    "ORF":
                        orf_name,

                    "Error":
                        str(
                            exc
                        ),
                }
            )


            alignment_rows.append(
                {
                    "Locus_ID":
                        locus_id,

                    "Primary_Subfamily":
                        row[
                            "Primary_Subfamily"
                        ],

                    "ORF":
                        orf_name,

                    "Script04_ORF_Status":
                        row[
                            f"{orf_name}"
                            "_Sequence_Status"
                        ],

                    "Script03_Alignment_QC":
                        row[
                            "Alignment_QC"
                        ],

                    "Alignment_Status":
                        "ERROR",

                    "Error_Message":
                        str(
                            exc
                        ),
                }
            )


    processed_loci += 1


    if (
        MODE == "TEST"
        or processed_loci % 100 == 0
    ):

        elapsed = (
            time.time()
            - start_time
        )

        print(
            f"  Processed "
            f"{processed_loci:,} / "
            f"{len(target_ids):,} loci "
            f"({elapsed/60:.2f} min)"
        )


# ============================================================
# STEP 6 — STREAM FASTA ZIP
# ============================================================

with zipfile.ZipFile(
    FASTA_ZIP,
    "r",
) as zf:

    fasta_members = [
        name
        for name
        in zf.namelist()
        if name.lower().endswith(
            (
                ".fa",
                ".fasta",
                ".fna",
                ".fas",
            )
        )
    ]


    if len(
        fasta_members
    ) != 1:

        raise ValueError(
            "Expected exactly one FASTA file "
            "inside the ZIP."
        )


    fasta_member = (
        fasta_members[0]
    )

    print(
        "\nFASTA inside ZIP:",
        fasta_member,
    )


    with zf.open(
        fasta_member,
        "r",
    ) as raw_handle:

        handle = (
            io.TextIOWrapper(
                raw_handle,
                encoding="utf-8",
                errors="replace",
            )
        )


        current_id = None
        sequence_parts = []


        for raw_line in handle:

            line = (
                raw_line.strip()
            )

            if not line:
                continue


            if line.startswith(">"):

                if (
                    current_id is not None
                    and current_id
                    in target_ids
                ):

                    process_locus(
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

                if (
                    current_id
                    in target_ids
                ):

                    sequence_parts.append(
                        line
                    )


        if (
            current_id is not None
            and current_id
            in target_ids
        ):

            process_locus(
                current_id,
                "".join(
                    sequence_parts
                ),
            )


# ============================================================
# CLOSE DETAILED STREAMS
# ============================================================

aa_gz_handle.close()
indel_gz_handle.close()

if aa_plain_handle:
    aa_plain_handle.close()

if indel_plain_handle:
    indel_plain_handle.close()


# ============================================================
# VERIFY PROCESSING
# ============================================================

print(
    f"\nRequested loci: "
    f"{len(target_ids):,}"
)

print(
    f"Processed loci: "
    f"{processed_loci:,}"
)


if (
    processed_loci
    != len(target_ids)
):

    raise ValueError(
        "Not all selected loci were found "
        "in the FASTA."
    )


# ============================================================
# SAVE ALIGNMENT RESULTS
# ============================================================

alignment_df = pd.DataFrame(
    alignment_rows
)


alignment_df.to_csv(
    ALIGNMENT_RESULTS_GZ,
    index=False,
    compression="gzip",
)


if (
    MODE == "TEST"
    and WRITE_TEST_PLAIN_CSV
):

    alignment_df.to_csv(
        ALIGNMENT_RESULTS_CSV,
        index=False,
    )


# ============================================================
# SUMMARY WORKBOOK
# ============================================================

status_counts = (
    alignment_df[
        "Alignment_Status"
    ]
    .value_counts(
        dropna=False
    )
    .rename_axis(
        "Alignment_Status"
    )
    .reset_index(
        name="ORF_Count"
    )
)


completed_df = (
    alignment_df.loc[
        alignment_df[
            "Alignment_Status"
        ]
        == "ALIGNMENT_COMPLETED"
    ]
    .copy()
)


summary_rows = [
    [
        "Script version",
        "05_V3_Biopython_semiglobal",
    ],
    [
        "Run mode",
        MODE,
    ],
    [
        "Loci selected",
        len(target_ids),
    ],
    [
        "Loci processed",
        processed_loci,
    ],
    [
        "ORF-specific alignments completed",
        completed_orfs,
    ],
    [
        "Partial Script 04 ORFs skipped",
        skipped_partial_orfs,
    ],
    [
        "Alignment errors",
        len(error_rows),
    ],
    [
        "ORF query flank bp",
        ORF_QUERY_FLANK_BP,
    ],
    [
        "Alignment engine",
        (
            "Biopython PairwiseAligner; "
            "BLASTN scoring; global DP "
            "with free query end overhangs"
        ),
    ],
    [
        "Difference terminology",
        (
            "Differences vs L1-RP; "
            "not automatically locus-specific mutations"
        ),
    ],
]


if not completed_df.empty:

    summary_rows.extend(
        [
            [
                "Median ORF-specific identity (%)",
                round(
                    completed_df[
                        "ORF_Specific_Identity_pct"
                    ]
                    .median(),
                    3,
                ),
            ],
            [
                "Completed ORFs with full L1-RP reference span",
                int(
                    completed_df[
                        "Reference_ORF_EndToEnd_Span"
                    ]
                    .fillna(False)
                    .sum()
                ),
            ],
        ]
    )


summary_df = pd.DataFrame(
    summary_rows,
    columns=[
        "Metric",
        "Value",
    ],
)


if not completed_df.empty:

    orf_summary = (
        completed_df
        .groupby(
            "ORF"
        )
        .agg(
            ORF_Count=(
                "Locus_ID",
                "count",
            ),

            Median_Identity_pct=(
                "ORF_Specific_Identity_pct",
                "median",
            ),

            Total_Nonsynonymous_AA_Differences=(
                "Nonsynonymous_AA_Differences",
                "sum",
            ),

            Total_Synonymous_Codon_Differences=(
                "Synonymous_Codon_Differences",
                "sum",
            ),

            ORFs_With_Frameshift=(
                "First_Frameshift_AA",
                lambda x:
                    x.notna().sum(),
            ),

            ORFs_With_Nonsense=(
                "First_Nonsense_AA",
                lambda x:
                    x.notna().sum(),
            ),
        )
        .reset_index()
    )


    subfamily_summary = (
        completed_df
        .groupby(
            [
                "Primary_Subfamily",
                "ORF",
            ],
            dropna=False,
        )
        .agg(
            ORF_Count=(
                "Locus_ID",
                "count",
            ),

            Median_Identity_pct=(
                "ORF_Specific_Identity_pct",
                "median",
            ),

            Total_Nonsynonymous_AA_Differences=(
                "Nonsynonymous_AA_Differences",
                "sum",
            ),
        )
        .reset_index()
    )


else:

    orf_summary = (
        pd.DataFrame()
    )

    subfamily_summary = (
        pd.DataFrame()
    )


priority_aa_df = pd.DataFrame(
    priority_aa_rows,
    columns=AA_FIELDS,
)


priority_indel_df = pd.DataFrame(
    priority_indel_rows,
    columns=INDEL_FIELDS,
)


error_df = pd.DataFrame(
    error_rows
)


with pd.ExcelWriter(
    SUMMARY_XLSX,
    engine="openpyxl",
) as writer:

    summary_df.to_excel(
        writer,
        sheet_name="Summary",
        index=False,
    )

    status_counts.to_excel(
        writer,
        sheet_name="Alignment_Status",
        index=False,
    )


    if not orf_summary.empty:

        orf_summary.to_excel(
            writer,
            sheet_name="ORF_Summary",
            index=False,
        )


    if not subfamily_summary.empty:

        subfamily_summary.to_excel(
            writer,
            sheet_name="Subfamily_Summary",
            index=False,
        )


    if MODE == "TEST":

        alignment_df.to_excel(
            writer,
            sheet_name="TEST_Details",
            index=False,
        )


    else:

        compact_columns = [
            column
            for column in [
                "Locus_ID",
                "Primary_Subfamily",
                "ORF",
                "Script04_ORF_Status",
                "Script03_Alignment_QC",
                "Alignment_Status",
                "ORF_Specific_Identity_pct",
                "Mismatch_bp_vs_L1RP",
                "Insertion_bp_vs_L1RP",
                "Deletion_bp_vs_L1RP",
                "First_Frameshift_AA",
                "First_Nonsense_AA",
                "Nonsynonymous_AA_Differences",
                "Synonymous_Codon_Differences",
                "Terminal_Stop_Status",
            ]
            if column
            in alignment_df.columns
        ]


        alignment_df[
            compact_columns
        ].to_excel(
            writer,
            sheet_name="Alignment_Compact",
            index=False,
        )


    if not priority_aa_df.empty:

        priority_aa_df.head(
            500_000
        ).to_excel(
            writer,
            sheet_name="Priority_AA_Differences",
            index=False,
        )


    if not priority_indel_df.empty:

        priority_indel_df.head(
            500_000
        ).to_excel(
            writer,
            sheet_name="Priority_ORF_Indels",
            index=False,
        )


    if not error_df.empty:

        error_df.to_excel(
            writer,
            sheet_name="Errors",
            index=False,
        )


# ============================================================
# FINAL REPORT
# ============================================================

elapsed = (
    time.time()
    - start_time
)


print("\n============================================================")
print("SCRIPT 05 V3 COMPLETE")
print("============================================================")

print(
    f"\nElapsed time: "
    f"{elapsed/60:.2f} minutes"
)

print(
    f"Loci processed: "
    f"{processed_loci:,}"
)

print(
    f"ORF alignments completed: "
    f"{completed_orfs:,}"
)

print(
    f"Partial ORFs skipped: "
    f"{skipped_partial_orfs:,}"
)

print(
    f"Errors: "
    f"{len(error_rows):,}"
)


print("\nOutputs:")

print(
    "  ",
    ALIGNMENT_RESULTS_GZ,
)

print(
    "  ",
    AA_DIFFERENCES_GZ,
)

print(
    "  ",
    INDEL_EVENTS_GZ,
)

print(
    "  ",
    SUMMARY_XLSX,
)


if MODE == "TEST":

    print(
        "\nSanity-check locus:"
    )

    print(
        "  hg38.RE.chr1.100199603-100206089.+"
    )

    print(
        "\nFor BLAST comparison, compare ORF1 to a "
        "BLAST alignment using the full genomic locus "
        "as query and the WT L1-RP ORF1 CDS as subject."
    )

    print(
        "\nIf TEST QC passes, change:"
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


print(
    "\nINTERPRETATION REMINDER:"
)

print(
    "AA calls are differences relative to "
    "L1-RP / AF148856.1."
)

print(
    "They are not automatically locus-specific mutations."
)
