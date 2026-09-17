# ============================================================
# LINE-1 FUNCTIONAL ANNOTATION PROJECT
#
# SCRIPT 03:
# Characterize LINE-1 locus structure relative to L1-RP
#
# WT/reference:
#   L1-RP
#   GenBank accession: AF148856.1
#
# WHAT THIS SCRIPT DOES:
#   1. Loads the prepared locus manifest from Script 02.
#   2. Selects candidate loci for detailed sequence alignment.
#   3. Downloads/caches the AF148856.1 GenBank record.
#   4. Reads ORF1 and ORF2 coordinates directly from GenBank.
#   5. Aligns each candidate locus to L1-RP in BOTH orientations.
#   6. Reports reference coverage of:
#        - 5' region before ORF1
#        - ORF1
#        - inter-ORF region
#        - ORF2
#        - 3' region after ORF2
#   7. Also reports coverage of selected ORF1/ORF2 protein domains.
#
# IMPORTANT:
#   - This is STRUCTURAL annotation only.
#   - It does NOT call a locus "active" or "inactive".
#   - It does NOT yet call point mutations.
#   - It does NOT treat every difference from L1-RP as a mutation.
#   - Shorter loci remain in the master Script 02 manifest and are
#     NOT deleted. This script performs detailed alignment on a
#     configurable subset to keep the analysis computationally
#     manageable.
# ============================================================


# ============================================================
# IMPORTS
# ============================================================

import csv
import io
import re
import sys
import time
import zipfile
from pathlib import Path

import pandas as pd
from Bio import Entrez, SeqIO
from Bio.Seq import Seq

try:
    import parasail
except ImportError:
    raise ImportError(
        "\nThe 'parasail' alignment package is not installed.\n\n"
        "Open Command Prompt and run:\n\n"
        "    pip install parasail\n\n"
        "Then run this script again.\n"
    )


# ============================================================
#                 EDIT THESE SETTINGS
# ============================================================

# NCBI requests an email address for Entrez access.
# Put your email inside the quotation marks.
Entrez.email = "manaal.dalwadi@cuanschutz.edu"


# ------------------------------------------------------------
# FIRST RUN:
# Keep RUN_MODE = "TEST".
#
# TEST mode aligns only a small number of loci so we can make
# sure Parasail and the reference download work correctly.
#
# After the TEST run works, change:
#
#     RUN_MODE = "FULL"
#
# ------------------------------------------------------------

RUN_MODE = "TEST"

# Number of candidate loci to align in TEST mode.
TEST_N = 25


# ------------------------------------------------------------
# DETAILED ALIGNMENT LENGTH THRESHOLD
# ------------------------------------------------------------
#
# Script 02 retains all ~983k loci.
#
# For this first structural-alignment pass, detailed alignment
# is performed on loci >= 5,000 bp by default.
#
# This is a TRIAGE threshold, not a functional threshold.
# Shorter loci are not discarded from the project.
#
# Later, we can:
#   - lower this threshold,
#   - analyze an expressed-locus subset,
#   - or analyze selected shorter loci separately.
#
MIN_LOCUS_LENGTH_FOR_ALIGNMENT = 5000


# ------------------------------------------------------------
# INPUT FILE LOCATIONS
# ------------------------------------------------------------

PREPARED_FOLDER_NAME = "02_prepared_LINE1_data"

PREPARED_MANIFEST_NAME = "hg38_LINE1_locus_manifest.csv.gz"

FASTA_ZIP_NAME = "hg38_RE_L1_seqs.zip"


# ------------------------------------------------------------
# OUTPUT FOLDER
# ------------------------------------------------------------

OUTPUT_FOLDER_NAME = "03_LINE1_structure"


# ------------------------------------------------------------
# ALIGNMENT SCORING
# ------------------------------------------------------------
#
# We use Parasail Smith-Waterman local DNA alignment.
# nuc44 is a nucleotide substitution matrix bundled with Parasail.
#
GAP_OPEN = 10
GAP_EXTEND = 2


# ============================================================
#             DO NOT NEED TO EDIT BELOW HERE
# ============================================================

SCRIPT_FOLDER = Path(__file__).resolve().parent

PREPARED_MANIFEST = (
    SCRIPT_FOLDER
    / PREPARED_FOLDER_NAME
    / PREPARED_MANIFEST_NAME
)

FASTA_ZIP = SCRIPT_FOLDER / FASTA_ZIP_NAME

OUTPUT_FOLDER = SCRIPT_FOLDER / OUTPUT_FOLDER_NAME
OUTPUT_FOLDER.mkdir(parents=True, exist_ok=True)

REFERENCE_GB = OUTPUT_FOLDER / "AF148856.1.gb"
REFERENCE_FASTA = OUTPUT_FOLDER / "AF148856.1.fasta"


if RUN_MODE.upper() == "TEST":
    OUTPUT_PREFIX = "TEST"
else:
    OUTPUT_PREFIX = "FULL"


ALIGNMENT_CSV = (
    OUTPUT_FOLDER
    / f"{OUTPUT_PREFIX}_LINE1_structure_alignment_results.csv"
)

ALIGNMENT_GZ = (
    OUTPUT_FOLDER
    / f"{OUTPUT_PREFIX}_LINE1_structure_alignment_results.csv.gz"
)

SUMMARY_XLSX = (
    OUTPUT_FOLDER
    / f"{OUTPUT_PREFIX}_LINE1_structure_summary.xlsx"
)


# ============================================================
# FILE CHECKS
# ============================================================

if not PREPARED_MANIFEST.exists():
    raise FileNotFoundError(
        f"\nCould not find the Script 02 manifest:\n"
        f"{PREPARED_MANIFEST}\n"
    )

if not FASTA_ZIP.exists():
    raise FileNotFoundError(
        f"\nCould not find the FASTA ZIP:\n"
        f"{FASTA_ZIP}\n"
    )

if Entrez.email == "YOUR_EMAIL_HERE@email.com":
    raise ValueError(
        "\nPlease edit Entrez.email near the top of this script "
        "before running it."
    )


# ============================================================
# HELPER FUNCTIONS
# ============================================================

def reverse_complement(sequence):
    return str(Seq(sequence).reverse_complement())


def clean_dna(sequence):
    """
    Convert sequence to uppercase.
    Non-IUPAC whitespace is removed by FASTA parsing.
    Parasail nuc44 supports common nucleotide ambiguity codes.
    """
    return sequence.upper().replace("U", "T")


def get_cigar_string(result):
    """
    Return a decoded SAM-style CIGAR string from a Parasail
    traceback-capable result.
    """
    cigar = result.cigar
    decoded = cigar.decode

    if isinstance(decoded, bytes):
        decoded = decoded.decode("ascii")

    return str(decoded)


def parse_cigar(cigar_string):
    """
    Convert e.g. '100M2I50M3D20M' into:
        [(100, 'M'), (2, 'I'), ...]
    """
    return [
        (int(length), op)
        for length, op
        in re.findall(r"(\d+)([MIDNSHP=X])", cigar_string)
    ]


def merge_intervals(intervals):
    """
    Merge 0-based half-open intervals.
    """
    if not intervals:
        return []

    intervals = sorted(intervals)

    merged = [list(intervals[0])]

    for start, end in intervals[1:]:

        last_start, last_end = merged[-1]

        if start <= last_end:
            merged[-1][1] = max(last_end, end)
        else:
            merged.append([start, end])

    return [
        (start, end)
        for start, end in merged
    ]


def interval_overlap_bp(intervals, feature_start_1based, feature_end_1based):
    """
    Count how many reference bases in 'intervals' overlap a
    1-based inclusive reference feature.

    intervals are 0-based half-open.
    """
    feature_start0 = feature_start_1based - 1
    feature_end0 = feature_end_1based

    total = 0

    for start0, end0 in intervals:
        overlap_start = max(start0, feature_start0)
        overlap_end = min(end0, feature_end0)

        if overlap_end > overlap_start:
            total += overlap_end - overlap_start

    return total


def feature_length(start_1based, end_1based):
    if end_1based < start_1based:
        return 0

    return end_1based - start_1based + 1


def coverage_percent(covered_bp, start_1based, end_1based):
    length = feature_length(start_1based, end_1based)

    if length <= 0:
        return None

    return round(
        100.0 * covered_bp / length,
        3
    )


def aa_to_reference_nt(
    orf_start_1based,
    aa_start,
    aa_end
):
    """
    Convert amino-acid coordinates within an ORF to nucleotide
    coordinates on the AF148856 reference.

    The returned coordinates are 1-based inclusive.
    """
    nt_start = (
        orf_start_1based
        + (aa_start - 1) * 3
    )

    nt_end = (
        orf_start_1based
        + aa_end * 3
        - 1
    )

    return nt_start, nt_end


def analyze_alignment(
    query_sequence,
    reference_sequence,
    result,
    orientation,
    features
):
    """
    Convert one Parasail local alignment into structural metrics.
    """

    cigar_string = get_cigar_string(result)
    cigar_ops = parse_cigar(cigar_string)

    # Parasail coordinates are 0-based.
    qpos = int(result.cigar.beg_query)
    rpos = int(result.cigar.beg_ref)

    q_start0 = qpos
    r_start0 = rpos

    covered_reference_intervals = []

    exact_matches = 0
    aligned_pair_columns = 0

    query_consumed = 0
    reference_consumed = 0

    insertions_bp = 0
    deletions_bp = 0

    for length, op in cigar_ops:

        # Alignment column consuming query + reference.
        if op in ("M", "=", "X"):

            q_segment = query_sequence[
                qpos:qpos + length
            ]

            r_segment = reference_sequence[
                rpos:rpos + length
            ]

            if op == "=":
                exact_matches += length

            elif op == "X":
                pass

            else:
                exact_matches += sum(
                    a == b
                    for a, b
                    in zip(q_segment, r_segment)
                )

            aligned_pair_columns += length

            covered_reference_intervals.append(
                (rpos, rpos + length)
            )

            qpos += length
            rpos += length

            query_consumed += length
            reference_consumed += length

        # Insertion relative to reference.
        elif op == "I":
            qpos += length
            query_consumed += length
            insertions_bp += length

        # Deletion relative to reference.
        elif op in ("D", "N"):
            rpos += length
            reference_consumed += length
            deletions_bp += length

        # Soft clip consumes query.
        elif op == "S":
            qpos += length

        # H and P consume neither sequence.
        elif op in ("H", "P"):
            pass


    covered_reference_intervals = merge_intervals(
        covered_reference_intervals
    )

    reference_aligned_bp = sum(
        end - start
        for start, end
        in covered_reference_intervals
    )

    identity_pct = None

    if aligned_pair_columns > 0:
        identity_pct = round(
            100.0
            * exact_matches
            / aligned_pair_columns,
            3
        )


    # Local-alignment bounds.
    ref_start_1based = (
        int(result.cigar.beg_ref) + 1
    )

    ref_end_1based = (
        int(result.end_ref) + 1
    )

    q_begin = int(result.cigar.beg_query)
    q_end = int(result.end_query)

    query_length = len(query_sequence)


    # Convert coordinates back to the ORIGINAL FASTA orientation.
    if orientation == "+":

        original_query_start_1based = q_begin + 1
        original_query_end_1based = q_end + 1

    else:

        original_query_start_1based = (
            query_length - q_end
        )

        original_query_end_1based = (
            query_length - q_begin
        )


    query_alignment_span_bp = (
        q_end - q_begin + 1
    )

    query_alignment_span_pct = round(
        100.0
        * query_alignment_span_bp
        / query_length,
        3
    )


    # Reference bases apparently absent before/after the aligned span.
    # These are APPROXIMATE because local alignment can stop in
    # highly diverged sequence even when sequence is physically present.
    approx_5prime_missing_ref_bp = max(
        0,
        ref_start_1based - 1
    )

    approx_3prime_missing_ref_bp = max(
        0,
        len(reference_sequence) - ref_end_1based
    )


    output = {
        "Alignment_Score": int(result.score),
        "FASTA_to_L1RP_Orientation": orientation,

        "Query_Alignment_Start_1based":
            original_query_start_1based,

        "Query_Alignment_End_1based":
            original_query_end_1based,

        "Query_Alignment_Span_bp":
            query_alignment_span_bp,

        "Query_Alignment_Span_pct":
            query_alignment_span_pct,

        "L1RP_Alignment_Start_1based":
            ref_start_1based,

        "L1RP_Alignment_End_1based":
            ref_end_1based,

        "Approx_5prime_Missing_L1RP_bp":
            approx_5prime_missing_ref_bp,

        "Approx_3prime_Missing_L1RP_bp":
            approx_3prime_missing_ref_bp,

        "Reference_Aligned_bp":
            reference_aligned_bp,

        "Alignment_Identity_pct":
            identity_pct,

        "Insertion_bp_vs_L1RP":
            insertions_bp,

        "Deletion_bp_vs_L1RP":
            deletions_bp,

        "CIGAR":
            cigar_string,
    }


    # Feature/domain coverage.
    for feature_name, feature_info in features.items():

        start = feature_info["start"]
        end = feature_info["end"]

        covered = interval_overlap_bp(
            covered_reference_intervals,
            start,
            end
        )

        output[
            f"{feature_name}_Coverage_bp"
        ] = covered

        output[
            f"{feature_name}_Coverage_pct"
        ] = coverage_percent(
            covered,
            start,
            end
        )


    return output


def align_one_orientation(
    query,
    reference,
    orientation,
    features
):
    """
    Smith-Waterman local alignment using Parasail.

    A traceback-capable vectorized implementation is used so we
    can recover the CIGAR and measure feature coverage.
    """

    result = parasail.sw_trace_striped_16(
        query,
        reference,
        GAP_OPEN,
        GAP_EXTEND,
        parasail.nuc44
    )

    metrics = analyze_alignment(
        query,
        reference,
        result,
        orientation,
        features
    )

    return metrics


def align_both_orientations(
    sequence,
    reference,
    features
):
    """
    Align the FASTA sequence in both orientations and retain the
    higher-scoring alignment.
    """

    forward = clean_dna(sequence)

    reverse = reverse_complement(
        forward
    )

    forward_metrics = align_one_orientation(
        forward,
        reference,
        "+",
        features
    )

    reverse_metrics = align_one_orientation(
        reverse,
        reference,
        "-",
        features
    )

    if (
        reverse_metrics["Alignment_Score"]
        > forward_metrics["Alignment_Score"]
    ):
        return reverse_metrics

    return forward_metrics


# ============================================================
# STEP 1 — LOAD PREPARED MANIFEST
# ============================================================

print("\n============================================================")
print("SCRIPT 03 — CHARACTERIZE LINE-1 STRUCTURE")
print("============================================================")

print("\nRun mode:", RUN_MODE.upper())
print(
    "Minimum locus length for detailed alignment:",
    f"{MIN_LOCUS_LENGTH_FOR_ALIGNMENT:,} bp"
)

print("\nLoading Script 02 manifest...")


manifest_columns = [
    "Locus_ID",
    "Chromosome",
    "Start_1based",
    "End_1based",
    "Strand",
    "Coordinate_Length",
    "Primary_Subfamily",
    "All_Subfamilies",
    "FASTA_QC",
    "Contained_Within_Other_L1",
    "Contains_Other_L1",
]

manifest = pd.read_csv(
    PREPARED_MANIFEST,
    compression="gzip",
    usecols=manifest_columns
)

print(
    f"Total loci in manifest: "
    f"{len(manifest):,}"
)


# Detailed candidate set.
candidate_df = manifest.loc[
    (manifest["FASTA_QC"] == "PASS")
    &
    (
        manifest["Coordinate_Length"]
        >= MIN_LOCUS_LENGTH_FOR_ALIGNMENT
    )
].copy()


print(
    f"Candidate loci meeting current alignment criteria: "
    f"{len(candidate_df):,}"
)


if RUN_MODE.upper() == "TEST":

    candidate_df = candidate_df.head(
        TEST_N
    ).copy()

    print(
        f"TEST mode: only the first "
        f"{len(candidate_df):,} candidates will be aligned."
    )


candidate_metadata = (
    candidate_df
    .set_index("Locus_ID")
    .to_dict(orient="index")
)

candidate_ids = set(
    candidate_metadata.keys()
)


# ============================================================
# STEP 2 — DOWNLOAD/CACHE AF148856.1 GENBANK
# ============================================================

print("\n============================================================")
print("STEP 2 — Loading AF148856.1 / L1-RP reference")
print("============================================================")


if not REFERENCE_GB.exists():

    print(
        "\nAF148856.1 GenBank record not cached."
    )

    print(
        "Downloading one record from NCBI..."
    )

    with Entrez.efetch(
        db="nuccore",
        id="AF148856.1",
        rettype="gb",
        retmode="text"
    ) as handle:

        reference_text = handle.read()

    REFERENCE_GB.write_text(
        reference_text,
        encoding="utf-8"
    )

    print(
        f"Saved reference cache:\n"
        f"{REFERENCE_GB}"
    )


reference_record = SeqIO.read(
    REFERENCE_GB,
    "genbank"
)

reference_sequence = (
    str(reference_record.seq)
    .upper()
)

reference_length = len(
    reference_sequence
)


SeqIO.write(
    reference_record,
    REFERENCE_FASTA,
    "fasta"
)


print(
    f"\nAF148856.1 length: "
    f"{reference_length:,} bp"
)


# ============================================================
# STEP 3 — PARSE ORF1 / ORF2 COORDINATES FROM GENBANK
# ============================================================

print("\n============================================================")
print("STEP 3 — Parsing L1-RP ORF coordinates")
print("============================================================")


orf_coordinates = {}


for feature in reference_record.features:

    if feature.type != "CDS":
        continue

    gene = (
        feature.qualifiers
        .get("gene", [""])[0]
        .upper()
    )

    if gene not in ("ORF1", "ORF2"):
        continue

    # BioPython feature start = 0-based.
    # feature end behaves as half-open.
    start_1based = int(
        feature.location.start
    ) + 1

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
        "protein_length_aa": len(translation),
    }


if (
    "ORF1" not in orf_coordinates
    or "ORF2" not in orf_coordinates
):
    raise ValueError(
        "\nCould not identify both ORF1 and ORF2 "
        "from the AF148856.1 GenBank CDS annotations."
    )


orf1_start = orf_coordinates[
    "ORF1"
]["start"]

orf1_end = orf_coordinates[
    "ORF1"
]["end"]

orf2_start = orf_coordinates[
    "ORF2"
]["start"]

orf2_end = orf_coordinates[
    "ORF2"
]["end"]


print(
    f"ORF1: {orf1_start:,}-{orf1_end:,} "
    f"({orf_coordinates['ORF1']['protein_length_aa']} aa)"
)

print(
    f"ORF2: {orf2_start:,}-{orf2_end:,} "
    f"({orf_coordinates['ORF2']['protein_length_aa']} aa)"
)


# ============================================================
# STEP 4 — DEFINE REFERENCE STRUCTURAL FEATURES
# ============================================================

print("\n============================================================")
print("STEP 4 — Defining reference architecture")
print("============================================================")


features = {
    "Pre_ORF1_5prime_Region": {
        "start": 1,
        "end": orf1_start - 1,
        "source": "Derived from AF148856.1 ORF1 CDS start",
    },

    "ORF1": {
        "start": orf1_start,
        "end": orf1_end,
        "source": "AF148856.1 GenBank CDS",
    },

    "Inter_ORF_Region": {
        "start": orf1_end + 1,
        "end": orf2_start - 1,
        "source": "Derived from AF148856.1 ORF1/ORF2 CDS coordinates",
    },

    "ORF2": {
        "start": orf2_start,
        "end": orf2_end,
        "source": "AF148856.1 GenBank CDS",
    },

    "Post_ORF2_3prime_Region": {
        "start": orf2_end + 1,
        "end": reference_length,
        "source": "Derived from AF148856.1 ORF2 CDS end",
    },
}


# ------------------------------------------------------------
# ORF1p domain coordinates
#
# Literature boundaries relative to L1-RP ORF1p:
#   CC   = aa 52-153
#   RRM  = aa 157-252
#   CTD  = aa 264-323
#
# We also define NTR = aa 1-51 as the region preceding CC.
# ------------------------------------------------------------

orf1_domains_aa = {
    "ORF1_NTR": (1, 51),
    "ORF1_CC": (52, 153),
    "ORF1_RRM": (157, 252),
    "ORF1_CTD": (264, 323),
}


for name, (aa_start, aa_end) in orf1_domains_aa.items():

    nt_start, nt_end = aa_to_reference_nt(
        orf1_start,
        aa_start,
        aa_end
    )

    features[name] = {
        "start": nt_start,
        "end": nt_end,
        "source":
            f"L1-RP ORF1p aa {aa_start}-{aa_end}",
    }


# ------------------------------------------------------------
# ORF2p domain coordinates
#
# Literature boundaries relative to L1-RP ORF2p:
#   EN    = aa 1-239
#   RT    = aa 453-883
#   CCHC  = aa 1096-1275
# ------------------------------------------------------------

orf2_domains_aa = {
    "ORF2_EN": (1, 239),
    "ORF2_RT": (453, 883),
    "ORF2_CCHC": (1096, 1275),
}


for name, (aa_start, aa_end) in orf2_domains_aa.items():

    nt_start, nt_end = aa_to_reference_nt(
        orf2_start,
        aa_start,
        aa_end
    )

    features[name] = {
        "start": nt_start,
        "end": nt_end,
        "source":
            f"L1-RP ORF2p aa {aa_start}-{aa_end}",
    }


for name, info in features.items():

    print(
        f"{name}: "
        f"{info['start']:,}-{info['end']:,}"
    )


# ============================================================
# STEP 5 — PREPARE OUTPUT COLUMNS
# ============================================================

metadata_columns = [
    "Locus_ID",
    "Chromosome",
    "Start_1based",
    "End_1based",
    "Strand",
    "Coordinate_Length",
    "Primary_Subfamily",
    "All_Subfamilies",
    "Contained_Within_Other_L1",
    "Contains_Other_L1",
]


alignment_columns = [
    "Alignment_Status",
    "Alignment_Error",

    "Alignment_Score",
    "FASTA_to_L1RP_Orientation",

    "Query_Alignment_Start_1based",
    "Query_Alignment_End_1based",
    "Query_Alignment_Span_bp",
    "Query_Alignment_Span_pct",

    "L1RP_Alignment_Start_1based",
    "L1RP_Alignment_End_1based",

    "Approx_5prime_Missing_L1RP_bp",
    "Approx_3prime_Missing_L1RP_bp",

    "Reference_Aligned_bp",
    "Alignment_Identity_pct",

    "Insertion_bp_vs_L1RP",
    "Deletion_bp_vs_L1RP",

    "CIGAR",
]


feature_columns = []

for feature_name in features:

    feature_columns.extend([
        f"{feature_name}_Coverage_bp",
        f"{feature_name}_Coverage_pct",
    ])


output_columns = (
    metadata_columns
    + alignment_columns
    + feature_columns
)


# ============================================================
# STEP 6 — STREAM FASTA AND ALIGN CANDIDATES
# ============================================================

print("\n============================================================")
print("STEP 6 — Aligning candidate loci to L1-RP")
print("============================================================")

print(
    "\nThe sequence FASTA is read directly from the ZIP."
)

print(
    "Each selected locus is aligned in both orientations."
)

print(
    "\nThis is the computationally expensive part of Script 03."
)


start_time = time.time()

processed = 0
errors = 0


with zipfile.ZipFile(
    FASTA_ZIP,
    "r"
) as zf:

    fasta_members = [
        name for name in zf.namelist()
        if name.lower().endswith(
            (".fa", ".fasta", ".fna", ".fas")
        )
    ]

    if len(fasta_members) != 1:
        raise ValueError(
            "\nExpected exactly one FASTA file in the ZIP.\n"
            f"Found: {fasta_members}"
        )

    fasta_member = fasta_members[0]

    print(
        f"\nFASTA inside ZIP: "
        f"{fasta_member}"
    )


    with open(
        ALIGNMENT_CSV,
        "w",
        newline="",
        encoding="utf-8"
    ) as output_handle:

        writer = csv.DictWriter(
            output_handle,
            fieldnames=output_columns
        )

        writer.writeheader()


        with zf.open(
            fasta_member,
            "r"
        ) as raw_handle:

            text_handle = io.TextIOWrapper(
                raw_handle,
                encoding="utf-8",
                errors="replace"
            )


            current_id = None
            sequence_parts = []


            def process_record(
                fasta_id,
                sequence
            ):
                global processed, errors

                if fasta_id not in candidate_ids:
                    return

                metadata = candidate_metadata[
                    fasta_id
                ]

                row = {
                    "Locus_ID": fasta_id,
                    **metadata,
                }

                # Locus_ID already supplied explicitly.
                row["Locus_ID"] = fasta_id


                try:

                    metrics = align_both_orientations(
                        sequence,
                        reference_sequence,
                        features
                    )

                    row.update(metrics)

                    row["Alignment_Status"] = (
                        "ALIGNMENT_COMPLETED"
                    )

                    row["Alignment_Error"] = ""


                except Exception as exc:

                    errors += 1

                    row["Alignment_Status"] = (
                        "ALIGNMENT_ERROR"
                    )

                    row["Alignment_Error"] = str(exc)


                writer.writerow(
                    {
                        column: row.get(
                            column,
                            ""
                        )
                        for column in output_columns
                    }
                )

                processed += 1


                if (
                    RUN_MODE.upper() == "TEST"
                    or processed % 100 == 0
                ):

                    elapsed = (
                        time.time()
                        - start_time
                    )

                    print(
                        f"  Processed "
                        f"{processed:,} / "
                        f"{len(candidate_ids):,} "
                        f"candidate loci "
                        f"({elapsed/60:.1f} min)"
                    )


            # FASTA parser.
            for raw_line in text_handle:

                line = raw_line.strip()

                if not line:
                    continue


                if line.startswith(">"):

                    if current_id is not None:

                        process_record(
                            current_id,
                            "".join(
                                sequence_parts
                            )
                        )


                    current_id = (
                        line[1:]
                        .split()[0]
                    )

                    sequence_parts = []


                else:

                    if (
                        current_id
                        in candidate_ids
                    ):

                        sequence_parts.append(
                            line
                        )


            # Final FASTA record.
            if current_id is not None:

                process_record(
                    current_id,
                    "".join(
                        sequence_parts
                    )
                )


# ============================================================
# STEP 7 — VERIFY ALL SELECTED CANDIDATES WERE FOUND
# ============================================================

print("\n============================================================")
print("STEP 7 — Checking candidate processing")
print("============================================================")

print(
    f"Candidates requested: "
    f"{len(candidate_ids):,}"
)

print(
    f"Candidates processed: "
    f"{processed:,}"
)

print(
    f"Alignment errors: "
    f"{errors:,}"
)


if processed != len(candidate_ids):

    print(
        "\nWARNING:"
        "\nNot every candidate locus was encountered in the FASTA."
    )


# ============================================================
# STEP 8 — COMPRESS FINAL ALIGNMENT TABLE
# ============================================================

print("\n============================================================")
print("STEP 8 — Compressing structural results")
print("============================================================")


results_df = pd.read_csv(
    ALIGNMENT_CSV
)


results_df.to_csv(
    ALIGNMENT_GZ,
    index=False,
    compression="gzip"
)


print(
    f"\nSaved compressed results:\n"
    f"{ALIGNMENT_GZ}"
)


# ============================================================
# STEP 9 — CREATE SUMMARY WORKBOOK
# ============================================================

print("\n============================================================")
print("STEP 9 — Creating structural summary workbook")
print("============================================================")


summary_rows = [
    ["Run mode", RUN_MODE.upper()],
    [
        "Detailed alignment minimum locus length (bp)",
        MIN_LOCUS_LENGTH_FOR_ALIGNMENT,
    ],
    [
        "Total loci in Script 02 manifest",
        len(manifest),
    ],
    [
        "Loci meeting full alignment criterion before TEST limit",
        int(
            (
                (manifest["FASTA_QC"] == "PASS")
                &
                (
                    manifest["Coordinate_Length"]
                    >= MIN_LOCUS_LENGTH_FOR_ALIGNMENT
                )
            ).sum()
        ),
    ],
    [
        "Loci requested in this run",
        len(candidate_ids),
    ],
    [
        "Loci processed",
        processed,
    ],
    [
        "Alignment errors",
        errors,
    ],
    [
        "L1-RP nucleotide reference",
        "AF148856.1",
    ],
    [
        "L1-RP length (bp)",
        reference_length,
    ],
    [
        "ORF1 coordinates",
        f"{orf1_start}-{orf1_end}",
    ],
    [
        "ORF1 protein length (aa)",
        orf_coordinates["ORF1"][
            "protein_length_aa"
        ],
    ],
    [
        "ORF2 coordinates",
        f"{orf2_start}-{orf2_end}",
    ],
    [
        "ORF2 protein length (aa)",
        orf_coordinates["ORF2"][
            "protein_length_aa"
        ],
    ],
]


summary_df = pd.DataFrame(
    summary_rows,
    columns=["Metric", "Value"]
)


feature_rows = []

for name, info in features.items():

    feature_rows.append(
        [
            name,
            info["start"],
            info["end"],
            feature_length(
                info["start"],
                info["end"]
            ),
            info["source"],
        ]
    )


feature_df = pd.DataFrame(
    feature_rows,
    columns=[
        "Feature",
        "L1RP_Start_1based",
        "L1RP_End_1based",
        "Length_bp",
        "Coordinate_Source",
    ]
)


source_df = pd.DataFrame(
    [
        [
            "AF148856.1 GenBank",
            "ORF1/ORF2 nucleotide coordinates and reference sequence",
            "https://www.ncbi.nlm.nih.gov/nuccore/AF148856.1",
        ],
        [
            "ORF1p domain boundaries",
            "CC aa 52-153; RRM aa 157-252; CTD aa 264-323",
            "https://pmc.ncbi.nlm.nih.gov/articles/PMC5700708/",
        ],
        [
            "ORF2p domain boundaries",
            "EN aa 1-239; RT aa 453-883; CCHC aa 1096-1275",
            "https://pmc.ncbi.nlm.nih.gov/articles/PMC8161598/",
        ],
        [
            "Detailed-alignment triage precedent",
            "Published intact-L1 workflows commonly first select >=5 kb candidate loci before ORF analysis",
            "https://link.springer.com/article/10.1038/s44318-023-00007-y",
        ],
    ],
    columns=[
        "Source",
        "Use_in_Script",
        "URL",
    ]
)


status_counts = (
    results_df["Alignment_Status"]
    .value_counts(
        dropna=False
    )
    .rename_axis(
        "Alignment_Status"
    )
    .reset_index(
        name="Locus_Count"
    )
)


subfamily_counts = (
    results_df["Primary_Subfamily"]
    .value_counts(
        dropna=False
    )
    .rename_axis(
        "Primary_Subfamily"
    )
    .reset_index(
        name="Locus_Count"
    )
)


with pd.ExcelWriter(
    SUMMARY_XLSX,
    engine="openpyxl"
) as writer:

    summary_df.to_excel(
        writer,
        sheet_name="Summary",
        index=False
    )

    feature_df.to_excel(
        writer,
        sheet_name="Reference_Features",
        index=False
    )

    source_df.to_excel(
        writer,
        sheet_name="Sources",
        index=False
    )

    status_counts.to_excel(
        writer,
        sheet_name="Alignment_Status",
        index=False
    )

    subfamily_counts.to_excel(
        writer,
        sheet_name="Subfamily_Counts",
        index=False
    )

    # TEST runs are intentionally small enough to place the
    # complete alignment output in the workbook for inspection.
    if RUN_MODE.upper() == "TEST":

        results_df.to_excel(
            writer,
            sheet_name="TEST_Results",
            index=False
        )


print(
    f"\nSaved summary workbook:\n"
    f"{SUMMARY_XLSX}"
)


# ============================================================
# FINAL REPORT
# ============================================================

elapsed = (
    time.time()
    - start_time
)


print("\n============================================================")
print("SCRIPT 03 COMPLETE")
print("============================================================")

print(
    f"\nElapsed time: "
    f"{elapsed/60:.2f} minutes"
)

print(
    f"\nProcessed loci: "
    f"{processed:,}"
)

print(
    f"Alignment errors: "
    f"{errors:,}"
)

print(
    "\nStructural results:"
)

print(
    ALIGNMENT_GZ
)

print(
    "\nSummary:"
)

print(
    SUMMARY_XLSX
)


if RUN_MODE.upper() == "TEST":

    print(
        "\nIMPORTANT:"
        "\nThis was a TEST run."
        "\nReview the TEST summary before changing:"
        '\n\n    RUN_MODE = "FULL"\n'
    )


print(
    "\nINTERPRETATION:"
    "\nCoverage columns describe how much of each L1-RP "
    "reference region is represented by the alignment."
)

print(
    "\nThey do NOT by themselves establish that a locus "
    "is transcriptionally active or retrotransposition-competent."
)
