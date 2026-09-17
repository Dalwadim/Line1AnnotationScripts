# ============================================================
# LINE-1 FUNCTIONAL ANNOTATION PROJECT
#
# SCRIPT 02:
# Prepare and QC the hg38 LINE-1 locus dataset
#
# INPUTS:
#   1. hg38_RE_only_L1.txt
#   2. hg38_RE_L1_seqs.zip
#
# IMPORTANT PROJECT ASSUMPTIONS:
#   - Each Locus_ID should be unique.
#   - Coordinates are 1-based and inclusive.
#   - Overlapping/nested/subset LINE-1 records are NOT merged.
#   - FASTA sequence is used here for QC/length only.
#   - This script does NOT classify loci as active/inactive.
# ============================================================

import gc
import io
import re
import zipfile
from pathlib import Path

import numpy as np
import pandas as pd


# ============================================================
#                  EDIT THESE SETTINGS
# ============================================================

# If these files are in the SAME folder as this script,
# you do not need to change these names.

LOCUS_TABLE_NAME = "hg38_RE_only_L1.txt"
FASTA_ZIP_NAME = "hg38_RE_L1_seqs.zip"

# Output folder. The script will create it automatically.
OUTPUT_FOLDER_NAME = "02_prepared_LINE1_data"

# Full prepared manifest.
# A compressed CSV is used because the dataset has ~1 million rows.
MANIFEST_NAME = "hg38_LINE1_locus_manifest.csv.gz"

# Small Excel workbook containing summary/QC information.
SUMMARY_NAME = "hg38_LINE1_preparation_summary.xlsx"


# ============================================================
#               DO NOT EDIT BELOW HERE
# ============================================================

SCRIPT_FOLDER = Path(__file__).resolve().parent

LOCUS_TABLE = SCRIPT_FOLDER / LOCUS_TABLE_NAME
FASTA_ZIP = SCRIPT_FOLDER / FASTA_ZIP_NAME

OUTPUT_FOLDER = SCRIPT_FOLDER / OUTPUT_FOLDER_NAME
OUTPUT_FOLDER.mkdir(parents=True, exist_ok=True)

MANIFEST_FILE = OUTPUT_FOLDER / MANIFEST_NAME
SUMMARY_FILE = OUTPUT_FOLDER / SUMMARY_NAME

MISSING_FASTA_FILE = OUTPUT_FOLDER / "LINE1_missing_FASTA_IDs.txt"
EXTRA_FASTA_FILE = OUTPUT_FOLDER / "LINE1_extra_FASTA_IDs.txt"
DUPLICATE_LOCUS_FILE = OUTPUT_FOLDER / "LINE1_duplicate_locus_IDs.txt"
DUPLICATE_FASTA_FILE = OUTPUT_FOLDER / "LINE1_duplicate_FASTA_IDs.txt"


# ============================================================
# BASIC FILE CHECKS
# ============================================================

if not LOCUS_TABLE.exists():
    raise FileNotFoundError(
        f"\nCould not find:\n{LOCUS_TABLE}\n\n"
        f"Place '{LOCUS_TABLE_NAME}' in the same folder as this script."
    )

if not FASTA_ZIP.exists():
    raise FileNotFoundError(
        f"\nCould not find:\n{FASTA_ZIP}\n\n"
        f"Place '{FASTA_ZIP_NAME}' in the same folder as this script."
    )


print("\n============================================================")
print("SCRIPT 02 — PREPARE hg38 LINE-1 LOCUS DATASET")
print("============================================================")

print("\nLocus table:")
print(LOCUS_TABLE)

print("\nFASTA ZIP:")
print(FASTA_ZIP)

print("\nCoordinate convention:")
print("1-based inclusive")


# ============================================================
# STEP 1 — LOAD THE LOCUS TABLE
# ============================================================

print("\n============================================================")
print("STEP 1 — Loading locus annotation table")
print("============================================================")

column_names = [
    "Locus_ID",
    "Chromosome",
    "Start_1based",
    "End_1based",
    "Strand",
    "Original_Length_Field",
    "Repeat_Annotation",
]

df = pd.read_csv(
    LOCUS_TABLE,
    sep="\t",
    header=None,
    names=column_names,
    dtype={
        "Locus_ID": "string",
        "Chromosome": "string",
        "Start_1based": "int64",
        "End_1based": "int64",
        "Strand": "string",
        "Original_Length_Field": "int64",
        "Repeat_Annotation": "string",
    },
    keep_default_na=False,
    na_filter=False,
)

print(f"Rows loaded: {len(df):,}")


# ============================================================
# STEP 2 — VALIDATE LOCUS IDs
# ============================================================

print("\n============================================================")
print("STEP 2 — Checking unique Locus_ID values")
print("============================================================")

duplicate_mask = df["Locus_ID"].duplicated(keep=False)
duplicate_locus_ids = (
    df.loc[duplicate_mask, "Locus_ID"]
    .drop_duplicates()
    .sort_values()
)

if len(duplicate_locus_ids) > 0:

    duplicate_locus_ids.to_csv(
        DUPLICATE_LOCUS_FILE,
        index=False,
        header=False
    )

    raise ValueError(
        f"\nSTOPPED: Found {len(duplicate_locus_ids):,} duplicated "
        f"Locus_ID values.\n"
        f"They were written to:\n{DUPLICATE_LOCUS_FILE}\n\n"
        f"No records were automatically merged or deleted."
    )

print("PASS: Every Locus_ID is unique.")


# ============================================================
# STEP 3 — COORDINATE QC
# ============================================================

print("\n============================================================")
print("STEP 3 — Coordinate QC")
print("============================================================")

# User confirmed that the coordinates are 1-based inclusive.
df["Coordinate_Length"] = (
    df["End_1based"] - df["Start_1based"] + 1
)

# Preserve the source length exactly as supplied.
# We do NOT overwrite it.
df["Coordinate_minus_Original_Length"] = (
    df["Coordinate_Length"] - df["Original_Length_Field"]
)

invalid_coordinate_mask = (
    (df["Start_1based"] < 1)
    | (df["End_1based"] < df["Start_1based"])
)

invalid_coordinate_count = int(invalid_coordinate_mask.sum())

print(f"Invalid coordinate rows: {invalid_coordinate_count:,}")

if invalid_coordinate_count > 0:
    print(
        "WARNING: Invalid coordinate rows were found. "
        "They will remain in the manifest for review."
    )


# ============================================================
# STEP 4 — PARSE LINE-1 SUBFAMILY INFORMATION
# ============================================================

print("\n============================================================")
print("STEP 4 — Parsing repeat/subfamily annotations")
print("============================================================")

# Handles standard annotations such as:
#   L1M2|L1|...
#
# and composite annotations such as:
#   L1PA4|LINE/L1 L1MB7|LINE/L1

subfamily_pattern = re.compile(
    r"([A-Za-z0-9_.-]+)\|(?:L1|LINE/L1)"
)


def parse_subfamilies(annotation):
    annotation = str(annotation)

    hits = subfamily_pattern.findall(annotation)

    # Preserve order while removing duplicate names.
    hits = list(dict.fromkeys(hits))

    if not hits:
        # Fallback: preserve the first annotation token.
        fallback = annotation.split("|")[0].strip()
        hits = [fallback] if fallback else ["UNKNOWN"]

    primary = hits[0]
    all_subfamilies = ";".join(hits)

    return primary, all_subfamilies, len(hits)


parsed = df["Repeat_Annotation"].apply(parse_subfamilies)

df["Primary_Subfamily"] = parsed.str[0]
df["All_Subfamilies"] = parsed.str[1]
df["Subfamily_Count"] = parsed.str[2].astype("int16")

df["Composite_Subfamily_Annotation"] = (
    df["Subfamily_Count"] > 1
)

print(
    "Composite/multi-subfamily annotations: "
    f"{int(df['Composite_Subfamily_Annotation'].sum()):,}"
)


# ============================================================
# STEP 5 — READ FASTA LENGTHS DIRECTLY FROM ZIP
# ============================================================

print("\n============================================================")
print("STEP 5 — Reading FASTA records directly from ZIP")
print("============================================================")


def read_fasta_lengths_from_zip(zip_path):
    """
    Stream the FASTA directly from the ZIP archive.

    We record only:
        FASTA ID -> sequence length

    The full nucleotide sequences are NOT stored in memory.
    """

    fasta_lengths = {}
    duplicate_fasta_ids = set()

    with zipfile.ZipFile(zip_path, "r") as zf:

        fasta_members = [
            name for name in zf.namelist()
            if name.lower().endswith(
                (".fa", ".fasta", ".fna", ".fas")
            )
        ]

        if len(fasta_members) == 0:
            raise ValueError(
                "No FASTA file was found inside the ZIP archive."
            )

        if len(fasta_members) > 1:
            raise ValueError(
                "More than one FASTA file was found inside the ZIP:\n"
                + "\n".join(fasta_members)
                + "\n\nPlease choose the intended FASTA explicitly."
            )

        fasta_member = fasta_members[0]

        print(f"FASTA inside ZIP: {fasta_member}")

        with zf.open(fasta_member, "r") as raw_handle:
            handle = io.TextIOWrapper(
                raw_handle,
                encoding="utf-8",
                errors="replace"
            )

            current_id = None
            current_length = 0
            record_count = 0

            for raw_line in handle:
                line = raw_line.strip()

                if not line:
                    continue

                if line.startswith(">"):

                    # Save previous record.
                    if current_id is not None:

                        if current_id in fasta_lengths:
                            duplicate_fasta_ids.add(current_id)
                        else:
                            fasta_lengths[current_id] = current_length

                        record_count += 1

                        if record_count % 100000 == 0:
                            print(
                                f"  FASTA records read: "
                                f"{record_count:,}"
                            )

                    # FASTA identifier = text after ">" up to first space.
                    current_id = line[1:].split()[0]
                    current_length = 0

                else:
                    current_length += len(line)

            # Save final FASTA record.
            if current_id is not None:

                if current_id in fasta_lengths:
                    duplicate_fasta_ids.add(current_id)
                else:
                    fasta_lengths[current_id] = current_length

                record_count += 1

    return fasta_lengths, duplicate_fasta_ids, fasta_member


fasta_lengths, duplicate_fasta_ids, fasta_member = (
    read_fasta_lengths_from_zip(FASTA_ZIP)
)

print(f"\nFASTA records loaded: {len(fasta_lengths):,}")


if duplicate_fasta_ids:

    with open(
        DUPLICATE_FASTA_FILE,
        "w",
        encoding="utf-8"
    ) as handle:

        for fasta_id in sorted(duplicate_fasta_ids):
            handle.write(f"{fasta_id}\n")

    raise ValueError(
        f"\nSTOPPED: Found {len(duplicate_fasta_ids):,} duplicated "
        f"FASTA IDs.\n"
        f"They were written to:\n{DUPLICATE_FASTA_FILE}"
    )

print("PASS: FASTA IDs are unique.")


# ============================================================
# STEP 6 — JOIN FASTA QC TO LOCUS TABLE
# ============================================================

print("\n============================================================")
print("STEP 6 — Matching FASTA records to locus IDs")
print("============================================================")

df["FASTA_Length"] = (
    df["Locus_ID"]
    .map(fasta_lengths)
    .astype("Int64")
)

df["FASTA_Present"] = df["FASTA_Length"].notna()

df["FASTA_Length_Match"] = (
    df["FASTA_Present"]
    & (
        df["FASTA_Length"].astype("Float64")
        == df["Coordinate_Length"].astype("Float64")
    )
)


def fasta_qc_status(row):

    if not row["FASTA_Present"]:
        return "MISSING_FASTA"

    if row["FASTA_Length_Match"]:
        return "PASS"

    return "LENGTH_MISMATCH"


df["FASTA_QC"] = df.apply(
    fasta_qc_status,
    axis=1
)


# Missing annotation -> FASTA relationships.
missing_fasta_ids = (
    df.loc[
        ~df["FASTA_Present"],
        "Locus_ID"
    ]
    .astype(str)
    .tolist()
)

with open(
    MISSING_FASTA_FILE,
    "w",
    encoding="utf-8"
) as handle:

    for locus_id in missing_fasta_ids:
        handle.write(f"{locus_id}\n")


# Identify FASTA IDs not represented in the annotation table.
# The set is created only for this QC comparison.
annotation_id_set = set(
    df["Locus_ID"].astype(str).tolist()
)

extra_fasta_ids = sorted(
    fasta_id
    for fasta_id in fasta_lengths.keys()
    if fasta_id not in annotation_id_set
)

with open(
    EXTRA_FASTA_FILE,
    "w",
    encoding="utf-8"
) as handle:

    for fasta_id in extra_fasta_ids:
        handle.write(f"{fasta_id}\n")


missing_count = len(missing_fasta_ids)
extra_count = len(extra_fasta_ids)

length_mismatch_count = int(
    (df["FASTA_QC"] == "LENGTH_MISMATCH").sum()
)

pass_count = int(
    (df["FASTA_QC"] == "PASS").sum()
)

print(f"FASTA PASS:             {pass_count:,}")
print(f"Missing FASTA records:  {missing_count:,}")
print(f"Extra FASTA records:    {extra_count:,}")
print(f"Length mismatches:      {length_mismatch_count:,}")


# We no longer need the ~1-million-entry FASTA dictionary.
del annotation_id_set
del fasta_lengths
gc.collect()


# ============================================================
# STEP 7 — FLAG NESTED / SUBSET LOCI
# ============================================================

print("\n============================================================")
print("STEP 7 — Flagging coordinate-contained LINE-1 records")
print("============================================================")

print(
    "NOTE: These are spatial containment flags only.\n"
    "Records are NOT merged, deleted, or assigned biological "
    "parent-child relationships."
)

n_rows = len(df)

contained_within = np.zeros(
    n_rows,
    dtype=bool
)

contains_other = np.zeros(
    n_rows,
    dtype=bool
)

starts = df["Start_1based"].to_numpy(dtype=np.int64)
ends = df["End_1based"].to_numpy(dtype=np.int64)

# Factorize chromosome values for efficient sorting.
chrom_codes, chrom_names = pd.factorize(
    df["Chromosome"],
    sort=True
)

# Sort by:
#   chromosome
#   start ascending
#   end descending
#
# End descending is important when multiple intervals share
# the same start: the larger containing interval comes first.
order = np.lexsort(
    (
        -ends,
        starts,
        chrom_codes
    )
)


# ---- Forward scan: Is this interval contained in an earlier one? ----

current_chrom = None
max_end = -1

for idx in order:

    chrom = chrom_codes[idx]
    end = ends[idx]

    if chrom != current_chrom:

        current_chrom = chrom
        max_end = end
        continue

    if end <= max_end:
        contained_within[idx] = True

    if end > max_end:
        max_end = end


# ---- Reverse scan: Does this interval contain a later one? ----

current_chrom = None
min_end = np.iinfo(np.int64).max

for idx in order[::-1]:

    chrom = chrom_codes[idx]
    end = ends[idx]

    if chrom != current_chrom:

        current_chrom = chrom
        min_end = end
        continue

    if end >= min_end:
        contains_other[idx] = True

    if end < min_end:
        min_end = end


df["Contained_Within_Other_L1"] = contained_within
df["Contains_Other_L1"] = contains_other

print(
    "Contained within another L1 interval: "
    f"{int(contained_within.sum()):,}"
)

print(
    "Contains at least one other L1 interval: "
    f"{int(contains_other.sum()):,}"
)


# ============================================================
# STEP 8 — ADD DESCRIPTIVE LENGTH BINS
# ============================================================

print("\n============================================================")
print("STEP 8 — Adding descriptive length bins")
print("============================================================")

# These are descriptive bins only.
# They do NOT imply functional activity.

df["Length_Bin"] = pd.cut(
    df["Coordinate_Length"],
    bins=[
        -np.inf,
        999,
        4999,
        5999,
        np.inf
    ],
    labels=[
        "<1 kb",
        "1-4.999 kb",
        "5-5.999 kb",
        ">=6 kb"
    ]
)


# ============================================================
# STEP 9 — REORDER OUTPUT COLUMNS
# ============================================================

output_columns = [
    "Locus_ID",
    "Chromosome",
    "Start_1based",
    "End_1based",
    "Strand",

    "Original_Length_Field",
    "Coordinate_Length",
    "Coordinate_minus_Original_Length",
    "Length_Bin",

    "Repeat_Annotation",
    "Primary_Subfamily",
    "All_Subfamilies",
    "Subfamily_Count",
    "Composite_Subfamily_Annotation",

    "FASTA_Present",
    "FASTA_Length",
    "FASTA_Length_Match",
    "FASTA_QC",

    "Contained_Within_Other_L1",
    "Contains_Other_L1",
]

df = df[output_columns]


# ============================================================
# STEP 10 — BUILD SUMMARY TABLES
# ============================================================

print("\n============================================================")
print("STEP 10 — Building QC summaries")
print("============================================================")

summary_rows = [
    ["Coordinate convention", "1-based inclusive"],
    ["Total annotation records", len(df)],
    ["Unique Locus_ID values", df["Locus_ID"].nunique()],
    ["FASTA file inside ZIP", fasta_member],
    ["FASTA records matched / PASS", pass_count],
    ["Missing FASTA records", missing_count],
    ["Extra FASTA records", extra_count],
    ["FASTA length mismatches", length_mismatch_count],
    ["Invalid coordinate rows", invalid_coordinate_count],
    [
        "Composite/multi-subfamily annotations",
        int(df["Composite_Subfamily_Annotation"].sum())
    ],
    [
        "Contained within another L1 interval",
        int(df["Contained_Within_Other_L1"].sum())
    ],
    [
        "Contains another L1 interval",
        int(df["Contains_Other_L1"].sum())
    ],
    [
        "Coordinate length >= 5 kb",
        int((df["Coordinate_Length"] >= 5000).sum())
    ],
    [
        "Coordinate length >= 6 kb",
        int((df["Coordinate_Length"] >= 6000).sum())
    ],
]

summary_df = pd.DataFrame(
    summary_rows,
    columns=["Metric", "Value"]
)


subfamily_counts = (
    df["Primary_Subfamily"]
    .value_counts(dropna=False)
    .rename_axis("Primary_Subfamily")
    .reset_index(name="Locus_Count")
)


fasta_qc_counts = (
    df["FASTA_QC"]
    .value_counts(dropna=False)
    .rename_axis("FASTA_QC")
    .reset_index(name="Locus_Count")
)


length_bin_counts = (
    df["Length_Bin"]
    .value_counts(sort=False, dropna=False)
    .rename_axis("Length_Bin")
    .reset_index(name="Locus_Count")
)


length_difference_counts = (
    df["Coordinate_minus_Original_Length"]
    .value_counts(dropna=False)
    .sort_index()
    .rename_axis("Coordinate_minus_Original_Length")
    .reset_index(name="Locus_Count")
)


fasta_issue_df = df.loc[
    df["FASTA_QC"] != "PASS",
    [
        "Locus_ID",
        "Chromosome",
        "Start_1based",
        "End_1based",
        "Strand",
        "Coordinate_Length",
        "FASTA_Length",
        "FASTA_QC",
        "Primary_Subfamily",
    ]
].copy()


# ============================================================
# STEP 11 — SAVE FULL MANIFEST
# ============================================================

print("\n============================================================")
print("STEP 11 — Saving prepared LINE-1 manifest")
print("============================================================")

print(
    "Writing compressed CSV.\n"
    "This can take a few minutes because the dataset contains "
    "nearly one million loci."
)

df.to_csv(
    MANIFEST_FILE,
    index=False,
    compression="gzip"
)

print(f"\nSaved full manifest:\n{MANIFEST_FILE}")


# ============================================================
# STEP 12 — SAVE SMALL EXCEL QC SUMMARY
# ============================================================

print("\nSaving QC summary workbook...")

with pd.ExcelWriter(
    SUMMARY_FILE,
    engine="openpyxl"
) as writer:

    summary_df.to_excel(
        writer,
        sheet_name="Summary",
        index=False
    )

    subfamily_counts.to_excel(
        writer,
        sheet_name="Subfamily_Counts",
        index=False
    )

    fasta_qc_counts.to_excel(
        writer,
        sheet_name="FASTA_QC_Counts",
        index=False
    )

    length_bin_counts.to_excel(
        writer,
        sheet_name="Length_Bins",
        index=False
    )

    length_difference_counts.to_excel(
        writer,
        sheet_name="Length_Field_QC",
        index=False
    )

    fasta_issue_df.to_excel(
        writer,
        sheet_name="FASTA_Issues",
        index=False
    )


print(f"Saved QC workbook:\n{SUMMARY_FILE}")


# ============================================================
# FINAL REPORT
# ============================================================

print("\n============================================================")
print("SCRIPT 02 COMPLETE")
print("============================================================")

print(f"\nTotal loci retained: {len(df):,}")
print("No overlapping or nested records were automatically merged.")

print("\nPrimary output:")
print(MANIFEST_FILE)

print("\nQC summary:")
print(SUMMARY_FILE)

print("\nMissing FASTA ID list:")
print(MISSING_FASTA_FILE)

print("\nExtra FASTA ID list:")
print(EXTRA_FASTA_FILE)

print(
    "\nNEXT PIPELINE STEP:\n"
    "Use this prepared manifest to decide which loci proceed "
    "to structural comparison against L1-RP."
)
