# Convert normalized AMR detection results into a common table format for comparison across StarAMR/ResFinder, RGI contig and RGI reads.
# Preserve all original source columns for detailed evidence reporting.

import re
import sys
from pathlib import Path

import pandas as pd

# Read the input/output paths and detection source labels from Snakemake.
input_file = Path(sys.argv[1])
output_file = Path(sys.argv[2])
tool = sys.argv[3]
analysis = sys.argv[4]

# Create the output directory if necessary.
output_file.parent.mkdir(parents=True, exist_ok=True)

# Read the normalized input table while preserving identifiers as strings and keeping empty fields as empty strings rather than NaN.
df = pd.read_csv(input_file, sep="\t", dtype=str, keep_default_na=False)

# Derive the sample identifier from the normalized input filename.
sample = input_file.name

# Remove the tool-specific filename suffix.
for suffix in [
    "_staramr_argnorm_resfinder.tsv",
    "_rgi_argnorm_contig.tsv",
    "_rgi_argnorm_reads.tsv",
]:
    if sample.endswith(suffix):
        sample = sample[:-len(suffix)]
        break

# Remove Illumina sample and lane suffixes for consistent report names.
# Example: Aci282_S51_L001 -> Aci282
sample = re.sub(r"_S\d+_L\d+$", "", sample)

# Identify a gene name from the available source columns.
# Use the first non-empty value in the listed priority order.
gene_columns = [
    "Gene",
    "ARO Term",
    "AMR Gene Family",
    "Best_Hit_ARO",
]

# Initialize empty gene names for all detections.
gene = pd.Series("", index=df.index, dtype=str)

# Fill missing gene names using alternative source columns.
for column in gene_columns:
    if column in df.columns:
        gene = gene.mask(gene.eq(""), df[column])

# Standardize ARO identifiers so non-empty accessions use the ARO: prefix.
if "ARO" in df.columns:
    aro = df["ARO"].copy()
else:
    aro = pd.Series("", index=df.index, dtype=str)

# Remove surrounding whitespace and add the prefix when necessary.
# Preserve empty values without creating invalid accessions.
aro = aro.str.strip()
aro = aro.where(aro.eq("") | aro.str.startswith("ARO:"), "ARO:" + aro)

# Retrieve an annotation column if present.
# Return empty values if the source does not contain that column.
def get_column(name):
    if name in df.columns:
        return df[name]
    return pd.Series("", index=df.index, dtype=str)


# Create the shared columns used by the cross-tool comparison workflow.
# Tool and Analysis identify which detection method produced each result.
common = pd.DataFrame({
    "Sample": sample,
    "Tool": tool,
    "Analysis": analysis,
    "ARO": aro,
    "Gene": gene,
    "ARO_name": get_column("ARO_name"),
    "Antibiotics": get_column("confers_resistance_to_names"),
    "Drug_classes": get_column("resistance_to_drug_classes_names"),
    "Detection": "Detected",
    "Source_file": input_file.name,
})

# Preserve every original column from the normalized input.
# Prefix source columns to distinguish them from common-format columns.
source = df.add_prefix("source_")

# Combine the shared comparison columns with the original detection evidence.
result = pd.concat(
    [common.reset_index(drop=True), source.reset_index(drop=True)],
    axis=1,
)

# Save the complete common-format table.
result.to_csv(output_file, sep="\t", index=False)

# Report the number of converted detections and the sample identifier.
print(
    f"Converted {len(result)} detections "
    f"from {input_file.name} for sample {sample}"
)