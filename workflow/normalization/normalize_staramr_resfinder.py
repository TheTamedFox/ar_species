# Normalize StarAMR/ResFinder detections using argNorm.
# Preserve the original StarAMR columns and append CARD ARO annotations.

import sys
from pathlib import Path
import subprocess
import tempfile

import pandas as pd

# Read input and output paths supplied by the Snakemake rule.
input_file = Path(sys.argv[1])
output_file = Path(sys.argv[2])

# Create the output directory if necessary.
output_file.parent.mkdir(parents=True, exist_ok=True)

# Load the original StarAMR ResFinder results.
df = pd.read_csv(input_file, sep="\t")

# Define the annotation columns expected from argNorm.
# These provide ARO identifiers, resistance associations and drug classes.
argnorm_columns = [
    "ARO",
    "ARO_name",
    "Cut_Off",
    "confers_resistance_to",
    "confers_resistance_to_names",
    "resistance_to_drug_classes",
    "resistance_to_drug_classes_names",
]

# Handle samples with no ResFinder detections.
# Preserve the original table structure and add empty annotation columns without calling argNorm on an empty input.
if df.empty:
    for column in argnorm_columns:
        df[column] = pd.Series(dtype="object")

    df.to_csv(output_file, sep="\t", index=False)
    print(f"No ResFinder detections: {input_file.name}")
    sys.exit(0)

# Rename StarAMR columns to the names required by argNorm's
# ResFinder input format, without modifying the original DataFrame.
argnorm_input = df.rename(columns={
    "Gene": "Resistance gene",
    "Accession": "Accession no.",
})

# Use temporary files for argNorm processing.
# The temporary directory and its contents are removed automatically.
with tempfile.TemporaryDirectory() as temp_dir:
    temp_input = Path(temp_dir) / "input.tsv"
    argnorm_output = Path(temp_dir) / "output.tsv"

    # Write the reformatted ResFinder table for argNorm.
    argnorm_input.to_csv(temp_input, sep="\t", index=False)

    # Run argNorm in its dedicated Conda environment.
    # check=True stops execution if the command fails.
    subprocess.run(
        [
            "conda", "run", "--no-capture-output",
            "-n", "argnorm",
            "argnorm", "resfinder",
            "-i", str(temp_input),
            "-o", str(argnorm_output),
        ],
        check=True,
    )

    # Load argNorm annotations, ignoring comment lines in its output.
    annotations = pd.read_csv(
        argnorm_output,
        sep="\t",
        comment="#",
    )

    # Verify that all expected annotation columns are present.
    for column in argnorm_columns:
        if column not in annotations.columns:
            raise ValueError(f"argNorm output is missing expected column: {column}")

    # Append argNorm annotations to the original StarAMR results.
    # IMPORTANT: This assumes argNorm preserves the input row order.
    result = pd.concat(
        [
            df.reset_index(drop=True),
            annotations[argnorm_columns].reset_index(drop=True),
        ],
        axis=1,
    )

    # Save the enriched results as a tab-separated file.
    result.to_csv(output_file, sep="\t", index=False)

# Report the number of processed detections.
print(f"Normalised {len(result)} detections: {output_file}")