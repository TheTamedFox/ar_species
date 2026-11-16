# Generate per-sample Excel workbooks comparing AMR detections across StarAMR/ResFinder, RGI contig and RGI read-based analyses.
# Preserve detailed source evidence alongside the comparison results.

import argparse
import re
from pathlib import Path

import pandas as pd
from openpyxl import load_workbook
from openpyxl.utils import get_column_letter


# Read the common-format input directory and report output directory supplied by the Snakemake workflow.
parser = argparse.ArgumentParser()
parser.add_argument("input_dir", type=Path)
parser.add_argument("output_dir", type=Path)
args = parser.parse_args()

# Create the main report output directory if necessary.
args.output_dir.mkdir(parents=True, exist_ok=True)

# Find all common-format TSV files in the input directory and subfolders.
files = sorted(args.input_dir.rglob("*.tsv"))

# Stop if no input TSV files are available.
if not files:
    raise SystemExit(f"No TSV files found in {args.input_dir}")

frames = []

# Load each TSV while preserving identifiers and empty fields as strings.
# Files containing no detection rows are excluded from the combined table.
for file in files:
    df = pd.read_csv(
        file,
        sep="\t",
        dtype=str,
        keep_default_na=False
    )

    if not df.empty:
        frames.append(df)

# Stop if all input tables contain zero detections.
if not frames:
    raise SystemExit("The input TSV files contain no detections.")

# Combine detections from all methods and samples into one evidence table.
evidence = pd.concat(frames, ignore_index=True, sort=False)

# Standardize ARO identifiers to ensure consistent comparisons.
# Preserve empty identifiers and add the ARO: prefix where necessary.
evidence["ARO"] = evidence["ARO"].astype(str).str.strip()
evidence["ARO"] = evidence["ARO"].apply(
    lambda x: x if not x or x.startswith("ARO:") else f"ARO:{x}"
)

# Create a descriptive label combining the detection tool and analysis type.
evidence["Method"] = (
    evidence["Tool"].astype(str)
    + " ("
    + evidence["Analysis"].astype(str)
    + ")"
)

# Define the preferred order of detection methods in comparison tables.
method_order = [
    "StarAMR (ResFinder)",
    "RGI (contig)",
    "RGI (reads)"
]

# Identify samples and generate reports separately for each one.
samples = sorted(evidence["Sample"].unique())

for sample in samples:

    # Select all detections belonging to the current sample.
    sample_evidence = evidence[
        evidence["Sample"] == sample
    ].copy()

    # Separate detections with valid non-empty ARO identifiers
    # from detections that could not be assigned an ARO.
    mapped = sample_evidence[
        sample_evidence["ARO"] != ""
    ].copy()

    unmapped = sample_evidence[
        sample_evidence["ARO"] == ""
    ].copy()

    # Include any additional detection methods beyond the three predefined methods, maintaining a consistent column order.
    additional_methods = sorted(
        set(mapped["Method"]) - set(method_order)
    )

    methods = method_order + additional_methods

    # Build one comparison row for each unique Sample + ARO combination.
    comparison_rows = []

    for (sample_name, aro), group in mapped.groupby(
        ["Sample", "ARO"],
        sort=True
    ):

        row = {
            "Sample": sample_name,
            "ARO": aro
        }

        # Collect distinct non-empty annotations reported for the same ARO.
        # Multiple values are combined into semicolon-separated strings.
        for column in [
            "Gene",
            "ARO_name",
            "Antibiotics",
            "Drug_classes"
        ]:
            values = sorted({
                value.strip()
                for value in group[column].astype(str)
                if value.strip()
            })

            row[column] = "; ".join(values)

        # Determine which methods reported the current ARO.
        detected_methods = set(group["Method"])

        # Mark each method as Detected or Not detected for this ARO.
        for method in methods:
            row[method] = (
                "Detected"
                if method in detected_methods
                else "Not detected"
            )

        comparison_rows.append(row)

    # Convert the ARO-level comparison records into a DataFrame.
    comparison = pd.DataFrame(comparison_rows)

    # Set a consistent column order for the comparison workbook.
    comparison_columns = [
        "Sample",
        "ARO",
        "Gene",
        "ARO_name",
        "Antibiotics",
        "Drug_classes"
    ] + methods

    comparison = comparison.reindex(
        columns=comparison_columns
    )

    # Remove the temporary Method column from the main evidence table.
    # Original tool and analysis labels remain available.
    sample_evidence = sample_evidence.drop(
        columns=["Method"]
    )

    # Replace filename-invalid characters in sample names to create valid output directory names.
    safe_sample = re.sub(
        r'[<>:"/\\|?*]',
        "_",
        str(sample)
    )

    # Create a separate report directory for the current sample.
    sample_dir = args.output_dir / safe_sample
    sample_dir.mkdir(parents=True, exist_ok=True)

    # Define the two Excel workbook paths.
    comparison_file = (
        sample_dir / "AMR_cross_tool_comparison.xlsx"
    )

    evidence_file = (
        sample_dir / "AMR_detection_evidence.xlsx"
    )

    # Write the ARO-level cross-tool comparison workbook.
    with pd.ExcelWriter(
        comparison_file,
        engine="openpyxl"
    ) as writer:

        comparison.to_excel(
            writer,
            sheet_name="Comparison",
            index=False
        )

    # Write the detailed evidence workbook.
    # The Evidence sheet retains all detections, while Unmapped contains detections without an ARO identifier.
    with pd.ExcelWriter(
        evidence_file,
        engine="openpyxl"
    ) as writer:

        sample_evidence.to_excel(
            writer,
            sheet_name="Evidence",
            index=False
        )

        unmapped.to_excel(
            writer,
            sheet_name="Unmapped",
            index=False
        )

    # Apply basic formatting to both Excel workbooks.
    for file in [comparison_file, evidence_file]:

        workbook = load_workbook(file)

        for sheet in workbook.worksheets:

            # Freeze the header row and enable column filtering.
            sheet.freeze_panes = "A2"
            sheet.auto_filter.ref = sheet.dimensions

            # Adjust column widths based on header lengths.
            # Limit widths to keep the workbooks manageable.
            for column in sheet.columns:

                letter = get_column_letter(
                    column[0].column
                )

                header = str(
                    column[0].value or ""
                )

                width = min(
                    max(len(header) + 2, 14),
                    35
                )

                sheet.column_dimensions[letter].width = width

        # Save the formatted Excel workbook.
        workbook.save(file)

    # Print a summary of detections and generated files.
    print(f"\nSample: {sample}")
    print(f"  Total detections: {len(sample_evidence)}")
    print(f"  Mapped detections: {len(mapped)}")
    print(f"  Unmapped detections: {len(unmapped)}")
    print(f"  Comparison rows: {len(comparison)}")
    print(f"  Comparison workbook: {comparison_file}")
    print(f"  Evidence workbook: {evidence_file}")

# Report how many samples received workbooks.
print(f"\nProcessed samples: {len(samples)}")