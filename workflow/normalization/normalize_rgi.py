
# Normalize CARD/RGI contig- and read-based results using argNorm's CARD ontology while preserving all original RGI detection evidence.

import sys
from pathlib import Path

import pandas as pd

from argnorm import lib
from argnorm.drug_categorization import (
    confers_resistance_to,
    drugs_to_drug_classes,
)

# Read the input path, output path and analysis type supplied by Snakemake.
input_file = Path(sys.argv[1])
output_file = Path(sys.argv[2])
result_type = sys.argv[3]

# Create the output directory if necessary.
output_file.parent.mkdir(parents=True, exist_ok=True)

# Load original RGI results as strings to preserve accession identifiers.
df = pd.read_csv(input_file, sep="\t", dtype=str)

# Identify the ARO accession column according to the RGI analysis type.
# Contig- and read-based RGI results use different column names.
if result_type == "contig":
    aro_column = "ARO"
elif result_type == "reads":
    aro_column = "ARO Accession"
else:
    raise ValueError(f"Unknown RGI result type: {result_type}")

# Stop if the expected accession column is missing.
if aro_column not in df.columns:
    raise ValueError(
        f"Expected column '{aro_column}' not found in {input_file.name}"
    )

# Preserve the original RGI accession for traceability.
df["ARO_original"] = df[aro_column]

# Standardize ARO identifiers to the common ARO:xxxxxxx format.
# Missing accessions remain missing; existing prefixes are preserved.
def standardise_aro(value):
    if pd.isna(value) or not str(value).strip():
        return pd.NA

    accession = str(value).strip()

    if accession.startswith("ARO:"):
        return accession

    return f"ARO:{accession}"

df["ARO"] = df["ARO_original"].apply(standardise_aro)

# Load the CARD Antibiotic Resistance Ontology through argNorm.
ontology = lib.get_aro_ontology()

# Cache annotations to avoid repeated ontology lookups when multiple detections share the same ARO accession.
annotation_cache = {}

# Retrieve gene names, resistance-associated antibiotics and drug classes from the CARD ontology using a standardized ARO accession.
def annotate_aro(accession):
    # Return missing annotations when no valid accession is available.
    if pd.isna(accession):
        return {
            "ARO_name": pd.NA,
            "confers_resistance_to": pd.NA,
            "confers_resistance_to_names": pd.NA,
            "resistance_to_drug_classes": pd.NA,
            "resistance_to_drug_classes_names": pd.NA,
        }

    # Reuse existing annotations for previously encountered AROs.
    if accession in annotation_cache:
        return annotation_cache[accession]

    try:
        # Retrieve the CARD ontology term associated with the ARO.
        term = ontology[accession]

        # Identify resistance-associated antibiotics and their drug classes.
        drugs = confers_resistance_to(accession)
        classes = drugs_to_drug_classes(drugs) if drugs else []

        # Store both ontology accessions and human-readable names.
        result = {
            "ARO_name": term.name,
            "confers_resistance_to": ",".join(drugs) if drugs else pd.NA,
            "confers_resistance_to_names": (
                ",".join(ontology[drug].name for drug in drugs)
                if drugs else pd.NA
            ),
            "resistance_to_drug_classes": (
                ",".join(classes) if classes else pd.NA
            ),
            "resistance_to_drug_classes_names": (
                ",".join(ontology[item].name for item in classes)
                if classes else pd.NA
            ),
        }

    # Warn about missing ontology entries without stopping processing.
    except KeyError:
        print(f"WARNING: ARO not found in ontology: {accession}")
        result = {
            "ARO_name": pd.NA,
            "confers_resistance_to": pd.NA,
            "confers_resistance_to_names": pd.NA,
            "resistance_to_drug_classes": pd.NA,
            "resistance_to_drug_classes_names": pd.NA,
        }

    # Cache the result and return it for the current detection.
    annotation_cache[accession] = result
    return result

# Define the ontology annotation fields to add to the original RGI table.
annotation_columns = [
    "ARO_name",
    "confers_resistance_to",
    "confers_resistance_to_names",
    "resistance_to_drug_classes",
    "resistance_to_drug_classes_names",
]

# Annotate each detection using its ARO identifier.
# Preserve the original DataFrame index to align annotations with rows.
annotations = pd.DataFrame(
    df["ARO"].apply(annotate_aro).tolist(),
    index=df.index,
    columns=annotation_columns,
)

# Append ontology annotations without removing original RGI evidence.
df = pd.concat([df, annotations], axis=1)

# Save the enriched results as a tab-separated file.
df.to_csv(output_file, sep="\t", index=False)

# Report the number of detections and unique annotated ARO accessions.
print(
    f"Annotated {len(df)} RGI {result_type} detections "
    f"({len(annotation_cache)} unique AROs): {output_file}"
)