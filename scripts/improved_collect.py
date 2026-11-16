#!/usr/bin/env python3

# Combine StarAMR/ResFinder detections with Kraken2 taxonomy and BLAST hits.
# This script preserves the existing identifiers used by the taxonomy branch.
import sys
import os
from typing import Dict


# Build lookup tables from the NCBI taxonomy database.
# These are used to translate Kraken2 tax IDs into names and ranks.
def load_taxonomy_data(taxonomy_dir):
    """
    Load taxonomy names, parents, and ranks from NCBI taxonomy files.

    Args:
        taxonomy_dir (str): Path to directory containing nodes.dmp and names.dmp

    Returns:
        tuple: Dictionaries for tax ID to name, parent, and rank
    """
    try:
        taxonomy_nodes = os.path.join(taxonomy_dir, "nodes.dmp")
        taxonomy_names = os.path.join(taxonomy_dir, "names.dmp")

        # Index taxonomy names, parent relationships and taxonomic ranks by tax ID.
        tax_id_to_name: Dict[str, str] = {}
        tax_id_to_parent: Dict[str, str] = {}
        tax_id_to_rank: Dict[str, str] = {}

        # Load taxon names from names.dmp.
        # Review later: this currently does not select only scientific names.
        # Load names
        with open(taxonomy_names, "r") as names_file:
            for line in names_file:
                parts = line.strip().split("\t|\t")
                tax_id = parts[0].strip()
                name = parts[1].strip()
                tax_id_to_name[tax_id] = name

        # Read parent relationships and ranks from nodes.dmp.
        # Load nodes
        with open(taxonomy_nodes, "r") as nodes_file:
            for line in nodes_file:
                parts = line.strip().split("\t|\t")
                tax_id = parts[0].strip()
                parent_id = parts[1].strip()
                rank = parts[2].strip()
                tax_id_to_parent[tax_id] = parent_id
                tax_id_to_rank[tax_id] = rank

        return tax_id_to_name, tax_id_to_parent, tax_id_to_rank

    except IOError as e:
        print(f"Error reading taxonomy files: {e}")
        sys.exit(1)


# Follow parent tax IDs upward to build a readable taxonomic lineage.
def get_full_lineage(tax_id,
                     tax_id_to_name,
                     tax_id_to_parent,
                     tax_id_to_rank):
    """
    Construct full taxonomic lineage for a given tax ID.

    Args:
        tax_id (str): Taxonomy ID to trace
        tax_id_to_name (Dict): Mapping of tax IDs to names
        tax_id_to_parent (Dict): Mapping of tax IDs to parent tax IDs
        tax_id_to_rank (Dict): Mapping of tax IDs to taxonomic ranks

    Returns:
        str: Taxonomic lineage string
    """
    # Restrict the displayed lineage to selected taxonomic ranks.
    RANKS_OF_INTEREST = ["species", "genus", "family", "order", "class", "phylum", "kingdom", "superkingdom"]

    lineage = []
    current_id = tax_id
    seen_ids = set()

    # Walk through parent nodes; track visited IDs to avoid a cycle.
    while current_id in tax_id_to_parent:
        if current_id in seen_ids:
            break
        seen_ids.add(current_id)

        rank = tax_id_to_rank.get(current_id, "Unknown")
        name = tax_id_to_name.get(current_id, "Unknown Taxon")

        # Current behavior: stop at root OR the first unranked node.
        # Review later: stopping at an unranked node can truncate lineage.
        if current_id == "1" or rank == "no rank":
            lineage.append(name)
            break

        if rank in RANKS_OF_INTEREST:
            lineage.append(name)

        current_id = tax_id_to_parent[current_id]

    return " > ".join(reversed(lineage)) if lineage else "Unknown Taxon"


# Convert Kraken2 sequence classifications into readable taxonomy records.
def translate_kraken(kraken2_output, kraken2_taxonomy):
    """
    Process Kraken2 output and translate taxonomy IDs.

    Args:
        kraken2_output (str): Path to Kraken2 output file
        kraken2_taxonomy (str): Path to taxonomy directory

    Returns:
        Dict: Processed taxonomy dictionary
    """
    # Load taxonomy data
    tax_id_to_name, tax_id_to_parent, tax_id_to_rank = load_taxonomy_data(kraken2_taxonomy)

    tax_dict = {}

    try:
        with open(kraken2_output, "r") as infile:
            for line in infile:
                # Kraken2 output provides status, sequence ID and taxonomy ID.
                parts = line.strip().split("\t", maxsplit=4)
                classification = "Classified" if parts[0] == "C" else "Unclassified"
                seq_id = parts[1]
                tax_id = parts[2]

                # Keep sequence identifiers unchanged so they match FASTA/BLAST IDs.
                if classification == "Classified":
                    rank = tax_id_to_rank.get(tax_id, "Unknown")
                    full_taxonomic_name = get_full_lineage(tax_id, tax_id_to_name, tax_id_to_parent, tax_id_to_rank)

                    tax_dict[seq_id] = [seq_id, classification, tax_id, full_taxonomic_name, rank]
                else:
                    tax_dict[seq_id] = [seq_id, classification, "-", "Unclassified", "-"]

    except IOError as e:
        print(f"Error reading Kraken2 output: {e}")
        sys.exit(1)

    return tax_dict


# Group BLAST hit descriptions by query sequence identifier.
def read_blast_result(blast_output):
    """
    Read and process BLAST results.

    Args:
        blast_output (str): Path to BLAST output file

    Returns:
        Dict: BLAST results with sequence IDs as keys
    """
    results = {}

    try:
        with open(blast_output, "r") as blast:
            for line in blast:
                # Extract the query ID and hit description from the configured BLAST output.
                seq_id = line.split("\t")[0]
                hit = line.split("\t")[2]

                if seq_id not in results:
                    results[seq_id] = [hit]
                else:
                    results[seq_id].append(hit)

    except IOError as e:
        print(f"Error reading BLAST output: {e}")
        sys.exit(1)

    return results


# Merge ResFinder detections with taxonomic classifications and BLAST hits.
def collect_results(resistances,
                    kraken2_output,
                    kraken2_taxonomy,
                    blast_output,
                    collected_results):
    """
    Collect and combine results from multiple sources.

    Args:
        resistances (str): Path to resistance gene file
        kraken2_output (str): Path to Kraken2 output
        kraken2_taxonomy (str): Path to taxonomy directory
        blast_output (str): Path to BLAST output
        collected_results (str): Path to output results file
    """
    # Store resistance detections and derive the same shortened sample prefix
    # used by the current scaffold-selection script.
    entries = {}
    name = os.path.basename(os.path.dirname(resistances)).split("_")[0]

    # Process resistance genes
    try:
        with open(resistances, "r") as res:
            for i, line in enumerate(res):
                if i == 0:  # Skip header
                    continue
                comp = line.split("\t")
                # Skip unknown phenotypes and keep the validated ResFinder fields:
                # sample, gene, phenotype, identity, overlap, contig, accession.
                if "unknown" not in comp[2]:
                    entries[str(i - 1)] = [name, comp[1], comp[2], comp[4], comp[5], comp[7], comp[10].rstrip()]

    except IOError as e:
        print(f"Error reading resistance file: {e}")
        sys.exit(1)

    # Process and write combined results
    try:
        # Read taxonomy and BLAST results, keyed by scaffold FASTA identifier.
        kraken_results = translate_kraken(kraken2_output, kraken2_taxonomy)
        blast_results = read_blast_result(blast_output)

        with open(collected_results, "w") as collection:
            # Write the combined report header.
            headers = ["Sample", "Gene", "Predicted phenotype", "%Identity", "%Overlap",
                       "Scaffold name", "Kraken2 Classified?", "Tax ID", "Lineage",
                       "Lowest identified rank", "Best 5 BLAST hits"]
            collection.write("\t".join(headers) + "\n")

            for key, entry in entries.items():
                # The first five fields are sample, gene, phenotype,
                # percent identity and percent overlap.
                base_info = entry[:5]
                # Construct the FASTA identifier used by Kraken2 and BLAST.
                # It combines the shortened sample name, gene and contig ID.
                unique_id = f"{entry[0]}_{entry[1]}_{entry[5]}"

                # Write base info
                collection.write("\t".join(base_info) + "\t")

                # Handle Kraken2 results
                # Append the Kraken2 record when present, or placeholder fields.
                # Review later: the length==1 branch does not match the
                # five-field records constructed by translate_kraken().
                if len(kraken_results.get(unique_id, [])) == 1:
                    collection.write(f"{unique_id}\tna\tna\tna\tna\t")
                else:
                    collection.write("\t".join(kraken_results.get(unique_id, ["na"] * 5)) + "\t")

                # Write BLAST hits
                # Append all BLAST hit descriptions found for this identifier.
                collection.write(", ".join(blast_results.get(unique_id, ["No hits"])))
                collection.write("\n")

    except Exception as e:
        print(f"Error processing results: {e}")
        sys.exit(1)


# Accept the five input/output paths provided by Snakemake.
def main():
    """Main entry point for the script."""
    if len(sys.argv) != 6:
        print(
            "Usage: python script.py <resistances> <kraken2_output> <kraken2_taxonomy> <blast_output> <collected_results>")
        sys.exit(1)

    collect_results(*sys.argv[1:])


# Run only when invoked as a script, not when imported as a module.
if __name__ == "__main__":
    main()
