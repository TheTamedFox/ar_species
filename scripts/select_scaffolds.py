# Extract resistance-associated scaffolds from a metaSPAdes assembly using the resistance gene detections reported by StarAMR/ResFinder.
# The resulting FASTA is used for downstream Kraken2 and BLAST analysis.

import os
import subprocess
import csv
from sys import argv
from Bio import SeqIO


def collect_scaffolds(metaspades, resfinder, res_scaffolds):

    # Store resistance detections and the corresponding scaffold sequences.
    entries = {}
    scaffolds = {}
    i = 0
    first_line = True

    # Derive a shortened sample identifier from the assembly directory.
    name = metaspades.split("/")[-2].split("_")[0]

    # Read the StarAMR/ResFinder results line by line.
    with open(resfinder, "r") as res:
        for line in res:

            # Skip the TSV header.
            if first_line:
                first_line = False
            else:

                # Split the tab-separated detection fields.
                comp = line.split("\t")
                print(comp)

                # Exclude entries marked as "unknown" in the third field.
                if "unknown" not in comp[2]:

                    # Preserve selected ResFinder information for each detection.
                    # comp[7] identifies the scaffold in the tested ResFinder format.
                    entries[str(i)] = [
                        name, comp[1], comp[2], comp[3],
                        comp[4], comp[7], comp[10].rstrip()
                    ]

                    # Register the scaffold for sequence extraction.
                    scaffolds[comp[7]] = ""

                i += 1

    # Read the assembled scaffolds and collect matching sequences.
    for record in SeqIO.parse(metaspades, "fasta"):
        if record.id in scaffolds:
            scaffolds[record.id] = str(record.seq)

    # Write a FASTA record for each retained resistance detection.
    with open(res_scaffolds, "w") as fasta:
        for key in entries:

            if "unknown" not in entries[key][2]:

                # Construct a header containing sample, gene and scaffold ID.
                fasta_string = (
                    ">" + entries[key][0]
                    + "_" + entries[key][1]
                    + "_" + entries[key][5]
                )

                fasta.write(fasta_string + "\n")

                # Retrieve the scaffold sequence corresponding to this hit.
                fasta.write(scaffolds[entries[key][-2]] + "\n")


# Read the three file paths passed by the Snakemake select rule.
# Run scaffold extraction and report the input/output paths.
print("collect_scaffolds " + "\t".join([argv[1], argv[2], argv[3]]))
collect_scaffolds(argv[1], argv[2], argv[3])