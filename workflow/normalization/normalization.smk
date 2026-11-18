# Project-local CARD ontology and its installation/verification script.
ARO_ONTOLOGY = os.path.join(
    workflow.basedir,
    "workflow", "normalization", "data", "aro.obo"
)

ARO_CHECK_SCRIPT = os.path.join(
    workflow.basedir,
    "workflow", "normalization", "install_aro_ontology.py"
)

# Normalize StarAMR/ResFinder resistance gene annotations using the custom normalization script for downstream cross-tool comparison.
rule normalize_staramr_resfinder:
    input:
        # Original ResFinder resistance gene detections from StarAMR.
        resfinder = "{output_folder}/{tag}_staramr/resfinder.tsv",

        # Validated CARD ontology used for resistance annotation.
        ontology = ARO_ONTOLOGY
    output:
        # Normalized StarAMR/ResFinder results for each sample.
        "{output_folder}/argnorm/staramr/{tag}_staramr_argnorm_resfinder.tsv"
    params:
        # Locate the StarAMR normalization script.
        script = os.path.join(
            workflow.basedir,
            "workflow",
            "normalization",
            "normalize_staramr_resfinder.py"
        ),

        # Locate the ontology verification script.
        ontology_check = ARO_CHECK_SCRIPT

    shell:
        """
        # Verify that argNorm uses the validated project ontology.
        conda run --no-capture-output -n argnorm python "{params.ontology_check}" check

        # Normalize the StarAMR/ResFinder results.
        python "{params.script}" "{input.resfinder}" "{output}"
        """

# Normalize contig-based CARD/RGI resistance gene detections for downstream comparison with StarAMR and read-based RGI results.
rule normalize_rgi_contig:
    input:
        # Original CARD/RGI results generated from assembled scaffolds.
        rgi = "{output_folder}/card/contig/{tag}_RGI.txt",

        # Validated CARD ontology used for resistance annotation.
        ontology = ARO_ONTOLOGY
    output:
        # Normalized contig-based RGI results for each sample.
        "{output_folder}/argnorm/rgi/contig/{tag}_rgi_argnorm_contig.tsv"
    params:
        # Locate the custom RGI normalization script.
        script = os.path.join(
            workflow.basedir,
            "workflow",
            "normalization",
            "normalize_rgi.py"
        ),

        # Locate the ontology verification script.
        ontology_check = ARO_CHECK_SCRIPT

    shell:
        """
        # Verify that argNorm uses the validated project ontology.
        conda run --no-capture-output -n argnorm python "{params.ontology_check}" check

        # Normalize contig-based RGI results.
        conda run --no-capture-output -n argnorm python "{params.script}" "{input.rgi}" "{output}" contig
        """

# Normalize read-based CARD/RGI resistance gene detections for downstream comparison with StarAMR and contig-based RGI results.
rule normalize_rgi_reads:
    input:
        # Original RGI gene-mapping results generated from paired-end reads.
        rgi = "{output_folder}/card/reads/{tag}_RGI.gene_mapping_data.txt",

        # Validated CARD ontology used for resistance annotation.
        ontology = ARO_ONTOLOGY
    output:
        # Normalized read-based RGI results for each sample.
        "{output_folder}/argnorm/rgi/reads/{tag}_rgi_argnorm_reads.tsv"
    params:
        # Locate the custom RGI normalization script.
        script = os.path.join(
            workflow.basedir,
            "workflow",
            "normalization",
            "normalize_rgi.py"
        ),

        # Locate the ontology verification script.
        ontology_check = ARO_CHECK_SCRIPT

    shell:
        """
        # Verify that argNorm uses the validated project ontology.
        conda run --no-capture-output -n argnorm python "{params.ontology_check}" check

        # Normalize read-based RGI results.
        conda run --no-capture-output -n argnorm python "{params.script}" "{input.rgi}" "{output}" reads
        """
# Convert normalized StarAMR/ResFinder results into the shared table format used for cross-tool AMR detection comparisons.
rule common_format_staramr:
    input:
        # StarAMR/ResFinder results processed by the normalization script.
        "{output_folder}/argnorm/staramr/{tag}_staramr_argnorm_resfinder.tsv"
    output:
        # Standardized comparison table for each sample.
        "{output_folder}/comparison/staramr/{tag}.tsv"
    params:
        # Locate the script that converts normalized results into the common comparison format.
        script = os.path.join(
            workflow.basedir,
            "workflow",
            "normalization",
            "to_common_format.py"
        )
    shell:
        """
        # Convert the results and label their tool and analysis source.
        python "{params.script}" "{input}" "{output}" "StarAMR" "ResFinder"
        """

# Convert normalized contig-based CARD/RGI results into the common table format used for cross-tool AMR comparisons.
rule common_format_rgi_contig:
    input:
        # Normalized RGI results obtained from assembled scaffolds.
        "{output_folder}/argnorm/rgi/contig/{tag}_rgi_argnorm_contig.tsv"
    output:
        # Common-format comparison table for contig-based RGI detections.
        "{output_folder}/comparison/rgi/contig/{tag}.tsv"
    params:
        # Locate the script responsible for converting normalized results into the common comparison format.
        script = os.path.join(
            workflow.basedir,
            "workflow",
            "normalization",
            "to_common_format.py"
        )
    shell:
        """
        # Convert the results and label their source as RGI contig analysis.
        python "{params.script}" "{input}" "{output}" "RGI" "contig"
        """

# Convert normalized read-based CARD/RGI results into the common table format used for cross-tool AMR comparisons.
rule common_format_rgi_reads:
    input:
        # Normalized RGI results obtained from paired-end read mapping.
        "{output_folder}/argnorm/rgi/reads/{tag}_rgi_argnorm_reads.tsv"
    output:
        # Common-format comparison table for read-based RGI detections.
        "{output_folder}/comparison/rgi/reads/{tag}.tsv"
    params:
        # Locate the script responsible for converting normalized results into the common comparison format.
        script = os.path.join(
            workflow.basedir,
            "workflow",
            "normalization",
            "to_common_format.py"
        )
    shell:
        """
        # Convert the results and label their source as RGI read-based analysis.
        python "{params.script}" "{input}" "{output}" "RGI" "reads"
        """


# Generate per-sample Excel reports comparing AMR gene detections from StarAMR/ResFinder, RGI contig and RGI read-based analyses.
rule build_comparison_reports:
    input:
        # Common-format StarAMR/ResFinder results for all samples.
        staramr=expand(
            output_folder + "/comparison/staramr/{tag}.tsv",
            tag=tags
        ),

        # Common-format CARD/RGI contig results for all samples.
        rgi_contig=expand(
            output_folder + "/comparison/rgi/contig/{tag}.tsv",
            tag=tags
        ),

        # Common-format CARD/RGI read results for all samples.
        rgi_reads=expand(
            output_folder + "/comparison/rgi/reads/{tag}.tsv",
            tag=tags
        )
    output:
        # Per-sample workbook comparing detections across analysis methods.
        comparison=expand(
            output_folder + "/comparison/reports/{sample}/AMR_cross_tool_comparison.xlsx",
            sample=sample_names
        ),

        # Per-sample workbook containing detailed detection evidence.
        evidence=expand(
            output_folder + "/comparison/reports/{sample}/AMR_detection_evidence.xlsx",
            sample=sample_names
        )
    params:
        # Locate the custom Excel report-generation script.
        script=os.path.join(
            workflow.basedir,
            "workflow",
            "normalization",
            "build_comparison_reports.py"
        ),

        # Directory containing common-format input tables.
        input_dir=output_folder + "/comparison",

        # Directory where the final Excel reports are saved.
        output_dir=output_folder + "/comparison/reports"
    shell:
        """
        # Generate the per-sample comparison and evidence workbooks.
        python "{params.script}" "{params.input_dir}" "{params.output_dir}"
        """