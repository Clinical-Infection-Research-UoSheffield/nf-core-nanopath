// Run-level QC report: collects every barcode's qc_<barcode>.json (emitted by GENERATE_REPORTS)
// into a single run sheet so the whole run can be quality-checked from the pipeline's point of
// view -- per barcode reads, clusters formed, clusters that failed to assemble (and their size),
// Racon-unpolished clusters, and the dropped-read fraction, with an expandable per-cluster table.
// Barcodes present in the samplesheet roster but missing a qc json (0 viable reads -> no report)
// are shown as failure rows. Pure-stdlib script, so it reuses the reports container.
process RUN_QC_REPORT {
    tag "run qc"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://mbdabrowska1/generate-reports:1.0' :
        'docker.io/mbdabrowska1/generate-reports:1.0' }"

    input:
    path qc_jsons          // collected qc_<barcode>.json (may be empty if nothing reported)
    path roster            // TSV 'barcode<TAB>status' for every samplesheet barcode
    tuple val(kit), val(run_id), val(seq_start)

    output:
    path "run_qc_report.html", emit: report
    path "run_qc_report.tsv",  emit: tsv
    path "versions.yml",       emit: versions

    script:
    def version = workflow.manifest.version ?: 'unknown'
    // Mirror GENERATE_REPORTS: fall back to the Nextflow run name when the run report wasn't found
    // (e.g. offline instrument) so the sheet still carries a run identifier.
    def run_identifier = (run_id && "${run_id}" != 'unknown') ? run_id : workflow.runName
    """
    run_qc_report.py \\
        --qc ${qc_jsons} \\
        --roster ${roster} \\
        --run_id "${run_identifier}" \\
        --kit "${kit}" \\
        --seq_start "${seq_start}" \\
        --version "${version}" \\
        --out run_qc_report

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version 2>&1 | cut -f2 -d ' ')
    END_VERSIONS
    """
}
