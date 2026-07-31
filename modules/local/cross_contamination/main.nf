// PROTOTYPE (contamination-clustering branch): sequence-level cross-contamination check.
// Pools every barcode's consensus for the run and flags negative-control consensus sequences that
// are near-identical to a sample's consensus (shared DNA => contamination), independent of how the
// classifier named them. Feed it the combined FASTA built in the workflow (headers encode
// barcode|status|cluster|species). See docs/cross_contamination.md.
process CONSENSUS_CROSS_CHECK {
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://mbdabrowska1/read-clustering:1.0' :
        'docker.io/mbdabrowska1/read-clustering:1.0' }"

    input:
    path(combined_consensus)

    output:
    path("cross_contamination.tsv"), emit: tsv
    path "versions.yml",             emit: versions

    script:
    def min_identity = params.contamination_identity ?: 0.99
    """
    consensus_cross_contamination.py \\
        --fasta ${combined_consensus} \\
        --min-identity ${min_identity} \\
        --out cross_contamination.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | cut -d ' ' -f2)
        edlib: \$(python -c "import edlib; print(getattr(edlib, '__version__', 'NA'))")
    END_VERSIONS
    """
}
