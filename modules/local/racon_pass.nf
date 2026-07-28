process RACON_PASS {
    tag "$meta.id"+ "_" + "$cluster"

    conda "bioconda::minimap2 racon"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://mbdabrowska1/racon-pass:1.0' :
        'docker.io/mbdabrowska1/racon-pass:1.0' }"

    input:
    tuple val(meta), path(draft_read), path(draft_log), path(corrected_reads), val(cluster)

    output:
    tuple val(meta), path('*_consensus.fasta'), path(draft_log), path(corrected_reads), val(cluster), env(SUCCESS),           emit: final_draft
    path "versions.yml",                                                                                                      emit: versions

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def cluster_id = "${cluster}"

    """
    SUCCESS=1

    # Perform read alignment
    minimap2 -ax map-ont --no-long-join -r100 -a "$draft_read" "$corrected_reads" -o "${prefix}_${cluster_id}_aligned.sam"

    # Check alignment success
    if ! racon --quality-threshold=9 -w 250 "$corrected_reads" "${prefix}_${cluster_id}_aligned.sam" "$draft_read" > "${prefix}_${cluster_id}_racon_consensus.fasta"; then
        SUCCESS=0
        cat "$draft_read" > "${prefix}_${cluster_id}_racon_consensus.fasta"
    fi

    # Check consensus file existence and size
    if [ ! -s "${prefix}_${cluster_id}_racon_consensus.fasta" ]; then
        SUCCESS=0
        cat "$draft_read" > "${prefix}_${cluster_id}_racon_consensus.fasta"
    fi


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        minimap2: \$(minimap2 --version | head -n1)
        racon: \$(racon --version | head -n1)
    END_VERSIONS
    """


}
