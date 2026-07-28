process MEDAKA_PASS {
    tag "$meta.id"+ "_" + "$cluster"

    conda "bioconda::medaka=1.4.4"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://quay.io/biocontainers/medaka:1.4.4--py38h130def0_0' :
        'biocontainers/medaka:1.4.4--py38h130def0_0' }"

    input:
    tuple val(meta), path(draft), path(draft_log), path(corrected_reads), val(cluster), val(success)

    output:
    tuple val(meta), path('*_consensus_medaka/consensus.fasta'), path(draft_log), val(cluster),          emit: consensus
    path "versions.yml",                                                                       emit: versions

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def cluster_id = "${cluster}"
    def success = "${success}"

    """
    if medaka_consensus -i $corrected_reads -d $draft -o ${prefix}_${cluster_id}_consensus_medaka -t 4 ; then
        echo "Command succeeded"
    else
        echo "Command failed. Copying draft reads to output."
        mkdir ${prefix}_${cluster_id}_consensus_medaka
        cat $draft > ${prefix}_${cluster_id}_consensus_medaka/consensus.fasta
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        medaka: \$(medaka --version)
    END_VERSIONS
    """


}
