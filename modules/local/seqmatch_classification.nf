process SEQMATCH_CLASSIFICATION {
    tag "$meta.id"+ "_" + "$cluster"

    conda "bioconda::rdptools=2.0.3"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/grep:3.4--hf43ccf4_4' :
        'biocontainers/grep:3.4--hf43ccf4_4' }"

    input:
    tuple val(meta), path(consensus), path(cluster_log), val(cluster)

    output:
    tuple val(meta), path('*_consensus_classification.csv'),      emit: classification
    tuple val(meta), path('*_classification.log'),                emit: log
    path "versions.yml",                                          emit: versions

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def cluster_id = "${cluster}"
    def seqmatch_db = params.seqmatch_db
    def seqmatch_accession = params.seqmatch_accession

    """
    echo "chosen classification: seqmatch"
    SequenceMatch seqmatch -k 5 ${seqmatch_db} $consensus | cut -f2,4 | sort | join -t \$'\t' -1 1 -2 1 -o 2.3,2.5,1.2 - ${seqmatch_accession} | sort -k3 -n -r -t '\t' | sed 's/\t/;/g' > ${prefix}_${cluster_id}_seqmatch_consensus_classification.csv
    cat $cluster_log > ${prefix}_${cluster_id}_classification.log
    echo -n ";" >> ${prefix}_${cluster_id}_classification.log
    SEQ_OUT=\$(head -n1 ${prefix}_${cluster_id}_seqmatch_consensus_classification.csv)
    echo \$SEQ_OUT >> ${prefix}_${cluster_id}_classification.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        seqmatch: 2.0.3
    END_VERSIONS
    """


}
