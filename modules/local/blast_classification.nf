process BLAST_CLASSIFICATION {
    tag "$meta.id"+ "_" + "$cluster"

    conda "bioconda::blast"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://mbdabrowska1/full-classification:1.0' :
        'docker.io/mbdabrowska1/full-classification:1.0' }"

    input:
    tuple val(meta), path(consensus), path(cluster_log), val(cluster)

    output:
    tuple val(meta), path('*_consensus_classification.csv'),      emit: classification
    tuple val(meta), path('*_classification.log'),                emit: log
    path "versions.yml",                                          emit: versions

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def cluster_id = "${cluster}"
    def blast_db = params.blast_db

    """
    export BLASTDB=\$(dirname ${blast_db})
    blastn -query $consensus -db \$(basename ${blast_db}) -task megablast -dust no -outfmt '10 sscinames staxids evalue length pident bitscore' -evalue 11 -max_hsps 50 -max_target_seqs 100 | 
    sort -t ',' -k5nr -k6nr |
    awk 'NR<=5' |
    sed 's/,/;/g' > ${prefix}_${cluster_id}_blastn_consensus_classification.csv

    if [ -s ${prefix}_${cluster_id}_blastn_consensus_classification.csv ]; then
        echo "success"
    else
        echo "unclassified;0;0" >> ${prefix}_${cluster_id}_blastn_consensus_classification.csv
    fi
    BLAST_OUT=\$(cut -d";" -f1,2,5 ${prefix}_${cluster_id}_blastn_consensus_classification.csv | head -n1)
    
    cat $cluster_log > ${prefix}_${cluster_id}_classification.log
    echo -n ";" >> ${prefix}_${cluster_id}_classification.log
    echo \$BLAST_OUT >> ${prefix}_${cluster_id}_classification.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        blastn: \$(blastn -version | head -n1 | cut -d' ' -f2)
    END_VERSIONS
    """


}