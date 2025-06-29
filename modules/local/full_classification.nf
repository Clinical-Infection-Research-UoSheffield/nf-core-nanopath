process FULL_CLASSIFICATION {
    tag "$meta.id"+ "_" + "$cluster"

    // If you intend to use Conda, ensure your environment.yml is correctly set up
    // and that the tools are available in the container or conda env.
    // If you prefer Docker-only (as implied by your current Dockerfile), you might remove conda.
    conda "bioconda::kraken2 rdptools blast"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://mbdabrowska1/full-classification:1.0' :
        'docker.io/mbdabrowska1/full-classification:1.0' }"

    input:
    tuple val(meta), path(consensus), path(cluster_log), val(cluster)
    path kraken2_db_dir
    path seqmatch_db_file
    path seqmatch_accession_file
    path blast_db_parent_dir
    val blast_db_prefix

    output:
    tuple val(meta), path('*_consensus_classification.csv'), emit: classification
    tuple val(meta), path('*_classification.log'),           emit: log
    tuple val(meta), path('*_classification_out.tsv'),       emit: tsv
    path "versions.yml",                                     emit: versions

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def cluster_id = "${cluster}"

    """
    echo "chosen classification: full"
    echo "classifying with kraken2"
    kraken2 --db ${kraken2_db_dir} --report ${prefix}_${cluster_id}_kraken2_consensus_classification.csv --output ${prefix}_${cluster_id}_kraken2_classification_out.tsv ${consensus}
    KR_OUT=\$(sed 's/\t/;/g' ${prefix}_${cluster_id}_kraken2_consensus_classification.csv | tr -s ' ' | sed 's/; /;/g' | cut -d ';' -f3,4,5,6 | grep -v '^0' | awk 'BEGIN {FS=";"; OFS=";"} {print \$4, \$3, \$2}')
    
    echo "classifying with seqmatch"
    SequenceMatch seqmatch -k 5 ${seqmatch_db_file} ${consensus} | cut -f2,4 | sort | join -t \$'\t' -1 1 -2 1 -o 2.3,2.5,1.2 - ${seqmatch_accession_file} | sort -k3 -n -r -t '\t' | sed 's/\t/;/g' > ${prefix}_${cluster_id}_seqmatch_consensus_classification.csv
    if [ -s ${prefix}_${cluster_id}_seqmatch_consensus_classification.csv ]; then
        echo "success"
    else
        echo "unclassified;0;0" >> ${prefix}_${cluster_id}_seqmatch_consensus_classification.csv
    fi
    SEQ_OUT=\$(head -n1 ${prefix}_${cluster_id}_seqmatch_consensus_classification.csv)

    echo "classifying with blastn"
    export BLASTDB=${blast_db_parent_dir}
    blastn -query ${consensus} -db ${blast_db_prefix} -task megablast -dust no -outfmt '10 sscinames staxids evalue length pident bitscore' -evalue 11 -max_hsps 50 -max_target_seqs 100 | sort -t ',' -k5nr -k6nr | head -n5 | sed 's/,/;/g' > ${prefix}_${cluster_id}_blastn_consensus_classification.csv
    if [ -s ${prefix}_${cluster_id}_blastn_consensus_classification.csv ]; then
        echo "success"
    else
        echo "unclassified;0;0" >> ${prefix}_${cluster_id}_blastn_consensus_classification.csv
    fi
    BLAST_OUT=\$(cut -d";" -f1,2,5 ${prefix}_${cluster_id}_blastn_consensus_classification.csv | head -n1)

    cat ${cluster_log} > ${prefix}_${cluster_id}_classification.log
    echo -n ";" >> ${prefix}_${cluster_id}_classification.log
    echo \$KR_OUT >> ${prefix}_${cluster_id}_classification.log
    echo \$SEQ_OUT >> ${prefix}_${cluster_id}_classification.log
    echo \$BLAST_OUT >> ${prefix}_${cluster_id}_classification.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        kraken2: \$(kraken2 --version | head -n1 | cut -d' ' -f3)
        blastn: \$(blastn -version | head -n1 | cut -d' ' -f2)
    END_VERSIONS
    """
}