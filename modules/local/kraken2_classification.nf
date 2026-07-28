process KRAKEN2_CLASSIFICATION {
    tag "$meta.id"+ "_" + "$cluster"

    conda "bioconda::kraken2=2.1.3 rdptools=2.0.3"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/grep:3.4--hf43ccf4_4' :
        'docker.io/mbdabrowska1/kraken2-classification:1.0' }"

    input:
    tuple val(meta), path(consensus), path(cluster_log), val(cluster)

    output:
    tuple val(meta), path('*_consensus_classification.csv'),      emit: classification
    tuple val(meta), path('*_classification.log'),                emit: log
    tuple val(meta), path('*_classification_out.tsv'),            emit: tsv
    path "versions.yml",                                          emit: versions

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def cluster_id = "${cluster}"
    def kraken2_db = params.kraken2_db

    if(params.reclassifyOnFail){
        seqmatch_accession=params.seqmatch_accession
        seqmatch_db = params.seqmatch_db
        """
        echo "chosen classification: kraken2"
        echo "classifying with kraken2"
        kraken2 --db ${kraken2_db} --report ${prefix}_${cluster_id}_kraken2_consensus_classification.csv --output ${prefix}_${cluster_id}_kraken2_classification_out.tsv $consensus
        CLASS_LVL=\$(cut -f4 ${prefix}_${cluster_id}_kraken2_consensus_classification.csv | tail -n1)
        echo \$CLASS_LVL
        if [[ \$CLASS_LVL != "S"* ]]
        then
            echo "reclassifying"
            SequenceMatch seqmatch -k 5 ${seqmatch_db} $consensus | cut -f2,4 | sort | join -t \$'\t' -1 1 -2 1 -o 2.3,2.5,1.2 - ${seqmatch_accession} | sort -k3 -n -r -t '\t' | sed 's/\t/;/g' > ${prefix}_${cluster_id}_seqmatch_consensus_classification.csv
            cat $cluster_log > ${prefix}_${cluster_id}_classification.log
            echo -n ";" >> ${prefix}_${cluster_id}_classification.log
            SEQ_OUT=\$(head -n1 ${prefix}_${cluster_id}_seqmatch_consensus_classification.csv)
            echo \$SEQ_OUT >> ${prefix}_${cluster_id}_classification.log
            echo ${params.classification}
        else
            cat $cluster_log > ${prefix}_${cluster_id}_classification.log
            echo -n ";" >> ${prefix}_${cluster_id}_classification.log
            KR_OUT=\$(sed 's/\t/;/g' ${prefix}_${cluster_id}_kraken2_consensus_classification.csv | tr -s ' ' | sed 's/; /;/g' | cut -d ';' -f3,4,5,6 | grep -v '^0' | awk 'BEGIN {FS=";"; OFS=";"} {print \$4, \$3, \$2}')
            echo \$KR_OUT >> ${prefix}_${cluster_id}_classification.log
        fi

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            kraken2: \$(kraken2 --version | head -n1 | cut -d' ' -f3)
        END_VERSIONS
        """
    }
    else {
        """
        echo "chosen classification: kraken2"
        kraken2 --db ${kraken2_db} --report ${prefix}_${cluster_id}_kraken2_consensus_classification.csv --output ${prefix}_${cluster_id}_kraken2_classification_out.tsv $consensus
        KR_OUT=\$(sed 's/\t/;/g' ${prefix}_${cluster_id}_kraken2_consensus_classification.csv | tr -s ' ' | sed 's/; /;/g' | cut -d ';' -f3,4,5,6 | grep -v '^0' | awk 'BEGIN {FS=";"; OFS=";"} {print \$4, \$3, \$2}')
        cat $cluster_log > ${prefix}_${cluster_id}_classification.log
        echo -n ";" >> ${prefix}_${cluster_id}_classification.log
        echo \$KR_OUT >> ${prefix}_${cluster_id}_classification.log

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            kraken2: \$(kraken2 --version | head -n1 | cut -d' ' -f3)
        END_VERSIONS
        """
    }
}
