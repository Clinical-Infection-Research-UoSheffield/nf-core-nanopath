process SPLIT_CLUSTERS {
    tag "$meta.id"
    label 'process_single'

    conda "bioconda::seqtk"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://quay.io/biocontainers/seqtk:1.4--he4a0461_1' :
        'biocontainers/seqtk:1.4--he4a0461_1' }"

    input:
    tuple val(meta), path(reads), path(clusters)

    output:
    tuple val(meta), path('*[0-9]*.fastq'), path('*[0-9]*.log'), env(CLUSTERS),      emit: reads
    path "versions.yml",                                                             emit: versions

    script:
    """
    sed 's/\\srunid.*//g' "$reads" > only_id_header_readfile.fastq

    CLUSTERS=\$(awk '(\$5 ~ /[0-9]/) {print \$5}' "$clusters" | sort | uniq)
    echo "CLUSTERS: \$CLUSTERS"

    for id in \$CLUSTERS; do
        echo \$id
        awk -v cluster="\$id" '(\$5 == cluster) {print \$1}' "$clusters" > cluster_\$id.txt
        seqtk subseq only_id_header_readfile.fastq -A cluster_\$id.txt > \$id.fastq
        READ_COUNT=\$(( \$(awk '{print \$1/4}' <(wc -l \$id.fastq)) ))
        echo -n "\$id;\$READ_COUNT" > \$id.log
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        seqtk: \$(seqtk | head -n1 | cut -f 2 -d ' ')
    END_VERSIONS
    """


}
