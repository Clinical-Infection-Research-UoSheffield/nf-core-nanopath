process DRAFT_SELECTION {
    tag "$meta.id"+ "_" + "$cluster"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://quay.io/biocontainers/fastani:1.34--h4dfc31f_1' :
        'biocontainers/fastani:1.34--h4dfc31f_1' }"

    input:
    tuple val(meta), path(corrected_reads), path(cluster_log), val(cluster)

    output:
    tuple val(meta), path('*_draft_read.fasta'), path('*_draft.log'), path(corrected_reads), val(cluster),      emit: draft
    path "versions.yml",                                                                                        emit: versions

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def cluster_id = "${cluster}"
    def frag_length = "${meta.assay}" == "16S" ? 1200 : "${meta.assay}" == "ITS2" ? 400 : 160
    """
    split -l 2 <(zcat $corrected_reads) split_reads
    find split_reads* > read_list.txt

    fastANI --ql read_list.txt --rl read_list.txt -o fastani_output.ani -t $task.cpus -k 16 --fragLen ${frag_length}

    DRAFT=\$(awk 'NR>1{name[\$1] = \$1; arr[\$1] += \$3; count[\$1] += 1}  END{for (a in arr) {print arr[a] / count[a], name[a] }}' fastani_output.ani | sort -rg | cut -d " " -f2 | head -n1)
    cat \$DRAFT > ${prefix}_cluster${cluster_id}_draft_read.fasta

    if [ -s ${prefix}_cluster${cluster_id}_draft_read.fasta ]; then
        ID=\$(head -n1 ${prefix}_cluster${cluster_id}_draft_read.fasta | sed 's/>//g')
        cat $cluster_log > ${prefix}_${cluster_id}_draft.log
        echo -n \$ID >> ${prefix}_${cluster_id}_draft.log
    else
        exit 73
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fastANI: \$(fastANI | head -n1 | cut -f 2 -d ' ')
    END_VERSIONS
    """


}
