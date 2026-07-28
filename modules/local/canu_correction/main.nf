process CANU_CORRECTION {
    tag "$meta.id"+ "_" + "$cluster"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/canu:2.2--ha47f30e_0':
        'biocontainers/canu:2.2--ha47f30e_0' }"

    input:
    tuple val(meta), path(reads), path(cluster_log), val(cluster)
    val mode
    val genomesize

    output:
    tuple val(meta), path("*.report")                                                         , emit: report
    tuple val(meta), path("*.correctedReads.fasta.gz"), path("*_racon.log"), val(cluster)     , emit: corrected_reads, optional: true
    path "versions.yml"                                                                       , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def valid_mode = ["-pacbio", "-nanopore", "-pacbio-hifi"]
    if ( !valid_mode.contains(mode) )  { error "Unrecognised mode to run Canu. Options: ${valid_mode.join(', ')}" }
    def count = params.polishing_reads
    def cluster_id = "${cluster}"
    """
    canu \\
        -correct \\
        -p ${prefix}_cluster${cluster_id} \\
        $mode \\
        genomeSize=$genomesize \\
        $args \\
        maxThreads=$task.cpus \\
        $reads \\

    gzip -d *.fasta.gz
    if grep "Found 0 reads" ${prefix}_cluster${cluster_id}.report
    then
        echo "Canu read correction has failed and the sample will be discontinued"
        exit 84
    fi

    READ_COUNT=\$(( \$(awk '{print \$1/2}' <(wc -l ${prefix}_cluster${cluster_id}.correctedReads.fasta)) ))
    cat $cluster_log > ${prefix}_cluster${cluster_id}_racon.log
    echo -n ";${count};\$READ_COUNT;" >> ${prefix}_cluster${cluster_id}_racon.log && cp ${prefix}_cluster${cluster_id}_racon.log ${prefix}_cluster${cluster_id}_racon_.log

    gzip *.fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        canu: \$(echo \$(canu --version 2>&1) | sed 's/^.*canu //; s/Using.*\$//' )
    END_VERSIONS
    """

}
