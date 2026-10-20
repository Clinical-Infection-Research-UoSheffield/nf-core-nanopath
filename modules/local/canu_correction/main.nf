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
    // After a few attempts, cap a very high-coverage cluster's reads. Canu is a genome assembler and
    // its correction step can fail stochastically when handed tens of thousands of reads of a single
    // short 16S amplicon (extreme coverage). Retries alone usually recover it (see errorStrategy);
    // if they don't, capping the input to a still-ample read count avoids the pathological coverage.
    def cap_reads = task.attempt > 3 ? 20000 : 0
    """
    INPUT_READS=${reads}
    CAP=${cap_reads}
    if [ "\$CAP" -gt 0 ] && head -c1 ${reads} | grep -q '>'; then
        NREC=\$(grep -c '^>' ${reads} || echo 0)
        if [ "\$NREC" -gt "\$CAP" ]; then
            awk -v cap="\$CAP" '/^>/{n++} n>cap{exit} {print}' ${reads} > canu_capped_input.fasta
            INPUT_READS=canu_capped_input.fasta
            echo "Attempt ${task.attempt}: capped \$NREC reads to \$CAP for canu correction"
        fi
    fi

    canu \\
        -correct \\
        -p ${prefix}_cluster${cluster_id} \\
        $mode \\
        genomeSize=$genomesize \\
        $args \\
        maxThreads=$task.cpus \\
        \$INPUT_READS \\

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
