process SUBSET_READS {
    tag "$meta.id"
    label 'process_single'

    conda "bioconda::grep=3.4"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/grep:3.4--hf43ccf4_4' :
        'biocontainers/grep:3.4--hf43ccf4_4' }"

    input:
    tuple val(meta), path(reads)
    val umap_set_size

    output:
    tuple val(meta), path("*_subset.fastq.gz"),  emit: subset
    path "versions.yml",                         emit: versions

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    zcat "$reads" | awk -v lines=\$((${umap_set_size}*4)) 'NR <= lines' > "${prefix}_subset.fastq"
    gzip "${prefix}_subset.fastq"


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cat: \$(cat --version | head -n1 | cut -f 4 -d ' ')
        gzip: \$(gzip --version | head -n1 | cut -f 2 -d ' ')
    END_VERSIONS
    """
}
