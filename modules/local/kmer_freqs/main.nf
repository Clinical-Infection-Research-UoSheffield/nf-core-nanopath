process KMER_FREQS {
    tag "$meta.id"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://mbdabrowska1/kmer_freqs:1.0' :
        'docker.io/mbdabrowska1/kmer_freqs:1.0' }"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*_freqs.txt"),        emit: freqs
    path "versions.yml",                         emit: versions

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def zipped = reads.toString().endsWith(".gz") ? '-z' : ''

    """
    echo ${zipped}
    kmer_freq.py -r $reads ${zipped} > ${prefix}_freqs.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | cut -d' ' -f2)
        biopython: \$(python -c 'import Bio; print(Bio.__version__)')
        tqdm: \$(python -c 'import tqdm; print(tqdm.__version__)')
    END_VERSIONS
    """
}
