process CAT_FASTQS {
    tag "$fastq_dir"
    label 'process_single'

    conda "bioconda::grep=3.4"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/grep:3.4--hf43ccf4_4' :
        'biocontainers/grep:3.4--hf43ccf4_4' }"

    input:
    path fastq_dir
    val outdir

    output:
    path "*.fastq.gz",                           emit: fastq
    env NEW_DIR,                                 emit: new_fastq_dir
    path "versions.yml",                         emit: versions

    script:
    """
    for i in ${fastq_dir}/barcode*; do
        barcode=\$(basename \$i)
        output_file="\${barcode}.fastq.gz"

        if ls \$i/*.fastq.gz >/dev/null 2>&1; then
            cat \$i/*.fastq.gz > \$output_file
        elif ls \$i/*.fastq >/dev/null 2>&1; then
            cat \$i/*.fastq | gzip -c > \$output_file
        fi
    done
    echo "${outdir}"
    NEW_DIR="${outdir}/cat"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cat: \$(cat --version | head -n1 | cut -f 4 -d ' ')
        grep: \$(grep --version | head -n1 | cut -f 4 -d ' ')
        find: \$(find --version | head -n1 | cut -f 4 -d ' ')
        gzip: \$(gzip --version | head -n1 | cut -f 2 -d ' ')
    END_VERSIONS
    """


}
