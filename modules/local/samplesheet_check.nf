process SAMPLESHEET_CHECK {
    tag "$samplesheet"
    label 'process_single'

    conda "conda-forge::python=3.8.3 pandas openpyxl pathlib"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://mbdabrowska1/samplesheet-check:1.0.2' :
        'docker.io/mbdabrowska1/samplesheet-check:1.0.2' }"

    input:
    path samplesheet
    val fastq_dir
    val clinical

    output:
    path '*.valid.csv'       , emit: csv
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script: 
    // Set fastqDir to fastq_dir path provided or leave as empty string
    def fastqDir = params.fastq_dir ? "--fastq_dir ${fastq_dir}" : ""
    def clinical = params.clinical ? "--clinical" : ""
    """
    python ${workflow.projectDir}/bin/check_samplesheet.py \\
        $samplesheet \\
        samplesheet.valid.csv \\
        $fastqDir \\
        $clinical \\

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
