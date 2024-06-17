process GENERATE_REPORTS {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://mbdabrowska1/generate-reports:1.0' :
        'docker.io/mbdabrowska1/generate-reports:1.0' }"

    input:
    tuple val(meta), path(sample_result), path(fastp_results)
    val(positive_control)
    val(negative_control)
    tuple val(kit), val(run_id), val(seq_start)
    path(samplesheet)

    output:
    tuple val(meta), path("patient_report*.html"), emit: report
    path "versions.yml",                           emit: versions

    script:
    def revision=workflow.revision
    def clustering_size=params.umap_set_size
    def report_template="$baseDir/assets/UoS_report_template.html"
    def logo="$baseDir/assets/UoS_white_logo.txt"
    def negative="${negative_control}" ? "${negative_control}" : "[None]"
    def positive="${positive_control}" ? "${positive_control}" : "[None]"
    """
    READS_COUNT=\$(zcat ${fastp_results} | grep 'runid' | wc -l)
    echo ${negative}
    echo ${positive}
    echo ${meta.id}
    echo ${sample_result}
    echo \$READS_COUNT
    results_report.py \
        --infile ${sample_result} \
        --output patient_report \
        --barcode ${meta.id} \
        --info ${samplesheet} \
        --demux 'Guppy 6.4.6' \
        --clustering_size ${clustering_size} \
        --negative ${negative} \
        --positive ${positive} \
        --reads_count \$READS_COUNT \
        --kit "${kit}" \
        --report_template ${report_template} \
        --logo ${logo} \
        --run_id "${run_id}" \
        --seq_start "${seq_start}"
        
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        grep: \$(grep --version | head -n1 | cut -f 4 -d ' ')
        sed: \$(sed --version | head -n1 | cut -f 4 -d ' ')
        cut: \$(cut --version | head -n1 | cut -f 4 -d ' ')
        awk: \$(awk --version | head -n1 | cut -f 3 -d ' ')
    END_VERSIONS
    """


}