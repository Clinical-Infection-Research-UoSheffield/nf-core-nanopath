process GET_ABUNDANCE {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://mbdabrowska1/get-abundance:1.0' :
        'docker.io/mbdabrowska1/get-abundance:1.0' }"

    input:
    tuple val(meta), path(result_table)

    output:
    tuple val(meta), path('rel_abundance*.csv'),                  emit: results
    tuple val(meta), path('rel_abundance*_S.csv'),                emit: species_results
    tuple val(meta), path('*_chosen_classifier.csv'),            emit: chosen
    tuple val(meta), path('*_unresolved_taxids.txt'),            emit: unresolved, optional: true
    path "versions.yml",                                          emit: versions

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    get_abundance.py --infile $result_table --prefix ${prefix}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | cut -d ' ' -f2)
        pandas: \$(python -c "import pandas; print(pandas.__version__)")
    END_VERSIONS
    """
}