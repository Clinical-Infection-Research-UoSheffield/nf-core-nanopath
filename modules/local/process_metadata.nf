process PROCESS_METADATA {
    tag "$report"
    label 'process_single'

    conda "bioconda::grep=3.4"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/grep:3.4--hf43ccf4_4' :
        'biocontainers/grep:3.4--hf43ccf4_4' }"

    input:
    tuple path(report), path(summary)

    output:
    tuple env(kit), env(run_id), env(seq_start), emit: metadata
    path "versions.yml",                         emit: versions

    script:
    """
    if grep -q '"Expansion kit", "value":' ${report}
    then
        kit=\$(grep '"Expansion kit", "value":' ${report} | grep -o -P 'Expansion.{0,35}' | cut -d '"' -f5)
    else
        kit=\$(grep '"Kit type", "value":' ${report} | grep -o -P 'Kit.{0,35}' | cut -d '"' -f5)

    fi
    run_id=\$(grep "protocol_run_id=" ${summary} | cut -d "=" -f2)
    seq_start=\$(grep 'started=' ${summary} | cut -d "=" -f2 | cut -d "." -f1 | sed 's/-/\\//g' | sed 's/T/ /g' | awk 'BEGIN{FS=OFS=" "} {split(\$1, a, /\\//); \$1 = a[3] "/" a[2] "/" a[1]} 1')
    echo \$kit
    echo \$run_id
    echo \$seq_start

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        grep: \$(grep --version | head -n1 | cut -f 4 -d ' ')
        sed: \$(sed --version | head -n1 | cut -f 4 -d ' ')
        cut: \$(cut --version | head -n1 | cut -f 4 -d ' ')
        awk: \$(awk --version | head -n1 | cut -f 3 -d ' ')
    END_VERSIONS
    """


}
