process JOIN_RESULTS {
    tag "$meta.id"
    label 'process_single'

    conda "bioconda::grep=3.4"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://ubuntu' :
        'docker.io/ubuntu' }"

    input:
    tuple val(meta), path(logs)
    path taxonomy

    output:
    tuple val(meta), path('*.nanoclust_out.txt'),                 emit: classification
    path "versions.yml",                                          emit: versions

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"

    if(params.classification=='full'){
        """
        echo "chosen classification: full"
        echo "id;reads_in_cluster;used_for_consensus;reads_after_corr;draft_id;kraken2_sciname;taxid;class_level;name;species;genus;family;order;seqmatch_sciname;taxid;class_level;name;species;genus;family;order;blast_sciname;taxid;class_level;name;species;genus;family;order;" > ${prefix}.nanoclust_out.txt
        for i in $logs; do
            while read line; do
                echo \$line
                TAXID=\$(echo \$line | awk -F';' '{print \$(NF-1)}')
                echo \$TAXID
                TAXinDB=\$(grep -w "^\${TAXID}" $taxonomy || [[ \$? == 1 ]])
                echo \$TAXinDB
                echo -n \$(echo \$line | tr -d '\n') >> ${prefix}.nanoclust_out.txt
                if [ "\$TAXID" != "0" ] | [ "\$TAXID" != "" ] | [ "\$TAXinDB" != "" ]; then
                    echo -n ";" >> ${prefix}.nanoclust_out.txt
                    TAXONOMY=\$(grep -w "^\${TAXID}" $taxonomy | tr -d '\t' | cut -d '|' -f2,3,4,5,6 --output-delimiter ';')
                    echo -n "\$TAXONOMY;" >> ${prefix}.nanoclust_out.txt
                else
                    echo -n ";;;;;;" >> ${prefix}.nanoclust_out.txt
                fi
            done <\$i
            echo -e "\n" >> ${prefix}.nanoclust_out.txt
        done
        sed -i 's/.\$//' ${prefix}.nanoclust_out.txt

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            cat: \$(cat --version | head -n1 | cut -f 4 -d ' ')
            grep: \$(grep --version | head -n1 | cut -f 4 -d ' ')
        END_VERSIONS
        """
    }
    else if(params.classification=='blast'){
        """
        echo "chosen classification: blast"
        echo "id;reads_in_cluster;used_for_consensus;reads_after_corr;draft_id;sciname;taxid;length;per_ident" > ${prefix}.nanoclust_out.txt

        for i in $logs; do
            cat \$i >> ${prefix}.nanoclust_out.txt
        done

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            cat: \$(cat --version | head -n1 | cut -f 4 -d ' ')
            grep: \$(grep --version | head -n1 | cut -f 4 -d ' ')
        END_VERSIONS
        """
    }
    else if(params.classification=='seqmatch'){
        """
        echo "chosen classification: seqmatch"
        echo "id;reads_in_cluster;used_for_consensus;reads_after_corr;draft_id;sciname;taxid;seqmatch_score;name;species;genus;family;order" > ${prefix}.nanoclust_out.txt

        for i in $logs; do
            TAXID=\$(cut -d ";" -f7 \$i)
            TAXinDB=\$(grep -w "^\${TAXID}" $taxonomy || [[ \$? == 1 ]])
            cat \$i | tr -d '\n' >> ${prefix}.nanoclust_out.txt
            if [ "\$TAXID" != "0" ] | [ "\$TAXID" != "" ] | [ "\$TAXinDB" != "" ]; then
                echo -n ";" >> ${prefix}.nanoclust_out.txt
                TAXONOMY=\$(grep -w "^\${TAXID}" $taxonomy | tr -d '\t' | cut -d '|' -f2,3,4,5,6 --output-delimiter ';')
                echo "\$TAXONOMY" >> ${prefix}.nanoclust_out.txt
            else
                echo ";;;;" >> ${prefix}.nanoclust_out.txt
            fi
        done

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            cat: \$(cat --version | head -n1 | cut -f 4 -d ' ')
            grep: \$(grep --version | head -n1 | cut -f 4 -d ' ')
        END_VERSIONS
        """
    }
    else if(params.classification=='kraken2'){
        """
        echo "chosen classification: kraken2"
        echo "id;reads_in_cluster;used_for_consensus;reads_after_corr;draft_id;sciname;taxid;class_level;name;species;genus;family;order" > ${prefix}.nanoclust_out.txt

        for i in $logs; do
            TAXID=\$(cut -d ";" -f7 \$i)
            TAXinDB=\$(grep -w "^\${TAXID}" $taxonomy || [[ \$? == 1 ]])
            cat \$i | tr -d '\n' >> ${prefix}.nanoclust_out.txt
            if [ "\$TAXID" != "0" ] | [ "\$TAXID" != "" ] | [ "\$TAXinDB" != "" ]; then
                echo -n ";" >> ${prefix}.nanoclust_out.txt
                TAXONOMY=\$(grep -w "^\${TAXID}" $taxonomy | tr -d '\t' | cut -d '|' -f2,3,4,5,6 --output-delimiter ';')
                echo "\$TAXONOMY" >> ${prefix}.nanoclust_out.txt
            else
                echo ";;;;" >> ${prefix}.nanoclust_out.txt
            fi
        done

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            cat: \$(cat --version | head -n1 | cut -f 4 -d ' ')
            grep: \$(grep --version | head -n1 | cut -f 4 -d ' ')
        END_VERSIONS
        """
    }
}