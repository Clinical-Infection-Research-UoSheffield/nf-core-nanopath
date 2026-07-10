#!/usr/bin/env nextflow
/*
 * Feed a ready-made consensus FASTA into the pipeline from the classification step onward,
 * skipping read clustering / Canu / Racon / Medaka. It reuses the real pipeline modules, so
 * identification, the QC report, etc. are identical to a normal run.
 *
 * Run from the repo root, e.g.:
 *   nextflow run classify_consensus.nf -profile docker \
 *       --classification full \            # REQUIRED: otherwise defaults to kraken2-only
 *       --consensus my_consensus.fasta \
 *       --barcode barcode01 \
 *       --kraken2_db /path/kraken2_db \
 *       --blast_db /path/16S_ribosomal_RNA \
 *       --seqmatch_db /path/seqmatch_trained \
 *       --seqmatch_accession /path/accession.txt \
 *       --taxonomy /path/rankedlineage.dmp \
 *       --outdir consensus_out \
 *       [--input experiment_info.xlsx]     # add this to also build the patient report
 *
 * The FASTA should contain a SINGLE consensus sequence (one cluster). To classify several,
 * run it once per FASTA (or extend the input channel below).
 */
nextflow.enable.dsl = 2

include { FULL_CLASSIFICATION } from './modules/local/full_classification'
include { JOIN_RESULTS        } from './modules/local/join_results'
include { GET_ABUNDANCE       } from './modules/local/get_abundance/main'
include { GENERATE_REPORTS    } from './modules/local/generate_reports/main'

// Recreate the per-cluster "cluster log" that SPLIT_CLUSTERS normally makes: "<id>;<reads>"
process STAGE_CONSENSUS {
    tag "$meta.id"
    input:
    tuple val(meta), path(consensus)
    output:
    tuple val(meta), path(consensus), path("cluster.log"), val(0)
    script:
    """
    printf '0;%s' "${params.fake_reads ?: 100}" > cluster.log
    """
}

// GENERATE_REPORTS counts reads by grepping 'runid' in a fastq.gz; supply a stand-in
process DUMMY_FASTQ {
    output:
    path "dummy.fastq.gz"
    script:
    """
    printf '@r runid=0\\nACGT\\n+\\n!!!!\\n' | gzip > dummy.fastq.gz
    """
}

workflow {
    if( !params.consensus ) { error "Provide --consensus <file.fasta>" }

    def meta = [ id            : params.barcode ?: 'barcode01',
                 single_end    : true,
                 assay         : params.assay  ?: '16S',
                 status        : params.status ?: 'sample',
                 specimen_number: params.specimen ?: (params.barcode ?: 'CONSENSUS') ]

    STAGE_CONSENSUS( Channel.of( [ meta, file(params.consensus) ] ) )

    def blast_parent = file(params.blast_db).parent
    def blast_name   = file(params.blast_db).getBaseName()

    // ---- identification (same modules as the full pipeline) ----
    FULL_CLASSIFICATION( STAGE_CONSENSUS.out,
                         params.kraken2_db, params.seqmatch_db, params.seqmatch_accession,
                         blast_parent, blast_name )

    JOIN_RESULTS( FULL_CLASSIFICATION.out.log.groupTuple(), params.taxonomy )
    GET_ABUNDANCE( JOIN_RESULTS.out.classification )

    // ---- optional: build the patient QC report (needs a samplesheet with this barcode) ----
    if( params.input ) {
        ch_hit = FULL_CLASSIFICATION.out.classification.groupTuple().map { m, f -> [ m, f.flatten() ] }
        DUMMY_FASTQ()

        ch_report = GET_ABUNDANCE.out.species_results
            .combine( DUMMY_FASTQ.out )                 // [meta, rel_abundance_S.csv, dummy.fastq.gz]
            .join( ch_hit, by: 0 )                       // + hit-detail csvs
            .join( GET_ABUNDANCE.out.chosen, by: 0 )     // + chosen_classifier.csv

        GENERATE_REPORTS(
            ch_report,
            Channel.value('[None]'),                     // positive control (none)
            Channel.value('[None]'),                     // negative control (none)
            Channel.value([ params.kit ?: 'unknown', params.run_id ?: 'unknown', params.seq_start ?: 'unknown' ]),
            file(params.input)
        )
    }
}
