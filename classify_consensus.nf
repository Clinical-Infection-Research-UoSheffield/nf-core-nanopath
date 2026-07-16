#!/usr/bin/env nextflow
/*
 * Feed a ready-made consensus FASTA into the pipeline from the classification step onward,
 * skipping read clustering / Canu / Racon / Medaka. It reuses the real pipeline modules, so
 * identification, the QC report, etc. are identical to a normal run.
 *
 * Run from the repo root, e.g. (use -profile singularity on an HPC, -profile docker on a GridION):
 *   nextflow run classify_consensus.nf -profile singularity \
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
 * --consensus accepts EITHER a single FASTA (one or more records) OR a quoted glob matching
 * several per-cluster FASTAs -- e.g. the separate consensus file Medaka emits for each cluster:
 *     --consensus 'medaka_out/*.fasta'      (quote it so Nextflow globs, not the shell)
 * Each record, across all matched files, becomes its own cluster (numbered 0..N-1) and all are
 * combined into a single sample report, exactly like a real run with several clusters. For
 * several INDEPENDENT samples (a separate report each), run the wrapper once per sample.
 */
nextflow.enable.dsl = 2

params.fake_reads = 100   // stand-in read count assigned to every cluster (-> even abundances)

include { FULL_CLASSIFICATION } from './modules/local/full_classification'
include { JOIN_RESULTS        } from './modules/local/join_results'
include { GET_ABUNDANCE       } from './modules/local/get_abundance/main'
include { GENERATE_REPORTS    } from './modules/local/generate_reports/main'

// Recreate the per-cluster "cluster log" that SPLIT_CLUSTERS normally makes: "<id>;<reads>".
// One instance per consensus record, so a multi-sequence FASTA becomes several clusters.
process STAGE_CONSENSUS {
    tag "${meta.id}_${cluster_id}"
    input:
    tuple val(meta), path(consensus), val(cluster_id)
    output:
    tuple val(meta), path(consensus), path("cluster.log"), val(cluster_id)
    script:
    """
    printf '%s;%s' "${cluster_id}" "${params.fake_reads ?: 100}" > cluster.log
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

    // one cluster per record: split the (possibly multi-sequence) FASTA and number them 0..N-1
    ch_clusters = Channel.fromPath( params.consensus )
        .splitFasta( file: true )
        .toList()
        .flatMap { recs -> recs.withIndex().collect { rec, i -> [ meta, rec, i ] } }

    STAGE_CONSENSUS( ch_clusters )

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
            file(params.input),
            Channel.value([                              // database names for the Run parameters table
                params.blast_db    ? file(params.blast_db).parent.name : 'n/a',
                params.kraken2_db  ? file(params.kraken2_db).name      : 'n/a',
                params.seqmatch_db ? file(params.seqmatch_db).name     : 'n/a',
                params.taxonomy    ? file(params.taxonomy).name        : 'n/a'
            ])
        )
    }
}
