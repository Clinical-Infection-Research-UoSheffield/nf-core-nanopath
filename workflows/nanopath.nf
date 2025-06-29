/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// Validate input parameters
WorkflowNanopath.initialise(params, log)

// TODO nf-core: Add all file path parameters for the pipeline to the list below
// Check input path parameters to see if they exist
def checkPathParamList = [ params.input, params.multiqc_config ]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Check mandatory parameters
if (params.input) { ch_input = file(params.input) } else { exit 1, 'Input samplesheet not specified!' }


if(params.onGridion){
        ch_metadata_files = Channel.fromPath(["$workflow.launchDir/../report*.html", "$workflow.launchDir/../final_summary*.txt"])
    }
else if(params.report_file && params.summary_file){
    ch_metadata_files = Channel.fromPath([params.report_file, params.summary_file])
}

//Create a boolean channel from params.clinical
ch_clinical = Channel.from([params.clinical])

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONFIG FILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

ch_multiqc_config          = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
ch_multiqc_custom_config   = params.multiqc_config ? Channel.fromPath( params.multiqc_config, checkIfExists: true ) : Channel.empty()
ch_multiqc_logo            = params.multiqc_logo   ? Channel.fromPath( params.multiqc_logo, checkIfExists: true ) : Channel.empty()
ch_multiqc_custom_methods_description = params.multiqc_methods_description ? file(params.multiqc_methods_description, checkIfExists: true) : file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
include { INPUT_CHECK } from '../subworkflows/local/input_check'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//
include { CAT_FASTQS                  } from '../modules/local/cat_fastqs'
include { PROCESS_METADATA            } from '../modules/local/process_metadata'
include { FASTP                       } from '../modules/nf-core/fastp/main'
include { FASTQC                      } from '../modules/nf-core/fastqc/main'
include { KRAKEN2_KRAKEN2             } from '../modules/nf-core/kraken2/kraken2/main'
include { SUBSET_READS                } from '../modules/local/subset_reads'
include { KMER_FREQS                  } from '../modules/local/kmer_freqs/main'
include { READ_CLUSTERING             } from '../modules/local/read_clustering/main'
include { SPLIT_CLUSTERS              } from '../modules/local/split_clusters'
include { CANU_CORRECTION             } from '../modules/local/canu_correction/main'
include { DRAFT_SELECTION             } from '../modules/local/draft_selection/main'
include { RACON_PASS                  } from '../modules/local/racon_pass'
include { MEDAKA_PASS                 } from '../modules/local/medaka_pass'
include { FULL_CLASSIFICATION         } from '../modules/local/full_classification'
include { BLAST_CLASSIFICATION        } from '../modules/local/blast_classification'
include { SEQMATCH_CLASSIFICATION     } from '../modules/local/seqmatch_classification'
include { KRAKEN2_CLASSIFICATION      } from '../modules/local/kraken2_classification'
include { JOIN_RESULTS                } from '../modules/local/join_results'
include { GET_ABUNDANCE               } from '../modules/local/get_abundance/main'
include { GENERATE_REPORTS            } from '../modules/local/generate_reports/main'
include { MULTIQC                     } from '../modules/nf-core/multiqc/main'
include { CUSTOM_DUMPSOFTWAREVERSIONS } from '../modules/nf-core/custom/dumpsoftwareversions/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Info required for completion email and summary
def multiqc_report = []

workflow NANOPATH {

    ch_versions = Channel.empty()

    if (params.fastq_dir && params.fastq_dir != false) {
        ch_fastq_dir = Channel.fromPath(params.fastq_dir)

        CAT_FASTQS (
            ch_fastq_dir,
            params.outdir,
        )

        ch_new_fastq_dir=CAT_FASTQS.out.new_fastq_dir
        ch_versions = ch_versions.mix(CAT_FASTQS.out.versions)
    } else {
        ch_new_fastq_dir = Channel.from([params.fastq_dir])
    }
    //
    // SUBWORKFLOW: Read in samplesheet, validate and stage input files
    //
    INPUT_CHECK (
        ch_input,
        ch_new_fastq_dir,
        ch_clinical
    )
    ch_versions = ch_versions.mix(INPUT_CHECK.out.versions)

    INPUT_CHECK.out.reads.branch {
        discontinued: it[0].status == "discontinued"
        samples: true
    }.set { ch_files }

    // if channel ch_metadata_files exists, then runn process_metadata
    if(params.onGridion || (params.report_file && params.summary_file)){
        PROCESS_METADATA (
            ch_metadata_files.collect()
        )
        ch_meta_final = PROCESS_METADATA.out.metadata
    } else {
        ch_meta_final = Channel.of([params.kit, params.run_id, params.seq_start])
    }

    FASTP (
        ch_files.samples,
        [],
        false,
        false
    )
    ch_versions = ch_versions.mix(FASTP.out.versions.first())

    // Filter out samples that failed FASTP (aka. no reads left due to fltering or quality control)
    FASTP.out.reads
        .join(FASTP.out.json)
        .map {
            meta, reads, json ->
                if (WorkflowNanopath.getFastpReadsAfterFiltering(json) > 0) {
                    [meta, reads]
                } else {
                    // Add a warning message to the console output
                    log.warn "${meta.id} has been discontinued due to a failed FASTP quality check."
                    null
                }
        }.set { ch_fastp_pass }


    FASTQC (
        INPUT_CHECK.out.reads
    )
    ch_versions = ch_versions.mix(FASTQC.out.versions.first())

    if(params.remove_unclassified){
        // Remove unclassified reads
        KRAKEN2_KRAKEN2 (
            ch_fastp_pass,
            params.kraken2_db,
            true,
            false
        )
        ch_for_subsetting = KRAKEN2_KRAKEN2.out.classified_reads_fastq
        ch_versions = ch_versions.mix(KRAKEN2_KRAKEN2.out.versions.first())
    } else {
        ch_for_subsetting = ch_fastp_pass
    }

    SUBSET_READS (
        ch_for_subsetting,
        params.umap_set_size
    )
    ch_versions = ch_versions.mix(SUBSET_READS.out.versions.first())

    KMER_FREQS (
        SUBSET_READS.out.subset
    )
    ch_versions = ch_versions.mix(KMER_FREQS.out.versions.first())

    READ_CLUSTERING (
        KMER_FREQS.out.freqs,
    )
    ch_versions = ch_versions.mix(READ_CLUSTERING.out.versions.first())

    SUBSET_READS.out.subset
        .join(READ_CLUSTERING.out.clusters, by: [0])
        .set{ ch_splitting }

    SPLIT_CLUSTERS (
        ch_splitting
    )
    ch_versions = ch_versions.mix(SPLIT_CLUSTERS.out.versions.first())

    SPLIT_CLUSTERS.out.reads.map {
        meta, reads, log, cluster_id ->
            clusters = cluster_id.split(" ").collect { it.toInteger() }
            [meta, reads, log, clusters]
    }.transpose().set {ch_split_cluster}

    CANU_CORRECTION (
        ch_split_cluster,
        "-nanopore",
        params.avg_amplicon_size
    )

    ch_versions = ch_versions.mix(CANU_CORRECTION.out.versions.first())

    DRAFT_SELECTION (
        CANU_CORRECTION.out.corrected_reads
    )

    ch_versions = ch_versions.mix(DRAFT_SELECTION.out.versions.first())

    RACON_PASS (
        DRAFT_SELECTION.out.draft
    )

    ch_versions = ch_versions.mix(RACON_PASS.out.versions.first())

    // Check for racon success and print warning if success is 0
    RACON_PASS.out.final_draft.map {
        meta, draft, log, corrected_reads, cluster_id, success ->
            if(success == "0"){
                log.warn "Sample ${meta.id} : Racon correction for cluster ${cluster_id} failed due to not enough overlaps. Taking draft read as consensus"
            }
    }

    MEDAKA_PASS (
        RACON_PASS.out.final_draft
    )

    ch_versions = ch_versions.mix(MEDAKA_PASS.out.versions.first())

    def blast_db_parent_path = file(params.blast_db).parent
    def blast_db_name = file(params.blast_db).getBaseName()

    // run FULL_CLASSIFICATION if params.classification is "full"
    if(params.classification == "full"){
        FULL_CLASSIFICATION (
            MEDAKA_PASS.out.consensus,
            params.kraken2_db,
            params.seqmatch_db,
            params.seqmatch_accession,
            blast_db_parent_path,
            blast_db_name
        )
        ch_join_results = FULL_CLASSIFICATION.out.log.groupTuple()
    } else if(params.classification == "blast"){
        BLAST_CLASSIFICATION (
            MEDAKA_PASS.out.consensus,
            params.blast_db
        )
        ch_join_results = BLAST_CLASSIFICATION.out.log.groupTuple()
    } else if(params.classification == "seqmatch"){
        SEQMATCH_CLASSIFICATION (
            MEDAKA_PASS.out.consensus,
            params.seqmatch_db,
            params.seqmatch_accession
        )
        ch_join_results = SEQMATCH_CLASSIFICATION.out.log.groupTuple()
    } else {
        KRAKEN2_CLASSIFICATION (
            MEDAKA_PASS.out.consensus,
            params.kraken2_db
        )
        ch_join_results = KRAKEN2_CLASSIFICATION.out.log.groupTuple()
    }

    JOIN_RESULTS (
        ch_join_results
    )

    GET_ABUNDANCE (
        JOIN_RESULTS.out.classification
    )
    
    if(!params.onGridion){
        Channel.from(params.kit, 'unknown').set{ch_barcoding_kit}
    }

    if(params.clinical && params.generateReports){

        GET_ABUNDANCE.out.species_results.branch{
            negative: it[0].status == "negative control"
                return it[1]
            positive: it[0].status == "positive control"
                return it[1]
        }.set { ch_controls }

        ch_reporting = GET_ABUNDANCE.out.species_results.join(FASTP.out.reads, by: [0])

        GENERATE_REPORTS(
            ch_reporting,
            ch_controls.positive.toList(),
            ch_controls.negative.toList(),
            ch_meta_final,
            ch_input
        )
    }
    

    
    

    CUSTOM_DUMPSOFTWAREVERSIONS (
        ch_versions.unique().collectFile(name: 'collated_versions.yml')
    )

    //
    // MODULE: MultiQC
    //
    workflow_summary    = WorkflowNanopath.paramsSummaryMultiqc(workflow, summary_params)
    ch_workflow_summary = Channel.value(workflow_summary)

    methods_description    = WorkflowNanopath.methodsDescriptionText(workflow, ch_multiqc_custom_methods_description)
    ch_methods_description = Channel.value(methods_description)

    ch_multiqc_files = Channel.empty()
    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect())
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]}.ifEmpty([]))

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList()
    )
    multiqc_report = MULTIQC.out.report.toList()
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COMPLETION EMAIL AND SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow.onComplete {
    if (params.email || params.email_on_fail) {
        NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report)
    }
    NfcoreTemplate.summary(workflow, params, log)
    if (params.hook_url) {
        NfcoreTemplate.IM_notification(workflow, params, summary_params, projectDir, log)
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
