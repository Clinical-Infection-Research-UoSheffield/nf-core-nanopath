// Opt-in (params.wait_for_run_report): let the pipeline be launched while a sequencing run is still
// going, and start processing only once the run report is written (i.e. the run has finished and the
// fastq_pass is complete). It polls the run folder for the ONT report ('report*.html' /
// 'final_summary*.txt') and exits as soon as one appears; the rest of the workflow is gated on it.
// Default behaviour (flag off) is unchanged, so already-completed runs -- whose report may be present
// or may have been deleted -- run immediately as before.
process AWAIT_RUN_REPORT {
    tag "await run report"
    label 'process_single'
    executor 'local'   // a long-lived poller belongs on the head node, not a submitted/queued job

    conda "bioconda::grep=3.4"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/grep:3.4--hf43ccf4_4' :
        'biocontainers/grep:3.4--hf43ccf4_4' }"

    input:
    path run_dir   // the run folder (parent of fastq_pass), staged live so new files are visible

    output:
    path "run_ready.flag", emit: ready

    script:
    def poll    = params.run_report_poll_interval ?: 60   // seconds between checks
    def timeout = params.run_report_timeout ?: 0          // seconds; 0 = wait indefinitely
    """
    echo "Waiting for the run report (report*.html / final_summary*.txt) in the run folder before starting; polling every ${poll}s..."
    waited=0
    until ls ${run_dir}/report*.html >/dev/null 2>&1 || ls ${run_dir}/final_summary*.txt >/dev/null 2>&1; do
        sleep ${poll}
        waited=\$(( waited + ${poll} ))
        if [ ${timeout} -gt 0 ] && [ \$waited -ge ${timeout} ]; then
            echo "ERROR: timed out after \${waited}s waiting for the run report in ${run_dir}" >&2
            exit 1
        fi
    done
    echo "Run report found after \${waited}s -- starting the pipeline."
    touch run_ready.flag
    """
}
