#!/usr/bin/env nextflow
/*
========================================================================================
                         nf-core/resolvinf
========================================================================================
 nf-core/resolvinf Analysis Pipeline.
 #### Homepage / Documentation
 https://github.com/nf-core/resolvinf
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl = 2

/*
 * DEFAULT PARAMETERS
 */
params.help = false
params.name = false
params.hostnames = false
params.tracedir = "${params.outdir}/pipeline_info"

def helpMessage() {
    log.info"""

    Usage:

    The typical command for running the pipeline is as follows:

    nextflow run nf-core/resolvinf --input samplesheet.csv -profile docker

    Mandatory arguments:
      --input [file]                  Path to comma-separated file containing information about the samples
      -profile [str]                  Configuration profile to use. Can use multiple (comma separated)
                                      Available: conda, docker, singularity, test, awsbatch, <institute> and more

    ResolVI options:
      --annotation_label [str]        Column name to use for cell type annotation (default: 'cell_type')
      --marker_genes [file]           Path to file containing marker genes (one per line) or comma-separated list
      --max_epochs [int]              Maximum training epochs (default: 100)
      --num_samples [int]             Number of posterior samples (default: 20)
      --da_comparisons [file]         JSON or CSV file specifying differential abundance comparisons

    scVIVA options:
      --scviva_max_epochs [int]       Maximum training epochs for scVIVA (default: 100)
      --scviva_comparisons [file]     JSON or CSV file specifying scVIVA niche-aware DE comparisons

    Validation options:
      --validate_inputs               Enable input validation (default: true)

    Other options:
      --outdir [file]                 The output directory where the results will be saved (default: './results')
      --skip_multiqc                  Skip MultiQC report generation
      --email [email]                 Set this parameter to your e-mail address to get a summary e-mail
      --email_on_fail [email]         Same as --email, except only send mail if the workflow is not successful
      -name [str]                     Name for the pipeline run

    AWSBatch options:
      --awsqueue [str]                The AWSBatch JobQueue that needs to be set when running on AWSBatch
      --awsregion [str]               The AWS Region for your AWS Batch job to run on
      --awscli [str]                  Path to the AWS CLI tool

    GPU Configuration:
      This pipeline REQUIRES GPU for training processes.
      Use -profile gpu for GPU execution (default)
      Use -profile cpu to force CPU-only mode (slower, not recommended)
    """.stripIndent()
}

// Show help message
if (params.help) {
    helpMessage()
    exit 0
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { RESOLVI_PREPROCESS } from './modules/local/resolvi_preprocess'
include { RESOLVI_TRAIN } from './modules/local/resolvi_train'
include { RESOLVI_ANALYZE } from './modules/local/resolvi_analyze'
include { RESOLVI_VISUALIZE } from './modules/local/resolvi_visualize'
include { SCVIVA_TRAIN } from './modules/local/scviva_train'
include { SCVIVA_ANALYZE } from './modules/local/scviva_analyze'
include { MULTIQC } from './modules/nf-core/multiqc/main'
include { VALIDATE_INPUTS } from './modules/local/validate_inputs'

/*
 * SET UP CONFIGURATION VARIABLES
 */

// Has the run name been specified by the user?
custom_runName = params.name
if (!(workflow.runName ==~ /[a-z]+_[a-z]+/)) {
    custom_runName = workflow.runName
}

if (workflow.profile.contains('awsbatch')) {
    // AWSBatch sanity checking
    if (!params.awsqueue || !params.awsregion) exit 1, "Specify correct --awsqueue and --awsregion parameters on AWSBatch!"
    if (!params.outdir.startsWith('s3:')) exit 1, "Outdir not on S3 - specify S3 Bucket to run on AWSBatch!"
    if (params.tracedir.startsWith('s3:')) exit 1, "Specify a local tracedir or run without trace! S3 cannot be used for tracefiles."
}

/*
 * Create a channel for input spatialdata files
 */
if (!params.input) {
    exit 1, "Input samplesheet not specified!"
}

// Simple GPU mode detection
def gpu_mode = workflow.profile.contains('cpu') ? false : true
def gpu_display = gpu_mode ? 'GPU (required)' : 'CPU-only (forced)'

// Header log info
log.info """
=======================================================
                    nf-core/resolvinf
=======================================================
""".stripIndent()

def summary = [:]
if (workflow.revision) summary['Pipeline Release'] = workflow.revision
summary['Run Name']         = custom_runName ?: workflow.runName
summary['Input']            = params.input
summary['Annotation Label'] = params.annotation_label ?: 'cell_type'
summary['Marker Genes']     = params.marker_genes ? (params.marker_genes instanceof List ? params.marker_genes.size() + " genes" : params.marker_genes) : 'Default set'
summary['Max Epochs']       = params.max_epochs
summary['Num Samples']      = params.num_samples
summary['GPU Mode']         = gpu_display
summary['DA Comparisons']   = params.da_comparisons ?: 'Default pairwise'
summary['scVIVA Max Epochs'] = params.scviva_max_epochs
summary['scVIVA Comparisons'] = params.scviva_comparisons ?: 'Default pairwise'
summary['Validate Inputs']  = params.validate_inputs
summary['Skip MultiQC']     = params.skip_multiqc
summary['Max Resources']    = "$params.max_memory memory, $params.max_cpus cpus, $params.max_time time per job"
if (workflow.containerEngine) summary['Container'] = "$workflow.containerEngine - $workflow.container"
summary['Output dir']       = params.outdir
summary['Launch dir']       = workflow.launchDir
summary['Working dir']      = workflow.workDir
summary['Script dir']       = workflow.projectDir
summary['User']             = workflow.userName
if (workflow.profile.contains('awsbatch')) {
    summary['AWS Region']   = params.awsregion
    summary['AWS Queue']    = params.awsqueue
    summary['AWS CLI']      = params.awscli
}
summary['Config Profile'] = workflow.profile
if (params.config_profile_description) summary['Config Description'] = params.config_profile_description
if (params.config_profile_contact)     summary['Config Contact']     = params.config_profile_contact
if (params.config_profile_url)         summary['Config URL']         = params.config_profile_url
if (params.email || params.email_on_fail) {
    summary['E-mail Address']    = params.email
    summary['E-mail on failure'] = params.email_on_fail
    summary['MultiQC maxsize']   = params.max_multiqc_email_size
}
log.info summary.collect { k,v -> "${k.padRight(20)}: $v" }.join("\n")
log.info "-\033[2m--------------------------------------------------\033[0m-"

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {
    // Create input channel
    ch_input = Channel
        .fromPath(params.input)
        .splitCsv(header: true)
        .map { row ->
            // Validate required fields
            if (!row.sample_id) error "Missing sample_id in row: ${row}"
            if (!row.zarr_path) error "Missing zarr_path in row: ${row}"

            def meta = [
                id: row.sample_id,
                sample: row.sample_id,
                condition: row.condition ?: 'unknown',
                batch: row.batch ?: 'batch1',
                single_cell: true
            ]
            [meta, file(row.zarr_path, checkIfExists: true)]
        }

    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()

    // Input validation
    if (params.validate_inputs) {
        VALIDATE_INPUTS(ch_input)
        ch_validated_input = VALIDATE_INPUTS.out.validated
        ch_versions = ch_versions.mix(VALIDATE_INPUTS.out.versions.first())
        ch_multiqc_files = ch_multiqc_files.mix(VALIDATE_INPUTS.out.report.collect())
    } else {
        ch_validated_input = ch_input
    }

    // Preprocessing
    RESOLVI_PREPROCESS(
        ch_validated_input,
        params.annotation_label
    )
    ch_versions = ch_versions.mix(RESOLVI_PREPROCESS.out.versions.first())

    // Training - GPU mode is determined by profile
    RESOLVI_TRAIN(
        RESOLVI_PREPROCESS.out.adata,
        params.annotation_label,
        params.max_epochs,
        params.num_samples,
        gpu_mode
    )
    ch_versions = ch_versions.mix(RESOLVI_TRAIN.out.versions.first())

    // ResolVI Analysis
    ch_resolvi_analysis = RESOLVI_TRAIN.out.model
        .join(RESOLVI_TRAIN.out.adata, by: 0)
        .map { meta, model_dir, adata -> 
            if (!model_dir || !adata) {
                error "Missing required files for sample ${meta.id}"
            }
            [meta, model_dir, adata] 
        }

    RESOLVI_ANALYZE(
        ch_resolvi_analysis,
        params.da_comparisons ?: "default_pairwise"
    )
    ch_versions = ch_versions.mix(RESOLVI_ANALYZE.out.versions.first())

    // scVIVA Training
    SCVIVA_TRAIN(
        RESOLVI_TRAIN.out.model,
        RESOLVI_TRAIN.out.adata,
        params.scviva_max_epochs,
        gpu_mode
    )
    ch_versions = ch_versions.mix(SCVIVA_TRAIN.out.versions.first())

    // scVIVA Analysis
    ch_scviva_analysis = SCVIVA_TRAIN.out.model_dir
        .join(SCVIVA_TRAIN.out.adata_trained, by: 0)
    
    SCVIVA_ANALYZE(
        ch_scviva_analysis,
        params.scviva_comparisons ?: "default_pairwise"
    )
    ch_versions = ch_versions.mix(SCVIVA_ANALYZE.out.versions.first())

    // Visualization
    RESOLVI_VISUALIZE(
        RESOLVI_TRAIN.out.adata,
        params.marker_genes
    )
    ch_versions = ch_versions.mix(RESOLVI_VISUALIZE.out.versions.first())

    // Collect MultiQC inputs
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_versions.unique().collectFile(name: 'collated_versions.yml')
    )

    // MultiQC
    if (!params.skip_multiqc) {
        // Handle optional MultiQC config files
        ch_multiqc_config = params.multiqc_config ? 
            Channel.fromPath(params.multiqc_config, checkIfExists: true) : 
            Channel.empty()
    
        ch_multiqc_logo = params.multiqc_logo ? 
            Channel.fromPath(params.multiqc_logo, checkIfExists: true) : 
            Channel.empty()
    
        MULTIQC(
            ch_multiqc_files.collect(),
            ch_multiqc_config.collect(),
            ch_multiqc_logo.collect(),
            Channel.empty(),
            Channel.empty(),
            Channel.empty()
        )
        ch_versions = ch_versions.mix(MULTIQC.out.versions)
    }

}

/*
 * Completion e-mail notification
 */
workflow.onComplete {
    // Set up the e-mail variables
    def subject = "[nf-core/resolvinf] Successful: $workflow.runName"
    if (!workflow.success) {
        subject = "[nf-core/resolvinf] FAILED: $workflow.runName"
    }
    def email_fields = [:]
    email_fields['version'] = workflow.manifest.version
    email_fields['runName'] = custom_runName ?: workflow.runName
    email_fields['success'] = workflow.success
    email_fields['dateComplete'] = workflow.complete
    email_fields['duration'] = workflow.duration
    email_fields['exitStatus'] = workflow.exitStatus
    email_fields['errorMessage'] = (workflow.errorMessage ?: 'None')
    email_fields['errorReport'] = (workflow.errorReport ?: 'None')
    email_fields['commandLine'] = workflow.commandLine
    email_fields['projectDir'] = workflow.projectDir
    email_fields['summary'] = summary
    email_fields['summary']['Date Started'] = workflow.start
    email_fields['summary']['Date Completed'] = workflow.complete
    email_fields['summary']['Pipeline script file path'] = workflow.scriptFile
    email_fields['summary']['Pipeline script hash ID'] = workflow.scriptId
    if (workflow.repository) email_fields['summary']['Pipeline repository Git URL'] = workflow.repository
    if (workflow.commitId) email_fields['summary']['Pipeline repository Git Commit'] = workflow.commitId
    if (workflow.revision) email_fields['summary']['Pipeline Git branch/tag'] = workflow.revision
    if (workflow.container) email_fields['summary']['Docker image'] = workflow.container

    // On completion, print out the summary
    log.info "========================================="
    log.info "Pipeline completed at: $workflow.complete"
    log.info "Execution status: ${ workflow.success ? 'OK' : 'failed' }"
    log.info "========================================="
}
