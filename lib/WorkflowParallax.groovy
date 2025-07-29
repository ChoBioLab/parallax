import nextflow.Nextflow
import groovy.json.JsonSlurper

class WorkflowParallax {

    //
    // Check and validate parameters
    //
    public static void initialise(params, log) {
        genomeExistsError(params, log)

        if (!params.input) {
            Nextflow.error "Please provide an input samplesheet to the pipeline e.g. '--input samplesheet.csv'"
        }

        // Validate GPU configuration
        validateGpuConfig(params, log)

        // Validate file inputs
        validateFileInputs(params, log)

        // Validate comparison files
        validateComparisonFiles(params, log)

        // Validate marker genes
        validateMarkerGenes(params, log)
    }

    //
    // Validate GPU configuration
    //
    public static void validateGpuConfig(params, log) {
        if (params.num_gpus != null) {
            if (params.num_gpus < -1) {
                Nextflow.error "Invalid num_gpus value: ${params.num_gpus}. Must be >= -1"
            }
            if (params.num_gpus > 8) {
                log.warn "num_gpus value ${params.num_gpus} is very high. Consider using fewer GPUs."
            }
        }

        if (params.gpu_memory_fraction != null) {
            if (params.gpu_memory_fraction <= 0 || params.gpu_memory_fraction > 1) {
                Nextflow.error "Invalid gpu_memory_fraction: ${params.gpu_memory_fraction}. Must be between 0 and 1"
            }
        }

        if (params.max_gpus_per_process != null) {
            if (params.max_gpus_per_process < 1) {
                Nextflow.error "Invalid max_gpus_per_process: ${params.max_gpus_per_process}. Must be >= 1"
            }
        }
    }

    //
    // Validate file inputs
    //
    public static void validateFileInputs(params, log) {
        // Check input samplesheet exists
        if (params.input && !file(params.input).exists()) {
            Nextflow.error "Input samplesheet not found: ${params.input}"
        }

        // Check marker genes file if provided as file path
        if (params.marker_genes && params.marker_genes.endsWith('.txt')) {
            if (!file(params.marker_genes).exists()) {
                Nextflow.error "Marker genes file not found: ${params.marker_genes}"
            }
        }

        // Check MultiQC config if provided
        if (params.multiqc_config && !file(params.multiqc_config).exists()) {
            log.warn "MultiQC config file not found: ${params.multiqc_config}"
        }
    }

    //
    // Validate comparison files
    //
    public static void validateComparisonFiles(params, log) {
        if (params.da_comparisons && params.da_comparisons != "default_pairwise") {
            if (!file(params.da_comparisons).exists()) {
                Nextflow.error "DA comparisons file not found: ${params.da_comparisons}"
            }
            validateComparisonFileFormat(params.da_comparisons, "DA", log)
        }

        if (params.scviva_comparisons && params.scviva_comparisons != "default_pairwise") {
            if (!file(params.scviva_comparisons).exists()) {
                Nextflow.error "scVIVA comparisons file not found: ${params.scviva_comparisons}"
            }
            validateComparisonFileFormat(params.scviva_comparisons, "scVIVA", log)
        }
    }

    //
    // Validate comparison file format
    //
    public static void validateComparisonFileFormat(filePath, type, log) {
        def file = file(filePath)

        if (filePath.endsWith('.json')) {
            try {
                def jsonSlurper = new JsonSlurper()
                def comparisons = jsonSlurper.parse(file)

                if (!(comparisons instanceof List)) {
                    Nextflow.error "${type} comparisons JSON must be a list of comparison objects"
                }

                comparisons.each { comparison ->
                    if (!comparison.name || !comparison.group1 || !comparison.group2) {
                        Nextflow.error "${type} comparison missing required fields: name, group1, group2"
                    }
                }

                log.info "✓ ${type} comparisons file validated: ${comparisons.size()} comparisons found"

            } catch (Exception e) {
                Nextflow.error "Invalid ${type} comparisons JSON format: ${e.message}"
            }
        } else if (filePath.endsWith('.csv')) {
            // Basic CSV validation
            def lines = file.readLines()
            if (lines.size() < 2) {
                Nextflow.error "${type} comparisons CSV must have header and at least one comparison"
            }

            def header = lines[0].split(',')
            if (!header.contains('name') || !header.contains('group1') || !header.contains('group2')) {
                Nextflow.error "${type} comparisons CSV must have columns: name, group1, group2"
            }

            log.info "✓ ${type} comparisons CSV validated: ${lines.size() - 1} comparisons found"
        } else {
            Nextflow.error "${type} comparisons file must be .json or .csv format"
        }
    }

    //
    // Validate marker genes
    //
    public static void validateMarkerGenes(params, log) {
        if (params.marker_genes) {
            if (params.marker_genes.endsWith('.txt')) {
                // File-based marker genes
                def file = file(params.marker_genes)
                def genes = file.readLines().findAll { it.trim() }
                if (genes.size() == 0) {
                    Nextflow.error "Marker genes file is empty: ${params.marker_genes}"
                }
                log.info "✓ Marker genes file validated: ${genes.size()} genes found"
            } else {
                // Comma-separated marker genes
                def genes = params.marker_genes.split(',').collect { it.trim() }.findAll { it }
                if (genes.size() == 0) {
                    Nextflow.error "No valid marker genes found in: ${params.marker_genes}"
                }
                log.info "✓ Marker genes validated: ${genes.size()} genes specified"
            }
        }
    }

    //
    // Validate epoch parameters
    //
    public static void validateEpochs(params, log) {
        if (params.max_epochs < 1 || params.max_epochs > 1000) {
            Nextflow.error "max_epochs must be between 1 and 1000, got: ${params.max_epochs}"
        }

        if (params.scviva_max_epochs < 1 || params.scviva_max_epochs > 1000) {
            Nextflow.error "scviva_max_epochs must be between 1 and 1000, got: ${params.scviva_max_epochs}"
        }

        if (params.num_samples < 1 || params.num_samples > 100) {
            Nextflow.error "num_samples must be between 1 and 100, got: ${params.num_samples}"
        }
    }

    //
    // Function to check if genome exists in config file
    //
    public static void genomeExistsError(params, log) {
        // This pipeline doesn't use reference genomes, so this is a placeholder
        // for consistency with nf-core template
    }

    //
    // Validate annotation label
    //
    public static void validateAnnotationLabel(params, log) {
        if (!params.annotation_label) {
            Nextflow.error "annotation_label parameter is required"
        }

        // Check for valid column name format
        if (params.annotation_label.contains(' ')) {
            log.warn "annotation_label contains spaces, this may cause issues. Consider using underscores."
        }
    }

    //
    // Generate summary of parameters
    //
    public static String paramsSummaryLog(workflow, params) {
        def summary_log = ''
        summary_log += "Core parameters:\n"
        summary_log += "  Input samplesheet : ${params.input}\n"
        summary_log += "  Output directory  : ${params.outdir}\n"
        summary_log += "  Annotation label  : ${params.annotation_label}\n"
        summary_log += "\nResolVI parameters:\n"
        summary_log += "  Max epochs        : ${params.max_epochs}\n"
        summary_log += "  Num samples       : ${params.num_samples}\n"
        summary_log += "  Marker genes      : ${params.marker_genes ?: 'All genes'}\n"
        summary_log += "\nscVIVA parameters:\n"
        summary_log += "  Max epochs        : ${params.scviva_max_epochs}\n"
        summary_log += "\nGPU configuration:\n"
        summary_log += "  Num GPUs          : ${params.num_gpus}\n"
        summary_log += "  Force CPU         : ${params.force_cpu}\n"
        summary_log += "  GPU memory frac   : ${params.gpu_memory_fraction}\n"
        summary_log += "\nResource limits:\n"
        summary_log += "  Max CPUs          : ${params.max_cpus}\n"
        summary_log += "  Max memory        : ${params.max_memory}\n"
        summary_log += "  Max time          : ${params.max_time}\n"

        return summary_log
    }
}

