//
// This file holds several functions specific to the main.nf workflow in the nf-core/qualitycontrol and alignassembly pipeline
//

class WorkflowMain {

    //
    // Citation string for pipeline
    //
    public static String citation(workflow) {
        return "If you use ${workflow.manifest.name} for your analysis please cite:\n\n" +
            // TODO nf-core: Add Zenodo DOI for pipeline after first release
            //"* The pipeline\n" +
            //"  https://doi.org/10.5281/zenodo.XXXXXXX\n\n" +
            "* The nf-core framework\n" +
            "  https://doi.org/10.1038/s41587-020-0439-x\n\n" +
            "* Software dependencies\n" +
            "  https://github.com/${workflow.manifest.name}/blob/master/CITATIONS.md"
    }

    //
    // Generate help string
    //
    public static String help(workflow, params, log) {
        def command = "nextflow run ${workflow.manifest.name} --input samplesheet.csv --genome GRCh37 -profile docker"
        def help_string = ''
        help_string += NfcoreTemplate.logo(workflow, params.monochrome_logs)
        help_string += NfcoreSchema.paramsHelp(workflow, params, command)
        help_string += '\n' + citation(workflow) + '\n'
        help_string += NfcoreTemplate.dashedLine(params.monochrome_logs)
        return help_string
    }

    //
    // Generate parameter summary log string
    //
    public static String paramsSummaryLog(workflow, params, log) {
        def summary_log = ''
        summary_log += NfcoreTemplate.logo(workflow, params.monochrome_logs)
        summary_log += NfcoreSchema.paramsSummaryLog(workflow, params)
        summary_log += '\n' + citation(workflow) + '\n'
        summary_log += NfcoreTemplate.dashedLine(params.monochrome_logs)
        return summary_log
    }

    //
    // Validate parameters and print summary to screen
    //
    public static void initialise(workflow, params, log) {
        // Print help to screen if required
        if (params.help) {
            log.info help(workflow, params, log)
            System.exit(0)
        }

        // Print workflow version and exit on --version
        if (params.version) {
            String workflow_version = NfcoreTemplate.version(workflow)
            log.info "${workflow.manifest.name} ${workflow_version}"
            System.exit(0)
        }

        // Print parameter summary log to screen
        log.info paramsSummaryLog(workflow, params, log)

        // Validate workflow parameters via the JSON schema
        if (params.validate_params) {
            NfcoreSchema.validateParameters(workflow, params, log)
        }

        // Check that a -profile or Nextflow config has been provided to run the pipeline
        NfcoreTemplate.checkConfigProvided(workflow, log)

        // Check that conda channels are set-up correctly
        if (params.enable_conda) {
            Utils.checkCondaChannels(log)
        }

        // Check AWS batch settings
        NfcoreTemplate.awsBatch(workflow, params)

        // Check input has been provided
        if (!params.input) {
            log.error "Please provide an input samplesheet to the pipeline e.g. '--input samplesheet.csv'"
            System.exit(1)
        }
    }
    //
    // Get attribute from genome config file e.g. fasta
    //
    public static Object getGenomeAttribute(params, attribute) {
        if (params.genomes && params.genome && params.genomes.containsKey(params.genome)) {
            if (params.genomes[ params.genome ].containsKey(attribute)) {
                return params.genomes[ params.genome ][ attribute ]
            }
        }
        return null
    }



    // Check if input ids are from the SRA
    public static Boolean isSraId(input, log) {
        def is_sra = false
        def total_ids = 0
        def no_match_ids = []
        def pattern = /^(((SR|ER|DR)[APRSX])|(SAM(N|EA|D))|(PRJ(NA|EB|DB))|(GS[EM]))(\d+)$/
        input.eachLine { line ->
            total_ids += 1
            if (!(line =~ pattern)) {
                no_match_ids << line
            }
        }

        def num_match = total_ids - no_match_ids.size()
        if (num_match > 0) {
            if (num_match == total_ids) {
                is_sra = true
            } else {
                log.error "Mixture of ids provided via --input: ${no_match_ids.join(', ')}\nPlease provide either SRA / ENA / DDBJ or Synapse ids!"
                System.exit(1)
            }
        }
        return is_sra
    }

    // Check if input ids are from the Synapse platform
    public static Boolean isSynapseId(input, log) {
        def is_synapse = false
        def total_ids = 0
        def no_match_ids = []
        def pattern = /^syn\d{8}$/
        input.eachLine { line ->
            total_ids += 1
            if (!(line =~ pattern)) {
                no_match_ids << line
            }
        }

        def num_match = total_ids - no_match_ids.size()
        if (num_match > 0) {
            if (num_match == total_ids) {
                is_synapse = true
            } else {
                log.error "Mixture of ids provided via --input: ${no_match_ids.join(', ')}\nPlease provide either SRA / ENA / DDBJ or Synapse ids!"
                System.exit(1)
            }
        }
        return is_synapse
    }



}
