/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// Validate input parameters
WorkflowWgsnano.initialise(params, log)

// TODO nf-core: Add all file path parameters for the pipeline to the list below
// Check input path parameters to see if they exist
def checkPathParamList = [ params.input, params.multiqc_config, params.fasta ]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Check mandatory parameters
if (params.input) { ch_input = file(params.input) } else { exit 1, 'Input samplesheet not specified!' }
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

include { MERGE_BASECALL as MERGE_BASECALL_ID           } from '../modules/local/MERGE_BASECALL'
include { MERGE_BASECALL as MERGE_BASECALL_SAMPLE       } from '../modules/local/MERGE_BASECALL'
include { PBMM2_ALIGNER                                 } from '../modules/local/PBMM2_ALIGNER'
include { SAMTOOLS_SORT                                 } from '../modules/local/SAMTOOLS_SORT'
include { DEEPVARIANT                                   } from '../modules/local/DEEPVARIANT'
include { SAMTOOLS_INDEX                                } from '../modules/local/SAMTOOLS_INDEX'
include { MOSDEPTH                                      } from '../modules/local/MOSDEPTH'
include { CUSTOM_DUMPSOFTWAREVERSIONS                   } from '../modules/nf-core/custom/dumpsoftwareversions/main'
include { MULTIQC                                       } from '../modules/local/MULTIQC'
include { WHATSHAP                                      } from '../modules/local/WHATSHAP.nf'
include { SAMTOOLS_STATS                                } from '../modules/local/SAMTOOLS_STATS.nf'
include { CLAIR3                                        } from '../modules/local/CLAIR3.nf'
include { DEEPSOMATIC                                   } from '../modules/local/DEEPSOMATIC.nf'
include { POST_DEEPSOMATIC                              } from '../modules/local/POST_DEEPSOMATIC.nf'
include { POST_CLAIR3                                   } from '../modules/local/POST_CLAIR3.nf'
include { PHASE                                         } from '../modules/local/PHASE.nf'
include { ANNOTATE_VARIANTS                             } from '../modules/local/ANNOTATE_VARIANTS.nf'
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Info required for completion email and summary
def multiqc_report = []

workflow WGSPACBIO {

    ch_versions = Channel.empty()

    //
    // SUBWORKFLOW: Read in samplesheet, validate and stage input files
    //
    INPUT_CHECK (
        ch_input
    )
    ch_versions = ch_versions.mix(INPUT_CHECK.out.versions)
    ch_phased_vcf = INPUT_CHECK.out.reads.map{ meta, files -> [[sample: meta.sample], meta.vcf, meta.vcf_tbi] }.dump(tag: "ch_phased_vcf")


if (params.reads_format == 'bam' ) {
    println 'bam\n'
    INPUT_CHECK
    .out
    .reads
    .flatMap { meta, bam_path -> 
        def bam_files = []
        if (file(bam_path).isDirectory()) {
            bam_files = file("${bam_path}/*.bam")
        } else if (bam_path.endsWith('.bam')) {
            bam_files = [file(bam_path)]
        }
        bam_files.collect { [[sample: meta.sample], it] }  // Create a list of [meta, file] pairs
    }
    .groupTuple(by: 0) // group bams by meta (i.e sample) which is zero-indexed
    // .dump(tag: 'basecall_sample', pretty: true)
    .set { ch_basecall_sample_merged_bams } // set channel name
}
    println ch_basecall_sample_merged_bams
    MERGE_BASECALL_SAMPLE (
        ch_basecall_sample_merged_bams
    )
    ch_versions = ch_versions.mix(MERGE_BASECALL_SAMPLE.out.versions)
    println 'merge_done\n'

    PBMM2_ALIGNER (
        MERGE_BASECALL_SAMPLE.out.merged_bam,
        file(params.fasta)
    )
    ch_versions = ch_versions.mix(PBMM2_ALIGNER.out.versions)


    //
    // MODULE: Samtools sort and indedx aligned bams
    //
    SAMTOOLS_SORT (
        PBMM2_ALIGNER.out.bam
    )
    ch_versions = ch_versions.mix(SAMTOOLS_SORT.out.versions)

    //
    // MODULE: CLAIR3
    //
    ch_dv_input = SAMTOOLS_SORT.out.bam.mix(SAMTOOLS_SORT.out.bai).groupTuple(size:2).map{ meta, files -> [ meta, files.flatten() ]}
    ch_dv_input.dump(tag: "clair3")
    CLAIR3 (
        ch_dv_input,
        file(params.fasta),
        file(params.fasta_index)
    )
    ch_versions = ch_versions.mix(CLAIR3.out.versions)

    // MODULE: DEEPSOMATIC
    //
    ch_dv_input = SAMTOOLS_SORT.out.bam.mix(SAMTOOLS_SORT.out.bai).groupTuple(size:2).map{ meta, files -> [ meta, files.flatten() ]}
    ch_dv_input.dump(tag: "deepsomatic")
    DEEPSOMATIC (
        ch_dv_input,
        file(params.fasta),
        file(params.fasta_index)
    )
    ch_versions = ch_versions.mix(DEEPSOMATIC.out.versions)

    // MODULE: POST_DEEPSOMATIC
    //
    POST_DEEPSOMATIC (
        DEEPSOMATIC.out.vcf,
        file(params.fasta),
        file(params.fasta_index)
    )
    ch_versions = ch_versions.mix(POST_DEEPSOMATIC.out.versions)

    // MODULE: POST_CLAIR3
    //
    POST_CLAIR3 (
        CLAIR3.out.vcf,
        file(params.fasta),
        file(params.fasta_index)
    )
    ch_versions = ch_versions.mix(POST_CLAIR3.out.versions)

    ch_phase_input = SAMTOOLS_SORT.out.bam.mix(SAMTOOLS_SORT.out.bai, POST_DEEPSOMATIC.out.fixed_vcf, POST_DEEPSOMATIC.out.fixed_vcf_tbi, POST_CLAIR3.out.fixed_vcf, POST_CLAIR3.out.fixed_vcf_tbi).groupTuple(size:6).map{ meta, files -> [ meta, files.flatten() ]}
    // MODULE: PHASE
    //
    PHASE (
        ch_phase_input,
        file(params.fasta),
        file(params.fasta_index)
    )
    ch_versions = ch_versions.mix(PHASE.out.versions)

    //ch_annotate_input = PHASE.out.phased_vcf.mix(PHASE.out.somatic_phased_vcf.
    // MODULE: ANNOTATE_VARIANTS
    //
    ANNOTATE_VARIANTS (
      PHASE.out.phased_vcf,
      file(params.fasta),
      file(params.fasta_index),
      file(params.vep_cache),
      file(params.cytobands),
      file(params.custom_annotation_vcf)
   )

    ANNOTATE_VARIANTS (
      PHASE.out.somatic_phased_vcf,
      file(params.fasta),
      file(params.fasta_index),
      file(params.vep_cache),
      file(params.cytobands),
      file(params.custom_annotation_vcf)
   )






//    //
//    // MODULE: DEEPVARIANT
//    //
//    ch_dv_input = SAMTOOLS_SORT.out.bam.mix(SAMTOOLS_SORT.out.bai).groupTuple(size:2).map{ meta, files -> [ meta, files.flatten() ]}
//    ch_dv_input.dump(tag: "deepvariant")
//    DEEPVARIANT (
//        ch_dv_input,
//        file(params.fasta),
//        file(params.fasta_index)
//    )
//    ch_versions = ch_versions.mix(DEEPVARIANT.out.versions)
//
//    
//    if (params.run_whatshap) {
//        //
//        // MODULE: Index PEPPER bam
//        //
//        ch_whatshap_input = SAMTOOLS_SORT.out.bam.mix(SAMTOOLS_SORT.out.bai,DEEPVARIANT.out.vcf).groupTuple(size:3).map{ meta, files -> [ meta, files.flatten() ]}
//        WHATSHAP (
//            ch_whatshap_input,
//            file(params.fasta),
//            file(params.fasta_index)
//        )
//        ch_versions = ch_versions.mix(WHATSHAP.out.versions)
//
//
//        //
//        // MODULE: MOSDEPTH for depth calculation
//        //
//        ch_mosdepth_input = WHATSHAP.out.bam.mix(WHATSHAP.out.bai).groupTuple(size:2).map{ meta, files -> [ meta, files.flatten() ]}
//        MOSDEPTH (
//            ch_mosdepth_input
//        )
//        ch_versions = ch_versions.mix(MOSDEPTH.out.versions)
//
//
//
//        SAMTOOLS_STATS (
//            WHATSHAP.out.bam
//        )
//        ch_versions = ch_versions.mix(SAMTOOLS_STATS.out.versions)
//    }
//
//    CUSTOM_DUMPSOFTWAREVERSIONS (
//        ch_versions.unique().collectFile(name: 'collated_versions.yml')
//    )
//
//    //
//    // MODULE: MultiQC
//    //
//    workflow_summary    = WorkflowWgsnano.paramsSummaryMultiqc(workflow, summary_params)
//    ch_workflow_summary = Channel.value(workflow_summary)
//
//    methods_description    = WorkflowWgsnano.methodsDescriptionText(workflow, ch_multiqc_custom_methods_description)
//    ch_methods_description = Channel.value(methods_description)
//
//    ch_multiqc_files = Channel.empty()
//    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
//    ch_multiqc_files = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml'))
//    ch_multiqc_files = ch_multiqc_files.mix(CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect())
//
//    if (params.run_whatshap) {
//        ch_multiqc_files = ch_multiqc_files.mix(MOSDEPTH.out.global_txt.collect{it[1]}.ifEmpty([]))
//        ch_multiqc_files = ch_multiqc_files.mix(MOSDEPTH.out.summary_txt.collect{it[1]}.ifEmpty([]))
//        ch_multiqc_files = ch_multiqc_files.mix(MOSDEPTH.out.regions_txt.collect{it[1]}.ifEmpty([]))
//        ch_multiqc_files = ch_multiqc_files.mix(MOSDEPTH.out.regions_bed.collect{it[1]}.ifEmpty([]))
//        ch_multiqc_files = ch_multiqc_files.mix(MOSDEPTH.out.regions_csi.collect{it[1]}.ifEmpty([]))
//        ch_multiqc_files = ch_multiqc_files.mix(MOSDEPTH.out.quantized_bed.collect{it[1]}.ifEmpty([]))
//        ch_multiqc_files = ch_multiqc_files.mix(MOSDEPTH.out.quantized_csi.collect{it[1]}.ifEmpty([]))
//    }
//
//    MULTIQC (
//        ch_multiqc_files.collect(),
//        ch_multiqc_config.toList(),
//        ch_multiqc_custom_config.toList(),
//        ch_multiqc_logo.toList()
//    )
//    multiqc_report = MULTIQC.out.report.toList()
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
