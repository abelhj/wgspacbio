process CPG_PILEUP {
    label 'process_high'
    conda '/opt/conda/envs/bio'
  
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'quay.io/pacbio/pb-cpg-tools@sha256:b95ff1c53bb16e53b8c24f0feaf625a4663973d80862518578437f44385f509b' :
        'quay.io/pacbio/pb-cpg-tools@sha256:b95ff1c53bb16e53b8c24f0feaf625a4663973d80862518578437f44385f509b' }"

    input:
        tuple val(meta), path(bam_bai_files)
        path (reference_fasta)
        path (reference_fasta_index)
        path (input_bam)
        path (phased_vcf)
        path (PON_tsv)
        path (trf_bed)

    output:
        tuple val(meta), path("./somatic_SVs/${meta.sample}.severus_somatic.vcf")    , emit: somatic_vcf
        tuple val(meta), path("./all_SVs/${meta.sample}.severus_all.vcf")  , emit: all_vcf
        tuple val(meta), path("./somatic_SVs/${meta.sample}.somatic.breakpoint_clusters_list.tsv")  , emit: somatic_clusters
        tuple val(meta), path("./all_SVs/${meta.sample}.all.breakpoint_clusters_list.tsv")  , emit: all_clusters
        path  ("versions.yml")                                      , emit: versions

    script:
    """

    aligned_bam_to_cpg_scores \
      --threads ${tasks.cpu} \
      --bam $bam \
      --ref $reference_fasta} \
      --output-prefix ${meta.sample} \
      --min-mapq 1 \
      --min-coverage 5 \
      --model "$PILEUP_MODEL_DIR"/pileup_calling_model.v1.tflite

    # gzip all bed
    gzip ~{output_prefix}.combined.bed
    # If hap1 and hap2 beds are present, gzip them
    if [ -f ~{output_prefix}.hap1.bed ]; then
      gzip ~{output_prefix}.hap1.bed
    fi
    if [ -f ~{output_prefix}.hap2.bed ]; then
      gzip ~{output_prefix}.hap2.bed
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        aligned_bam_to_cpg_scores: \$(aligned_bam_to_cpg_scores --version  )
    END_VERSIONS
    """

}
