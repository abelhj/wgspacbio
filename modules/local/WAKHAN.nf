process WAKHAN {
    label 'process_high'
    conda '/opt/conda/envs/bio'
  
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'quay.io/biocontainers/severus@sha256:df04bec4a0ae9c55104ff91d6063d8c7c58b05145db04ad481f5d3d9527a9b7d' :
        'quay.io/biocontainers/severus@sha256:df04bec4a0ae9c55104ff91d6063d8c7c58b05145db04ad481f5d3d9527a9b7d' }"

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

    wakhan \
      --threads ${task.cpus} \
      --target-bam $input_bam \
      --reference $refence_fasta \
      --genome-name ${meta.sample} \
      --out-dir ./ \
      --breakpoints $severus_sv_vcf \
      --loh-enable \
      --use-sv-haplotypes \
      --ploidy-range 1.0-6.0 \
      --purity-range 0.2-0.99 \
      --centromere /opt/wakhan/src/annotations/grch38.cen_coord.curated.bed \
      --tumor-phased-vcf $tumor_germline_vcf  --hets-ratio 0.25



    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        wakhan: 0.0  )
    END_VERSIONS
    """

}
