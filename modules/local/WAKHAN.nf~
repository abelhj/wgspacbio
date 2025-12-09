process SEVERUS {
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


    severus \
      --target-bam $input_bam \
      --phasing-vcf $phased_vcf \
      --PON $PON_tsv \
      --out-dir . \
      -t ${task.cpus} \
      --vntr-bed ${trf_bed} \
      --min-support 2 \
      --resolve-overlaps \
      --between-junction-ins \
      --single-bp

    # Compress SVs plots HTML inside somatic_SVs/plots directory
    # Check if the directory exists first
    if [[ -d ./somatic_SVs/plots ]]
      then tar -czvf ${meta.sample}.plots.tar.gz ./somatic_SVs/plots
    fi

    mv ./somatic_SVs/severus_somatic.vcf ./somatic_SVs/${meta.sample}.severus_somatic.vcf
    mv ./all_SVs/severus_all.vcf ./all_SVs/${meta.sample}.severus_all.vcf
    mv ./somatic_SVs/breakpoint_clusters_list.tsv ./somatic_SVs/${meta.sample}.somatic.breakpoint_clusters_list.tsv
    mv ./all_SVs/breakpoint_clusters_list.tsv ./all_SVs/${meta.sample}.all.breakpoint_clusters_list.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        severus: \$(severus --version  )
    END_VERSIONS
    """

}
