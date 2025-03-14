process DEEPVARIANT {
    label 'process_high'
    conda '/opt/conda/envs/bio'
  
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'abelhj/dv_fix' :
        'abelhj/dv_fix' }"

    input:
        tuple val(meta), path(bam_bai_files)
        path (reference_fasta)
        path (reference_fasta_index)

    output:
        tuple val(meta), path("${meta.sample}*.deepvariant.vcf.gz")    , emit: vcf
        tuple val(meta), path("${meta.sample}*.deepvariant.g.vcf.gz")  , emit: gvcf
        path  ("versions.yml")                                      , emit: versions

    script:
    """
        export PATH=/opt/deepvariant/bin/:$PATH
        run_deepvariant --model_type PACBIO \\
        --ref $reference_fasta \\
        --reads ${meta.sample}.sorted.bam \\
        --output_gvcf ${meta.sample}.deepvariant.g.vcf.gz \\
        --output_vcf ${meta.sample}.deepvariant.vcf.gz \\
        --num_shards ${task.cpus} \\

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        deepvariant: \$(run_deepvariant --version | grep version )
    END_VERSIONS
    """

}
