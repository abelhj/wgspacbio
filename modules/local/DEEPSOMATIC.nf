process DEEPSOMATIC {
    label 'process_high'
    conda '/opt/conda/envs/bio'
  
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'abelhj/deepsomatic:1.9.0' :
        'abelhj/deepsomatic:1.9.0' }"

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
        /opt/deepvariant/bin/deepsomatic/run_deepsomatic \
        --model_type=PACBIO_TUMOR_ONLY \
        --use_default_pon_filtering=true \
        --ref=$reference_fasta \
        --reads_tumor=${meta.sample}.sorted.bam \
        --output_vcf=${meta.sample}.deepsomatic.vcf.gz \
        --output_gvcf=${meta.sample}.deepsomatic.g.vcf.gz \
        --num_shards=${task.cpus} 

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        deepsomatic: \$(/opt/deepvariant/bin/deepsomatic/run_deepsomatic --version | grep version )
    END_VERSIONS
    """

}
