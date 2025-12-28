process WHATSHAP {
    label 'process_high'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'abelhj/whatshap:1.0.0' :
        'abelhj/whatshap:1.0.0' }"

    input:
        tuple val(meta), path(bam_bai_vcf_files)
        path(reference_fasta)
        path(reference_fasta_index)

    output:
        tuple val(meta), path("${meta.sample}*.haplotagged.bam")     , emit: bam
        tuple val(meta), path("${meta.sample}*.haplotagged.bam.bai") , emit: bai
        tuple val(meta), path("${meta.sample}*.phased.vcf.gz")       , emit: phased_vcf
        tuple val(meta), path("${meta.sample}*.phased.vcf.gz.tbi")   , emit: phased_vcf_tbi    
        path  ("versions.yml")                                       , emit: versions

    script:
    """
    
    echo ""

    whatshap phase --output ${meta.sample}.phased.vcf.gz \\
    --reference $reference_fasta \\
    ${meta.sample}.deepvariant.vcf.gz \\
    ${meta.sample}.sorted.bam

    bcftools index -t ${meta.sample}.phased.vcf.gz

    whatshap haplotag --output ${meta.sample}.haplotagged.bam \\
    --reference $reference_fasta \\
    ${meta.sample}.phased.vcf.gz \\
    ${meta.sample}.sorted.bam

    samtools index ${meta.sample}.haplotagged.bam 


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        whatshap: \$(whatshap --version)
    END_VERSIONS
    """
}
