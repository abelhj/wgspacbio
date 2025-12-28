process POST_DEEPSOMATIC {
    label 'process_high'
    conda '/opt/conda/envs/bio'
  
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'quay.io/biocontainers/bcftools:1.17--h3cc50cf_1' :
        'quay.io/biocontainers/bcftools:1.17--h3cc50cf_1' }"

    input:
        tuple val(meta), path(input_vcf)
        path (reference_fasta)
        path (reference_fasta_index)


    output:
	tuple val(meta), path("${meta.sample}*.deepsomatic.PASS.vcf.gz")    , emit: pass_vcf
	tuple val(meta), path("${meta.sample}*.deepsomatic.PASS.fixed.vcf.gz")    , emit: fixed_vcf
	tuple val(meta), path("${meta.sample}*.deepsomatic.PASS.fixed.vcf.gz.tbi")    , emit: fixed_vcf_tbi
        path  ("versions.yml")                                      , emit: versions

    script:
    """
       

        bcftools view \
          -f PASS -Oz \
          -o ${meta.sample}.deepsomatic.PASS.vcf.gz \
	${meta.sample}.deepsomatic.vcf.gz

        echo ${meta.sample} > tmp.txt

        # Set GT to 0/1 for somatic variants as DeepSomatic set 1/1 which will not be phased by HiPhase
        bcftools +setGT  ${meta.sample}.deepsomatic.PASS.vcf.gz -- -t q -n c:"0/1" -i 'FMT/GT="1/1"' \
        | bcftools reheader -s tmp.txt \
        | bcftools sort -Oz -o ${meta.sample}.deepsomatic.PASS.fixed.vcf.gz

        tabix -p vcf ${meta.sample}.deepsomatic.PASS.vcf.gz 
        tabix -p vcf ${meta.sample}.deepsomatic.PASS.fixed.vcf.gz


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version  )
    END_VERSIONS
    """

}
