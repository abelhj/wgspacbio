process LONGPHASE {
    label 'process_high'
    conda '/opt/conda/envs/bio'
  
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'abelhj/longphase:2.0' :
        'abelhj/longphase:2.0' }"

    input:
        tuple val(meta), path(input_files)
        path (reference_fasta)
        path (reference_fasta_index)

    output:
	tuple val(meta), path("${meta.sample}.longphase.vcf.gz")    , emit: phased_vcf
        path  ("versions.yml")                                      , emit: versions

    script:
    """

       /opt/longphase/longphase phase -s ${meta.sample}.clair3.correct.vcf.gz \
       -b ${meta.sample}.sorted.bam \
       -r ${reference_fasta} \
       -t ${task.cpus} \
       -o ${meta.sample}.longphase \
       --pb 
       bgzip ${meta.sample}.longphase.vcf
       tabix -p vcf ${meta.sample}.longphase.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        longphase: \$(/opt/longphase/longphase -v | grep Version  )
    END_VERSIONS
    """

}
