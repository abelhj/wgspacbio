process HAPLOTAG {
    label 'process_high'
    conda '/opt/conda/envs/bio'


    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'mkolmogo/whatshap:2.3' :
        'mkolmogo/whatshap:2.3' }"

    input:
        tuple val(meta), path(input_files)
        path (reference_fasta)
        path (reference_fasta_index)

    output:
        tuple val(meta), path("${meta.sample}.haplotagged.bam"), emit: bam
        tuple val(meta), path("${meta.sample}.haplotagged.bam.bai"), emit: bai
        path  ("versions.yml")                                      , emit: versions

    script:
    """

        whatshap haplotag --reference ${reference_fasta} \
        ${meta.sample}.rephased.vcf.gz ${meta.sample}.sorted.bam \
        -o ${meta.sample}.haplotagged.bam \
        --ignore-read-groups \
        --tag-supplementary --skip-missing-contigs --output-threads 4
        samtools index -@8 ${meta.sample}.haplotagged.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        whatshap: \$(whatshap --version)
    END_VERSIONS
    """

}
