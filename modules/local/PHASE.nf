process PHASE {
    label 'process_high'
    conda '/opt/conda/envs/bio'
  
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'quay.io/pacbio/hiphase@sha256:353b4ffdae4281bdd5daf5a73ea3bb26ea742ef2c36e9980cb1f1ed524a07482' :
        'quay.io/pacbio/hiphase@sha256:353b4ffdae4281bdd5daf5a73ea3bb26ea742ef2c36e9980cb1f1ed524a07482' }"

    input:
        tuple val(meta), path(input_files)
        path (reference_fasta)
        path (reference_fasta_index)

    output:
	tuple val(meta), path("${meta.sample}.hiphase.bam")    , emit: phased_bam
	tuple val(meta), path("${meta.sample}.hiphase.vcf.gz")    , emit: phased_vcf
	tuple val(meta), path("${meta.sample}.somatic.hiphase.vcf.gz")    , emit: phased_somatic_vcf
	tuple val(meta), path("${meta.sample}.hiphase.*.t*"), emit: phasing_files
        path  ("versions.yml")                                      , emit: versions

    script:
    """

       hiphase --bam  ${meta.sample}.sorted.bam \
        -t ${task.cpus} \
        --output-bam ${meta.sample}.hiphase.bam \
        --vcf ${meta.sample}.clair3.correct.vcf.gz \
        --output-vcf ${meta.sample}.hiphase.vcf.gz \
        --vcf ${meta.sample}.deepsomatic.PASS.fixed.vcf.gz \
        --output-vcf ${meta.sample}.somatic.hiphase.vcf.gz \
        -r ${reference_fasta} \
        --stats-file ${meta.sample}.hiphase.stats.txt \
        --summary-file ${meta.sample}.hiphase.summary.tsv \
        --blocks-file ${meta.sample}.hiphase.blocks.tsv \
        --ignore-read-groups


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        hiphase: \$(hiphase --version  )
    END_VERSIONS
    """

}
