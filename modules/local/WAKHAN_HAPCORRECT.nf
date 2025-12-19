process WAKHAN_HAPCORRECT {
    label 'process_high'
    conda '/opt/conda/envs/bio'
  
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'mkolmogo/wakhan:dev_14f2192' :
        'mkolmogo/wakhan:dev_14f2192' }"

    input:
        tuple val(meta), path(input_files)
        path (reference_fasta)
        path (reference_fasta_index)

    output:
        tuple val(meta), path("${meta.sample}.wakhan_out"),  emit: wakhanHPOutput
        tuple val(meta), path("${meta.sample}.rephased.vcf.gz"), emit: rephased_vcf
        tuple val(meta), path("${meta.sample}.rephased.vcf.gz.tbi"), emit: rephased_vcf_tbi
        path  ("versions.yml")                                      , emit: versions

    script:
    """

        wakhan hapcorrect --threads ${task.cpus} \
        --reference ${reference_fasta} \
        --target-bam ${meta.sample}.sorted.bam \
        --tumor-phased-vcf ${meta.sample}.longphase.vcf.gz \
        --genome-name ${meta.sample} \
        --out-dir-plots ${meta.sample}.wakhan_out \
        --bin-size 50000 \
        --phaseblocks-enable --contigs chr1-22,chrX,chrY \
        --copynumbers-subclonal-enable

    if [ ! -e ${meta.sample}.wakhan_out/phasing_output/${meta.sample}.rephased.vcf.gz ]; then
      cp ${meta.sample}.longphase.vcf.gz ${meta.sample}.rephased.vcf.gz
    else
      cp ${meta.sample}.wakhan_out/phasing_output/${meta.sample}.rephased.vcf.gz ${meta.sample}.rephased.vcf.gz
    fi

    tabix -p vcf ${meta.sample}.rephased.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        wakhan: wakhan:dev_14f2192
    END_VERSIONS
    """

}
