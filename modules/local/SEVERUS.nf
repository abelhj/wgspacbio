process SEVERUS {
    label 'process_high'
    conda '/opt/conda/envs/bio'
  
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'gokcekeskus/severus:v1_6' :
        'gokcekeskus/severus:v1_6' }"

    input:
        tuple val(meta), path(input_files)
        path (reference_fasta)
        path (reference_fasta_index)
        path (trf)
        path (pon)

    output:
        tuple val(meta), path("${meta.sample}*_severus_out")    , emit: severus_output
        tuple val(meta), path("${meta.sample}*_severus_out/somatic_SVs/severus_somatic.vcf"), emit: severus_vcf
        path  ("versions.yml")                                      , emit: versions

    script:
    """
        echo "hi"
        severus --target-bam ${meta.sample}.haplotagged.bam \
        --out-dir ${meta.sample}_severus_out \
        -t ${task.cpus} \
        --phasing-vcf ${meta.sample}.rephased.vcf.gz \
        --vntr-bed ${trf} \
        --PON ${pon} \
        --output-read-ids \
        --min-reference-flank 0 \
        --single-bp \
        --resolve-overlaps \
        --max-unmapped-seq 7000 \
        --between-junction-ins 


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        severus: \$(severus --version)
    END_VERSIONS
    """

}
