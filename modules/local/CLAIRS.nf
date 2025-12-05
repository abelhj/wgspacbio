process CLAIRS {
    label 'process_high'
    conda '/opt/conda/envs/bio'
  
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker pull hkubal/clair3@sha256:857af16c759b0893fc757511a17c1efdfe253cbb64dffbcc8eecac0d33a60f60' :
        'docker pull hkubal/clair3@sha256:857af16c759b0893fc757511a17c1efdfe253cbb64dffbcc8eecac0d33a60f60' }"

    input:
        tuple val(meta), path(bam_bai_files)
        path (reference_fasta)
        path (reference_fasta_index)

    output:
        tuple val(meta), path("clair3_${meta.sample}*.deepvariant.vcf.gz")    , emit: vcf
        path  ("versions.yml")                                      , emit: versions

    script:
    """
        /opt/bin/run_clair3.sh \
        --bam_fn=$meta.sample}.sorted.bam \
        --ref_fn=$reference_fasta \
        --threads=$threads \
        --platform="hifi" \
        --model_path="/opt/models/hifi_revio" \
        --output=clair3_${meta.sample} \
        --sample_name=${meta.sample}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        clairs: \$(run_deepvariant --version | grep version )
    END_VERSIONS
    """

}
