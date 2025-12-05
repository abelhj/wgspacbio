process PBMM2_ALIGNER {
    label 'process_high'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'quay.io/biocontainers/pbmm2@sha256:c1ec77296850cbdb02621bca1addcc25e510aacdabbad753ab3b0b8ba43ccd52' :
        'quay.io/biocontainers/pbmm2@sha256:c1ec77296850cbdb02621bca1addcc25e510aacdabbad753ab3b0b8ba43ccd52' }"

    input:

        tuple val(meta), path (reads_paths) 
        path (index)

    output:
        tuple val(meta), path ("*.bam")       , emit: bam
        path "versions.yml"                   , emit: versions

    script:
        def args = task.ext.args ?: ''
        """        
        pbmm2 align --preset CCS -j ${task.cpus}  \\
        --unmapped --log-level INFO \\
        ${index} ${reads_paths} ${meta.sample}.bam
        
        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            pbmm2: \$(pbmm2 --version 2>&1)
        END_VERSIONS
        """
}
