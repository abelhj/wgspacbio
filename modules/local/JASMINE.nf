process JASMINE {
    label 'process_xtra_high'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'mdivr/basic-tools:v11' :
        'mdivr/basic-tools:v11' }"

    input:

        tuple val(meta), path (reads_paths) 

    output:
        tuple val(meta), path ("*.bam")       , emit: bam
        path "versions.yml"                   , emit: versions

    script:
        def args = task.ext.args ?: ''
        """        
        /opt/jasmine/bin/jasmine -j ${task.cpus} --6mA=False \\
        ${reads_paths} ${meta.sample}.bam
        
        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            jasmine: \$(/opt/jasmine/bin/jasmine --version 2>&1)
        END_VERSIONS
        """
}
