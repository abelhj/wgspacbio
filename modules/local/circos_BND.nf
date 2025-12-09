process circos_BND {
    label 'process_high'
    conda '/opt/conda/envs/bio'
  
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'quay.io/pacbio/somatic_general_tools@sha256:99159e099d9044c52debabdc9491b168487aaa37534c1a748809bc69f169812a' :
        'quay.io/pacbio/somatic_general_tools@sha256:99159e099d9044c52debabdc9491b168487aaa37534c1a748809bc69f169812a' }"

    input:
        tuple val(meta), path(bam_bai_files)
        path (sv_vcf)

    output:
        tuple val(meta), path("${meta.sample}.circos.pdf")    , emit: circos_pdf
        tuple val(meta), path("${meta.sample}.circos.png")    , emit: circos_png
        tuple val(meta), path("${meta.sample}.circos.mitelman_fusions.tsv")    , emit: mitelman
        tuple val(meta), path("${meta.sample}.circos.known_gene_pairs.tsv")    , emit: known_genes
        path  ("versions.yml")                                      , emit: versions

    script:
    """

     python /app/plot_circos.py \
      $sv_vcf \
      ${meta.sample}.circos \
      /app/MCGENE.TXT.DATA \
      --max-gene-labels ${max_gene_labels}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        circos: 1.0
    END_VERSIONS
    """

}
