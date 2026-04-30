process FILTER_CNV {
    tag "filter_sv"
    label 'process_single'

    conda "conda-forge::python=3.8.3"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'ghcr.io/dhslab/docker-cleutils:240129' :
        'ghcr.io/dhslab/docker-cleutils:240129' }"

    input:
    tuple val(meta), path(vcfs)
    path(sv_bed)
    path(genes_bed)
    path(cytobands)
    path(centromeres)


    output:
    tuple val(meta), path("${meta.sample}*.filtered.txt"), emit: txt
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script: // This script is bundled with the pipeline, in nf-core/wgsnano/bin/
    """


    clonal_vcf=${meta.sample}_wakhan_cna_out/*/vcf_output/*wakhan_cna_integers.vcf
    subclonal_vcf=${meta.sample}_wakhan_cna_out/*/vcf_output/*wakhan_cna_subclonals.vcf

    python /storage2/fs1/dspencer/Active/spencerlab/abelhj/pacbio/wgspacbio/bin/filter_cnv.py -v \$clonal_vcf -b \$subclonal_vcf \
    --genebed $genes_bed --svgenebed $sv_bed \
    --outfile ${meta.sample}.cnv.filtered.txt --cb $cytobands \
    --min_cnv_size 100000 -c $centromeres



    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
