process FILTER_SV {
    tag "filter_sv"
    label 'process_single'

    conda "conda-forge::python=3.8.3"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'ghcr.io/dhslab/docker-cleutils:240129' :
        'ghcr.io/dhslab/docker-cleutils:240129' }"

    input:
    tuple val(meta), path(vcf)
    path(sv_bed)
    path(genes_bed)
    path(recurrent_sv)
    path(hprc_ccdg_ins)


    output:
    tuple val(meta), path("${meta.sample}*.sv.combined.vcf.gz"), emit: vcxf
    tuple val(meta), path("${meta.sample}*.sv.combined.filtered.vcf.gz"), emit: filtered_vcf
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script: // This script is bundled with the pipeline, in nf-core/wgsnano/bin/
    """

    cat ${meta.sample}.severus.svpack.vcf \
      | grep -v "^#" \
      | python /storage2/fs1/dspencer/Active/spencerlab/abelhj/pacbio/wgspacbio/bin/zjoin3  \
        -v -a stdin -b <(zcat  ${meta.sample}.severus.svpack.vep.vcf.gz | grep -v "^#" ) -1 3 -2 3 \
      | cat  <(zcat ${meta.sample}.severus.svpack.vep.vcf.gz | grep -v "^#") - \
      | sort -k1,1V -k2,2n \
      | cat  <(zcat ${meta.sample}.severus.svpack.vep.vcf.gz | grep "^#") - \
      | bgzip -c >  ${meta.sample}.sv.combined.vcf.gz

    tabix -p vcf ${meta.sample}.sv.combined.vcf.gz


    python  /storage2/fs1/dspencer/Active/spencerlab/abelhj/pacbio/wgspacbio/bin/filter_sv_lr.3.py \
       -b $sv_bed \
       -g $recurrent_sv \
       -c $genes_bed    \
       -n $hprc_ccdg_ins \
       --caller severus \
       -o ${meta.sample}.sv.combined.filtered.vcf \
       ${meta.sample}.sv.combined.vcf.gz

    bgzip ${meta.sample}.sv.combined.filtered.vcf
    tabix -p vcf ${meta.sample}.sv.combined.filtered.vcf.gz


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
