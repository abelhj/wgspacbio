process ANNOTATE_SV_VEP {
    tag "${meta.sample}"
    label "process_low"

    container "ghcr.io/dhslab/docker-vep_release113:250810"

    input:
    tuple val(meta), path(vcf)
    path(reference_fasta)
    path(reference_fasta_index)
    path(vep_cache)
    path(cytobands)
    path(sv_bed)


    output:
    tuple val(meta), path("*.vep.vcf.gz*"), emit: vcf
    path("versions.yml")                        , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """

    distance=20000

    tabix -p bed $cytobands
    tabix -p bed $sv_bed

    /usr/bin/perl -I /opt/lib/perl/VEP/Plugins /opt/vep/src/ensembl-vep/vep --format vcf --vcf \
    --fasta ${reference_fasta} --flag_pick --symbol --distance \$distance \
    --term SO -o ${meta.sample}.severus.svpack.vep.vcf \
    -i ${meta.sample}.severus.svpack.vcf \
    --force_overwrite \
    --custom ${cytobands},cytobands,bed --custom ${sv_bed},KnownSvGenes,bed \
    --offline --cache --dir ${vep_cache}

    cat ${meta.sample}.severus.svpack.vep.vcf | grep -v "^#" \
    | sort -k1,1V -k2,2n \
    | cat <(cat ${meta.sample}.severus.svpack.vep.vcf | grep "^#" ) - \
    | bgzip -c > ${meta.sample}.severus.svpack.vep.vcf.gz

    tabix -p vcf ${meta.sample}.severus.svpack.vep.vcf.gz


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vep: \$(/opt/vep/src/ensembl-vep/vep 2>&1 | grep ensembl-vep | awk -F ': ' '{print \$NF}')
    END_VERSIONS
    """

}