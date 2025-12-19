process WAKHAN_CNA {
    label 'process_high'
    conda '/opt/conda/envs/bio'
  
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'abelhj/wakhan' :
        'abelhj/wakhan' }"

    input:
        tuple val(meta), path(input_files)
        path (reference_fasta)
        path (reference_fasta_index)
        path (cnv_genes)

    output:
        tuple val(meta), path("${meta.sample}*_wakhan_cna_out")    , emit: wakhan_output
        path  ("versions.yml")                                      , emit: versions

    script:
    """

        cp -r ${meta.sample}.wakhan_out/* .



        python3 /opt/wakhan/Wakhan/wakhan.py  cna --threads ${task.cpus} \
        --reference ${reference_fasta}\
        --target-bam ${meta.sample}.haplotagged.bam \
        --tumor-phased-vcf ${meta.sample}.rephased.vcf.gz \
        --genome-name ${meta.sample} \
        --cancer-genes ${cnv_genes} \
        --use-sv-haplotypes \
        --out-dir-plots . \
        --bin-size 50000  \
        --breakpoints severus_somatic.vcf \
        --phaseblocks-enable \
        --contigs 'chr1-22,chrX' \
        --copynumbers-subclonal-enable \
        --loh-enable

        mkdir -p ${meta.sample}_wakhan_cna_out
        find . -mindepth 1 -maxdepth 1 -type d ! -name "${meta.sample}_wakhan_cna_out" -print0 | xargs -0 -I {} mv "{}" ${meta.sample}_wakhan_cna_out/
        find . -maxdepth 1 -type f -name "*.html" -print0 | xargs -0 -I {} mv "{}" ${meta.sample}_wakhan_cna_out/

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        wakhan: wakhan:abelhj
    END_VERSIONS
    """

}
