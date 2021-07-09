#!/usr/bin/env nextflow

process prepare_vcf_file {
    // NOTE: the script prepare_vcf_file.py has optional arguments --bcftools_path, --bgzip_path and --tabix_path that are not used when using the docker container

    input:
        file(vcf_in)
        val(outprefix)
        path(tmp_dir)

    publishDir "${tmp_dir}"

    // NOTE: the output processed VCF file will be used by vcf_to_plink_files and vcf_to_pastml_matrix
    output:
        // file("${outprefix}.reformatted.vcf"), emit: vcf_reformatted
        file("${outprefix}.reformatted.vcf")
        file("${outprefix}.reformatted.vcf.gz")
        file("${outprefix}.reformatted.vcf.gz.tbi")

    """
    prepare_vcf_file.py --input_vcf "$vcf_in" --output_vcf "$outprefix".reformatted.vcf --to_haploid_format
    """
}


process vcf_to_plink_files {

    // NOTE: the script vcf_to_plink_files.py has optional arguments --plink_path and --bcftools_path which are not used when using the docker container

    input:
        file(vcf_rf_in)
  	val(outprefix)
        path(tmp_dir)

    publishDir "${tmp_dir}"

    output:
         file "${outprefix}.bed"
         file "${outprefix}.bim"
         file "${outprefix}.fam"
         file "${outprefix}.log"
         file "${outprefix}.nosex"

    """
    vcf_to_plink_files.py --input_vcf "$vcf_rf_in" --output_prefix "$outprefix"
    """
}


process vcf_to_pastml_matrix {

    input:
        file(vcf_rf_in)
  	val(outprefix)
        path(tmp_dir)

    publishDir "${tmp_dir}"

    output:
        file "${outprefix}.vcf.pastml.csv"

    """
    vcf_to_pastml_matrix.py --input_vcf "$vcf_rf_in" --output_table "$outprefix".vcf.pastml.csv
    """
}



process roary_to_plink_files {

    // NOTE: the script roary_to_plink_files.py optional argument --plink_path not needed when using the docker container
    // NOTE: the script roary_to_plink_files.py optional argument --samples_list not used in this process
    // NOTE: Roary file format .Rtab assumed in this process

    input:
        file(pg_in)
  	val(outprefix)
        path(tmp_dir)

    publishDir "${tmp_dir}"

    output:
         file "${outprefix}.bed"
         file "${outprefix}.bim"
         file "${outprefix}.fam"
         file "${outprefix}.log"
         file "${outprefix}.nosex"


    """
    roary_to_plink_files.py --gene_presence_absence "$pg_in" --input_format 2 --output_prefix "$outprefix"
    """
}


process roary_to_pastml_matrix {

    // NOTE: the script roary_to_pastml_matrix.py optional argument --samples_list not used in this process

    input:
        file(pg_in)
  	val(outprefix)
        path(tmp_dir)

    publishDir "${tmp_dir}"
    output:
        file "${outprefix}.pastml.csv"

    """
    roary_to_pastml_matrix.py --gene_presence_absence "$pg_in" --input_format 2 --output_table "$outprefix".pastml.csv
    """
}



process annotate_nodes {

    input:
        file(tree_in)
  	val(outprefix)
        path(tmp_dir)

    publishDir "${tmp_dir}"

    output:
        file "${outprefix}.tree.annotated.nwk"

    // NOTE: the single quotes introduced by script annotate_nodes_newick.py must be removed before running PastML

    """
    annotate_nodes_newick.py --input_tree "$tree_in" --node_prefix node --output_tree "$outprefix".tree.annotated.nwk
    """
}




