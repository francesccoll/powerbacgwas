#!/usr/bin/env nextflow

process ancestral_state_reconstruction_vcf {

    input:
        file(vcf_pastml_in)
        file(tree_in)
        val(outprefix)
        path(tmp_dir)
        val(num_cores)

    // cpus "${num_cores}"

    publishDir "${tmp_dir}"

    output:
        file "${outprefix}.vcf.ancestral.csv"
        file "${outprefix}.vcf.ancestral_steps.csv"

    """
    ancestral_state_reconstruction.py --input_vcf_table "$vcf_pastml_in" --input_tree "$tree_in" --output_table "$outprefix".vcf.ancestral.csv --output_steps "$outprefix".vcf.ancestral_steps.csv --processes "$num_cores" --pastml_dir "$outprefix"_tmp_pastml_dir/ --rm_tmp
    """
}


process ancestral_state_reconstruction_pg {

    input:
        file(pg_pastml_in)
        file(tree_in)
        val(outprefix)
        path(tmp_dir)
        val(num_cores)

    // cpus "${num_cores}"

    publishDir "${tmp_dir}"

    output:
        file "${outprefix}.ancestral.csv"
        file "${outprefix}.ancestral_steps.csv"

    """
    ancestral_state_reconstruction_roary.py --input_pastml_table "$pg_pastml_in" --input_tree "$tree_in" --output_table "$outprefix".ancestral.csv --output_steps "$outprefix".ancestral_steps.csv --processes "$num_cores" --pastml_dir "$outprefix"_tmp_pastml_dir/ --rm_tmp
    """
}

// NOTE: the pull path to calculate_changes.py needs to be specified

process calculate_changes_per_region {

    input:
        file(vcf_anc_in)
        file(tree_in)
        file(regions_in)
        val(outprefix)
        path(tmp_dir)
        val(num_cores)

    // cpus "${num_cores}"

    publishDir "${tmp_dir}"

    output:
        file "${outprefix}.regions.ancestral_steps.csv"

    """
    if [ ! -d "$tmp_dir"/"$outprefix"_tmp_pastml_dir ]
    then
        mkdir "$tmp_dir"/"$outprefix"_tmp_pastml_dir
    fi
    calculate_changes_per_region.py --pastml_output_table "$vcf_anc_in" --input_tree "$tree_in" --regions_file "$regions_in" --processes "$num_cores" --output_steps "$outprefix".regions.ancestral_steps.csv --tmp_dir "$tmp_dir"/"$outprefix"_tmp_pastml_dir --calculate_changes_path /usr/local/bin/calculate_changes.py
    rm -r "$tmp_dir"/"$outprefix"_tmp_pastml_dir
    """
}
