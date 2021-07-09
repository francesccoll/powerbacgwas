#!/usr/bin/env nextflow

process phylogeny_distance {

    input:
        file(tree_in)
        val(outprefix)
        path(tmp_dir)

    publishDir "${tmp_dir}"

    output:
        file "${outprefix}.tree_distances.csv"

    """
    phylogeny_distance.py --calc-C "$tree_in" > "$outprefix".tree_distances.csv
    """
}

process prepare_gwas_runs_subsampling_pg {

    input:
        file(pg_in)
        file(similarity_in)
        file(parameters_in)
        file(causal_in)
        file(phenotype_in)
        val(outprefix)
        path(tmp_dir)
        path(gwas_tmp_dir)
        val(pyseer_model)
        val(num_cores)

    publishDir "${tmp_dir}"

    output:
       val true
       path "${outprefix}.gwas_runs.csv"
       path "${outprefix}.gwas_runs.sh"

    """
    prepare_gwas_runs_subsampling_roary.py --roary_table "$pg_in" --input_format 2 --similarity "$similarity_in" --parameters_file "$parameters_in" --causal-loci "$causal_in" --phenotype "$phenotype_in" --output_dir "$gwas_tmp_dir" --output_prefix "$outprefix" --pyseer_model "$pyseer_model" --cpu "$num_cores"

    """
}

process prepare_gwas_runs_subsampling_vcf {

    input:
        file(vcf_in)
        file(vcf_tbi_in)
        file(similarity_in)
        file(parameters_in)
        file(causal_in)
        file(phenotype_in)
        val(outprefix)
        path(tmp_dir)
        path(gwas_tmp_dir)
        val(pyseer_model)
        val(num_cores)

    publishDir "${tmp_dir}"

    output:
       val true
       path "${outprefix}.gwas_runs.csv"
       path "${outprefix}.gwas_runs.sh"

    """
    prepare_gwas_runs_subsampling.py --input_vcf "$vcf_in" --similarity "$similarity_in" --parameters_file "$parameters_in" --causal-loci "$causal_in" --phenotype "$phenotype_in" --output_dir "$gwas_tmp_dir" --output_prefix "$outprefix" --pyseer_model "$pyseer_model" --cpu "$num_cores"

    """
}

process prepare_gwas_runs_subsampling_burden {

    input:
        file(vcf_in)
        file(vcf_tbi_in)
        file(similarity_in)
        file(parameters_in)
        file(causal_in)
        file(phenotype_in)
        file(regions_in)
        val(outprefix)
        path(tmp_dir)
        path(gwas_tmp_dir)
        val(pyseer_model)
        val(num_cores)

    publishDir "${tmp_dir}"

    output:
       val true
       path "${outprefix}.gwas_runs.csv"
       path "${outprefix}.gwas_runs.sh"

    """
    prepare_gwas_runs_subsampling.py --input_vcf "$vcf_in" --similarity "$similarity_in" --parameters_file "$parameters_in" --causal-loci "$causal_in" --phenotype "$phenotype_in" --burden "$regions_in" --output_dir "$gwas_tmp_dir" --output_prefix "$outprefix" --pyseer_model "$pyseer_model" --cpu "$num_cores"

    """
}


process prepare_gwas_runs_phensim_pg {

    input:
        file(pg_in)
        file(similarity_in)
        file(parameters_in)
        file(pastml_steps_file)
        val(outprefix)
        file(plink_bed)
        file(plink_bim)
        file(plink_fam)
        path(tmp_dir)
        path(gwas_tmp_dir)
        val(pyseer_model)
        val(phenotype_type)
        val(num_cores)

    publishDir "${tmp_dir}"

    output:
       val true
       path "${outprefix}.gwas_runs.csv"
       path "${outprefix}.gwas_runs.sh"

    """
    prepare_gwas_runs_roary.py --roary_table "$pg_in" --input_format 2 --plink_prefix "$outprefix" --similarity "$similarity_in" --pastml_steps_file "$pastml_steps_file" --parameters_file "$parameters_in" --output_dir "$gwas_tmp_dir" --output_prefix "$outprefix" --pyseer_model "$pyseer_model" --phenotype_type "$phenotype_type" --cpu "$num_cores"

    """
}


process prepare_gwas_runs_phensim_vcf {

    input:
        file(vcf_in)
        file(vcf_tbi_in)
        file(similarity_in)
        file(parameters_in)
        file(pastml_steps_file)
        val(outprefix)
        file(plink_bed)
        file(plink_bim)
        file(plink_fam)
        path(tmp_dir)
        path(gwas_tmp_dir)
        val(pyseer_model)
        val(phenotype_type)
        val(num_cores)


    publishDir "${tmp_dir}"

    output:
       val true
       path "${outprefix}.gwas_runs.csv"
       path "${outprefix}.gwas_runs.sh"

    """
    prepare_gwas_runs.py --input_vcf "$vcf_in" --plink_prefix "$outprefix" --similarity "$similarity_in" --pastml_steps_file "$pastml_steps_file" --parameters_file "$parameters_in" --output_dir "$gwas_tmp_dir" --output_prefix "$outprefix" --pyseer_model "$pyseer_model" --phenotype_type "$phenotype_type" --cpu "$num_cores"

    """
}

process prepare_gwas_runs_phensim_burden {
 input:
        file(vcf_in)
        file(vcf_tbi_in)
        file(similarity_in)
        file(parameters_in)
        file(pastml_steps_file)
        val(outprefix)
        file(plink_bed)
        file(plink_bim)
        file(plink_fam)
        file(regions_in)
        path(tmp_dir)
        path(gwas_tmp_dir)
        val(pyseer_model)
        val(phenotype_type)
        val(num_cores)

    publishDir "${tmp_dir}"

    output:
       val true
       path "${outprefix}.gwas_runs.csv"
       path "${outprefix}.gwas_runs.sh"

    """
    prepare_gwas_runs.py --input_vcf "$vcf_in" --burden "$regions_in" --plink_prefix "$outprefix" --similarity "$similarity_in" --pastml_steps_file "$pastml_steps_file" --parameters_file "$parameters_in" --output_dir "$gwas_tmp_dir" --output_prefix "$outprefix" --pyseer_model "$pyseer_model" --phenotype_type "$phenotype_type" --cpu "$num_cores"

    """

}


// NOTE: https://github.com/nextflow-io/patterns/blob/master/docs/process-per-file-output.adoc

process collect_gwas_scripts {
    input:
        val(flag)
        path(gwas_tmp_dir)

    publishDir "${gwas_tmp_dir}"

    output:
        file "${gwas_tmp_dir}/*.gwas_run.sh"

    """
    ls "$gwas_tmp_dir"/*.gwas_run.sh

    """
}


// See https://bioinformatics.stackexchange.com/questions/15633/nextflow-run-sbatch-slurm-job-script-inside-a-process
// NOTE: because each process is run in a separate unique temporary working directory, all input files used by the gwas bash script must be defined
// NOTE: because pyseer may fail in certain conditions, errors in this process are ignored (https://github.com/nextflow-io/patterns/blob/master/docs/ignore-failing-process.adoc)
// NOTE: it was not needed to define the output files of gwas_bash_script in the process output (due to soft link to directory gwas_tmp_dir?)
// NOTE: if run_gwas_bash_script process failed, the consecutive process (i.e. using run_gwas_bash_script.out) did not launch
// NOTE: to do: design generic run_gwas_bash_script processes (depending on pipeline branch)

process run_gwas_bash_script_subsampling_pg {

    input:
        file(gwas_bash_script)
        path(gwas_tmp_dir)
        path(variant_in)
        path(similarity_in)
        path(causal_in)
        path(phenotype_in)

    output:
        val true
        // file "*.pyseer.phen"
        // file "*_*.pyseer_output.csv"

    errorStrategy 'ignore'

    script:
    """
    bash "$gwas_bash_script"
    """

}


process run_gwas_bash_script_phensim_pg {

    input:
        file(gwas_bash_script)
        path(gwas_tmp_dir)
        path(pangenome_in)
        path(similarity_in)
        path(pastml_steps)
        file(plink_bed)
        file(plink_bim)
        file(plink_fam)

    output:
        val true
        // file "*.pyseer.phen"
        // file "*_*.pyseer_output.csv"

    errorStrategy 'ignore'

    script:
    """
    bash "$gwas_bash_script"
    """

}
// here

process run_gwas_bash_script_phensim_vcf {

    input:
        file(gwas_bash_script)
        path(gwas_tmp_dir)
        path(vcf_in)
        path(vcf_tbi_in)
        path(similarity_in)
        path(pastml_steps)
        file(plink_bed)
        file(plink_bim)
        file(plink_fam)

    output:
        val true
        // file "*.pyseer.phen"
        // file "*_*.pyseer_output.csv"

    errorStrategy 'ignore'

    script:
    """
    bash "$gwas_bash_script"
    """

}

process run_gwas_bash_script_phensim_burden {

    input:
        file(gwas_bash_script)
        path(gwas_tmp_dir)
        path(vcf_in)
        path(vcf_tbi_in)
        path(regions_in)
        path(similarity_in)
        path(pastml_steps)
        file(plink_bed)
        file(plink_bim)
        file(plink_fam)

    output:
        val true
        // file "*.pyseer.phen"
        // file "*_*.pyseer_output.csv"

    errorStrategy 'ignore'

    script:
    """
    bash "$gwas_bash_script"
    """

}

process run_gwas_bash_script_subsampling_burden {

    input:
        file(gwas_bash_script)
        path(gwas_tmp_dir)
        path(vcf_in)
        path(vcf_tbi_in)
        path(regions_in)
        path(similarity_in)
        path(causal_in)
        path(phenotype_in)

    output:
        val true
        // file "*.pyseer.phen"
        // file "*_*.pyseer_output.csv"

    errorStrategy 'ignore'

    script:
    """
    bash "$gwas_bash_script"
    """

}


process run_gwas_bash_script_subsampling_vcf {

    input:
        file(gwas_bash_script)
        path(gwas_tmp_dir)
        path(vcf_in)
        path(vcf_tbi_in)
        path(similarity_in)
        path(causal_in)
        path(phenotype_in)

    output:
        val true
        // file "*.pyseer.phen"
        // file "*_*.pyseer_output.csv"

    errorStrategy 'ignore'

    script:
    """
    bash "$gwas_bash_script"
    """

}




process process_gwas_runs_subsampling {

    input:
        val(flag)
        path(gwas_runs_in_table)
        path(gwas_tmp_dir)
        val(variant_type)
        path(causal_variants)
        val(outprefix)
        path(tmp_dir)

    publishDir "${tmp_dir}"

    output:
        file "${outprefix}.gwas_runs.results.csv"

    script:
    """
    process_gwas_runs.py --gwas_runs_in_table "$gwas_runs_in_table" --gwas_runs_out_table "$outprefix".gwas_runs.results.csv --output_dir "$gwas_tmp_dir" --variant_type "$variant_type" --causal-loci "$causal_variants"
    """
}


process process_gwas_runs {

    input:
        val(flag)
        path(gwas_runs_in_table)
        path(gwas_tmp_dir)
        val(variant_type)
        val(outprefix)
        path(tmp_dir)

    publishDir "${tmp_dir}"

    output:
        file "${outprefix}.gwas_runs.results.csv"

    script:
    """
    process_gwas_runs.py --gwas_runs_in_table "$gwas_runs_in_table" --gwas_runs_out_table "$outprefix".gwas_runs.results.csv --output_dir "$gwas_tmp_dir" --variant_type "$variant_type"
    """
}



// Process to calculate number of unique pan-genome patterns from Roary .Rtab file
process pan_genome_patterns {

    input:
        path(pg_in)

    output:
        stdout

    script:
    """
    cat "$pg_in" | cut -f2- | sed 1d | sort | uniq | wc -l
    """

}


// Process to calculate number of unique VCF patterns
process vcf_patterns {

    input:
        path(vcf_in)
        path(vcf_tbi_in)


    output:
        stdout

    script:
    """
    bcftools view "$vcf_in" | bcftools query -f '[%GT]\n' | sort | uniq | wc -l
    """

}

process number_of_genes {

    input:
        path(regions_in)

    output:
        stdout

    script:
    """
    cat "regions_in" | wc -l
    """
}

// NOTE: plot_gwas_runs_subsampling plots both average p-value and power

process plot_gwas_runs_subsampling {

    input:
        path(gwas_runs_out_table)
        file(parameters_file)
        val(variants_tested)
        val(significance_level)
        val(outprefix)
        path(tmp_dir)
        val(pyseer_model)


    publishDir "${tmp_dir}"

    output:
        val true
        file "${outprefix}.subsampling.power.pdf"
        file "${outprefix}.subsampling.avg_pval.pdf"

    script:
    """
    plot_gwas_runs_subsampling.R --input_table "$gwas_runs_out_table" --parameters_file "$parameters_file" --variants_tested "$variants_tested" --significance_level "$significance_level" --output_plot "$outprefix".subsampling.power.pdf --pyseer_model "$pyseer_model" --plot_metrics 1
    plot_gwas_runs_subsampling.R --input_table "$gwas_runs_out_table" --parameters_file "$parameters_file" --variants_tested "$variants_tested" --significance_level "$significance_level" --output_plot "$outprefix".subsampling.avg_pval.pdf --pyseer_model "$pyseer_model" --plot_metrics 2
    """
}



process plot_gwas_runs {

    input:
        path(gwas_runs_out_table)
        file(parameters_file)
        val(variants_tested)
        val(significance_level)
        val(plot_type)
        val(outprefix)
        path(tmp_dir)
        val(pyseer_model)


    publishDir "${tmp_dir}"

    output:
        val true
        file "${outprefix}.phensim.plot_type${plot_type}.power.pdf"
        file "${outprefix}.phensim.plot_type${plot_type}.avg_pval.pdf"

    script:
    """
    plot_gwas_runs.R --input_table "$gwas_runs_out_table" --parameters_file "$parameters_file" --variants_tested "$variants_tested" --significance_level "$significance_level" --output_plot "$outprefix".phensim.plot_type"$plot_type".power.pdf --plot_type "$plot_type" --pyseer_model "$pyseer_model" --plot_metrics 1
    plot_gwas_runs.R --input_table "$gwas_runs_out_table" --parameters_file "$parameters_file" --variants_tested "$variants_tested" --significance_level "$significance_level" --output_plot "$outprefix".phensim.plot_type"$plot_type".avg_pval.pdf --plot_type "$plot_type" --pyseer_model "$pyseer_model" --plot_metrics 2
    """
}





