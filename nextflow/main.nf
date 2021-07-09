/*
 * Nextflow pipeline for PowerBacGWAS: Power calculations for Bacterial GWAS
 *
 */

// Enable DSL 2
nextflow.enable.dsl=2


// Import modules
include {printHelp} from './modules/help.nf'
include {prepare_vcf_file; vcf_to_plink_files; vcf_to_pastml_matrix; roary_to_plink_files; roary_to_pastml_matrix; annotate_nodes} from './modules/preprocessing.nf'
include {ancestral_state_reconstruction_vcf; ancestral_state_reconstruction_pg; calculate_changes_per_region} from './modules/asr.nf'
include {phylogeny_distance; prepare_gwas_runs_subsampling_pg; collect_gwas_scripts; run_gwas_bash_script_subsampling_pg; process_gwas_runs_subsampling; pan_genome_patterns; plot_gwas_runs_subsampling} from './modules/gwas.nf'
include {prepare_gwas_runs_subsampling_burden; run_gwas_bash_script_subsampling_burden; prepare_gwas_runs_subsampling_vcf; run_gwas_bash_script_subsampling_vcf} from './modules/gwas.nf'
include {prepare_gwas_runs_phensim_pg; run_gwas_bash_script_phensim_pg; process_gwas_runs; plot_gwas_runs; prepare_gwas_runs_phensim_vcf; run_gwas_bash_script_phensim_vcf; vcf_patterns; prepare_gwas_runs_phensim_burden; run_gwas_bash_script_phensim_burden; number_of_genes } from './modules/gwas.nf'
include {copy_files} from './modules/other.nf'




// Help message
if (params.help){
    printHelp()
    exit 0
}

// Make sure required arguments are specified
if (params.outprefix == ""){
    println("Please specify an output prefix with --outprefix.")
    println("Print help with --help")
    System.exit(1)
}

// Create directories if they don't already exist
// NOTE: default params values are defined in nextflow.config
tmp_dir = file(params.tmp_dir)
tmp_dir.mkdir()

gwas_tmp_dir = file(params.gwas_tmp_dir)
gwas_tmp_dir.mkdir()

results_dir = file(params.results_dir)
results_dir.mkdir()

// NOTE: the code below makes sure the right input files and/or intermediate files are present before running user-selected workflows

if (params.tree){
    Channel.fromPath( params.tree, checkIfExists: true )
        .set { tree_ch }
}

if (params.run_prepare_vcf){
    if (params.vcf == ""){
        println("Please specify a VCF file with --vcf.")
        println("Print help with --help")
        System.exit(1)
    }
}

if (params.vcf){
    // Create read pairs channel
    Channel.fromPath( params.vcf, checkIfExists: true )
        .set { vcf_ch }
}

if (params.run_prepare_pangenome){
    if (params.pangenome == ""){
        println("Please specify a Pan-genome file (Roary's gene_presense_absence.Rtab) with --pangenome.")
        println("Print help with --help")
        System.exit(1)
    }
}

if (params.pangenome){
    // Create read pairs channel
    Channel.fromPath( params.pangenome, checkIfExists: true )
        .set { pg_ch }
}

if (params.run_vcf_regions_ast | params.run_subsampling_burden){
    if (params.regions == ""){
        println("Please specify a regions file with --regions.")
        println("Print help with --help")
        System.exit(1)
    }
}

if (params.run_subsampling_burden){
    if (params.causal_regions == ""){
        println("Please specify a causal regions file with --causal_regions.")
        println("Print help with --help")
        System.exit(1)
    }
}



if (params.run_vcf_ast | params.run_pangenome_ast | params.run_vcf_regions_ast){
    if (params.tree == ""){
        println("Please specify a tree in newick format with --tree to run an ancestral state reconstruction pipeline")
        println("Print help with --help")
        System.exit(1)
    }
}


if (params.regions){
    // Create read pairs channel
    Channel.fromPath( params.regions, checkIfExists: true )
        .set { regions_ch }
}

if (params.run_subsampling_pangenome | params.run_phensim_pangenome){
    if (params.pangenome == ""){
        println("Please specify a pan-genome file with --pangenome to run --run_subsampling_pangenome or --run_phensim_pangenome")
        println("Print help with --help")
        System.exit(1)
    }
}

if (params.run_subsampling_vcf | params.run_phensim_vcf){
    if (params.vcf== ""){
        println("Please specify a VCF file with --vcf to run --run_subsampling_vcf or --run_phensim_vcf.")
        println("Print help with --help")
        System.exit(1)
    }
}

if (params.run_subsampling_burden | params.run_phensim_burden){
    if (params.vcf== ""){
        println("Please specify a VCF file with --vcf to run --run_subsampling_burden or --run_phensim_burden.")
        println("Print help with --help")
        System.exit(1)
    }
    if (params.regions== ""){
        println("Please specify a regions file with --regions to run --run_subsampling_burden or --run_phensim_burden.")
        println("Print help with --help")
        System.exit(1)
    }
}

if (params.run_subsampling_pangenome | params.run_subsampling_vcf | params.run_subsampling_burden){
    if (params.causal_variants == ""){
        println("Please specify a causal variants file with --causal_variants to run a sub-sampling pipeline")
        println("Print help with --help")
        System.exit(1)
    }
    if (params.phenotype == ""){
        println("Please specify a phenotype file with --phenotype to run a sub-sampling pipeline")
        println("Print help with --help")
        System.exit(1)
    }
    if (params.parameters == ""){
        println("Please specify a parameters file with --parameters to run a sub-sampling pipeline")
        println("Print help with --help")
        System.exit(1)
    }
    if (params.tree == ""){
        println("Please specify a tree in newick format --tree to run a sub-sampling pipeline")
        println("Print help with --help")
        System.exit(1)
    }
}

if (params.run_phensim_pangenome | params.run_phensim_vcf | params.run_phensim_burden){
    if (params.parameters == ""){
        println("Please specify a parameters file with --parameters to run a phenotype simulation pipeline")
        println("Print help with --help")
        System.exit(1)
    }
    if (params.tree == ""){
        println("Please specify a tree in newick format --tree to run a phenotype simulation pipeline")
        println("Print help with --help")
        System.exit(1)
    }
}

// Check Pyseer model selected is correct
if (params.pyseer_model){
    if (params.pyseer_model != "LMM"){
        if (params.pyseer_model != "ENET"){
            println("Please select a valid Pyseer model with --pyseer_model: 'LMM' or 'ENET'")
            println("Print help with --help")
            System.exit(1)
	}
    }
}
// Check phenotype type selected is correct
if (params.phenotype_type){
    if (params.phenotype_type != "binary"){
        if (params.phenotype_type != "quantitative"){
            println("Please select a valid phenotype type with --phenotype_type: either 'binary' or 'quantitative'")
            println("Print help with --help")
            System.exit(1)
	}
    }
}




// If --num_cpus specified
println "${params.num_cores} cores will be used for multiprocessing tasks"


// Workflow to prepare VCF files
workflow PREPARE_VCF {

    take:
        vcf_ch

    main:
        prepare_vcf_file(vcf_ch, params.outprefix, tmp_dir)

        vcf_to_plink_files(prepare_vcf_file.out[0], params.outprefix, tmp_dir)

	vcf_to_pastml_matrix(prepare_vcf_file.out[0], params.outprefix, tmp_dir)

    emit:
	// not sure what process output to emit here --> to check
        vcf_to_pastml_matrix.out
}

// Workflow to prepare pan-genome file
workflow PREPARE_PG {

    take:
        pg_ch

    main:
        roary_to_plink_files(pg_ch, params.outprefix, tmp_dir)

	roary_to_pastml_matrix(pg_ch, params.outprefix, tmp_dir)

    emit:
	// not sure what process output to emit here --> to check
        roary_to_pastml_matrix.out
}


// Main Workflow
workflow {

    main:

        if (params.run_prepare_vcf){
            PREPARE_VCF(vcf_ch)
	}

        if (params.run_prepare_pangenome){
            PREPARE_PG(pg_ch)
	}

	// The internal nodes of the phylogenetic tree need to be annotated before AST
        if (params.run_vcf_ast | params.run_pangenome_ast | params.run_vcf_regions_ast){
		annotate_nodes(tree_ch, params.outprefix, tmp_dir)
	}

        if (params.run_vcf_ast){

    		vcf_pastml_input = file("${tmp_dir}/${params.outprefix}.vcf.pastml.csv")

    		if( ! vcf_pastml_input.isFile() ){
	  		println "File ${params.outprefix}.vcf.pastml.csv not found. Run --run_prepare_vcf before running --run_vcf_ast"
           		System.exit(1)
    		}

            	ancestral_state_reconstruction_vcf(vcf_pastml_input, annotate_nodes.out, params.outprefix, tmp_dir, params.num_cores)
	}

	if (params.run_vcf_regions_ast){

    		vcf_anc_input = file("${tmp_dir}/${params.outprefix}.vcf.ancestral.csv")
    		if( ! vcf_anc_input.isFile() ){
	   		println "File ${params.outprefix}.vcf.ancestral.csv not found. Run --run_vcf_ast before running --run_vcf_regions_ast"
           		System.exit(1)
    		}

         	vcf_anc_input = file("${tmp_dir}/${params.outprefix}.vcf.ancestral.csv")

		calculate_changes_per_region(vcf_anc_input, annotate_nodes.out, regions_ch, params.outprefix, tmp_dir, params.num_cores)

	}

        if (params.run_pangenome_ast){

    		pg_pastml_input = file("${tmp_dir}/${params.outprefix}.pastml.csv")
    		if( ! pg_pastml_input.isFile() ){
	   		println "File ${params.outprefix}.pastml.csv not found. Run --run_prepare_pangenome before running --run_pangenome_ast"
           		System.exit(1)
    		}

            	ancestral_state_reconstruction_pg(pg_pastml_input, annotate_nodes.out, params.outprefix, tmp_dir, params.num_cores)
	}

	if (params.run_subsampling_vcf | params.run_subsampling_burden | params.run_phensim_vcf | params.run_phensim_burden){
		vcf_input = file("${tmp_dir}/${params.outprefix}.reformatted.vcf.gz")
		if( ! vcf_input.isFile() ){
	   		println "File ${params.outprefix}.reformatted.vcf.gz not found. Run --run_prepare_vcf"
           		System.exit(1)
    		}
	}

	// to do: update when Elastic Net option added
	if (params.run_subsampling_pangenome | params.run_subsampling_vcf | params.run_subsampling_burden | params.run_phensim_pangenome | params.run_phensim_vcf | params.run_phensim_burden  ){
		phylogeny_distance(tree_ch, params.outprefix, tmp_dir)
	}

	if (params.run_subsampling_pangenome){

		// NOTE: it is very important that only channel collect_gwas_scripts.out is passed into run_gwas_bash_script() process to allow processing all bash scripts
		// That's why all other arguments are passed as file paths

        		parameters = file(params.parameters, checkIfExists: true)
        		causal_variants = file(params.causal_variants, checkIfExists: true)
        		phenotype = file(params.phenotype, checkIfExists: true)
        		pangenome = file(params.pangenome, checkIfExists: true)

		prepare_gwas_runs_subsampling_pg(pg_ch, phylogeny_distance.out, parameters, causal_variants, phenotype, params.outprefix, tmp_dir, gwas_tmp_dir, params.pyseer_model, params.num_cores)

		collect_gwas_scripts(prepare_gwas_runs_subsampling_pg.out[0], gwas_tmp_dir)

		distance = file("${tmp_dir}/${params.outprefix}.tree_distances.csv")
		run_gwas_bash_script_subsampling_pg(collect_gwas_scripts.out.flatten(), gwas_tmp_dir, pangenome, distance, causal_variants, phenotype)

		gwas_runs_in_table = file("${tmp_dir}/${params.outprefix}.gwas_runs.csv")
		process_gwas_runs_subsampling(run_gwas_bash_script_subsampling_pg.out.last(), gwas_runs_in_table, gwas_tmp_dir, "r", causal_variants, params.outprefix, tmp_dir)

		pan_genome_patterns(pg_ch)

		plot_gwas_runs_subsampling(process_gwas_runs_subsampling.out, parameters, pan_genome_patterns.out, params.siglev, params.outprefix, tmp_dir, params.pyseer_model)

		final_files = [ file("${tmp_dir}/${params.outprefix}.gwas_runs.results.csv"), file("${tmp_dir}/${params.outprefix}.subsampling.power.pdf"), file("${tmp_dir}/${params.outprefix}.subsampling.avg_pval.pdf")]
		copy_files(plot_gwas_runs_subsampling.out[0], final_files, results_dir)
	}

	if (params.run_subsampling_vcf){
		
		parameters = file(params.parameters, checkIfExists: true)
        		causal_variants = file(params.causal_variants, checkIfExists: true)
        		phenotype = file(params.phenotype, checkIfExists: true)
        		vcf_input = file("${tmp_dir}/${params.outprefix}.reformatted.vcf.gz")
        		vcf_input_tbi = file("${tmp_dir}/${params.outprefix}.reformatted.vcf.gz.tbi")

		prepare_gwas_runs_subsampling_vcf(vcf_input, vcf_input_tbi, phylogeny_distance.out, parameters, causal_variants, phenotype, params.outprefix, tmp_dir, gwas_tmp_dir, params.pyseer_model, params.num_cores)

		collect_gwas_scripts(prepare_gwas_runs_subsampling_vcf.out[0], gwas_tmp_dir)

		distance = file("${tmp_dir}/${params.outprefix}.tree_distances.csv")
		run_gwas_bash_script_subsampling_vcf(collect_gwas_scripts.out.flatten(), gwas_tmp_dir, vcf_input, vcf_input_tbi, distance, causal_variants, phenotype)

		gwas_runs_in_table = file("${tmp_dir}/${params.outprefix}.gwas_runs.csv")
		process_gwas_runs_subsampling(run_gwas_bash_script_subsampling_vcf.out.last(), gwas_runs_in_table, gwas_tmp_dir, "v", causal_variants, params.outprefix, tmp_dir)

		vcf_patterns(vcf_input, vcf_input_tbi)

		plot_gwas_runs_subsampling(process_gwas_runs_subsampling.out, parameters, vcf_patterns.out, params.siglev, params.outprefix, tmp_dir, params.pyseer_model)

		final_files = [ file("${tmp_dir}/${params.outprefix}.gwas_runs.results.csv"), file("${tmp_dir}/${params.outprefix}.subsampling.power.pdf"), file("${tmp_dir}/${params.outprefix}.subsampling.avg_pval.pdf")]
		copy_files(plot_gwas_runs_subsampling.out[0], final_files, results_dir)
	}

	if (params.run_subsampling_burden){

                 parameters = file(params.parameters, checkIfExists: true)
        		causal_variants = file(params.causal_variants, checkIfExists: true)
		causal_regions = file(params.causal_regions, checkIfExists: true)
        		phenotype = file(params.phenotype, checkIfExists: true)
        		regions = file(params.regions, checkIfExists: true)
        		vcf_input = file("${tmp_dir}/${params.outprefix}.reformatted.vcf.gz")
        		vcf_input_tbi = file("${tmp_dir}/${params.outprefix}.reformatted.vcf.gz.tbi")

		prepare_gwas_runs_subsampling_burden(vcf_input, vcf_input_tbi, phylogeny_distance.out, parameters, causal_variants, phenotype, regions_ch, params.outprefix, tmp_dir, gwas_tmp_dir, params.pyseer_model, params.num_cores)

		collect_gwas_scripts(prepare_gwas_runs_subsampling_burden.out[0], gwas_tmp_dir)

		distance = file("${tmp_dir}/${params.outprefix}.tree_distances.csv")
		run_gwas_bash_script_subsampling_burden(collect_gwas_scripts.out.flatten(), gwas_tmp_dir, vcf_input, vcf_input_tbi, regions, distance, causal_variants, phenotype)

		gwas_runs_in_table = file("${tmp_dir}/${params.outprefix}.gwas_runs.csv")
		process_gwas_runs_subsampling(run_gwas_bash_script_subsampling_burden.out.last(), gwas_runs_in_table, gwas_tmp_dir, "b", causal_regions, params.outprefix, tmp_dir)

		number_of_genes(regions)

		plot_gwas_runs_subsampling(process_gwas_runs_subsampling.out, parameters, number_of_genes.out, params.siglev, params.outprefix, tmp_dir, params.pyseer_model)

		final_files = [ file("${tmp_dir}/${params.outprefix}.gwas_runs.results.csv"), file("${tmp_dir}/${params.outprefix}.subsampling.power.pdf"), file("${tmp_dir}/${params.outprefix}.subsampling.avg_pval.pdf")]
		copy_files(plot_gwas_runs_subsampling.out[0], final_files, results_dir)
	}


	if (params.run_phensim_pangenome){
		pastml_steps = file("${tmp_dir}/${params.outprefix}.ancestral_steps.csv", checkIfExists: true)
		if( ! pastml_steps.isFile() ){
	   		println "File ${tmp_dir}/${params.outprefix}.ancestral_steps.csv not found. Run ---run_pangenome_ast first"
           		System.exit(1)
    		}
	}

	if (params.run_phensim_vcf){
		pastml_steps = file("${tmp_dir}/${params.outprefix}.vcf.ancestral_steps.csv", checkIfExists: true)
		if( ! pastml_steps.isFile() ){
	   		println "File ${tmp_dir}/${params.outprefix}.vcf.ancestral_steps.csv not found. Run ---run_vcf_ast first"
           		System.exit(1)
    		}
	}

	if (params.run_phensim_burden){
		pastml_steps = file("${tmp_dir}/${params.outprefix}.regions.ancestral_steps.csv", checkIfExists: true)
		if( ! pastml_steps.isFile() ){
	   		println "File ${tmp_dir}/${params.outprefix}.regions.ancestral_steps.csv not found. Run ---run_vcf_regions_ast first"
           		System.exit(1)
    		}
	}


	if (params.run_phensim_pangenome){

		// NOTE: it is very important that only channel collect_gwas_scripts.out is passed into run_gwas_bash_script() process to allow processing all bash scripts
		// That's why all other arguments are passed as file paths

        		parameters = file(params.parameters, checkIfExists: true)
        		pangenome = file(params.pangenome, checkIfExists: true)
		pastml_steps = file("${tmp_dir}/${params.outprefix}.ancestral_steps.csv", checkIfExists: true)
		plink_bed = file("${tmp_dir}/${params.outprefix}.bed", checkIfExists: true)
		plink_bim = file("${tmp_dir}/${params.outprefix}.bim", checkIfExists: true)
		plink_fam = file("${tmp_dir}/${params.outprefix}.fam", checkIfExists: true)
		
		prepare_gwas_runs_phensim_pg(pg_ch, phylogeny_distance.out, parameters, pastml_steps, params.outprefix, plink_bed, plink_bim, plink_fam, tmp_dir, gwas_tmp_dir, params.pyseer_model, params.phenotype_type, params.num_cores)

		collect_gwas_scripts(prepare_gwas_runs_phensim_pg.out[0], gwas_tmp_dir)

		distance = file("${tmp_dir}/${params.outprefix}.tree_distances.csv")
		run_gwas_bash_script_phensim_pg(collect_gwas_scripts.out.flatten(), gwas_tmp_dir, pangenome, distance, pastml_steps, plink_bed, plink_bim, plink_fam)

		gwas_runs_in_table = file("${tmp_dir}/${params.outprefix}.gwas_runs.csv")
		process_gwas_runs(run_gwas_bash_script_phensim_pg.out.last(), gwas_runs_in_table, gwas_tmp_dir, "r", params.outprefix, tmp_dir)

		pan_genome_patterns(pg_ch)

		// to do: add plot_type option
                 plot_type = "1"
		plot_gwas_runs(process_gwas_runs.out, parameters, pan_genome_patterns.out, params.siglev, plot_type, params.outprefix, tmp_dir, params.pyseer_model)

		final_files = [ file("${tmp_dir}/${params.outprefix}.gwas_runs.results.csv"), file("${tmp_dir}/${params.outprefix}.phensim.plot_type${plot_type}.power.pdf"), file("${tmp_dir}/${params.outprefix}.phensim.plot_type${plot_type}.avg_pval.pdf")]
		copy_files(plot_gwas_runs.out[0], final_files, results_dir)
	}


	if (params.run_phensim_vcf){

		// NOTE: it is very important that only channel collect_gwas_scripts.out is passed into run_gwas_bash_script() process to allow processing all bash scripts
		// That's why all other arguments are passed as file paths

        		parameters = file(params.parameters, checkIfExists: true)
		vcf_input = file("${tmp_dir}/${params.outprefix}.reformatted.vcf.gz", checkIfExists: true)
        		vcf_input_tbi = file("${tmp_dir}/${params.outprefix}.reformatted.vcf.gz.tbi", checkIfExists: true)
		pastml_steps = file("${tmp_dir}/${params.outprefix}.vcf.ancestral_steps.csv", checkIfExists: true)
		plink_bed = file("${tmp_dir}/${params.outprefix}.bed", checkIfExists: true)
		plink_bim = file("${tmp_dir}/${params.outprefix}.bim", checkIfExists: true)
		plink_fam = file("${tmp_dir}/${params.outprefix}.fam", checkIfExists: true)

		prepare_gwas_runs_phensim_vcf(vcf_input, vcf_input_tbi, phylogeny_distance.out, parameters, pastml_steps, params.outprefix, plink_bed, plink_bim, plink_fam, tmp_dir, gwas_tmp_dir, params.pyseer_model, params.phenotype_type, params.num_cores)

		collect_gwas_scripts(prepare_gwas_runs_phensim_vcf.out[0], gwas_tmp_dir)

		distance = file("${tmp_dir}/${params.outprefix}.tree_distances.csv")
		run_gwas_bash_script_phensim_vcf(collect_gwas_scripts.out.flatten(), gwas_tmp_dir, vcf_input, vcf_input_tbi, distance, pastml_steps, plink_bed, plink_bim, plink_fam)

		gwas_runs_in_table = file("${tmp_dir}/${params.outprefix}.gwas_runs.csv")
		process_gwas_runs(run_gwas_bash_script_phensim_vcf.out.last(), gwas_runs_in_table, gwas_tmp_dir, "v", params.outprefix, tmp_dir)

		vcf_patterns(vcf_input, vcf_input_tbi)

		// to do: add plot_type option
                 plot_type = "1"
		plot_gwas_runs(process_gwas_runs.out, parameters, vcf_patterns.out, params.siglev, plot_type, params.outprefix, tmp_dir, params.pyseer_model)

		final_files = [ file("${tmp_dir}/${params.outprefix}.gwas_runs.results.csv"), file("${tmp_dir}/${params.outprefix}.phensim.plot_type${plot_type}.power.pdf"), file("${tmp_dir}/${params.outprefix}.phensim.plot_type${plot_type}.avg_pval.pdf")]
		copy_files(plot_gwas_runs.out[0], final_files, results_dir)
	}

	if (params.run_phensim_burden){

		// NOTE: it is very important that only channel collect_gwas_scripts.out is passed into run_gwas_bash_script() process to allow processing all bash scripts
		// That's why all other arguments are passed as file paths

        		parameters = file(params.parameters, checkIfExists: true)
		vcf_input = file("${tmp_dir}/${params.outprefix}.reformatted.vcf.gz", checkIfExists: true)
        		vcf_input_tbi = file("${tmp_dir}/${params.outprefix}.reformatted.vcf.gz.tbi", checkIfExists: true)
		pastml_steps = file("${tmp_dir}/${params.outprefix}.regions.ancestral_steps.csv", checkIfExists: true)
		plink_bed = file("${tmp_dir}/${params.outprefix}.bed", checkIfExists: true)
		plink_bim = file("${tmp_dir}/${params.outprefix}.bim", checkIfExists: true)
		plink_fam = file("${tmp_dir}/${params.outprefix}.fam", checkIfExists: true)
		regions = file(params.regions, checkIfExists: true)

		prepare_gwas_runs_phensim_burden(vcf_input, vcf_input_tbi, phylogeny_distance.out, parameters, pastml_steps, params.outprefix, plink_bed, plink_bim, plink_fam, regions, tmp_dir, gwas_tmp_dir, params.pyseer_model, params.phenotype_type, params.num_cores)

		collect_gwas_scripts(prepare_gwas_runs_phensim_burden.out[0], gwas_tmp_dir)

		distance = file("${tmp_dir}/${params.outprefix}.tree_distances.csv")
		run_gwas_bash_script_phensim_burden(collect_gwas_scripts.out.flatten(), gwas_tmp_dir, vcf_input, vcf_input_tbi, regions, distance, pastml_steps, plink_bed, plink_bim, plink_fam)

		gwas_runs_in_table = file("${tmp_dir}/${params.outprefix}.gwas_runs.csv")
		process_gwas_runs(run_gwas_bash_script_phensim_burden.out.last(), gwas_runs_in_table, gwas_tmp_dir, "b", params.outprefix, tmp_dir)

		number_of_genes(regions)

		// to do: add plot_type option
                 plot_type = "1"
		plot_gwas_runs(process_gwas_runs.out, parameters, number_of_genes.out, params.siglev, plot_type, params.outprefix, tmp_dir, params.pyseer_model)

		final_files = [ file("${tmp_dir}/${params.outprefix}.gwas_runs.results.csv"), file("${tmp_dir}/${params.outprefix}.phensim.plot_type${plot_type}.power.pdf"), file("${tmp_dir}/${params.outprefix}.phensim.plot_type${plot_type}.avg_pval.pdf")]
		copy_files(plot_gwas_runs.out[0], final_files, results_dir)
	}






}


