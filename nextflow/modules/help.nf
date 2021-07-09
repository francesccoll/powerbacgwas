
def printHelp() {
  log.info"""
  Usage:
    nextflow run main.nf --tree [phylogenetic tree] --outprefix [output prefix]

  Description:
    PowerBacGWAS: Power calculations for Bacterial GWAS

  Additional options:

    Inputs:
      --vcf                         multi-sample VCF file needed to run variant or burden test GWAS
      --pangenome                   pan-genome file in Roary's gene_presence_absence.Rtab format
      --regions                     The regions file can be used to specify gene coordinates. It must contain one region per line, with their name and the bcftools style region co-ordinates 
                                    delimited by tab (locus_name chr_name:start-end).

    Outputs:
      --results_dir                 Results directory for output files. (Default: './results')
      --tmp_dir                     Temporary directory for intermediate files. (Default: './tmp')

    Pipeline Options:
      --run_prepare_vcf             Run the pre-processings steps needed to prepare the VCF for downstream analyses. Must also specify
                                    --vcf input. (Default: false)
      --run_prepare_pangenome       Run the pre-processings steps needed to prepare Roary's gene_presence_absence.Rtab file for downstream analyses.
                                    Must also specify --pangenome input. (Default: false)
      --run_vcf_ast                 Run ancestral state reconstruction for VCF. --run_prepare_vcf must be run first. Must also specify
                                    --vcf input. (Default: false)
      --run_vcf_regions_ast         Calculate number of independent variant acquisitions (homoplasies) at the region level (e.g. gene). Needed for burden testing power calculations.
                                    Must also specify --vcf and --regions inputs.
      --run_pangenome_ast           Run ancestral state reconstruction for pan-genome. --run_prepare_pangenome must be run first. Must also specify
                                    --pangenome input. (Default: false)
      --run_subsampling_vcf         Run the sub-sampling approach for individual VCF variants. Must also specify --pangenome and --tree inputs.
      --run_subsampling_pangenome   Run the sub-sampling approach for the pan-genome. Must also specify --vcf and --tree inputs.
      --run_subsampling_burden	  Run the sub-sampling approach for burden testing. Must also specify --vcf, --tree and --regions inputs.
      --run_phensim_vcf             Run the phenotype simulation approach for individual VCF variants. Must also specify --pangenome and --tree inputs.
      --run_phensim_pangenome       Run the phenotype simulation approach for the pan-genome. Must also specify --vcf and --tree inputs.
      --run_phensim_burden	  Run the phenotype simulation approach for burden testing. Must also specify --vcf, --tree and --regions inputs.


    Pyseer options:
      --pyseer_model                Pyseer model to run. To choose between: "LMM" (Linear Mixed Models) or "ENET" (Elastic Net) (Default: LMM)
      --phenotype_type              Type of phenotype to simulate. To choose between: "binary" or "quantitative" (Default: binary)


    Required input files for sub-sampling pipelines (when running --run_subsampling_vcf, --run_subsampling_pangenome or --run_subsampling_burden):
      --causal_variants             A file with causal variants (a GCTA-formatted file with a list of variants IDs)
      --phenotype                   A Pyseer-formatted phenotype file, where the variants in the causal variants file are known to determine the trait in the phenotype file. Only binary phenotypes supported.
      --parameters                  A parameters file specifying the range of allele frequencies, sample sizes and effect sizes to test.

    Required input files for phenotype simulation pipelines (when running --run_phensim_vcf, --run_phensim_pangenome or --run_phensim_burden):
      --parameters                  A parameters file specifying the range of allele frequencies, sample sizes, effect sizes, heritability and homoplasy level to test.


    Process running Options:
      --num_cores                   Number of CPUs to be used for multiprocessing tasks (Default: 8).
                                    NOTE: if using Docker Desktop, make sure CPUs in Resources is equal or greater than --num_cores



  """.stripIndent()
}
