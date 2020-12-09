#!/usr/bin/env Rscript

# This R script is used to plot GWAS power calculations results produced by script process_gwas_runs.py

#============================================================================================================
#						Making sure required R packages are installed
#============================================================================================================

if(!"ggplot2" %in% installed.packages()[,"Package"]) {
  cat("ERROR: Package ggplot2 could not be found. Run the following command to install it: install.packages(\"ggplot2\")\n"); quit();
}
if(!"optparse" %in% installed.packages()[,"Package"]) {
  cat("ERROR: Package optparse could not be found. Run the following command to install it: install.packages(\"optparse\")\n"); quit();
}

#============================================================================================================
#						R functions
#============================================================================================================

check_file_exists = function(file, file_name)
{
  if(file == ""){
    print(paste("ERROR: required input file ", file_name, " not specified.", sep = "")); quit();
  } else {
    if (!file.exists(file)) {
      print(paste("ERROR: required input ", file_name, " ", file, " could not be found.", sep = "")); quit();
    }
  }
}

check_vector_not_empty = function(vector, vector_name)
{
  if(length(vector)==0){
    print(paste("ERROR: vector ", vector_name, " is empty. Check --parameters_file format", sep = "")); quit();
  }
}

is_integer = function(x)
{
  is_int = FALSE
  if(is.numeric(x)){
    if(x == round(x)){
      is_int = TRUE
    }
  }
  return(is_int)
}

#============================================================================================================
#						Parsing script arguments
#============================================================================================================

library(optparse)

option_list = list(
  make_option(c("-i", "--input_table"), type="character", default="", 
              help="Table with GWAS runs results created by script process_gwas_runs.py", metavar="character"),
  make_option(c("-p", "--parameters_file"), type="character", default="", 
              help="File with GWAS power calculations parameters used by prepare_gwas_runs.py", metavar="character"),
  make_option(c("-t", "--plot_type"), type="integer", default=1, 
              help="Select 1, 2, 3 and 4 to choose type of plot: 1. Across sample sizes and effect sizes. 2. Across sample sizes and heritability values. 3. Across sample sizes and homoplasy levels (fixed allele frequency)", metavar="integer"),
  make_option(c("-f", "--allele_frequency"), type="character", default=NA, 
              help="Select allele frequency range to keep (e.g. 0.05,0.10) from --input_table. Required for plot type 3", metavar="character"),
  make_option(c("-a", "--average_pval_across"), type="character", default="both", 
              help="Choose to average p-values across: variant 'sampling' replicates, 'phenotype' simulation replicates or 'both' (default)", metavar="character"),
  make_option(c("-v", "--variants_tested"), type="integer", default=NA, 
              help="Number of variants/genes tested in GWAS. Used to calculate and plot genome-wide significance level.", metavar="integer"),
  make_option(c("-l", "--significance_level"), type="double", default=0.05, 
              help="Significance level. Used to calculate and plot genome-wide significance level.", metavar="double"),
  make_option(c("-e", "--add_error_bars"), type="character", default="no", 
              help="whether to add -log10(p-value) error bars. To choose from 'yes' or 'no'.", metavar="character"),
  make_option(c("-o", "--output_plot"), type="character", default="out.pdf", 
              help="output file name [default= %default]", metavar="character"),
  make_option(c("-x", "--plot_width"), type="integer", default="21", 
              help="Plot width in cm", metavar="integer"),
  make_option(c("-y", "--plot_height"), type="integer", default="10", 
              help="Plot height in cm", metavar="integer"),
  make_option(c("-z", "--ylim_to"), type="integer", default=NA, 
              help="Maximum significance level (in -log10(p-value) units) to show on plot", metavar="integer")
  );

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);


#============================================================================================================
#						Checks related to user-defined options
#============================================================================================================

if(opt$plot_type != 1 & opt$plot_type != 2 & opt$plot_type != 3 & opt$plot_type != 4){
  cat("ERROR: plot_type option provided (", opt$plot_type, ") must be either 1, 2, 3 or 4\n"); quit();
}
if(opt$average_pval_across != "sampling" & opt$average_pval_across != "phenotype" & opt$average_pval_across != "both"){
  cat("ERROR: average_pval_across option provided (", opt$average_pval_across, ") must be either 'sampling', 'phenotype' or 'both'\n"); quit();
}
if(opt$add_error_bars != "yes" & opt$add_error_bars != "no"){
  cat("ERROR: add_error_bars option provided (", opt$add_error_bars, ") must be either 'yes' or 'no'\n"); quit();
}
if(opt$plot_type == 3){
  if(is.na(opt$allele_frequency)){
    cat("ERROR: Option --allele_frequency must be specified if --plot_type 3 is to be created\n"); quit();
  }
}
if(!is.na(opt$variants_tested)){
  if(!is_integer(opt$variants_tested)){
    cat("ERROR: Option --variants_tested (", opt$variants_tested, ") must be an integer\n"); quit();
  }
}
if(!is.na(opt$ylim_to)){
  if(!is_integer(opt$ylim_to)){
    cat("ERROR: Option --ylim_to (", opt$ylim_to, ") must be an integer\n"); quit();
  }
}
if(!is.numeric(opt$significance_level)){
    cat("ERROR: Option --significance_level (", opt$significance_level, ") must be numeric\n"); quit();
}

#============================================================================================================
#						Making sure input files exist
#============================================================================================================

check_file_exists(opt$input_table, "input_table")
check_file_exists(opt$parameters_file, "parameters_file")


#============================================================================================================
#						Extracting GWAS variables from parameters_file
#============================================================================================================

cat("Extracting GWAS variables from parameters_file\n");

allele_frequency_from = vector(); allele_frequency_to = vector(); effect_size = vector();
homoplasy_steps_from = vector(); homoplasy_steps_to = vector(); sample_size = vector();
heritability = vector();

lines = as.vector(as.matrix(read.delim(opt$parameters_file, sep = "\n", header = F)))
for(l in 1:length(lines)){
  if(!startsWith(lines[l], "#")){
    if(startsWith(lines[l], "allele_frequency_from")){
      tmp = unlist(strsplit(lines[l], " ")); allele_frequency_from = unlist(strsplit(tmp[2], ","));
    }
    if(startsWith(lines[l], "allele_frequency_to")){
      tmp = unlist(strsplit(lines[l], " ")); allele_frequency_to = unlist(strsplit(tmp[2], ","));
    }
    if(startsWith(lines[l], "effect_size")){
      tmp = unlist(strsplit(lines[l], " ")); effect_size = unlist(strsplit(tmp[2], ","));
    }
    if(startsWith(lines[l], "homoplasy_steps_from")){
      tmp = unlist(strsplit(lines[l], " ")); homoplasy_steps_from = unlist(strsplit(tmp[2], ","));
    }
    if(startsWith(lines[l], "homoplasy_steps_to")){
      tmp = unlist(strsplit(lines[l], " ")); homoplasy_steps_to = unlist(strsplit(tmp[2], ","));
    }
    if(startsWith(lines[l], "sample_size")){
      tmp = unlist(strsplit(lines[l], " ")); sample_size = unlist(strsplit(tmp[2], ","));
    }
    if(startsWith(lines[l], "heritability")){
      tmp = unlist(strsplit(lines[l], " ")); heritability = unlist(strsplit(tmp[2], ","));
    }
  }
}

check_vector_not_empty(allele_frequency_from, "allele_frequency_from"); check_vector_not_empty(allele_frequency_to, "allele_frequency_to");
check_vector_not_empty(effect_size, "effect_size"); check_vector_not_empty(sample_size, "sample_size");
check_vector_not_empty(homoplasy_steps_from, "homoplasy_steps_from"); check_vector_not_empty(homoplasy_steps_to, "homoplasy_steps_to");
check_vector_not_empty(heritability, "heritability");

allele_frequency = paste(allele_frequency_from, allele_frequency_to, sep = ",")
cat("\tallele_frequency", paste(allele_frequency, sep = " "), "\n")
cat("\teffect_size", paste(effect_size, sep = " "), "\n")
homoplasy_steps = paste(homoplasy_steps_from, homoplasy_steps_to, sep = ",")
cat("\thomoplasy_steps", paste(homoplasy_steps, sep = " "), "\n")
cat("\tsample_size", paste(sample_size, sep = " "), "\n")
cat("\theritability", paste(heritability, sep = " "), "\n")

if(opt$plot_type == 3){
  if(!(length(homoplasy_steps)>1)){
    cat("ERROR: --plot_type 3 option cannot be used as homoplasy_steps parameter in",opt$parameters_file,"must include multiple levels\n"); quit();
  }
}
# if(length(homoplasy_steps)==1 & (homoplasy_steps[1] == '0,0')){ homoplasy_steps[1] = "1,1000"; }  # if selection of variants by degree of homoplasy disabled


#============================================================================================================
#						Loading input data
#============================================================================================================

cat("Loading input table with GWAS runs results", opt$input_table, "\n");
table = read.delim(opt$input_table, header = T, sep = '\t')

# Making sure opt$input_table has expected number of columns (n=23)
if(ncol(table) != 23){
  cat(paste("ERROR: Input file ",opt$input_table," has ",ncol(table)," number of columns, different from expected (n=23))\n",sep="")); quit();
}

if(opt$plot_type == 3){
  if(!opt$allele_frequency %in% table$allele_frequency){
    cat(paste("ERROR: Selected --allele_frequency value (",opt$allele_frequency,") not found in ",opt$input_table," table\n",sep="")); quit();
  }
}
  
#============================================================================================================
#						Keeping entries matching variables provided in parameters_file
#============================================================================================================

cat("Keeping entries in", opt$input_table,"matching variables provided in", opt$parameters_file, "\n");
# make sure expected column names exist
if(!is.element('allele_frequency', colnames(table))){ cat(paste("ERROR: Input file ",opt$input_table," has no column named 'allele_frequency'\n",sep="")); quit(); }
if(!is.element('effect_size', colnames(table))){ cat(paste("ERROR: Input file ",opt$input_table," has no column named 'effect_size'\n",sep="")); quit(); }
if(!is.element('homoplasy_steps', colnames(table))){ cat(paste("ERROR: Input file ",opt$input_table," has no column named 'homoplasy_steps'\n",sep="")); quit(); }
if(!is.element('sample_size', colnames(table))){ cat(paste("ERROR: Input file ",opt$input_table," has no column named 'sample_size'\n",sep="")); quit(); }
if(!is.element('heritability', colnames(table))){ cat(paste("ERROR: Input file ",opt$input_table," has no column named 'heritability'\n",sep="")); quit(); }

table = table[is.element(table$allele_frequency, allele_frequency),]
table = table[is.element(table$effect_size, effect_size),]
table = table[is.element(table$homoplasy_steps, homoplasy_steps),]
table = table[is.element(table$sample_size, sample_size),]
table = table[is.element(table$heritability, heritability),]

if(opt$plot_type == 3){
  table = table[is.element(table$allele_frequency, opt$allele_frequency),]
}

#============================================================================================================
#						Averaging p-values across repetitions
#============================================================================================================

cat("Averaging p-values across replicates for each unique parameter combination\n");
# print(dim(table))

# For each unique combination of allele frequency, effect size, homoplasy steps and sample size, p-values must be averaged
# creating unique parameter combinations
parameters_id = paste(table$allele_frequency, table$effect_size, table$homoplasy_steps, table$sample_size, table$heritability, sep = ";")
table$parameters_id = parameters_id
parameters_id = unique(parameters_id)

# -log10(p-value)
table$log10_pval = -log10(table$lrt.pvalue)

# Before averaging p-values for each parameter combination across replicates, decide whether to keep 
# variant 'sampling' replicates, 'phenotype' simulation replicates or 'both'
if(opt$average_pval_across == "sampling"){ table = subset(table, simulation_repetition == 1); }
if(opt$average_pval_across == "phenotype"){ table = subset(table, sample_repetition == 1); }
# print(dim(table))

# creating new table of averaged p-values (same dimentions as table, but with two extra columns, average and sd p-values)
avg_table = mat.or.vec(1, ncol(table) + 2)
colnames(avg_table) = c(colnames(table), "pval_avg", "pval_sd")
for(pc in 1:length(parameters_id))
{
  tmp = which(table$parameters_id == parameters_id[pc])
  pval_avg = mean(table$log10_pval[tmp])
  pval_sd = sd(table$log10_pval[tmp])
  newrow = table[tmp[1],]; newrow$pval_avg = pval_avg; newrow$pval_sd = pval_sd;
  avg_table = rbind(avg_table, newrow)
}
avg_table = avg_table[-1,]

#============================================================================================================
#						PLOTTING DATA
#============================================================================================================

cat("Creating GWAS power calculations plot and saving it as ", opt$output_plot, "\n");
library(ggplot2)

# changing dataframe values (needed to use aes shape=)
avg_table = transform(avg_table, effect_size = as.character(effect_size))
avg_table = transform(avg_table, heritability = as.character(heritability))

# plotting variables
sample_sizes = unique(sort(as.numeric(as.vector(avg_table$sample_size))));
sample_sizes_range = c(sample_sizes[1], sample_sizes[length(sample_sizes)]);
genome_wide_sig_level = -log10((opt$significance_level)/(opt$variants_tested)); # Bonferroni correction
axis_lines_width = 0.5
font = "Times";

# Plot 1 - across sample sizes and effect sizes
if(opt$plot_type == 1){
  avg_table$group1 = paste(avg_table$allele_frequency, avg_table$effect_size, sep = " ")
  plot = ggplot(data = avg_table, mapping = aes(x = sample_size, y = pval_avg, group = group1))
  plot = plot + theme_light()
  plot = plot + geom_line(aes(color=allele_frequency))
  plot = plot + geom_point(aes(color=allele_frequency, shape=effect_size))
  plot = plot + labs(title="GWAS power calculations", subtitle = "For allele frequency/effect size combinations across sample sizes", x = "Sample size", y = "-log10(p-value)")
  plot = plot + scale_x_continuous(breaks=sample_sizes, labels=as.character(sample_sizes))
  plot = plot + geom_hline(yintercept=genome_wide_sig_level, linetype="dashed", size = axis_lines_width)
  # plot = plot + theme(axis.text.x = element_text(size = 7))             # x-axis text size
  if(!is.na(opt$ylim_to)){
    plot = plot + coord_cartesian(ylim=c(0, as.integer(opt$ylim_to)))
  }
  if(opt$add_error_bars == "yes"){ 
    plot = plot + geom_errorbar(aes(ymin=pval_avg-pval_sd, ymax=pval_avg + pval_sd), width=.2, position=position_dodge(0.05));
    }
}

# Plot 2 - Across sample sizes and heritability values
if(opt$plot_type == 2){
  avg_table$group1 = paste(avg_table$allele_frequency, avg_table$heritability, sep = " ")
  plot = ggplot(data = avg_table, mapping = aes(x = sample_size, y = pval_avg, group = group1))
  plot = plot + theme_light()
  plot = plot + geom_line(aes(color=allele_frequency))
  plot = plot + geom_point(aes(color=allele_frequency, shape=heritability))
  plot = plot + labs(title="GWAS power calculations", subtitle = "For allele frequency/heritability combinations across sample sizes", x = "Sample size", y = "-log10(p-value)")
  plot = plot + scale_x_continuous(breaks=sample_sizes, labels=as.character(sample_sizes))
  plot = plot + geom_hline(yintercept=genome_wide_sig_level, linetype="dashed", size = axis_lines_width)
  if(!is.na(opt$ylim_to)){
    plot = plot + coord_cartesian(ylim=c(0, as.integer(opt$ylim_to)))
  }
  if(opt$add_error_bars == "yes"){ 
    plot = plot + geom_errorbar(aes(ymin=pval_avg-pval_sd, ymax=pval_avg + pval_sd), width=.2, position=position_dodge(0.05));
  }
}

# Plot 2 - For a given sample size
if(opt$plot_type == 3){
  avg_table$group1 = paste(avg_table$homoplasy_steps, avg_table$effect_size, sep = " ")
  plot = ggplot(data = avg_table, mapping = aes(x = sample_size, y = pval_avg, group = group1))
  plot = plot + theme_light()
  plot = plot + geom_line(aes(color=effect_size))
  plot = plot + geom_point(aes(color=effect_size, shape=homoplasy_steps))
  plot = plot + labs(title="GWAS power calculations", subtitle = paste("For homoplasy level/effect size combinations across sample sizes (allele frequency: ",opt$allele_frequency,")",sep=""), x = "Sample size", y = "-log10(p-value)")
  plot = plot + scale_x_continuous(breaks=sample_sizes, labels=as.character(sample_sizes))
  plot = plot + geom_hline(yintercept=genome_wide_sig_level, linetype="dashed", size = axis_lines_width)
  plot = plot + theme(axis.text.x = element_text(size = 7))             # x-axis text size
  if(!is.na(opt$ylim_to)){
    plot = plot + coord_cartesian(ylim=c(0, as.integer(opt$ylim_to)))
  }
  if(opt$add_error_bars == "yes"){ 
    plot = plot + geom_errorbar(aes(ymin=pval_avg-pval_sd, ymax=pval_avg + pval_sd), width=.2, position=position_dodge(0.05));
  }
}

ggsave(opt$output_plot, plot = plot, device = "pdf", width = opt$plot_width, height = opt$plot_height, dpi = 300, units = "cm")






