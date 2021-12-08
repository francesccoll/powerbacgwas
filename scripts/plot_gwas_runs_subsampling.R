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
  make_option(c("-d", "--pyseer_model"), type="character", default="LMM", 
              help="Pyseer model used. To choose between: LMM (Linear Mixed Models) or ENET (Elastic Net)", metavar="character"),
  make_option(c("-m", "--plot_metrics"), type="character", default=1, 
              help="Select 1 or 2 to choose what metrics to plot: 1. power (percentage of GWAS runs with p-value above genome-wide significance) or 2. average p-value across GWAS runs (default: 1)", metavar="integer"),
  make_option(c("-p", "--parameters_file"), type="character", default="", 
              help="File with GWAS power calculations parameters used by prepare_gwas_runs_subsampling[_roary].py", metavar="character"),
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

if(opt$pyseer_model != "LMM" & opt$pyseer_model != "ENET"){
  cat("ERROR: pyseer_model option provided (", opt$pyseer_model, ") must be either LMM or ENET\n"); quit();
}
if(opt$plot_metrics != 1 & opt$plot_metrics != 2){
  cat("ERROR: plot_metrics option provided (", opt$plot_metrics, ") must be either 1 or 2\n"); quit();
}
if(opt$add_error_bars != "yes" & opt$add_error_bars != "no"){
  cat("ERROR: add_error_bars option provided (", opt$add_error_bars, ") must be either 'yes' or 'no'\n"); quit();
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

allele_frequency = vector(); effect_size = vector(); sample_size = vector();

lines = as.vector(as.matrix(read.delim(opt$parameters_file, sep = "\n", header = F)))
for(l in 1:length(lines)){
  if(!startsWith(lines[l], "#")){
    if(startsWith(lines[l], "allele_frequency")){
      tmp = unlist(strsplit(lines[l], " ")); allele_frequency = unlist(strsplit(tmp[2], ","));
    }
    if(startsWith(lines[l], "effect_size")){
      tmp = unlist(strsplit(lines[l], " ")); effect_size = unlist(strsplit(tmp[2], ","));
    }
    if(startsWith(lines[l], "sample_size")){
      tmp = unlist(strsplit(lines[l], " ")); sample_size = unlist(strsplit(tmp[2], ","));
    }
  }
}

check_vector_not_empty(allele_frequency, "allele_frequency");
check_vector_not_empty(effect_size, "effect_size");
check_vector_not_empty(sample_size, "sample_size");

cat("\tallele_frequency", paste(allele_frequency, sep = " "), "\n")
cat("\teffect_size", paste(effect_size, sep = " "), "\n")
cat("\tsample_size", paste(sample_size, sep = " "), "\n")

#============================================================================================================
#						Loading input data
#============================================================================================================

cat("Loading input table with GWAS runs results", opt$input_table, "\n");
table = read.delim(opt$input_table, header = T, sep = '\t')

# Making sure opt$input_table has expected number of columns (n=23)
if(ncol(table) != 18){
  cat(paste("ERROR: Input file ",opt$input_table," has ",ncol(table)," number of columns, different from expected (n=18))\n",sep="")); quit();
}

#============================================================================================================
#						Keeping entries matching variables provided in parameters_file
#============================================================================================================

cat("Keeping entries in", opt$input_table,"matching variables provided in", opt$parameters_file, "\n");
# make sure expected column names exist
if(!is.element('allele_frequency', colnames(table))){ cat(paste("ERROR: Input file ",opt$input_table," has no column named 'allele_frequency'\n",sep="")); quit(); }
if(!is.element('effect_size', colnames(table))){ cat(paste("ERROR: Input file ",opt$input_table," has no column named 'effect_size'\n",sep="")); quit(); }
if(!is.element('sample_size', colnames(table))){ cat(paste("ERROR: Input file ",opt$input_table," has no column named 'sample_size'\n",sep="")); quit(); }

table = table[is.element(table$allele_frequency, as.numeric(allele_frequency)),]
table = table[is.element(table$effect_size, as.numeric(effect_size)),]
table = table[is.element(table$sample_size, as.numeric(sample_size)),]

#============================================================================================================
#						Averaging p-values across repetitions
#============================================================================================================

cat("Averaging p-values across replicates for each unique parameter combination\n");
# print(dim(table))

# For each unique combination of allele frequency, effect size, homoplasy steps and sample size, p-values must be averaged
# creating unique parameter combinations
parameters_id = paste(table$allele_frequency, table$effect_size, table$sample_size, sep = ";")
table$parameters_id = parameters_id
parameters_id = unique(parameters_id)

# -log10(p-value)
# NOTE: when using ENET, lrt.pvalue (the p-value of association, adjusted for population structure.) is not calculated, thus filter.pvalue is chosen instead
table$log10_pval = -log10(table$lrt.pvalue)
if(opt$pyseer_model == "ENET"){ table$log10_pval = -log10(table$filter.pvalue); }
genome_wide_sig_level = -log10((opt$significance_level)/(opt$variants_tested)); # Bonferroni correction

# creating new table of averaged p-values (same dimentions as table, but with two extra columns, average and sd p-values)
avg_table = mat.or.vec(1, ncol(table) + 3)
colnames(avg_table) = c(colnames(table), "pval_avg", "pval_sd","power")
for(pc in 1:length(parameters_id))
{
  tmp = which(table$parameters_id == parameters_id[pc])
  pval_avg = mean(table$log10_pval[tmp])
  pval_sd = sd(table$log10_pval[tmp])
  # power is calculated as the percentage of runs with a p-value above genome-wide significance level
  power = (length(which(table$log10_pval[tmp] > genome_wide_sig_level))/length(tmp))*100
  newrow = table[tmp[1],]; newrow$pval_avg = pval_avg; newrow$pval_sd = pval_sd; newrow$power = power;
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
avg_table = transform(avg_table, allele_frequency = as.character(allele_frequency))

# plotting variables
sample_sizes = unique(sort(as.numeric(as.vector(avg_table$sample_size))));
sample_sizes_range = c(sample_sizes[1], sample_sizes[length(sample_sizes)]);
axis_lines_width = 0.5
font = "Times";
y_label = "NA";
if(opt$plot_metrics == 1){ y_label = "power"; }
if(opt$plot_metrics == 2){ y_label = "-log10(p-value)"; }


# Plot 1 - across sample sizes and effect sizes
avg_table$group1 = paste(avg_table$allele_frequency, avg_table$effect_size, sep = " ")
if(opt$plot_metrics == 1){ plot = ggplot(data = avg_table, mapping = aes(x = sample_size, y = power, group = group1)); }
if(opt$plot_metrics == 2){ plot = ggplot(data = avg_table, mapping = aes(x = sample_size, y = pval_avg, group = group1)); }
plot = plot + theme_light()
# plot = plot + scale_color_manual(values=c("#525252", "#969696", "#525252"))
plot = plot + geom_line(aes(color=allele_frequency))
plot = plot + geom_point(aes(color=allele_frequency, shape=effect_size))
plot = plot + labs(title="Sub-sampling GWAS power calculations", subtitle = "For allele frequency/effect size combinations across sample sizes", x = "Sample size", y = y_label)
plot = plot + scale_x_continuous(breaks=sample_sizes, labels=as.character(sample_sizes))
if(opt$plot_metrics == 2){ plot = plot + geom_hline(yintercept=genome_wide_sig_level, linetype="dashed", size = axis_lines_width); }
if(opt$plot_metrics == 1){ plot = plot + geom_hline(yintercept=80, linetype="dashed", size = axis_lines_width); }
if(!is.na(opt$ylim_to)){
  plot = plot + coord_cartesian(ylim=c(0, as.integer(opt$ylim_to)))
}
if(opt$add_error_bars == "yes"){ 
  plot = plot + geom_errorbar(aes(ymin=pval_avg-pval_sd, ymax=pval_avg + pval_sd), width=.2, position=position_dodge(0.05));
}

ggsave(opt$output_plot, plot = plot, device = "pdf", width = opt$plot_width, height = opt$plot_height, dpi = 300, units = "cm")






