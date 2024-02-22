#!/usr/bin/env Rscript
options( warn = -1 )

suppressPackageStartupMessages( if(!require("pacman")) install.packages ("pacman") )
suppressPackageStartupMessages(p_load("optparse"))

# # parse the directory this file is located (THIS DOESN'T WORK. using opt$libdir in option_list)
# script.dir <- commandArgs()[4]
# script.dir <- sub('^(.*(/|\\\\)).*', '\\1', sub('.*?\\=','', script.dir))
# # this doesn't seem to be reliably working on Windows OS
# # if called from the directory is 'ssgsea-cli.R' resides in
# # the line below attempts to fix this
# script.dir <- ifelse(dir.exists(script.dir), script.dir, '.')

# specify command line arguments
option_list <- list(
  make_option( c("-i", "--input"), action='store', type='character',  dest='input_ds', help='Path to input GCT file.'),
  make_option( c("-o", "--output"), action='store', type='character',  dest='output_prefix', help='File prefix for output files.', default='out'),
  make_option( c("-d", "--db"), action='store', type='character',  dest='gene_set_databases', help='Path to gene set database (GMT format).'),
  make_option( c("-n", "--norm"), action='store', type='character',  dest='sample_norm_type', help='Sample normalization: "rank", "log", "log.rank" or "none".', default = 'rank'),
  make_option( c("-w", "--weight"), action='store', type='character',  dest='weight', help='When weight==0, all genes have the same weight; if weight>0 actual values matter and can change the resulting score.', default = 0.75),
  make_option( c("-c", "--correl"), action='store', type='character',  dest='correl_type', help='Correlation type: "rank", "z.score", "symm.rank".', default = 'z.score'),
  make_option( c("-t", "--test"), action='store', type='character',  dest='statistic', help='Test statistic: "area.under.RES", "Kolmogorov-Smirnov"', default = 'area.under.RES'),
  make_option( c("-s", "--score"), action='store', type='character',  dest='output_score_type', help='Score type: "ES" - enrichment score,  "NES" - normalized ES', default = 'NES'),
  make_option( c("-p", "--perm"), action='store', type='character',  dest='nperm', help='Number of permutations', default = 1000),
  make_option( c("-m", "--minoverlap"), action='store', type='character',  dest='min_overlap', help='Minimal overlap between signature and data set.', default = 10),
  make_option( c("-q", "--tolerate_min_overlap_error"), action='store', type='logical',  dest='tolerate_min_overlap_error', help='Boolean to allow ssGSEA2() to tolerate input_ds having < min_overlap with gene_set_databases, without propogating error to shell.', default = FALSE),
  make_option( c("-x", "--extendedoutput"), action='store', type='character',  dest='extended_output', help='If TRUE additional stats on signature coverage etc. will be included as row annotations in the GCT results files.', default = TRUE),
  make_option( c("-e", "--export"), action='store', type='character',  dest='export_signat_gct', help='For each signature export expression GCT files.', default = TRUE),
  make_option( c("-g", "--globalfdr"), action='store', type='character',  dest='global_fdr', help='If TRUE global FDR across all data columns is calculated.', default = FALSE),
  make_option( c("-l", "--lightspeed"), action='store', type='character',  dest='multi_core', help='If TRUE processing will be parallized across gene sets. (I ran out of single letters to define parameters...)', default = TRUE),
  make_option( c("-y", "--yaml"), action='store', type='character',  dest='yaml_file', help='Parameter file (.yaml)', default = NA),
  make_option( c("-z", "--scrdir"), action='store', type='character',  dest='script.dir', help="Folder where 'ssgsea-cli.R' script is located.", default = '.')
)

## #####################################
# parse script-directory straight from command line inputs. or use '.' by default.
script.dir = parse_args( OptionParser(option_list=option_list) )$script.dir # create script.dir variable for backwards compatibility

## source the actual script
source(file.path(script.dir, 'src', 'ssGSEA2.0.R'))
source(file.path(script.dir, 'src', 'parse_yaml_ssgsea.R'))

# parse command line parameters
opt <- parse_param_ssgsea(option_list) # reparse args with our special yaml overwrite function
 
# hard-coded parameters
spare.cores <- 0 # use all available cpus
log.file <- paste(opt$output_prefix, '_ssgsea.log.txt', sep='')


## ######################################################################################################
##
##                   run ssGSEA
##
## ######################################################################################################
res <- tryCatch(ssGSEA2(
  input.ds=opt$input_ds,
  output.prefix=opt$output_prefix,
  gene.set.databases=opt$gene_set_databases,
  sample.norm.type=opt$sample_norm_type,
  weight=opt$weight,
  statistic=opt$statistic,
  output.score.type=opt$output_score_type,
  nperm=opt$nperm,
  min.overlap=opt$min_overlap,
  correl.type=opt$correl_type,
  export.signat.gct=opt$export_signat_gct,
  extended.output=opt$extended_output,
  global.fdr=opt$global_fdr,
  par=opt$multi_core,
  spare.cores=spare.cores,
  log.file=log.file
), error = function(e) {
  if (grepl("does not meet minimum-overlap", e) && # if we have a minimum-overlap error
      opt$tolerate_min_overlap_error) { # AND we have chosen to tolerate minimum-overlap errors
    message(paste0("\n### WARNING\n",e)) # print minimum overlap error as warning, but do not stop()
  } else stop(e) # otherwise, print error as normal and stop
} )




