#!/usr/bin/env Rscript
#
# Copyright (c) 2020 The Broad Institute, Inc. All rights reserved.
#
suppressPackageStartupMessages(library("pacman"))
suppressPackageStartupMessages(p_load("optparse"))
suppressPackageStartupMessages(p_load("glue"))

options( warn = -1, stringsAsFactors=F )

# specify command line arguments
option_list <- list(
  make_option( c("-i", "--input"), action='store', type='character',  dest='gct_str', help='Path to input GCT file.'),
  make_option( c("-l", "--level"), action='store', type='character',  dest='level', help="Mode of report, 'ssc' - single-site-centric, 'gc' - gene-centric, 'gcr' - gene-centric-redundant.", default='ssc'),
  make_option( c("-t", "--id_type"), action='store', type='character',  dest='id_type', help="Notation of site-ids: 'sm' - Spectrum Mill; 'wg' - Web Gestalt; 'ph' - Philosopher", default='sm'),
  make_option( c("-o", "--id_type_out"), action='store', type='character',  dest='id_type_out', help="Type of site id for output: 'uniprot', 'refseq', 'seqwin', 'psp' (psp not implemented yet).", default='uniprot'),
  make_option( c("-a", "--acc_type_in"), action='store', type='character',  dest='acc_type_in', help="Type of accession number in 'rid' object in GCT file (uniprot, refseq, symbol).", default='refseq'),
  make_option( c("-s", "--seqwin_column"), action="store", type='character', dest='seqwin_column', help="Column containing flanking sequences, separated by '|'. Only relevant if '--id_type_out' = 'seqwin'", default='VMsiteFlanks'),
  make_option( c("-g", "--gene_symbol_column"), action='store', type='character',  dest='gene_symbol_column', help="Name of column listing gene names; used for gene centric reports.", default='geneSymbol'),
  make_option( c("-k", "--humanize_gene_symbol"), action='store', type='character',  dest='humanize_gene_names', help="if TRUE, gene symbols will be capitalized (for e.g. mouse or rat).", default=FALSE),
  make_option( c("-v", "--sgt_column"), action='store', type='character',  dest='sgt_column', help="Column used to collpase subgroup-top (SGT) reports.", default='subgroupNum'),
  make_option( c("-d", "--localized"), action='store', type='character',  dest='localized', help="CAUTION: it is NOT RECOMMENDED to set this flag to FALSE. If TRUE only fully localized sites will be considered." , default=TRUE),
  make_option( c("-m", "--mode"), action='store', type='character',  dest='mode', help="Determines how multiple sites per gene will be combined: sd - most variable (standard deviation) across sample columns; SGT - subgroup top: first subgroup in protein group (Spectrum Mill); abs.max - for log-transformed, signed p-values" , default='median'),
  make_option( c("-r", "--residue"), action='store', type='character',  dest='residue', help='Modified residues, e.g. "S|T|Y" or "K".', default='S|T|Y'),
  make_option( c("-p", "--ptm"), action='store', type='character',  dest='ptm', help='Type of modification, e.g "p" or "ac".', default = 'p'),
  make_option( c("-u", "--preprocess_gct"), action='store', type='character',  dest='preprocess_gct', help='If FALSE nothing will be done; probably needed for to make this step optional in a FireCLoud WDL.', default = FALSE),
  make_option( c("-z", "--libdir"), action='store', type='character',  dest='libdir', help='Folder to source from.', default = 'c:/Users/karsten/Dropbox/Devel/PANOPLY/src/panoply_ssgsea/'),
  make_option( c("-y", "--yaml_file"), action='store', type='character',  dest='yaml_file', help='yaml parameter file.', default = 'NA')

  
  )

# parse command line parameters
opt <- parse_args( OptionParser(option_list=option_list) )

source(glue("{opt$libdir}/preprocess_gct.R"))
source(glue("{opt$libdir}/parse_yaml_preprocess_gct.R"))

# parse command line parameters
opt <-parse_param_preprocess_gct(option_list) 


###################################################
##       run the function
out <- preprocessGCT(gct.str = opt$gct_str,
              level=opt$level,
              id.type=opt$id_type,
              id.type.out=opt$id_type_out,
              acc.type=opt$acc_type_in,
              seqwin.col=opt$seqwin_column,
              gene.col=opt$gene_symbol_column,
              humanize.gene.names=opt$humanize_gene_names,
              SGT.col=opt$sgt_column,
              loc=opt$localized,
              mode=opt$mode,
              mod.res = opt$residue,
              mod.type=opt$ptm,
              appenddim=F,
              preprocess.gct=opt$preprocess_gct
)
writeLines(out, con="fn.out")
