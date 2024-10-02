#
# Copyright (c) 2020 The Broad Institute, Inc. All rights reserved.
#
################################################
## funtion to parse parse and update parameters
## - cmd line
## - yaml file
## parameters in yaml file will be updated with 
## parameters specified on cmd
parse_param_preprocess_gct <- function(cmd_option_list, yaml_section='panoply_preprocess_gct'){
  
  ## #########################################################
  # parse command line parameters
  opt_cmd <- parse_args( OptionParser(option_list=cmd_option_list) )
  
  ##save(opt_cmd, file='debug.RData')
  
  ############################################################
  ## parse yaml file
  if(!is.na(opt_cmd$yaml_file)) {
    
    if(file.exists(opt_cmd$yaml_file)){
    
      p_load(yaml)
      
      ## import yaml
      opt_yaml <- read_yaml(opt_cmd$yaml_file)
      
      ## extract relevant section
      opt_yaml <- opt_yaml[[yaml_section]]
      
      if(length(opt_yaml) > 0){
        
          ## parse cmd params
          cat('\n\nparsing command line parameters:\n', paste0(rep('-', 60), collapse = ''), '\n', sep='')
          for(x in names(opt_cmd))
            cat('---', x, opt_cmd[[x]], '; prefer yaml?', opt_cmd[[x]] == 'NA','\n')
          
          ## update yaml with parameters specified on cmd line 
          cat('\n\nUpdating parameter file with command line parameters:\n', paste0(rep('-', 60), collapse = ''), '\n', sep='')
          
          ## cmd parameters
          cmd_not_null <- which( !sapply(opt_cmd, function(x) x == 'NA' ) )
          cmd_to_update <- intersect( names(opt_cmd)[ cmd_not_null], names(opt_yaml) )
          cmd_to_add <- setdiff( names(opt_cmd), names(opt_yaml) )
          
          
          sapply(cmd_to_update, function(x) cat(x, ':', opt_yaml[[x]], '->', opt_cmd[[x]], '\n'))
         
          ## update yaml by cmd 
          opt_yaml[cmd_to_update] <- opt_cmd[cmd_to_update]
          
          cat(paste0(rep('-', 60), collapse = ''), '\n\n')
          
          ## add parameters only specified on cmd
          if(length(cmd_to_add) > 0){
            opt_cmd_to_add <- opt_cmd[cmd_to_add]
            opt_yaml <- append(opt_yaml, opt_cmd_to_add)
          }
          
          ## updated params
          opt <- opt_yaml
       } ## end if(length(opt_yaml) > 0)
     } ## end if(file.exists(opt_cmd$yaml_file))   
    } else {
    ## no yaml file
    opt <- opt_cmd
  }
  
  #########################################
  ## force correct mode
  opt$level <- as.character(opt$level)
  opt$id_type <- as.character(opt$id_type)
  opt$id_type_out <- as.character(opt$id_type_out)
  opt$acc_type_in <- as.character(opt$acc_type_in)
  opt$seqwin_column <- as.character(opt$seqwin_column)
  opt$gene_symbol_column <- as.character(opt$gene_symbol_column)
  opt$sgt_column <- as.character(opt$sgt_column)
  opt$localized <- as.logical(opt$localized)
  opt$mode <- as.character(opt$mode)
  opt$ptm <- as.character(opt$ptm)
  opt$preprocess_gct <- as.logical(opt$preprocess_gct)
  opt$libdir <- as.character(opt$libdir)
  
  return(opt)
} 
