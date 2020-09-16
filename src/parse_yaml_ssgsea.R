################################################
## funtion to parse parse and update parameters
## - cmd line
## - yaml file
## parameters in yaml file will be updated with 
## parameters specified on cmd
parse_param_ssgsea <- function(cmd_option_list, yaml_section='panoply_ssgsea'){
  
  ## #########################################################
  # parse command line parameters
  opt_cmd <- parse_args( OptionParser(option_list=cmd_option_list) )
 
  ############################################################
  ## parse yaml file
  if(!is.na(opt_cmd$yaml_file)) {
    
    if(file.exists(opt_cmd$yaml_file)){
    
      p_load(yaml)
      
      ## import yaml
      opt_yaml <- read_yaml(opt_cmd$yaml_file)
      
      ## extract relevant section
      opt_yaml <- opt_yaml[[yaml_section]]
      
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
  
      }    
    } else {
    ## no yaml file
    opt <- opt_cmd
  }
  
  #########################################
  ## force correct mode
  opt$nperm <- as.integer(opt$nperm)
  opt$weight <- as.numeric(opt$weight)
  opt$sample_norm_type <- as.character(opt$sample_norm_type)
  opt$correl_type <- as.character(opt$correl_type)
  opt$statistic <- as.character(opt$statistic)
  opt$output_score_type <- as.character(opt$output_score_type)
  opt$output_prefix <- as.character(opt$output_prefix)
  opt$min_overlap <- as.integer(opt$min_overlap)
  opt$extended_output <- as.logical(opt$extended_output)
  opt$export_signat_gct <- as.logical(opt$export_signat_gct)
  opt$global_fdr <- as.logical(opt$global_fdr)
  opt$multi_core <- as.logical(opt$multi_core)
  
  return(opt)
} 
