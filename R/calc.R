# Functions for calculating


#' Calculate the difference between random pairs of flux samples
#' 
#' This function uses sdata to calculate random subtractions between ko and wt sampling results.
#' Metabolite names should be in the form of "EX_xxxx(e)".
#' 
#' @param sdata An sdata list containing ko and wt dataframes
#' @param seed A seed to set for reproducibility
#' @importFrom magrittr %>%
#' @importFrom stringr str_replace_all
#' @importFrom dplyr select bind_cols rename_with mutate group_by summarise
#' @importFrom tidyr gather
#' @returns A "long" dataframe containing all diff samples for all metabolites
#' @export
calc_diff = function(sdata, seed = NA) {
  # Get metabs that are in the results
  metabs = data.frame(par = colnames(sdata$ko)[str_replace_all(str_replace_all(colnames(sdata$ko),"\\(e\\)", ""), "EX_",
                                                               "")
                                               %in% str_replace_all(metab_names$Metab.Recon2v2.IDs, "_e", "")])
  # Store both naming conventions
  metabs$us = str_replace_all(str_replace_all(metabs$par,"(?<!\\(e\\))$", "_e"), "\\(e\\)", "_e")
  
  # Calculate the difference between each random pair of fluxes (wt vs ko)
  n = length(metabs$par)
  n1 = substr(n, 1, 1)
  n2 = substr(n, 2, 2)
  
  if (!(is.na(seed))){
    set.seed = seed
  }
  
  diff = sdata$ko %>%
    select(metabs$par) %>%
    bind_cols(sdata$wt %>%
                select(metabs$par)) %>%
    # Some regex to deal with the duplicate column names --> wty and mut
    rename_with( ~ gsub(paste0("\\.\\.\\.([1-9]$)|\\.\\.\\.([1-",as.numeric(n1)-1,"][0-9]$)|\\.\\.\\.(",n1,"[0-",n2,"])$"), "_mut", .x, fixed = F)) %>%
    rename_with( ~ gsub(paste0("\\.\\.\\.(",n1,"[",as.numeric(n2)+1,"-9]$)|\\.\\.\\.([",as.numeric(n1)+1,"-9][0-9]$)|\\.\\.\\.(1[0-9][0-9]$)"), "_wty", .x, fixed = F)) %>%
    gather(key=Metab, value=Value) %>%
    mutate("Type" = substr(Metab, nchar(Metab)-2, nchar(Metab))) %>%
    mutate(Metab = substr(Metab, 1, nchar(Metab)-4)) %>%
    group_by(Metab) %>%
    # Difference: mut - wty
    summarise(Value = sample(Value[Type=="mut"], 100000) - sample(Value[Type=="wty"], 100000))
  
  return(diff)
}

#' Calculate a pseudo-z-score of a difference distribution
#' 
#' This function uses a calculated diff (calc_diff) to calculate a pseudo-z-score for the 
#' distribution.
#' 
#' @param diff A "long" dataframe containing Metabs and diff values
#' @param disease.code The 3-letter disease code for the current disease
#' @returns A zscore dataframe containing zscores for each metabolite
#' @importFrom magrittr %>%
#' @export
calc_zscore = function(diff, disease.code){
  zscore = diff %>%
    group_by(Metab) %>%
    summarise(!!disease.code := mean(Value) / sd(Value))
  return(zscore)
}


#' Calculate a monte carlo p-value estimation
#' 
#' This function uses a z-score prediction table to estimate p-values using monte carlo estimation. It calculates
#' the significance of a z-score compared to the metabolite's "normal behaviour" over multiple conditions. 
#' 
#' @param diff A dataframe containing Metabs as rows and conditions as columns
#' @returns A dataframe containing estimated p-values for each predicted value
#' @importFrom magrittr %>%
#' @export
calc_mc_pval = function(zscore.df){
  zscore.df = read.csv("~/these/code/R/samdata/tests/test_results/pred.csv",check.names = F)
  # Transpose the df so that metabs are columns
  zscore.df.t = zscore.df %>% tibble::column_to_rownames("Metab") %>% t() %>% as.data.frame()
  
  pred.mc = lapply(zscore.df.t, FUN=function(x) {  # Each column (metab) is x
    names(x) = rownames(zscore.df.t)
    x = lapply(x, function(y) {  # Each value of the column is y
      ifelse(y < 0, log2((sum(x<y)+1)/(length(x)+1)), -log2((sum(x>y)+1)/(length(x)+1)) )
    })
    x
  })
  # Unlist and turn into data frame:
  temp = as.data.frame(apply(do.call(cbind, pred.mc), 2, unlist))
  # Transpose back to the correct format
  pred = temp  %>% t() %>% as.data.frame() %>% tibble::rownames_to_column("Metab") 

  return(pred)
}





