# Functions for plots

#' Plot a heatmap representation of WT vs KO sampling results for Recon 2
#' 
#' This functions plots a heatmap using previously calculated scores. It uses an existing observed dataset in order to 
#' reproduce Fig 3 from the Recon 2 article.
#' 
#' @param pred The dataframe with Metabs as rows and IEMs as columns, with Metabs as column 1
#' @param title The title to be used for the plot
#' @param outfile The output path and filename for the plot
#' @param squash Whether to squash extreme values (squash>0) or not (squash=0)
#' @importFrom magrittr %>%
#' @importFrom dplyr left_join
#' @importFrom tidyr gather
#' @import ggplot2
#' @export
plot_scale_heatmap_recon2 = function(pred, title, outfile, squash=0){
  
  # Turn metab IDs into metab names
  pred$Metab = metab_names$Metab.names[match(pred$Metab, metab_names$Metab.EX.IDs)]
  # Make sure the IEM names are the correct form
  colnames(pred) = c("Metab",iem_names$IEM.full.upper[match(names(pred[,-1]),
                                                            iem_names$IEM.code)])
  # Re-order colnames based on obs
  pred = pred[,colnames(obs)]

  # Squash extreme values
  if (squash > 0){
    pred[,-1][pred[,-1] >squash] = squash
    pred[,-1][pred[,-1] < -squash] = -squash
  }
  
  obs.pred.merged = pred %>%
    gather(key = "IEM", value = "pred", 2:50) %>%
    left_join(., obs %>% gather(key = "IEM", value = "obs", 2:50), by=c("Metab", "IEM"))
  
  # Rename IEM names to a shorter version
  obs.pred.merged$IEM = iem_names$IEM.short[match(obs.pred.merged$IEM, iem_names$IEM.full.upper)]
  
  # Plot the heatmap
  hm = obs.pred.merged %>%
    ggplot(aes(x=factor(Metab, levels=metab_names$Metab.names), 
               y=forcats::fct_rev(factor(IEM, levels=iem_names$IEM.short)), color=as.factor(as.numeric(obs))))+
    # Sampling predicted values:
    geom_tile(color="black", aes(fill=pred))+  # Make tiles and fill them based on value
    scale_fill_gradient2(name="Predicted diff", low = "brown1", mid = "white", high = "cornflowerblue", na.value = "grey")+
    # Observed values:
    # Add + and - as geom text:
    geom_text(aes(label = ifelse(obs=="+1" , "+", ifelse(obs == "-1", "-", ""))),show.legend = F,fontface = "bold", size = 6) +
    # Create a dummy geom_point for the legend:
    geom_point(size = NA) +
    # Colour based on value:
    scale_color_manual(name="Observed",values = c("navy", "brown4", "white"), labels=c("Increased", "Decreased"), 
                       breaks = c("1", "-1"))+
    # Override Observed legend to add "+" and "-":
    guides(colour = guide_legend(override.aes = list(size = 4,shape = c("+", "-")))) + 
    coord_equal()+  # Make it square
    theme(axis.text.x=element_text(angle=90,hjust=0,vjust=0.5), axis.text.y = element_text(hjust=0))+
    scale_x_discrete(position = "top") +
    ggtitle(title)+
    xlab("Metabolites")+
    ylab("IEMs")
  
  ggsave(plot = hm, filename = outfile, width = 17, height = 10)
}


#' Plot a generic scaled heatmap representation
#' 
#' Compare predicted and observed for any set of metabolites and KOs. You must provide both the pred and obs datasets.
#' The first column of both obs and pred must be the "Metab" column.
#' The other columns are the condition columns. The column names must match between obs and pred.
#' Obs must contain "+1", "-1" and "00". Pred must be numeric.
#' Conditions (column names) must match the ones in condition_order.
#' Metabolites (first column) must match the ones in metab_order.
#' 
#' @param pred The predicted dataframe with Metabs as rows and conditions as columns, with Metabs as column 1
#' @param obs The observed dataframe with Metabs as rows and conditions as columns, with Metabs as column 1
#' @param title The title to be used for the plot
#' @param metab_order The order in which to plot the metabolites (must match metabolites in Metabs)
#' @param condition_order The order in which to plot the onditions (must match column names)
#' @param squash Whether to squash extreme values (squash>0) or not (squash=0)
#' @importFrom magrittr %>%
#' @importFrom dplyr left_join
#' @importFrom tidyr gather
#' @import ggplot2
#' @export
plot_scale_heatmap = function(pred, obs, title, metab_order, condition_order, squash){

  # TODO: remove plot test stuff ####
  # pred = read.csv("~/these/code/R/samdata/tests/test_results/pred.csv",check.names = F)
  # pred$Metab = metab_names$Metab.names[match(pred$Metab, metab_names$Metab.EX.IDs)]
  # colnames(pred) = c("Metab",iem_names$IEM.full.upper[match(names(pred[,-1]),
  #                                                          iem_names$IEM.code)])
  # metab_order = metab_names$Metab.names
  # condition_order = iem_names$IEM.short
  # title= ""
  # squash = 2
  
  # TODO: test with other datasets ####
  # TODO: check if obs needs to be +1 -1 00 ####
  
  # Squash extreme values
  if (squash > 0){
    pred[,-1][pred[,-1] >squash] = squash
    pred[,-1][pred[,-1] < -squash] = -squash
  }
  
  obs.pred.merged = pred %>%
    gather(key = "Condition", value = "pred", 2:length(pred)) %>%
    left_join(., obs %>% gather(key = "Condition", value = "obs", 2:length(pred)), by=c("Metab", "Condition"))
  
  # obs.pred.merged$Condition = iem_names$IEM.short[match(obs.pred.merged$Condition, iem_names$IEM.full.upper)]
  
  hm = obs.pred.merged %>%
    ggplot(aes(x=factor(Metab, levels=metab_order), 
               y=forcats::fct_rev(factor(Condition, levels=condition_order)), color=as.factor(as.numeric(obs))))+
    # Sampling predicted values:
    geom_tile(color="black", aes(fill=pred))+  # Make tiles and fill them based on value
    scale_fill_gradient2(name="Predicted diff", low = "brown1", mid = "white", high = "cornflowerblue", na.value = "grey")+
    # Observed values:
    # Add + and - as geom text:
    geom_text(aes(label = ifelse(obs=="+1" , "+", ifelse(obs == "-1", "-", ""))),show.legend = F,fontface = "bold", size = 6) +
    # Create a dummy geom_point for the legend:
    geom_point(size = NA) +
    # Colour based on value:
    scale_color_manual(name="Observed",values = c("navy", "brown4", "white"), labels=c("Increased", "Decreased"), 
                       breaks = c("1", "-1"))+
    # Override Observed legend to add "+" and "-":
    guides(colour = guide_legend(override.aes = list(size = 4,shape = c("+", "-")))) + 
    coord_equal()+  # Make it square
    theme(axis.text.x=element_text(angle=90,hjust=0,vjust=0.5), axis.text.y = element_text(hjust=0))+
    scale_x_discrete(position = "top") +
    ggtitle(title)+
    xlab("Metabolites")+
    ylab("IEMs")
  
  return(hm)
}