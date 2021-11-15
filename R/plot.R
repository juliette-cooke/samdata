# Functions for plots

#' Plot a heatmap representation of WT vs KO sampling results
#' 
#' This functions plots a heatmap using previously calculated scores. It uses an existing observed dataset in order to 
#' reproduce Fig 3 from the Recon 2 article.
#' 
#' @param pred The dataframe with Metabs as rows and IEMs as columns, with Metabs as column 1
#' @param title The title to be used for the plot
#' @param outfile The output path and filename for the plot
#' @param squash Whether to squash extreme values (squash>0) or not (squash=0)
#' @importFrom magrittr %>%
#' @import ggplot2
#' @export
plot_scale_heatmap = function(pred, title, outfile, squash=0){
  load(file = "data/metab_names.rda")
  load(file="data/iem_names.rda")
  load(file="data/obs.rda")
  
  # Turn metab IDs into metab names
  pred$Metab = metab_names$Metab.names[match(pred$Metab, paste0("EX_", metab_names$Metab.Recon2v2.IDs))]
  # Make sure the IEM names are the correct form
  colnames(pred) = c("Metab",iem_names$IEM.full.upper[match(trimws(tolower(names(pred[,-1]))),
                                                            trimws(tolower(iem_names$IEM.full.upper)))])
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