# Functions for calculating statistics on predicted vs observed data

#' Calculate precision and recall for predicted and observed results
#' 
#' @param pred.num A data frame containing numeric prediction values with metabolites as rows and conditions as columns
#' @param obs An data frame with metabolites as rows and conditions as columns
#' @returns A named vector containing precision and recall
#' @importFrom magrittr %>%
#' @importFrom tibble column_to_rownames rownames_to_column
#' @importFrom dplyr mutate select
#' @export
stat_counts = function(pred.num, obs) {
  # TODO: make obs generic ####
  # pred.num=pred.temp
  
  # Turn metab IDs into metab names
  pred.num$Metab = metab_names$Metab.names[match(pred.num$Metab, metab_names$Metab.EX.IDs)]
  # Make sure the IEM names are the correct form
  colnames(pred.num) = c("Metab",iem_names$IEM.short[match(names(pred.num[,-1]),
                                                       iem_names$IEM.code)])
  pred.num = pred.num %>% column_to_rownames("Metab")
  
  obs_f = obs
  colnames(obs_f) = c("Metab",iem_names$IEM.short[match(names(obs_f[,-1]),
                                                        iem_names$IEM.full.upper)])
  obs_f = obs_f %>% column_to_rownames("Metab")
  # Re-order colnames & rownames based on obs
  pred.num = pred.num[,colnames(obs_f)]
  pred.num = pred.num[rownames(obs_f),]
  
  pred.char = pred.num
  pred.char[pred.num < 0] = "-1"
  pred.char[pred.num > 0] = "+1"
  pred.char[pred.num == 0] = "00"
  
  
  TP_pos = length(pred.char[pred.char == obs_f & pred.char == "+1" & pred.char != "00" & !is.na(pred.char) & !is.na(obs_f)])
  TP_neg = length(pred.char[pred.char == obs_f & pred.char == "-1" & pred.char != "00" & !is.na(pred.char) & !is.na(obs_f)])
  FN_neg = length(obs_f[obs_f == "-1" & pred.char == "00" & !is.na(pred.char) & !is.na(obs_f)])
  FN_pos = length(obs_f[obs_f == "+1" & pred.char == "00" & !is.na(pred.char) & !is.na(obs_f)])
  FP_pos = length(pred.char[pred.char == "+1" & obs_f == "00" & !is.na(pred.char) & !is.na(obs_f)])
  FP_neg = length(pred.char[pred.char == "-1" & obs_f == "00" & !is.na(pred.char) & !is.na(obs_f)])
  TN = length(pred.char[pred.char == obs_f & pred.char == "00" & !is.na(pred.char) & !is.na(obs_f)]) 
  MP_pos = length(pred.char[pred.char == "-1" & obs_f == "+1" & !is.na(pred.char) & !is.na(obs_f)])
  MP_neg = length(pred.char[pred.char == "+1" & obs_f == "-1" & !is.na(pred.char) & !is.na(obs_f)])
  NA_count = length(pred.char[(pred.char == obs_f & is.na(pred.char)) | (pred.char != obs_f & is.na(pred.char)) 
                              | (pred.char != obs_f & is.na(obs_f))])
  stats = data.frame("Increase" = c(TP_pos, MP_neg,FP_pos),
                     "Decrease" = c(MP_pos, TP_neg,FP_neg),
                     "No change" = c(FN_pos,FN_neg,TN),
                     row.names = c("Increase", "Decrease", "No change"), check.names = F)
  
  counts = data.frame("Obs inc" = c(TP_pos, MP_pos, FN_pos),
                            "Obs dec" = c(MP_neg, TP_neg, FN_neg),
                            "Obs no change" = c(FP_pos, FP_neg, TN),
                            row.names = c("Pred inc", "Pred dec", "Pred no change"), check.names = F)
  return(counts)
}


#' @export
stat_pr_curve = function(pred, obs, min.thr=0, max.thr=NA, thr.step=0.1, method = c("merge_pos", "separate_dir", "only_pos")) {
  method = match.arg(method)
  # method = "separate_dir"
  # min.thr = 0
  # max.thr = max(pred[,-1])
  # thr.step = 0.1
  if (is.na(max.thr)) {
    max.thr = max(pred[,-1])
  }

  if (method == "merge_pos"){
    results = data.frame(matrix(ncol=3,nrow=0, dimnames=list(NULL, c("precision", "recall", "fpr"))))
  }else if (method == "separate_dir") {
    results = data.frame(matrix(ncol=6,nrow=0, dimnames=list(NULL, c("precision_inc", "recall_inc",
                                                                     "precision_dec", "recall_dec", 
                                                                     "fpr_inc", "fpr_dec"))))
  }else if (method == "only_pos"){
    results = data.frame(matrix(ncol=1,nrow=0, dimnames=list(NULL, c("accuracy"))))
  }
  
  # thr = 0
  for (thr in seq(from = min.thr, to = max.thr, by = thr.step)){
    pred.temp = pred
    if (method == "merge_pos"){
      pred.temp[,-1][(pred.temp[,-1] < thr & pred.temp[,-1] > -thr)] = 0
      
      counts = stat_counts(pred.temp, obs)
      # Merging TP with mis predictions => no matter the change direction
      TP = counts[1,1] + counts[2,2] + counts[1,2] + counts[2,1]
      FP = counts[1,3] + counts[2,3]
      FN = counts[3,1] + counts[3,2]
      TN = counts[3,3]
      
      precision = TP / (TP + FP)
      recall = TP / (TP + FN)
      fpr = FP / (FP + TN)
      
      stats = list(precision, recall, fpr)
      names(stats) = c("precision", "recall", "fpr")
    }else if (method == "separate_dir") {
      # Increases
      pred.temp[,-1][(pred.temp[,-1] <= 0)] = NA
      pred.temp[,-1][(pred.temp[,-1] > 0 & pred.temp[,-1] < thr)] = 0
      
      counts = stat_counts(pred.temp, obs)
      
      TP_inc = counts[1,1]
      FP_inc = counts[1,2] + counts[1,3]
      FN_inc = counts[2,1] + counts[3,1]
      TN_inc = counts[3,3] + counts[2,2] + counts[2,3] + counts[3,2]

      # Decreases
      pred.temp = pred
      pred.temp[,-1][(pred.temp[,-1] >= 0)] = NA
      pred.temp[,-1][(pred.temp[,-1] < 0 & pred.temp[,-1] > -thr)] = 0
      
      counts = stat_counts(pred.temp, obs)
      
      TP_dec = counts[2,2]
      FP_dec = counts[2,1] + counts[2,3]
      FN_dec = counts[1,2] + counts[3,2]
      TN_dec = counts[3,3] + counts[1,1] + counts[1,3] + counts[3,1]
      
      precision_inc = TP_inc / (TP_inc + FP_inc)
      recall_inc = TP_inc / (TP_inc + FN_inc)
      fpr_inc = FP_inc / (FP_inc + TN_inc)
      precision_dec = TP_dec / (TP_dec + FP_dec)
      recall_dec = TP_dec / (TP_dec + FN_dec)
      fpr_dec = FP_dec / (FP_dec + TN_dec)
      
      stats = list(precision_inc, recall_inc, precision_dec, recall_dec, fpr_inc, fpr_dec)
      names(stats) = c("precision_inc", "recall_inc","precision_dec", "recall_dec", "fpr_inc", "fpr_dec")
    }else if (method == "only_pos"){
      TP = counts[1,1] + counts[2,2]
      FP = counts[1,2] + counts[2,1]
      
      stats = list(TP / (TP + FP))
      names(stats) = c("accuracy")
    }
    
    # Add as a row
    if (thr == min.thr){
      # For the first row, it's a bit more complicated as the df is empty
      results = do.call(data.frame,setNames(stats, names(results)))
    }else{
      results = rbind(results, stats)
    }
  }
  rownames(results) = seq(from = min.thr, to = max.thr, by = thr.step)
  
  return(results)
}


#' @export
stat_roc_curve = function(pred.num, obs, min.thr=0, max.thr=NA, thr.step=0.1){
  if (is.na(max.thr)) {
    max.thr = max(pred.num[,-1])
  }
  
  roc = data.frame(matrix(ncol=2,nrow=0, dimnames=list(NULL, c("FPR", "TPR"))))
  
  for (thr in seq(from = min.thr, to = max.thr, by = thr.step)){
    pred.temp = pred.num
    pred.temp[,-1][(pred.temp[,-1] < thr & pred.temp[,-1] > -thr)] = 0
    counts = stat_pr(pred.temp, obs)
    TPR = as.numeric(stats$stats["Recall"])
    FPR = stats$counts[1,2] / (stats$counts[1,2] + stats$counts[2,2])
    # Add as a row
    if (thr == min.thr){
      # For the first row, it's a bit more complicated as the df is empty
      roc = do.call(data.frame,setNames(as.list(c(FPR, TPR)), names(roc)))
    }else{
      roc = rbind(roc, c(FPR, TPR))
    }
  }
  rownames(roc) = seq(from = min.thr, to = max.thr, by = thr.step)
  
  return(roc)
}

