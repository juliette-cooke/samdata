# Functions for loading all types of data


#' Read in KO and WT sampling data
#' 
#' This function reads in KO and WT sampling data. It assumes that the results files are prefixed
#' using a 3 letter code (xaa, xab..) and that a code-to-disease dictionary is provided. It supports
#' using a single WT file or a specific WT file. Files should be in a /ko/ and /wt/ directories.
#' 
#' @param indir Path to the input dir, should contain a ko folder and wt folder if not using single.wt
#' @param disease.code Three-letter code designating the requested disease KO to load
#' @param single.wt If True, uses a non-disease-specific WT (provide the file using wt.file)
#' @param wt.file Use this to provide the wt path and filename when using single.wt=T
#' @returns A list containing ko and wt dataframes of all samples
#' @export
load_sampling_results = function(indir, disease.code, single.wt = F, wt.file = NA) {
  #disease.code="xaa"
  #indir="/home/juliette/these/data/IEM/sampling/results/wtopt_0_1_gzip_100000"
  # KO import
  ko.dir = paste0(indir,"/ko/")
  # List all files in the dir
  ko.files = list.files(ko.dir)
  # Get the file corresponding to the input disease code
  ko.file = ko.files[stringr::str_detect(ko.files, disease.code)]
  # Get the corresponding disease name using the disease code
  disease = iem_names$IEM.full.upper[match(stringr::str_sub(ko.file, start = 1, end = 3),iem_names$IEM.code)]
  
  # Read in the file
  ko = data.table::fread(paste0(ko.dir, ko.file))
  
  if (single.wt){
    # Only using one WT file
    if (is.na(wt.file)){
      stop("Provide the path and filename to the WT file when using single.wt = T")
    }else{
      wt = data.table::fread(wt.file)
    }
  }else {
    # WT import
    wt.dir = paste0(indir,"/wt/")
    # List all files in the dir
    wt.files = list.files(wt.dir)
    # Get the file corresponding to the input disease code
    wt.file = wt.files[stringr::str_detect(wt.files, disease.code)]
    
    wt = data.table::fread(paste0(wt.dir, wt.file))
  }
  
  sdata = list(ko,wt)
  names(sdata) = c("ko","wt")
  sdata = lapply(sdata, as.data.frame)
  
  # Fix an issue with how certain metabolites are named
  sdata = lapply(sdata, function(x) {
    names(x)[names(x) =="EX_tetdece1crn"] = "EX_tetdece1crn(e)"
    names(x)[names(x) =="EX_c101crn"] = "EX_c101crn(e)"
    names(x)[names(x) =="EX_tetdec2crn"] = "EX_tetdec2crn(e)"
    names(x)[names(x) =="EX_glc(e)"] = "EX_glc_D(e)"
    x
    })

  return(sdata)
}

