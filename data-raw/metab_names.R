metab_names = read.table("inst/extdata/metab_names_ids.tsv", sep = "\t",
                         col.names = c("Metab.names", "Metab.Recon1.IDs", "Metab.Recon2v2.IDs"))
usethis::use_data(metab_names)
