
### SET UP --------
## clear workspace and set directories ----
rm(list = ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
## libraries ----
library(maplet)
library(tidyverse)
library(magrittr)
source('custom_functions.R')

# datasets
datasets <- c("plasma_metabo", "plasma_proteo", 'plasma_lipids')

#### MAIN --------

# medications effect ----
# medication data
med_file <- 'input/medication_data.xlsx'
medication_columns <- 'input/medication_columns.xlsx'
# medication correction ----
# loop over omics datasets
for (dataset in datasets) {
  outfile <- sprintf("results/supplementary_table_%s_medication_covid19.xlsx", dataset)
  # load processed data 
  D <- mt_load_se_xls(file=paste0('input/', dataset, '_processed.xlsx')) %>%
    # flag that data is logged
    mt_load_flag_logged() %>%
    # add additional medication information
    mt_anno_xls(file=med_file, sheet=1, anno_type = 'samples', anno_id_col = 'subject_id',
                data_id_col = 'Subject_ID') %>%
    # filter for covid samples
    mt_modify_filter_samples(filter=Group=="Co19-ARDS")
  
Dmc <- D %>% # converting all medication columns to factors
  mt_anno_class(anno_type = 'samples', file=medication_columns, sheet=1) %>%
  # medication correction for medications which were administered on more than 4 samples
  mt_pre_confounding_correction_stepaic(cols_to_correct = names(which(apply(Y, 2, FUN=function(x) length(which(x>0)))>4)), 
                                        cols_to_exclude = NULL)

# medication effect
df_medcor <- metadata(Dmc)$results[[grep('pre_confounding_correction_stepaic', names(metadata(Dmc)$results))]]$output %>% 
  data.frame() %>% bind_cols(rowData(Dmc) %>% data.frame() %>% select(name)) %>%
  select(-feature) %>% select(name, covariates, model.rsq, model.pvalue)%>%
  dplyr::rename(medications=covariates,
                model_rsq = model.rsq,
                model_pval = model.pvalue)

# write out medication effects i.e supplementary table 3
out <- 'medication_correction'; this_res <- df_medcor
wb <- openxlsx::createWorkbook()
# creat worksheet
openxlsx::addWorksheet(wb,sprintf('%s', out))
# write data
openxlsx::writeData(wb, sprintf('%s', out), this_res, rowNames = F, colNames = T)
# create and add a style to the column headers
headerStyle <- createStyle(fontName = 'Arial', fontSize = 12, halign = 'center', valign = 'center', textDecoration = 'bold')
# style for body
bodyStyle <- createStyle(fontName = 'Arial', fontSize = 12, halign = 'center', valign = 'center')
# apply style
addStyle(wb, sheet = sprintf('%s', out), bodyStyle, rows = 1:nrow(this_res), cols = 1:ncol(this_res), gridExpand = TRUE)
addStyle(wb, sheet = sprintf('%s', out), headerStyle, rows = 1, cols = 1:ncol(this_res), gridExpand = TRUE)

# write out
openxlsx::saveWorkbook (wb, file=outfile, overwrite=TRUE)
}
## finished
print("Done! estimating mediation effect completed.") 
print("Generated excel files with supplementary tables in results folder!") 
