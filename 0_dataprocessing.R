#### data processing script
# The processed data used in this study can be downloaded at https://doi.org/10.6084/m9.figshare.19775359
# clear workspace
rm(list = ls())
# set working directory to location of source code
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

#### init ----
library(maplet)
library(tidyverse)
library(magrittr)
# run preprocessing for datasets
data_preprocessing <- function(D, dataset) {
  
  D <- D %>%
    mt_pre_filter_missingness(samp_max=0.5) %>%
    mt_pre_filter_missingness(feat_max=0.25) %>%
    {.}
  # normalization
  if(grepl("_proteo", dataset)){
    # undo log transformation of data for quot norm
    D %<>% mt_pre_trans_exp()
  } 
  # normalization, log, imputation
  D  %<>%
    mt_pre_norm_quot(feat_max = 0) %>%
    # transformation & imputation
    mt_pre_trans_log() %>% 
    mt_pre_impute_knn() %>% 
    {.}
  
  # average duplicate molecules --> specifically needed for proteins
  if (grepl("proteo", dataset)){
    D %<>% mt_modify_avg_features(group_col = 'name')  
  }
  # final dataset info
  D %<>%  # dataset info
    mt_reporting_data()
  
}

# make sure result directory exists
dir.create("input/", showWarnings = F, recursive = T)
datasets <- c('plasma_lipids', "plasma_metabo", "plasma_proteo") #

##### start analysis ----

# loop over datasets ----
for (dataset in datasets) {
  # load and preprocess
  D <-
    # load data
    mt_load_se_xls(file=sprintf('input/%s_raw.xlsx', dataset)) %>%
    # adjust the age variable with ">" sign
    mt_anno_mutate(anno_type = "samples", col_name='age', term=case_when(age%in%">90" ~ "90", TRUE~age)) %>%
    # make sure age is numeric
    mt_anno_mutate(anno_type = "samples", col_name = "age", term = as.numeric(as.matrix(age))) %>%
    # preprocess
    data_preprocessing(dataset=dataset)
  # select columns needed
  Y <- D %>% colData() %>% data.frame() %>% select(Subject_ID, sex, age, aki, platelet, pf, death, Group)
  D1 <- SummarizedExperiment(assay= assay(D), rowData=rowData(D), colData=Y)
  # save processed data
  D1 %<>% mt_write_se_xls(file = paste0("input/", dataset, '_processed.xlsx'))
  
} # closing dataset loop

## finished
print("Done! per omics data processing!") 
print("Generated excel files with processed data in input folder!") 

