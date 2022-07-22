#### internal data processing script
# clear workspace
rm(list = setdiff(ls(),c("codes.makepath","data.makepath")))
# set working directory to location of source code
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

#### init ----
library(maplet)
library(tidyverse)
library(magrittr)

## define all files
# clinical
file_clin <- file.path(data.makepath("COVID19/Choi_COVID19_data/ICU/clinical_lab_annotations/ICU_plasma_annotations.xlsx"))
file_clin_sheet <- "annotations"
file_clin_checksum <- "ee7231ed4495840736e84f8fc77c5c72"

# plasma metabolomics
file_plasma_metabo <- file.path(data.makepath("COVID19/Choi_COVID19_data/ICU/metabolomics/plasma_metabolite_intensity.xlsx"))
file_plasma_metabo_checksum <- "269691828d3e9159970a83b27c58d920"
# plasma lipidomics
file_plasma_lipid <- file.path(data.makepath("COVID19/Choi_COVID19_data/ICU/metabolomics/plasma_lipids.xlsx"))
file_plasma_lipid_checksum <- "4b4f484245be4c642cedef4446d858be"
# plasma proteomics
file_plasma_prot <- file.path(data.makepath("COVID19/Choi_COVID19_data/ICU/olink/plasma_proteomics.xlsx"))
file_plasma_prot_checksum <- "66af3506af1e249d6acef8259bb2555e"
file_plasma_prot_anno <- file.path(data.makepath("COVID19/Choi_COVID19_data/ICU/olink/proteomics_sampleinfo.xlsx"))
file_plasma_prot_anno_checksum <- "b4cd02a4b926c19e227d15e50b6c8d11"
# For loading any of the ICU datasets
# @returns D: a summarized experiment object
data_loader <- function(dataset) {
  
  if(!(dataset %in% c("plasma_metabo","plasma_lipids", "plasma_proteo")))
    stop("dataset can only be plasma_metabo or plasma_lipids or plasma_proteo")
  
  if (dataset=="plasma_metabo") {
    # input data
    # load data
    D <- mt_load_metabolon_v1(file=file_plasma_metabo, sheet = "OrigScale") 
    # load sample annotations
    D <- D %>%
      # verify data and annotations are unchanged
      mt_load_checksum(file = file_plasma_metabo, checksum = file_plasma_metabo_checksum) %>%
      mt_load_checksum(file = file_clin, checksum = file_clin_checksum) %>%
      # load clinical data
      mt_anno_xls(file=file_clin, sheet = file_clin_sheet, anno_type = "samples", anno_id = "Sample_ID", data_id= "CLIENT SAMPLE ID") %>%
      # make sure age is numeric
      mt_anno_mutate(anno_type = "samples", col_name = "age", term = as.numeric(as.matrix(age))) %>%
      # the two ards
      mt_anno_mutate(anno_type = "samples", col_name = "ARDSTypes", term = Group %>% recode(`Bact-Seps`="B", `Infl-ARDS`="none", `Co19-ARDS`="C", `Cont`="none", `Influ-Pneu`="none", `Co19-Pneu`="none"))  %>%
      
      # filter out samples without ards information
      mt_modify_filter_samples(filter=ARDSTypes !='none') %>%
      
      {.}
    
  } else if (dataset=="plasma_lipids"){
    sheet_list <- c("Species Concentrations", "Fatty Acid Concentrations")
    
    # load data
    D <- mt_load_metabolon_lipidomics(file = file_plasma_lipid, sheet_list = sheet_list)
    # load sample annotations
    D <- D %>%
      # verify data and annotations are unchanged
      mt_load_checksum(file = file_plasma_lipid, checksum = file_plasma_lipid_checksum) %>%
      mt_load_checksum(file = file_clin, checksum = file_clin_checksum) %>%
      # load clinical data
      mt_anno_xls(file=file_clin, sheet = file_clin_sheet, anno_type = "samples", anno_id = "Sample_ID", data_id= "CLIENT_SAMPLE_ID") %>%
      # make sure age is numeric
      mt_anno_mutate(anno_type = "samples", col_name = "age", term = as.numeric(as.matrix(age))) %>%
      # the two ards
      mt_anno_mutate(anno_type = "samples", col_name = "ARDSTypes", term = Group %>% recode(`Bact-Seps`="B", `Infl-ARDS`="none", `Co19-ARDS`="C", `Cont`="none", `Influ-Pneu`="none", `Co19-Pneu`="none"))  %>%
      # filter out samples without ards information
      mt_modify_filter_samples(filter=ARDSTypes !='none') %>%
      {.}
    
  } else if (dataset=="plasma_proteo"){
    # inpute data
    anno_sheet <- "ID matching plasma"
    # load data
    D <- mt_load_olink(file = file_plasma_prot)
    
    # load sample annotations
    D <- D %>%
      # verify data and annotations are unchanged
      mt_load_checksum(file = file_plasma_prot, checksum = file_plasma_prot_checksum) %>%
      mt_load_checksum(file = file_clin, checksum = file_clin_checksum ) %>%
      mt_load_checksum(file = file_plasma_prot_anno, checksum = file_plasma_prot_anno_checksum) %>%
      # load clinical data
      mt_anno_xls(file = file_plasma_prot_anno, sheet = anno_sheet, anno_type = "samples", anno_id = "sample id WCM-Q", data_id = "sample_id")%>%
      mt_anno_xls(file=file_clin, sheet = file_clin_sheet, anno_type = "samples", anno_id = "Prot_ID", data_id= "subject.id.WCM.NY")%>%
      
      # make sure age is numeric
      mt_anno_mutate(anno_type = "samples", col_name = "age", term = as.numeric(as.matrix(age))) %>%
      # filter out samples without group information
      mt_modify_filter_samples(filter = !is.na(Group)) %>%
      # same names as in plasma dataset
      mt_anno_mutate(anno_type = "samples", col_name = "Group", term = Group %>% recode(`sepsis_ARDS`="Bact-Seps", `Influenza_ARDS`="Infl-ARDS", `COVID19_ARDS`="Co19-ARDS", `CNTRL`="Cont", `Influenza_PNA`="Influ-Pneu", `COVID19_PNA`="Co19-Pneu")) %>%
      # the two ards
      mt_anno_mutate(anno_type = "samples", col_name = "ARDSTypes", term = Group %>% recode(`Bact-Seps`="B", `Infl-ARDS`="none", `Co19-ARDS`="C", `Cont`="none", `Influ-Pneu`="none", `Co19-Pneu`="none"))  %>%
      
      # filter out samples without ards information
      mt_modify_filter_samples(filter=ARDSTypes !='none') %>%
      
      {.}
    
  }
  
  # return 
  D
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
    data_loader(dataset=dataset) #%>%
    # preprocess
   # data_preprocessing(dataset=dataset)
  # select columns needed
  Y <- D %>% colData() %>% data.frame() %>% select(Subject_ID, sex, age, aki, platelet, pf, death, Group)
  D1 <- SummarizedExperiment(assay= assay(D), rowData=rowData(D), colData=Y)
  # save processed data
  D1 %<>% mt_write_se_xls(file = paste0("input/", dataset, '_raw.xlsx'))

} # closing dataset loop

## finished
print("Done! per omics data processing!") 
print("Generated excel files with raw data in input folder!") 

