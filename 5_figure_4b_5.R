# Script to generate network of Figure 4b and subnetworks of Figure 5
#
# !! NOTE: Cytoscape needs to be running and allow incoming connections to
#          generate GGM networks (see also RCy3 package). !!
#
# Once the networks are printed in cytoscape we remove self-loops: Edit --> remove self-loops
# Extract biggest connected component in AKI and Platelet subnetworks: Tools --> Network analyzer --> Subnetwork creation ->  Extract connected components -> from the poplist choose the biggest one
# Then save it using File --> save as --> 'supplementary_file_1.cys'

#### SET UP --------
### clear workspace and set directories ------
rm(list = ls())
# set working directory to location of source code
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
### libraries ------
library(stringi) # string match
library(RColorBrewer) # colors
source('custom_functions.R') # customized functions
### input variables ------
# my color palette
annocols <- brewer.pal(n = 9, name = "Pastel1")
annocols[9] <- '#999999'; annocols[7] <- '#80DEEA'
annocols[3:4] <- c('#80CBC4', '#C5E1A5'); annocols[5] <- '#E6EE9C' 
# datasets
datasets <- c("plasma_metabo", "plasma_proteo", 'plasma_lipids')
# clinical manifestations
outcomes <- data.frame(outcome = c('aki', 'platelet', 'pf', 'death'),
                       outcome_type = c(rep('numeric', 3), 'binary'),
                       outcome_mode = c(rep('numeric', 3), 'character'))
# groups to compare
complst <- list(comps = list(c("Group","Co19-ARDS","Bact-Seps")), 
                check_groups = c("Bact-Seps", "Co19-ARDS"))
group_names <- sub('-', '_', complst$check_groups)




#### MAIN --------
### compile 'within' stats result files ------
compiled_stats <- get_compiled_stats(
  sub_set=paste0('results/tmp_', datasets, '_within_stats.xlsx'))

### get ggm edge list ------
ggm_edges <- get_multiomics_network()

### build a node attribute data frame ------
node_attributes <- data.frame(node_name=unique(c(as.matrix(ggm_edges$source), as.matrix(ggm_edges$target))), 
                              node_type="Lipid") %>% 
  dplyr::mutate(node_type = case_when(node_name%in%(compiled_stats %>% 
                                 filter(mol_type=='Proteins') %>% 
                                 pull(mol_uname)) ~ "Protein", 
                               node_name%in%(compiled_stats%>% 
                                 filter(mol_type=='Metabolites')  %>% 
                                 pull(mol_uname))~ "Metabolite",
                               TRUE~"Lipid")) %>% 
  dplyr::left_join(compiled_stats, by=c("node_name"= "mol_uname"))

### create Figure 4b ------
elist_to_cytoscape(node_attributes=node_attributes,
                   edge_list=ggm_edges,
                   outcome_name='overall', 
                   group_names=group_names,
                   collection_name='ggm_bon_05',
                   whole_net=T)

### create Figure 5a ------
# subgraph with neighbours that are significant at 10% FDR
get_manifestation_subnetwork (ggm_edges=ggm_edges,
                              neigh_order=0,
                              sig_col_suf = 'fdr10',
                              node_attributes= node_attributes %>% 
                                mutate(aki_Co19_ARDS_pval_score = case_when(aki_Co19_ARDS_pval_score<1.3 ~0, TRUE~aki_Co19_ARDS_pval_score),
                                       aki_Bact_Seps_pval_score = case_when(aki_Bact_Seps_pval_score<1.3 ~0, TRUE~aki_Bact_Seps_pval_score)), 
                              outcome_name='aki',
                              group_names=group_names, 
                              collection_name='group_outcome_subnet_n1') 

### create Figure 5b ------
# subgraph with neighbours that are significant at 10% FDR
get_manifestation_subnetwork (ggm_edges=ggm_edges,
                              neigh_order=0,
                              sig_col_suf = 'fdr10',
                              node_attributes= node_attributes %>% 
                                mutate(platelet_Co19_ARDS_pval_score = case_when(platelet_Co19_ARDS_pval_score<1.3 ~0, TRUE~platelet_Co19_ARDS_pval_score),
                                       platelet_Bact_Seps_pval_score = case_when(platelet_Bact_Seps_pval_score<1.3 ~0, TRUE~platelet_Bact_Seps_pval_score)), 
                              outcome_name='platelet',
                              group_names=group_names, 
                              collection_name='group_outcome_subnet2') 

### finished ------
print("Done!") 
print("Generated panels of figure 4b and figure 5 in cytoscape interface!") 
