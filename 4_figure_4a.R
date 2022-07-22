# Script to generate panles of Figure 4a
# Once the PDFs are generated the figure was assembled in illustrator

#### SET UP --------
### clear workspace and set directories ------
rm(list = ls())
# set working directory to location of source code
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

### input variables ------
# my color palette
annocols <- c('#E6EE9C', '#C5E1A5', '#80CBC4', '#999999')
# datasets
datasets <- c("plasma_metabo", "plasma_proteo", 'plasma_lipids')
# clinical manifestations
outcomes <- data.frame(outcome = c('aki', 'platelet', 'pf', 'death'),
                       outcome_type = c(rep('numeric', 3), 'binary'),
                       outcome_mode = c(rep('numeric', 3), 'character'))
# groups to compare
complst <- list(comps = list(c("Group","Co19-ARDS","Bact-Seps")), 
                check_groups = c("Bact-Seps", "Co19-ARDS"))
# adjusted p-value cutoff
pcut <- 0.05

### libraries ------
library(RColorBrewer) # colors
source('custom_functions.R') # customized functions

#### MAIN --------
### compile 'within' stats result files ------
compiled_stats <- get_compiled_stats(
  sub_set=paste0('results/tmp_', datasets, '_within_stats.xlsx'))

### generate 'within' stats bar plots ------
fdr5_plist <- within_ards_barplots(compiled_stats=compiled_stats,
                                   outcomes=outcomes, annocols = annocols,
                                   sig_col_suf='fdr5',
                                   group_names=sub('-', '_', complst$check_groups))

### save to pdf ------
pdf('results/Figure4a_plasma_clinicals_fdr5.pdf', height=5, width=6)
grid.arrange(grobs=fdr5_plist, ncol=2)
dev.off()

###  finished ------
print("Done!") 
print("Generated PDF file with figure 4a in results folder!") 
