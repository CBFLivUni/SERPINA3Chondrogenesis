#load the quarto library
library(rmarkdown)

#run the analysis for both datasets
rmarkdown::render('Rmds/analysis_workflow.Rmd',params=list(study_id = "DW_E36_120123_Osteogenesis"))
rmarkdown::render('Rmds/analysis_workflow.Rmd',params=list(study_id = "DW_E36_120123_Chondrogenesis"))