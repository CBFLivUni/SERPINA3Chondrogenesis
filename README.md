# SERPINA3 is a marker of cartilage differentiation and drives extracellular matrix production in chondrogenesis
Comparative analysis of SERPINA3 knock-down in Chondrogenesis and Osteogenesis RNA-Seq time course data.


## Running the analysis

### With Singularity
Requires singularity to be [installed](https://apptainer.org/docs/admin/main/installation.html) 
```
#download the data and code
git clone https://github.com/CBFLivUni/SERPINA3Chondrogenesis.git

#move the downloaded repositry 
cd d.wilk_serpins

#run the analysis using the singularity container with all the needed dependencies
singularity run https://pgb.liv.ac.uk/~jsoul/DW/runAnalysis.img

```
### Manual installation and execution
```
#or install the needed R libraries
Rscript ./install/installRLibs.R

#run the analysis
Rscript runAnalysis.R

```

## Directories & Files

* ```[DW_E36_120123_Osteogenesis,DW_E36_120123_Chondrogenesis].[DATE].R_session_info.txt``` contains the R session info for the run of the pipeline on that date.
* ```DW_E36_120123_Osteogenesis``` and ```DW_E36_120123_Chondrogenesis``` contains:
	
	* ```data``` contains ```expr```, ```annot``` and ```metadata```:
		
		* ```annot``` contains ID conversions from biomart and REACTOME pathway-wise gene sets.
		* ```expr``` contains .rds and .tsv counts/expression matrices from the SkeletalViz pipeline (https://github.com/soulj/SkeletalVis-Pipeline).
		* ```metadata``` contains metadata for the Osteogenesis and Chondrogenesis experiments, respectively. These are important for the running the .Rmd documents.

	* ```results``` contains ```deseq2```, ```geva```, ```qc``` and ```pathway```:

		* ```deseq2``` contains the results tables .tsv's for DESeq2.
		* ```geva``` contains the results tables .tsv's for GEVA (https://github.com/sbcblab/geva).
		* ```qc``` contains the quality control analyses, including the PCA score plots.
		* ```pathway``` contains the GSEA pathway analysis results.
 
* ```Rmds``` contains the R markdown (.Rmd) file used to analyses the gene expression data (post-Nextflow pipeline).
* ```extra_scripts``` contains helper functions in R scripts for the .Rmd document(s).
* ```install``` provides the installation script for the needed R libraries
