library("seqinr")
library("Biostrings")
library("MASS")
library("GenomicRanges")
library("dndscv")

## Optional: install the latest version of dndscv from GitHub if not already installed
library(devtools)
devtools::install_github("im3sanger/dndscv")

### Driver discovery (positive selection) in cancer exomes/genomes

## INPUT format:
#sampleID	chr	pos	ref	mut
#P01	1	16959851	A	C
#P02	1	31836951	G	T

#Load input mutation data
pul <- read.csv(file = "samples_dndscv.csv", header = T, sep = "\t")
head(pul)

# Ensure the position column is numeric
pul$pos = as.numeric(as.character(pul$pos))

# Run dNdScv analysis to infer positive selection
dndsout = dndscv(pul)

# Extract results from the global positive selection model
sel_cv = dndsout$sel_cv
print(head(sel_cv), digits = 3)

# Here you can filter genes
signif_genes = sel_cv[sel_cv$qglobal_cv<0.2, c("gene_name","qglobal_cv")]
rownames(signif_genes) = NULL
print(signif_genes)

# Save to file
write.csv(sel_cv, "dndscv_output.csv", row.names=FALSE, quote=FALSE) 

# Print global dN/dS ratios for mutation classes
print(dndsout$globaldnds)

# Print the theta parameter used in the negative binomial regression
print(dndsout$nbreg$theta)

# Identify genes under positive selection according to the local model (more sensitive)
signif_genes_localmodel = as.vector(dndsout$sel_loc$gene_name[dndsout$sel_loc$qall_loc<0.1])
print(signif_genes_localmodel)

