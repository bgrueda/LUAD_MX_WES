# Whole Exome Sequencing of Mexican Lung Adenocarcinoma
This is a workflow to analyze WES (*Whole exome sequencing*) data in therms of the genomic variants identification.

## Motivation

Lung cancer is the leading cause of cancer incidence and mortality worldwide (Globocan 2022). In Mexico, while it ranks ninth in new cases, it is the third leading cause of cancer-related deaths.

Despite major advances in cancer genomics, most large-scale studies have focused on European or Asian populations, leaving Latin American populations underrepresented. This lack of representation can lead to biased conclusions that overlook population-specific mechanisms, and potentially misguide clinical decisions.

Latin American populations, including Mexicans, have a rich and complex genetic background shaped by migration, cultural diversity, and unique sociodemographic factors. Understanding how these differences influence cancer biology is crucial for achieving more equitable precision medicine.

One of the main challenges in countries like Mexico is the limited availability of large, well-annotated cohorts. In this study, we applied a bootstrapping-based strategy to extract robust insights from a small but valuable set of samples.

This repository is part of a broader effort to advance the genomic characterization of underrepresented populations, and to contribute to a more inclusive and accurate understanding of lung cancer biology.


## Workflow
This image below shows the steps of the Workflow for the identification and analysis of genomic variants and the scripts needed for this purpose.

![IMAGE](https://github.com/bgrueda/WES_MX_LUAD/blob/main/figures/BioinformaticWorkflow.png)

## Content
This repository has the following organization:

```
|--- bin
|   |-- 1_Quality.sh
|   |-- 2_Mapping.sh
|   |-- 3_Preprocessing.sh
|   |-- 4_Germline_ancestry.sh
|   |-- 5_Germline_variants.sh
|   |-- 6_Somatic_variants.sh
|   |__ 6.1_Signatures.py
|   |__ 6.2_Oncodrive-Mutsig.md
|   |__ 6.2_dNdScv.R
|   |__ 6.3_Bootstrap_analysis.R
|
|--- data
|   |-- data.md
|  
|--- figures

```

**[bin](/bin)**. All the scripts needed for the analysis.

**[data](/data)**. This contains information about the links and type of data to use this workflow. 

**[figures](/figures)**. This directory contains all the images that are related or resulted from this analysis.


## Data specifications
**Data files**
1. The data that are needed for this kind of analysis are **fastq files** which are the result of sequencing of your samples (raw data).
2. A reference genome version (hg19). This can be found in one of the links in *data* folder.

**Sample preparation requirements**
+ **DNA** samples
+ *Tumor and/or normal* tissue
+ *Exon enrichment*
+ Sequenced in a NGS platform (Illumina Nextseq, paired ends)
+ Calculating an average depth of **30X** approximately.  

## Requirements and software versions
- [FastQC (v.0.11.9)](https://www.bioinformatics.babraham.ac.uk/projects/download.html#fastqc)
- [Trimmomatic (v.0.39)](http://www.usadellab.org/cms/?page=trimmomatic)
- [BWA (v.0.7.12-r1039)](https://sourceforge.net/projects/bio-bwa/)
- [GATK (v.4.1.9.0)](https://github.com/broadinstitute/gatk/releases)
- [Samtools (v.1.10)](https://sourceforge.net/projects/samtools/)
- [R (v.3.4.2)](https://cran.r-project.org)
    + [MAFtoos (v.2.6.0)](https://bioconductor.org/packages/release/bioc/html/maftools.html)

## Contributors
+ Bertha Rueda-Zarazua
+ Humberto García-Ortiz

## Publication
{in revision}

## References
+ Sung H, Ferlay J, Siegel RL, Laversanne M, Soerjomataram I, Jemal A, Bray F. Global cancer statistics 2020:
GLOBOCAN estimates of incidence and mortality worldwide for 36 cancers in 185 countries. CA Cancer J Clin.
2021 Feb 4. doi: 10.3322/caac.21660. Epub ahead of print. PMID: 33538338.

## Useful sites: 
- https://gatk.broadinstitute.org/hc/en-us
- https://www.genepattern.org/modules/docs/MutSigCV/1
- https://bbglab.irbbarcelona.org/oncodrivefml/home
- https://cancer.sanger.ac.uk/signatures/tools/
- https://grch37.ensembl.org/Homo_sapiens/Tools/VEP