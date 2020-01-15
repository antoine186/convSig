# convSig ![Rpack](https://img.shields.io/badge/R-Package-brightgreen) ![vers](https://img.shields.io/badge/Version-1-blue) ![gpllic](https://img.shields.io/badge/GPL-3%20and%20Higher-lightgrey) ![CancGen](https://img.shields.io/badge/Made%204-Cancer%20Genomics-red) ![Bioinf](https://img.shields.io/badge/Bio-Informatics-orange)

Convolutional Filtering Combined with Expectation Maximisation for the Detection of Mutational Processes in the Cancer Genome

Implementation Work Based on Novel Method Developed by [Yun Feng](https://github.com/yun-feng)

# Installation

## Windows System

To install convSig on your Windows system, please download the convSig package folder from the main convSig GitHub page. Once downloaded, running the following code should complete the installation.

```
install.packages("convSig-master", repos = NULL, type="source")
```

Please note that you may encounter an issue with updating/installing a dependency package stringi during installation. In such case, please manually download [![stringi](https://img.shields.io/badge/R-stringi-blue)](https://cran.r-project.org/web/packages/stringi/index.html) and install from source prior to convSig installation.

## Mac OS

To install convSig on your Mac OS system, please run the following command:

```
library(devtools)

install_github("antoine186/convSig")
```

Please note that you may encounter an issue with updating/installing a dependency package stringi during installation. In such case, please manually download [![stringi](https://img.shields.io/badge/R-stringi-blue)](https://cran.r-project.org/web/packages/stringi/index.html) and install from source prior to convSig installation. It is advisable to update all of your R packages prior to convSig installation. Finally, while convSig is building on your system and when prompted to update dependencies, you should choose the option *None*.

# Package Overview

This package implements an entirely novel approach for the detection and extraction of mutational processes that is based on Yun Feng’s work (convSig project, Oxford University PhD).
Every cancer-affected patient has specific mutational patterns in their genome, hence detecting, extracting, and characterising such patterns can provide a useful source of information for the personalisation of cancer treatment. The methods implemented in this package make use of the convolutional filtering techniques found in computer vision. Our model is centred around the concept of each mutational process being represented as a filter, which slides over the entire genome of interest and probabilistically induces mutations in genetic sequences that it has a strong rapport with. Our methods extract and characterise these processes using two variations of expectation maximisation (i.e. EM).

## ICGC Pipeline Quick Start 

This is a quick guide for the ICGC mutation pipeline. For more detailed information, please visit the convSig vignette.

Let us start with the ICGC mutation input data:

```
datapath <- "simple_somatic_mutation.open.COCA-CN.tsv"
```

The following function call provides a good example of the correct ICGC file format (note: we cannot use the example for the pipeline as it is not multi-sampled, for more information please visit coonvSig vignette):

```
loadICGCexample()
```

The pre-processing steps required for your input data is followed by the application of the following commands:

```
mut_file <- icgc2mut(datapath, assembly = "GRCh37", Seq = "WGS")
cur_mut_file <- icgc_curate(mut_file, remove.nonSNP = TRUE)
```

Once pre-processing of data is completed, we can generate the assets needed for convSig's expectation maximisation operations:

```
assembly <- "Homo_sapiens.GRCh37.dna.primary_assembly.fa"
EMu_prepped <- mut_count(assembly, cur_mut_file, five = TRUE)
```

Please note that the `five` parameter allows us to choose between 3-base and 5-base signatures extraction. If set to **FALSE**, 3-base signatures are extracted, otherwise 5-base signatures are extracted.

## VCF Pipeline Quick Start 

This is a quick guide for the VCF pipeline. For more detailed information, please visit the convSig vignette.

```
assembly <- "Homo_sapiens.GRCh37.dna.primary_assembly.fa"
mut_file <- "multisampled_mutations.vcf"
```

Running the following function call will pre-process our VCF input file and generate assets required for EM simultaneously.

```
EMu_prepped <- vcf2mut(mut_file, geno = "GT", assembly, five = TRUE)
```

### Expectation Maximisation via ReLU or Exponential Transformation

convSig provides us with two different approaches for applying EM. The ReLU transform is applied using the below command:

```
relu_res <- relu_transform(EMu_prepped, five = TRUE, K = 5)
```

Followed by the Exponential transform which is applied using the below command here:

```
exp_res <- exp_transform(EMu_prepped, five = TRUE, K = 5)
```
