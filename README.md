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

The following function call provides a good example of the correct ICGC file format:

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

## Expectation Maximisation Using ReLU or the Exponential Transformation

### Expectation Maximisation Parameters and Results

The ReLU and the Exponential transformations both enable the learning of key parameters via a novel expectation maximisation technique (EM) subject to a negative log-likelihood loss function. The parameters whose estimation are sought for are a Feature matrix (contains the weights of the mutational filters), a P matrix (contains the mutational rate for each mutational process in each sample), and a mat matrix (contains the probabilities for all mutation types; Take T mutating to C as an instance of a mutation type).
The Feature matrix has (numbase×4)−2 rows. This number covers all possible bases at each base index in a particular trimer fragment - 4 possible bases that pertain to a single base index are represented as 4 contiguous rows in the matrix. The Feature matrix also has the same number of columns as the number of mutational processes as indicated by the user.
The P matrix has as many rows as the number of samples present in your input data. It also contains the same number of columns as the number of mutational processes specified by the user.
The mat result is a 3 dimensional matrix with dimensions (2,3,number of mutational processes). The first dimension equals 2 because there are only two distinct original bases possible pre-mutation; That is because of the reverse complementarity nature of the genome. The second dimension equals 3 because there are only three possible mutated bases per each of the 2 original bases.

### ReLU or Exponential Transformation

The below function calls require the data output by the previous function call [i.e. mut_count()], a parameter specifying whether the transformation should consider a 3-base long window vs a 5-base long window, and a parameter indicating the number of mutational processes to extract from the pre-processed input data.

```
relu_res <- relu_transform(EMu_prepped, five = TRUE, K = 5)
```

```
exp_res <- exp_transform(EMu_prepped, five = TRUE, K = 5)
```
