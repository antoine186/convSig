# convSig R Package ![vers](https://img.shields.io/badge/Version-1-blue) ![CancGen](https://img.shields.io/badge/Made%204-Cancer%20Genomics-red)

Convolutional Filtering Combined with Expectation Maximisation for the Detection of Mutational Processes in the Cancer Genome

Implementation Work Based on Novel Method Developed by [Yun Feng](https://github.com/yun-feng)

# Installation

## Windows System

To install convSig on your Windows system, please download the convSig package folder from the main convSig GitHub page. Once downloaded, running the following code should complete the installation.

```
install.packages("convSig-master", repos = NULL, type="source")
```

Please note that you may encounter an issue with updating/installing a dependency package stringi during installation. In such case, please manually download [stringi](https://cran.r-project.org/web/packages/stringi/index.html) and install from source prior to convSig installation.

## Mac OS

To install convSig on your Mac OS system, please run the following command:

```
library(devtools)

install_github("antoine186/convSig")
```

Please note that you may encounter an issue with updating/installing a dependency package stringi during installation. In such case, please manually download [stringi](https://cran.r-project.org/web/packages/stringi/index.html) and install from source prior to convSig installation. It is advisable to update all of your R packages prior to convSig installation. Finally, while convSig is building on your system and when prompted to update dependencies, you should choose the option *None*.

# Package Overview

This package implements an entirely novel approach for the detection and extraction of mutational processes that is based on Yun Feng’s work (convSig project, Oxford University PhD).
Every cancer-affected patient has specific mutational patterns in their genome, hence detecting, extracting, and characterising such patterns can provide a useful source of information for the personalisation of cancer treatment. The methods implemented in this package make use of the convolutional filtering techniques found in computer vision. Our model is centred around the concept of each mutational process being represented as a filter, which slides over the entire genome of interest and probabilistically induces mutations in genetic sequences that it has a strong rapport with. Our methods extract and characterise these processes using two variations of expectation maximisation (i.e. EM).

## Data Preparation for Expectation Maximisation

The EM methods (ReLU and exponential transformations) both require multiple input data that conform to specific data format, therefore the data prepration outlined in this section will certainly be a requisite step in most usages of this analysis pipeline.
The earliest upstream data input file is an ICGC (i.e. International Cancer Genome Consortium) mutation file either in the form of a .csv file or a .tsv file. Please note that there currently exists more than one ICGC format (e.g. donor format), however the one which this pipeline is concerned with is the ICGC mutation file format. To find a typical example, please visit the online ICGC data portal (https://dcc.icgc.org/), navigate to DCC DATA RELEASE, select a project and download any file which contains “simple_somatic_mutation” in its name. Or alternatively you can run the below line to load a suitable data input example into your global R environment:

```
loadICGCexample()
```

This input ICGC mutation data can either be a path leading to a mutation file on your workstation or data previously loaded into your R global environment. In the latter case, the data has to either be in the form of a data.frame or a matrix. This data preparation stage imposes an additional format requirement for this data. Namely, the following headers must be present in the input file (Note: this requirement is easily fulfilled if you have an input file that conforms to the ICGC format):
- icgc_sample_id
- chromosome
- chromosome_start
- chromosome_end
- mutated_from_allele
- mutated_to_allele
- assembly_version
- sequencing_strategy

The below block of code designates the reference genome file that corresponds to the relevant experiment and the sequencing strategy used to obtain it.
This function call skips all of the X and Y chromosomes (i.e. sex chromosomes), gathers all of the mutations that pertain to experiments run with your specified assembly version and sequencing strategy, and applies a base to number conversion to help filter out invalid types of mutations (e.g. insertions or deletions).

```
datapath <- "simple_somatic_mutation.open.COCA-CN.tsv"
mut_file <- icgc2mut(datapath, "GRCh37", "WGS")
```

The below function, in a first instance, sorts your input in ascending order by chromosome number first and then by the mutation chromosome start position. This sorting is required as the downstream convolution operation can only scan the reference genome in an ordered manner.
In a second instance (optional: depending on user parameter specification), it removes all non-single nucleotide polymorphisms as our methods are only valid for single nucleotide substitutions.
In a third instance, it removes all duplicate mutations present in the input mutation data.

```
cur_mut_file <- icgc_curate(mut_file, remove.nonSNP = TRUE)
```

The mut_count() function constructs two essential count tables: the mutation and the genome background count tables.
The mutation count table is a matrix, which records the frequency of each mutation types for each sample present in your supplied mutation input data. Each of its rows therefore represents a sample and each of its columns represents a distinct mutation type/possibility. The number of possible mutation types is calculated in accordance to the length of the mutation signature window (i.e. 3 or 5 bases only for our methods). The formula used is: 4<sup>numbase−1</sup>×2×3 (where numbase is equal to the length of the window). We first perform 4<sup>numbase−1</sup>×2 because all flanking bases in the window (i.e. non-middle bases) have exactly 4 base possibilities and the middle base only has 2 possibilities. The latter is due to the reverse complementarity nature of the genome; This is to avoid double counting the mutations that have actually occurred. We furthermore multiply the intermediate result by 3 here because of the three possibilities for the middle base substitution.
The genome background count table is a vector with length equal to the total possible number of mutation types as calculated above. It also records the frequency that we could expect from each mutation under circumstances of equal probability for the occurrence of such mutations.
Please note that this function call requires you to supply it with the correct corresponding reference genome.
The function call can also take a parameter that helps indicate whether you wish to construct the count tables for a 3-base long or a 5-base long mutation signature window. By default, the call is set to run for a 3-base long window.

```
assembly <- "Homo_sapiens.GRCh37.dna.primary_assembly.fa"
EMu_prepped <- mut_count(assembly, cur_mut_file, five = FALSE)
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
