# convSig R Package

Convolutional Filtering Combined with Expectation Maximisation for the Detection of Mutational Processes in the Cancer Genome

Implementation Work Based on Novel Method Developed by [Yun Feng](https://github.com/yun-feng)

## Package Overview

This package implements an entirely novel approach for the detection and extraction of mutational processes that is based on Yun Feng’s work (convSig project, Oxford University PhD).
Every cancer-affected patient has specific mutational patterns in their genome, hence detecting, extracting, and characterising such patterns can provide a useful source of information for the personalisation of cancer treatment. The methods implemented in this package make use of the convolutional filtering techniques found in computer vision. Our model is centred around the concept of each mutational process being represented as a filter, which slides over the entire genome of interest and probabilistically induces mutations in genetic sequences that it has a strong rapport with. Our methods extract and characterise these processes using two variations of expectation maximisation (i.e. EM).

## Data Preparation for Expectation Maximisation

The EM methods (ReLU and exponential transformations) both require multiple input data that conform to specific data format, therefore the data prepration outlined in this section will certainly be a requisite step in most usages of this analysis pipeline.
The earliest upstream data input file is an ICGC (i.e. International Cancer Genome Consortium) mutation file either in the form of a .csv file or a .tsv file. Please note that there currently exists more than one ICGC format (e.g. donor format), however the one which this pipeline is concerned with is the ICGC mutation file format. To find a typical example, please visit the online ICGC data portal (https://dcc.icgc.org/), navigate to DCC DATA RELEASE, select a project and download any file which contains “simple_somatic_mutation” in its name. Or alternatively you can run the below line to load a suitable data input example into your global R environment:

```
loadICGCexample()
```

This input ICGC mutation data can either be a path leading to a mutation file on your workstation or data previously loaded into your R global environment. In the latter case, the data has to either be in the form of a data.frame or a matrix. This data preparation stage imposes an additional format requirement for this data. Namely, the following headers must be present in the input file (Note: this requirement is easily fulfilled if you have an input file that conforms to the ICGC format):
