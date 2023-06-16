# Monocytes_leishmania
Repository under progress

Update 16/06/23: Repository will be updated shortly after publication to include the download of the count matrix from GEO.


This repository contains scripts for the scRNAseq analysis of monocytes infected by Leishmania


## Data Accessibility
For now, the count matrix can be downloaded from the [Helmholtz Cloud](https://nubes.helmholtz-berlin.de/).
The count matrix should be manually download and place in the data folder of the repository.

## Repository Structure
**data:** Directory containing required raw data. The count matrix should be manually placed in this directory too. A metadata table describing relevant experimental information for each single-cell library is also included (Info_experiment.csv)

**analysis:** Directory containing the R scripts.

**outputs:** This directory is created when running the scripts. It will contain the processed data and different results including differentially expressed gene tables, annotated Seurat objects, plots...

## Code Execution
#### 1- Download repository
Option 1: Download manually the repository as a ZIP archive and extract it locally on your computer

Option 2: Clone the repository
```shell
git clone https://github.com/saliba-lab/Monocytes_leishmania.git
cd Monocytes_leishmania/analysis
```


#### 2- Install R and python dependencies 
See Dependencies section.


#### 3- Run sequentially the scripts in the **analysis** directory 
Make sure to set the **analysis** directory as the working directory when running the scripts.

## Dependencies
#### Required R libraries
- R                       4.0.3
- dplyr                   1.0.10
- stringr                 1.4.0
- biomaRt                 2.46.3
- tibble                  3.1.8
- ggplot2                 3.4.0
- cowplot                 1.1.1
- scales                  1.2.1
- gridExtra               2.3
- SCnorm                  1.12.1
- SingleCellExperiment    1.12.0
- Seurat                  4.0.0
- harmony                 0.1.1
- igraph                  1.2.6
- leidenbase              0.1.3
- tidyr                   1.2.1
- RColorBrewer            1.1-3
