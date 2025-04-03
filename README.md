# Investigating genomic landscape of influenza virus A in Swiss wastewater

Influenza A virus poses significant public health challenges, causing seasonal outbreaks and pandemics. Its rapid evolution motivates continuous monitoring of circulating influenza genomes to inform vaccine and antiviral development. Wastewater-based surveillance offers an unbiased, cost-effective approach for genomic surveillance. We developed a tiling amplicon primer panel for wastewater-based genomic surveillance of influenza A virus, targeting segments of the surface proteins HA, NA, and M of subtypes H1N1 and H3N2. Using this panel, we sequenced nucleic acid extracts from 59 Swiss wastewater samples collected at four locations during the 2022/2023 and 2023/2024 winter seasons.

**In this repo we present the bioinformatic tools and configurations needed to generate the results presented in our publication xxx. Furthur we provide the scripts needed to generate the figures presented in the publication.**

## Main analysis steps

### 1. V-pipe

Raw sequencing reads were processed using the bioinformatics pipeline [**V-pipe** ](https://github.com/cbg-ethz/V-pipe) (Version 3.0). The pipeline included the main steps of primer trimming using samtools (Version 1.19), alignment using bwa (Version 0.7.17), mutation calling using LoFreq (Version 2.1.3) and, where applicable, local haplotype reconstruction using VILOCA (Version 1.0.0). 

<ins>**1.1 Set-up**</ins>
> - For each segment subtype a v-pipe working directory is set up ``work-IA_XX/v-pipe``

<ins>**1.2 Inputs**</ins>

**i)** Configuration files
> - specifies settings for different steps within pipeline
> - for each data set/ analysis discussed in the publication we provide a configuration file template (see work ``work-IA_XX/v-pipe/v-pipe_configFiles``)
  
**ii)** Raw fastq files
> - In our study, we sequenced using MinION Mk1C, basecalling was performed using the Dorado pipeline (Version 7.3.11)
> - Raw data needs to be foramted as ``Raw_data_dirName/barcode_xx/batch_yy/raw_data/raw_readsName.fastq.gz`` (``barcode_xx`` = physical sample, ``batch_yy``= batch sample belongs to, more infos under [**V-pipe** ](https://github.com/cbg-ethz/V-pipe))
> - Path to raw data needs to be specified in configuration file

**iii)** Samples file
> - Specifies which samples (``barcode_xx/batch_yy``) to process (see, ``work-IA_XX/v-pipe/samples_zz``)

<ins>**1.3 Outputs**</ins>

For each data set and the corresponding analysis, a results directory with the structure ``work-IA_XX/v-pipe/results/barcode_xx/batch_yy`` is generated. Per sample (``barcode_xx``) they contain various (intermediate) outputs of the pipeline. As the entire output can become memory heavy we only provide selected example.


### 2. Clade abundance estimation

Using the output of the V-pipe alignment step, LolliPop (Version 0.5.1) was applied to estimate the relative abundance of H1N1 clades, in the time-series wastewater dataset.

### 3. Converting Nucleotides to amino acids


