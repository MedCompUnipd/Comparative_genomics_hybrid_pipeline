## 🧬Hybrid Comparative Genomics Pipeline

This repository provides an automated pipeline for comparative genomic analysis using a hybrid approach integrating:
- Oxford Nanopore Technologies (ONT) long-read DNA sequencing
- Illumina RNA-seq short-read transcriptome data

It streamlines analysis from raw signal to polished genome assemblies, including comparison and evaluation metrics.

## 💻 Features
This pipeline includes scripts and configurations for:

- **Conversion**: fast5 to pod5 using pod5tools

- **Basecalling**: performed using Dorado

- **Taxonomic analysis**: using Kraken2

- **Quality Control**: Nanoplot, FastQC, Qualimap

- **Trimming & Filtering**: fastp, seqtk

- **De novo Assembly**: Flye, Canu, wtdbg2, NextDenovo, Hifiasm
  
- **Comparison & Evaluation**: QUAST, MUMmer, Minimap2, BLAST, GLSEARCH, BCFtools
    
- **Polishing**: Racon, Medaka, BWA, Pilon


## ✅ Requirements
- Bash (Linux environment)
- Conda or Miniconda installed
- Compatible bioinformatics tools installed and accessible in your $PATH (**read the SetupEnvironments file**)


## 🚀 Running the Pipeline
The pipeline is controlled by a master bash script which runs all analysis steps and logs progress:
bash 0_master.sh

## 📊 Outputs
- Polished genome assemblies and their stats and reports;
- Alignment to reference genes and RBH summary Excel and TSV files;
- Nucleotidic and aminoacidic sequences from denovo assembly FASTA file;
- Alignment logs and summary reports.

## 📚 Additional Notes
For reference genomes and annotation data, download FASTA files from NCBI and GFF3 files EMBL-EBI:
- NCBI Genomes
- Use tools like wget or Entrez Direct for batch downloads.
- Visualize your results with genome browsers like IGV or Artemis.

## 🧠 Citation
If you use this pipeline in your research, please cite the relevant tools listed above and credit this repository.

## 📫 Contact
For questions or contributions, please open an issue or submit a pull request.
