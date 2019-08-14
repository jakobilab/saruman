# SARUMAN

## Note: this repository serves archival purposes

SARUMAN was developed in 2009 with CUDA version 4. The source code is now on GitHub, however, while I was able to compile it with recent CUDA versions (>=10), I cannot give any guaranty that the software a) works or b) produces still correct results. 

### Introduction

Using GPU programming for short read mapping  
  
Since the introduction of next generation sequencing technologies like Solexa, 454, and SOLiD the amount of generated data rises with each new technology upgrade. As the application scenarios especially of the short read techniques include the re-sequencing of known genomes or sequencing of closely related strains, new software tools are needed for the fast mapping of sequencing reads against a reference genome.  

Currently, there are several tools available, but most of them are limited either in speed or accuracy. Limitations in accuracy lead to non detected mappings, which could become important in post processing steps like SNP calling. Because of those limitations our goal was to develop an exact and complete mapping algorithm with equivalent running time compared to available heuristic implementations.  
  

The result is SARUMAN (**S**emiglobal **A**lignment of Short **R**eads using C**U**DA and Needle**man**-Wunsch). SARUMAN uses a qgram index based filter algorithm followed by a modified Needleman-Wunsch alignment. To speed up he normally time-consuming alignment step all alignments are processed on a NVIDIA graphics card to exploit the massive parallel architecture of new graphics processing units (GPUs). Based on this technique, depending on the input read length, SARUMAN is able to process hundreds of thousands of alignments in just a few seconds. As a result of this alignment strategy SARUMAN not only detects mismatches, but also allows to detect and handle all insertions and deletions correctly. The mapping algorithm is exact and complete, it identifies all possible matching positions for a given error threshold and always returns the optimal local alignment.

### Requirements

Before downloading and testing SARUMAN please make sure that your own system meets the following requirements:

- OS:	64bit Linux system (SARUMAN was tested on Ubuntu & Gentoo) 
- Hardware: 4-8GB RAM, dual core CPU, GPU with at least 512MB for reasonable performance 
- [**CUDA compatible graphics card**](http://www.nvidia.com/object/cuda_gpus.html) 
- [**CUDA capable driver**](http://developer.nvidia.com/object/cuda_download.html) for your card 
- [**CUDA runtime environment**](http://developer.nvidia.com/object/cuda_download.html) for proper functioning of the CUDA module 
- Two additional libraries, both available as installation package for almost all linux distributions:  
- [**argtable2**](http://argtable.sourceforge.net/)  
- [**uthash.h**](http://uthash.sourceforge.net/) 
- A working [**BioPerl**](http://www.bioperl.org/) installation for converting the SARUMAN output into SAM format 

### Downloads

The download package contains the SARUMAN Linux binaries tested on Ubuntu and Gentoo Linux, a Perl script for converting the original SARUMAN output into SAM format and a short documentation with installation instructions and commandline options.

[**Download SARUMAN**](ftp://ftp.cebitec.uni-bielefeld.de/pub/software/saruman/saruman-current.tar.bz2)   
  
[**E.Coli K12 artificial sample data**](ftp://ftp.cebitec.uni-bielefeld.de/pub/software/saruman/saruman_sample_data_artificial.tar.bz2)

[**Corynebacterium Glutamicum Illumina data**](ftp://ftp.cebitec.uni-bielefeld.de/pub/software/saruman/saruman_sample_data_CG.tar.bz2)

### Publication

If you use SARUMAN please cite the following publication:

Exact and complete short read alignment to microbial genomes using GPU programming 
Jochen Blom, Tobias Jakobi, Daniel Doppmeier, Sebastian Jaenicke, Jorn Kalinowski, Jens Stoye, and Alexander Goesmann

Bioinformatics published 30 March 2011, 10.1093/bioinformatics/btr151 [**http://bioinformatics.oxfordjournals.org/cgi/content/abstract/btr151v1?papetoc**](http://bioinformatics.oxfordjournals.org/cgi/content/abstract/btr151v1?papetoc) 

[**Download preprint manuscript**](ftp://ftp.cebitec.uni-bielefeld.de/pub/software/saruman/publication.pdf)

### Licence

SARUMAN is free for non-commercial use.

Commercial users: please contact [**tobias@jako.bi**](mailto:tobias@jako.bi).

