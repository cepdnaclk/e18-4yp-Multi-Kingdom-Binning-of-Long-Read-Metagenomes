---
layout: home
permalink: index.html

# Please update this with your repository name and title
repository-name: eYY-4yp-project-template
title:
---

[comment]: # "This is the standard layout for the project, but you can clean this and use your own template"

# Long-reads Binning For Microbial Metagenomics Considering Multi-kingdoms


#### Team

- E/18/030, Aththanayaka A.M.S., [e18030@pdn.ac.lk](mailto:name@email.com)
- E/18/282, Ranasinghe R.A.N.S., [e18282@pdn.ac.lk](mailto:name@email.com)
- E/18/283, Ranasinghe R.D.J.M., [e18283@pdn.ac.lk](mailto:name@email.com)

#### Supervisors

- Dr. Damayanthi Herath, [damayanthiherath@eng.pdn.ac.lk](mailto:name@eng.pdn.ac.lk)
- Dr. Vijini Mallawarachchi, [vijini.mallawaarachchi@flinders.edu.au](mailto:name@eng.pdn.ac.lk)

#### Table of content

1. [Abstract](#abstract)
2. [Background](#background)
3. [Related works](#related-works)
4. [Methodology](#methodology)
5. [Experiment Setup and Implementation](#experiment-setup-and-implementation)
6. [Results and Analysis](#results-and-analysis)
7. [Conclusion](#conclusion)
8. [Publications](#publications)
9. [Links](#links)


## Abstract

## Background 

<p align="center">
  <img src="https://useruploads.socratic.org/puDVmAVgSBy1dqrzo38g_cellsToDNA.gif" alt="Alt text" width="500" height="320">
</p>

<p style="text-align: justify">
Microorganisms thrive in a multitude of environments worldwide, fulfilling critical roles in human health, agriculture, food production, climate regulation, and numerous other processes. Every living organism consists of tiny units known as cells, serving dual functions of providing structure and facilitating various biological processes. Enclosed within the nucleus of each cell lies the genome, a comprehensive blueprint encompassing instructions for the construction and sustenance of the entire organism, including its distinct characteristics and behaviors. This genetic blueprint resides within slender, thread-like structures called chromosomes, composed of DNA and proteins. DNA, the carrier of genetic information, adopts a double helix structure comprising two intertwined strands. Comprised of nucleotides, each denoted by specific letters—A (Adenine), C (Cytosine), G (Guanine), or T (Thymine)—DNA serves as the foundation for genetic coding. Genes, the fundamental units of heredity, constitute segments of DNA containing instructions for synthesizing proteins or functional RNA molecules. Serving as conduits of hereditary information, genes perpetuate traits across generations, thereby ensuring the perpetuation of life.
</p>

## Related works
<p style="text-align: justify">
<b>Early long-reads Binning Tools<br></b>
Megan-LR stands out as one of the earliest tools, employing a reference database. Megan-LR utilizes a protein-alignment-based approach and introduces two algorithms; one for taxonomic binning (based on Lowest Common Ancestor) and another for functional binning (based on an Interval-tree algorithm).
</p>
<p style="text-align: justify">
Two other noteworthy reference-independent tools, MetaProb and BusyBee Web, significantly contributed to the domain of unsupervised metagenomic binning. BusyBee Web, in particular, includes a web-based interface, offering additional visual insights into the binning process. However, despite their respective strengths, both MetaProb and BusyBee Web faced challenges related to scalability as input dataset sizes increased, impeding their ability to bin entire datasets in a single iteration.
</p>

<p align="justify">
<b>MetaBCC-LR<br></b>
MetaBCC-LR, a reference-free binning tool, utilizes composition and coverage as read features, relying on trinucleotide frequency vectors for composition and k-mer coverage histograms for coverage. The tool initially clusters reads based on coverage information, which will be re-clustered using composition information. 
Only a sample of reads is utilized for this process, contributing to computational efficiency. At the final stage, it creates statistical models for each cluster and bin the remaining reads. Despite its high accuracy, it may suffer from potential misclassification issues, particularly for low-abundance species, as well as the need for subsampling large datasets.
</p>

<p style="text-align: justify">
<b>LRBinner<br></b>
LRBinner adopts an innovative approach to reference-free binning by concurrently computing composition and coverage information for the entire dataset. It merges these features through a variational autoencoder, eliminating the need for subsampling and improving overall binning accuracy. It uses tetranucleotide frequency vectors for composition and k-mer coverage vectors as coverage information of reads. However, the tool faces challenges in distinguishing long reads from similar regions shared between different species.
</p>

<p style="text-align: justify">
<b>OBLR<br></b>
OBLR introduces a novel strategy in reference-free binning, leveraging read overlap graphs to estimate coverages and improve binning outcomes. It then employs the HDBSCAN hierarchical density-based clustering algorithm for read clustering. Additionally, it uses a sample of reads for initial clustering sampled using a probabilistic downsampling strategy. This results in clusters with similar sizes and fewer isolated points. OBLR then utilizes inductive learning with the GraphSAGE neural network architecture to assign bins to remaining reads.
</p>

<div align="center">
  <img src="./images/metabcc.png" alt="Workflow MetaBCC-LR" width="30%" title="Workflow MetaBCC-LR">
  <img src="https://media.springernature.com/full/springer-static/image/art%3A10.1186%2Fs13015-022-00221-z/MediaObjects/13015_2022_221_Fig1_HTML.png?as=webp " alt="Workflow LRBinner" width="30%" title="Workflow LRBinner">
  <img src="https://media.springernature.com/lw685/springer-static/image/chp%3A10.1007%2F978-3-031-06220-9_15/MediaObjects/526061_1_En_15_Fig1_HTML.png" alt="Workflow OBLR" width="30%" title="Workflow LRBinner">
</div>

## Proposed Work
We have identified the following as the challenges in existing tools.

- Mainly focus on composition and coverage as primary features. However, marker genes-based kingdom-level information can enhance the binning process.

- Existing long reads binning tools overlook differential abundance in multiple samples. Considering species abundance across samples could enhance binning accuracy.

- Lack of binning refinements for long reads binning tools. Introducing refining mechanisms could improve the precision of bin assignments.

Therefore, this project aims to develop a method to bin long reads from multiple metagenomic samples while being aware of the underlying microbial kingdoms. Specifically, it will be a Python-based command-line tool addressing the scalability issues with massive datasets.

## Methodology

## Experiment Setup and Implementation

## Results and Analysis

## Conclusion

## Publications
[//]: # "Note: Uncomment each once you uploaded the files to the repository"

<!-- 1. [Semester 7 report](./) -->
<!-- 2. [Semester 7 slides](./) -->
<!-- 3. [Semester 8 report](./) -->
<!-- 4. [Semester 8 slides](./) -->
<!-- 5. Author 1, Author 2 and Author 3 "Research paper title" (2021). [PDF](./). -->


## Links

[//]: # ( NOTE: EDIT THIS LINKS WITH YOUR REPO DETAILS )

- [Project Repository](https://github.com/cepdnaclk/repository-name)
- [Project Page](https://cepdnaclk.github.io/repository-name)
- [Department of Computer Engineering](http://www.ce.pdn.ac.lk/)
- [University of Peradeniya](https://eng.pdn.ac.lk/)

[//]: # "Please refer this to learn more about Markdown syntax"
[//]: # "https://github.com/adam-p/markdown-here/wiki/Markdown-Cheatsheet"
