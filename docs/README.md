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
<p style="text-align: justify">
DNA metagenomics, which analyzes the entire genetic pool of an environmental sample, offers powerful insights into microbial communities. Traditionally, short-read sequencing technology dominated metagenomic analysis. As sequencing technology advanced, long-read sequencing emerged, generating significantly longer reads. Then several binning tools have developed enabling reconstruction of more complete genomes. Most of these tools have used coverage and composition features for binning procedure and have achieved good accuracy. 


This research introduces [G11 refiner], a novel long-read binning refiner designed to address further additional read features like kingdom level information of microorganisms to enhance the accuracy and work along with long read binning tools like OBLR, MetaBCC-LR. By incorporating these advancements, [G11 refiner] aims to significantly improve the accuracy and efficiency of binning long-reads by using multi-kingdom data. 

</p>

## Background 

<p align="center">
  <img src="https://useruploads.socratic.org/puDVmAVgSBy1dqrzo38g_cellsToDNA.gif" alt="Alt text" width="500" height="320">
</p>

<p style="text-align: justify">
Our planet harbors a hidden universe teeming with life forms invisible to the naked eye – microorganisms. These diverse single-celled organisms, including bacteria, archaea, protists, and fungi, exist in a multitude of environments, from the scorching deserts to the frigid depths of the oceans. Despite their minute size, microorganisms play a critical role in the intricate tapestry of life on Earth.


Every living organism, from towering trees to microscopic bacteria, is built from fundamental units called cells. These microscopic marvels serve a dual purpose: providing structure and carrying out the essential chemical reactions that sustain life. Tucked away within the cell's nucleus lies the blueprint for the entire organism – its genome. This blueprint dictates everything from physical appearance to specialized functions.


The genome itself is meticulously organized into thread-like structures called chromosomes, composed of DNA and proteins. DNA, the molecule of heredity, boasts a unique double helix structure. Each rung of this ladder is formed by a pair of molecules called nucleotides, identified by specific letters: A, C, G, and T. The specific sequence of these nucleotides along the DNA strand acts as a coded message that dictates the production of proteins, the workhorses of the cell. These proteins perform a multitude of tasks, ensuring the proper functioning of the organism. Genes, discrete segments along the DNA strand, house the instructions for specific proteins or functional RNA molecules. These genes are the cornerstone of heredity, faithfully passing traits from parents to offspring, guaranteeing the continuation of life in all its magnificent diversity.  


While most microorganisms are single-celled entities, viruses occupy a unique niche. These infectious agents are much smaller and simpler than cells, lacking the complex machinery for independent life. Viruses consist of genetic material (either DNA or RNA) enclosed in a protein coat. They rely on hijacking the cellular machinery of host organisms to replicate and spread.  Despite their parasitic nature, viruses play a significant role in ecosystems, influencing the evolution of their hosts and participating in nutrient cycling.   
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

<img src="./images/metabcc.png" alt="Workflow MetaBCC-LR" width="50%" title="Workflow MetaBCC-LR"><br>
<p style="text-align: justify">
<b>LRBinner<br></b>
LRBinner adopts an innovative approach to reference-free binning by concurrently computing composition and coverage information for the entire dataset. It merges these features through a variational autoencoder, eliminating the need for subsampling and improving overall binning accuracy. It uses tetranucleotide frequency vectors for composition and k-mer coverage vectors as coverage information of reads. However, the tool faces challenges in distinguishing long reads from similar regions shared between different species.
</p>

<img src="https://media.springernature.com/full/springer-static/image/art%3A10.1186%2Fs13015-022-00221-z/MediaObjects/13015_2022_221_Fig1_HTML.png?as=webp " alt="Workflow LRBinner" width="50%" title="Workflow LRBinner"><br>

<p style="text-align: justify">
<b>OBLR<br></b>
OBLR introduces a novel strategy in reference-free binning, leveraging read overlap graphs to estimate coverages and improve binning outcomes. It then employs the HDBSCAN hierarchical density-based clustering algorithm for read clustering. Additionally, it uses a sample of reads for initial clustering sampled using a probabilistic downsampling strategy. This results in clusters with similar sizes and fewer isolated points. OBLR then utilizes inductive learning with the GraphSAGE neural network architecture to assign bins to remaining reads.
</p>

<img src="https://media.springernature.com/lw685/springer-static/image/chp%3A10.1007%2F978-3-031-06220-9_15/MediaObjects/526061_1_En_15_Fig1_HTML.png" alt="Workflow OBLR" width="50%" title="Workflow OBLR"><br>

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
