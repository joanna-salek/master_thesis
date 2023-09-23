# Master thesis - Searching for evolutionary relationships based on the  chaos game representation
## This is a project for Master Thesis in Bioinformatics

### 1. The aim of this project is to find better solutions for biological sequence comparinsion 

One of the main interests of bioinformatics research is the effective comparison of the biological sequences. This allows studying sequences similarity, and thus forms the basis of phylogenetics. There are several popular sequence comparison approaches which are usually based on dynamic programming techniques. Despite the high accuracy of sequence alignment, it is ineffective in analysis of large set of long sequences. Nowadays, we observe the increase in the number of biological data, mainly because of next generation sequencing methods. Therefore, “free-alignment” methods have recently gained a lot of interest. Some of them transform sequences into a numerical series.

### 2. Methods used: 
##### First step was to translate DNA sequence to sequence of numbers (euclidean coordinates) based on chaos game representation and then searching for relationships between sequences.
##### Second step was to construct phylogenetic trees based on this relationships and compering them to classic phylogenteics methods. 

In this study, three types of nucleotide sequences were analyzed: the SPARC gene of vertebrate representatives, the neuraminidase gene of influenza viruses and the genomes of coronaviruses. The Hurst exponent, Discrete Fourier Transform (DFT), descriptors of the moment of inertia and the mixed method were used to analyze the numerical series, which represent sequences. Based on the results, phylogenetic trees were constructed using the Unweighted Pair Group Method with Arithmetic mean (UPGMA). To verify the reliability of these trees, literature data and classical approaches were employed, i.e. based on alignment (Clustal Omega), and k-mers and Lempel-Ziv complexity.

### 3. Conclutions:
The Discrete Fourier Transform provided the most reliable phylogenetic trees, especially when it comes to coronaviruses and influenza viruses. Remaining methods turned out to be insufficient, which was revealed during phylogenetic trees inspection. Using reference methods the results remained in consistent with those from DFT. 
