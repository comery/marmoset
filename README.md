## Project 

The common or white-tufted-ear marmoset Callithrix jacchus is a New World primate that inhabits the coastal rainforests of eastern Brazil. Its genome sequence will be useful for comparative genomics studies of primates. The C. jacchus genome is distributed over 22 pairs of autosomes and two sex chromosomes.

This assembly represents the principal haplotype of the diploid genome Callithrix jacchus. This is the VGP trio assembly merged haplotype data from the male and female. The alternate maternal and paternal haplotype sequences are in WGS projects JAALXQ000000000 and JAALXR000000000 respectively.

We applied the trio-binning approach to sequence all three individuals of a trio. The parental samples were sequenced with illumina short reads. And the male offspring was sequenced with PacBio long read, 10X Genomics linked-reads, Bionano optical molecules, and Hi-C reads. We first binned the long-reads into maternal or paternal reads based on the unique k-mer sets from maternal and paternal genomic data. Then the haplotypes were assembled independently using VGP standard genome assembly pipeline which I am not going to explain in detail here. After several rounds of polishing and scaffolding using 10x, Hi-C, and Bionano data, we produce the haploid assemble for the maternal and paternal genomes separately. 

In traditional genome sequencing efforts, heterozygosity was normally estimated by mapping the sequencing reads onto a mosaic reference genome, resulting in limited phase information of the heterozygous variants. Our chromosome-level diploid assemblies allow us to directly compare the two parentally inherited genomes and identify the full spectrum of genetic variants between the parental alleles in de novo way. These genetic variations include not only the single nucleotide variations (SNVs), insertion/deletions (indels), but also other large structure variations (SVs). 


## Analysis material

This repository collects most customized scripts for genetic diveristy analysis for marmset project, including genetic diversity analysis, sequencing & polishing error measurement, run of homozygosity, and chimerism analysis.

- Heterozygosity spectrum between paternal and maternal genomes
	- [Mummer alignment, SNV and Small indels](Mummer_alignment/README.md)
	- [10x_and_Pacbio_mapping](10x_and_Pacbio_mapping/README.md)
	- [Sequencing and polishing error measurement](Sequencing_and_polishing_error/README.md)
	- [SV detection](Structure_variations/README.md)
	- [PCR validation for SNVs and small indels](Experimental_validation/README.md)
	- [Chimerism](Chimerism/README.md)
	- [ROH](ROH/README.txt)

- Marmoset biology
	- [Positive analysis](positive_selection/README.md)







