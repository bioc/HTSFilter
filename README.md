# HTSFilter: Data-based filtering for replicated transcriptome sequencing experiments

Authors: Andrea Rau, M&eacute;lina Gallopin, Gilles Celeux, Florence Jaffr&eacute;zic
 
High-throughput sequencing (HTS) data, such as RNA-sequencing (RNA-seq) data, are often used to conduct differential analyses, in which statistical tests are performed for each biological feature (e.g., a gene, transcript, exon) in order to identify those whose expression levels show systematic covariation with a particular condition, such as a treatment or phenotype of interest. Because hypothesis tests are performed for gene-by-gene differential analyses, the obtained *p*-values must be adjusted to correct for multiple testing. However, procedures to adjust *p*-values to control the number of detected false positives often lead to a loss of power to detect truly differentially expressed (DE) genes due to the large number of hypothesis tests performed. To reduce the impact of such procedures, independent data filters are often used to identify and remove genes that appear to generate an uninformative signal; this in turn moderates the correction needed to adjust for multiple testing. 

The *HTSFilter* package (DOI: 10.18129/B9.bioc.HTSFilter) implements a novel data-based filtering procedure based on the calculation of a similarity index among biological replicates for read counts arising from replicated transcriptome sequencing (RNA-seq) data. This technique provides an intuitive data-driven way to filter high-throughput transcriptome sequencing data and to effectively remove genes with low, constant expression levels without incorrectly removing those that would otherwise have been identified as DE. 

### Reference

Rau, A., Gallopin, M., Celeux, G., and Jaffrézic, F. (2013). Data-based filtering for replicated high-throughput transcriptome sequencing experiments. *Bioinformatics* 29(17): 2146-2152. 

### License

The *HTSFilter* package is free software; you can copy or redistribute it under the terms of the Artistic-2.0 License. This program is distributed in the hope that it will be useful, but without any warranty; without even the implied 
warranty of merchantability or fitness for a particular purpose. See the Artistic-2.0 License for more details.

A copy of the Artistic License, version 2.0, is available at https://opensource.org/licenses/Artistic-2.0.