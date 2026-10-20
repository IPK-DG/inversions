# Detection of large chromosomal inversions with chromsome conformation capture sequencing (Hi-C) data — a usage example

## Usage example

1. Load the package `data.table`. This code has been tested with version 1.9.6 and 1.10.4 of `data.table`. More information on `data.table` can be found [here](https://cran.r-project.org/web/packages/data.table/).

    ```r
    library(data.table)
    ```

2. Read the R object holding the Hi-C link information for Morex and Barke from a file. This file can be downloaded from the [Plant Genomics and Phenomics Research Data Repository](https://doi.ipk-gatersleben.de/DOI/c7305c69-6b95-46a0-97ca-2ad6df8cc343/3f9f60f5-10b1-4b0a-a26c-a7684947a4c5/2/1847940088).

    ```r
    readRDS('hic_links_morex_barke.Rds') -> hic_links_morex_barke
    ```

3. Read the pseudomolecules lengths of the reference genome sequence assembly of barley cultivar Morex. This file is available in a public [Bitbucket repository](https://bitbucket.org/inversions/inversions.bitbucket.io/src).

    ```r
    fread('barley_chromosome_lengths.txt') -> chromosome_lengths
    ```

4. Read the function `findInversions()` from a file. This file is available in a public [Bitbucket repository](https://bitbucket.org/inversions/inversions.bitbucket.io/src).

    ```r
    source('findInversions.R')
    ```

5. Run the inversion detection by calling `findInversions()`.

    ```r
    morex_barke_inversions <- findInversions(links=hic_links_morex_barke,
                                             genotypes=c("Morex", "Barke"),
                                             chromosome_lengths=chromosome_lengths,
                                             maxdist=1e8, multiplier=1e8, cab=c(1,2),
                                             rsquare=0.4, coefficient=-2, dcab=1.8)
    ```

    `findInversions()` takes the following parameters:

    - **links**: A data.table in [long format](http://www.cookbook-r.com/Manipulating_data/Converting_data_between_wide_and_long_format/) with the columns as in the example `hic_links_morex_barke`.

    - **genotypes**: A character vector with exactly two elements listing the genotypes to compare. The genotypes should occurs in `genotype` column of `links`.

    - **chromosome_lengths**: A two-column data.table indicating the length (in bp) of each chromosome.

    - **maxdist**: Hi-C links between genomic bins separated by more than `maxdist` bp will be ignored.

    - **multiplier**: Effectively a depth threshold: a multiplier of `1e8` (100 million) means that pairs of genomic bins with less than 1 link per 100 total million links will be discarded.

    - **cab**: Numeric vector with two elements indicating the `CAB` values of genotypes 1 and 2 at boundaries of candidate inversions as described in Himmelbach *et al.* 2017. Smaller values and of the first element or higher values of the second element increase stringency.

    - **rsquare**: Minimum adjusted `r2` value of the linear model described in Himmelbach *et al.*, 2017. Higher values will result in more stringent filtering.

    - **coefficient**: Maximum coefficient of the linear model described in Himmelbach *et al.*, 2017. Smaller values impose a more stringent filter.

    - **dcab**: Minimum absolute difference between `CAB` of genotypes 1 and 2 at boundaries of inversions. Higher values increase stringency.

6. Plot the directionality biases along the genome and highlight putative inversions.

    ```r
    pdf("Inversions_between_Morex_and_Barke.pdf", height=10, width=7)
    par(mfrow=c(7,1))
    par(oma=c(4,0,0,0))
    par(mar=c(0,4,5,1))
    ylim <- 300
    with(morex_barke_inversions, lapply(1:7, function(i){
     directionality_bias[chr == i, plot(bin/1e6, directionality_bias,
                                        pch=20, las=1, col=1, ylim=c(-ylim, ylim),
                                        bty='l', ylab="D", xlab="")]
     if(nrow(putative_inversions[chr == i]) > 0)
      putative_inversions[chr == i, rect(col="#FF000011", border='red', start/1e6, -2*ylim, end/1e6, 2*ylim)]
     title(paste0(i, "H"), line=1.4)
    }))
    title(outer=T, xlab="genomic position (Mb)", xpd=NA)
    dev.off()
    ```
