# McDonald-Kreitman Test for *Melospiza* melanogenesis genes 
Code for running the McDonald-Kreitman Test (McDonald & Kreitman, 1991) on Song Sparrow melanogenesis genes for Oliver Brown et al (in prep). 

This script is custom-written for this specific instance (2 *Melospiza melodia maxima* individuals and 1 outgroup *Melospiza georgiana*), but if there is sufficient interest, I can rewrite the script for more general use. 

The core math for this script was adapted from the approach implemented by Li and Fu (2003) in NeutralityTest. 

This script requires codon-aligned multi-sequence FASTA files with 3 entries (named maxima, maxima_UAM, and georgiana). 

To run the McDonald-Kreitman Test, use the function:

```
calc_melospiza_MK(fastafile)
```

This function generates a data frame with the following variables
* Number of synonymous divergent sites (Ds)
* Number of non-synonymous divergent sites (Dn)
* Number of synonymous polymorphic sites (Ps)
* Number of non-synonymous polymorphic sites (Pn)
* The neutrality index
* Alpha (proportion of substitutions driven by positive selection)
* P-value of the Fisher's Exact Test

The final code block at the end runs the MK test for a defined set of genes and spits out a CSV file with the raw MK test data. 

# References
* Li, H., & Fu, Y.-X. (2003). NeutralityTest: novel software for performing tests of neutrality. https://www.picb.ac.cn/evolgen/softwares/download/NeutralityTest/MKtest.htm
* McDonald, J. H., & Kreitman, M. (1991). Adaptive protein evolution at the Adh locus in Drosophila. Nature, 351(6328), 652–654.
