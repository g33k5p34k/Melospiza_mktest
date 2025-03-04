# McDonald-Kreitman Test for Melospiza melanogenesis genes 
Code for running the McDonald-Kreitman Test on Song Sparrow melanogenesis genes for Oliver Brown et al (in prep). 

This script is custom-written for this specific instance (2 Melospiza melodia maxima individuals and 1 outgroup Melospiza georgiana), but if there is sufficient interest, I can rewrite the script for more general use. 

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
