# CasPhi-guide-targeting
This script takes a list gene fasta sequences "Genes.fasta" and generates:

1. A list of possible spacers for targeting each gene based on the "TBN" Pam sequence of CasPhi.The spacer positions are picked such that the distance of each spacer from the preceding spacer is > 10% of the gene length. The GC% of the spacer region is also required to be 40-80%.
2. For each spacer it creates variants of target sequence by varying the PAM sequence or by creating mismatches between the spacer (Guide RNA) and the target.
3. A final list of oligos that contains contant primer binding sites on 5' and 3' ends that flank the repeat, spacer, barcode and the target sequence (WT, PAM or mismatch variants)
