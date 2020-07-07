# Case study

Reference: CP028309.1
Generate by badreads:
- version: v0.1.5
- param: `simulate --reference {ref path} --quantity 25x  --identity 95,100,5 --start_adapter 0,0 --end_adapter 0,0 --junk_reads 0 --random_reads 0 --chimeras 0 --glitches 0,0,0`

Correction perform with kmer from reference

Read used in file one.fasta.

In this file 150 first base, 200 last base are removed

## Correction with greedy only

First line reads second ref

```
ATG-CGCGATCGGCGCGAGGTTCAACGATCCTGGCGTAACGGTCGCTGACGAATACACCGACACCGTAAACGA--AGACAATCAGCGTCAGCATCGCTTTATCACCGAGCAGACTGCGCAACTCTTTGATACCCAGATTAAAAATATTGCGTAAATGGCGCATCATCCCTCCTGTTTTTTCAGCAGCAGGATACTTAAGCCCATCACCAGCGGGATGGCTATCAGCAACGGGATAAAAAGTTGCCACAAATCAGTCAGATCCAGCGCTTTCGAGAACGTTCCGCGGGCGATAGTCAGAAAATGACTGGTCGGGTAAACCTCGCCGATCCCAACGTCCAGGCCCTTCCAGCGAAGCTACCGGATCGATC-TCCCGGAAA-CTGTGTCGCCGGGATCAACGTGATAATCGCCGTTCCGAAAATGGCGGCGA
||| ||  |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||  |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||| |||||||||||||||||||||||||||||||||||||| ||||||||| |||||||||||||||||||||||||||||||||||||||||||||||| |
ATGGCG--ATCGGCGCGAGGTTCAACGATCCTGGCGTAACGGTCGCTGACGAATACACCGACACCGTAAACGAGAAGACAATCAGCGTCAGCATCGCTTTATCACCGAGCAGACTGCGCAACTCTTTGATACCCAGATTAAAAATATTGCGTAAATGGCGCATCATCCCTCCTGTTTTTTCAGCAGCAGGATACTTAAGCCCATCACCAGCGGGATGGCTATCAGCAACGGGATAAAAAGTTGCCACAAATCAGTCAGATCCAGCGCTTTCGAGAACGTTCCGCGGGCGATAGTCAGAAAATGACTGGTCGGGTAAACCTCGCCGATCC-AACGTCCAGGCCCTTCCAGCGAAGCTACCGGATCGATCATCCCGGAAAACTGTGTCGCCGGGATCAACGTGATAATCGCCGTTCCGAAAATGGCGGC-A
   |  |                                                                  |                                                                                                                                                                                                                                                               |                                      |         |                                                |
   ---|-> 154 Alts with A or C                                           -> 223 One alts A but no valid correction scenario                                                                                                                                                                                                              -> 477 Alts with A or G                ----------|--> 515 corrected                               -> 573 corrected
      |                                                                                                                                                                                                                                                                                                                                                                                   |
      -> Didn't correct because previous error not correct error didn't see                                                                                                                                                                                                                                                                                                               -> 524 corrected
```

bwa mem CIGAR string:
- raw read: 153M1I66M2D254M1I40M1D6M1D203M
- cor read: 153M1I66M2D254M1I251M

## Correction with graph only

First line reads second ref

```
ATG-CGCGATCGGCGCGAGGTTCAACGATCCTGGCGTAACGGTCGCTGACGAATACACCGACACCGTAAACGA--AGACAATCAGCGTCAGCATCGCTTTATCACCGAGCAGACTGCGCAACTCTTTGATACCCAGATTAAAAATATTGCGTAAATGGCGCATCATCCCTCCTGTTTTTTCAGCAGCAGGATACTTAAGCCCATCACCAGCGGGATGGCTATCAGCAACGGGATAAAAAGTTGCCACAAATCAGTCAGATCCAGCGCTTTCGAGAACGTTCCGCGGGCGATAGTCAGAAAATGACTGGTCGGGTAAACCTCGCCGATCCCAACGTCCAGGCCCTTCCAGCGAAGCTACCGGATCGATC-TCCCGGAAA-CTGTGTCGCCGGGATCAACGTGATAATCGCCGTTCCGAAAATGGCGGCGA
||| ||  |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||  |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||| |||||||||||||||||||||||||||||||||||||| ||||||||| |||||||||||||||||||||||||||||||||||||||||||||||| |
ATGGCG--ATCGGCGCGAGGTTCAACGATCCTGGCGTAACGGTCGCTGACGAATACACCGACACCGTAAACGAGAAGACAATCAGCGTCAGCATCGCTTTATCACCGAGCAGACTGCGCAACTCTTTGATACCCAGATTAAAAATATTGCGTAAATGGCGCATCATCCCTCCTGTTTTTTCAGCAGCAGGATACTTAAGCCCATCACCAGCGGGATGGCTATCAGCAACGGGATAAAAAGTTGCCACAAATCAGTCAGATCCAGCGCTTTCGAGAACGTTCCGCGGGCGATAGTCAGAAAATGACTGGTCGGGTAAACCTCGCCGATCC-AACGTCCAGGCCCTTCCAGCGAAGCTACCGGATCGATCATCCCGGAAAACTGTGTCGCCGGGATCAACGTGATAATCGCCGTTCCGAAAATGGCGGC-A
   |  |                                                                  |                                                                                                                                                                                                                                                               |                                      |         |                                                |
   ---|-> 154 not correct k=15 correct k=17                              -> 223 not correct k=15 correct k=17                                                                                                                                                                                                                            -> 477 not corr k=15 corr k=17         ----------|--> 515 not correct k=15 correct k=17           -> 573 not corr k=15 not corr k=17
      |                                                                                                                                                                                                                                                                                                                                                                                   |
      -> correct with 154 ?                                                                                                                                                                                                                                                                                                                                                               -> 524 correct with 515 ?
```

bwa mem CIGAR string:
- raw read: 153M1I66M2D254M1I40M1D6M1D203M
- cor read k=15: 153M1I66M2D254M1I40M1D6M1D203M
- cor read k=17: 726M
