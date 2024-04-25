# PGS Pipeline
Functionality will support multiple PGS methods


1. Summary statistics preparation
SNP_ID_harmonizer will detect and harmonize the ID between genotype data and GWAS summary stats. the output includes the following columns:

SNPname: original SNP name from the summary stats.


P: p-value.


CHR: chromosome.


BP: position.


A1: Allele1 according to summary statistics.


A2: Allele2 according to summary statistics.


BETA: Effect size (Log(OR), if OR is indicated as the effect size).


NCase: If -Ntype and -Ncase are specified or (Sample size if -Ntype=="quantitative").


CHR_POS: chromosome_position name in case is needed for reference files.


rA1: Allele 1 as specified by the SNP list and harmonized in the name.


rA2: Allele 2 as specified by the SNP list and harmonized in the name.


SNP: SNP name in format CHR_POS_rA1_rA2. <NA> if no matching SNP was found.


SNPchange: If SNP alleles needed to be swapped.*

*Note that these alleles are swapped only for the SNP column (name). Effect size direction is still dictated by A1 and A2 columns.

If -Ntype option is used, the PRS_ready file is also produced
