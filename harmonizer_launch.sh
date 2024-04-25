wkdir="/u/project/loes/elopera/genotypes/QCedv3/" ### directory with plink files
scriptdir="/u/project/loes/elopera/scripts/" ##directory with the scripts
GWASdir="/u/project/loes/elopera/pubGWAS_sumstats/processed" ### directory with the GWAS summary statistics
name="HEIGHT" ## Phenotype
target="5.HEIGHT_lifted.txt" ## summary statistics file
Ntype="quantitative" ## "binary" or "quantitative" according to phenotype. will decide if sample size is one or 2 columns. leave empty if unknown or not necessary
Ncase=2197673 ### number if GWAS has a single sample size or Column name if each SNP has individual sample sizes. leave empty if unknown or not necessary
Ncontrol="Nco" ### 

#### ### create SNP list from the plink files
cat ${wkdir}/*bim > ${wkdir}/SNPs.list
head ${wkdir}/SNPs.list
## add header to the all snps lits
sed -i '1i\CHROM\tID\tCM\tPOS\tREF\tALT' ${wkdir}/SNPs.list

cd $GWASdir

Rscript ${scriptdir}/SNP_ID_harmonizer.R \
--input="${GWASdir}/${target}" \
--refinput="${wkdir}/SNPs.list" \
--CHR="CHR" --POS="POS" --A1="A1" --A2="A2" --PVAL="P" --SNP="SNP" --effect="OR" \
--rCHR="CHROM" --rPOS="POS" --rA2="ALT" --rA1="REF" \
--Ntype="$Ntype" --Ncase="$Ncase" \
--prefix="$name"

#done
