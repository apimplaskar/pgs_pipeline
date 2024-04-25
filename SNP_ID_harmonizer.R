#!/usr/bin/env Rscript
################################## 
### Function: SNP ID harmonizer for GWAS summary statistics
### date created: 25-Oct-2023
### version 1.2.0
### author: EALM (ealopera@gmail.com)
##################################
## Notes
# This script takes GWAS summary statistics and a list of SNPs as reference
# both need to contain these columns: 
# CHR, POS, ALT, REF
# the target file needs to contain additional Pval and optional SNPname columns
# The harmonizer will create SNP ID column harmonized with the reference list
# and output a file ready to be used by plink
# the harmonize ID will be in format CHR_POS_A1_A2
# Req. columns for input: rCHR, rPOS , rA1, rA2, CHR, POS , A1, A2,
# indicated by column names.
### New
## 17-nov-2023
# added snippet to output the BETA column in the harmonized files
## 24-Apr-2024
# add an output for PRS-estimation software
# add option and variables of sample size | effective sample size | or defined

##################################
#### set up ####
options(stringsAsFactors=F)
Sys.setlocale("LC_CTYPE", "C.UTF-8")

req_packages <- c("optparse", "data.table","tidyverse")
for (pack in req_packages) {
  if(!require(pack, character.only = TRUE)) {
    #install.packages(pack, repos = "https://cloud.r-project.org")
    install.packages(pack, repos='http://cran.us.r-project.org')
  }
}
## bash parser
option_list <- list(
  make_option(c("-r", "--refinput"), type="character", default="",
              help="full path to reference input file"),
  make_option(c("-i","--input"), type="character", default="",
              help="input file to be harmonized"),
  make_option("--rCHR", type="character", default="rchrom",
              help="colnames for chromosome in the reference"),
  make_option("--rA1", type="character", default="rA1",
              help="colnames for allele1 in the reference"),
  make_option("--rA2", type="character", default="rA2",
              help="colnames for allele2 in the reference"),
  make_option("--rPOS", type="character", default="rpos",
              help="colnames for genome position in the reference"),
  make_option("--CHR", type="character", default="chrom",
              help="colnames for chromosome"),
  make_option("--A1", type="character", default="A1",
              help="colnames for allele1"),
  make_option("--A2", type="character", default="A2",
              help="colnames for allele2"),
  make_option("--POS", type="character", default="pos",
              help="colnames for genome position"),
  make_option("--PVAL", type="character", default="P",
              help="colnames for p value of target file"),
    make_option("--effect", type="character", default="BETA",
              help="colnames for effect (beta/BETA or OR/or) "),
  make_option("--secol", type="character", default="SE",
              help="colnames for effect (beta/BETA or OR/or) "),
  make_option("--SNP", type="character", default="SNP",
              help="colnames for SNP identifier  (rs code) in the target file"),
  make_option(c("-p", "--prefix"), type="character", default="",
              help="prefix for output files"),
  make_option("--Ntype", type="character", default="quantitative",
              help="type of phenotype of the summary statistics: binary or quantitative,
              if binary then the Ncase and Ncontrol must be defined"),
  make_option(c("--Ncase","-N"), type="character", default="defined",
              help="sample size or sample size column, if individual
              sample size input number otherwise character"),
  make_option("--Ncontrol", type="character", default="defined",
              help="Control sample size or control sample size column, 
              if individual sample size input number otherwise character")
)
opt_parser  <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)
print(opt)
# test & developing data
#  opt<-list()
# opt$CHR<-"CHR"
# opt$POS<-"POS"
# opt$A2<-"A2"
# opt$A1<-"A1"
# opt$rCHR<-"CHROM"
# opt$rPOS<-"POS"
# opt$rA2<-"ALT"
# opt$rA1<-"REF"
# opt$input<-"/u/project/loes/elopera/pubGWAS_sumstats/processed/1.daner_MDDwoBP_20201001_2015iR15iex_HRC_MDDwoBP_iPSYCH2015i_lifted.txt"
# opt$refinput<-"/u/project/loes/elopera/genotypes/QCedv3/SNPs.list"
# opt$prefix<-"MDD"
# opt$PVAL<-"P"
# opt$SNP<-"SNP"
# opt$effect<-"OR"
##################################
#### main ####
##################################
## define variables
chrcol <- opt$CHR
poscol <- opt$POS
A1col <- opt$A1
A2col <- opt$A2
rchrcol <- opt$rCHR
rposcol <- opt$rPOS
rA1col <- opt$rA1
rA2col <- opt$rA2
prefix<-opt$prefix
input<-opt$input
refinput<-opt$refinput
pvalcol<-opt$PVAL
snpcol<-opt$SNP
effectcol<-opt$effect
secol<-opt$secol
selecv<-c(snpcol,pvalcol,chrcol, poscol, A1col, A2col, effectcol,secol) ### vector of columns to extract
namev<-c("SNPname","P","CHR", "BP", "A1", "A2","BETA","SE") ### names of these columns

##### setting the PRS file #####
#opt<-c()
Ntype<-opt$Ntype
Ncase<-opt$Ncase
Ncontrol<-opt$Ncontrol
suppressWarnings(
  if (!is.null(Ntype)){
    if (is.na(as.numeric(Ncase))){ ###if Ncase is not a number then it will interpreted as column name
      Ncase<-as.character(Ncase)
      selecv<-c(selecv,Ncase)
      namev<-c(namev,"Sample_size")
      if (Ntype=="binary"){
        Ncontrol<-as.character(Ncontrol)
        selecv<-c(selecv,Ncontrol)
        namev[length(namev)]<-"Ncase"
        namev<-c(namev,"Ncontrol")
      }
    } else{
      Ncase<-as.numeric(Ncase)
      if (Ntype=="binary"){
        Ncontrol<-as.numeric(Ncontrol)
      }
    }
  }
  
)


##### Format and retrieve data ######
outfile <- paste0(gsub(".+/([^/]+$)","\\1",prefix),"_SNPharmo")
file_result <- paste0(outfile,".txt")
file_report <- paste0(outfile,"_report.txt")

## import objective file
if( grepl(".gz$",input) | grepl(".bgz$",input) ) {
  dat1_1 = fread(cmd=paste0("gunzip -c ", input), header=T, select=selecv)
} else {
  dat1_1 <- fread(input, header=T, select=selecv)
}
## import reference file
if( grepl(".gz$",refinput) | grepl(".bgz$",refinput) ) {
  datref <-fread(cmd=paste0("gunzip -c ", refinput), header=T, select=c(rchrcol, rposcol, rA1col, rA2col))
} else {
  datref<- fread(refinput, header=T, select=c(rchrcol, rposcol, rA1col, rA2col))
}

print("[INFO] ref and input files imported,  with columns:")
print(names(dat1_1))
print("[INFO] ref and input files imported,  with columns:")
print(names(datref))

## format SNP data
dat1 = dat1_1[, selecv, with=FALSE]
setnames(dat1, selecv,namev)

######## fill the data with sample size if needed

if (!is.null(Ntype)){
  if (!is.na(as.numeric(Ncase))){
    dat1$Sample_size<-Ncase
    if (Ntype=="binary"){
      Ncontrol<-as.numeric(Ncontrol)
      dat1$Ncontrol<-Ncontrol
      names(dat1)[which(names(dat1)=="Sample_size")]<-"Ncase"
    }}
}   
head(dat1)


setnames(datref, c(rchrcol, rposcol, rA1col, rA2col), c("rCHR", "rBP", "rA1", "rA2"))

dat1$A1 = toupper(dat1$A1)
dat1$A2 = toupper(dat1$A2)
datref$rA1 = toupper(datref$rA1)
datref$rA2 = toupper(datref$rA2)
dat1$CHR <- gsub("chr","",dat1$CHR)
datref$CHR <- gsub("chr","",datref$CHR)
dat1$CHR <- as.numeric(gsub("X","23",dat1$CHR))
##format effect column
if (effectcol %like% "OR"| effectcol %like% "or" 
    | effectcol %like% "odds" |effectcol %like% "atio" ) {
  dat1$BETA<-as.numeric(dat1$BETA)
  dat1$BETA<-log(dat1$BETA)
  print(paste0("Odds ratio column ", effectcol, " was transformed into effect size BETA "))
}
#build equal positions column 
dat1$CHR_POS<-paste0(dat1$CHR,"_",dat1$BP)
datref$CHR_POS<-paste0(datref$rCHR,"_",datref$rBP)
#cut reference to the SNPs contained in data only
inref<-which(datref$CHR_POS %in% dat1$CHR_POS )
datref<-datref[inref,]
##remove SNPs without position in original data
dat1<-dat1[which(!is.na(dat1$BP)),]

# bring data from ref
dat1$rA1<-datref$rA1[match(dat1$CHR_POS,datref$CHR_POS)]
dat1$rA2<-datref$rA2[match(dat1$CHR_POS,datref$CHR_POS)]
dat1$SNP<-case_when(dat1$rA1==dat1$A1 & dat1$rA2==dat1$A2~ paste0(dat1$CHR_POS,"_",dat1$A1,"_",dat1$A2),
                        dat1$rA1==dat1$A2 & dat1$rA2==dat1$A1~ paste0(dat1$CHR_POS,"_",dat1$A2,"_",dat1$A1),
                        dat1$rA1==dat1$A1 & dat1$rA2!=dat1$A2~ paste0(dat1$CHR_POS,"_",dat1$A1,"_",dat1$A2,"b"),
                        dat1$rA2==dat1$A1 & dat1$rA1!=dat1$A2~ paste0(dat1$CHR_POS,"_",dat1$A1,"_",dat1$A2,"b"))

dat1$SNPchange<-case_when(dat1$rA1==dat1$A1 & dat1$rA2==dat1$A2~ "no_change",
                        dat1$rA1==dat1$A2 & dat1$rA2==dat1$A1~ "swapped",
                        dat1$rA1==dat1$A1 & dat1$rA2!=dat1$A2~ "possible_triallelic",
                        dat1$rA2==dat1$A1 & dat1$rA1!=dat1$A2~ "possible_triallelic")

report<-table(dat1$SNPchange)
print("finished harmonization, found:")
print(report)
write.table(dat1,file_result,quote=F,sep='\t',row.names = F)
write.table(report,file_report,quote=F,sep='\t',row.names = F)

if (!is.null(Ntype)){
  file_PRSready <- paste0(outfile,"_PRSready.txt")
  if (Ntype=="binary"){
    PRSdat<-dat1[which(!is.na(SNP)),c("SNPname","SNP","CHR","BP","A1","A2", "P","BETA","SE","Ncase","Ncontrol")]
  } else {
    PRSdat<-dat1[which(!is.na(SNP)),c("SNPname","SNP","CHR","BP","A1","A2","P","BETA","SE","Sample_size")]
  }
  write.table(PRSdat,file_PRSready,quote=F,sep='\t',row.names = F)
}
   



print("[INFO] Harmonized file saved in current directory")

###done###
