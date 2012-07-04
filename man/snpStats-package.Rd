\name{snpStats-package}
\alias{snpStats-package}
\alias{snpStats}
\docType{package}
\title{
SnpMatrix and XSnpMatrix classes and methods
}
\description{
Classes and statistical methods for large SNP association studies, extending the snpMatrix package (now removed)
}
\details{
\tabular{ll}{
Package: \tab snpStats\cr
Version: \tab 1.7.2\cr
Date: \tab 2012-07-04\cr
Depends: \tab R(>= 2.10.0), survival, methods, Matrix\cr
Imports: \tab graphics, grDevices, methods, stats, survival, utils, Matrix\cr
Suggests: \tab hexbin\cr
License: \tab GPL-3\cr
Collate: \tab ss.R contingency.table.R convert.R compare.R glm-test.R imputation.R indata.R long.R misc.R ld.R mvtests.R pedfile.R outdata.R plink.R qc.R qq-chisq.R single.R tdt-single.R structure.R xstuff.R zzz.R\cr
LazyLoad: \tab yes\cr
biocViews: \tab Microarray, SNP, GeneticVariability\cr
Packaged: \tab 2012-07-04\cr
Built: \tab R 2.13.1; i686-pc-linux-gnu; 2012-07-04 11:19:01 UTC; unix\cr
}

Index:
\preformatted{
Fst                     Calculate fixation indices
GlmEstimates-class      Class "GlmEstimates"
GlmTests-class          Classes "GlmTests" and "GlmTestsScore"
ImputationRules-class   Class "ImputationRules"
SingleSnpTests-class    Classes "SingleSnpTests" and
                        "SingleSnpTestsScore"
SnpMatrix-class         Class "SnpMatrix"
XSnpMatrix-class        Class "XSnpMatrix"
chi.squared             Extract test statistics and p-values
convert.snpMatrix       Convert 'snpMatrix' objects to 'snpStats'
                        objects
example-new             An example of intensity data for SNP genotyping
families                Test data for family association tests
filter.rules            Filter a set of imputation rules
for.exercise            Data for exercise in use of the snpStats
                        package
glm.test.control        Set up control object for GLM computations
ibsCount                Count alleles identical by state
ibsDist                 Distance matrix based on identity by state
                        (IBS)
imputation.maf          Extract statistics from imputation rules
impute.snps             Impute snps
ld                      Pairwise linkage disequilibrium measures
ld.example              Datasets to illustrate calculation of linkage
                        disequilibrium statistics
mean2g                  Raw coding of posterior probabilities of SNP
                        genotype
misinherits             Find non-Mendelian inheritances in family data
mvtests                 Multivariate SNP tests
plotUncertainty         Plot posterior probabilities of genotype
                        assignment
pool                    Pool test results from several studies or
                        sub-studies
pool2                   Pool results of tests from two independent
                        datasets
pp                      Unpack posterior probabilities from one-byte
                        codes
qq.chisq                Quantile-quantile plot for chi-squared tests
read.beagle             Read genotypes imputed by the BEAGLE program
read.impute             Read genotypes imputed by the IMPUTE2 program
read.long               Read SNP genotype data in long format
read.mach               Read genotypes imputed by the MACH program
read.pedfile            Read a pedfile as '"SnpMatrix"' object
read.plink              Read a PLINK binary data file as a SnpMatrix
read.snps.long          Read SNP data in long format (deprecated)
row.summary             Summarize rows or columns of a snp matrix
sample.info             Sample datasets to illustrate data input
single.snp.tests        1-df and 2-df tests for genetic associations
                        with SNPs (or imputed SNPs)
sm.compare              Compare two SnpMatrix objects
snp.cbind               Bind together two or more SnpMatrix objects
snp.cor                 Correlations with columns of a SnpMatrix
snp.imputation          Calculate imputation rules
snp.lhs.estimates       Logistic regression with SNP genotypes as
                        dependent variable
snp.lhs.tests           Score tests with SNP genotypes as dependent
                        variable
snp.pre.multiply        Pre- or post-multiply a SnpMatrix object by a
                        general matrix
snp.rhs.estimates       Fit GLMs with SNP genotypes as independent
                        variable(s)
snp.rhs.tests           Score tests with SNP genotypes as independent
                        variable
snpStats-package        SnpMatrix and XSnpMatrix classes and methods
switch.alleles          Switch alleles in columns of a SnpMatrix or in
                        test results
tdt.snp                 1-df and 2-df tests for genetic associations
                        with SNPs (or imputed SNPs) in family data
test.allele.switch      Test for switch of alleles between two
                        collections
testdata                Test data for the snpStats package
write.SnpMatrix         Write a SnpMatrix object as a text file
write.plink             Write files for analysis in the PLINK toolset
xxt                     X.X-transpose for a standardized SnpMatrix
}

Further information is available in the following vignettes:
\tabular{ll}{
\code{data-input-vignette} \tab Data input (source)\cr
\code{differences} \tab snpMatrix-differences (source)\cr
\code{imputation-vignette} \tab Imputation and meta-analysis (source)\cr
\code{ld-vignette} \tab LD statistics (source)\cr
\code{pca-vignette} \tab Principal components analysis (source)\cr
\code{snpStats-vignette} \tab snpStats introduction (source)\cr
\code{tdt-vignette} \tab TDT tests (source)\cr
}
}
\author{
David Clayton <dc208@cam.ac.uk>

Maintainer: David Clayton <dc208@cam.ac.uk>
}
\keyword{ package }
