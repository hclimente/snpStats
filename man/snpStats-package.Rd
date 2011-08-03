\name{snpStats-package}
\alias{snpStats-package}
\alias{snpStats}
\docType{package}
\title{
SnpMatrix and XSnpMatrix classes and methods
}
\description{
Classes and statistical methods for large SNP association studies, extending the snpMatrix package
}
\details{
This package extends the \code{snpMatrix} package allowing for uncertainty
in genotype assignments, for example, as  consequence of
imputation. It provides most of the functionality of \code{snpMatrix}, although
it currently lacks a few of the functions. It is intended that the most 
important of these will be provided in the near future. It also
includes some novel commands. Note that the \code{snpMatrix} classes
have been renamed and, in two cases, changed slightly. For more
details, see the "differences" vignette.

Further development of the snpMatrix project by this developer will
take place entirely within \code{snpStats}.

\tabular{ll}{
Package: \tab snpStats\cr
Version: \tab 1.3.4\cr
Date: \tab 2011-08-03\cr
Depends: \tab R(>= 2.3.0), survival, methods, Matrix\cr
Imports: \tab graphics, grDevices, methods, stats, survival, utils, Matrix\cr
Suggests: \tab hexbin\cr
License: \tab GPL-3\cr
URL: \tab http://www-gene.cimr.cam.ac.uk/clayton\cr
Collate: \tab ss.R contingency.table.R convert.R compare.R glm-test.R imputation.R indata.R long.R misc.R ld.R mvtests.R pedfile.R outdata.R plink.R qc.R qq-chisq.R single.R tdt-single.R structure.R xstuff.R zzz.R\cr
LazyLoad: \tab yes\cr
biocViews: \tab Microarray, SNP, GeneticVariability\cr
Packaged: \tab 2011-08-03 10:21:40 UTC; david\cr
Built: \tab R 2.13.1; x86_64-unknown-linux-gnu; 2011-08-03 08:39:51 UTC; unix\cr
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
glm.test.control        Set up control object for GLM tests
ibsCount                Count alleles identical by state
ibsDist                 Distance matrix based on identity by state
                        (IBS)
imputation.maf          Extract statistics from imputation rules
impute.snps             Impute snps
ld                      Pairwise linkage disequilibrium measures
ld.example              Datasets to illustrate calculation of linkage
                        disequilibrium statistics
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
\code{data-input-vignette} \tab Data input (source, pdf)\cr
\code{differences} \tab snpMatrix-differences (source, pdf)\cr
\code{imputation-vignette} \tab Imputation and meta-analysis (source, pdf)\cr
\code{ld-vignette} \tab LD statistics (source, pdf)\cr
\code{pca-vignette} \tab Principal components analysis (source, pdf)\cr
\code{snpStats-vignette} \tab snpStats introduction (source, pdf)\cr
\code{tdt-vignette} \tab TDT tests (source, pdf)\cr
}
}
\author{
David Clayton \email{david.clayton@cimr.cam.ac.uk}

Maintainer: David Clayton \email{david.clayton@cimr.cam.ac.uk}
}
\keyword{ package }
