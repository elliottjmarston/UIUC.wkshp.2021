#### UIUC WORKSHOP ON APPLIED QUANTITATIVE GENETICS FOR PLANT BREEDERS ####
#### SET UP AND NOTES ####

### DEPENDENCIES ----

# set wd

setwd("C:/Users/Elliott/R/UIUC.wkshp/UIUC.wkshp.2021")

# libraries

library(lme4)
library(ggplot2)
library(rrBLUP)
library(qqman)
library(arm)
library(MASS)
library(BiocManager)
library(multtest) # not available for current version of Rstudio
library(gplots)

### DAY 1 NOTES ----
# Dr. Jessica Rutkoski

## LECTURE 1 ----
# MODELS OF INHERITANCE AND EFFECTS OF SINGLE LOCI

# OBJECTIVE: introducce/refresh basic population genetics concepts

# terminology: locus, gene, allele, genotype, segregation, population

# Mendelian model of Inheritance
#   - states that traits can be explained by single genes
#   - Ex. seed shattering, pea color
#   - mendel's law of inheritance: law of dominance, law of
#     segregation, law of independent assortment
#   - did not seem to explain continuous variation

# Infinitesimal model of inheritance
#   - traits can be explained by an impossibly large number
#      of unlinked additive alleles
#   - Ex. human height, skin color
#   - Explains that continuous trait variation is due to the 
#      segregation of alleles at a very large number of loci
#   - reconciles mendelian inheritance with continuous trait
#     variation

# Infinitesimal model assumptions
#   - an infinite number of unlinked loci affect the trait,
#     and each locus has a very small effect
#   - each locus has the same effects and frequencies
#   - selection will not affect allele frequencies, nor
#     genetic variance

# Most traits are infinitesimal
#   - most of the time, traits have been subject to selection
#     follow an infinitesimal model of inheritance
#   - why?
#      - large effect favorable alleles rapidly become fixed
#         and only small effect alleles remain segregating
#         in the population

# Allele effects
#   - additive: allele effect does not depend on other allele states
#       - addition of alleles has fixed additive effect on
#         phenotype
#   - non-additive: allele effect depends on other allele states
#       - dominance: some alleles may be masked
#           - phenotype dependent on presence of dominant allele
#       - epistatic: depends on genotype at other locus

# Breeding Value
#   - additive genotypic effect at a single locus, or summed
#     across multiple loci
#   - the genotypic effect that is passed on to progeny upon mating
#   - key for selecting breeding parents because it tells us
#     the expected genetic value of the progeny

# Breeding Value of a Genotype
#   - additive allele substitution effect x genotype =
#     breeding value of a genotype at a single locus
#   EXample:
#   additive allele sub effect = 0.001 = alpha
#   coded genotype = 1
#   breeding value = 0.001 * 1 = 0.001
#   Example: multi-locus genotype
#   Locus X breeding value = 0.001 x 1
#   Locus Y breeding value = 0.04 * -1
#   Total breeding value = 0.001 - 0.04 = -0.39

##  the additive allele effect = 
##  measured effect of allele on observed phenotype

# Additivity explains most effects
#   - dominance and epistasis are expressed as 
#     deviations from additivity
#   - even at a completely dominant locus, an additive model
#     can explain a large part of the variation

## LECTURE 2 ----
# PHENOTYPES, GENOTYPES, AND BREEDING VALUES #

# OBJECTIVES:
#   - make connection between genotypes and phenotypes
#   - understand phenotypes and breeding values at 
#     the population level

# Phenotypes
#   - phenotype = environment + genetics + residual
#   y-i = mew-i + g-i + sigma-i
#       environmental effect "permanent evironment", 
#       genetic value, residual "specific environment"

# Genetic Value
#   - can be broken down into different components
#   - g-i = a-i (additive genetic = breeding value) +
#     d-i (dominance) + e-i (epistasis)

# Breeding value
#   - of individuals
#       - the portion of genetic value that is transmitted
#         from parent to progeny
#       - value as a parent
#       - sum total of all additive effects of individual's alleles

# Breeding value prediction based on phenotype
#   - a-i = a + bxa(x-i - x)
#   average breeding value of population + 
#   regression coefficient (phenotypic value of individual 
#   - population mean phenotype)

# Transmitting ability
#   - average effect of a random sample of half of an individual's alleles
#   - 1/2 of an individual's total breeding value

# Breeding value of parents and progeny
#   - each parent contributes half of its alleles to the progeny
#   - average breeding balue of progeny is the average breeding
#     value of the two parents
#   E(a) = 1/2 ap1 + 1/2 ap2

# Variation due to mendelian sampling
#   - breeding valu eof progeny vary due to random sampling
#     of alleles
#   ai = 1/2 ap1 + 1/2 ap2 + mi (mendelian sampling)

# Mendelian sampling term
#   - how much an individual's breeding value deviates from
#     the average breeding value of its parents
#   - mendelian sampling allows us to ID individuals that are
#     superior to their parents

# The normal distribution
#   - mean: measure of center
#   - variance: measure of spread
#   - standard deviation: square root of variance, same units as trait

# Variance: mean of squared deviations
#   - Total phenotypic variance = genetic variance + error variance
#   - additive genetic variance: variance of the true breeding values

# Covariance: measures how two variables change together
# Cov(x,y) = theta^2xy = sigma(xi-x)(yi-y) / n

# Correlation
#   - normalized version of covariance, ranges -1 to 1
#   - indicates the strength of a linear relationship between
#     two variables
#   rxy = Cov(x,y) / theta-x * theta-y

# Correlation and regression coefficient
#   - from standard regression theory, the regression coefficient
#     for the regression of y on x is 
#   byx = theta-xy/theta^2x = rxy * theta-y / theta-x

# Heritability
#   - regression of the breeding value on the phenotypic value
#   - an individual's breeding value is the mean of its progeny's
#     phenotypic value
#   - regression coefficient of the breeding value on the phenotypic value
#   - square root of heritability
#       - accuracy of selection based on phenotypic values
#       - correlation coefficient between breeding value and
#         phenotypic value

# Heritability and the base population
#   - it is the proportion of the total phenotypic variance
#     that is due to additive genetic factors in the base population
#   - phenotypes are on individual non-inbred plants

# Alternative definitions of heritability
#   - broad-sense heritability: non-additive variance is included
#   - line-mean heritability: broad sense heritability
#       - square of the selection accuracy based on a mean
#         of multiple phenotypes. Referred to as "reliability"
#         in animal breeding.

## LECTURE 3 ----
# BASIC PRINCIPLES OF RESPONSE TO SELECTION
















