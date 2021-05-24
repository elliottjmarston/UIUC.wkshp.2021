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

# response to selection = genetic gain
# R = gt - gt-1
# R = delta g

# Prediction of selection reponse hinges on 2 basic principles
# Principle 1
#   - breeding value can be predicted from phenotypic value
#     using the linear regression equation
#   - linear regression: breeding value vs phenotypic value
#   - same principle applies to groups of individuals as individuals
# Principle 2
#   - average true or predicted breeding value of two parents
#     = average breeding value of their progeny

# HOW CAN WE PREDICT SELECTION RESPONSE USING THESE PRINCIPLES?

# breeding values of population 0, truncate selection
#   - select everything w/ phenotype above value in normal distribution

# regression coefficient for the regression of breeding value
#   on phenotype or selection index values

# The Breeder's equation
#   R = rxa * theta-a * k
# response per breeding cycle = selection accurary *
#   additive genetic standard deviation * selection intensity

# Selection intensity = k = z / p
# z = height of the distribution at the truncation point
# p = proportion severed (truncated portion)

# Expected gain per unit time
# R / L
# L = time (usually years) required to complete one breeding cycle
#   - expected slope of tehr egression of breeding value on time

# Assumptions
#   - discrete generations
#   - constant variance (infinitesimal model), accuracy, intensity
#   - truncation selection

# Key points
#   - breeder's equation can be derived from two basic principles
#   - gain from selection per unit time depends on selection
#     accuracy, additive genetic standard deviation, selection
#     intensity, and length of breeding cycle
#   - breeders' equation makes several assumptions that are
#     often violated in reality

# Exercise
# population mean = 10, standard dev = 9
# mean of selected = 17, truncation point = 15
# percent selected = 20%

# R = pop mean - mean of selected
# R = 17 - 10 = 7
# R = 7?

## LECTURE 4 ----
# GENOMIC SELECTION

# OBJECTIVES
#   - learn what is genomic selection and its application to breeding
#   - learn basic guidelines for implementing genomic selection
#   - understand GBLUP and RR-BLUP

# Genomic Selection
#   - GS is using genome-wide marker data to help predict
#     breeding value for the purpose of selection
#   - GS is useful for traits conferred by at least 3 loci
#   - GS is a disruptive breeding technology
#       - it will change how the breeding program operates if
#         done correctly

# Why Gs?
#   - increase selection accuracy
#   - reduce the breeding cycle time bc parents can be selected
#     sooner

# GS de-couples population improvement and variety development
#   - generating "evaluation units", lines, hybrids, or clones
#     is required for releasing varieties and also to generate
#     data used in GS model

# GS selection model training set
#   - the model training set refers to the set of genotypes
#     used in the GS model that have phenotypes for the primary
#     traits of interest
#   - model training set also referred to as training pop,
#     calibration set, reference set
#   - the size and composition of the training set is the most
#     important factor affecting GS accuracy

# Training sets and Selection Candidates
#   - model training set can include your selection candidates
#   - model traiing set may not include your selection candidates

# What makes a good training set?
#   - data should be relevant: from target pop of environments
#     on primary traits of interest at minimum
#   - data should be meaningful, collected consistently w/
#     best practices
#   - germplasm used for model training set should be closely
#     related to your selection candidates

# Updating the training set
#   - germplasm used for model training set should be closely
#     related to selection candidates
#   - with each breeding cycle, traning set will become less
#     predictive as genetic realationship to selection candidates
#     decreases, therefore the training set must be continually
#     updated

# Training Set Size
#   - incresaing size will increase selection accuracy
#   - helps maintain accuracy over cycles of selection
#   - number of individuals needed to achieve high GS accuracy
#     depends on effective pop size and heritability of trait
#     of interest

# Expected accuracy based on training set size 
#   and population history
#   - Me is the # of independent chrosome segments
#   - Np is the training set size

# Number of genome-wide markers required
#   - depends on rate of LD decay in the population
#   - need at least one marker per segregating 
#     segment of the genome
#   - 10 NeL markers required to achieve high accuracies,
#     where L is the genome size in Morgans
#   - NeL markers could be used to achieve moderate to high accuracies

# Choice of Prediction Model
#   - genomic selection w/ Genome-wide markers = Quant
#       - RR/G-BLUP
#   - vs. Marker Assisted Selection w/ QTL linked markers = Mendelian
#       - Multiple linear regression
#   - Bayesian Lasso Bayes mid point between the two

# GBLUP or RR-BLUP work well in most cases
#   - GBLUP is equivalent to RR-BLUP and assumes 
#     infinitesimal model of inheritance
#   - across range of datasets GBLUP has been shown to perform well

# GBLUP
#   - Genomic Relationship Matrix (G)
#       - will not be analogous to a pedigree relationship matrix
#       - inbreeding coefficients cannot be obtained from
#         diagonal of B and estimates of genetic variance
#         will not tell us the additive genetic variance in base pop

# Model Fitting in Two Steps
#   - reduces computational complexity
#   Step 1: generate means by environment or across all
#       environments using a mixed model assuming genotypes
#       are not related
#   Step 2: Fit the gS model using the means generated in
#       STep 1 and be sure to account for heterogenous error
#       variances

# Implementing GS
#   - Fitting the model is the easiest part
#   1. Use a database
#   2. Start by genotyping everything that is being phenotyped
#       and selected based on a GS model including all available
#       data
#   3. Avoid selection w/o recording phenotypic data
#   4. Begin genotyping and selecting individuals prior to
#       phenotyping when enough genotypic and phenotypic data
#       has been accumulated

# True or False?
# The goal of GS is to predict phenotypic value?
# FALSE, predict BREEDING value


## LECTURE 5 ----
# APPLICATION OF GS AND FACTORS AFFECTING ITS SUCCESS

#OBJECTIVES
#   - learn different ways GS can be used to improve rates
#       of genetic gain
#   - Understand what are the critical factors affecting
#     success of a GS strategy

# Application of GS
# 2 types
#   - rapid cycle recurrent selection
#     apply GS among non-inbred selection candidates
#     challenging, greater potential to accelerate genetic gain
#     requires fast, cheap genotyping and a very large training set
#   - GS among lines
#     apply GS after creating a line
#     less potential to accelerate genetic gain, easier to
#     implement

# Can do both rapid cycle GS + GS among lines
#   - L = 0.5 years min
#   - all during line advancement you're using rapid GS techniques
#   - potential benefits:
#       - Faster levels of genetic gain w/i breeding programs

# GS among lines: strategies
#   - genotype more lines than are phenotyped
#       - phenotype using a conventional testing strategy
#       - phenotyping using a sparse testing strategy
#   - genotype and phenotype all lines
#       - phenotype using a conventional testing strategy
#       - phenotype using a sparse testing strategy

# GS among lines: optimal testing strategy depends on costs
#   - low genotyping cost relative to phenotyping
#     better to genotype more lines, use less reps on each,
#     and possibly not phenotype some lines
#   - similar genotyping costs relative to phenotyping
#     better to have some reps, fewer lines, and phenotype all
#     lines that have been genotyped

# Factors affecting success of a GS program
#   - strategy
#   - operations (is it well equipped to carry out strategy)

# What makes a Good GS Strategy?
#   - GS applied to all traits of interest, especially low
#     heritability traits
#   - involves phenotyping a large # of individuals from the
#     breeding program every year for model training and updating
#   - uses genomic prediction to reduce the breeding cycle time
#       - parent selection done based on GEBVs
#   - Allocates testing resources appropriately 
#       - possibly using sparese testing
#   - involves cheap genotyping w/ an adequate number of markers

# What are critical operational factors?
#   - phenotyping: as accurate as possible given allocated resources
#   - data management: relational database used for all phenotypic
#       and marker data
#   - germplasm management: the right seed gets in the right envelope
#   - fast turn around time of inexpensive marker data

## Example ----

load("Breeding Population Data.RData")















