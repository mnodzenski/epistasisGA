#' Genotypes for the mothers of case-parent triads.
#'
#' A simulated dataset containing the counts of the alternate allele
#' for 100 SNPs for the mothers in 1000 simulated case-parent
#' triads. Columns represent SNPs, rows are individuals. SNPs
#' in columns 51, 52, 76, and 77 represent a true risk pathway.
#'
#' @format A data frame with 1000 rows and 100 variables
#' @usage data(mom)
"mom"

#' Genotypes for the fathers of case-parent triads.
#'
#' A simulated dataset containing the counts of the alternate allele
#' for 100 SNPs for the fathers in 1000 simulated case-parent
#' triads. Columns represent SNPs, rows are individuals. SNPs
#' in columns 51, 52, 76, and 77 represent a true risk pathway.
#'
#' @format A data frame with 1000 rows and 100 variables
#' @usage data(dad)
"dad"

#' Genotypes for the affected children of case-parent triads.
#'
#' A simulated dataset containing the counts of the alternate allele
#' for 100 SNPs for the affected child in 1000 simulated case-parent
#' triads. Columns represent SNPs, rows are individuals. SNPs
#' in columns 51, 52, 76, and 77 represent a true risk pathway.
#'
#' @format A data frame with 1000 rows and 100 variables
#' @usage data(case)
"case"

#' RSID, REF, and ALT annotations for example dataset SNPs
#'
#' A data.frame containing the RSID, REF allele and ALT allele
#' for each SNP in the example datasets. The SNPs are in the same
#' order as they appear in the example dataset columns.
#'
#' @format A data frame with 100 rows and 3 variables
#' @usage data(snp.annotations)
"snp.annotations"

#' Genotypes for the mothers of case-parent triads with a simulated gene
#' environment interaction.
#'
#' A simulated dataset containing the counts of the alternate allele
#' for 100 SNPs for the mothers in 1000 simulated case-parent
#' triads. Columns represent SNPs, rows are individuals. SNPs
#' in columns 25, 50, and 75 represent a simulated risk pathway,
#' where at least one copy of the alternate allele for each path SNP
#' in addition to exposure 2 confers increased disease risk.
#' .
#'
#' @format A data frame with 1000 rows and 100 variables
#' @usage data(mom.gxe)
"mom.gxe"

#' Genotypes for the fathers of case-parent triads with a simulated gene
#' environment interaction.
#'
#' A simulated dataset containing the counts of the alternate allele
#' for 100 SNPs for the fathers in 1000 simulated case-parent
#' triads. Columns represent SNPs, rows are individuals. SNPs
#' in columns 25, 50, and 75 represent a simulated risk pathway,
#' where at least one copy of the alternate allele for each path SNP
#' in addition to exposure 2 confers increased disease risk.
#' .
#'
#' @format A data frame with 1000 rows and 100 variables
#' @usage data(dad.gxe)
"dad.gxe"

#' Genotypes for the cases of case-parent triads with a simulated gene
#' environment interaction.
#'
#' A simulated dataset containing the counts of the alternate allele
#' for 100 SNPs for the cases in 1000 simulated case-parent
#' triads. Columns represent SNPs, rows are individuals. SNPs
#' in columns 25, 50, and 75 represent a simulated risk pathway,
#' where at least one copy of the alternate allele for each path SNP
#' in addition to exposure 2 confers increased disease risk.
#' .
#'
#' @format A data frame with 1000 rows and 100 variables
#' @usage data(case.gxe)
"case.gxe"

#' Exposures for the cases of case-parent triads with a simulated gene
#' environment interaction.
#'
#' A simulated vector containing the exposure status for each case
#' of the case-parent triads data.
#' .
#'
#' @format A data frame with 1000 rows and 100 variables
#' @usage data(exposure)
"exposure"

