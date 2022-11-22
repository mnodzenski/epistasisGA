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
#' for each SNP in the example datasets \code{case}, \code{mom}, and \code{dad}.
#' The SNPs are in the same order as they appear in the example dataset columns.
#'
#' @format A data frame with 100 rows and 3 variables
#' @usage data(snp.annotations)
"snp.annotations"

#' Genotypes for the mothers of case-parent triads with a simulated
#' maternal-fetal interaction.
#'
#' A simulated dataset containing the counts of the alternate allele
#' for 24 SNPs for the mothers in 1000 simulated case-parent
#' triads. Columns represent SNPs, rows are individuals. The SNP
#' in column 6 interacts with the SNPs in columns 12 and 18
#' of dataset \code{case.mci} to increase risk of disease in the child,
#' where at least one copy of the alternate allele (genotype 1 or 2)
#' is required at each implicated locus.
#' .
#'
#' @format A matrix with 1000 rows and 24 variables
#' @usage data(mom.mci)
"mom.mci"

#' Genotypes for the fathers of case-parent triads with a simulated
#' maternal-fetal interaction.
#'
#' A simulated dataset containing the counts of the alternate allele
#' for 24 SNPs for the fathers in 1000 simulated case-parent
#' triads. Columns represent SNPs, rows are individuals. The SNP
#' in column 6 of the corresponding maternal dataset \code{mom.mci}
#' interacts with the SNPs in columns 12 and 18
#' of the corresponding child dataset \code{case.mci} to increase risk of
#' disease in the child, where at least one copy of the alternate allele
#' (genotype 1 or 2) is required at each implicated locus.
#' .
#'
#' @format A matrix with 1000 rows and 24 variables
#' @usage data(dad.mci)
"dad.mci"

#' Genotypes for the affected cases of case-parent triads with a simulated
#' maternal-fetal interaction.
#'
#' A simulated dataset containing the counts of the alternate allele
#' for 24 SNPs for the cases in 1000 simulated case-parent
#' triads. Columns represent SNPs, rows are individuals. The SNP
#' in column 6 of the corresponding maternal dataset \code{mom.mci}
#' interacts with the SNPs in columns 12 and 18
#' of \code{case.mci} to increase risk of
#' disease in the child, where at least one copy of the alternate allele
#' (genotype 1 or 2) is required at each implicated locus.
#' .
#'
#' @format A matrix with 1000 rows and 24 variables
#' @usage data(case.mci)
"case.mci"

#' RSID, REF, and ALT annotations for example dataset SNPs
#'
#' A data.frame containing the RSID, REF allele and ALT allele
#' for each SNP in the example datasets \code{case.mci}, \code{mom.mci},
#' \code{dad.mci}, \code{case.gxe}, \code{mom.gxe}, and \code{dad.gxe}.
#' The SNPs are in the same order as they appear in the example dataset columns.
#'
#' @format A matrix with 24 rows and 3 variables
#' @usage data(snp.annotations.mci)
"snp.annotations.mci"

#' Genotypes for the fathers of case-parent triads with a simulated gene
#' environment interaction.
#'
#' A simulated dataset containing the counts of the alternate allele
#' for 24 SNPs for the fathers in 1000 simulated case-parent
#' triads. Columns represent SNPs, rows are individuals. SNPs
#' in columns 6, 12, and 18 represent a simulated risk pathway,
#' where, in the child, at least one copy of the alternate allele
#' for each path SNP in addition to exposure 1 confers increased
#' disease risk.
#' .
#'
#' @format A data frame with 1000 rows and 24 variables
#' @usage data(dad.gxe)
"dad.gxe"

#' Genotypes for the cases of case-parent triads with a simulated gene
#' environment interaction.
#'
#' A simulated dataset containing the counts of the alternate allele
#' for 24 SNPs for the cases in 1000 simulated case-parent
#' triads. Columns represent SNPs, rows are individuals. SNPs
#' in columns 6, 12, and 18 represent a simulated risk pathway,
#' where, in the child, at least one copy of the alternate allele
#' for each path SNP in addition to exposure 1 confers increased
#' disease risk.
#' .
#'
#' @format A data frame with 1000 rows and 24 variables
#' @usage data(case.gxe)
"case.gxe"

#' Genotypes for the mothers of case-parent triads with a simulated gene
#' environment interaction.
#'
#' A simulated dataset containing the counts of the alternate allele
#' for 24 SNPs for the mothers in 1000 simulated case-parent
#' triads. Columns represent SNPs, rows are individuals. SNPs
#' in columns 6, 12, and 18 represent a simulated risk pathway,
#' where, in the child, at least one copy of the alternate allele
#' for each path SNP in addition to exposure 1 confers increased
#' disease risk.
#' .
#'
#' @format A data frame with 1000 rows and 24 variables
#' @usage data(mom.gxe)
"mom.gxe"

#' Genotypes for the fathers of case-parent triads with a simulated gene
#' environment interaction.
#'
#' A simulated dataset containing the counts of the alternate allele
#' for 24 SNPs for the fathers in 1000 simulated case-parent
#' triads. Columns represent SNPs, rows are individuals. SNPs
#' in columns 6, 12, and 18 represent a simulated risk pathway,
#' where, in the child, at least one copy of the alternate allele
#' for each path SNP in addition to exposure 1 confers increased
#' disease risk.
#' .
#'
#' @format A data frame with 1000 rows and 24 variables
#' @usage data(dad.gxe)
"dad.gxe"

#' Genotypes for the cases of case-parent triads with a simulated gene
#' environment interaction.
#'
#' A simulated dataset containing the counts of the alternate allele
#' for 24 SNPs for the cases in 1000 simulated case-parent
#' triads. Columns represent SNPs, rows are individuals. SNPs
#' in columns 6, 12, and 18 represent a simulated risk pathway,
#' where, in the child, at least one copy of the alternate allele
#' for each path SNP in addition to exposure 1 confers increased
#' disease risk.
#' .
#'
#' @format A data frame with 1000 rows and 24 variables
#' @usage data(case.gxe)
"case.gxe"

#' Exposures for the cases of case-parent triads with a simulated gene
#' environment interaction.
#'
#' A data.frame containing simulated exposure status for each case
#' of the case-parent triads data. Rows correspond to different
#' families. The single column represents a binary exposure, where
#' in combination with the relevant risk-associated alleles (columns
#' 6, 12, and 18 in data set case.gxe), is associated with increased risk.
#' .
#'
#' @format A data frame with 1000 rows and 1 variables
#' @usage data(exposure)
"exposure"

