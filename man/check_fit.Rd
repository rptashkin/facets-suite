% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/check-fit.R
\name{check_fit}
\alias{check_fit}
\title{Sample QC}
\usage{
check_fit(
  facets_output,
  genome = c("hg19", "hg18", "hg38"),
  algorithm = c("em", "cncf"),
  maf = NULL
)
}
\arguments{
\item{facets_output}{Output from \code{run_facets}.}

\item{genome}{Genome build.}

\item{algorithm}{Choose assessing the fit from the \code{em} or \code{cncf} algorithm.}

\item{maf}{Optional: mutation calls for assessed samples, should only include mutations for given sample.}
}
\value{
A list object with the following items:
\itemize{
      \item{\code{dipLogR_flag}:} {Boolean indicating extreme dipLogR value.}
      \item{\code{n_alternative_dipLogR}:} {Number of alternative dipLogR values.}
      \item{\code{n_dip_bal_segs}, \code{frac_dip_bal_segs}:} {Number of balanced segments at dipLogR and the fraction of genome they represent.}
      \item{\code{n_dip_imbal_segs}, \code{frac_dip_imbal_segs}:} {Number of imbalanced segments at dipLogR and the fraction of genome they represent.}
      \item{\code{n_amp}:} {Number of segments at total copy number >= 10.}
      \item{\code{n_homdels}:} {Number of homozygously deleted segments (total copy number = 0).}
      \item{\code{n_homdels_clonal}, \code{frac_homdels_clonal}:} {Number of clonal homdel segments and the fraction of the genome they represent.}
      \item{\code{n_cn_states}:} {Number of unique copy-number states (i.e. combinations of major and minor copy number).}
      \item{\code{n_segs}:} {Number of segments.}
      \item{\code{n_cnlr_clusters}:} {Number of copy-number log-ratio clusters}
      \item{\code{n_lcn_na}:} {Number of segments where no minor copy number was inferred (lcn is NA).}
      \item{\code{n_loh}, \code{n_loh}:} {Number of segments where the minor copy number is 0 and the fraction of the genome they represent.}
      \item{\code{n_snps}:} {Number of SNPs used for segmentation.}
      \item{\code{n_het_snps}, \code{frac_het_snps}:} {Number of heterozyous SNPs used for segmentation and their fraction of the total.}
      \item{\code{n_het_snps_hom_in_tumor_1pct}, \code{frac_het_snps_hom_in_tumor_1pct}:} {Number of heterozyous SNPs where the tumor allele frequency is <0.01/>0.99 their fraction of the total.}
      \item{\code{n_het_snps_hom_in_tumor_5pct}, \code{frac_het_snps_hom_in_tumor_5pct}:} {Number of heterozyous SNPs where the tumor allele frequency is <0.05/>0.95 their fraction of the total.}
      \item{\code{mean_cnlr_residual}, \code{sd_cnlr_residual}:} {Mean and standard deviation of SNPs' log-ratio from their segments copy-number log-ratio.}
      \item{\code{n_segs_discordant_tcn}, \code{frac_segs_discordant_tcn}:} {Number of segments where the naïve and EM algorithm estimates of the total copy number are discordant and the fraction of the genome they represent.}
      \item{\code{n_segs_discordant_lcn}, \code{frac_segs_discordant_lcn}:} {Number of segments where the naïve and EM algorithm estimates of the minor copy number are discordant and the fraction of the genome they represent.}
      \item{\code{n_segs_discordant_both}, \code{frac_segs_discordant_both}:} {Number of segments where the naïve and EM algorithm estimates of the both copy numbers are discordant and the fraction of the genome they represent.}
      \item{\code{n_segs_icn_cnlor_discordant}, \code{frac_icn_cnlor_discordant}:} {Number of clonal segments where the log-ratio shows balance but the copy-number solution does not, and the reverse, and the fraction of the genome they represent.}
      \item{\code{dip_median_vaf}:} {If MAF input: median tumor VAF of somatic mutations on clonal segments with total copy number 2 and allelic balance.}
      \item{\code{n_homdel_muts}:} {If MAF input: number of somatic mutations in homozygously deleted segments.}
      \item{\code{median_vaf_homdel_muts}:} {If MAF input: Median tumor VAF of somatic  mutations homozygously deleted segments.}
   }
}
\description{
Generate QC metrics for sample. Can use mutation calls in MAF file.
}
