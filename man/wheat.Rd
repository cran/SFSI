\name{wheat}
\alias{wheatHTP}
\alias{Y}
\alias{X}
\alias{WL}
\alias{VI}
\alias{wheatHTP.E1}
\alias{wheatHTP.E2}
\alias{wheatHTP.E3}
\alias{wheatHTP.E4}
\title{Wheat dataset}
\description{

The data set consists of 1,092 inbred wheat lines grouped into 39 trials and grown during the 2013-2014 season at the Norman Borlaug experimental research station in Ciudad Obregon, Sonora, Mexico. Each trial consisted of 28 breeding lines that were arranged in an alpha-lattice design with three replicates and six sub-blocks. The trials were grown in four different environments:
\itemize{
\item E1: Flat-Drought (sowing in flat with irrigation of 180 mm through drip system)
\item E2: Bed-2IR (sowing in bed with 2 irrigations approximately 250 mm)
\item E3: Bed-5IR (bed sowing with 5 normal irrigations)
\item E4: Bed-EHeat (bed sowing 30 days before optimal planting date with 5 normal irrigations approximately 500 mm)
}

\ifelse{html}{\out{<b>Phenotypic data.</b>}}{\eqn{\textbf{Phenotypic data.}}{Phenotypic data.}}
Measurements of grain yield (YLD) were reported as the total plot yield after maturity. Records for YLD are reported as adjusted means from which trial, replicate and sub-block effects were removed. Measurements for days to heading (DTH), days to maturity (DTM), and plant height (PH) were recorded only in the first replicate at each trial and thus no phenotype adjustment was made.

\ifelse{html}{\out{<b>Reflectance data.</b>}}{\eqn{\textbf{Reflectance data.}}{Reflectance data.}}
Reflectance data was collected from the fields using both infrared and hyper-spectral cameras mounted on an aircraft on 9 different dates (time-points) between January 10 and March 27th, 2014. During each flight, data from 250 wavelengths ranging from 392 to 850 nm were collected for each pixel in the pictures. The average reflectance of all the pixels for each wavelength was calculated from each of the geo-referenced trial plots and reported as each line reflectance. Data for reflectance and Green NDVI and Red NDVI are reported as adjusted phenotypes from which trial, replicate and sub-block effects were removed.

\ifelse{html}{\out{<b>Marker data.</b>}}{\eqn{\textbf{Marker data.}}{Marker data.}}
Lines were sequenced using GBS technology at 192-plexing on Illumina HiSeq2000 or HiSeq2500 with 1 x 100 bp reads. SNPs were called across all lines anchored to the genome assembly of Chinese Spring (International Wheat Genome Sequencing Consortium 2014). Next, SNP were extracted and filtered so that lines >50\% missing data were removed. Markers were recoded as –1, 0, and 1, corresponding to homozygous for the minor allele, heterozygous, and homozygous for the major allele, respectively. Next, markers with a minor allele frequency >0.05 and >15\% of missing data were removed. Remaining SNPs with missing values were imputed using the mean of the observed marker genotypes at a given locus. Only data for 967 lines evaluated in 35 trials with marker information for 3,274 SNPs remained after QC.

Marker information, phenotypic and reflectance for all four environments (wheatHTP.E1,...,wheatHTP.E4 data sets) are provided through the GitHub (development version) of the SFSI R-package available at https://github.com/MarcooLopez/SFSI.

The CRAN version of the package includes data from the Bed-5IR environment only and reflectance data for the latest time-point only (wheatHTP data set). This data set contains the phenotypic and reflectance adjusted means for 967 lines evaluated in 35 trials for which marker information on 3,274 SNPs is available.
}

\usage{
  data(wheatHTP)
  # data(wheatHTP.E1) # GitHub version
  # data(wheatHTP.E2) # GitHub version
  # data(wheatHTP.E3) # GitHub version
  # data(wheatHTP.E4) # GitHub version
}

\format{
\itemize{
\item \code{Y}: (matrix) phenotypic data for YLD, DTH, DTM, and PH; and the trial in which each genotype was tested. 
\item \code{X}: (matrix) marker data with SNPs in columns (only in wheatHTP data set).
\item \code{WL[[i]]}: (list) reflectance data for each time-point (\code{i=1,...,9}).
\item \code{WL}: (matrix) reflectance data for time-point 9 (only in wheatHTP data set).
\item \code{VI[[i]]}: (list) green and red NDVI for each time-point (\code{i=1,...,9}). 
\item \code{VI}: (matrix) green and red NDVI for each time-point 9 (only in wheatHTP data set). 
}
}

\source{
  International Maize and Wheat Improvement Center (CIMMYT), Mexico.
}

\keyword{datasets}