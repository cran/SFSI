### Version 1.1.0 (Mar-10-2022)

**New features**

  - Some problems were fixed in the documentation structure as required by CRAN.
  - Functions 'lars2', 'SSI_CV', 'plotNet', 'plotPath' changed their names to 'LARS', 'SSI.CV', 'net.plot', and 'path.plot', respectively.
  - Some arguments' functions changed their names to a more informative name (e.g., 'minLambda' => 'lambda.min')
  - More functionalities added to 'net.plot' function


### Version 1.0.1 (Jan-26-2022)

**New features**

  - Functions 'SSI' and 'SSI_CV' allow providing either 'theta' (residual/genetic variances ratio) or the 'h2' (heritability)

**Bug fixes**

  - C-based routine associated to the 'readBinary' function now uses the 'Rf_allocMatrix' method to handle matrices whose length (number of rows x number of columns) exceed 2^31-1 = 2147483647


### Version 1.0.0 (Sep-30-2021)

**New features**

  - Function 'solveEN' allows early stop when a user-provided number of non-zero predictors (at a given value of lambda) is reached (argument 'maxDF')
  - Functions 'solveEN' and 'lars2' return object 'beta' as matrix with predictors in rows (rather than in columns)
  - Function 'cov2cor2' allows multiplying the resulting correlation matrix times a constant 'a' (default is 'a=1')
  - Provided 'wheatHTP' dataset includes now an array of 4-folds partitions ('CV' column in object 'Y') and calculations of genetic and residual covariances between YLD and each of the wavelengths ('genCOV_xy' and 'resCOV_xy' objects), and among YLD from each environment ('genCOV_yy' object). Residuals covariances among YLD from each environment ('resCOV_yy' object) is also included

**Bug fixes**

  - Function 'fitBLUP' performs the new checking varU <= 2*var(y) to declare a possible error if FALSE
  - Function reshape2::melt is used instead of reshape::melt


### Version 0.4.0 (May-12-2021)

**New features**

  - More detailed functions' documentation
  - Function 'fitBLUP' performs a quality control for very small or negatives eigenvalues
  - Function 'saveBinary' does not save columns' nor rows' names anymore
  - Function 'SSI' uses now a C-based routine called 'add2diag' created to add a numeric value to the diagonal of a symmetric matrix (single or double precision). This routine is not at the user level
  - Function 'getGenCov' has the argument 'warn' to whether show warnings from 'BLUP' analyses

**Bug fixes**

  - All C-based routines: a 'long long' variable type, instead of an 'int' type, was used for indexing arrays (matrices). This change allows dealing with matrices whose length (number of rows x number of columns) exceed 2^31-1 = 2147483647 (e.g., a matrix of 46341 x 46341)  


### Version 0.3.0 (April-29-2021)

**Features**

- First released version
- Function 'solveMixed' (from GitHub version) was renamed to 'fitBLUP'
