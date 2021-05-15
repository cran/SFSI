# Version 0.4.0 (May-12-2021)

## New features

  - More detailed functions' documentation
  - Function 'fitBLUP' performs a quality control for very small or negatives eigenvalues
  - Function 'saveBinary' does not save columns' nor rows' names anymore
  - Function 'SSI' uses now a C-based routine called 'add2diag' created to add a numeric value to the diagonal of a symmetric matrix (single or double precision). This routine is not at the user level
  - Function 'getGenCov' has the argument 'warn' to whether show warnings from 'BLUP' analyses

## Bug fixes

  - All C-based routines: a 'long long' variable type, instead of an 'int' type, was used for indexing arrays (matrices). This change allows dealing with matrices whose length (number of rows x number of columns) exceed 2^31-1 = 2147483647 (e.g., a matrix of 46341 x 46341)  


# Version 0.3.0 (April-29-2021)

## Features

- First released version
