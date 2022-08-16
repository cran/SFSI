#include <R.h>
#include <stdio.h>
#include <Rmath.h>
#include <Rinternals.h>
#include <Rdefines.h>
#include <string.h>
#include <R_ext/Lapack.h>
#include <float/float32.h>

//#define FLOAT(x) ((float*) INTEGER(x))

// ----------------------------------------------------------
// Update betas until convergence, i.e., the difference
// between two consecutive solutions is smaller than a threshold
// for a set of decreasing lambda: l0>l1>...>lk
// Using the whole matrix XtX as input stored in object XL passed as vector
//
//      n:         Number of beta parameters
//      XtX:       Crossprod matrix X'X in vector form, stacked by columns
//      Xty:       Crossprod vector X'y in vector form
//      q:         Number of lambdas
//      lambda:    Vector of lambdas for which the betas will be calculated
//      alpha:     alpha value in the Elastic-net problem
//      tol:       Maximum value between two consecutive solutions for beta to be accepted for convergence
//      maxiter:   Number of iterations to run before the updating stops
// ----------------------------------------------------------
SEXP updatebeta(SEXP n_, SEXP XtX_, SEXP Xty_, SEXP q_, SEXP lambda_, SEXP alpha_, SEXP tol_,
            SEXP maxiter_, SEXP dfmax_, SEXP isfloat_, SEXP scale_, SEXP sd_, SEXP verbose_)
{
    double *lambda, *sd, alpha, tol;
    double L1, L2, error;
    int i, j, k, n, q, maxiter, iter, verbose, scale, isfloat, dfmax;
    int *df;
    int inc1=1;
    double delta, bOLS, bNew;
    double *B, *b, *currfit;
    double value;
    SEXP list, df_, B_;

    n=INTEGER_VALUE(n_);
    q=INTEGER_VALUE(q_);
    maxiter=INTEGER_VALUE(maxiter_);
    dfmax=INTEGER_VALUE(dfmax_);
    verbose=asLogical(verbose_);
    scale=asLogical(scale_);
    alpha=NUMERIC_VALUE(alpha_);
    tol=NUMERIC_VALUE(tol_);
    isfloat=asLogical(isfloat_);

    PROTECT(lambda_=AS_NUMERIC(lambda_));
    lambda=NUMERIC_POINTER(lambda_);

    PROTECT(sd_=AS_NUMERIC(sd_));
    sd=NUMERIC_POINTER(sd_);

    df_=PROTECT(allocVector(INTSXP, q));
    df=INTEGER_POINTER(df_);

    B_ = PROTECT(allocMatrix(REALSXP, n, q)); // dimension N predictors x N lambdas
    B=NUMERIC_POINTER(B_);

    b=(double *) R_alloc(n, sizeof(double));
    currfit=(double *) R_alloc(n, sizeof(double));

    memset(currfit, 0, sizeof(double)*n);  // Initialize all currentfit to zero
    memset(b, 0, sizeof(double)*n);  // Initialize all coefficients to zero
    for(j=0; j<q; j++) df[j]=n;

    //double eps = DBL_EPSILON;

    if(isfloat){
      PROTECT(XtX_=AS_INTEGER(XtX_));
      float *XtX=FLOAT(XtX_);

      PROTECT(Xty_=AS_INTEGER(Xty_));
      float *Xty=FLOAT(Xty_);

      for(k=0; k<q; k++)
      {
          L1=alpha*lambda[k];
          L2=(1-alpha)*lambda[k];
          iter=0;
          error=INFINITY;    // Set to INF to enter to the WHILE
          while(iter<maxiter && error>tol)
          {
              iter++;
              error=0;
              for(j=0; j<n; j++)
              {
                  //bOLS=(Xty[j] - currfit)/XtX[n*j + j];
                  bOLS=Xty[j] - currfit[j];
                  if(fabs(bOLS) > L1){
                      bNew=sign(bOLS)*(fabs(bOLS)-L1)/(1+L2); // (XtX[n*j + j]+L2)
                  }else{
                      bNew=0;
                  }

                  delta = bNew-b[j];
                  if(fabs(delta)>0){
                      // update the current fit for all variables: cf=cf+XtX[j]*(bNew-b[j])
                    	for(i=0; i<n; i++){
                        currfit[i] = currfit[i] + delta*XtX[(long long)j*(long long)n + (long long)i];
                    	}

                      currfit[j] -= delta; //delta*XtX[n*j + j]
                      if(fabs(delta)>error){
                          error=fabs(delta);
                      }
                  }
                  b[j]=bNew;
              }
          }
          if(verbose){
              Rprintf(" lambda[%d]=%f.\t  nIters=%d.\t  Error=%G\n",k+1,lambda[k],iter,error);
              if(error>tol){
                Rprintf("    Warning: The process did not converge after %d iterations for lambda[%d]=%f\n",maxiter,k+1,lambda[k]);
              }
          }

          F77_NAME(dcopy)(&n, b, &inc1, B + k*n, &inc1);
          df[k]=0;
          for(j=0; j<n; j++){
            if(fabs(B[k*n + j])>0) df[k]++;
          }

          if(dfmax<n && df[k]>=dfmax) break;
      }
    }else{
      PROTECT(XtX_=AS_NUMERIC(XtX_));
      double *XtX=NUMERIC_POINTER(XtX_);

      PROTECT(Xty_=AS_NUMERIC(Xty_));
      double *Xty=NUMERIC_POINTER(Xty_);

      for(k=0; k<q; k++)
      {
          L1=alpha*lambda[k];
          L2=(1-alpha)*lambda[k];
          iter=0;
          error=INFINITY;    // Set to INF to enter to the WHILE
          while(iter<maxiter && error>tol)
          {
              iter++;
              error=0;
              for(j=0; j<n; j++)
              {
                  //bOLS=(Xty[j] - currfit)/XtX[n*j + j];
                  bOLS=Xty[j] - currfit[j];
                  if(fabs(bOLS) > L1){
                      bNew=sign(bOLS)*(fabs(bOLS)-L1)/(1+L2); // (XtX[n*j + j]+L2)
                  }else{
                      bNew=0;
                  }

                  delta = bNew-b[j];
                  if(fabs(delta)>0){
                      // update the current fit for all variables: cf=cf+XtX[j]*(bNew-b[j])
                      F77_NAME(daxpy)(&n, &delta, XtX + (long long)j*(long long)n, &inc1, currfit, &inc1);

                      currfit[j] -= delta; //delta*XtX[n*j + j]
                      if(fabs(delta)>error){
                          error=fabs(delta);
                      }
                  }
                  b[j]=bNew;
              }
          }
          if(verbose){
              Rprintf(" lambda[%d]=%f.\t  nIters=%d.\t  Error=%G\n",k+1,lambda[k],iter,error);
              if(error>tol){
                Rprintf("    Warning: The process did not converge after %d iterations for lambda[%d]=%f\n",maxiter,k+1,lambda[k]);
              }
          }

          F77_NAME(dcopy)(&n, b, &inc1, B + k*n, &inc1);
          df[k]=0;
          for(j=0; j<n; j++){
            if(fabs(B[k*n + j])>0) df[k]++;
          }

          if(dfmax<n && df[k]>=dfmax) break;
      }
    }
    if(scale){
      for(j=0; j<n; j++){
        value=1/sd[j];
        F77_NAME(dscal)(&q, &value, B + j, &n);
      }
    }

    // Creating a list with 1 vector elements:
    PROTECT(list = allocVector(VECSXP, 2));
    SET_VECTOR_ELT(list, 0, B_);
    SET_VECTOR_ELT(list, 1, df_);

    UNPROTECT(7);

    return(list);
}

// ----------------------------------------------------------
// Transform a covariance matrix to a (squared) distance matrix
// The distance between variables i and j is
// d^2_ij = v_ii + v_jj -2*v_ij
// where v_ii is the variance (diagonal value) of the variable i.
//
//      n:       Number of variables (columns in V)
//      V:       Covariance matrix
// ----------------------------------------------------------
SEXP cov2distance(SEXP n_, SEXP X_, SEXP isfloat_)
{
    double *di; // diagonal values
    int n;
    int i, j, isfloat;
    long long pos;

    n=INTEGER_VALUE(n_);
    isfloat=asLogical(isfloat_);

    di=(double *) R_alloc(n, sizeof(double));

    if(isfloat){
      PROTECT(X_=AS_INTEGER(X_));
      float *X=FLOAT(X_);

      // Diagonal values
      for(j=0; j<n; j++){
          di[j]=X[(long long)n*(long long)j + (long long)j];
          X[(long long)n*(long long)j + (long long)j]=0;
      }
      // Calculate distance  V[i,j] <- V[i,i] + V[j,j] -2*V[i,j]
      for(j=0; j<n-1; j++){
        for(i=j+1; i<n; i++)
        {
            pos=(long long)n*(long long)j + (long long)i;
            X[pos]=di[i] + di[j] -2*X[pos];
            pos=(long long)n*(long long)i + (long long)j;
            X[pos]=di[i] + di[j] -2*X[pos];
        }
      }

    }else{
      double *one;
      double two=-2;
      int inc0=0, inc1=1;

      PROTECT(X_=AS_NUMERIC(X_));
      double *X=NUMERIC_POINTER(X_);

      one=(double *) R_alloc(1, sizeof(double));
      one[0] =1;

      // Diagonal values
      for(j=0; j<n; j++){
          di[j]=X[(long long)n*(long long)j + (long long)j];
      }
      // Calculate distance  V[i,j] <- V[i,i] + V[j,j] -2*V[i,j]
      for(j=0; j<n; j++){
        F77_NAME(dscal)(&n, &two, X + (long long)n*(long long)j, &inc1);
        F77_NAME(daxpy)(&n, di + j, one, &inc0, X + (long long)n*(long long)j, &inc1);
        F77_NAME(daxpy)(&n, one, di, &inc1, X + (long long)n*(long long)j, &inc1);
        X[(long long)n*(long long)j + (long long)j]=0;
      }
    }

    // Creating a NULL variable with 1 elements:
    SEXP out = PROTECT(allocVector(NILSXP, 1));

    UNPROTECT(2);

    return(out);
}

// ----------------------------------------------------------
// Transform a covariance matrix to a correlation matrix
// The correlation between variables i and j is
// r_ij = v_ij/(sqrt(v_ii)*sqrt(v_jj))
// where sqrt(v_ii) is the SD of the variable i.
//
//      n:     Number of variables (columns in X)
//      X:     Covariance matrix
//      a:     Numeric value to multiply the resulting matrix
// ----------------------------------------------------------
SEXP cov2correlation(SEXP n_, SEXP X_, SEXP isfloat_, SEXP a_)
{
    double *sd;  // (Inverse) standard deviation (diagonal values)
    double a;
    double value;
    int n, nOK;
    int inc1=1;
    int i, j, isfloat;
    long long pos;
    SEXP list;

    n=INTEGER_VALUE(n_);
    isfloat=asLogical(isfloat_);
    a=NUMERIC_VALUE(a_);

    nOK=0;
    sd=(double *) R_alloc(n, sizeof(double));

    if(isfloat){
      PROTECT(X_=AS_INTEGER(X_));
      float *X=FLOAT(X_);

      for(j=0; j<n; j++)
      {
         sd[j]=sqrt(1/X[(long long)n*(long long)j + (long long)j]);
         X[(long long)n*(long long)j + (long long)j]=a*1;
         nOK=nOK+isfinite(sd[j]);
      }
      for(j=0; j<n-1; j++){
         for(i=j+1; i<n; i++){
           pos=(long long)n*(long long)j + (long long)i;
           X[pos]=a*X[pos]*sd[j]*sd[i];

           pos=(long long)n*(long long)i + (long long)j;
           X[pos]=a*X[pos]*sd[j]*sd[i];
         }
      }

    }else{
      PROTECT(X_=AS_NUMERIC(X_));
      double *X=NUMERIC_POINTER(X_);

      for(j=0; j<n; j++)
      {
         sd[j]=sqrt(1/X[(long long)n*(long long)j + (long long)j]);
         value=a*sd[j];
         F77_NAME(dscal)(&n, &value, X + (long long)n*(long long)j, &inc1);
         F77_NAME(dscal)(&n, sd+j, X + j, &n);
         X[(long long)n*(long long)j + (long long)j]=a*1;
         nOK=nOK+isfinite(sd[j]);
      }
    }

    PROTECT(list = allocVector(VECSXP, 1));

    SET_VECTOR_ELT(list, 0, ScalarInteger(nOK));

    UNPROTECT(2);

    return(list);
}

// ----------------------------------------------------------
// Add a numeric value to the diagonal of a squared matrix
//        n:  Number of variables (columns in X)
//    value:  Numeric value to be added
//        X:  squared matrix
// ----------------------------------------------------------
SEXP addvalue2diag(SEXP n_, SEXP X_, SEXP value_, SEXP isfloat_)
{
    double value;
    double *one;
    int inc0, inc1;
    int n;
    int i, isfloat;
    long long pos;

    n=INTEGER_VALUE(n_);
    value=NUMERIC_VALUE(value_);
    isfloat=asLogical(isfloat_);

    one=(double *) R_alloc(1, sizeof(double));
    one[0] =1;

    if(isfloat){
      PROTECT(X_=AS_INTEGER(X_));
      float *X=FLOAT(X_);

      for(i=0; i<n; i++)
      {
        pos=(long long)n*(long long)i + (long long)i;
        X[pos]=X[pos]+value;
      }
    }else{
      PROTECT(X_=AS_NUMERIC(X_));
      double *X=NUMERIC_POINTER(X_);

      inc0=0;
      inc1 = n+1;
      F77_NAME(daxpy)(&n, &value, one, &inc0, X, &inc1);
    }

    // Creating a NULL variable with 1 elements:
    SEXP out = PROTECT(allocVector(NILSXP, 1));

    UNPROTECT(2);

    return(out);
}

// ----------------------------------------------------------
// The p x p-1 matrix R has been formed from a
// p x p upper-triangular matrix by deleting column k
// integer p,k,n,nz; double precision R[p,p-1], z[p,nz]
// ----------------------------------------------------------
SEXP delete_col(SEXP R, SEXP p0, SEXP k0, SEXP z, SEXP nz0, SEXP flagfloat)
{
    double *pR2, *pz2;
    float *pR1, *pz1;
    int p, k, nz;
    int i,j, isFloat;
    double a, b, c, s, tau;
    //float a1, b1, c1, s1, tau1;
    long long pos1, pos2;
    SEXP list;

    p=INTEGER_VALUE(p0);
    k=INTEGER_VALUE(k0);
    nz=INTEGER_VALUE(nz0);
    isFloat=asLogical(flagfloat);

    if(isFloat){
      PROTECT(R=AS_INTEGER(R));
      pR1=FLOAT(R);
      pR2=(double *) R_alloc(0, sizeof(double));   // will not be used

      PROTECT(z=AS_INTEGER(z));
      pz1=FLOAT(z);
      pz2=(double *) R_alloc(0, sizeof(double));   // will not be used
    }else{
      PROTECT(R=AS_NUMERIC(R));
      pR2=NUMERIC_POINTER(R);
      pR1=(float *) R_alloc(0, sizeof(float));   // will not be used

      PROTECT(z=AS_NUMERIC(z));
      pz2=NUMERIC_POINTER(z);
      pz1=(float *) R_alloc(0, sizeof(float));   // will not be used
    }

    for(i=k-1; i<p-1; i++) // loop thats goes j=i+1,...,p-1
    {
        pos1=(long long)p*(long long)i + (long long)i;
        if(isFloat){
          a=pR1[pos1];
          b=pR1[pos1 + 1];
        }else{
          a=pR2[pos1];
          b=pR2[pos1 + 1];
        }
        if(b!=0.0f)
        {
          // Compute the rotation
	         if(fabs(b)>fabs(a)){
             tau = -a/b;
             s = 1/sqrt(1 + tau*tau);
             c = s * tau;
           }else{
             tau = -b/a;
             c = 1/sqrt(1 + tau*tau);
             s = c * tau;
           }

           // update r and z
           if(isFloat){
             pR1[pos1]=c*a - s*b;
             pR1[pos1 + 1]=s*a + c*b;
           }else{
             pR2[pos1]=c*a - s*b;
             pR2[pos1 + 1]=s*a + c*b;
           }

           for(j=i+1; j<p-1; j++)  // loop thats goes j=i+1,...,p-1
           {
             pos2=(long long)p*(long long)j + (long long)i;
             if(isFloat){
               a=pR1[pos2];
               b=pR1[pos2 + 1];
               pR1[pos2]=c*a - s*b;
               pR1[pos2 + 1]=s*a + c*b;
             }else{
               a=pR2[pos2];
               b=pR2[pos2 + 1];
               pR2[pos2]=c*a - s*b;
               pR2[pos2 + 1]=s*a + c*b;
             }
           }
           for(j=0; j<nz; j++)  // loop thats goes j=1,...,nz
           {
             pos2=(long long)p*(long long)j + (long long)i;
              if(isFloat){
               a=pz1[pos2];
               b=pz1[pos2 + 1];
               pz1[pos2]=c*a - s*b;
               pz1[pos2 + 1]=s*a + c*b;
             }else{
               a=pz2[pos2];
               b=pz2[pos2 + 1];
               pz2[pos2]=c*a - s*b;
               pz2[pos2 + 1]=s*a + c*b;
             }
           }
        }
    }

    // Creating a list with 1 vector elements:
    PROTECT(list = allocVector(VECSXP, 2));
    SET_VECTOR_ELT(list, 0, R);
    SET_VECTOR_ELT(list, 1, z);

    UNPROTECT(3);

    return(list);
}

// ----------------------------------------------------------
// vartype: 1:integer, 2: logical, 3: double
// ----------------------------------------------------------
SEXP writeBinFile(SEXP filename, SEXP n_, SEXP p_, SEXP X_, SEXP isfloat_, SEXP precision_)
{
    FILE *f=NULL;
    int i, j;
    int n, p, varsize, isfloat, precision, vartype;
    float valuefloat;
    SEXP list;

    n=INTEGER_VALUE(n_);
    p=INTEGER_VALUE(p_);
    isfloat=INTEGER_VALUE(isfloat_);
    precision=INTEGER_VALUE(precision_);

    varsize=0;  // Initialize value
    vartype=0;  // Initialize value

    f=fopen(CHAR(STRING_ELT(filename,0)),"wb");
    fwrite(&n,4, 1, f);
    fwrite(&p,4, 1, f);
    fwrite(&isfloat, 4, 1, f);

    if(TYPEOF(X_) ==  INTSXP || TYPEOF(X_) ==  LGLSXP)
    {
      vartype=TYPEOF(X_) ==  INTSXP ? 1 : 2;

      PROTECT(X_=AS_INTEGER(X_));
      if(isfloat){
        float *X=FLOAT(X_);
        varsize=sizeof(X[0]);
        fwrite(&vartype, 4, 1, f);
        fwrite(&varsize, 4, 1, f);

        for(j=0; j<p; j++){
          fwrite(X+(long long)n*(long long)j, varsize, n, f);
        }
      }else{
        int *X=INTEGER_POINTER(X_);
        varsize=sizeof(X[0]);
        fwrite(&vartype, 4, 1, f);
        fwrite(&varsize, 4, 1, f);

        for(j=0; j<p; j++){
          fwrite(X+(long long)n*(long long)j, varsize, n, f);
        }
      }

    }else{
      if(TYPEOF(X_) ==  REALSXP)
      {
        vartype=3;
        PROTECT(X_=AS_NUMERIC(X_));
        double *X=NUMERIC_POINTER(X_);

        if(precision==1)
        {
          varsize=sizeof(float);
          fwrite(&vartype, 4, 1, f);
          fwrite(&varsize, 4, 1, f);

          for(j=0; j<p; j++)
          {
            for(i=0; i<n; i++){
              valuefloat = X[(long long)n*(long long)j + (long long)i];
              fwrite(&valuefloat, varsize, 1, f);
            }
          }
        }else{
          varsize=sizeof(X[0]);
          fwrite(&vartype, 4, 1, f);
          fwrite(&varsize, 4, 1, f);

          for(j=0; j<p; j++){
            fwrite(X+(long long)n*(long long)j, varsize, n, f);
          }
        }

      }else{
        Rprintf("    File can not be saved with the current type format\n");
      }
    }

    fclose(f);

    PROTECT(list = allocVector(VECSXP, 5));

    SET_VECTOR_ELT(list, 0, ScalarInteger(n));
    SET_VECTOR_ELT(list, 1, ScalarInteger(p));
    SET_VECTOR_ELT(list, 2, ScalarInteger(isfloat));
    SET_VECTOR_ELT(list, 3, ScalarInteger(vartype));
    SET_VECTOR_ELT(list, 4, ScalarInteger(varsize));

    UNPROTECT(2);

    return(list);
}

// ----------------------------------------------------------
// vartype: 1:integer, 2: logical, 3: double
// ----------------------------------------------------------
SEXP readBinFile(SEXP filename, SEXP nsetrow_, SEXP nsetcol_,
                 SEXP maxsetrow_, SEXP maxsetcol_, SEXP setrow, SEXP setcol)
{
    FILE *f=NULL;
    int i, j;
    off_t file_length;
    off_t offset;
    int *psetrow, *psetcol;
    int nrow, ncol, varsize, vartype, nsetrow, nsetcol, isfloat;
    int maxsetcol, maxsetrow;
    int n, p, intvalue, nerror; // number of elements read

    SEXP list;
    SEXP X_=NULL;

    nsetrow=INTEGER_VALUE(nsetrow_);
    nsetcol=INTEGER_VALUE(nsetcol_);
    maxsetrow=INTEGER_VALUE(maxsetrow_);
    maxsetcol=INTEGER_VALUE(maxsetcol_);

    PROTECT(setrow=AS_INTEGER(setrow));
    psetrow=INTEGER_POINTER(setrow);

    PROTECT(setcol=AS_INTEGER(setcol));
    psetcol=INTEGER_POINTER(setcol);

    nerror=0;
    f=fopen(CHAR(STRING_ELT(filename,0)),"rb");

    intvalue=fread(&nrow, 4, 1, f);
    intvalue+=fread(&ncol, 4, 1, f);
    intvalue+=fread(&isfloat, 4, 1, f);
    intvalue+=fread(&vartype, 4, 1, f);
    intvalue+=fread(&varsize, 4, 1, f);

    if(intvalue < 5){
      Rprintf("    Error: The function failed to read data information\n");
      nerror++;
    }

    // Check if any index is larger than n or p
    if(nsetrow > 0 && maxsetrow>nrow){
       Rprintf("    Error in reading row %d: file contains only %d rows \n",maxsetrow,nrow);
       nerror++;
    }
    if(nsetcol > 0 && maxsetcol>ncol){
       Rprintf("    Error in reading column %d: file contains only %d columns \n",maxsetcol,ncol);
       nerror++;
    }

    n=nsetrow > 0 ? nsetrow : nrow;
    p=nsetcol > 0 ? nsetcol : ncol;

    fseeko(f, 0, SEEK_END);
    file_length=ftello(f);
    //Rprintf("    file_length= %lld \n",file_length);
    offset=(long long)nrow*(long long)ncol*(long long)varsize + 20;
    //Rprintf("    Offset= %lld \n",offset);

    if(nerror == 0 && offset==file_length)
    {
      fseeko(f, 20, SEEK_SET);
      if(varsize == sizeof(int) || varsize == sizeof(float)){ // NUMERIC SINGLE or INTEGER|LOGICAL
        X_=PROTECT(Rf_allocMatrix(INTSXP, n, p));

        if(isfloat || (vartype==3 && varsize == sizeof(float))){
          float *line;
          float *X=FLOAT(X_);
          if(nsetrow>0){
           line=(float *) R_alloc(nrow, sizeof(float));
          }

          for(j=0; j<p; j++)
          {
            if(nsetcol > 0){
              fseeko(f, 20 + (long long)nrow*(long long)varsize*((long long)psetcol[j]-1), SEEK_SET);
            }

            if(nsetrow>0){
               intvalue=fread(line, varsize, nrow, f);
               for(i=0; i<n; i++){
                 X[(long long)n*(long long)j + (long long)i]=line[psetrow[i]-1];
               }
            }else{
               intvalue=fread(X + (long long)n*(long long)j, varsize, nrow, f);
            }
            if(intvalue<nrow){
               Rprintf("    Error: The function failed to read data at col %d \n",j+1);
               nerror++;
            }
            if(nerror>0){
              break;
            }
          }

        }else{
          int *line;
          int *X=INTEGER_POINTER(X_);
          if(nsetrow>0){
           line=(int *) R_alloc(nrow, sizeof(int));
          }

          for(j=0; j<p; j++)
          {
            if(nsetcol > 0){
              fseeko(f, 20 + (long long)nrow*(long long)varsize*((long long)psetcol[j]-1), SEEK_SET);
            }

            if(nsetrow>0){
               intvalue=fread(line, varsize, nrow, f);
               for(i=0; i<n; i++){
                 X[(long long)n*(long long)j + (long long)i]=line[psetrow[i]-1];
               }
            }else{
               intvalue=fread(X + (long long)n*(long long)j, varsize, nrow, f);
            }
            if(intvalue<nrow){
               Rprintf("    Error: The function failed to read data at col %d \n",j+1);
               nerror++;
            }
            if(nerror>0){
              break;
            }
          }
        }

      }else{
        if(varsize == sizeof(double)){ // DOUBLE
          double *line;
          X_=PROTECT(Rf_allocMatrix(REALSXP, n, p));
          double *X=NUMERIC_POINTER(X_);

          if(nsetrow>0){
           line=(double *) R_alloc(nrow, sizeof(double));
          }

          for(j=0; j<p; j++)
          {
            if(nsetcol > 0){
              fseeko(f, 20 + (long long)nrow*(long long)varsize*((long long)psetcol[j]-1), SEEK_SET);
            }

            if(nsetrow>0){
               intvalue=fread(line, varsize, nrow, f);
               for(i=0; i<n; i++){
                 X[(long long)n*(long long)j + (long long)i]=line[psetrow[i]-1];
               }
            }else{
               intvalue=fread(X + (long long)n*(long long)j, varsize, nrow, f);
            }
            if(intvalue<nrow){
               Rprintf("    Error: The function failed to read data at col %d \n",j+1);
               nerror++;
            }
            if(nerror>0){
              break;
            }
          }
        }else{
          Rprintf("    Error: File can not be read with the current type format\n");
          nerror++;
        }
      }
    }else{
      Rprintf("    Error: The function failed to read data from file \n");
      nerror++;
    }

    fclose(f);

    PROTECT(list = allocVector(VECSXP, 7));

    SET_VECTOR_ELT(list, 0, ScalarInteger(n));
    SET_VECTOR_ELT(list, 1, ScalarInteger(p));
    SET_VECTOR_ELT(list, 2, ScalarInteger(isfloat));
    SET_VECTOR_ELT(list, 3, ScalarInteger(vartype));
    SET_VECTOR_ELT(list, 4, ScalarInteger(varsize));
    SET_VECTOR_ELT(list, 5, ScalarInteger(nerror));
    SET_VECTOR_ELT(list, 6, X_);

    UNPROTECT(4);
    return(list);
}
