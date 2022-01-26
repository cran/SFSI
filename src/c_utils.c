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
//      a:         alpha value in the Elastic-net problem
//      maxtole:   Maximum value between two consecutive solutions for beta to be accepted for convergence
//      maxsteps:  Number of iterations to run before the updating stops
// ----------------------------------------------------------
SEXP updatebeta(SEXP n, SEXP XtX, SEXP Xty, SEXP q, SEXP lambda, SEXP a, SEXP maxtole, SEXP maxsteps, SEXP maxDF, SEXP echo, SEXP flagfloat)
{
    float *pXtX1, *pXty1;
    double *pXtX2, *pXty2, *plambda, alpha, maxTol, eps;
    double L1, L2, maxdiff;
    int i, j, k, np, maxIter, iter, nlambda, verbose, isFloat, maxdf;
    int *pDF;
    //int inc=1;
    double delta, bOLS, bNew;
    double *pB, *b, *currfit;
    long long pos1;
    SEXP list, DF, B;

    np=INTEGER_VALUE(n);
    nlambda=INTEGER_VALUE(q);
    maxIter=INTEGER_VALUE(maxsteps);
    maxdf=INTEGER_VALUE(maxDF);
    verbose=asLogical(echo);
    alpha=NUMERIC_VALUE(a);
    maxTol=NUMERIC_VALUE(maxtole);
    isFloat=asLogical(flagfloat);

    if(isFloat){
      PROTECT(XtX=AS_INTEGER(XtX));
      pXtX1=FLOAT(XtX);
      pXtX2=(double *) R_alloc(0, sizeof(double)); // will not be used

      PROTECT(Xty=AS_INTEGER(Xty));
      pXty1=FLOAT(Xty);
      pXty2=(double *) R_alloc(0, sizeof(double)); // will not be used
    }else{
      PROTECT(XtX=AS_NUMERIC(XtX));
      pXtX2=NUMERIC_POINTER(XtX);
      pXtX1=(float *) R_alloc(0, sizeof(float)); // will not be used

      PROTECT(Xty=AS_NUMERIC(Xty));
      pXty2=NUMERIC_POINTER(Xty);
      pXty1=(float *) R_alloc(0, sizeof(float)); // will not be used
	  }

    PROTECT(lambda=AS_NUMERIC(lambda));
    plambda=NUMERIC_POINTER(lambda);

    DF = PROTECT(allocVector(INTSXP, nlambda));
    pDF=INTEGER_POINTER(DF);

    B = PROTECT(allocMatrix(REALSXP, np, nlambda)); // dimension N predictors x N lambdas
    pB=NUMERIC_POINTER(B);

    b=(double *) R_alloc(np, sizeof(double));
    currfit=(double *) R_alloc(np, sizeof(double));

    memset(currfit,0, sizeof(double)*np);  // Initialize all currentfit to zero
    memset(b,0, sizeof(double)*np);  // Initialize all coefficients to zero
    for(j=0; j<nlambda; j++) pDF[j]=np;

    eps = DBL_EPSILON;

    for(k=0; k<nlambda; k++)
    {
        L1=alpha*plambda[k];
        L2=(1-alpha)*plambda[k];
        iter=0;
        maxdiff=INFINITY;    // Set to INF to enter to the WHILE
        while(iter<maxIter && maxdiff>maxTol-eps)
        {
            iter++;
            maxdiff=0;
            for(j=0; j<np; j++)
            {
                //bOLS=(pXty[j] - currfit)/pXtX[np*j + j];
                if(isFloat){
                	bOLS=pXty1[j] - currfit[j];
                }else{
                	bOLS=pXty2[j] - currfit[j];
                }
                if(fabs(bOLS) > L1){
                    bNew=sign(bOLS)*(fabs(bOLS)-L1)/(1+L2); // (pXtX[np*j + j]+L2)
                }else{
                    bNew=0;
                }

                delta = bNew-b[j];
                if(fabs(delta)>eps){
                    // update the current fit for all variables: cf=cf+XtX[j]*(bNew-b[j])
                    //F77_NAME(daxpy)(&np, &delta, pXtX + j*np, &inc, currfit, &inc);
                    if(isFloat){
                    	for(i=0; i<np; i++)
                    	{ pos1=(long long)j*(long long)np + (long long)i;
                        currfit[i] = currfit[i] + delta*pXtX1[pos1];
                    	}
                    }else{
                    	for(i=0; i<np; i++)
                    	{ pos1=(long long)j*(long long)np + (long long)i;
                        currfit[i] = currfit[i] + delta*pXtX2[pos1];
                    	}
                    }
                    currfit[j] -= delta; //delta*pXtX[np*j + j]
                    if(fabs(delta)>maxdiff){
                        //Rprintf(" iter=%d. var=%d. abs(delta)=%g.  \t maxdiff=%g.  \t tol=%g\n",iter,j+1,fabs(delta),maxdiff,maxTol);
                        maxdiff=fabs(delta);
                    }
                }
                b[j]=bNew;
            }
        }
        if(verbose){
            Rprintf(" lambda[%d]=%f.\t  nIters=%d.\t  maxError=%g\n",k+1,plambda[k],iter,maxdiff);
            if(maxdiff>maxTol){
              Rprintf("    Warning: The process did not converge after %d iterations for lambda[%d]=%f\n",maxIter,k+1,plambda[k]);
            }
        }
        //F77_NAME(dcopy)(&np, b, &inc, pB+k, &nlambda);
        pDF[k]=0;
        for(j=0; j<np; j++){
          if(fabs(b[j])>eps) pDF[k]++;
          pB[k*np+j]=b[j];
          //pB[k+j*nlambda]=b[j];
        }
        if(maxdf<np && pDF[k]>=maxdf) break;
    }

    // Creating a list with 1 vector elements:
    PROTECT(list = allocVector(VECSXP, 2));
    SET_VECTOR_ELT(list, 0, B);
    SET_VECTOR_ELT(list, 1, DF);

    UNPROTECT(6);

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
SEXP cov2distance(SEXP n, SEXP V, SEXP flagfloat)
{
    float *pV1, *pd1;
    double *pV2, *pd2;
    int np;
    int i, j, isFloat;
    long long pos1;

    np=INTEGER_VALUE(n);
    isFloat=asLogical(flagfloat);

    if(isFloat){
      PROTECT(V=AS_INTEGER(V));
      pV1=FLOAT(V);
      pd1=(float *) R_alloc(np, sizeof(float));

      pV2=(double *) R_alloc(0, sizeof(double)); // will not be used
      pd2=(double *) R_alloc(0, sizeof(double)); // will not be used
    }else{
      PROTECT(V=AS_NUMERIC(V));
      pV2=NUMERIC_POINTER(V);
      pd2=(double *) R_alloc(np, sizeof(double));

      pV1=(float *) R_alloc(0, sizeof(float));   // will not be used
      pd1=(float *) R_alloc(0, sizeof(float));   // will not be used
    }

    // Diagonal values
    for(j=0; j<np; j++)
    {
        pos1=(long long)np*(long long)j + (long long)j;
        if(isFloat){
          pd1[j]=pV1[pos1];
          pV1[pos1]=0;
        }else{
          pd2[j]=pV2[pos1];
          pV2[pos1]=0;
        }
    }

    for(j=0; j<np-1; j++)
    {
        for(i=j+1; i<np; i++)
        {
          if(isFloat){
            pos1=(long long)np*(long long)j + (long long)i;
            pV1[pos1]=pd1[i] + pd1[j] -2*pV1[pos1];
            pos1=(long long)np*(long long)i + (long long)j;
            pV1[pos1]=pd1[i] + pd1[j] -2*pV1[pos1];
          }else{
            pos1=(long long)np*(long long)j + (long long)i;
            pV2[pos1]=pd2[i] + pd2[j] -2*pV2[pos1];
            pos1=(long long)np*(long long)i + (long long)j;
            pV2[pos1]=pd2[i] + pd2[j] -2*pV2[pos1];
          }
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
//      n:       Number of variables (columns in V)
//      V:       Covariance matrix
// ----------------------------------------------------------
SEXP cov2correlation(SEXP n, SEXP V, SEXP flagfloat, SEXP a)
{
    float *pV1, *psdx1;
    double *pV2, *psdx2;
    double a0;
    int np, nOK;
    int i, j, isFloat;
    long long pos1;
    SEXP list;

    np=INTEGER_VALUE(n);
    isFloat=asLogical(flagfloat);
    a0=NUMERIC_VALUE(a);

    if(isFloat){
      PROTECT(V=AS_INTEGER(V));
      pV1=FLOAT(V);
      psdx1=(float *) R_alloc(np, sizeof(float));

      pV2=(double *) R_alloc(0, sizeof(double)); // will not be used
      psdx2=(double *) R_alloc(0, sizeof(double)); // will not be used
    }else{
      PROTECT(V=AS_NUMERIC(V));
      pV2=NUMERIC_POINTER(V);
      psdx2=(double *) R_alloc(np, sizeof(double));

      pV1=(float *) R_alloc(0, sizeof(float));   // will not be used
      psdx1=(float *) R_alloc(0, sizeof(float));   // will not be used
    }

    // Get standard deviations
    nOK=0;
    for(i=0; i<np; i++)
    {
        pos1=(long long)np*(long long)i + (long long)i;
        if(isFloat){
          psdx1[i]=sqrt(pV1[pos1]);
          pV1[pos1]=a0*1;
          nOK=nOK+isfinite(1/psdx1[i]);
        }else{
          psdx2[i]=sqrt(pV2[pos1]);
          pV2[pos1]=a0*1;
          nOK=nOK+isfinite(1/psdx2[i]);
        }
    }

    for(j=0; j<np-1; j++)
    {
        for(i=j+1; i<np; i++)
        {
          if(isFloat){
            pos1=(long long)np*(long long)j + (long long)i;
            pV1[pos1]=a0*pV1[pos1]/(psdx1[j]*psdx1[i]);
            pos1=(long long)np*(long long)i + (long long)j;
            pV1[pos1]=a0*pV1[pos1]/(psdx1[j]*psdx1[i]);
          }else{
            pos1=(long long)np*(long long)j + (long long)i;
            pV2[pos1]=a0*pV2[pos1]/(psdx2[j]*psdx2[i]);
            pos1=(long long)np*(long long)i + (long long)j;
            pV2[pos1]=a0*pV2[pos1]/(psdx2[j]*psdx2[i]);
          }
        }
    }

    PROTECT(list = allocVector(VECSXP, 1));

    SET_VECTOR_ELT(list, 0, ScalarInteger(nOK));

    UNPROTECT(2);

    return(list);
}

// ----------------------------------------------------------
// Add a numeric value to the diagonal of a squared matrix
// a: numeric value to be added
// V: squared matrix
// ----------------------------------------------------------
SEXP addvalue2diag(SEXP n, SEXP V, SEXP a, SEXP flagfloat)
{
    float *pV1;
    double *pV2;
    double value;
    int np;
    int i, isFloat;
    long long pos1;

    np=INTEGER_VALUE(n);
    value=NUMERIC_VALUE(a);
    isFloat=asLogical(flagfloat);

    if(isFloat){
      PROTECT(V=AS_INTEGER(V));
      pV1=FLOAT(V);
      pV2=(double *) R_alloc(0, sizeof(double));   // will not be used
    }else{
      PROTECT(V=AS_NUMERIC(V));
      pV2=NUMERIC_POINTER(V);
      pV1=(float *) R_alloc(0, sizeof(float));   // will not be used
    }

    // Diagonal values
    for(i=0; i<np; i++)
    {
        pos1=(long long)np*(long long)i + (long long)i;
        if(isFloat){
          pV1[pos1]=pV1[pos1]+value;
        }else{
          pV2[pos1]=pV2[pos1]+value;
        }
    }

    // Creating a NULL variable with 1 elements:
    SEXP out = PROTECT(allocVector(NILSXP, 1));

    UNPROTECT(2);

    return(out);
}

// ----------------------------------------------------------
// Return the index of variables that are very correlated
// output a 0-1 vector indicating whether the variable has
// a correlation of at most a threshold
//
//      COV:         Covariance matrix, stacked by columns
//      maxVal:      Maximum correlation allowed
// ----------------------------------------------------------
SEXP getCorrelated(SEXP n, SEXP COV, SEXP maxVal)
{
    double *pCOV, *psdx;
    double maxCor, corre;
    int *pout;
    int np,isok, cont;
    int i,j;
    long long pos1;
    SEXP list;

    np=INTEGER_VALUE(n);
    maxCor=NUMERIC_VALUE(maxVal);

    PROTECT(COV=AS_NUMERIC(COV));
    pCOV=NUMERIC_POINTER(COV);

    SEXP out = PROTECT(allocVector(INTSXP, np));
    pout=INTEGER_POINTER(out);

    psdx=(double *) R_alloc(np, sizeof(double));

    // Get standard deviations
    for(i=0; i<np; i++)
    {
        pos1=(long long)np*(long long)i + (long long)i;
        psdx[i]=sqrt(pCOV[pos1]);
    }

    cont=0;
    for(j=0; j<np-1; j++)
    {
        if(psdx[j]>0)
        {
            isok=1;
            for(i=j+1; i<np; i++)
            {
                if(psdx[i]>0){
                  pos1=(long long)np*(long long)j + (long long)i;
                    corre=pCOV[pos1]/(psdx[j]*psdx[i]);
                    if(corre>maxCor){
                        isok=0;
                        break;
                    }
                }
            }

            if(isok==1){
              pout[cont]=j+1;
              cont++;
            }
        }
    }

    if(psdx[np-1]>0){
        pout[cont]=np;
        cont++;
    }

    SEXP le = PROTECT(ScalarInteger(cont));

    // Creating a list with 1 vector elements:
    PROTECT(list = allocVector(VECSXP, 2));
    // Attaching outputs to list:
    SET_VECTOR_ELT(list, 0, out);
    SET_VECTOR_ELT(list, 1, le);

    UNPROTECT(4);

    return(list);
}

// ----------------------------------------------------------
// the p x p-1 matrix R has been formed from a
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
// ----------------------------------------------------------
SEXP writeBinFileFloat(SEXP filename, SEXP n, SEXP p, SEXP size, SEXP X, SEXP flagfloat)
{
    FILE *f=NULL;
    int i, j;
    long long pos1;
    int nrows, ncols, sizevar, isFloat;
    //int inc=1;
    float *pX1;
    double *pX2;
    //double *linedouble;
    float valuefloat;
    SEXP list;

    nrows=INTEGER_VALUE(n);
    ncols=INTEGER_VALUE(p);
    sizevar=INTEGER_VALUE(size);
    isFloat=asLogical(flagfloat);

    if(isFloat){
      PROTECT(X=AS_INTEGER(X));
      pX1=FLOAT(X);
      pX2=(double *) R_alloc(0, sizeof(double));   // will not be used
    }else{
      PROTECT(X=AS_NUMERIC(X));
      pX2=NUMERIC_POINTER(X);
      pX1=(float *) R_alloc(0, sizeof(float));   // will not be used
    }

    f=fopen(CHAR(STRING_ELT(filename,0)),"wb");
    fwrite(&nrows,4, 1 , f);
    fwrite(&ncols,4, 1 , f);
    fwrite(&sizevar,4, 1 , f);
    fwrite(&isFloat,4, 1 , f);

    // Write lines
    for(i=0; i<nrows; i++)
    {
      for(j=0; j<ncols; j++)
      {
        pos1=(long long)nrows*(long long)j + (long long)i;
        if(sizevar==4)
        {
          if(isFloat){
            fwrite(pX1+pos1,sizevar, 1 , f);
          }else{
            valuefloat = pX2[pos1];
            fwrite(&valuefloat,sizevar, 1 , f);
          }
        }else{
          fwrite(pX2+pos1,sizevar, 1 , f);
        }
      }
    }

    fclose(f);

    PROTECT(list = allocVector(VECSXP, 4));

    SET_VECTOR_ELT(list, 0, ScalarInteger(nrows));
    SET_VECTOR_ELT(list, 1, ScalarInteger(ncols));
    SET_VECTOR_ELT(list, 2, ScalarInteger(sizevar));
    SET_VECTOR_ELT(list, 3, ScalarInteger(isFloat));

    UNPROTECT(2);

    return(list);
}

// ----------------------------------------------------------
// ----------------------------------------------------------
SEXP readBinFileFloat(SEXP filename, SEXP nsetRow, SEXP nsetCol, SEXP setRow, SEXP setCol)
{
    FILE *f=NULL;
    int i, j;
    //int k;
    off_t file_length;
    off_t offset;
    int *psetRow, *psetCol;
    int nrows, ncols, sizevar, lsrow, lscol, isFloat;
    int n, p, intvalue, nerror; // number of elements read
    float *pX1;
    float  *linefloat;
    double *pX2;
    double *linedouble;
    long long pos1;
    SEXP list, X;

    lsrow=INTEGER_VALUE(nsetRow);
    lscol=INTEGER_VALUE(nsetCol);

    PROTECT(setRow=AS_INTEGER(setRow));
    psetRow=INTEGER_POINTER(setRow);

    PROTECT(setCol=AS_INTEGER(setCol));
    psetCol=INTEGER_POINTER(setCol);

    nerror=0;
    f=fopen(CHAR(STRING_ELT(filename,0)),"rb");
    intvalue=fread(&nrows, 4, 1, f);
    if(intvalue<1){
      Rprintf("    Error: The function failed to read information on the number of rows\n");
      nerror++;
    }
    intvalue=fread(&ncols, 4, 1, f);
    if(intvalue<1){
      Rprintf("    Error: The function failed to read information on the number of columns\n");
      nerror++;
    }
    intvalue=fread(&sizevar, 4, 1, f);
    if(intvalue<1){
      Rprintf("    Error: The function failed to read information on the size in bytes\n");
      nerror++;
    }
    intvalue=fread(&isFloat, 4, 1, f);
    if(intvalue<1){
      Rprintf("    Error: The function failed to read information on variable type\n");
      nerror++;
    }

    n=lsrow > 0 ? lsrow : nrows;
    p=lscol > 0 ? lscol : ncols;

    if(isFloat || sizevar==4){
      linedouble=(double *) R_alloc(0, sizeof(double));   // will not be used
      linefloat=(float *) R_alloc(ncols, sizeof(float));

      X=PROTECT(Rf_allocMatrix(INTSXP, n, p));
      pX1=FLOAT(X);
      pX2=(double *) R_alloc(0, sizeof(double));   // will not be used
    }else{
      linedouble=(double *) R_alloc(ncols, sizeof(double));
      linefloat=(float *) R_alloc(0, sizeof(float));   // will not be used

      X=PROTECT(Rf_allocMatrix(REALSXP, n, p));
      pX2=NUMERIC_POINTER(X);
      pX1=(float *) R_alloc(0, sizeof(float));   // will not be used
    }


    fseeko(f, 0, SEEK_END);
    file_length=ftello(f);
    //Rprintf("    file_length= %lld \n",file_length);
    offset=(long long)nrows*(long long)ncols*(long long)sizevar + 16;
    //Rprintf("    Offset= %lld \n",offset);

    if(offset==file_length)
    {
      //Rprintf("    Reading lines... \n");
      fseeko(f, 16, SEEK_SET);
      for(i=0; i<n; i++)
      {
          if(lsrow > 0){
              // Move to the line indicated by setRow
            if(psetRow[i]<=nrows)
            {
              offset=16 + (long long)ncols*(long long)sizevar*((long long)psetRow[i]-1);
              intvalue=fseeko(f, offset, SEEK_SET);
              if(intvalue != 0){
                 Rprintf("    Error in line %d: fseek failed at offset=%lld \n",psetRow[i],offset);
                 nerror++;
              }
            }else{
              Rprintf("    Error in reading row %d: file contains only %d rows \n",psetRow[i],nrows);
              nerror++;
            }
          }

          if(isFloat || sizevar==4){
              intvalue=fread(linefloat,sizevar,ncols,f);
          }else{
              intvalue=fread(linedouble,sizevar,ncols,f);
          }
          if(intvalue<ncols){
            Rprintf("    Error: The function failed to read data at row %d \n",i+1);
            nerror++;
          }
          // Read columns
          for(j=0; j<p; j++)
          {
            pos1=(long long)n*(long long)j + (long long)i;
            if(isFloat || sizevar==4){
              if(lscol>0)
              {
                  if(psetCol[j]>ncols){
                    Rprintf("    Error in reading column %d: file contains only %d columns \n",psetCol[j],ncols);
                    nerror++;
                  }
                  pX1[pos1]=linefloat[psetCol[j]-1];
              }else{
                  pX1[pos1]=linefloat[j];
              }
            }else{
                if(lscol>0)
                {
                  if(psetCol[j]>ncols){
                    Rprintf("    Error in reading column %d: file contains only %d columns \n",psetCol[j],ncols);
                    nerror++;
                  }
                  pX2[pos1]=linedouble[psetCol[j]-1];
                }else{
                  pX2[pos1]=linedouble[j];
                }
            }
          }
          if(nerror>0){
            break;
          }
      }
    }else{
      Rprintf("    Error: The function failed to read data from file \n");
      nerror++;
    }

    fclose(f);

    PROTECT(list = allocVector(VECSXP, 6));
    // Attaching outputs to list:
    SET_VECTOR_ELT(list, 0, ScalarInteger(n));
    SET_VECTOR_ELT(list, 1, ScalarInteger(p));
    SET_VECTOR_ELT(list, 2, ScalarInteger(sizevar));
    SET_VECTOR_ELT(list, 3, ScalarInteger(isFloat));
    SET_VECTOR_ELT(list, 4, ScalarInteger(nerror));
    SET_VECTOR_ELT(list, 5, X);

    UNPROTECT(4);
    return(list);
}
