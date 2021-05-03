#include <R.h>
#include <stdio.h>
#include <Rmath.h>
#include <Rinternals.h>
#include <Rdefines.h>
#include <string.h>
#include <R_ext/Lapack.h>
#include <float/float32.h>

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
SEXP updatebeta(SEXP n, SEXP XtX, SEXP Xty, SEXP q, SEXP lambda, SEXP a, SEXP maxtole, SEXP maxsteps, SEXP echo, SEXP flagfloat)
{
    float *pXtX1, *pXty1;
    double *pXtX2, *pXty2, *plambda, alpha, maxTol, eps;
    double L1, L2, maxdiff;
    int i, j, k, np, maxIter, iter, nlambda, verbose, isFloat;
    //int inc=1;
    double delta, bOLS, bNew;
    double *pB, *b, *currfit;
    SEXP list, B;

    np=INTEGER_VALUE(n);
    nlambda=INTEGER_VALUE(q);
    maxIter=INTEGER_VALUE(maxsteps);
    verbose=asLogical(echo);
    alpha=NUMERIC_VALUE(a);
    maxTol=NUMERIC_VALUE(maxtole);
    isFloat=asLogical(flagfloat);

    if(isFloat){
      PROTECT(XtX=AS_INTEGER(XtX));
      pXtX1=FLOAT(XtX);

      PROTECT(Xty=AS_INTEGER(Xty));
      pXty1=FLOAT(Xty);
    }else{
      PROTECT(XtX=AS_NUMERIC(XtX));
      pXtX2=NUMERIC_POINTER(XtX);

      PROTECT(Xty=AS_NUMERIC(Xty));
      pXty2=NUMERIC_POINTER(Xty);
	}

    PROTECT(lambda=AS_NUMERIC(lambda));
    plambda=NUMERIC_POINTER(lambda);

    B = PROTECT(allocMatrix(REALSXP, nlambda, np));
    pB=NUMERIC_POINTER(B);

    b=(double *) R_alloc(np, sizeof(double));
    currfit=(double *) R_alloc(np, sizeof(double));

    memset(currfit,0, sizeof(double)*np);  // Initialize all currentfit to zero
    memset(b,0, sizeof(double)*np);  // Initialize all coefficients to zero

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
                    	{   currfit[i] = currfit[i] + delta*pXtX1[j*np+i];
                    	}
                    }else{
                    	for(i=0; i<np; i++)
                    	{   currfit[i] = currfit[i] + delta*pXtX2[j*np+i];
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
        for(i=0; i<np; i++){
          pB[k+i*nlambda]=b[i];
        }
    }

    // Creating a list with 1 vector elements:
    PROTECT(list = allocVector(VECSXP, 1));
    SET_VECTOR_ELT(list, 0, B);

    UNPROTECT(5);

    return(list);
}

// ----------------------------------------------------------
// ----------------------------------------------------------
SEXP cov2distance(SEXP n, SEXP V, SEXP flagfloat)
{
    float *pV1, *pd1;
    double *pV2, *pd2;
    int np;
    int i, j, isFloat ;

    np=INTEGER_VALUE(n);
    isFloat=asLogical(flagfloat);

    if(isFloat){
      PROTECT(V=AS_INTEGER(V));
      pV1=FLOAT(V);

      pd1=(float *) R_alloc(np, sizeof(float));
      pd2=(double *) R_alloc(0, sizeof(double)); // will not be used
    }else{
      PROTECT(V=AS_NUMERIC(V));
      pV2=NUMERIC_POINTER(V);

      pd1=(float *) R_alloc(0, sizeof(float));   // will not be used
      pd2=(double *) R_alloc(np, sizeof(double));
    }

    // Diagonal values
    for(j=0; j<np; j++)
    {
        if(isFloat){
          pd1[j]=pV1[np*j + j];
          pV1[np*j + j]=0;
        }else{
          pd2[j]=pV2[np*j + j];
          pV2[np*j + j]=0;
        }
    }

    for(j=0; j<np-1; j++)
    {
        for(i=j+1; i<np; i++)
        {
          if(isFloat){
            pV1[np*j + i]=pd1[i] + pd1[j] -2*pV1[np*j + i];
            pV1[np*i + j]=pd1[i] + pd1[j] -2*pV1[np*i + j];
          }else{
            pV2[np*j + i]=pd2[i] + pd2[j] -2*pV2[np*j + i];
            pV2[np*i + j]=pd2[i] + pd2[j] -2*pV2[np*i + j];
          }
        }
    }

    // Creating a NULL variable with 1 elements:
    SEXP out = PROTECT(allocVector(NILSXP, 1));

    UNPROTECT(2);

    return(out);
}

// ----------------------------------------------------------
// Scale XtX matrix to have all diagonal elements equal to one.
// Scaling is perform by dividing each entry x_ij by
// x_ij = x_ij/(sqrt(x_ii)*sqrt(x_jj))
// where sqrt(x_ii) is the SD of the variable i. These values are returned
//
//      n:         Number of variable (columns in XtX)
//      XtX:       Crossprod matrix X'X in vector form, stacked by columns
// ----------------------------------------------------------
SEXP cov2correlation(SEXP n, SEXP V, SEXP flagfloat)
{
    float *pV1, *psdx1;
    double *pV2, *psdx2;
    int np, nOK;
    int i, j, isFloat;
    SEXP list;

    np=INTEGER_VALUE(n);
    isFloat=asLogical(flagfloat);

    if(isFloat){
      PROTECT(V=AS_INTEGER(V));
      pV1=FLOAT(V);

      psdx1=(float *) R_alloc(np, sizeof(float));
      psdx2=(double *) R_alloc(0, sizeof(double)); // will not be used
    }else{
      PROTECT(V=AS_NUMERIC(V));
      pV2=NUMERIC_POINTER(V);

      psdx1=(float *) R_alloc(0, sizeof(float));   // will not be used
      psdx2=(double *) R_alloc(np, sizeof(double));
    }

    // Get standard deviations
    nOK=0;
    for(i=0; i<np; i++)
    {
        if(isFloat){
          psdx1[i]=sqrt(pV1[np*i + i]);
          pV1[np*i + i]=1;
          nOK=nOK+isfinite(1/psdx1[i]);
        }else{
          psdx2[i]=sqrt(pV2[np*i + i]);
          pV2[np*i + i]=1;
          nOK=nOK+isfinite(1/psdx2[i]);
        }
    }

    for(j=0; j<np-1; j++)
    {
        for(i=j+1; i<np; i++)
        {
          if(isFloat){
            pV1[np*j + i]=pV1[np*j + i]/(psdx1[j]*psdx1[i]);
            pV1[np*i + j]=pV1[np*i + j]/(psdx1[j]*psdx1[i]);
          }else{
            pV2[np*j + i]=pV2[np*j + i]/(psdx2[j]*psdx2[i]);
            pV2[np*i + j]=pV2[np*i + j]/(psdx2[j]*psdx2[i]);
          }
        }
    }

    PROTECT(list = allocVector(VECSXP, 1));

    SET_VECTOR_ELT(list, 0, ScalarInteger(nOK));

    UNPROTECT(2);

    return(list);
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
        psdx[i]=sqrt(pCOV[np*i + i]);
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
                    corre=pCOV[np*j + i]/(psdx[j]*psdx[i]);
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
    SEXP list;

    p=INTEGER_VALUE(p0);
    k=INTEGER_VALUE(k0);
    nz=INTEGER_VALUE(nz0);
    isFloat=asLogical(flagfloat);

    if(isFloat){
      PROTECT(R=AS_INTEGER(R));
      pR1=FLOAT(R);

      PROTECT(z=AS_INTEGER(z));
      pz1=FLOAT(z);
    }else{
      PROTECT(R=AS_NUMERIC(R));
      pR2=NUMERIC_POINTER(R);

      PROTECT(z=AS_NUMERIC(z));
      pz2=NUMERIC_POINTER(z);
    }

    for(i=k-1; i<p-1; i++) // loop thats goes j=i+1,...,p-1
    {
        if(isFloat){
          a=pR1[p*i + i];
          b=pR1[p*i + i + 1];
        }else{
          a=pR2[p*i + i];
          b=pR2[p*i + i + 1];
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
             pR1[p*i + i]=c*a - s*b;
             pR1[p*i + i + 1]=s*a + c*b;
           }else{
             pR2[p*i + i]=c*a - s*b;
             pR2[p*i + i + 1]=s*a + c*b;
           }

           for(j=i+1; j<p-1; j++)  // loop thats goes j=i+1,...,p-1
           {
             if(isFloat){
               a=pR1[p*j + i];
               b=pR1[p*j + i + 1];
               pR1[p*j + i]=c*a - s*b;
               pR1[p*j + i + 1]=s*a + c*b;
             }else{
               a=pR2[p*j + i];
               b=pR2[p*j + i + 1];
               pR2[p*j + i]=c*a - s*b;
               pR2[p*j + i + 1]=s*a + c*b;
             }
           }
           for(j=0; j<nz; j++)  // loop thats goes j=1,...,nz
           {
              if(isFloat){
               a=pz1[p*j + i];
               b=pz1[p*j + i + 1];
               pz1[p*j + i]=c*a - s*b;
               pz1[p*j + i + 1]=s*a + c*b;
             }else{
               a=pz2[p*j + i];
               b=pz2[p*j + i + 1];
               pz2[p*j + i]=c*a - s*b;
               pz2[p*j + i + 1]=s*a + c*b;
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

SEXP writeBinFileFloat(SEXP filename, SEXP n, SEXP p, SEXP size, SEXP X,
   SEXP nRowNames, SEXP nColNames, SEXP rowNames, SEXP colNames, SEXP flagfloat)
{
    FILE *f=NULL;
    int i, j;
    int nrows, ncols, sizevar, sizeRowNames, sizeColNames, isFloat;
    //int inc=1;
    float *pX1;
    double *pX2;
    //double *linedouble;
    const char *s;
    char space[]="\0";
    float valuefloat;
    SEXP list;

    nrows=INTEGER_VALUE(n);
    ncols=INTEGER_VALUE(p);
    sizevar=INTEGER_VALUE(size);
    sizeRowNames=INTEGER_VALUE(nRowNames);
    sizeColNames=INTEGER_VALUE(nColNames);
    isFloat=asLogical(flagfloat);

    if(isFloat){
      PROTECT(X=AS_INTEGER(X));
      pX1=FLOAT(X);
    }else{
      PROTECT(X=AS_NUMERIC(X));
      pX2=NUMERIC_POINTER(X);
    }

    f=fopen(CHAR(STRING_ELT(filename,0)),"wb");
    fwrite(&nrows,4, 1 , f);
    fwrite(&ncols,4, 1 , f);
    fwrite(&sizevar,4, 1 , f);
    fwrite(&isFloat,4, 1 , f);
    fwrite(&sizeRowNames,4, 1 , f);
    fwrite(&sizeColNames,4, 1 , f);

    // Write rownames and colnames
    if(sizeRowNames > 0){
      for(i=0; i<nrows; i++){
         s=CHAR(STRING_ELT(rowNames, i));
         fwrite(s, 1, strlen(s), f);
         fwrite(&space, 1, 1, f);
      }
    }
    if(sizeColNames > 0){
      for(j=0; j<ncols; j++){
          s=CHAR(STRING_ELT(colNames, j));
          fwrite(s, 1, strlen(s), f);
          fwrite(&space, 1, 1, f);
      }
    }
    // Write lines
    for(i=0; i<nrows; i++)
    {
      for(j=0; j<ncols; j++)
      {
        if(sizevar==4)
        {
          if(isFloat){
            fwrite(pX1+nrows*j + i,sizevar, 1 , f);
          }else{
            valuefloat = pX2[nrows*j + i];
            fwrite(&valuefloat,sizevar, 1 , f);
          }
        }else{
          fwrite(pX2+nrows*j + i,sizevar, 1 , f);
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
    int i, j, k;
    int *psetRow, *psetCol;
    int nrows, ncols, sizevar, lsrow, lscol, sizeRowNames, sizeColNames, isFloat;
    int n, p, nr; // number of elements read
    float *pX1;
    float  *linefloat;
    double *pX2;
    double *linedouble;
    SEXP list, rowNames, colNames, X;

    lsrow=INTEGER_VALUE(nsetRow);
    lscol=INTEGER_VALUE(nsetCol);

    PROTECT(setRow=AS_INTEGER(setRow));
    psetRow=INTEGER_POINTER(setRow);

    PROTECT(setCol=AS_INTEGER(setCol));
    psetCol=INTEGER_POINTER(setCol);

    f=fopen(CHAR(STRING_ELT(filename,0)),"rb");
    nr=fread(&nrows, 4, 1, f);
    if(nr<1){
      Rprintf("    Warning: The function failed to read information on the number of rows\n");
    }
    nr=fread(&ncols, 4, 1, f);
    if(nr<1){
      Rprintf("    Warning: The function failed to read information on the number of columns\n");
    }
    nr=fread(&sizevar, 4, 1, f);
    if(nr<1){
      Rprintf("    Warning: The function failed to read information on the size in bytes\n");
    }
    nr=fread(&isFloat, 4, 1, f);
    if(nr<1){
      Rprintf("    Warning: The function failed to read information whether is 'float32'\n");
    }
    nr=fread(&sizeRowNames, 4, 1, f);
    if(nr<1){
      Rprintf("    Warning: The function failed to read information on rows names\n");
    }
    nr=fread(&sizeColNames, 4, 1, f);
    if(nr<1){
      Rprintf("    Warning: The function failed to read information on columns names\n");
    }

    n=lsrow > 0 ? lsrow : nrows;
    p=lscol > 0 ? lscol : ncols;

    if(isFloat || sizevar==4){
      linedouble=(double *) R_alloc(0, sizeof(double));   // will not be used
      linefloat=(float *) R_alloc(ncols, sizeof(float));

      X=PROTECT(allocMatrix(INTSXP, n, p));
      pX1=FLOAT(X);
    }else{
      linedouble=(double *) R_alloc(ncols, sizeof(double));
      linefloat=(float *) R_alloc(0, sizeof(float));   // will not be used

      X=PROTECT(allocMatrix(REALSXP, n, p));
      pX2=NUMERIC_POINTER(X);
    }

    PROTECT(rowNames = allocVector(STRSXP, nrows));
    char lineRowNames[sizeRowNames];

    PROTECT(colNames = allocVector(STRSXP, ncols));
    char lineColNames[sizeColNames];

    const char *s[100];

    // Read rownames and colnames
    if(sizeRowNames > 0){
      nr=fread(&lineRowNames,sizeRowNames,1,f);
      if(nr<1){
        Rprintf("    Warning: The function failed to read rows names of the file\n");
      }
      i=0;
      j=0;
      memset(s,'\0',sizeof(s));
      for(k=0; k<sizeRowNames; k++)
      {
        if(lineRowNames[k]=='\0'){
            SET_STRING_ELT(rowNames, i, mkChar(*s));
            j=0;
            i++;
            memset(s,'\0',sizeof(s));
        }else{
            s[j]=&lineRowNames[k];
            j++;
        }
      }
    }
    if(sizeColNames > 0){
      nr=fread(&lineColNames,sizeColNames,1,f);
      if(nr<1){
        Rprintf("    Warning: The function failed to read columns names of the file\n");
      }
      i=0;
      j=0;
      memset(s,'\0',sizeof(s));
      for(k=0; k<sizeColNames; k++)
      {
        if(lineColNames[k]=='\0'){
            SET_STRING_ELT(colNames, i, mkChar(*s));
            j=0;
            i++;
            memset(s,'\0',sizeof(s));
        }else{
            s[j]=&lineColNames[k];
            j++;
        }
      }
    }
    // Read lines
    for(i=0; i<n; i++)
    {
        if(lsrow > 0){
            // Move to the line indicated by setRow
            fseek(f, 24 + sizeRowNames + sizeColNames + ncols*sizevar*(psetRow[i]-1), SEEK_SET);
        }

        if(isFloat || sizevar==4){
            nr=fread(linefloat,sizevar,ncols,f);
        }else{
            nr=fread(linedouble,sizevar,ncols,f);
        }
        if(nr<ncols){
          Rprintf("    Warning: The function failed to read items at line %d of the file\n",i+1);
        }

        for(j=0; j<p; j++)
        {
          if(isFloat || sizevar==4){
            if(lscol>0)
            {
                pX1[n*j + i]=linefloat[psetCol[j]-1];
            }else{
                pX1[n*j + i]=linefloat[j];
            }
          }else{
              if(lscol>0)
              {
                  pX2[n*j + i]=linedouble[psetCol[j]-1];
              }else{
                  pX2[n*j + i]=linedouble[j];
              }
          }
        }
    }

    fclose(f);

    PROTECT(list = allocVector(VECSXP, 9));
    // Attaching outputs to list:
    SET_VECTOR_ELT(list, 0, ScalarInteger(n));
    SET_VECTOR_ELT(list, 1, ScalarInteger(p));
    SET_VECTOR_ELT(list, 2, ScalarInteger(sizevar));
    SET_VECTOR_ELT(list, 3, ScalarInteger(isFloat));
    SET_VECTOR_ELT(list, 4, ScalarInteger(sizeRowNames));
    SET_VECTOR_ELT(list, 5, ScalarInteger(sizeColNames));
    SET_VECTOR_ELT(list, 6, rowNames);
    SET_VECTOR_ELT(list, 7, colNames);
    SET_VECTOR_ELT(list, 8, X);

    UNPROTECT(6);
    return(list);
}
