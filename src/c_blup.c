#include "SFSI.h"
//#include "utils.c"

//====================================================================
// Recursive quantities for the GEMMA algorithm:
// a'Pib, a'PiPib, and a'PiPiPib, for any a and b
//====================================================================
double atPb(int i, int n, double *a, double *b, double *w, double *dbar){
  if(i == 0){  // atP0b - atP0w*btP0w/w1tP0w
    return(ddot3(n,a,dbar,b) - ddot3(n,a,dbar,w + n*i)*ddot3(n,b,dbar,w + n*i)/ddot3(n,w + n*i,dbar,w + n*i));

  }else{
    return(atPb(i-1,n,a,b,w,dbar) -
           atPb(i-1,n,a,w+n*i,w,dbar)*atPb(i-1,n,b,w+n*i,w,dbar)/atPb(i-1,n,w+n*i,w+n*i,w,dbar));
  }
}

//==========================================

double atPPb(int i, int n, double *a, double *b, double *w, double *dbar){
  double aPw, bPw, wPw;
  if(i == 0){  // atP0P0b + (atP0w)(btP0w)/(wtP0P0w)^2 + ...
    aPw = ddot3(n, a, dbar, w + n*i);
    bPw = ddot3(n, b, dbar, w + n*i);
    wPw = ddot3(n, w + n*i, dbar, w + n*i);
    return(ddot4(n,a,dbar,dbar,b) + aPw*bPw*ddot4(n,w + n*i,dbar,dbar,w + n*i)/pow(wPw,2) -
           aPw*ddot4(n,b,dbar,dbar,w + n*i)/wPw - bPw*ddot4(n,a,dbar,dbar,w + n*i)/wPw);

  }else{
    aPw = atPb(i-1,n,a,w+n*i,w,dbar);
    bPw = atPb(i-1,n,b,w+n*i,w,dbar);
    wPw = atPb(i-1,n,w+n*i,w+n*i,w,dbar);
    return(atPPb(i-1,n,a,b,w,dbar) + aPw*bPw*atPPb(i-1,n,w+n*i,w+n*i,w,dbar)/pow(wPw,2) -
           aPw*atPPb(i-1,n,b,w+n*i,w,dbar)/wPw - bPw*atPPb(i-1,n,a,w+n*i,w,dbar)/wPw);
  }
}

//==========================================

double atPPPb(int i, int n, double *a, double *b, double *w, double *dbar){
  double aPw, bPw, wPw, aPPw, bPPw, wPPw;
  if(i == 0){
    aPw = ddot3(n, a, dbar, w + n*i);
    bPw = ddot3(n, b, dbar, w + n*i);
    wPw = ddot3(n, w + n*i, dbar, w + n*i);
    aPPw = ddot4(n, a, dbar, dbar, w + n*i);
    bPPw = ddot4(n, b, dbar, dbar, w + n*i);
    wPPw = ddot4(n, w + n*i, dbar, dbar, w + n*i);
    return(ddot5(n,a,dbar,dbar,dbar,b) - aPw*bPw*pow(wPPw,2)/pow(wPw,3) -
           aPw*ddot5(n,b,dbar,dbar,dbar,w + n*i)/wPw - bPw*ddot5(n,a,dbar,dbar,dbar,w + n*i)/wPw -
           aPPw*bPPw/wPw + aPw*bPPw*wPPw/pow(wPw,2) + bPw*aPPw*wPPw/pow(wPw,2) +
           aPw*bPw*ddot5(n,w + n*i,dbar,dbar,dbar,w + n*i)/pow(wPw,2)
          );

  }else{
    aPw = atPb(i-1,n,a,w+n*i,w,dbar);
    bPw = atPb(i-1,n,b,w+n*i,w,dbar);
    wPw = atPb(i-1,n,w+n*i,w+n*i,w,dbar);
    aPPw = atPPb(i-1,n,a,w+n*i,w,dbar);
    bPPw = atPPb(i-1,n,b,w+n*i,w,dbar);
    wPPw = atPPb(i-1,n,w+n*i,w+n*i,w,dbar);
    return(atPPPb(i-1,n,a,b,w,dbar) - aPw*bPw*pow(wPPw,2)/pow(wPw,3) -
           aPw*atPPPb(i-1,n,b,w+n*i,w,dbar)/wPw - bPw*atPPPb(i-1,n,a,w+n*i,w,dbar)/wPw -
           aPPw*bPPw/wPw + aPw*bPPw*wPPw/pow(wPw,2) + bPw*aPPw*wPPw/pow(wPw,2) +
           aPw*bPw*atPPPb(i-1,n,w+n*i,w+n*i,w,dbar)/pow(wPw,2)
          );
  }
}

//==========================================

double tr_P(int i, int n, double *w, double *dbar){
  if(i == 0){  //  Tr_P0 - wtP0P0w/wP0w
    return(dsum(n,dbar) - ddot4(n,w + n*i,dbar,dbar,w + n*i)/ddot3(n,w + n*i,dbar,w + n*i));

  }else{
    return(tr_P(i-1,n,w,dbar) - atPPb(i-1,n,w+n*i,w+n*i,w,dbar)/atPb(i-1,n,w+n*i,w+n*i,w,dbar));
  }
}

//==========================================

double tr_PP(int i, int n, double *w, double *dbar){
  double wPw;
  if(i == 0){
    int inc1 = 1;
    wPw = ddot3(n, w + n*i, dbar, w + n*i);
    return(F77_NAME(ddot)(&n, dbar, &inc1, dbar, &inc1) +
           pow(ddot4(n,w + n*i,dbar,dbar,w + n*i),2)/pow(wPw,2) -
           2*ddot5(n,w + n*i,dbar,dbar,dbar,w + n*i)/wPw);

  }else{
    wPw = atPb(i-1,n,w+n*i,w+n*i,w,dbar);
    return(tr_PP(i-1,n,w,dbar) +
           pow(atPPb(i-1,n,w+n*i,w+n*i,w,dbar),2)/pow(wPw,2) -
           2*atPPPb(i-1,n,w+n*i,w+n*i,w,dbar)/wPw
          );
  }
}

//==========================================

double det_WtHinvW(int i, int n, double *w, double *dbar){
  // Obtained from Leibniz formula
  if(i == 0){  // det(W0'H^{-1}W0)*w1'P0w1 = w1'H^{-1}w1 = wPw
    return(ddot3(n, w + n*i, dbar, w + n*i)); // wPw

  }else{
    return(det_WtHinvW(i-1,n,w,dbar)*atPb(i-1,n,w+n*i,w+n*i,w,dbar));
  }
}

//====================================================================
// Log Likelihood (ML) and Restricted Log Likelihood (REML)
//====================================================================
double logLik(double ratio, int n, int p, double *Uty, double *UtX,
              double *d, double pi, double *dbar)
{
  int k;
  double logdetH = 0;   // log(det(H)) = log(prod(ratio*d+1))

  for(k=0; k<n; k++){
    dbar[k] = 1/(ratio*d[k] + 1);
    logdetH += log(ratio*d[k] + 1);
  }

  // Remove the factor 0.5: 0.5*n*log(0.5*n/pi) - 0.5*n - 0.5*logdetH - 0.5*n*log(ytPy)
  return(n*log(0.5*n/pi) - n - logdetH - n*log(atPb(p-1,n,Uty,Uty,UtX,dbar)));
}

//====================================================================
// Zhou & Stephen (2012) use degrees of freedom in REML equal to
// df = n-p-1 because they adjust single marker regression as well
// thus substracting 1 df. Here we do not adjust SMR so df = n-p

double logResLik(double ratio, int n, int p, double *Uty, double *UtX,
                 double *d, double pi, double *dbar)
{
  int k;
  double logdetH = 0;   // log(det(H)) = log(prod(ratio*d+1))

  for(k=0; k<n; k++){
    dbar[k] = 1/(ratio*d[k] + 1);
    logdetH += log(ratio*d[k] + 1);
  }

  // Remove the factor 0.5 and the term +log(detXtX)
  return((n-p)*log(0.5*(n-p)/pi) - (n-p) -
        logdetH - log(det_WtHinvW(p-1,n,UtX,dbar)) -
        (n-p)*log(atPb(p-1,n,Uty,Uty,UtX,dbar)));
}

//====================================================================
// First derivative of the Log Likelihood (ML)
//====================================================================
double dlogLik(double ratio, int n, int p, double *Uty, double *UtX,
               double *d, double pi, double *dbar)
{
  int k;
  for(k=0; k<n; k++){
    dbar[k] = 1/(ratio*d[k] + 1);
  }

  double Tr_HinvG = (n - dsum(n,dbar))/ratio;  // sum(diag(Hinv%*%G)) = (n-sumd)/ratio
  double ytPy = atPb(p-1, n, Uty, Uty, UtX, dbar);
  double ytPGPy = (ytPy - atPPb(p-1,n,Uty,Uty,UtX,dbar))/ratio;   // y'Px G Px y = (ytPy-ytPPy)/ratio

  // Remove the factor 0.5: -0.5*Tr_HinvG+0.5*n*ytPGPy/ytPy
  return(-Tr_HinvG + n*ytPGPy/ytPy);
}

//====================================================================
// First derivative of the Log-restricted Likelihood (REML)
//====================================================================
double dlogResLik(double ratio, int n, int p, double *Uty, double *UtX,
                  double *d, double pi, double *dbar)
{
  int k;

  for(k=0; k<n; k++){
    dbar[k] = 1/(ratio*d[k] + 1);
  }

  double Tr_PG = (n-p-tr_P(p-1,n,UtX,dbar))/ratio;  // sum(diag(Px%*%G))
  double ytPy = atPb(p-1, n, Uty, Uty, UtX, dbar);  // t(y)%*%Px%*%y
  double ytPGPy = (ytPy - atPPb(p-1,n,Uty,Uty,UtX,dbar))/ratio;  // y' Px G Px y = (ytPy-ytPPy)/ratio

  // Remove the factor 0.5: -0.5*Tr_PG+0.5*(n-p)*ytPGPy/ytPy
  return(-Tr_PG + (n-p)*ytPGPy/ytPy);
}

//====================================================================
// TRUE if x1*x2 negative
int RootBracketed(double x1,double x2) {
  int result;

  if((sign(x1)*sign(x2)) <= 0){
    result = 1;
  }else{
    result = 0;
  }

  return result;
}

// Define fun and dfun that will be later either REML or ML
double (*fun)(double, int, int, double*, double*, double*, double, double*) = NULL;
double (*dfun)(double, int, int, double*, double*,double*, double, double*) = NULL;

/*******************************************************
*              Brent Method Function                   *
* ---------------------------------------------------- *
* Reference:  BORLAND MATHEMATICAL LIBRARY             *
*                                                      *
*                C++ version by J-P Moreau, Paris.     *
*                       (www.jpmoreau.fr)              *
* ---------------------------------------------------- *
* The purpose is to find a real root of a real         *
* function f(x) using Brent's method.                  *
*                                                      *
* INPUTS:  x1,x2     : interval of root                *
*          Tolerance : desired accuracy for root       *
*          maxIter   : maximum number of iterations    *
*                                                      *
* OUTPUTS: The function returns the root value         *
*          ValueAtRoot : value of f(root)              *
*          niter    : number of done iterations        *
*          error    : =0, all OK                       *
*                   : =1, no root found in interval    *
*                   : =2, no more iterations !         *
*******************************************************/
double BrentRoots(double x1, double x2,
                  double Tolerance, int maxIterations,
                  double *valueAtRoot, int *niter, int *error,
                  int n, int p, double *Uty, double *UtX,
                  double *d, double pi, double *dbar){

  double FPP = DBL_EPSILON; // 1e-11;
  double nearzero = DBL_EPSILON/10; //1e-20;
  double result, AA, BB, CC, DD, EE, FA, FB, FC, Tol1, PP, QQ, RR, SS, xm;
  int i, done;

  i = 0;
  done = 0;
  *error = 0;
  result = NA_REAL;

  AA = x1;  BB = x2;
  FA = dfun(AA, n, p, Uty, UtX, d, pi, dbar);
  FB = dfun(BB, n, p, Uty, UtX, d, pi, dbar);
  if(RootBracketed(FA,FB)){
    FC = FB;
    do{
      if(!(RootBracketed(FC,FB))){
        CC = AA; FC = FA; DD = BB - AA; EE = DD;
      }
      if(fabs(FC) < fabs(FB)) {
        AA = BB; BB = CC; CC = AA;
        FA = FB; FB = FC; FC = FA;
      }
      Tol1 = 2.0 * FPP * fabs(BB) + 0.5 * Tolerance;
      xm = 0.5 * (CC-BB);
      if((fabs(xm) <= Tol1) || (fabs(FA) < nearzero)){
        result = BB;
        done = 1;
        *valueAtRoot = dfun(result, n, p, Uty, UtX, d, pi, dbar);
        //Rprintf("A root has been found\n");
      }else{
        if((fabs(EE) >= Tol1) && (fabs(FA) > fabs(FB))){
          SS = FB/ FA;
          if(fabs(AA - CC) < nearzero){
            PP = 2.0 * xm * SS;
            QQ = 1.0 - SS;
          }else{
            QQ = FA/FC;
            RR = FB /FC;
            PP = SS * (2.0 * xm * QQ * (QQ - RR) - (BB-AA) * (RR - 1.0));
            QQ = (QQ - 1.0) * (RR - 1.0) * (SS - 1.0);
          }
          if(PP > nearzero) QQ = -QQ;
          PP = fabs(PP);
          if((2.0 * PP) < fmin2(3.0*xm *QQ-fabs(Tol1 * QQ), fabs(EE * QQ))){
            EE = DD;  DD = PP/QQ;
          }else{
            DD = xm;   EE = DD;
          }
        }else{
          DD = xm;
          EE = DD;
        }
        AA = BB;
        FA = FB;
        if(fabs(DD) > Tol1){
          BB = BB + DD;
        }else{
          if(xm > 0){
             BB = BB + fabs(Tol1);
          }else{
             BB = BB - fabs(Tol1);
          }
        }
        FB = dfun(BB, n, p, Uty, UtX, d, pi, dbar);
        i++;
      }
	  }
    while((!done) && (i < maxIterations));
    if(i >= maxIterations) *error = 2;
  }else{
    //Rprintf("Values f(x) at bounds (%f,%f) are NOT of opposite signs\n",x1,x2);
    *error = 1;
  }
  *niter = i;
  return result;
}

SEXP R_solve_mixed(SEXP N_, SEXP ratio_, SEXP trn_, SEXP ytrn_,
            SEXP X_, SEXP Z_, SEXP K_, SEXP U_,
            SEXP d_, SEXP bounds_, SEXP tol_,
            SEXP maxiter_, SEXP dmin_, SEXP isREML_,
            SEXP isEigen_, SEXP BLUE_, SEXP BLUP_)
{
    int j;
    int *trn;
    double x1, x2, value, ratio, varE, varU, h2;
    double *ytrn, *X, *d, *U, *work, *dbar, *tmp1, *tmp2, *tmp3;
    double *Uty, *UtX, *yHat;
    double *bounds;
    double one = 1;
    int inc1 = 1;
    int intval;
    SEXP list, uHat_, yHat_, bHat_;
    int nprotect = 8;

    double tol=NUMERIC_VALUE(tol_);
    double dmin=NUMERIC_VALUE(dmin_);
    int maxiter=INTEGER_VALUE(maxiter_);
    int N=INTEGER_VALUE(N_);
    int ntrn=Rf_length(trn_);
    int p=Rf_ncols(X_);
    int isREML=asLogical(isREML_);
    int isEigen=asLogical(isEigen_);
    int BLUE=asLogical(BLUE_);
    int BLUP=asLogical(BLUP_);
    int nintervals=Rf_length(bounds_)-1;

    PROTECT(trn_=AS_INTEGER(trn_));
    trn=INTEGER_POINTER(trn_);

    PROTECT(ytrn_=AS_NUMERIC(ytrn_));
    ytrn=NUMERIC_POINTER(ytrn_);

    PROTECT(X_=AS_NUMERIC(X_));
    X=NUMERIC_POINTER(X_);

    PROTECT(d_=AS_NUMERIC(d_));
    d=NUMERIC_POINTER(d_);

    PROTECT(U_=AS_NUMERIC(U_));
    U=NUMERIC_POINTER(U_);

    PROTECT(bounds_=AS_NUMERIC(bounds_));
    bounds=NUMERIC_POINTER(bounds_);

    double dfx0 = NA_REAL;
    int convergence = NA_INTEGER;
    int status = NA_INTEGER;
    int ndsmall = NA_INTEGER;
    int niter = NA_INTEGER;
    int error = NA_INTEGER;
    double pi = M_PI, eps = DBL_EPSILON;

    Uty = (double *) R_alloc(ntrn, sizeof(double));
    UtX = (double *) R_alloc(ntrn*p, sizeof(double));

    yHat_ = PROTECT(Rf_allocVector(REALSXP, N));
    yHat = NUMERIC_POINTER(yHat_);

    // Check for small eigenvalues
    ndsmall=0;
    for(j=0; j<ntrn; j++){
      if(d[j] < dmin){
        ndsmall++;
        d[j] = 0;
      }
    }

    intval = N>(p*p) ? N : p*p;
    dbar=(double *) R_alloc(ntrn, sizeof(double));
    work=(double *) R_alloc(intval, sizeof(double));
    tmp1=(double *) R_alloc(intval, sizeof(double));
    tmp2=(double *) R_alloc(intval, sizeof(double));
    tmp3=(double *) R_alloc(intval, sizeof(double));

    //Rprintf(" Calculating UtX and Uty...\n");
    for(j=0; j<p; j++){
      slice_matrix(N,X,tmp1,ntrn,trn,j,2); // X[trn,j]
      matrix_vector_product(ntrn,ntrn,&one,U,tmp1,1,UtX+ntrn*j,1); // UtX[,j]=U'X[trn,j]
    }
    matrix_vector_product(ntrn,ntrn,&one,U,ytrn,1,Uty,1); // U'ytrn

    if(Rf_isNull(ratio_)){ // Perform likelihood maximization
      // Define function for REML or ML
      if(isREML){
        fun = &logResLik;
        dfun = &dlogResLik;
      }else{
        fun = &logLik;
        dfun = &dlogLik;
      }

      double fx1, fx2, fx0, slope, dfx1, dfx2;
      error = 1;
      status = 0;
      convergence = 0;

      // Evaluate the likelihood at the bounds and get slope
      x1 = bounds[0];
      x2 = bounds[nintervals];
      fx1 = fun(x1, ntrn, p, Uty, UtX, d, pi, dbar);
      fx2 = fun(x2, ntrn, p, Uty, UtX, d, pi, dbar);
      slope = (fx2-fx1)/(x2-x1);

      //Rprintf(" Starting searching of optimal value ...\n");
      if(fabs(slope) < eps){
        ratio = x1;
        convergence = 0;
        status = 1;
      }else{
        //Rprintf("Starting Brent's algorithm (within interval)...\n");
        for(j=0; j<nintervals; j++){
          dfx1 = dfun(bounds[j], ntrn, p, Uty, UtX, d, pi, dbar);
          dfx2 = dfun(bounds[j+1], ntrn, p, Uty, UtX, d, pi, dbar);

          if(sign(dfx1)*sign(dfx2) <= 0){
            value = BrentRoots(bounds[j], bounds[j+1], tol, maxiter,
                             &dfx0, &niter, &error,
                             ntrn, p, Uty, UtX, d, pi, dbar);
            if(error != 1){
              break;
            }
          }
        }

        if(error == 1){
          //Rprintf("No solution was found in the interval\n");
          //Rprintf("The coordinate bound yielding the larger likelihood is returned\n");
          if(fx1 > fx2){
            ratio = x1;
            status = 3;
          }else{
            ratio = x2;
            status = 4;
          }
          convergence = 0;
        }else{
          ratio = value;
          fx0 = fun(ratio, ntrn, p, Uty, UtX, d, pi, dbar);
          if(fx0 < fmax2(fx1, fx2)){  // Degenerated solution
            //Rprintf("A degenerated solution was found: x=%f and f(x)=%f\n",ratio,fx0);
            if(fx1 > fx2){
              ratio = x1;
              status = 3;
            }else{
              ratio = x2;
              status = 4;
            }
            convergence = 0;
          }else{
            //Rprintf("A good solution was found: x=%f and f(x)=%f\n",ratio,fx0);
            convergence = (error ==  0 ? 1 : 0);
            if(!convergence){
              status = 2;
            }
          }
        }
      }
    }else{
      //Rprintf("No likelihood optimization was performed\n");
      ratio = NUMERIC_VALUE(ratio_);
    }

    // Variance components and fixed effects
    for(j=0; j<ntrn; j++){
      dbar[j] = 1/(ratio*d[j] + 1);
      tmp1[j] = Uty[j]*dbar[j];
    }

    if(BLUE){
      double *bHat;

      bHat_ = PROTECT(Rf_allocVector(REALSXP, p));
      bHat = NUMERIC_POINTER(bHat_);
      nprotect++;

      crossproduct(ntrn,1,p,tmp1,UtX,tmp2);   // tmp2=t(Uty*dbar)%*%UtX
      crossproduct_scale(ntrn,p,p,UtX,dbar,UtX,tmp1,work); // tmp1 = UtX' dbar UtX
      invert_matrix(p,tmp1,tmp3,&eps,work);  // tmp3=solve(UtX'*dbar*UtX)
      matrix_vector_product(p,p,&one,tmp3,tmp2,1,bHat,0); // b = solve(UtX'*dbar*UtX)*tmp2
      matrix_vector_product(N,p,&one,X,bHat,1,yHat,0);    // yHat = Xb

      for(j=0; j<ntrn; j++){  // ytrn* = ytrn - X[trn,]b
        ytrn[j]-=yHat[trn[j]];
      }
    }else{
      bHat_ = R_NilValue;
      memset(yHat, 0, N*sizeof(double));
    }

    varE = atPb(p-1,ntrn,Uty,Uty,UtX,dbar)/(isREML ? ntrn-p : ntrn);
    varU = ratio*varE;
    h2 = varU/(varU + varE);

    if(BLUP){
      double *Z, *K, *M, *uHat;
      int nu = Rf_isNull(Z_) ? N : Rf_ncols(Z_);

      uHat_ = PROTECT(Rf_allocVector(REALSXP, nu));
      uHat = NUMERIC_POINTER(uHat_);
      nprotect++;

      M = (double *) R_alloc(N*(ntrn>nu ? ntrn : nu), sizeof(double));

      if(!Rf_isNull(Z_)){
        PROTECT(Z_=AS_NUMERIC(Z_));
        Z = NUMERIC_POINTER(Z_);
        nprotect++;
      }
      if(!Rf_isNull(K_)){
        PROTECT(K_=AS_NUMERIC(K_));
        K = NUMERIC_POINTER(K_);
        nprotect++;
      }

      memset(M, 0, ntrn*ntrn*sizeof(double));
      if(isEigen){
        // Matrix B = KZ'Hinv = U diag{d*ratio/(d*ratio+1)} U'
        for(j=0; j<ntrn; j++){
          tmp1[j] = d[j]*ratio/(ratio*d[j] + 1);
        }
        tcrossproduct_scale(ntrn,ntrn,ntrn,U,tmp1,U,M,work);
        matrix_vector_product(ntrn,ntrn,&one,M,ytrn,1,uHat,0);
      }else{
        //  Hinv = U diag{1/(theta+d)} U' = U diag{ratio/(ratio*d+1)} = UDU'
        for(j=0; j<ntrn; j++){
          tmp1[j] = ratio/(ratio*d[j] + 1);
        }

        tcrossproduct_scale(ntrn,ntrn,ntrn,U,tmp1,U,M,work); // Hinv = UDU'
        matrix_vector_product(ntrn,ntrn,&one,M,ytrn,1,tmp1,0); // tmp1=Hinv%*%(y-Xb)

        if(Rf_isNull(Z_) && Rf_isNull(K_)){
          memset(uHat, 0, nu*sizeof(double));
          for(j=0; j<ntrn; j++){
            uHat[trn[j]] = tmp1[j];
          }
        }else{
          if(Rf_isNull(Z_)){   // Z=NULL, K=K:   u = K[,trn]*Hinv*(y-Xb)
            matrix_vector_product_subset(N,N,K,tmp1,uHat,0,NULL,ntrn,trn,0,work);
          }else{
            if(Rf_isNull(K_)){  // Z=Z, K=NULL:   u = Z[trn,]'Hinv*(y-Xb)
              matrix_vector_product_subset(N,nu,Z,tmp1,uHat,ntrn,trn,0,NULL,1,work);

            }else{  // u = KZ[trn,]'Hinv*(y-Xb)
              matrix_vector_product_subset(N,nu,Z,tmp1,tmp2,ntrn,trn,0,NULL,1,work);
              matrix_vector_product(nu,nu,&one,K,tmp2,1,uHat,0);
            }
          }
        }
      }

      if(Rf_isNull(Z_)){
        F77_NAME(daxpy)(&N,&one,uHat,&inc1,yHat,&inc1);   // yHat = Xb + u
      }else{
        matrix_vector_product(N,nu,&one,Z,uHat,1,tmp1,0); // Zu
        F77_NAME(daxpy)(&N,&one,tmp1,&inc1,yHat,&inc1);   // yHat = Xb + Zu
      }
    }else{
      uHat_ = R_NilValue;
    }

    PROTECT(list = Rf_allocVector(VECSXP, 12));
    SET_VECTOR_ELT(list, 0, ScalarReal(ratio));
    SET_VECTOR_ELT(list, 1, ScalarReal(dfx0));
    SET_VECTOR_ELT(list, 2, ScalarInteger(niter));
    SET_VECTOR_ELT(list, 3, ScalarInteger(status));
    SET_VECTOR_ELT(list, 4, ScalarInteger(convergence));
    SET_VECTOR_ELT(list, 5, ScalarInteger(ndsmall));
    SET_VECTOR_ELT(list, 6, ScalarReal(varU));
    SET_VECTOR_ELT(list, 7, ScalarReal(varE));
    SET_VECTOR_ELT(list, 8, ScalarReal(h2));
    SET_VECTOR_ELT(list, 9, bHat_);
    SET_VECTOR_ELT(list, 10, yHat_);
    SET_VECTOR_ELT(list, 11, uHat_);

    UNPROTECT(nprotect);

    return(list);
}
