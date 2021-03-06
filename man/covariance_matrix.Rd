
\name{covariance_matrix}
\alias{cov2dist}
\alias{cov2cor2}
\title{Conversion of a covariance matrix to a distance or correlation matrix}
\usage{
cov2dist(V, void = FALSE)

cov2cor2(V, a = 1, void = FALSE)
}
\arguments{
\item{V}{(numeric matrix) Symmetric variance-covariance matrix among \eqn{p}{p} variables. It can be of the "float32" type as per the 'float' R-package}

\item{void}{\code{TRUE} or \code{FALSE} to whether return or not return the output.
When \code{FALSE} no result is displayed but the input \code{V} is modified. Default \code{void=FALSE}}

\item{a}{(numeric) A number to multiply the whole resulting matrix by. Default \code{a=1}}
}
\value{
Function 'cov2dist' returns a matrix containing the (square) Euclidean distances. Function 'cov2cor2' returns a correlation matrix
}
\description{
Computes a correlation matrix or a Euclidean distance matrix from a covariance matrix
}
\details{
For any variables
\ifelse{html}{\out{X<sub>i</sub>}}{\eqn{X_i}{X_i}} and
\ifelse{html}{\out{X<sub>j</sub>}}{\eqn{X_j}{X_j}} with mean zero and with sample vectors
\ifelse{html}{\out{<b>x</b><sub>i</sub> = (x<sub>i1</sub>,...,x<sub>in</sub>)'}}{\eqn{\textbf{x}_i=(x_{i1},...,x_{in})'}{x_i=(x_i1,...,x_in)'}} and
\ifelse{html}{\out{<b>x</b><sub>j</sub> = (x<sub>j1</sub>,...,x<sub>jn</sub>)'}}{\eqn{\textbf{x}_j=(x_{j1},...,x_{jn})'}{x_j=(x_j1,...,x_jn)'}}
, their (sample) variances are equal (up-to a constant) to their cross-products, this is,
\ifelse{html}{\out{var(X<sub>i</sub>) = <b>x</b>'<sub>i</sub><b>x</b><sub>i</sub>}}{\eqn{var(X_i)=\textbf{x}_i'\textbf{x}_i}{var(X_i) = x'_i x_i}} and
\ifelse{html}{\out{var(X<sub>j</sub>) = <b>x</b>'<sub>j</sub><b>x</b><sub>j</sub>}}{\eqn{var(X_j)=\textbf{x}_j'\textbf{x}_j}{var(X_j) = x'_j x_j}}.
Likewise, the covariance is
\ifelse{html}{\out{cov(X<sub>i</sub>,X<sub>j</sub>) = <b>x</b>'<sub>i</sub><b>x</b><sub>j</sub>}}{\eqn{cov(X_i,X_j)=\textbf{x}_i'\textbf{x}_j}{cov(X_i,X_j) = x'_i x_j}}.

\ifelse{html}{\out{<b>Distance.</b>}}{\eqn{\textbf{Distance.}}{Distance.}}
The square of the distance
\ifelse{html}{\out{d(X<sub>i</sub>,X<sub>j</sub>)}}{\eqn{d(X_i,X_j)}{d(X_i,X_j)}}
between the variables expressed in terms of cross-products is

\ifelse{html}{\out{<p style='text-align:center'>d<sup>2</sup>(X<sub>i</sub>,X<sub>j</sub>) = <b>x</b>'<sub>i</sub><b>x</b><sub>i</sub> + <b>x</b>'<sub>j</sub><b>x</b><sub>j</sub> - 2<b>x</b>'<sub>i</sub><b>x</b><sub>j</sub></p>}}{\deqn{d^2(X_i,X_j) = \textbf{x}_i'\textbf{x}_i + \textbf{x}_j'\textbf{x}_j - 2\textbf{x}_i'\textbf{x}_j}{d^2(X_i,X_j) = x'_i x_i + x'_j x_j - 2x'_i x_j}}

Therefore, the output (square) distance matrix will contain as off-diagonal entries

\ifelse{html}{\out{<p style='text-align:center'>d<sup>2</sup>(X<sub>i</sub>,X<sub>j</sub>) = var(X<sub>i</sub>) + var(X<sub>j</sub>) - 2cov(X<sub>i</sub>,X<sub>j</sub>)</p>}}{\deqn{d^2(X_i,X_j) = var(X_i)+var(X_j)- 2cov(X_i,X_j)}{d^2(X_i,X_j) = var(X_i) + var(X_j) - 2cov(X_i,X_j)}}

while in the diagonal, the distance between one variable with itself is
\ifelse{html}{\out{d<sup>2</sup>(X<sub>i</sub>,X<sub>i</sub>) = 0}}{\eqn{d^2(X_i,X_i)=0}{d^2(X_i,X_i) = 0}}

\ifelse{html}{\out{<b>Correlation.</b>}}{\eqn{\textbf{Correlation.}}{Correlation.}}
The correlation between the variables is obtained from variances and covariances as

\ifelse{html}{\out{<p style='text-align:center'>cor(X<sub>i</sub>,X<sub>j</sub>) = cov(X<sub>i</sub>,X<sub>j</sub>)/(sd(X<sub>i</sub>)sd(X<sub>j</sub>))</p>}}{\deqn{cor(X_i,X_j) = cov(X_i,X_j)/(sd(X_i)sd(X_j))}{cor(X_i,X_j) = cov(X_i,X_j)/(sd(X_i)sd(X_j))}}

where \ifelse{html}{\out{sd(X<sub>i</sub>)=sqrt(var(X<sub>i</sub>))}}{\eqn{sd(X_i)=\sqrt{var(X_i)}}{sd(X_i) = sqrt(var(X_i))}}; while in the diagonal, the correlation between one variable with itself is
\ifelse{html}{\out{cor(X<sub>i</sub>,X<sub>i</sub>) = 1}}{\eqn{cor(X_i,X_i)=1}{cor(X_i,X_i) = 1}}

Variances are obtained from the diagonal values while covariances are obtained from the out-diagonal.
}
\examples{
  require(SFSI)
  data(wheatHTP)
  
  X = scale(Y[,4:7])
  (V = crossprod(X))             # Covariance matrix 
  
  # Covariance matrix to distance matrix
  (D1 = cov2dist(V))
  # it must equal (but faster) to:
  D0 = as.matrix(dist(t(X)))^2
  max(abs(D0-D1))
  
  # Covariance to a correlation matrix
  (R1 = cov2cor2(V))
  # it must equal (but faster) to:
  R0 = cov2cor(V)
  max(abs(R0-R1))
  
  if(requireNamespace("float")){
   # Using a 'float' type variable
   V2 = float::fl(V)
   D2 = cov2dist(V2)
   max(abs(D1-D2))   # discrepancy with previous matrix
   R2 = cov2cor2(V2)
   max(abs(R1-R2))   # discrepancy with previous matrix
  }
  
  # Using void=TRUE
  cov2dist(V,void=TRUE)
  V       # notice that V was modified
  cov2dist(V2,void=TRUE)
  V2       # notice that V2 was modified
  
}
\author{
Marco Lopez-Cruz (\email{maraloc@gmail.com}) and Gustavo de los Campos
}
