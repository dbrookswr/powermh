\name{pwAnova}
\alias{pwAnova}
\title{
Estimates power for oneway ANOVA.}
\description{
Using simulation methods, the power for a oneway ANOVA with multiple comparisons/hypotheses is estimated. The function allows control over which hypotheses are considered and how many should be significant for the study to be successful. Further details are in the vignette.}
\usage{
pwAnova(dataset, replics, k = ifelse(is.matrix(dataset) || is.data.frame(dataset), length(unique(dataset[, 1])), length(dataset)), dcon = 0, pcon = 1, signcon = 0, numcon = k - 1, conmat = contr.treatment(k, 1), adjustcon = "none", radjust = FALSE, r2size = 0, r2p = 1, extrasuccess = NULL, rseed = FALSE, varyef = 0, varyn = 0, plotit = TRUE, varycolumn = NULL, dfv = 10, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataset}{
The data can be entered in one of two formats. First, it may be entered as a \texttt{list} of functions that define data for \texttt{k} conditions. The values \texttt{ef} and \texttt{n} can be used to vary any parameters used in constructing these functions, probably from a combination of \textsf{R}'s built-in random number generators (e.g., \texttt{rnorm}) and \texttt{sample}. The second format is as a $n \times 2$ matrix of data where the first column is the group variable and the second column is the response variable. There is no default.
}
  \item{replics}{
The number of replications to have. There is no default.
}
  \item{k}{
This is the number of conditions. The default is determined from the number of conditions defined by \texttt{dataset}. The value \texttt{k} is used within the function. The user may enter others values, but if a value other than the number of conditions in \texttt{dataset} is entered then there may be unintended consequences. For most purposes the user should not enter any value for this as the default is appropriate. 
}
  \item{dcon}{
A vector of length $k-1$ for the effect size of each contrast, defined as $\beta / sd$ where $sd$ is the residual standard error reported by \texttt{summary.lm}. If only a single value is used this is assumed for all contrasts. The default is 0, and 0 should be used for any contrast for which the effect size is not used for defining `success'. 
}
  \item{pcon}{
A vector of length $k-1$ for the $p$ value needed for each contrast. If only a single value is used this is assumed for all contrasts. The default is 1, and 1 should be used for any contrast for which the significance of the contrast is not used for defining `success'.
}
  \item{signcon}{
A vector of length $k-1$ for the sign of the $\beta$s. The elements of the vector should be -1, 0, or 1, to match the output of the \textsf{R} function \texttt{sign}. If a single value is given it is assumed to be the value for all contrasts. The default is 0, and 0 should be used for any contrast for which the direction of the effect is not used for defining `success'. In the rare case that $\hat\beta = 0$, this counts as a `success' for this contrast for the direction (though it would be non-significant and have a zero effect size so would not pass these hurdles if they are included).
}
  \item{numcon}{
The number of contrasts that need to meet criteria to have `success'. Defaults to $k-1$.
}
  \item{conmat}{
A contrast matrix of the appropriate dimension ($k$ by $k-1$) or reserved word (e.g., the \texttt{contr.treatment(3,2)}) appropriate for the \textsf{R} \texttt{contrasts} function's \texttt{contrasts} slot. The default for this is \texttt{contr.treatment(k-1,1)}.
}
  \item{adjustcon}{
The value entered into the \textsf{R} function \texttt{p.adjust} (see \texttt{p.adjust.methods}). The number of $p$ values used is the number of elements of \texttt{pcon} not equal to 1. The default is \texttt{"none"}, which corresponds to no adjustment.
}
  \item{radjust}{
Whether to use the adjusted or unadjusted $R^2$ produced by \textsf{R}'s \texttt{summary.lm}. Default is \texttt{FALSE} for unadjusted. 
\item{\texttt{r2size}} The critical value for the $R^2$ (adjusted or unadjusted), defaults to 0. If the \texttt{radjust = TRUE} and adjusted $R^2 < 0$ this will be below the default size for \texttt{r2size}. This means that if \texttt{r2size} has not been changed this will produce a `failure', which is appropriate.
}
  \item{r2p}{
The critical $p$ value for the overall fit of the model. Defaults to 1. 
}
  \item{extrasuccess}{
A list of functions appropriate for \textsf{R} that can be \texttt{TRUE} or \texttt{FALSE} with variables active when this part of the function is evaluated.  For example,\\
\texttt{extrasuccess <- list(function(...) \\
\hspace{1cm} pi < 4,function(...) mean(dd[,2])>1)} \\
will work because \texttt{dd} is the name of the \texttt{data.frame} used within the function and it is in the global environment when this function begins evaluating `success'. Thus, this option may require some knowledge of the internal workings of \texttt{pwAnova} (though for many procedures the user only needs to know about \texttt{dd}). Examples are included later in this paper. Defaults to \texttt{NULL}. 
}
  \item{rseed}{
This controls the seed to allow exact replication of the simulated results by using any integer (non-integers are truncated). The default (\texttt{FALSE}) is not to change the seed. The value \texttt{NULL} can also be used which causes the package to re-set the seed based on the time and ID process.
}
  \item{varyef}{
This can be entered as a matrix, a vector (or one-dimensional matrix), or a single value. If entered as a matrix it is a $r \times k$ matrix where $r$ is the number of values which are evaluated and $k$ is the number of conditions. If a vector is input with length greater than \texttt{k} it is assumed that its length is $r$ and the values are replicated across the \texttt{k} columns. If the length is equal to \texttt{k} then it is assumed $r = 1$ and each value of \texttt{varyef} corresponds to one of the \emph{k} conditions. If a single numeric value is given, this value is replicated to create a $1 \times k$ matrix.
}
  \item{varyn}{
This can be entered as a matrix, a vector (or one-dimensional matrix), or a single value. If entered as a matrix it is a $r \times k$ matrix where $r$ is the number of values which are evaluated and $k$ is the number of conditions. If a vector is input with length greater than \texttt{k} it is assumed that its length is $r$ and the values are divided by \texttt{k} (i.e., it is assumed the value is the total sample size) and replicated across the \texttt{k} columns. If the length is equal to \texttt{k} then it is assumed $r = 1$ and each value of \texttt{varyn} corresponds to one of the \emph{k} conditions. If a single numeric value is given, this value is divided by \texttt{k} replicated to create a $1 \times k$ matrix.
}
  \item{plotit}{
A binary (\texttt{TRUE}/\texttt{FALSE}) variable for whether to include a plot. The default is \texttt{TRUE}.
}
  \item{varycolumn}{
This tells the function which column to plot on the \emph{x-axis} of the plot. The default is to use the variable with the greatest variance for all the conditions. The numbers $1 \dotso k$ are for the effect sizes (\texttt{varyef}) for each of the \texttt{k} conditions and the numbers $k+1 \dotso 2k$ for the sample sizes (\texttt{varyn}) for each of the \texttt{k} conditions.
}
  \item{dfv}{
Controls the degree of freedom of the spline in the plot. The default is 10 if the number of unique values in \texttt{varycolumn} is at least 12, otherwise it defaults to 3 less than the number of unique values. The default is higher than will be appropriate for many uses.
}
  \item{\dots}{

}
}
\details{
See vignette.
}
\value{A matrix is returned, and if plotit = TRUE a plot is made. The matrix has the number of rows for the number of effects/sample sizes being examined. The first column is rep (the condition number), then the proportion of successes out of the number of replications. For the k groups, there are then k columns for effects and k columns for condition size.
}
\references{
}
\author{Daniel B. Wright}
\examples{
reps <- 1000
par(mfrow=c(2,2))
dd1ab <- cbind(gr <- sample(1:5,500,replace=TRUE),rnorm(500))
options(warn=-1)
eg1a <- pwAnova(dd1ab,replics=reps,varyn=seq(100,500,20),
                pcon=.05,dfv=3)
eg1b <- pwAnova(dd1ab,replics=reps,varyn=seq(100,500,20),
               pcon=.05,numcon=1,dfv=5)
dd1cd <- cbind(gr <- sample(1:5,1000,replace=TRUE),
               rnorm(1000)+.1*(gr>1)+.2*(gr>3))
eg1c <- pwAnova(dd1cd,replics=reps,varyn=seq(100,500,20),
               pcon=.05,dfv=3)
eg1d <- pwAnova(dd1cd,replics=reps,varyn=seq(100,500,20),
               numcon=1,pcon=.05,dfv=5)
abline(h=.8,lty=3)
}
\keyword{ ANOVA } 
\keyword{ Multiple Comparisons } 
\keyword{ Power } 
