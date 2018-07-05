pwAnova <-
function(dataset,replics,
   k = ifelse(is.matrix(dataset) || is.data.frame(dataset),
              length(unique(dataset[,1])),length(dataset)),
   dcon=0,pcon=1,signcon=0,numcon=k-1,conmat=contr.treatment(k,1),
   adjustcon="none",radjust=FALSE,r2size=0,r2p=1,extrasuccess=NULL,
   rseed=FALSE,varyef=0,varyn=0,plotit=TRUE,
   varycolumn=NULL,dfv = 10, ...)
  {
  if (rseed != FALSE) set.seed(rseed)
  if (plotit) {
     if (requireNamespace("splines", quietly = TRUE)==FALSE){ 
       plotit <- FALSE
       warning("you need splines to do plots")
     }}
  
  if (is.data.frame(dataset)) dataset <- as.matrix(dataset)
  if (is.matrix(dataset)==FALSE && is.list(dataset)==FALSE)
     stop("dataset not of proper type")
  if (is.matrix(dataset) && ncol(dataset) != 2)
     stop("your data set needs two columns")
  if (numcon > k-1) warning("It appears you are requiring more contrasts 
     to be successful than you have contrasts. Computation will proceed.
     Hope you know what you are doing.")
  if (numcon < 1) warning("It appears you are requiring less than one 
     contrast to be successful. There are other ways to performance this.
     Hope you know what you are doing.")
  vinput <- list(dcon,pcon,signcon)
  vnames <- c("dcon","pcon","signcon")
  for (q in seq_along(vinput))
    if (length(vinput[[q]]) != 1 && length(vinput[[q]]) != k-1) 
      stop(paste(vnames[q], "not proper length"))
  if (is.matrix(conmat) && (nrow(conmat) != k || ncol(conmat) != k-1)) 
    stop("conmat not the right dimensions")

  r <- max(NROW(varyn),NROW(varyef))
  if (r<=k) r <- 1
  if (NROW(varyef) > k && NROW(varyn) > k && plotit) {
        plotit <- FALSE
        warning("You are varying by both n and ef. The plotting is not
           designed for this and is off. Also, this function
           assumes these go equal numbers of trials. Be cautious.")}
  if (is.matrix(dataset) && (length(varyef) == k || NCOL(varyef) == k))
    warning("Looks like you are varying the effect size with the dataset
            as a data matrix. This may not produce what you want.
            Just the initial value in each row is used.")
  if ((length(varyn) == 1 && varyn%%k != 0) ||
      (is.vector(varyn) && length(varyn) > k && any(varyn%%k!=0)))
      warning("varyn not divible by k, values have be truncated.
             Enter as a matrix if you want don't want truncation.")
  if (length(varyn) == 1) n <- matrix(trunc(varyn/k),ncol=k,nrow=r)
  if (is.vector(varyn) && length(varyn) > k) 
    n <- replicate(k,trunc(varyn/3))
  if (length(varyn) == k) n <- t(replicate(r,varyn))
  if (is.matrix(varyn)) n <- varyn
  
  if (length(varyef) == 1) ef <- matrix(varyef,ncol=k,nrow=r)
  if (length(varyef) > k) ef <- replicate(k,varyef)
  if (length(varyef) == k) ef <- t(replicate(r,varyef))
  if (is.matrix(varyef)) ef <- varyef
  if (identical(dim(ef),dim(n))==FALSE) 
      stop("ef and n are not same dimensions. If you vary both
            they must have the same r values.")
  if (any(n < 3)) 
    stop("Some population group sizes too small. Minimum is set to 3.")
  pwOutput <- matrix(nrow=nrow(ef),ncol=2+2*k)
  for (i in 1:r){
    numsucc <- 0
    nv <- n[i,]
    efv <- ef[i,]
    for (j in 1:replics){
      if (is.matrix(dataset)==FALSE) {
         dv <- {}
         for (x in 1:k) dv <- c(dv,dataset[[x]](n=nv[x],ef=efv[x])) 
         dd <- as.data.frame(cbind(rep(1:k,nv),dv))}
      if (is.matrix(dataset)==TRUE && ncol(dataset)==2){ 
         dd1 <- rep(unique(dataset[,1]),nv)
         dd2 <- {}
         for (w in 1:k){ 
              subd <- dataset[dataset[,1]==unique(dataset[,1])[w],]
              dd2 <- c(dd2,
                subd[sample(1:nrow(subd),size=nv[w],replace=TRUE),][,2])}
         dd <- as.data.frame(cbind(dd1,dd2))
         }
      if (min(table(dd[,1])) < 3){
        warning("A group in a sample had a group n of 0, 1, or 2. 
          The trial was assigned failure. It is possible varyn is 
          not what you intended.")
        next()
      }
  
      
  assign("dd",dd,envir=.GlobalEnv)
  succ <- TRUE
  groupv <- as.factor(dd[,1])
  contrasts(groupv) <- conmat
  lmout <- lm(dd[,2]~groupv)
  if (summary(lmout)$r.square < r2size && radjust == FALSE) succ <- FALSE  
  if (summary(lmout)$adj.r.square < r2size && radjust == TRUE) succ<-FALSE
  if (pf(summary(lmout)$fstat[1],summary(lmout)$fstat[2],
    summary(lmout)$fstat[3],lower.tail=FALSE) > r2p) succ <- FALSE
  if (any(signcon != 0)) 
   {if (length(signcon)==1) signcon <- rep(signcon,k-1)
    signson <- sum(sign(coef(lmout)[2:k])*signcon >= 0)
    if (signson < numcon) succ <- FALSE}
  if (any(pcon!= 1)) 
   {if (length(pcon)==1) pcon <- rep(pcon,k-1)
    pon <- sum(summary(lmout)$coef[2:k,4] < pcon)
    if (pon < numcon) succ <- FALSE
        }
  if (any(dcon!= 0)) 
   {if (length(dcon)==1) dcon <- rep(dcon,k-1)
    coefstan <- summary(lmout)$coef[2:k,1]/summary(lmout)$sigma
    don <- sum(coefstan > dcon)
    if (don < numcon) succ <- FALSE}
    if (is.null(extrasuccess)==FALSE)
    for (p in seq_along(extrasuccess))
       if (extrasuccess[[p]]()==FALSE) succ <- FALSE
    numsucc <- numsucc +succ
    }
     pwOutput[i,] <- c(i,numsucc/replics,ef[i,],n[i,])
      }
    colnames(pwOutput) <- c("rep","propsucc",
          paste("ef",1:k,sep=""),paste("n",1:k,sep=""))

if (plotit && r > k) {
    if (is.null(varycolumn))
      varycolumn <- 2+ which.max(c(apply(ef,2,sd),c(apply(n,2,sd))))
    x <- pwOutput[,varycolumn]
    x2 <- seq(min(x),max(x),length.out=200)
    if (r < 12 && dfv > r-3) {
      dfv <- r-3
      print("The df for bs you entered has been changed")
         }
    if (dfv < 1) {plotit <- FALSE
       warning("r to small for plotting")}
    ml <- predict(glm(cbind(replics*pwOutput[,2],replics*(1-pwOutput[,2]))
                      ~bs(x,df=dfv),family=binomial),
                  type="response",newdata=data.frame(x=x2),se.fit=TRUE)
    plot(x,pwOutput[,2],ylim=c(max(min(ml$fit-2*ml$se,pwOutput[,2]),0),
                               min(max(ml$fit+2*ml$se,pwOutput[,2]),1)),
       xlab=paste("Variable in",colnames(pwOutput)[varycolumn]),
       ylab="Proportion success")
    polygon(c(x2,rev(x2)),
       c(ml$fit-2*ml$se.fit,rev(ml$fit+2*ml$se.fit)),
          border=NA,col="grey90")
whichmx <- which.max(ml$fit)
powmax <- ml$fit[whichmx]
lines(x,pwOutput[,2],col="grey65") 
lines(x2,ml$fit)
lines(rep(x2[whichmx],2),c(0,ml$fit[whichmx]),lty=3)
points(x,pwOutput[,2],col="grey65")
}
  return(pwOutput)  
}
