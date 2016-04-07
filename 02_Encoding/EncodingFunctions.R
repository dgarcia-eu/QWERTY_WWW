
rndLinModels <- function(df, dsname, nrep=10000, verbose=FALSE, robust=FALSE)
{
  if (robust)
  { require(robust) }
  rests1 <- NULL
  rests2 <- NULL

  for (i in seq(1,nrep))
  {  
    if (verbose) { print(i) }
    ord <- sample(length(df$R)) 
    rdata <- data.frame(V=df$V, R = df$R[ord], L=df$L[ord])
    if (! robust)
      { rmodel <- lm(V~R*L, data=rdata, model=FALSE) }
    if (robust)
      { rmodel <- lmrob(V~R*L, data=rdata, model=FALSE) }
    rests1[i] <- rmodel$coefficients[2]
   	rests2[i] <- rmodel$coefficients[3]
  }
  dir.create(dsname, showWarnings = FALSE)  
  save(rests1, file=paste(dsname,"RndEsts1.RData",sep="/"))
  save(rests2, file=paste(dsname,"RndEsts2.RData",sep="/"))
  
}


bootLinModels <- function(df, dsname, nrep=10000, verbose=FALSE, robust=FALSE)
{
  if (robust)
  { require(robust) }

  bests1 <- NULL
  bests2 <- NULL
  
  for (i in seq(1,nrep))
  {  
    if (verbose) { print(i) }
    binds <- sample(length(df$R),replace=T)
    bR <- df$R[binds] 
    bL <- df$L[binds] 
    bV <- df$V[binds]
    bdata <- data.frame(V=bV, R = bR, L=bL)
    if (! robust)
      { rmodel <- lm(V~R*L, data=bdata, model=FALSE) }
    if (robust)
      { rmodel <- lmrob(V~R*L, data=bdata, model=FALSE) }
    bests1[i] <- rmodel$coefficients[2]
    bests2[i] <- rmodel$coefficients[3]
  }
  dir.create(dsname, showWarnings = FALSE)  
  save(bests1, file=paste(dsname, "Boot1.RData", sep="/"))  
  save(bests2, file=paste(dsname, "Boot2.RData", sep="/"))  
}

