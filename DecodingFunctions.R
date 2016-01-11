

linModels <- function(df, dsname, summary=FALSE)
{
	suppressMessages(require(robust))
    # df is a data frame like df <- data.frame(V=avgR, nameLen, avgWordFreq, avgLetterFreq, RSR, type, salesrank, nR)
    # dsnames is a prefix to save the files 
    linModel <- lm(V ~ RSR, data=df, model=FALSE)
    linModelRobust <- lmrob(V ~ RSR, data=df)
  
    V.norm <- (df$V-mean(df$V))/sd(df$V)
    RSR.norm <- (df$RSR-mean(df$RSR))/sd(df$RSR)
    df.norm <- data.frame(V = V.norm, RSR = RSR.norm)
    linModelNorm <- lm(V ~ RSR, data = df.norm, model=FALSE)
  
    dir.create(dsname, showWarnings = FALSE)  
    save(linModel, file=paste(dsname,"LinModel.RData",sep="/"))
    save(linModelRobust, file=paste(dsname,"LinModelRobust.RData",sep="/"))
    save(linModelNorm, file=paste(dsname,"LinModelNorm.RData",sep="/"))

    if (summary)
    {

    	summaryLinModel <- summary(linModel)
   	 	save(summaryLinModel, file=paste(dsname,"SummaryLinModel.RData",sep="/"))

    	summaryLinModelNorm <- summary(linModelNorm)
	    save(summaryLinModelNorm, file=paste(dsname,"SummaryLinModelNorm.RData",sep="/"))
    }
    
}


rndLinModels <- function(df, dsname, nrep=10000, verbose=FALSE)
{
  rests <- NULL
  rintercepts <- NULL

  for (i in seq(1,nrep))
  {  
    if (verbose) { print(i) }
    rRSR <- sample(df$RSR) 
    rdata <- data.frame(V=df$V, RSR = rRSR)
    rmodel <- lm(V ~ RSR, data=rdata, model=FALSE)
    rests[i] <- rmodel$coefficients[2]
    rintercepts[i] <- rmodel$coefficients[1]
  }
  dir.create(dsname, showWarnings = FALSE)  
  save(rests, file=paste(dsname,"RndEsts.RData",sep="/"))
  save(rintercepts, file=paste(dsname,"RndIntercepts.RData",sep="/"))
  
}


bootLinModels <- function(df, dsname, nrep=10000, verbose=FALSE)
{
  bests <- NULL
  rhoests <- NULL
  
  for (i in seq(1,nrep))
  {  
    if (verbose) { print(i) }
    binds <- sample(length(df$RSR),replace=T)
    bRSR <- df$RSR[binds] 
    bV <- df$V[binds]
    bdata <- data.frame(V=bV, RSR = bRSR)
    rmodel <- lm(V ~ RSR, data=bdata, model=FALSE)
    bests[i] <- rmodel$coefficients[2]
    rhoests[i] <- cor.test(bdata$V,bdata$RSR, method="spearman")$estimate
  }
  dir.create(dsname, showWarnings = FALSE)  
  save(bests, file=paste(dsname, "Boot.RData", sep="/"))  
  save(rhoests, file=paste(dsname, "SpearmanBoot.RData", sep="/"))  
}


linModelsControl <- function(V, RSR, dfcontrols, dsname)
{
  # df is a data frame like df <- data.frame(V=avgR, nameLen, avgWordFreq, avgLetterFreq, RSR, type, salesrank, nR)
  # dsnames is a prefix to save the files 

  data <- cbind(RSR,cbind(V, dfcontrols))
  formula <- as.formula(paste("V ~ ", paste(names(dfcontrols), collapse= "+")))
  preModelControl <- lm(formula, data=data, model=FALSE)
  
  rV <- preModelControl$residuals
  extradata <- data.frame(rV, RSR)
  linModelControl <- lm(rV ~ RSR, data=extradata, model=FALSE)
  
  rV.norm <- (rV - mean(rV)) / sd(rV)
  RSR.norm <- (RSR - mean(RSR)) / sd(RSR)
  extradata.norm <- data.frame(rV=rV.norm, RSR=RSR.norm)
  linModelControlNorm <- lm(rV ~ RSR, data = extradata.norm, model=FALSE)
  
  dir.create(dsname, showWarnings = FALSE)  

  save(preModelControl, file=paste(dsname,"PreModelControl.RData",sep="/"))
  save(linModelControl, file=paste(dsname,"LinModelControl.RData",sep="/"))
  save(linModelControlNorm, file=paste(dsname,"LinModelControlNorm.RData",sep="/"))
}



plotLM <- function(plotData, rests, bests, xlab, ylab, booted=TRUE, coords=c(0,0))
{
  attach(plotData)
  par(mar=c(4,4,1,0.5))
  if (booted)
  { par(mfrow=c(1,2))}

  ylims <- c(min(plotData$lowerPred), max(plotData$upperPred))
  plot(plotData$RSRs, plotData$mnPred, xlim=range(plotData$RSRs), ylim=ylims, lwd=2, type="l", xlab="", ylab="", cex.axis=1.2)
  mtext(side=1,line=3,xlab, cex=1.7)
  mtext(side=2,line=2.5,ylab, cex=1.7)

  polygon(x=c(plotData$RSRs, rev(plotData$RSRs)), y=c(plotData$upperPred, rev(plotData$lowerPred)), col="lightgray", lty=0)
  lines(plotData$RSRs, plotData$mnPred, type="l", lwd=2)
	abline(v=coords[1], col="darkgray", lty=2)
	abline(h=coords[2], col="darkgray", lty=2)

  if (booted)
  {
    h1 <- hist(rests, plot=F, breaks=30)
    h2 <- hist(bests, plot=F, breaks=30)
  
    
	suppressMessages(require(Hmisc))
  
    xlims <- range(c(rests,bests))
	ylims <- c(max(c(h1$density,h2$density))*0.01,1.1*max(c(h1$density,h2$density)))
	plot(h1$mids, h1$density, col="blue", type="l", xlim=xlims,ylim=ylims, xlab="", ylab="", lwd=2, cex.axis=1.2) 
	polygon(x=c(0,h1$mids,0),  y = c(0,h1$density ,0), col=rgb(0,0,1,0.3), border="blue", lwd=2)
	polygon(x=c(0,h2$mids,0),  y = c(0,h2$density ,0), col=rgb(1,0,0,0.3), border="red", lwd=2)
	abline(h=0, lwd=2)
	abline(v=0, lty=1, col="black", lwd=1)
	abline(v=mean(bests), col="red", lwd=2, lty=2)
	abline(v=mean(rests), col="blue", lty=2, lwd=2)
	mtext(side=1, line=3, expression(hat(b)), cex=1.7)
	mtext(side=2, line=2.5, "Density", cex=1.7)
	box()
  }
}

statsTable <- function(dsname, linModel, bests, rests, rhoests, linModelRobust, msg=dsname)
{
  
  suppressMessages(require(xtable))
  sm <- summary(linModel)
  t <- sm$coefficients[2,3]
  tp <- sm$coefficients[2,4]

  est <- mean(bests)  
  bp <- 1-sum(bests>0)/length(bests)
  rp <- sum(rests>=est)/length(rests)  
  rhop <- 1-sum(rhoests>0)/length(rhoests)
  
  est.robust <- linModelRobust$coefficients[2]
  srobust <- summary(linModelRobust)
  p.robust <- srobust$coefficients[2,4]

  dftable <- data.frame(tp, bp, rp, rhop, p.robust)
  colnames(dftable) <- c("t test p-value", "bootstrap p-value", "permutation p-value", "Spearman p-value", "robust p-value")
  row.names(dftable) <- ""
  xtable(x=dftable, caption=paste("Tests for",msg), digits=10, type="html", comment=FALSE, table.placement="h!")
}



calcPlotDataSimple <- function(df, dsname)
{
  load(file=paste(dsname,"LinModel.RData",sep="/"))  
  w <- 0.01
  RSRs <- seq(min(df$RSR), max(df$RSR)+w, by=w)
  err <- stats::predict(linModel, newdata=data.frame(RSR=RSRs), interval = "confidence")
  preds <- err[,1]
  ucl <- err[,3]
  lcl <- err[,2]
  
  plotData <- data.frame(RSRs=RSRs, mnPred=preds, upperPred=ucl, lowerPred=lcl)
  
  dir.create(dsname, showWarnings = FALSE)  
  save(plotData, file=paste(dsname, "PlotDataSimple.RData", sep="/"))  
  
}


linModelsControl <- function(V, RSR, dfcontrols, dsname)
{
  # df is a data frame like df <- data.frame(V=avgR, nameLen, avgWordFreq, avgLetterFreq, RSA, type, salesrank, nR)
  # dsnames is a prefix to save the files 

  data <- cbind(RSR,cbind(V, dfcontrols))
  formula <- as.formula(paste("V ~ ", paste(names(dfcontrols), collapse= "+")))
  preModelControl <- lm(formula, data=data, model=FALSE)
  
  rV <- preModelControl$residuals
  extradata <- data.frame(rV, RSR)
  linModelControl <- lm(rV ~ RSR, data=extradata, model=FALSE)
  
  rV.norm <- (rV - mean(rV)) / sd(rV)
  RSR.norm <- (RSR - mean(RSR)) / sd(RSR)
  extradata.norm <- data.frame(rV=rV.norm, RSA=RSR.norm)
  linModelControlNorm <- lm(rV ~ RSR, data = extradata.norm, model=FALSE)
  
  dir.create(dsname, showWarnings = FALSE)  

  save(preModelControl, file=paste(dsname,"PreModelControl.RData",sep="/"))
  save(linModelControl, file=paste(dsname,"LinModelControl.RData",sep="/"))
  save(linModelControlNorm, file=paste(dsname,"LinModelControlNorm.RData",sep="/"))
}