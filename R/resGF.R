# In general, try to group together related functions into the same .R file
# https://tinyheero.github.io/jekyll/update/2015/07/26/making-your-first-R-package.html

# tag above your function to indicate this function to be “exposed” to users to use.


# check this pacakge - https://github.com/cran/randomForest


#' Generate resistance surface from GF analysis
#'
#' This function generate a resistance surface based on gradient forest analysis.
#'
#' @param obj A GF object
#' @param raster_stack A raster stack oject containing the same variable used in the GF analysis
#' @param save.image Save individual conversion of landscape variables into resistance surfaces
#' @param results_dir Directory to save the transformation images
#' @return A resistance surface
#' @export
resGF <- function(obj,
                  raster_stack,
                  save.image = TRUE,
                  results_dir) {
  #GF_R_Total <- GF_total_importance(obj) # simplify this?
  mylength <- raster::ncell(raster_stack[[1]])
  myvector <- vector(mode="numeric", length=mylength)
  Env_pc <- 0
  s <- raster::stack()
  for (i in names(raster_stack)) {
    importance.df <- obj$res[obj$res$var==i,] # me
    mycumu <- getCU(importance.df, importance.randomForest(obj, "Weighted")[i], i, obj)
    CU <- cumimp(obj, i)
    ymax <- max(CU$y)
    imp <- importance(obj)[i]
    resA <- obj$res[obj$res$var == i, ]
    splits <- resA$split
    w <- pmax(resA$improve.norm, 0)
    X <- na.omit(obj$X[,i])
    rX <- range(X)
    dX <- diff(rX)

    # importance density I(x)
    dImp <- density(splits, weight = w/sum(w), from = rX[1], to = rX[2])
    if((dX/dImp$bw) > 50)
      dImp <- density(splits, weight = w/sum(w), from = rX[1], to = rX[2], bw=dX/50)
    dImpNorm <- tryCatch( { normalize.density(dImp,imp,integrate=T)}, error = function(e){dImpNorm <- normalize.histogram(dImp,imp)} )

    # data density d(x)
    dObs <- density(X, from = rX[1], to = rX[2])
    if((dX/dObs$bw) > 50)
      dObs <- density(X, from = rX[1], to = rX[2], bw=dX/50)
    dObs <- whiten(dObs,lambda=0.9)
    # dObsNorm <- normalize.density(dObs,imp,integrate=T)
    dObsNorm <-tryCatch( { normalize.density(dObs,imp,integrate=T)}, error = function(e){dObsNorm <- normalize.histogram(dObs,imp)} )

    # standardized density f(x) = I(x)/d(x)
    dStd <- dImp
    dStd$y <- dImp$y/dObs$y
    # dStdNorm <- tryCatch( { normalize.density(dStd,imp,integrate=T)}, error = function(e){dStdNorm <- dStd} )
    dStdNorm <- try(normalize.density(dStd,imp,integrate=T), silent=TRUE) # this sometimes does not converge, added the silence
    if (class(dStdNorm) == "try-error")
      dStdNorm <- normalize.histogram(dStd,imp)

    # getting the y axis normalised to get the individual contribution of each variables
    dStd2 <- dStd
    dStd2$y <- dStd2$y * importance.randomForest(obj, "Weighted")[i] / max(dStd2$y)
    dImpNorm2 <- dImpNorm
    dImpNorm2$y <- dImpNorm$y * importance.randomForest(obj, "Weighted")[i] / max(dImpNorm$y)
    dObsNorm2 <- dObsNorm
    dObsNorm2$y <- dObsNorm$y * importance.randomForest(obj, "Weighted")[i] / max(dObsNorm$y)
    dStdNorm2 <- dStdNorm
    dStdNorm2$y <- dStdNorm$y * importance.randomForest(obj, "Weighted")[i] / max(dStdNorm$y)

    # dinv <- dStdNorm
    dinv <- dStdNorm2
    dinv$call <- "fp(x)"
    r <- raster_stack[[i]]
    #print(r)
    v0 <- raster::values(r)
    v1 <- as.matrix(v0)
    df1 <- as.data.frame(v0)
    df1$ID <- seq.int(nrow(df1))
    v1b <- v1[!rowSums(!is.finite(v1)),]
    uniq <- unique(v1b)

    # Transofomation
    v2 <- approx(dinv$x, dinv$y, uniq, rule=2)$y
    df2 <- data.frame(uniq,v2)
    names(df2)[names(df2) == "uniq"] <- "v0"
    myderivDF5 <- merge(df1, df2, by='v0', all = TRUE)
    myderivDF5 <- myderivDF5[order(myderivDF5$ID),]
    # Getting individual raster
    r_GF_tempmin  <- raster::raster()
    names(r_GF_tempmin)=i
    raster::extent(r_GF_tempmin) <- raster::extent(r)
    dim(r_GF_tempmin) <- dim(r)
    raster::crs(r_GF_tempmin) <- raster_stack
    # myderivDF5$v2[myderivDF5$v2<0] <- 0 # change $y into $v2
    # stdvalues <- myderivDF5$v2 * gf$overall.imp[v] / GF_R_Total # method 1
    stdvalues <- myderivDF5$v2  # method 2
    # stdvalues <- myderivDF5$v2 * imp / sum(extendedForest::importance(obj)) # method 2
    # stdvalues <- myderivDF5$v2 * extendedForest::importance(obj, "Weighted") / extendedForest::importance(obj, "Weighted")[i]

    myvector <- myvector + stdvalues # change $y into $v2
    raster::values(r_GF_tempmin) <- myderivDF5$v2 # change $y into $v2 - non-standardized values

    variable_name <- paste0(i, ".pdf")
    #Plotting
    bin=F
    nbin=101
    leg.panel=1
    barwidth=1
    leg.posn="topright"
    cex.legend=0.8
    line.ylab=1.5
    ci <- cumimp(obj,i,standardize=T, standardize_after=T)
    ci$y <- diff(c(0,ci$y))
    nice.names <- i
    # ci <- tryCatch( {normalize.histogram(ci, imp, bin=bin | !is.binned(obj), nbin=nbin)}, error = function(e){ci <- ci} )
    ci <-normalize.histogram(ci, imp, bin=bin | !is.binned(obj), nbin=nbin)

    ci$y <- ci$y * importance.randomForest(obj, "Weighted")[i] / max(ci$y)


    if (save.image) {
      pdf(paste0(results_dir, variable_name)) # to get it save to file
      par(mfrow=c(2,3), oma=c(0,0,2,0))
      # plot(mycumu$x, mycumu$y, main= "Cummulative values")
      # lines(mycumu, col=2)
      plot(cumimp(obj,i,standardize=FALSE),type='n',main="",ylab="Cum. importance",xlab=i)
      # lines(gradientForest::cumimp(gf,i),type='l',col="black")
      # lines(gradientforest::cumimp(gf,i,standardize_after=TRUE),type='l',col="blue")
      lines(cumimp(obj,i,standardize=FALSE),type='l',col="black")
      # legend(par("usr")[1],par("usr")[4],legend=c("Standardize then Normalize (default)",
      #                                             "Normalize then Standardize","Normalize only"),col=c("black","blue","red"),lty=1,cex=0.8)
      raster::plot(r, main= "original raster")
      hist(v1, breaks=50, col = 'grey', xlab = i, main = "intial values" ) # value of the initial raster ; main = paste("intial values for" , v)
      plot(ci, type='h', col="grey60", xlim=range(splits), lwd=barwidth,
           ylim = c(0, max(dStdNorm2$y)*1.1), lend=2, xlab = i, ylab = "density")
      lines(dStdNorm2, col = "blue", lwd = 2)
      abline(h = mean(dStdNorm2$y)/mean(dStd2$y), lty = 2, col = "blue")
      # legend(leg.posn, legend = c("Density of splits", "Density of data", "Ratio of densities", "Ratio=1"), lty = c(1, 1, 1, 2), col = c("black", "red", "blue", "blue"), cex = cex.legend, bty = "n", lwd = 1)
      raster::plot(r_GF_tempmin, main= "resistance surface")
      hist(myderivDF5$v2, breaks=50, col = 'grey', xlab = i,  main = "transformed values") # value of the initial raster
      mtext(i, line=0, side=3, outer=TRUE, cex=1.2)
      dev.off() # to get it save to file

    }
    # percent = round(gf$overall.imp[v] / GF_R_Total * 100, 2)
    percent = round(imp / sum(importance.randomForest(obj)) * 100, 2)
    Env_pc <- Env_pc + percent
    par(mfrow=c(2,3), oma=c(0,0,2,0))
    # plot(mycumu$x, mycumu$y, main= "Cummulative values")
    # lines(mycumu, col=2)
    plot(cumimp(obj,i,standardize=FALSE),type='n',main="",ylab="Cum. importance",xlab=i)
    # lines(gradientForest::cumimp(gf,i),type='l',col="black")
    # lines(gradientForest::cumimp(gf,i,standardize_after=TRUE),type='l',col="blue")
    lines(cumimp(obj,i,standardize=FALSE),type='l',col="black")
    raster::plot(r, main= "original raster")
    hist(v1, breaks=50, col = 'grey', xlab = i, main = "intial values" ) # value of the initial raster ; main = paste("intial values for" , v)
    plot(ci, type='h', col="grey60", xlim=range(splits), lwd=barwidth,
         ylim = c(0, max(dStdNorm2$y)*1.1), lend=2, xlab = i, ylab = "density")
    lines(dStdNorm2, col = "blue", lwd = 2)
    abline(h = mean(dStdNorm2$y)/mean(dStd2$y), lty = 2, col = "blue")
    raster::plot(r_GF_tempmin, main= "resistance surface")
    hist(myderivDF5$v2, breaks=50, col = 'grey', xlab = i,  main = "transformed values") # value of the initial raster
    mtext(paste0(i,': ',percent,'%'), line=0, side=3, outer=TRUE, cex=1.2)
    s <- raster::stack(s, r_GF_tempmin)
    print(paste0(i,': ',percent,'%'))
  }
  # creating a single resistant raster
  par(mfrow=c(1,1))
  single_r  <- raster::raster()
  names(single_r)='final_resistance'
  raster::extent(single_r) <- raster::extent(r)
  dim(single_r) <- dim(r)
  raster::values(single_r) <- myvector
  raster::crs(single_r) <- raster::crs(raster_stack)
  single_r <- climateStability::rescale0to1(single_r)
  # raster::plot(single_r)
  print(paste0('overall predicators: ',Env_pc,'%'))
  print(paste0('positive snp: ', obj$species.pos.rsq))
  return(single_r)
}



#' Return gradient forest importance
#'
#' This function generate a resistance surface based on gradient forest analysis.
#'
#' @param obj A GF object
#' @param raster_stack A raster stack oject containing the same variable used in the GF analysis
#' @return importance of the different variables
#' @export
get_importance <- function(obj, raster_stack) {
  Env_pc <- 0
  myimp <- NULL;
  for (i in names(raster_stack)) {
    importance.df <- obj$res[obj$res$var==i,] # me
    # https://r-forge.r-project.org/scm/viewvc.php/pkg/gradientForest/R/whiten.R?view=markup&revision=2&root=gradientforest&pathrev=16
    imp <- importance.randomForest(obj)[i]
    percent = round(imp / sum(importance.randomForest(obj)) * 100, 2)
    print(paste0(i,': ',percent,'%'))
    myobj <- c(i, percent)
    myimp <- rbind(myimp, myobj)
    Env_pc <- Env_pc + percent
  }
  print(paste0('overall predicators: ',Env_pc,'%'))
  myobj <- c('overall', Env_pc)
  myimp <- rbind(myimp, myobj)
  rownames(myimp) <- myimp[,1]
  myimp <- as.data.frame(myimp)
  myimp[,1] <- NULL
  return(myimp)
}

######################
### Other function
######################
lower <- function(matrix) {
  if (is.vector(matrix) == TRUE ||
      dim(matrix)[1] != dim(matrix)[2]) {
    warning("Must provide square distance matrix with no column or row names")
  }
  lm <- matrix[lower.tri(matrix)]
  return(lm)
}

getCU <- function(importance.df, Rsq, i, gf) {  # ME - made using both standardization
  agg <- with(importance.df, agg.sum(improve.norm, list(split), sort.it=TRUE))
  cum.split <- agg[,1]
  height <- agg[,2]
  dinv <- normalize(inverse(gf$dens[[i]])) # crucial to normalize this case
  dinv.vals <- approx(dinv$x, dinv$y, cum.split, rule=2)$y
  #par(mfrow=c(1,2), oma=c(0,0,2,0))
  #plot(dinv)
  #plot(dinv.vals)
  par(mfrow = c(1, 1))
  if (any(bad <- is.na(height))) {
    cum.split <- cum.split[!bad]
    height <- height[!bad]
    dinv.vals <- dinv.vals[!bad]
  }
  # height <- height/sum(height)*Rsq # if (standardize & !standardize_after
  height <- height * dinv.vals
  res <- list(x=cum.split, y=cumsum(height))
}

# aggregate
agg.sum <- function(x, by, sort.it = F)
{
  if(!is.data.frame(by))
    by <- data.frame(by)
  if(sort.it) {
    ord <- do.call("order", unname(by))
    x <- x[ord]
    by <- by[ord,  , drop=F]
  }
  logical.diff <- function(group)
    group[-1] != group[ - length(group)]
  change <- logical.diff(by[[1]])
  for(i in seq(along = by)[-1])
    change <- change | logical.diff(by[[i]])
  by <- by[c(T, change),  , drop = F]
  by$x <- diff(c(0, cumsum(x)[c(change, T)]))
  by
}

normalize <- function(f) {
  integral <- integrate(approxfun(f,rule=2),lower=min(f$x),upper=max(f$x),  stop.on.error = FALSE)$value
  f$y <- f$y/integral*diff(range(f$x));
  f
}

inverse <- function(dens) {dens$y <- 1/dens$y; dens}

normalize.histogram <- function(ci,integral=1,bin=F,nbin=101) {
  # scale y values so that histogram integrates to integral
  # optionally aggregate the y's into binned x ranges
  if (bin) {
    brks <- seq(min(ci$x),max(ci$x),len=nbin)
    xx <- cut(ci$x,breaks=brks,inc=T)
    yy <- tapply(ci$y,xx,sum)
    yy[is.na(yy)] <- 0
    ci <- list(x=0.5*(brks[-1]+brks[-nbin]),y=yy)
  }
  dx <- min(diff(ci$x));
  Id <- sum(ci$y*dx);
  ci$y <- ci$y/Id*integral;
  ci
}

normalize.density <- function(d,integral=1,integrate=T) {
  # scale y values so that density integrates to integral
  Id <- if(integrate) integrate.density(d) else 1;
  d$y <- d$y/Id*integral;
  d
}

integrate.density <- function(d) {
  integrate(approxfun(d,rule=2),lower=min(d$x),upper=max(d$x))$value
}




scale.density <- function(d,scale=1/mean(d$y)) {
  d$y <- d$y*scale
  d
}

whiten <-
  function(dens, lambda=0.9)
  {
    # add a small uniform value to the density to avoid zeroes when taking inverse
    dens$y <- lambda*dens$y + (1-lambda)/diff(range(dens$x))
    dens
  }

is.binned <- function(obj) {
  compact <- obj$call$compact
  if (is.null(compact))
    FALSE
  else
    eval(compact)
}


gf_importance <- function(obj) {
  Env_pc <- 0
  myimp <- NULL;
  variable <- NULL;
  for (i in names(obj$overall.imp) ) {
    imp <- importance.randomForest(obj)[i]
    percent = round(imp / sum(importance.randomForest(obj)) * 100, 2)
    #print(paste0(i,': ',percent,'%'))
    myimp <- rbind(myimp, percent)
    variable <- rbind(variable, i)
    Env_pc = Env_pc + imp
  }
  #print(paste0('overall predicators: ',Env_pc,'%'))
  myimp <- rbind(myimp, Env_pc)
  variable <- rbind(variable, 'overall predicators')
  rownames(myimp) <- variable
  colnames(myimp) <- "imp"
  return(myimp)
}



#' Create gradientforest objects
#'
#' gradientForest uses an extended version of the package randomForest (Liaw and Wiener 2002), extendedForest which retains all of the split values and fit improvements
#' for each of the response variables (species catches in our case) for further analysis. gradientForest collates the numerous split values along each gradient and their `
#' associated fit improvements for each species that were retained by extendedForest, for each predictor in each tree and each forest. Details on the method are given in
#' Ellis et al. (2012) and applications are described in Pitcher et al. (2012).
#' https://rdrr.io/rforge/gradientForest/man/gradientForest.html
#'
#' @param data data.frame containing where rows identify sites and columns contain response variables (usually species catch (numbers or weight) or predictor variables such as physical or chemical, variables. Column names identify species or specific predictor variable. If the species are numeric variables, a regression forest is calculated. If the species are factor variables, a classification forest is calculated.
#' @param predictor.vars vector identifying which columns containing predictor variables (e.g., physical variables) are to be used in the randomForest analysis. This vector can contain column names (as a character) or column number.
#' @param response.vars vector identifying which species are to be used in the randomForest analysis. This vector can contain column names (as a character) or column number.
#' @param ntree number of bootstrapped trees to be generated by randomForest. Default set to 10 trees.
#' @param mtry number of predictor variables randomly sampled as candidates at each split. Setting to NULL accepts default values. Note that the default values are different for classification (sqrt(p) where p is number of variables in x) and regression (p/3).
#' @param transform a function defining a transformation to be applied the species data. For example, a square-root transformation would be entered as transform=function(x){sqrt(x)}. Default set to no transformation.
#' @param maxLevel if maxLevel == 0, compute importance from marginal permutation distribution of each variable (the default). If maxLevel > 0, compute importance from conditional permutation distribution of each variable, permuted within 2^{maxLevel} partitions of correlated variables.
#' @param corr.threshold if maxLevel > 0, OOB permuting is conditioned on partitions of variables having absolute correlation > corr.threshold.
#' @param compact logical variable to choose standard method or compact method for aggregating importance measures across species. Compact=TRUE to be chosen when memory problems cause a crash in this function. Still experimental.
#' @param nbin number of bins for compact option. Default set to 101.
#' @param trace if TRUE show the progress. Default FALSE.
#' @param check.names if TRUE then ensure that all predictor and response vars are syntactically valid col names in R and throw an error if they are not. gradientForest should still work with invalid col names, but with potentially more bugs. Default TRUE.
#' @return A resistance surface
#' @export
gradientForest <- function (data, predictor.vars, response.vars, ntree = 10, mtry = NULL,
          transform = NULL, maxLevel = 0, corr.threshold = 0.5, compact = FALSE,
          nbin = 101, trace = FALSE) {
  if (!inherits(data, "data.frame"))
    stop("'data' must be a data.frame")
  X <- data[predictor.vars]
  Y <- data[response.vars]
  if (compact) {
    bins <- do.call("cbind", lapply(X, function(x) bin(x,
                                                       nbin = nbin)))
  }
  if (!is.null(transform)) {
    Y <- apply(Y, 2, transform)
  }
  imp <- matrix(0, 0, 2, dimnames = list(NULL, c("%IncMSE",
                                                 "IncNodePurity")))
  if (is.null(mtry))
    fitcmd <- quote(randomForest::randomForest(Species ~ rhs, data = cbind(Y,
                                                             X), maxLevel = maxLevel, keep.forest = TRUE, importance = TRUE,
                                 ntree = ntree, keep.group = TRUE, keep.inbag = TRUE,
                                 corr.threshold = corr.threshold, na.action = na.omit))
  else fitcmd <- quote(randomForest::randomForest(Species ~ rhs, data = cbind(Y,
                                                                X), maxLevel = maxLevel, keep.forest = TRUE, importance = TRUE,
                                    ntree = ntree, mtry = mtry, keep.group = TRUE, keep.inbag = TRUE,
                                    corr.threshold = corr.threshold, na.action = na.omit))
  result <- list()
  species.pos.rsq <- 0
  form.rhs <- as.formula(paste("Y ~ ", paste(predictor.vars,
                                             collapse = "+")))
  if (trace) {
    spcount <- 0
    cat("Calculating forests for", length(response.vars),
        "species\n")
  }
  for (spec in response.vars) {
    if (trace)
      cat(if ((spcount <- spcount + 1)%%options("width")$width ==
              0)
        "\n."
        else ".")
    try({
      thisfitcmd <- do.call("substitute", list(fitcmd,
                                               list(Species = as.name(spec), SpeciesName = spec,
                                                    ntree = ntree, rhs = form.rhs[[3]])))
      fit <- eval(thisfitcmd)
      if (fit$type == "regression") {
        if (!is.na(fit$rsq[fit$ntree])) {
          if (fit$rsq[fit$ntree] > 0) {
            species.pos.rsq <- species.pos.rsq + 1
            if (compact) {
              result[[spec]] <- getSplitImproveCompact(fit,
                                                       bins)
            }
            else {
              result[[spec]] <- getSplitImprove(fit,
                                                X)
            }
            imp <- rbind(imp, fit$importance)
          }
        }
      }
      else if (fit$type == "classification") {
        if (!is.na(fit$err.rate[fit$ntree, "OOB"])) {
          p <- sum(Y[[spec]] == levels(Y[[spec]])[1])/length(Y[[spec]])
          err0 <- 2 * p * (1 - p)
          if (fit$err.rate[fit$ntree, "OOB"] < 2 * p *
              (1 - p)) {
            species.pos.rsq <- species.pos.rsq + 1
            if (compact) {
              result[[spec]] <- getSplitImproveClassCompact(fit,
                                                            bins, err0)
            }
            else {
              result[[spec]] <- getSplitImproveClass(fit,
                                                     X, err0)
            }
            nclass <- length(levels(Y[[spec]]))
            imp <- rbind(imp, fit$importance[, -(1:nclass)])
          }
        }
      }
      else stop(paste("unknown randomForest type:", fit$type))
    }, silent = FALSE)
  }
  if (!length(result)) {
    warning("No species models provided a positive R^2. \nThe gradient forest is empty")
    return(NULL)
  }
  rsq <- sapply(result, function(x) x$rsq[1])
  imp.rsq <- matrix(imp[, 1], length(predictor.vars), dimnames = list(predictor.vars,
                                                                      names(result)))
  imp.rsq[imp.rsq < 0] <- 0
  imp.rsq <- sweep(imp.rsq, 2, colSums(imp.rsq, na.rm = T),
                   "/")
  imp.rsq <- sweep(imp.rsq, 2, rsq, "*")
  overall.imp <- tapply(imp[, 1], dimnames(imp)[[1]], mean,
                        na.rm = T)
  overall.imp2 <- tapply(imp[, 2], dimnames(imp)[[1]], mean,
                         na.rm = T)
  out1 <- list(X = X, Y = Y, result = result, overall.imp = overall.imp,
               overall.imp2 = overall.imp2, ntree = ntree, imp.rsq = imp.rsq,
               species.pos.rsq = species.pos.rsq, ranForest.type = fit$type)
  out2 <- Impurity.based.measures(out1)
  out1$result <- rsq
  out <- c(out1, out2, call = match.call())
  class(out) <- c("gradientForest", "list")
  out
}




importance <- function(x, ...)  UseMethod("importance")

importance.default <- function(x, ...)
  stop("No method implemented for this class of object")


#' Extract variable importance measure
#' This is the extractor function for variable importance measures as produced by randomForest.
#'
#' Here are the definitions of the variable importance measures. For each tree, the prediction accuracy on the out-of-bag portion of the data is recorded. Then the same is done after permuting each predictor variable. The difference between the two accuracies are then averaged over all trees, and normalized by the standard error. For regression, the MSE is computed on the out-of-bag data for each tree, and then the same computed after permuting a variable. The differences are averaged and normalized by the standard error. If the standard error is equal to 0 for a variable, the division is not done (but the measure is almost always equal to 0 in that case).
#'
#'The second measure is the total decrease in node impurities from splitting on the variable, averaged over all trees. For classification, the node impurity is measured by the Gini index. For regression, it is measured by residual sum of squares.
#'
#'
#' @param x an object of class randomForest
#' @param type either 1 or 2, specifying the type of importance measure (1=mean decrease in accuracy, 2=mean decrease in node impurity).
#' @param class for classification problem, which class-specific measure to return.
#' @return A (named) vector of importance measure, one for each predictor variable.
#' @examples
#' set.seed(4543)
#' data(mtcars)
#' mtcars.rf <- randomForest(mpg ~ ., data=mtcars, ntree=1000, keep.forest=FALSE, importance=TRUE)
#' importance(mtcars.rf)
#' importance(mtcars.rf, type=1)
#' @export
importance.randomForest <- function(x, type=NULL, class=NULL, scale=TRUE,
                                    ...) {
  if (!inherits(x, "randomForest"))
    stop("x is not of class randomForest")
  classRF <- x$type != "regression"
  hasImp <- !is.null(dim(x$importance)) || ncol(x$importance) == 1
  hasType <- !is.null(type)
  if (hasType && type == 1 && !hasImp)
    stop("That measure has not been computed")
  allImp <- is.null(type) && hasImp
  if (hasType) {
    if (!(type %in% 1:2)) stop("Wrong type specified")
    if (type == 2 && !is.null(class))
      stop("No class-specific measure for that type")
  }

  imp <- x$importance
  if (hasType && type == 2) {
    if (hasImp) imp <- imp[, ncol(imp), drop=FALSE]
  } else {
    if (scale) {
      SD <- x$importanceSD
      imp[, -ncol(imp)] <-
        imp[, -ncol(imp), drop=FALSE] /
        ifelse(SD < .Machine$double.eps, 1, SD)
    }
    if (!allImp) {
      if (is.null(class)) {
        ## The average decrease in accuracy measure:
        imp <- imp[, ncol(imp) - 1, drop=FALSE]
      } else {
        whichCol <- if (classRF) match(class, colnames(imp)) else 1
        if (is.na(whichCol)) stop(paste("Class", class, "not found."))
        imp <- imp[, whichCol, drop=FALSE]
      }
    }
  }
  imp
}



getSplitImproveCompact <- function(fit, bins) {
  #   Return a data-frame: var name, rsq, split value, improvement
  #   Compact the splits into bins defined by bins matrix
  #   The i'th bin for predictor p is the interval (bin[i,p],bin[i+1,p])
  #   Every predictor is split into the same number of bins (nrow(bins)-1)

  #   extract all trees to a matrix and select for splits with some improvement
  trees <- lapply(1:fit$ntree, function(k) try(getTree(fit, k),silent=TRUE)) #Nick Ellis 10/12/2009
  ok <- sapply(trees, class) != "try-error"
  tmp <- do.call("rbind", lapply((1:fit$ntree)[ok], function(k) cbind(tree = k, trees[[k]])))
  tmp <- tmp[tmp[,"status"]==-3 & zapsmall(tmp[,"improve"]) > 0,c("split var","split point","improve")]
  colnames(tmp) <- c("var_n","split","improve")
  rownames(tmp) <- NULL

  #   assign the split to the appropriate bin and aggregate importance in each bin
  Xnames <- colnames(bins)
  tmp <- data.frame(var=Xnames[tmp[,"var_n"]], tmp, bin=rep(0,nrow(tmp)))
  for(p in Xnames) {
    sub <- with(tmp,var==p)
    tmp$bin[sub] <- as.numeric(cut(tmp$split[sub], bins[,p], include=TRUE, ordered=TRUE))
  }
  tmp <- with(tmp[tmp$bin>0,],agg.sum(improve,list(var,bin),sort.it=TRUE))
  names(tmp) <- c("var","bin","improve")

  #   Set the split value to the bin centre, but retain the bin number in case
  #   the bin centre is not appropriate value
  tmp <- cbind(tmp,split=rep(NA,nrow(tmp)),rsq=fit$rsq[fit$ntree])
  midpoints <- function(x) 0.5*(x[-1]+x[-length(x)]) # points between equally spaced points
  for(p in Xnames) {
    sub <- with(tmp,var==p)
    tmp$split[sub] <- midpoints(bins[,p])[tmp$bin[sub]]
  }
  tmp[,c("var","rsq","split","improve","bin")]
}


getSplitImprove <-function(fit, X) {
  #   return a data-frame: var name, rsq, var number, split value, improvement
  trees <- lapply(1:fit$ntree, function(k) try(getTree(fit, k),silent=TRUE)) #Nick Ellis 10/12/2009
  ok <- sapply(trees, class) != "try-error"
  tmp <- do.call("rbind", lapply((1:fit$ntree)[ok], function(k) cbind(tree = k, trees[[k]])))
  tmp <- tmp[tmp[,"status"]==-3 & zapsmall(tmp[,"improve"]) > 0,c("split var","split point","improve")]
  colnames(tmp) <- c("var_n","split","improve")
  rownames(tmp)<-NULL     #S.J. Smith 11/05/2009
  res <- cbind(data.frame(var=names(X)[tmp[,"var_n"]],rsq=rep(fit$rsq[fit$ntree],nrow(tmp))),tmp)
  ok <- zapsmall(res[,"improve"]) > 0
  res[ok,]
}


getSplitImproveClassCompact <- function(fit, bins, err0) {
  #   Return a data-frame: var name, rsq, split value, improvement
  #   Compact the splits into bins defined by bins matrix
  #   The i'th bin for predictor p is the interval (bin[i,p],bin[i+1,p])
  #   Every predictor is split into the same number of bins (nrow(bins)-1)

  #   extract all trees to a matrix and select for splits with some improvement
  trees <- lapply(1:fit$ntree, function(k) try(getTree(fit, k),silent=TRUE)) #Nick Ellis 10/12/2009
  ok <- sapply(trees, class) != "try-error"
  tmp <- do.call("rbind", lapply((1:fit$ntree)[ok], function(k) cbind(tree = k, trees[[k]])))
  tmp <- tmp[tmp[,"status"]== 1 & zapsmall(tmp[,"improve"]) > 0,c("split var","split point","improve")]
  colnames(tmp) <- c("var_n","split","improve")
  rownames(tmp) <- NULL

  #   assign the split to the appropriate bin and aggregate importance in each bin
  Xnames <- colnames(bins)
  tmp <- data.frame(var=Xnames[tmp[,"var_n"]], tmp, bin=rep(0,nrow(tmp)))
  for(p in Xnames) {
    sub <- with(tmp,var==p)
    tmp$bin[sub] <- as.numeric(cut(tmp$split[sub], bins[,p], include=TRUE, ordered=TRUE))
  }
  tmp <- with(tmp[tmp$bin>0,],agg.sum(improve,list(var,bin),sort.it=TRUE))
  names(tmp) <- c("var","bin","improve")

  #   Set the split value to the bin centre, but retain the bin number in case
  #   the bin centre is not appropriate value
  tmp <- cbind(tmp,split=rep(NA,nrow(tmp)),rsq=rep((err0-fit$err.rate[fit$ntree, "OOB"])/err0, nrow(tmp)))
  for(p in Xnames) {
    sub <- with(tmp,var==p)
    tmp$split[sub] <- midpoints(bins[,p])[tmp$bin[sub]]
  }
  tmp[,c("var","rsq","split","improve","bin")]
}

getSplitImproveClass <-
  function(fit, X, err0)
  {
    #   return a data-frame: var name, rsq, var number, split value, improvement
    trees <- lapply(1:fit$ntree, function(k) try(getTree(fit, k),silent=TRUE)) #Nick Ellis 10/12/2009
    ok <- sapply(trees, class) != "try-error"
    tmp <- do.call("rbind", lapply((1:fit$ntree)[ok], function(k) cbind(tree = k, trees[[k]])))
    tmp <- tmp[tmp[,"status"]==1,c("split var","split point","improve")]
    dimnames(tmp) <- list(NULL,c("var_n","split","improve"))
    res<-cbind(data.frame(var=names(X)[tmp[,"var_n"]],rsq=rep((err0-fit$err.rate[fit$ntree,"OOB"])/err0,nrow(tmp))),tmp)
    res
}


Impurity.based.measures <-function(obj)
{
  #becomes an internal function not usually used by users
  #Modified 07/10/2009 by SJS re: NE changes for classification trees.
  dens <- lapply(names(obj$X), function(i) density(na.omit(obj$X[,i]),from=min(na.omit(obj$X[,i])),to=max(na.omit(obj$X[,i]))))
  dens <- lapply(dens,whiten,lambda=0.90) # hard-coded whitening
  names(dens) <- names(obj$X)
  res <- do.call("rbind", lapply(names(obj$result), function(spec) cbind(spec=spec,obj$result[[spec]]))) #added by Smith 13/05/2009
  res$improve <- pmax(0,res$improve)
  res$rsq <- pmax(0,res$rsq)   #added by Ellis 12/05/2009
  res$improve.tot <- tapply(res$improve,res$spec,sum)[res$spec]
  res$improve.tot.var <- tapply(res$improve,interaction(res$spec,res$var),sum)[interaction(res$spec,res$var)]
  res$improve.norm <- with(res,improve/improve.tot*rsq)
  nodup <- !duplicated(res[,1:2])
  res.u <- res[nodup, c("spec","var","rsq","improve.tot","improve.tot.var")]
  res.u$rsq <- with(res.u, ifelse(is.na(rsq), 0, rsq))
  res.u$rsq.var <- with(res.u,rsq*improve.tot.var/improve.tot)
  list(res=res,res.u=res.u,dens=dens)
}


cumimp <- function (x, ...)
  UseMethod("cumimp")

cumimp.gradientForest <- function (x, predictor, type=c("Overall","Species")[1], standardize=TRUE, standardize_after=FALSE, ...)
  {
    if (!inherits(x,"gradientForest"))
      stop(paste("'x' must be a gradientForest object"))
    if (length(predictor) != 1)
      stop(paste("'predictor' must be a single string"))
    if (!is.element(predictor,levels(x$res$var)))
      stop(paste("Predictor",predictor,"does not belong to gradientForest object"))
    if (is.na(option <- pmatch(type,c("Overall","Species"))))
      stop(paste('Unmatched type "',type,'". Expecting one of "Overall" or "Species"',sep=''))

    # convert density to its inverse
    inverse <- function(dens) {dens$y <- 1/dens$y; dens}

    # crude integral
    crude.integrate <- function(f) sum(f$y)*diff(f$x)[1]

    # normalize f(x) to f(x)/fbar
    normalize <- function(f) {
      integral <- try(integrate(approxfun(f,rule=2),lower=min(f$x),upper=max(f$x))$value)
      if (class(integral)=="try-error") integral <- crude.integrate(f)
      f$y <- f$y/integral*diff(range(f$x));
      f
    }

    getCU <- function(importance.df, Rsq) {
      if (nrow(importance.df) == 0) {
        return( list(x=0, y=0))
      }
      agg <- with(importance.df, agg.sum(improve.norm, list(split), sort.it=TRUE))
      cum.split <- agg[,1]
      height <- agg[,2]
      if (standardize & standardize_after) # crucial to normalize this case
        dinv <- normalize(inverse(x$dens[[predictor]]))
      else dinv <- inverse(x$dens[[predictor]]) #
      dinv.vals <- approx(dinv$x, dinv$y, cum.split, rule=2)$y
      if (any(bad <- is.na(height))) {
        cum.split <- cum.split[!bad]
        height <- height[!bad]
        dinv.vals <- dinv.vals[!bad]
      }
      if (standardize & !standardize_after) height <- height * dinv.vals
      height <- height/sum(height)*Rsq
      if (standardize & standardize_after) height <- height * dinv.vals
      res <- list(x=cum.split, y=cumsum(height))
    }

    importance.df <- x$res[x$res$var==predictor,]
    if (option==1) {
      res <- getCU(importance.df, importance(x, "Weighted")[predictor])
    } else {
      species <- names(x$result)
      res <- lapply(namenames(species), function(sp)
        getCU(subset(importance.df, spec==sp), x$imp.rsq[predictor,sp]))
    }
    res
  }

cumimp.combinedGradientForest <-
  function (x, predictor, weight=c("uniform","species","rsq.total","rsq.mean")[3], gear, ...)
  {
    if (!inherits(x,"combinedGradientForest"))
      stop(paste("'x' must be a gradientForest object"))
    if (length(predictor) != 1)
      stop(paste("'predictor' must be a single string"))
    if (!is.element(predictor,names(x$X)[-1]))
      stop(paste("Predictor",predictor,"does not belong to combinedGradientForest object"))
    if (is.na(option <- pmatch(weight,c("uniform","species","rsq.total","rsq.mean"))))
      stop(paste('Unmatched weight "',weight,'". Expecting one of "uniform", "species", "rsq.total" or "rsq.mean"',sep=""))

    if (missing(gear)) {
      res <- x$CU[[predictor]][[paste("combined",weight,sep=".")]]
    } else {
      res <- x$CU[[predictor]][[gear]]
    }
    res
  }




