# In general, try to group together related functions into the same .R file
# https://tinyheero.github.io/jekyll/update/2015/07/26/making-your-first-R-package.html

# tag above your function to indicate this function to be “exposed” to users to use.


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
    fitcmd <- quote(randomForest(Species ~ rhs, data = cbind(Y,
                                                             X), maxLevel = maxLevel, keep.forest = TRUE, importance = TRUE,
                                 ntree = ntree, keep.group = TRUE, keep.inbag = TRUE,
                                 corr.threshold = corr.threshold, na.action = na.omit))
  else fitcmd <- quote(randomForest(Species ~ rhs, data = cbind(Y,
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



randomForest <- function (x, ...) UseMethod("randomForest")

## mylevels() returns levels if given a factor, otherwise 0.
mylevels <- function(x) if (is.factor(x)) levels(x) else 0

#' randomForest: Classification and Regression with Random Forest
#' randomForest implements Breiman's random forest algorithm (based on Breiman and Cutler's original Fortran code) for classification and regression.
#' It can also be used in unsupervised mode for assessing proximities among data points.
#' #'
#'
#' @param data an optional data frame containing the variables in the model. By default the variables are taken from the environment which randomForest is called from.
#' @param subset an index vector indicating which rows should be used. (NOTE: If given, this argument must be named.)
#' @param na.action A function to specify the action to be taken if NAs are found. (NOTE: If given, this argument must be named.)
#' @param x formula, a data frame or a matrix of predictors, or a formula describing the model to be fitted (for the print method, an randomForest object).
#' @param y A response vector. If a factor, classification is assumed, otherwise regression is assumed. If omitted, randomForest will run in unsupervised mode.
#' @param xtest a data frame or matrix (like x) containing predictors for the test set.
#' @param ytest response for the test set.
#' @param ntree Number of trees to grow. This should not be set to too small a number, to ensure that every input row gets predicted at least a few times.
#' @param mtry Number of variables randomly sampled as candidates at each split. Note that the default values are different for classification (sqrt(p) where p is number of variables in x) and regression (p/3)
#' @param replace Should sampling of cases be done with or without replacement?
#' @param classwt Priors of the classes. Need not add up to one. Ignored for regression.
#' @param cutoff (Classification only) A vector of length equal to number of classes. The ‘winning’ class for an observation is the one with the maximum ratio of proportion of votes to cutoff. Default is 1/k where k is the number of classes (i.e., majority vote wins).
#' @param strata A (factor) variable that is used for stratified sampling.
#' @param sampsize Size(s) of sample to draw. For classification, if sampsize is a vector of the length the number of strata, then sampling is stratified by strata, and the elements of sampsize indicate the numbers to be drawn from the strata.
#' @param nodesize Minimum size of terminal nodes. Setting this number larger causes smaller trees to be grown (and thus take less time). Note that the default values are different for classification (1) and regression (5).
#' @param importance Should importance of predictors be assessed?
#' @param localImp Should casewise importance measure be computed? (Setting this to TRUE will override importance.)
#' @param nPerm Number of times the OOB data are permuted per tree for assessing variable importance. Number larger than 1 gives slightly more stable estimate, but not very effective. Currently only implemented for regression.
#' @param proximity Should proximity measure among the rows be calculated?
#' @param oob.prox Should proximity be calculated only on “out-of-bag” data?
#' @param norm.votes If TRUE (default), the final result of votes are expressed as fractions. If FALSE, raw vote counts are returned (useful for combining results from different runs). Ignored for regression.
#' @param do.trace If set to TRUE, give a more verbose output as randomForest is run. If set to some integer, then running output is printed for every do.trace trees.
#' @param keep.forest If set to FALSE, the forest will not be retained in the output object. If xtest is given, defaults to FALSE.
#' @param corr.bias perform bias correction for regression? Note: Experimental. Use at your own risk.
#' @param keep.inbag Should an n by ntree matrix be returned that keeps track of which samples are “in-bag” in which trees (but not how many times, if sampling with replacement)
#' @param maxLevel If maxLevel == 0, compute importance from marginal permutation distribution of each variable (the default). If maxLevel > 0, compute importance from conditional permutation distribution of each variable, permuted within 2^maxLevel partitions of correlated variables.
#' @param corr.threshold If maxLevel > 0, OOB permuting is conditioned on partitions of variables having absolute correlation > corr.threshold.
#' @param corr.method Method for computing correlation between variables. Default "pearson".
#' @param keep.group If TRUE keep diagnostic information on the partitioning for the importance calculation.
#' @param ... optional parameters to be passed to the low level function randomForest.default.
#'
#'
#' @return data.frame of predictor variables.
#' @examples
#' data(CoMLsimulation)
#' preds <- colnames(Xsimulation)
#' specs <- colnames(Ysimulation)
#' f1 <- gradientForest(data.frame(Ysimulation,Xsimulation), preds, specs, ntree=10)
#' f1
#' @export
randomForest.default <- function (x, y = NULL, xtest = NULL, ytest = NULL, ntree = 500,
                                  mtry = if (!is.null(y) && !is.factor(y)) max(floor(ncol(x)/3),
                                                                               1) else floor(sqrt(ncol(x))), replace = TRUE, classwt = NULL,
                                  cutoff, strata, sampsize = if (replace) nrow(x) else ceiling(0.632 *
                                                                                                 nrow(x)), nodesize = if (!is.null(y) && !is.factor(y)) 5 else 1,
                                  importance = FALSE, localImp = FALSE, nPerm = 1, proximity,
                                  oob.prox = proximity, norm.votes = TRUE, do.trace = FALSE,
                                  keep.forest = !is.null(y) && is.null(xtest), corr.bias = FALSE,
                                  keep.inbag = FALSE, maxLevel = 0, keep.group = FALSE, corr.threshold = 1,
                                  corr.method = "pearson", ...) {
  addclass <- is.null(y)
  classRF <- addclass || is.factor(y)
  if (!classRF && length(unique(y)) <= 5) {
    warning("The response has five or fewer unique values.  Are you sure you want to do regression?")
  }
  if (classRF && !addclass && length(unique(y)) < 2)
    stop("Need at least two classes to do classification.")
  n <- nrow(x)
  p <- ncol(x)
  if (n == 0)
    stop("data (x) has 0 rows")
  x.row.names <- rownames(x)
  x.col.names <- if (is.null(colnames(x)))
    1:ncol(x)
  else colnames(x)
  keep.forest <- keep.forest
  testdat <- !is.null(xtest)
  if (testdat) {
    if (ncol(x) != ncol(xtest))
      stop("x and xtest must have same number of columns")
    ntest <- nrow(xtest)
    xts.row.names <- rownames(xtest)
  }
  if (mtry < 1 || mtry > p)
    warning("invalid mtry: reset to within valid range")
  mtry <- max(1, min(p, round(mtry)))
  if (!is.null(y)) {
    if (length(y) != n)
      stop("length of response must be the same as predictors")
    addclass <- FALSE
  }
  else {
    if (!addclass)
      addclass <- TRUE
    y <- factor(c(rep(1, n), rep(2, n)))
    x <- rbind(x, x)
  }
  if (any(is.na(x)))
    stop("NA not permitted in predictors")
  if (testdat && any(is.na(xtest)))
    stop("NA not permitted in xtest")
  if (any(is.na(y)))
    stop("NA not permitted in response")
  if (!is.null(ytest) && any(is.na(ytest)))
    stop("NA not permitted in ytest")
  if (is.data.frame(x)) {
    xc <- as.data.frame(lapply(x, function(y) if (is.numeric(y))
      y
      else if (is.ordered(y))
        as.numeric(y)
      else rep(NA, length(y))))
    corr <- cor(xc, method = corr.method)
    corr[is.na(corr)] <- 0
    xlevels <- lapply(x, mylevels)
    ncat <- sapply(xlevels, length)
    ncat <- ifelse(sapply(x, is.ordered), 1, ncat)
    x <- data.matrix(x)
    if (testdat) {
      if (!is.data.frame(xtest))
        stop("xtest must be data frame if x is")
      xfactor <- which(sapply(xtest, is.factor))
      if (length(xfactor) > 0) {
        for (i in xfactor) {
          if (any(!levels(xtest[[i]]) %in% xlevels[[i]]))
            stop("New factor levels in xtest not present in x")
          xtest[[i]] <- factor(xlevels[[i]][match(xtest[[i]],
                                                  xlevels[[i]])], levels = xlevels[[i]])
        }
      }
      xtest <- data.matrix(xtest)
    }
  }
  else {
    corr <- cor(x, method = corr.method)
    ncat <- rep(1, p)
    xlevels <- as.list(rep(0, p))
  }
  maxcat <- max(ncat)
  if (maxcat > 32)
    stop("Can not handle categorical predictors with more than 32 categories.")
  if (classRF) {
    nclass <- length(levels(y))
    if (any(table(y) == 0))
      stop("Can't have empty classes in y.")
    if (!is.null(ytest)) {
      if (!is.factor(ytest))
        stop("ytest must be a factor")
      if (!all(levels(y) == levels(ytest)))
        stop("y and ytest must have the same levels")
    }
    if (missing(cutoff)) {
      cutoff <- rep(1/nclass, nclass)
    }
    else {
      if (sum(cutoff) > 1 || sum(cutoff) < 0 || !all(cutoff >
                                                     0) || length(cutoff) != nclass) {
        stop("Incorrect cutoff specified.")
      }
      if (!is.null(names(cutoff))) {
        if (!all(names(cutoff) %in% levels(y))) {
          stop("Wrong name(s) for cutoff")
        }
        cutoff <- cutoff[levels(y)]
      }
    }
    if (!is.null(classwt)) {
      if (length(classwt) != nclass)
        stop("length of classwt not equal to number of classes")
      if (!is.null(names(classwt))) {
        if (!all(names(classwt) %in% levels(y))) {
          stop("Wrong name(s) for classwt")
        }
        classwt <- classwt[levels(y)]
      }
      if (any(classwt <= 0))
        stop("classwt must be positive")
      ipi <- 1
    }
    else {
      classwt <- rep(1, nclass)
      ipi <- 0
    }
  }
  else addclass <- FALSE
  if (missing(proximity))
    proximity <- addclass
  if (proximity) {
    prox <- matrix(0, n, n)
    proxts <- if (testdat)
      matrix(0, ntest, ntest + n)
    else double(1)
  }
  else {
    prox <- proxts <- double(1)
  }
  if (localImp) {
    importance <- TRUE
    impmat <- matrix(0, p, n)
  }
  else impmat <- double(1)
  if (importance) {
    if (nPerm < 1)
      nPerm <- as.integer(1)
    else nPerm <- as.integer(nPerm)
    if (classRF) {
      impout <- matrix(0, p, nclass + 2)
      impSD <- matrix(0, p, nclass + 1)
    }
    else {
      impout <- matrix(0, p, 2)
      impSD <- double(p)
      names(impSD) <- x.col.names
    }
  }
  else {
    impout <- double(p)
    impSD <- double(1)
  }
  nsample <- if (addclass)
    2 * n
  else n
  Stratify <- length(sampsize) > 1
  if ((!Stratify) && sampsize > nrow(x))
    stop("sampsize too large")
  if (Stratify && (!classRF))
    stop("sampsize should be of length one")
  if (classRF) {
    if (Stratify) {
      if (missing(strata))
        strata <- y
      if (!is.factor(strata))
        strata <- as.factor(strata)
      nsum <- sum(sampsize)
      if (length(sampsize) > nlevels(strata))
        stop("sampsize has too many elements.")
      if (any(sampsize <= 0) || nsum == 0)
        stop("Bad sampsize specification")
      if (!is.null(names(sampsize))) {
        sampsize <- sampsize[levels(strata)]
      }
      if (any(sampsize > table(strata)))
        stop("sampsize can not be larger than class frequency")
    }
    else {
      nsum <- sampsize
    }
    nrnodes <- 2 * trunc(nsum/nodesize) + 1
  }
  else {
    nrnodes <- 2 * trunc(sampsize/max(1, nodesize - 4)) +
      1
  }
  x <- t(x)
  storage.mode(x) <- "double"
  if (testdat) {
    xtest <- t(xtest)
    storage.mode(xtest) <- "double"
    if (is.null(ytest)) {
      ytest <- labelts <- 0
    }
    else {
      labelts <- TRUE
    }
  }
  else {
    xtest <- double(1)
    ytest <- double(1)
    ntest <- 1
    labelts <- FALSE
  }
  nt <- if (keep.forest)
    ntree
  else 1
  if (classRF) {
    error.test <- if (labelts)
      double((nclass + 1) * ntree)
    else double(1)
    rfout <- .C("classRF", x = x, xdim = as.integer(c(p,
                                                      n)), y = as.integer(y), nclass = as.integer(nclass),
                ncat = as.integer(ncat), maxcat = as.integer(maxcat),
                sampsize = as.integer(sampsize), strata = if (Stratify) as.integer(strata) else integer(1),
                Options = as.integer(c(addclass, importance, localImp,
                                       proximity, oob.prox, do.trace, keep.forest, replace,
                                       Stratify, keep.inbag, keep.group)), ntree = as.integer(ntree),
                mtry = as.integer(mtry), ipi = as.integer(ipi), classwt = as.double(classwt),
                cutoff = as.double(cutoff), nodesize = as.integer(nodesize),
                outcl = integer(nsample), counttr = integer(nclass *
                                                              nsample), prox = prox, impout = impout, impSD = impSD,
                impmat = impmat, nrnodes = as.integer(nrnodes), ndbigtree = integer(ntree),
                nodestatus = integer(nt * nrnodes), bestvar = integer(nt *
                                                                        nrnodes), treemap = integer(nt * 2 * nrnodes),
                nodepred = integer(nt * nrnodes), xbestsplit = double(nt *
                                                                        nrnodes), errtr = double((nclass + 1) * ntree),
                testdat = as.integer(testdat), xts = as.double(xtest),
                clts = as.integer(ytest), nts = as.integer(ntest),
                countts = double(nclass * ntest), outclts = as.integer(numeric(ntest)),
                labelts = as.integer(labelts), proxts = proxts, errts = error.test,
                inbag = if (keep.inbag) matrix(integer(n * ntree),
                                               n) else integer(n), nodeImprove = double(nt *
                                                                                          nrnodes), as.integer(maxLevel), group = if (keep.group) matrix(integer(n *
                                                                                                                                                                   p), n) else integer(1), permX = if (keep.group) matrix(double(n *
                                                                                                                                                                                                                                   p), n) else double(1), as.integer(abs(corr) >
                                                                                                                                                                                                                                                                       corr.threshold), DUP = FALSE, PACKAGE = "extendedForest")[-1]
    if (keep.forest) {
      max.nodes <- max(rfout$ndbigtree)
      treemap <- aperm(array(rfout$treemap, dim = c(2,
                                                    nrnodes, ntree)), c(2, 1, 3))[1:max.nodes, ,
                                                                                  , drop = FALSE]
    }
    if (!addclass) {
      out.class <- factor(rfout$outcl, levels = 1:nclass,
                          labels = levels(y))
      names(out.class) <- x.row.names
      con <- table(observed = y, predicted = out.class)[levels(y),
                                                        levels(y)]
      con <- cbind(con, class.error = 1 - diag(con)/rowSums(con))
    }
    out.votes <- t(matrix(rfout$counttr, nclass, nsample))[1:n,
                                                           ]
    oob.times <- rowSums(out.votes)
    if (norm.votes)
      out.votes <- t(apply(out.votes, 1, function(x) x/sum(x)))
    dimnames(out.votes) <- list(x.row.names, levels(y))
    if (testdat) {
      out.class.ts <- factor(rfout$outclts, levels = 1:nclass,
                             labels = levels(y))
      names(out.class.ts) <- xts.row.names
      out.votes.ts <- t(matrix(rfout$countts, nclass, ntest))
      dimnames(out.votes.ts) <- list(xts.row.names, levels(y))
      if (norm.votes)
        out.votes.ts <- t(apply(out.votes.ts, 1, function(x) x/sum(x)))
      if (labelts) {
        testcon <- table(observed = ytest, predicted = out.class.ts)[levels(y),
                                                                     levels(y)]
        testcon <- cbind(testcon, class.error = 1 - diag(testcon)/rowSums(testcon))
      }
    }
    cl <- match.call()
    cl[[1]] <- as.name("randomForest")
    out <- list(call = cl, type = if (addclass) "unsupervised" else "classification",
                predicted = if (addclass) NULL else out.class, err.rate = if (addclass) NULL else t(matrix(rfout$errtr,
                                                                                                           nclass + 1, ntree, dimnames = list(c("OOB", levels(y)),
                                                                                                                                              NULL))), confusion = if (addclass) NULL else con,
                votes = out.votes, oob.times = oob.times, classes = levels(y),
                importance = if (importance) matrix(rfout$impout,
                                                    p, nclass + 2, dimnames = list(x.col.names, c(levels(y),
                                                                                                  "MeanDecreaseAccuracy", "MeanDecreaseGini"))) else matrix(rfout$impout,
                                                                                                                                                            ncol = 1, dimnames = list(x.col.names, "MeanDecreaseGini")),
                importanceSD = if (importance) matrix(rfout$impSD,
                                                      p, nclass + 1, dimnames = list(x.col.names, c(levels(y),
                                                                                                    "MeanDecreaseAccuracy"))) else NULL, localImportance = if (localImp) matrix(rfout$impmat,
                                                                                                                                                                                p, n, dimnames = list(x.col.names, x.row.names)) else NULL,
                proximity = if (proximity) matrix(rfout$prox, n,
                                                  n, dimnames = list(x.row.names, x.row.names)) else NULL,
                ntree = ntree, mtry = mtry, forest = if (!keep.forest) NULL else {
                  list(ndbigtree = rfout$ndbigtree, nodestatus = matrix(rfout$nodestatus,
                                                                        ncol = ntree)[1:max.nodes, , drop = FALSE],
                       bestvar = matrix(rfout$bestvar, ncol = ntree)[1:max.nodes,
                                                                     , drop = FALSE], treemap = treemap, nodepred = matrix(rfout$nodepred,
                                                                                                                           ncol = ntree)[1:max.nodes, , drop = FALSE],
                       xbestsplit = matrix(rfout$xbestsplit, ncol = ntree)[1:max.nodes,
                                                                           , drop = FALSE], pid = rfout$classwt, cutoff = cutoff,
                       ncat = ncat, maxcat = maxcat, nrnodes = max.nodes,
                       ntree = ntree, nclass = nclass, xlevels = xlevels,
                       nodeImprove = matrix(rfout$nodeImprove, ncol = ntree)[1:max.nodes,
                                                                             , drop = FALSE])
                }, y = if (addclass) NULL else y, test = if (!testdat) NULL else list(predicted = out.class.ts,
                                                                                      err.rate = if (labelts) t(matrix(rfout$errts,
                                                                                                                       nclass + 1, ntree, dimnames = list(c("Test",
                                                                                                                                                            levels(y)), NULL))) else NULL, confusion = if (labelts) testcon else NULL,
                                                                                      votes = out.votes.ts, proximity = if (proximity) matrix(rfout$proxts,
                                                                                                                                              nrow = ntest, dimnames = list(xts.row.names,
                                                                                                                                                                            c(xts.row.names, x.row.names))) else NULL),
                inbag = if (keep.inbag) rfout$inbag else NULL, group = if (keep.group) matrix(rfout$group,
                                                                                              nrow(rfout$group), dimnames = list(x.row.names,
                                                                                                                                 x.col.names)) else NULL, permX = if (keep.group) matrix(rfout$permX,
                                                                                                                                                                                         nrow(rfout$permX), dimnames = list(x.row.names,
                                                                                                                                                                                                                            x.col.names)) else NULL)
  }
  else {
    rfout <- .C("regRF", x, as.double(y), as.integer(c(n, p)), as.integer(sampsize), as.integer(nodesize),
                as.integer(nrnodes), as.integer(ntree), as.integer(mtry),
                as.integer(c(importance, localImp, nPerm)), as.integer(ncat),
                as.integer(maxcat), as.integer(do.trace), as.integer(proximity),
                as.integer(oob.prox), as.integer(corr.bias), ypred = double(n),
                impout = impout, impmat = impmat, impSD = impSD,
                prox = prox, ndbigtree = integer(ntree), nodestatus = matrix(integer(nrnodes * nt), ncol = nt), leftDaughter = matrix(integer(nrnodes * nt), ncol = nt), rightDaughter = matrix(integer(nrnodes *  nt), ncol = nt), nodepred = matrix(double(nrnodes *  nt), ncol = nt), bestvar = matrix(integer(nrnodes *  nt), ncol = nt), xbestsplit = matrix(double(nrnodes * nt), ncol = nt), mse = double(ntree), keep = as.integer(c(keep.forest,
                                                                                                                                                                                                                                                                                                                                                                                          keep.inbag, keep.group)), replace = as.integer(replace),
                testdat = as.integer(testdat), xts = xtest, ntest = as.integer(ntest),
                yts = as.double(ytest), labelts = as.integer(labelts),
                ytestpred = double(ntest), proxts = proxts, msets = double(if (labelts) ntree else 1),
                coef = double(2), oob.times = integer(n), inbag = if (keep.inbag) matrix(integer(n *   ntree), n) else integer(1), nodeSS = matrix(double(nrnodes *  nt), ncol = nt), nodeImprove = matrix(double(nrnodes * nt), ncol = nt), as.integer(maxLevel), group = if (keep.group) matrix(integer(n *   p), n) else integer(1), permX = if (keep.group) matrix(double(n * p), n) else double(1), as.integer(abs(corr) > corr.threshold), DUP = FALSE, PACKAGE = "extendedForest")[c(16:28, 36:43, 45:46)]
    if (keep.forest) {
      max.nodes <- max(rfout$ndbigtree)
      rfout$nodestatus <- rfout$nodestatus[1:max.nodes,
                                           , drop = FALSE]
      rfout$bestvar <- rfout$bestvar[1:max.nodes, , drop = FALSE]
      rfout$nodepred <- rfout$nodepred[1:max.nodes, , drop = FALSE]
      rfout$nodeSS <- rfout$nodeSS[1:max.nodes, , drop = FALSE]
      rfout$nodeImprove <- rfout$nodeImprove[1:max.nodes,
                                             , drop = FALSE]
      rfout$xbestsplit <- rfout$xbestsplit[1:max.nodes,
                                           , drop = FALSE]
      rfout$leftDaughter <- rfout$leftDaughter[1:max.nodes,
                                               , drop = FALSE]
      rfout$rightDaughter <- rfout$rightDaughter[1:max.nodes,
                                                 , drop = FALSE]
    }
    cl <- match.call()
    cl[[1]] <- as.name("randomForest")
    out <- list(call = cl, type = "regression", predicted = structure(rfout$ypred,
                                                                      names = x.row.names), mse = rfout$mse, rsq = 1 -
                  rfout$mse/(var(y) * (n - 1)/n), oob.times = rfout$oob.times,
                importance = if (importance) matrix(rfout$impout,
                                                    p, 2, dimnames = list(x.col.names, c("%IncMSE",
                                                                                         "IncNodePurity"))) else matrix(rfout$impout,
                                                                                                                        ncol = 1, dimnames = list(x.col.names, "IncNodePurity")),
                importanceSD = if (importance) rfout$impSD else NULL,
                localImportance = if (localImp) matrix(rfout$impmat,
                                                       p, n, dimnames = list(x.col.names, x.row.names)) else NULL,
                proximity = if (proximity) matrix(rfout$prox, n,
                                                  n, dimnames = list(x.row.names, x.row.names)) else NULL,
                ntree = ntree, mtry = mtry, forest = if (keep.forest) c(rfout[c("ndbigtree",
                                                                                "nodestatus", "leftDaughter", "rightDaughter",
                                                                                "nodepred", "bestvar", "xbestsplit", "nodeSS",
                                                                                "nodeImprove")], list(ncat = ncat), list(nrnodes = max.nodes),
                                                                        list(ntree = ntree), list(xlevels = xlevels)) else NULL,
                coefs = if (corr.bias) rfout$coef else NULL, y = y,
                test = if (testdat) {
                  list(predicted = structure(rfout$ytestpred, names = xts.row.names),
                       mse = if (labelts) rfout$msets else NULL, rsq = if (labelts) 1 -
                         rfout$msets/(var(ytest) * (n - 1)/n) else NULL,
                       proximity = if (proximity) matrix(rfout$proxts/ntree,
                                                         nrow = ntest, dimnames = list(xts.row.names,
                                                                                       c(xts.row.names, x.row.names))) else NULL)
                } else NULL, inbag = if (keep.inbag) matrix(rfout$inbag,
                                                            nrow(rfout$inbag), dimnames = list(x.row.names,
                                                                                               NULL)) else NULL, group = if (keep.group) matrix(rfout$group,
                                                                                                                                                nrow(rfout$group), dimnames = list(x.row.names,
                                                                                                                                                                                   x.col.names)) else NULL, permX = if (keep.group) matrix(rfout$permX,
                                                                                                                                                                                                                                           nrow(rfout$permX), dimnames = list(x.row.names,
                                                                                                                                                                                                                                                                              x.col.names)) else NULL)
  }
  class(out) <- "randomForest"
  return(out)
}



randomForest.formula <-
  function(formula, data = NULL, ..., subset, na.action = na.fail) {
    ### formula interface for randomForest.
    ### code gratefully stolen from svm.formula (package e1071).
    ###
    if (!inherits(formula, "formula"))
      stop("method is only for formula objects")
    m <- match.call(expand = FALSE)
    ## Catch xtest and ytest in arguments.
    if (any(c("xtest", "ytest") %in% names(m)))
      stop("xtest/ytest not supported through the formula interface")
    names(m)[2] <- "formula"
    if (is.matrix(eval(m$data, parent.frame())))
      m$data <- as.data.frame(data)
    m$... <- NULL
    m$na.action <- na.action
    m[[1]] <- as.name("model.frame")
    m <- eval(m, parent.frame())
    y <- model.response(m)
    Terms <- attr(m, "terms")
    attr(Terms, "intercept") <- 0
    ## Drop any "negative" terms in the formula.
    ## test with:
    ## randomForest(Fertility~.-Catholic+I(Catholic<50),data=swiss,mtry=2)
    m <- model.frame(terms(reformulate(attributes(Terms)$term.labels)),
                     data.frame(m))
    ## if (!is.null(y)) m <- m[, -1, drop=FALSE]
    for (i in seq(along=ncol(m))) {
      if (is.ordered(m[[i]])) m[[i]] <- as.numeric(m[[i]])
    }
    ret <- randomForest(m, y, ...)
    cl <- match.call()
    cl[[1]] <- as.name("randomForest")
    ret$call <- cl
    ret$terms <- Terms
    if (!is.null(attr(m, "na.action")))
      ret$na.action <- attr(m, "na.action")
    class(ret) <- c("randomForest.formula", "randomForest")
    return(ret)
  }


rfImpute <- function(x, ...)
  UseMethod("rfImpute")

rfImpute.formula <- function(x, data, ..., subset) {
  if (!inherits(x, "formula"))
    stop("method is only for formula objects")
  call <- match.call()
  m <- match.call(expand.dots = FALSE)
  names(m)[2] <- "formula"
  if (is.matrix(eval(m$data, parent.frame())))
    m$data <- as.data.frame(data)
  m$... <- NULL
  m$na.action <- as.name("na.pass")
  m[[1]] <- as.name("model.frame")
  m <- eval(m, parent.frame())
  Terms <- attr(m, "terms")
  attr(Terms, "intercept") <- 0
  y <- model.response(m)
  if (!is.null(y)) m <- m[,-1]
  for (i in seq(along=ncol(m))) {
    if(is.ordered(m[[i]])) m[[i]] <- as.numeric(m[[i]])
  }
  ret <- rfImpute.default(m, y, ...)
  names(ret)[1] <- deparse(as.list(x)[[2]])
  ret
}

rfImpute.default <- function(x, y, iter=5, ntree=300, ...) {
  if (any(is.na(y))) stop("Can't have NAs in", deparse(substitute(y)))
  if (!any(is.na(x))) stop("No NAs found in ", deparse(substitute(x)))
  xf <- na.roughfix(x)
  hasNA <- which(apply(x, 2, function(x) any(is.na(x))))
  if (is.data.frame(x)) {
    isfac <- sapply(x, is.factor)
  } else {
    isfac <- rep(FALSE, ncol(x))
  }

  for (i in 1:iter) {
    prox <- randomForest(xf, y, ntree=ntree, ..., do.trace=ntree,
                         proximity=TRUE)$proximity
    for (j in hasNA) {
      miss <- which(is.na(x[, j]))
      if (isfac[j]) {
        lvl <- levels(x[[j]])
        catprox <- apply(prox[-miss, miss, drop=FALSE], 2,
                         function(v) lvl[which.max(tapply(v, x[[j]][-miss], mean))])
        xf[miss, j] <- catprox
      } else {
        sumprox <- colSums(prox[-miss, miss, drop=FALSE])
        xf[miss, j] <- (prox[miss, -miss, drop=FALSE] %*% xf[,j][-miss]) / sumprox
      }
      NULL
    }
  }
  xf <- cbind(y, xf)
  names(xf)[1] <- deparse(substitute(y))
  xf
}


classCenter <- function(x, label, prox, nNbr = min(table(label))-1) {
  ## nPrototype=rep(3, length(unique(label))), ...) {
  label <- as.character(label)
  clsLabel <- unique(label)
  ## Find the nearest nNbr neighbors of each case
  ## (including the case itself).
  idx <- t(apply(prox, 1, order, decreasing=TRUE)[1:nNbr,])
  ## Find the class labels of the neighbors.
  cls <- label[idx]
  dim(cls) <- dim(idx)
  ## Count the number of neighbors in each class for each case.
  ncls <- sapply(clsLabel, function(x) rowSums(cls == x))
  ## For each class, find the case(s) with most neighbors in that class.
  clsMode <- max.col(t(ncls))
  ## Identify the neighbors of the class modes that are of the target class.
  nbrList <- mapply(function(cls, m) idx[m,][label[idx[m,]] == cls],
                    clsLabel, clsMode, SIMPLIFY=FALSE)
  ## Get the X data for the neighbors of the class `modes'.
  xdat <- t(sapply(nbrList, function(i) apply(x[i,,drop=FALSE], 2,
                                              median)))
  xdat
}

MDSplot <- function(rf, fac, k=2, palette=NULL, pch=20, ...) {
  if (!inherits(rf, "randomForest"))
    stop(deparse(substitute(rf)), " must be a randomForest object")
  if(is.null(rf$proximity))
    stop(deparse(substitute(rf)), " does not contain a proximity matrix")
  op <- par(pty="s")
  on.exit(par(op))
  rf.mds <- stats:::cmdscale(1 - rf$proximity, eig=TRUE, k=k)
  colnames(rf.mds$points) <- paste("Dim", 1:k)
  nlevs <- nlevels(fac)
  if (is.null(palette)) {
    palette <- if (require(RColorBrewer) && nlevs < 12)
      brewer.pal(nlevs, "Set1") else rainbow(nlevs)
  }
  if (k <= 2) {
    plot(rf.mds$points, col=palette[as.numeric(fac)], pch=pch, ...)
  } else {
    pairs(rf.mds$points, col=palette[as.numeric(fac)], pch=pch, ...)
  }
  invisible(rf.mds)
}


combine <- function(...) {
  pad0 <- function(x, len) c(x, rep(0, len-length(x)))
  padm0 <- function(x, len) rbind(x, matrix(0, nrow=len-nrow(x),
                                            ncol=ncol(x)))
  rflist <- list(...)
  areForest <- sapply(rflist, function(x) inherits(x, "randomForest"))
  if (any(!areForest)) stop("Argument must be a list of randomForest objects")
  ## Use the first component as a template
  rf <- rflist[[1]]
  classRF <- rf$type == "classification"
  trees <- sapply(rflist, function(x) x$ntree)
  ntree <- sum(trees)
  rf$ntree <- ntree
  nforest <- length(rflist)
  haveTest <- ! any(sapply(rflist, function(x) is.null(x$test)))
  ## Check if predictor variables are identical.
  vlist <- lapply(rflist, function(x) rownames(importance(x)))
  numvars <- sapply(vlist, length)
  if (! all(numvars[1] == numvars[-1]))
    stop("Unequal number of predictor variables in the randomForest objects.")
  for (i in seq_along(vlist)) {
    if (! all(vlist[[i]] == vlist[[1]]))
      stop("Predictor variables are different in the randomForest objects.")
  }
  ## Combine the forest component, if any
  haveForest <- sapply(rflist, function(x) !is.null(x$forest))
  if (all(haveForest)) {
    nrnodes <- max(sapply(rflist, function(x) x$forest$nrnodes))
    rf$forest$nrnodes <- nrnodes
    rf$forest$ndbigtree <-
      unlist(sapply(rflist, function(x) x$forest$ndbigtree))
    rf$forest$nodestatus <-
      do.call("cbind", lapply(rflist, function(x)
        padm0(x$forest$nodestatus, nrnodes)))
    rf$forest $bestvar <-
      do.call("cbind",
              lapply(rflist, function(x)
                padm0(x$forest$bestvar, nrnodes)))
    rf$forest$xbestsplit <-
      do.call("cbind",
              lapply(rflist, function(x)
                padm0(x$forest$xbestsplit, nrnodes)))
    rf$forest$nodepred <-
      do.call("cbind", lapply(rflist, function(x)
        padm0(x$forest$nodepred, nrnodes)))
    tree.dim <- dim(rf$forest$treemap)
    if (classRF) {
      rf$forest$treemap <-
        array(unlist(lapply(rflist, function(x) apply(x$forest$treemap, 2:3,
                                                      pad0, nrnodes))),
              c(nrnodes, 2, ntree))
    } else {
      rf$forest$leftDaughter <-
        do.call("cbind",
                lapply(rflist, function(x)
                  padm0(x$forest$leftDaughter, nrnodes)))
      rf$forest$rightDaughter <-
        do.call("cbind",
                lapply(rflist, function(x)
                  padm0(x$forest$rightDaughter, nrnodes)))
    }
    rf$forest$ntree <- ntree
    if (classRF) rf$forest$cutoff <- rflist[[1]]$forest$cutoff
  } else {
    rf$forest <- NULL
  }

  if (classRF) {
    ## Combine the votes matrix:
    rf$votes <- 0
    rf$oob.times <- 0
    areVotes <- all(sapply(rflist, function(x) any(x$votes > 1, na.rf=TRUE)))
    if (areVotes) {
      for(i in 1:nforest) {
        rf$oob.times <- rf$oob.times + rflist[[i]]$oob.times
        rf$votes <- rf$votes +
          ifelse(is.na(rflist[[i]]$votes), 0, rflist[[i]]$votes)
      }
    } else {
      for(i in 1:nforest) {
        rf$oob.times <- rf$oob.times + rflist[[i]]$oob.times
        rf$votes <- rf$votes +
          ifelse(is.na(rflist[[i]]$votes), 0, rflist[[i]]$votes) *
          rflist[[i]]$oob.times
      }
      rf$votes <- rf$votes / rf$oob.times
    }
    rf$predicted <- factor(colnames(rf$votes)[max.col(rf$votes)],
                           levels=levels(rf$predicted))
    if(haveTest) {
      rf$test$votes <- 0
      if (any(rf$test$votes > 1)) {
        for(i in 1:nforest)
          rf$test$votes <- rf$test$votes + rflist[[i]]$test$votes
      } else {
        for (i in 1:nforest)
          rf$test$votes <- rf$test$votes +
            rflist[[i]]$test$votes * rflist[[i]]$ntree
      }
      rf$test$predicted <-
        factor(colnames(rf$test$votes)[max.col(rf$test$votes)],
               levels=levels(rf$test$predicted))
    }
  } else {
    rf$predicted <- 0
    for (i in 1:nforest) rf$predicted <- rf$predicted +
        rflist[[i]]$predicted * rflist[[i]]$ntree
    rf$predicted <- rf$predicted / ntree
    if (haveTest) {
      rf$test$predicted <- 0
      for (i in 1:nforest) rf$test$predicted <- rf$test$predicted +
          rflist[[i]]$test$predicted * rflist[[i]]$ntree
      rf$test$predicted <- rf$test$predicted / ntree
    }
  }

  ## If variable importance is in all of them, compute the average
  ## (weighted by the number of trees in each forest)
  have.imp <- !any(sapply(rflist, function(x) is.null(x$importance)))
  if (have.imp) {
    rf$importance <- rf$importanceSD <- 0
    for(i in 1:nforest) {
      rf$importance <- rf$importance +
        rflist[[i]]$importance * rflist[[i]]$ntree
      ## Do the same thing with SD of importance, though that's not
      ## exactly right...
      rf$importanceSD <- rf$importanceSD +
        rflist[[i]]$importanceSD^2 * rflist[[i]]$ntree
    }
    rf$importance <- rf$importance / ntree
    rf$importanceSD <- sqrt(rf$importanceSD / ntree)
    haveCaseImp <- !any(sapply(rflist, function(x)
      is.null(x$localImportance)))
    ## Average casewise importance
    if (haveCaseImp) {
      rf$localImportance <- 0
      for (i in 1:nforest) {
        rf$localImportance <- rf$localImportance +
          rflist[[i]]$localImportance * rflist[[i]]$ntree
      }
      rf$localImportance <- rf$localImportance / ntree
    }
  }

  ## If proximity is in all of them, compute the average
  ## (weighted by the number of trees in each forest)
  have.prox <- !any(sapply(rflist, function(x) is.null(x$proximity)))
  if(have.prox) {
    rf$proximity <- 0
    for(i in 1:nforest)
      rf$proximity <- rf$proximity + rflist[[i]]$proximity * rflist[[i]]$ntree
    rf$proximity <- rf$proximity / ntree
  }

  ## Set confusion matrix and error rates to NULL
  if(classRF) {
    rf$confusion <- NULL
    rf$err.rate <- NULL
    if(haveTest) {
      rf$test$confusion <- NULL
      rf$err.rate <- NULL
    }
  } else {
    rf$mse <- rf$rsq <- NULL
    if(haveTest) rf$test$mse <- rf$test$rsq <- NULL
  }
  rf
}

getTree <- function(rfobj, k=1, labelVar=FALSE) {
  if (is.null(rfobj$forest)) {
    stop("No forest component in ", deparse(substitute(rfobj)))
  }
  if (k > rfobj$ntree) {
    stop("There are fewer than ", k, "trees in the forest")
  }
  if (rfobj$type == "regression") {
    tree <- cbind(rfobj$forest$leftDaughter[,k],
                  rfobj$forest$rightDaughter[,k],
                  rfobj$forest$bestvar[,k],
                  rfobj$forest$xbestsplit[,k],
                  rfobj$forest$nodestatus[,k],
                  rfobj$forest$nodepred[,k],
                  rfobj$forest$nodeSS[,k],
                  rfobj$forest$nodeImprove[,k])[1:rfobj$forest$ndbigtree[k],,drop=F]
    dimnames(tree) <- list(1:nrow(tree), c("left daughter", "right daughter",
                                           "split var", "split point", "status",
                                           "prediction", "SS", "improve"))
  } else {
    tree <- cbind(rfobj$forest$treemap[,,k],
                  rfobj$forest$rightDaughter[,k],
                  rfobj$forest$bestvar[,k],
                  rfobj$forest$xbestsplit[,k],
                  rfobj$forest$nodestatus[,k],
                  rfobj$forest$nodepred[,k],
                  rfobj$forest$nodeImprove[,k])[1:rfobj$forest$ndbigtree[k],,drop=F]
    dimnames(tree) <- list(1:nrow(tree), c("left daughter", "right daughter",
                                           "split var", "split point", "status",
                                           "prediction", "improve"))
  }


  if (labelVar) {
    tree <- as.data.frame(tree)
    v <- tree[[3]]
    v[v == 0] <- NA
    tree[[3]] <- factor(rownames(rfobj$importance)[v])
    if (rfobj$type == "classification") {
      v <- tree[[6]]
      v[! v %in% 1:nlevels(rfobj$y)] <- NA
      tree[[6]] <- levels(rfobj$y)[v]
    }
  }
  tree
}

grow <- function(x, ...) UseMethod("grow")

grow.default <- function(x, ...)
  stop("grow has not been implemented for this class of object")

grow.randomForest <- function(x, how.many, ...) {
  y <- update(x, ntree=how.many)
  combine(x, y)
}


margin <- function(rf, observed) {
  if( !inherits(rf, "randomForest") ) {
    stop("margin defined for Random Forests")
  }
  if( is.null(rf$votes) ) {
    stop("margin is only defined if votes are present")
  }
  if( !is.factor(observed) ) {
    stop(deparse(substitute(observed)), " is not a factor")
  }
  augD <- rf$votes
  if( any(augD > 1) ) {
    augD <- sweep(augD, 1, rowSums(augD), "/")
  }
  augD <- data.frame(augD, observed)
  names(augD) <- c(dimnames(rf$votes)[[2]], "observed")
  nlev <- length(levels(observed))

  ans<- apply(augD, 1, function(x) { pos <- match(x[nlev+1], names(x));
  t1 <- as.numeric(x[pos]);
  t2 <- max(as.numeric(x[-c(pos, nlev+1)]));
  t1 - t2 }
  )
  names(ans) <- observed
  class(ans) <- "margin"
  ans
}

plot.margin <- function(x, sort=TRUE, ...) {
  if (sort) x <- sort(x)
  nF <- factor(names(x))
  nlevs <- length(levels(nF))
  if ( require(RColorBrewer) && nlevs < 12) {
    pal <- brewer.pal(nlevs,"Set1")
  } else {
    pal <- rainbow(nlevs)
  }
  plot.default(x, col=pal[as.numeric(nF)], pch=20, ... )
}


na.roughfix <- function(object, ...)
  UseMethod("na.roughfix")

na.roughfix.data.frame <- function(object, ...) {
  ##  n <- length(object)
  ##  vars <- seq(length = n)
  isfac <- sapply(object, is.factor)
  isnum <- sapply(object, is.numeric)
  hasNA <- which(sapply(object, function(x) any(is.na(x))))
  if (any(!(isfac | isnum)))
    stop("na.roughfix only works for numeric or factor")
  for (j in hasNA) {
    if (isfac[j]) {
      freq <- table(object[[j]])
      xmode <- names(freq)[max.col(rbind(freq))]
      object[[j]][is.na(object[[j]])] <- xmode
    } else {
      xmed <- median(object[[j]], na.rm=TRUE)
      object[[j]][is.na(object[[j]])] <- xmed
    }
  }
  object
}

na.roughfix.default <- function(object, ...) {
  if (!is.atomic(object))
    return(object)
  d <- dim(object)
  if (length(d) > 2)
    stop("can't handle objects with more than two dimensions")
  if (all(!is.na(object)))
    return(object)
  if (!is.numeric(object))
    stop("roughfix can only deal with numeric data.")
  if (length(d) == 2) {
    hasNA <- which(apply(object, 2, function(x) any(is.na(x))))
    for (j in hasNA)
      object[is.na(object[, j]), j] <- median(object[, j], na.rm=TRUE)
  } else {
    object[is.na(object)] <- median(object, na.rm=TRUE)
  }
  object
}


outlier <- function(x, ...) UseMethod("outlier")

outlier.randomForest <- function(x, ...) {
  if (!inherits(x, "randomForest")) stop("x is not a randomForest object")
  if (x$type == "regression") stop("no outlier measure for regression")
  if (is.null(x$proximity)) stop("no proximity measures available")
  outlier.default(x$proximity, x$y)
}

outlier.default <- function(x, cls=NULL, ...) {
  if (nrow(x) != ncol(x)) stop ("x must be a square matrix")
  n <- nrow(x)
  if (is.null(cls)) cls <- rep(1, n)
  cls <- factor(cls)
  lvl <- levels(cls)
  cls.n <- table(cls)[lvl]
  id <- if (is.null(rownames(x))) 1:n else rownames(x)
  outlier <- structure(rep(NA, n), names=id)
  for (i in lvl) {
    out <- rowSums(x[cls == i, cls == i]^2)
    out <- n / ifelse(out == 0, 1, out)
    out <- (out - median(out)) / mad(out)
    outlier[names(out)] <- out
  }
  outlier
}

partialPlot <- function(x, ...) UseMethod("partialPlot")

partialPlot.default <- function(x, ...)
  stop("partial dependence plot not implemented for this class of objects.\n")

partialPlot.randomForest <-
  function (x, pred.data, x.var, which.class, w, plot=TRUE, add=FALSE,
            n.pt = min(length(unique(pred.data[, xname])), 51), rug = TRUE,
            xlab=deparse(substitute(x.var)), ylab="",
            main=paste("Partial Dependence on", deparse(substitute(x.var))),
            ...)
  {
    classRF <- x$type != "regression"
    if (is.null(x$forest))
      stop("The randomForest object must contain the forest.\n")
    x.var <- substitute(x.var)
    xname <- if (is.character(x.var)) x.var else {
      if (is.name(x.var)) deparse(x.var) else {
        eval(x.var)
      }
    }
    xv <- pred.data[, xname]
    n <- nrow(pred.data)
    if (missing(w)) w <- rep(1, n)
    if (classRF) {
      if (missing(which.class)) {
        focus <- 1
      }
      else {
        focus <- charmatch(which.class, colnames(x$votes))
        if (is.na(focus))
          stop(which.class, "is not one of the class labels.")
      }
    }
    if (is.factor(xv) && !is.ordered(xv)) {
      x.pt <- levels(xv)
      y.pt <- numeric(length(x.pt))
      for (i in seq(along = x.pt)) {
        x.data <- pred.data
        x.data[, xname] <- factor(rep(x.pt[i], n), levels = x.pt)
        if (classRF) {
          pr <- predict(x, x.data, type = "prob")
          y.pt[i] <- weighted.mean(log(ifelse(pr[, focus] > 0,
                                              pr[, focus], 1)) -
                                     rowMeans(log(ifelse(pr > 0, pr, 1))),
                                   w, na.rm=TRUE)
        } else y.pt[i] <- weighted.mean(predict(x, x.data), w, na.rm=TRUE)

      }
      if (add) {
        points(1:length(x.pt), y.pt, type="h", lwd=2, ...)
      } else {
        if (plot) barplot(y.pt, width=rep(1, length(y.pt)), col="blue",
                          xlab = xlab, ylab = ylab, main=main,
                          names.arg=x.pt, ...)
      }
    } else {
      if (is.ordered(xv))
        xv <- as.numeric(xv)
      x.pt <- seq(min(xv), max(xv), length = n.pt)
      y.pt <- numeric(length(x.pt))
      for (i in seq(along = x.pt)) {
        x.data <- pred.data
        x.data[, xname] <- rep(x.pt[i], n)
        if (classRF) {
          pr <- predict(x, x.data, type = "prob")
          y.pt[i] <- weighted.mean(log(ifelse(pr[, focus] == 0, 1, pr[, focus]))
                                   - rowMeans(log(ifelse(pr == 0, 1, pr))),
                                   w, na.rm=TRUE)
        } else {
          y.pt[i] <- weighted.mean(predict(x, x.data), w, na.rm=TRUE)
        }
      }
      if (add) {
        lines(x.pt, y.pt, ...)
      } else {
        if (plot) plot(x.pt, y.pt, type = "l", xlab=xlab, ylab=ylab,
                       main = main, ...)
      }
      if (rug && plot) {
        if (n.pt > 10) {
          rug(quantile(xv, seq(0.1, 0.9, by = 0.1)), side = 1)
        } else {
          rug(unique(xv, side = 1))
        }
      }
    }
    invisible(list(x = x.pt, y = y.pt))
  }


plot.randomForest <- function(x, type="l", main=deparse(substitute(x)), ...) {
  if(x$type == "unsupervised")
    stop("No plot for unsupervised randomForest.")
  test <- !(is.null(x$test$mse) || is.null(x$test$err.rate))
  if(x$type == "regression") {
    err <- x$mse
    if(test) err <- cbind(err, x$test$mse)
  } else {
    err <- x$err.rate
    if(test) err <- cbind(err, x$test$err.rate)
  }
  if(test) {
    colnames(err) <- c("OOB", "Test")
    matplot(1:x$ntree, err, type = type, xlab="trees", ylab="Error",
            main=main, ...)
  } else {
    matplot(1:x$ntree, err, type = type, xlab="trees", ylab="Error",
            main=main, ...)
  }
  invisible(err)
}

predict.randomForest <-
  function (object, newdata, type = "response", norm.votes = TRUE,
            predict.all=FALSE, proximity = FALSE, nodes=FALSE, cutoff, ...)
  {
    if (!inherits(object, "randomForest"))
      stop("object not of class randomForest")
    if (is.null(object$forest)) stop("No forest component in the object")
    out.type <- charmatch(tolower(type),
                          c("response", "prob", "vote", "class"))
    if (is.na(out.type))
      stop("type must be one of 'response', 'prob', 'vote'")
    if (out.type == 4) out.type <- 1
    if (out.type != 1 && object$type == "regression")
      stop("'prob' or 'vote' not meaningful for regression")
    if (out.type == 2)
      norm.votes <- TRUE
    if (missing(newdata)) {
      if (object$type == "regression") return(object$predicted)
      if (proximity & is.null(object$proximity))
        warning("cannot return proximity without new data if random forest object does not already have proximity")
      if (out.type == 1) {
        if (proximity) {
          return(list(pred = object$predicted,
                      proximity = object$proximity))
        } else return(object$predicted)
      }
      if (norm.votes) {
        t1 <- t(apply(object$votes, 1, function(x) { x/sum(x) }))
        if(proximity) return(list(pred = t1, proximity = object$proximity))
        else return(t1)
      } else {
        if (proximity) return(list(pred = object$votes, proximity = object$proximity))
        else return(object$votes)
      }
    }
    if (missing(cutoff)) {
      cutoff <- object$forest$cutoff
    } else {
      if (sum(cutoff) > 1 || sum(cutoff) < 0 || !all(cutoff > 0) ||
          length(cutoff) != length(object$classes)) {
        stop("Incorrect cutoff specified.")
      }
      if (!is.null(names(cutoff))) {
        if (!all(names(cutoff) %in% object$classes)) {
          stop("Wrong name(s) for cutoff")
        }
        cutoff <- cutoff[object$classes]
      }
    }

    if (object$type == "unsupervised")
      stop("Can't predict unsupervised forest.")

    if (inherits(object, "randomForest.formula")) {
      newdata <- as.data.frame(newdata)
      rn <- row.names(newdata)
      Terms <- delete.response(object$terms)
      x <- model.frame(Terms, newdata, na.action = na.omit)
      keep <- match(row.names(x), rn)
    } else {
      if (is.null(dim(newdata)))
        dim(newdata) <- c(1, length(newdata))
      x <- newdata
      if (nrow(x) == 0)
        stop("newdata has 0 rows")
      if (any(is.na(x)))
        stop("missing values in newdata")
      keep <- 1:nrow(x)
      rn <- rownames(x)
      if (is.null(rn)) rn <- keep
    }
    vname <- if (is.null(dim(object$importance))) {
      names(object$importance)
    } else {
      rownames(object$importance)
    }
    if (is.null(colnames(x))) {
      if (ncol(x) != length(vname)) {
        stop("number of variables in newdata does not match that in the training data")
      }
    } else {
      if (any(! vname %in% colnames(x)))
        stop("variables in the training data missing in newdata")
      x <- x[, vname, drop=FALSE]
    }
    if (is.data.frame(x)) {
      xfactor <- which(sapply(x, is.factor))
      if (length(xfactor) > 0 && "xlevels" %in% names(object$forest)) {
        for (i in xfactor) {
          if (any(! levels(x[[i]]) %in% object$forest$xlevels[[i]]))
            stop("New factor levels not present in the training data")
          x[[i]] <-
            factor(x[[i]],
                   levels=levels(x[[i]])[match(levels(x[[i]]), object$forest$xlevels[[i]])])
        }
      }
      cat.new <- sapply(x, function(x) if (is.factor(x) && !is.ordered(x))
        length(levels(x)) else 1)
      if (!all(object$forest$ncat == cat.new))
        stop("Type of predictors in new data do not match that of the training data.")
    }
    mdim <- ncol(x)
    ntest <- nrow(x)
    ntree <- object$forest$ntree
    maxcat <- max(object$forest$ncat)
    nclass <- object$forest$nclass
    nrnodes <- object$forest$nrnodes
    ## get rid of warning:
    op <- options(warn=-1)
    on.exit(options(op))
    x <- t(data.matrix(x))

    if (predict.all) {
      treepred <- if (object$type == "regression") {
        matrix(double(ntest * ntree), ncol=ntree)
      } else {
        matrix(integer(ntest * ntree), ncol=ntree)
      }
    } else {
      treepred <- numeric(ntest)
    }
    proxmatrix <- if (proximity) matrix(0, ntest, ntest) else numeric(1)
    nodexts <- if (nodes) integer(ntest * ntree) else integer(ntest)

    if (object$type == "regression") {
      if (!is.null(object$forest$treemap)) {
        object$forest$leftDaughter <-
          object$forest$treemap[,1,, drop=FALSE]
        object$forest$rightDaughter <-
          object$forest$treemap[,2,, drop=FALSE]
        object$forest$treemap <- NULL
      }

      keepIndex <- "ypred"
      if (predict.all) keepIndex <- c(keepIndex, "treepred")
      if (proximity) keepIndex <- c(keepIndex, "proximity")
      if (nodes) keepIndex <- c(keepIndex, "nodexts")
      ans <- .C("regForest",
                as.double(x),
                ypred = double(ntest),
                as.integer(mdim),
                as.integer(ntest),
                as.integer(ntree),
                as.integer(object$forest$leftDaughter),
                as.integer(object$forest$rightDaughter),
                as.integer(object$forest$nodestatus),
                as.integer(object$forest$nrnodes),
                as.double(object$forest$xbestsplit),
                as.double(object$forest$nodepred),
                as.integer(object$forest$bestvar),
                as.integer(object$forest$ndbigtree),
                as.integer(object$forest$ncat),
                as.integer(maxcat),
                as.integer(predict.all),
                treepred = as.double(treepred),
                as.integer(proximity),
                proximity = as.double(proxmatrix),
                nodes = as.integer(nodes),
                nodexts = as.integer(nodexts),
                DUP=FALSE,
                PACKAGE = "extendedForest")[keepIndex]
      ## Apply bias correction if needed.
      if (!is.null(object$coefs)) {
        yhat <- object$coefs[1] + object$coefs[2] * ans$ypred
      } else {
        yhat <- ans$ypred
      }
      if (predict.all) {
        treepred <- matrix(ans$treepred, length(keep),
                           dimnames=list(rn[keep], NULL))
      }
      if (!proximity) {
        res <- if (predict.all)
          list(aggregate=yhat, individual=treepred) else yhat
      } else {
        res <- list(predicted = yhat,
                    proximity = structure(ans$proximity,
                                          dim=c(ntest, ntest), dimnames=list(rn, rn)))
      }
      if (nodes) {
        attr(res, "nodes") <- matrix(ans$nodexts, ntest, ntree,
                                     dimnames=list(rn[keep], 1:ntree))
      }
    } else {
      countts <- matrix(0, ntest, nclass)
      t1 <- .C("classForest",
               mdim = as.integer(mdim),
               ntest = as.integer(ntest),
               nclass = as.integer(object$forest$nclass),
               maxcat = as.integer(maxcat),
               nrnodes = as.integer(nrnodes),
               jbt = as.integer(ntree),
               xts = as.double(x),
               xbestsplit = as.double(object$forest$xbestsplit),
               pid = as.double(object$forest$pid),
               cutoff = as.double(cutoff),
               countts = as.double(countts),
               treemap = as.integer(aperm(object$forest$treemap,
                                          c(2, 1, 3))),
               nodestatus = as.integer(object$forest$nodestatus),
               cat = as.integer(object$forest$ncat),
               nodepred = as.integer(object$forest$nodepred),
               treepred = as.integer(treepred),
               jet = as.integer(numeric(ntest)),
               bestvar = as.integer(object$forest$bestvar),
               nodexts = nodexts,
               ndbigtree = as.integer(object$forest$ndbigtree),
               predict.all = as.integer(predict.all),
               prox = as.integer(proximity),
               proxmatrix = as.double(proxmatrix),
               nodes = as.integer(nodes),
               DUP=TRUE,
               PACKAGE = "extendedForest")
      if (out.type > 1) {
        out.class.votes <- t(matrix(t1$countts, nrow = nclass, ncol = ntest))
        if (norm.votes)
          out.class.votes <-
            sweep(out.class.votes, 1, rowSums(out.class.votes), "/")
        z <- matrix(NA, length(rn), nclass,
                    dimnames=list(rn, object$classes))
        z[keep, ] <- out.class.votes
        res <- z
      } else {
        out.class <- factor(rep(NA, length(rn)),
                            levels=1:length(object$classes),
                            labels=object$classes)
        out.class[keep] <- object$classes[t1$jet]
        names(out.class[keep]) <- rn[keep]
        res <- out.class
      }
      if (predict.all) {
        treepred <- matrix(object$classes[t1$treepred],
                           nrow=length(keep), dimnames=list(rn[keep], NULL))
        res <- list(aggregate=res, individual=treepred)
      }
      if (proximity)
        res <- list(predicted = res, proximity = structure(t1$proxmatrix,
                                                           dim = c(ntest, ntest),
                                                           dimnames = list(rn[keep], rn[keep])))
      if (nodes) attr(res, "nodes") <- matrix(t1$nodexts, ntest, ntree,
                                              dimnames=list(rn[keep], 1:ntree))
    }
    res
  }



print.randomForest <-
  function(x, ...) {
    cat("\nCall:\n", deparse(x$call), "\n")
    cat("               Type of random forest: ", x$type, "\n", sep="")
    cat("                     Number of trees: ", x$ntree, "\n",sep="")
    cat("No. of variables tried at each split: ", x$mtry, "\n\n", sep="")
    if(x$type == "classification") {
      if(!is.null(x$confusion)) {
        cat("        OOB estimate of  error rate: ",
            round(x$err.rate[x$ntree, "OOB"]*100, dig=2), "%\n", sep="")
        cat("Confusion matrix:\n")
        print(x$confusion)
        if(!is.null(x$test$err.rate)) {
          cat("                Test set error rate: ",
              round(x$test$err.rate[x$ntree, "Test"]*100, dig=2), "%\n",
              sep="")
          cat("Confusion matrix:\n")
          print(x$test$confusion)
        }
      }
    }
    if(x$type == "regression") {
      if(!is.null(x$mse)) {
        cat("          Mean of squared residuals: ", x$mse[length(x$mse)],
            "\n", sep="")
        cat("                    % Var explained: ",
            round(100*x$rsq[length(x$rsq)], dig=2), "\n", sep="")
        if(!is.null(x$test$mse)) {
          cat("                       Test set MSE: ",
              round(x$test$mse[length(x$test$mse)], dig=2), "\n", sep="")
          cat("                    % Var explained: ",
              round(100*x$test$rsq[length(x$test$rsq)], dig=2), "\n", sep="")
        }
      }
      if (!is.null(x$coefs)) {
        cat("  Bias correction applied:\n")
        cat("  Intercept: ", x$coefs[1], "\n")
        cat("      Slope: ", x$coefs[2], "\n")
      }
    }
  }

varImpPlot <- function(x, sort=TRUE,
                       n.var=min(30, nrow(x$importance)),
                       type=NULL, class=NULL, scale=TRUE,
                       main=deparse(substitute(x)), ...) {
  if (!inherits(x, "randomForest"))
    stop("This function only works for objects of class `randomForest'")
  imp <- importance(x, class=class, scale=scale, type=type, ...)
  ## If there are more than two columns, just use the last two columns.
  if (ncol(imp) > 2) imp <- imp[, -(1:(ncol(imp) - 2))]
  nmeas <- ncol(imp)
  if (nmeas > 1) {
    op <- par(mfrow=c(1, 2), mar=c(4, 5, 4, 1), mgp=c(2, .8, 0),
              oma=c(0, 0, 2, 0), no.readonly=TRUE)
    on.exit(par(op))
  }
  for (i in 1:nmeas) {
    ord <- if (sort) rev(order(imp[,i],
                               decreasing=TRUE)[1:n.var]) else 1:n.var
    xmin <- if (colnames(imp)[i] %in%
                c("IncNodePurity", "MeanDecreaseGini")) 0 else min(imp[ord, i])
    dotchart(imp[ord,i], xlab=colnames(imp)[i], ylab="",
             main=if (nmeas == 1) main else NULL,
             xlim=c(xmin, max(imp[,i])), ...)
  }
  if (nmeas > 1) mtext(outer=TRUE, side=3, text=main, cex=1.2)
  invisible(imp)
}


varUsed <- function(x, by.tree=FALSE, count=TRUE) {
  if (!inherits(x, "randomForest"))
    stop(deparse(substitute(x)), "is not a randomForest object")
  if (is.null(x$forest))
    stop(deparse(substitute(x)), "does not contain forest")

  p <- length(x$forest$ncat)  # Total number of variables.
  if (count) {
    if (by.tree) {
      v <- apply(x$forest$bestvar, 2, function(x) {
        xx <- numeric(p)
        y <- table(x[x>0])
        xx[as.numeric(names(y))] <- y
        xx
      })
    } else {
      v <- numeric(p)
      vv <- table(x$forest$bestvar[x$forest$bestvar > 0])
      v[as.numeric(names(vv))] <- vv
    }
  } else {
    v <- apply(x$forest$bestvar, 2, function(x) sort(unique(x[x>0])))
    if(!by.tree) v <- sort(unique(unlist(v)))
  }
  v
}


.onAttach <- function(libname, pkgname) {
  RFver <- read.dcf(file=system.file("DESCRIPTION", package=pkgname),
                    fields="Version")
  message(paste(pkgname, RFver))
  message("Type rfNews() to see new features/changes/bug fixes.")
}




