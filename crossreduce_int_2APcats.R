##########################################################
### crossreduce_int function for interaction analysis 
### To get the estimates by AP_cat (2 categories)
### Modified from the crossreduce() in the dlnm package
### Kai Chen. August 2017. Helmholtz Zentrum Muenchen 
##########################################################

crossreduce_int_2APcats <- function (basis, model = NULL, type = "overall", value = NULL, 
                             coef = NULL, vcov = NULL, model.link = NULL, at = NULL, from = NULL, 
                             to = NULL, by = NULL, lag, bylag = 1, cen = NULL, ci.level = 0.95) 
{
  if (all(class(basis) != "crossbasis")) {
    stop("the first argument must be an object of class 'crossbasis'")
  }
  name <- deparse(substitute(basis))
  attr <- attributes(basis)
  if (ncol(basis) == 1) 
    cond <- name
  if (is.null(model) && (is.null(coef) || is.null(vcov))) {
    stop("At least 'model' or 'coef'-'vcov' must be provided")
  }
  type <- match.arg(type, c("overall", "var", "lag"))
  if (type != "overall") {
    if (is.null(value)) 
      stop("'value' must be provided for type 'var' or 'lag'")
    else if (!is.numeric(value) || length(value) > 1) {
      stop("'value' must be a numeric scalar")
    }
    if (type == "lag" && (any(value < attr$lag[1]) || any(value > 
                                                          attr$lag[2]))) {
      stop("'value' of lag-specific effects must be within the lag range")
    }
  }
  else value <- NULL
  lag <- if (missing(lag)) 
    attr$lag
  else mklag(lag)
  if (lag != attr$lag && attr$arglag$fun == "integer") 
    stop("prediction for lag sub-period not allowed for type 'integer'")
  if (!is.numeric(ci.level) || ci.level >= 1 || ci.level <= 
      0) {
    stop("'ci.level' must be numeric and between 0 and 1")
  }
  cond <- if (ncol(basis) == 1L) 
    name
  else paste(name, "[[:print:]]*v[0-9]{1,2}\\.l[0-9]{1,2}", 
             sep = "")
  cond.AP_cat2 <-  paste(name, "[[:print:]]*v[0-9]{1,2}\\.l[0-9]{1,2}\\:AP_cat2", 
                         sep = "")
  
  if (!is.null(model)) {
    model.class <- class(model)
    ##For ozone, Error occurs due to singluarity for cb and AP_cat, some coefs for cb are NA
    # coef <- summary(model)$coefficients[,1] ##change from 
    coef <- dlnm:::getcoef(model, model.class)
    ind.all <- grep(cond, names(coef))
    ind.AP_cat2 <- grep(cond.AP_cat2, names(coef))
    ind.main <- ind.all[ind.all != ind.AP_cat2]
    
    ###Extract the coef and vcov from the interaction model for AP categories
    coef.main <- coef[ind.main]
    coef.int_cat2 <- coef[ind.AP_cat2]
    ##vcov for AP categories
    vcov.all <- dlnm:::getvcov(model, model.class)
    vcov.main <- dlnm:::getvcov(model, model.class)[ind.main, ind.main, drop = FALSE]
    vcov.int_cat2 <- dlnm:::getvcov(model, model.class)[ind.AP_cat2, ind.AP_cat2, drop = FALSE]

    #cat=1
    coef_APcat1 <- coef.main
    vcov_APcat1 <- vcov.main
    #cat=2
    coef_APcat2 <- coef_APcat1+coef.int_cat2
    ####Important!! note that for interaction analysis, vcov(b1*b2)=var(b1)+var(b2)+2cov(b1,b2)
    ####This is only correct for cov(b1, b2) == cov(b2, b1); otherwise(like here), must using cov(b1, b2) + cov(b2,b1)
    vcov_APcat2 <- vcov_APcat1+vcov.int_cat2+dlnm:::getvcov(model, model.class)[ind.main, ind.AP_cat2, drop=FALSE]+
      dlnm:::getvcov(model, model.class)[ind.AP_cat2, ind.main, drop=FALSE]
     
    #model.link
    model.link <- dlnm:::getlink(model, model.class)
  }
  else model.class <- NA
  npar <- ncol(basis)
  range <- attr$range
  at <- dlnm:::mkat(at, from, to, by, range, lag, bylag)
  cen <- dlnm:::mkcen(cen, type = "cb", basis, range)
  attributes(basis)$argvar$cen <- attr$argvar$cen <- NULL
  if (type == "overall") {
    lagbasis <- do.call("onebasis", c(list(x = dlnm:::seqlag(lag)), 
                                      attr$arglag))
    M <- diag(ncol(basis)/ncol(lagbasis)) %x% (t(rep(1, diff(lag) + 
                                                       1)) %*% lagbasis)
    newbasis <- do.call("onebasis", c(list(x = at), attr$argvar))
    if (!is.null(cen)) {
      basiscen <- do.call("onebasis", c(list(x = cen), 
                                        attr$argvar))
      newbasis <- scale(newbasis, center = basiscen, scale = FALSE)
    }
  }
  else if (type == "lag") {
    lagbasis <- do.call("onebasis", c(list(x = value), attr$arglag))
    M <- diag(ncol(basis)/ncol(lagbasis)) %x% lagbasis
    newbasis <- do.call("onebasis", c(list(x = at), attr$argvar))
    if (!is.null(cen)) {
      basiscen <- do.call("onebasis", c(list(x = cen), 
                                        attr$argvar))
      newbasis <- scale(newbasis, center = basiscen, scale = FALSE)
    }
  }
  else if (type == "var") {
    varbasis <- do.call("onebasis", c(list(x = value), attr$argvar))
    if (!is.null(cen)) {
      basiscen <- do.call("onebasis", c(list(x = cen), 
                                        attr$argvar))
      varbasis <- scale(varbasis, center = basiscen, scale = FALSE)
    }
    M <- varbasis %x% diag(ncol(basis)/ncol(varbasis))
    newbasis <- do.call("onebasis", c(list(x = seqlag(lag, 
                                                      bylag)), attr$arglag))
  }
  dimnames(newbasis) <- list(seq(nrow(newbasis)), paste0("b", 
                                                         seq(ncol(newbasis))))
  ##cat=1
  newcoef_APcat1 <- as.vector(M %*% coef_APcat1)
  names(newcoef_APcat1) <- colnames(newbasis)
  newvcov_APcat1 <- M %*% vcov_APcat1 %*% t(M)
  dimnames(newvcov_APcat1) <- list(colnames(newbasis), colnames(newbasis))
  fit_APcat1 <- as.vector(newbasis %*% newcoef_APcat1)
  se_APcat1 <- sqrt(pmax(0, rowSums((newbasis %*% newvcov_APcat1) * newbasis)))
  
  if (type == "var") {
    names(fit_APcat1) <- names(se_APcat1) <- outer("lag", seqlag(lag, bylag), 
                                                   paste, sep = "")
  }
  else names(fit_APcat1) <- names(se_APcat1) <- at
  
  ##cat=2
  newcoef_APcat2 <- as.vector(M %*% coef_APcat2)
  names(newcoef_APcat2) <- colnames(newbasis)
  newvcov_APcat2 <- M %*% vcov_APcat2 %*% t(M)
  dimnames(newvcov_APcat2) <- list(colnames(newbasis), colnames(newbasis))
  fit_APcat2 <- as.vector(newbasis %*% newcoef_APcat2)
  se_APcat2 <- sqrt(pmax(0, rowSums((newbasis %*% newvcov_APcat2) * newbasis)))
  
  if (type == "var") {
    names(fit_APcat2) <- names(se_APcat2) <- outer("lag", seqlag(lag, bylag), 
                                                   paste, sep = "")
  }
  else names(fit_APcat2) <- names(se_APcat2) <- at
  
  
  ##result list
  list <- list(coef_APcat1 = newcoef_APcat1, vcov_APcat1 = newvcov_APcat1, 
               coef_APcat2 = newcoef_APcat2, vcov_APcat2 = newvcov_APcat2, 
               basis = newbasis, type = type, value = value)
  if (type != "var") 
    list$predvar <- at
  if (!is.null(cen)) 
    list$cen <- cen
  list <- c(list, list(lag = lag, bylag = bylag, fit_APcat1 = fit_APcat1, se_APcat1 = se_APcat1,
                       fit_APcat2 = fit_APcat2, se_APcat2 = se_APcat2))
  
  z <- qnorm(1 - (1 - ci.level)/2)
  if (model.link %in% c("log", "logit")) {
    #cat=1
    list$RRfit_APcat1 <- exp(fit_APcat1)
    list$RRlow_APcat1 <- exp(fit_APcat1 - z * se_APcat1)
    names(list$RRlow_APcat1) <- names(fit_APcat1)
    list$RRhigh_APcat1 <- exp(fit_APcat1 + z * se_APcat1)
    names(list$RRhigh_APcat1) <- names(fit_APcat1)
    #cat=2
    list$RRfit_APcat2 <- exp(fit_APcat2)
    list$RRlow_APcat2 <- exp(fit_APcat2 - z * se_APcat2)
    names(list$RRlow_APcat2) <- names(fit_APcat2)
    list$RRhigh_APcat2 <- exp(fit_APcat2 + z * se_APcat2)
    names(list$RRhigh_APcat2) <- names(fit_APcat2)
  }
  else {
    #cat1
    list$low_APcat1 <- fit_APcat1 - z * se_APcat1
    names(list$low_APcat1) <- names(fit_APcat1)
    list$high_APcat1 <- fit_APcat1 + z * se_APcat1
    names(list$high_APcat1) <- names(fit_APcat1)
    #cat2
    list$low_APcat2 <- fit_APcat2 - z * se_APcat2
    names(list$low_APcat2) <- names(fit_APcat2)
    list$high_APcat2 <- fit_APcat2 + z * se_APcat2
    names(list$high_APcat2) <- names(fit_APcat2)
     }
  list$ci.level <- ci.level
  list$model.class <- model.class
  list$model.link <- model.link
  class(list) <- "crossreduce"
  return(list)
}
