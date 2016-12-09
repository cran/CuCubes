#' Build MultiDimensional Feature Selector from IGs
#'
#' @param IGs max conditional information gains
#' @param dimensions number of dimensions
#' @param divisions number of divisions
#' @param response_divisions number of response divisions (i.e. categories-1)
#' @param IG_bits input is in binary log (as opposed to natural log)
#' @param IG_doubled input is doubled (to follow the chi-squared distribution)
#' @param ignore_lowest number of variables with the lowest IG to ignore (ignored if computed)
#' @param variable_number number of irrelevant variables (ignored if computed)
#' @param calc_variable_number whether to compute the number of neglected and irrelevant variables
#' @param mode_1D "exp" - exponential distribution, "lin" - linear function of chi-squared, "raw" - raw chi-squared
#' @param min_variable_number minimum number of irrelevant variables
#' @param max_ignore_lowest maximum number of ignored variables
#' @param max_iterations maximum number of iterations in variable number calculation
#' @param acceptable_error acceptable error level for distribution parameter
#' @return MDFS (list-based S3 class object) with the following named elements:
#'   "IGs" is a vector of information gains (input copy)
#'   "order" is a vector of ordinal numbers (order of variables by decreasing score)
#'   "chi.squared" is a vector of chi-squared p-values
#'   "p.values" is a vector of eventual p-values
#'   "scores" is a list of two vectors FDR and FWER with FDR and FWER scores respectively
#'   "lo.sq.dev." is a vector of square deviations used to calculate the number of ignored variables
#'   "hi.sq.dev." is a vector of square deviations used to calculate the number of irrelevant variables
#'   "ign.lowest" is a number of ignored variables
#'   "var.number" is a number of irrelevant variables
#'   "dist.param." is an exponential distribution parameter or linear coefficient
#'   "err.param." is a square error of the parameter
#' @importFrom stats pchisq
#' @export
MDFS <- function(
 IGs,
 dimensions,
 divisions,
 response_divisions=1,
 IG_bits=TRUE,
 IG_doubled=FALSE,
 ignore_lowest=length(IGs)%/%10,
 variable_number=length(IGs),
 calc_variable_number=TRUE,
 mode_1D="exp",
 min_variable_number=variable_number%/%2,
 max_ignore_lowest=variable_number%/%3,
 max_iterations=20,
 acceptable_error=0.05) {
 #check the reasonability of input
 stopifnot(dimensions>0,divisions>0,response_divisions>0,
           length(IGs)>0,
           variable_number>0,variable_number<=length(IGs),
           ignore_lowest>=0,ignore_lowest<variable_number,
           min_variable_number<length(IGs),max_ignore_lowest>0)

 #bits to nats and the factor 2
 if (IG_bits) IGs<-log(2)*IGs
 if (!IG_doubled) IGs<-2*IGs

 #order variables
 ord<-order(IGs)
 ord_IGs<-IGs[ord]

 #degrees of freedom
 df<-response_divisions*divisions*(divisions+1)^(dimensions-1)
 n_vars<-length(IGs)

 #compute p-values
 pvalues<-pchisq(IGs,df,lower.tail=FALSE)
 ord_pvalues<-pvalues[ord]
 log_pvalues<-pchisq(IGs,df,log.p=TRUE)[ord]

 #calculate number of irrelevant variables and distribution parameter
 nv<-variable_number
 n0<-ignore_lowest
 n0p<-(-1)
 iteration<-1

 k<-1:n_vars

 if (dimensions>1 || mode_1D=="exp") {
  lk<-log(k)
  w<-(k-1)/(n_vars-k+1) #weights
  lp<-log_pvalues[k]
 } else if (mode_1D=="lin") {
  w<-1/k/(n_vars-k+1) #weights
  lp<-ord_pvalues[k]
  lk<-k
 } else {
  w<-1/k/(n_vars-k+1) #weights
  lp<-1-ord_pvalues[k]
  lk<-k
 }

 Skk<-cumsum(w*lk^2)
 Sk<-cumsum(w*lk)
 S<-cumsum(w)
 Spp<-cumsum(w*lp^2)
 Spk<-cumsum(w*lp*lk)
 Sp<-cumsum(w*lp)

 S0<-cumsum(ord_pvalues[k])

 #repeat until n0 and nv stabilize
 while(n0!=n0p && iteration<max_iterations) {
  iteration<-iteration+1

  skk<-Skk-Skk[n0]
  spk<-Spk-Spk[n0]
  spp<-Spp-Spp[n0]
  Nsk<-lk*(Sk-Sk[n0])
  NNs<-lk^2*(S-S[n0])
  Nsp<-lk*(Sp-Sp[n0])

  #to raw only
  Nspk<-lk*spk
  NNspp<-lk^2*spp

  if (dimensions==1 && mode_1D=="raw") {
   sqdev<- (skk-2*Nspk+NNspp)/spp/(k-n0)
  } else {
   sqdev<-((skk-2*Nsk+NNs)/spp-(spk-Nsp)^2/spp^2)/(k-n0)
  }

  if (calc_variable_number) {
   nv<-min_variable_number-1+which.min(sqdev[min_variable_number:length(sqdev)]) #optimum number of irrelevant variables
  }

  #calculate least square estimate of distribution parameter
  alpha<-(spk[nv]-Nsp[nv])/spp[nv]
  d_alpha<-sqrt(sqdev[nv])

  #calculate number of ignored variables that minimize square error; additional 0th value of n0 was needed
  if (dimensions>1 || mode_1D=="exp") {
   s0<-c(S0[nv],S0[nv]-S0)
   sqdev0<-(alpha-(nv-c(0,k)+1)/s0)^2
  } else {
   skk<-c(Skk[nv],Skk[nv]-Skk)
   spk<-c(Spk[nv],Spk[nv]-Spk)
   spp<-c(Spp[nv],Spp[nv]-Spp)
   Nsk<-lk[nv]*c(Sk[nv],Sk[nv]-Sk)
   NNs<-lk[nv]^2*c(S[nv],S[nv]-S)
   Nsp<-lk[nv]*c(Sp[nv],Sp[nv]-Sp)

   #to raw only
   Nspk<-lk[nv]*spk
   NNspp<-lk[nv]^2*spp

   if (mode_1D=="raw") {
    sqdev0<- (skk-2*Nspk+NNspp)/spp/(nv-c(0,k))
   } else {
    sqdev0<-((skk-2*Nsk+NNs)/spp-(spk-Nsp)^2/spp^2)/(nv-c(0,k))
   }

  }
  n0p<-n0
  if (calc_variable_number) { n0<-which.min(sqdev0[1:max_ignore_lowest]-1) }
 }
 if (calc_variable_number) {
  ignore_lowest<-n0
  variable_number<-nv
 }

 if (dimensions>1 || mode_1D=="exp") {
  var_pvalues<-1-exp(alpha*log_pvalues)
 } else if (mode_1D=="lin") {
  alpha<-alpha/(n0-nv)
  d_alpha<-d_alpha/(nv-n0)
  var_pvalues<-alpha*ord_pvalues
 } else {
  alpha<-1
  d_alpha<-d_alpha/(nv-n0)
  var_pvalues<-ord_pvalues
 }

 sqdev[1:ignore_lowest+1]<-NA
 sqdev0[variable_number:length(sqdev0)]<-NA

 if (iteration>=max_iterations) {warning("Possibly not convergent")}
 if (ignore_lowest>=max_ignore_lowest || variable_number<=min_variable_number) {warning("Border value reached for variable number")}
 if (d_alpha>acceptable_error*alpha) {warning("Too big error value for distribution parameter")}

 var_scores <- list(FDR = (n_vars * var_pvalues / (1 + n_vars - seq(var_pvalues)))[order(ord)],
                    FWER = (seq(var_pvalues) * var_pvalues)[order(ord)])

 result <- list(IGs = IGs,
                order = rev(ord),
                chi.squared = pvalues,
                p.values = var_pvalues[order(ord)],
                scores = var_scores,
                lo.sq.dev. = sqdev0[order(ord)],
                hi.sq.dev. = sqdev[order(ord)],
                ign.lowest = ignore_lowest,
                var.number = variable_number,
                dist.param. = alpha,
                err.param. = d_alpha)
 class(result) <- 'MDFS'

 return(result)
}

#' Find indices of relevant variables
#'
#' @param fs feature selector
#' @param ... arguments passed to methods
#' @return indices of important variables
#' @export
RelevantVariables <- function(fs, ...) {
  UseMethod('RelevantVariables')
}

#' Find indices of relevant variables from MDFS
#'
#' @param fs an MDFS object
#' @param level statistical significance level
#' @param score score to use
#' @param ... ignored
#' @return indices of relevant variables
#' @export
RelevantVariables.MDFS <- function(fs, level=0.05, score='FDR', ...) {
  stopifnot(any(score == c('FDR', 'FWER')))

  ord <- fs[['order']]
  ind <- seq(ord)
  scores <- fs[['scores']][[score]]
  count <- min(ind[scores[ord] > level]) - 1

  if (count == 0)
    return(numeric())
  else
    return(ord[1:count])
}

#' Plot MDFS details
#'
#' @param x an MDFS object
#' @param plots plots to plot (I for max IG, p for p-values, FDR, FWER)
#' @param ... ignored
#' @importFrom graphics plot
#' @export
plot.MDFS <- function(x, plots=c('I', 'p', 'FDR'), ...) {
  ord <- x[['order']]
  for (plt in plots) {
    switch(plt,
           I = plot(seq(ord), x[['IGs']][ord], xlab="index", ylab=expression('I'[max]*'(X)')),
           p = plot(x[['p.values']][ord], seq(ord)/351, xlab="exponential p-value", ylab="experimental p-value"),
           FDR = plot(seq(ord), x[['scores']][['FDR']][ord], xlab="index", ylab='FDR'),
           FWER = plot(seq(ord), x[['scores']][['FWER']][ord], xlab="index", ylab='FWER'),
           stop(paste('I don\'t now how to plot plot', plt)))
  }
}
