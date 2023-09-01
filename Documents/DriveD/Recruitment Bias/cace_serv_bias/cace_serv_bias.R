# Call library
library("lme4")       # For multilevel models
library("estimatr")

# Define function
cace_serv_bias <- function(data_csv,serv_var,y_var,trt_var,xtlogit,xc_logit,
                           x_wls,clus_var,wgt,marg,se_ipw,df_wls,out_regr,out_est) {

# Read in csv data

mydat1 <- read.csv(data_csv)

sum_dat1  <- summary(mydat1)
nobs_dat1 <- NROW(mydat1)

# Define defaults

if ((clus_var == '0') | (clus_var == "")) {
    clus_var <- c(0)
    }

if ((wgt == '0') | (wgt == "")) {
    wgt <- c(0)
    }

if ((marg == '1') | (marg == "")) {
    marg <- c(1)
    } else if (marg == '0') {
      marg <- c(0)
      }

if ((se_ipw == '1') | (se_ipw == "")) {
  se_ipw <- c(1)
  } else if (se_ipw == '0') {
    se_ipw <- c(0)
    }

if ((df_wls == '1') | (df_wls == "")) {
  df_wls <- c(1)
} else if (df_wls == '0') {
  df_wls <- c(0)
}

if (out_regr == "") {
  out_regr  <- c('cace regr log.txt')
  }

if (out_est == "") {
  out_est  <- c('cace estimates.txt')
  }

# If there is no clus_var than must use the marginal model

if ((clus_var == 0) | (clus_var == "")) {
  marg <- c(1)
  }

# Add the weight onto mydat

if (wgt == 0) {
  mydat1['ww'] <- 1
  } else {
      mydat1['ww'] <- mydat1[wgt]
      }

# Add on service variable and treatment variable

mydat1['serv_var']  <- mydat1[serv_var]
mydat1['trt_var']   <- mydat1[trt_var]

# Add on cluster variable, which is set to 1 for non-clustered designs for now

if (clus_var == 0) {
  mydat1['clusv'] <- 1
  } else {
    mydat1['clusv'] <- mydat1[clus_var]
    }

mobs_dat1 <- NROW(unique(mydat1$clusv))

# Keep only cases with treatment and service receipt = 0,1 and >0 wgts

mydat1 <- mydat1[which(((mydat1$trt_var==0)  | (mydat1$trt_var==1))),]
mydat1 <- mydat1[which(((mydat1$serv_var==0) | (mydat1$serv_var==1))),]
mydat1 <- mydat1[which(mydat1$ww > 0),]

# Check whether x_wls has elements

hasx <- 0
if (any(is.element(x_wls,'0')) == FALSE | any(is.element(x_wls,0)) == FALSE) {
  hasx <- 1
  }

#
# Conduct a complete case analysis so add on other X and Y variables
#

compl  <- c('serv_var','trt_var','ww','clusv')

# Add on xt_logit variables

compla <- c("")
for (i in 1:NROW(xt_logit)) {
  vname <- paste0("zzxt", i)
  mydat1[[vname]] <- as.numeric(unlist(mydat1[xt_logit[[i]]]))

  if (i == 1) {
    compla <- c('zzxt1')
    } else {
      compla <- c(compla,paste0("zzxt", i))
      }
    }

compl <- c(compl,compla)

# Add on xc_logit variables

compla <- c("")
for (i in 1:NROW(xc_logit)) {
  vname <- paste0("zzxc", i)
  mydat1[[vname]] <- as.numeric(unlist(mydat1[xc_logit[[i]]]))

  if (i == 1) {
    compla <- c('zzxc1')
    } else {
      compla <- c(compla,paste0("zzxc", i))
      }
    }

compl <- c(compl,compla)

# Add x_wls if specified

if (hasx == 1) {

  compla <- c("")
  for (i in 1:NROW(x_wls)) {
    vname <- paste0("zzx", i)
    mydat1[[vname]] <- as.numeric(unlist(mydat1[x_wls[[i]]]))

    if (i == 1) {
      compla <- c('zzx1')
      } else {
        compla <- c(compla,paste0("zzx", i))
        }
      }

  compl <- c(compl,compla)
  }

# DO OVER Y

for (iy in 1:NROW(y_var)) {

yname <- y_var[[iy]]
vname <- paste0("zzy", iy)
mydat1[[vname]] <- as.numeric(unlist(mydat1[y_var[[iy]]]))

compla <- c(paste0("zzy", iy))

complg <- c(compl,compla)

###
# Take complete cases only
###

mydat <- mydat1[complete.cases(mydat1[ , complg]), ]

nobs_dat <- NROW(mydat)
mobs_dat <- NROW(unique(mydat$clusv))

# Define clusv for non-clustered designs to row number

if (clus_var == 0) {
  mydat['clusv'] <- seq_len(NROW(mydat['serv_var']))
  }

# Normalize the ww weights to nt and nc

mydat['ttt'] <- mydat[trt_var]

normt <- sum(mydat$ww[mydat$ttt == 1])
n_wwt <- NROW(mydat$ww[mydat$ttt == 1])

mydat$ww[mydat$ttt == 1] <- mydat$ww[mydat$ttt == 1]*n_wwt/normt

normc <- sum(mydat$ww[mydat$ttt == 0])
n_wwc <- NROW(mydat$ww[mydat$ttt == 0])

mydat$ww[mydat$ttt == 0] <- mydat$ww[mydat$ttt == 0]*n_wwc/normc

# Create X covariates for logit regressions

xtdat <- mydat[xt_logit]
xcdat <- mydat[xc_logit]

# Check whether x_wls has elements and create X covariates for WLS regression

if (hasx == 1) {
  xdat <- mydat[x_wls]
  }

###
# RUN LOGIT MODEL FOR THE TREATMENT GROUP
###

# Create ppt variable for treatment logit model
mydat$ppt <- mydat$serv_var
mydat$ppt[mydat[trt_var]==0] <- NA

# Create logit data for treatment model
tldat <- cbind(mydat['ppt'],xtdat)

# Add on weight variable
tldat <- cbind(tldat,mydat['ww'])

# Add on cluster variable
tldat_clus   <- cbind(tldat,mydat['clusv'])

# Create formulas for treatment logit model

formt <- as.formula(paste("ppt ~ ", paste(colnames(xtdat), collapse = " + ")))
formt_clus <- as.formula(paste("ppt ~ ", paste(colnames(xtdat), collapse = " + "),"+ (1|clusv)"))

# Fit logistic regression model for the treatment group

if (marg == 1) {
  logt <- suppressWarnings(glm(formt,
              family="binomial",
              weights=ww,
              data=tldat))

  # Generate predicted probs for everyone since logt$fitted.values
  # is only for those with the ppt dependent variable (treatments)

  predt <- predict(logt, newdata = tldat, type = 'response')

  } else if ((marg == 0) | (marg == '0')) {
    logt <- suppressWarnings(glmer(formt_clus,
                  family="binomial",
                  weights=ww,
                  nAGQ=10,
                  data=tldat_clus))
    predt <- predict(logt, newdata = tldat, re.form=NA, type = 'response')
    }

#view model summary
logt_res <- summary(logt)

# Add predicted probs back to the main dataset and calculate residuals
mydat$predt <- predt
mydat$resid_lt <- mydat[serv_var] - mydat$predt

# NOW DO THE SAME FOR THE CONTROL GROUP

# Create ppc variable for control logit model
mydat$ppc <- mydat$serv_var
mydat$ppc[mydat[trt_var]==1] <- NA

# Create logit data for control model
cldat      <- cbind(mydat['ppc'],xcdat)
cldat      <- cbind(cldat,mydat['ww'])
cldat_clus <- cbind(cldat,mydat['clusv'])

# Create formulas for treatment model

formc <- as.formula(paste("ppc ~ ", paste(colnames(xcdat), collapse = " + ")))
formc_clus <- as.formula(paste("ppc ~ ", paste(colnames(xcdat), collapse = " + "),"+ (1|clusv)"))

# Fit the logistic regression model for the control group

if (marg == 1) {
  logc <- suppressWarnings(glm(formc,
              family="binomial",
              weights=ww,
              data=cldat))
  predc <- predict(logc, newdata = cldat, type = 'response')

  } else if ((marg == 0) | (marg == '0')) {
    logc <- suppressWarnings(glmer(formc_clus,
                  family="binomial",
                  weights=ww,
                  nAGQ=10,
                  data=cldat_clus))
    predc <- predict(logc, newdata = cldat, re.form=NA, type = 'response')
    }

#view model summary
logc_res <- summary(logc)

# Add predicted probs back to the main dataset and calculate residuals
mydat$predc    <- predc
mydat$resid_lc <- mydat[serv_var] - mydat$predc

# Construct weights

rrr <- mydat[serv_var]
ttt <- mydat[trt_var]

mydat['wgt_tc'] <- rrr  + (1-rrr)*(ttt*predc +(1-ttt)*predt)
mydat['wgt_t']  <- ttt*rrr + (1-ttt)*predt
mydat['wgt_11'] <- rrr*(ttt*predc + (1-ttt)*predt)

# Multiply weight by sample weight but first normalize ww to sum to

mydat$wgt_tc <- mydat$wgt_tc*mydat$ww
mydat$wgt_t  <- mydat$wgt_t*mydat$ww
mydat$wgt_11 <- mydat$wgt_11*mydat$ww

# Run WLS regressions

# Create formula

if (hasx == 1) {
  form_wls <- as.formula(paste("yyy ~ ttt +", paste(colnames(xcdat), collapse = " + ")))
  } else if (hasx == 0) {
    form_wls <- as.formula(paste("yyy ~ ttt"))
    }

# Create WLS data

mydat['yyy'] <- mydat[yname]
#mydat['ttt'] <- mydat[trt_var]
wlsdat <- cbind(mydat['yyy'],mydat['ttt'])
wlsdat <- cbind(wlsdat,mydat['ww'])
wlsdat <- cbind(wlsdat,mydat['clusv'])

if (hasx == 1) {
  wlsdat <- cbind(wlsdat,xdat)
  }

wlsdat_tc <- cbind(wlsdat,mydat['wgt_tc'])
wlsdat_t  <- cbind(wlsdat,mydat['wgt_t'])
wlsdat_11 <- cbind(wlsdat,mydat['wgt_11'])

# Calculate prop estimates - first for Ts and then for Cs

# CACE_TC for IV
swgt_tc_t <- sum(wlsdat_tc$wgt_tc[wlsdat_tc$ttt == 1])
sww_t     <- sum(wlsdat_tc$ww[wlsdat_tc$ttt == 1])

prop_tc_t <- swgt_tc_t/sww_t

swgt_tc_c <- sum(wlsdat_tc$wgt_tc[wlsdat_tc$ttt == 0])
sww_c     <- sum(wlsdat_tc$ww[wlsdat_tc$ttt == 0])

prop_tc_c <- swgt_tc_c/sww_c

# For CACE_11
swgt_11_t <- sum(wlsdat_11$wgt_11[wlsdat_11$ttt == 1])
sww_11_t  <- sum(wlsdat_11$ww[wlsdat_11$ttt == 1])

prop_11_t <- swgt_11_t/sww_11_t

swgt_11_c <- sum(wlsdat_11$wgt_11[wlsdat_11$ttt == 0])
sww_11_c  <- sum(wlsdat_11$ww[wlsdat_11$ttt == 0])

prop_11_c <- swgt_11_c/sww_11_c

# For CACE_T can use service receipt proportion in treatment group

sww_t_t   <- sum(wlsdat_t$ww[wlsdat_t$ttt == 1])
serv_wt   <- mydat$serv_var[mydat$ttt == 1]*mydat$ww[mydat$ttt == 1]
prop_t_t  <- sum(serv_wt)/sww_t_t

swgt_t_c  <- sum(wlsdat_t$wgt_t[wlsdat_t$ttt == 0])
sww_t_c   <- sum(wlsdat_t$ww[wlsdat_t$ttt == 0])

prop_t_c  <- swgt_t_c/sww_t_c

# Use the treatment group value for the IV estimator

prop_tc   <- prop_tc_t

###
# RUN REGRESSIONS
###

###
# DO THE REGRESSIONS IF DO NOT WANT ADJ IPW STANDARD ERRORS - USE LM_ROBUST EXISTING PACKAGES
###

if (se_ipw == 0) {
  if (clus_var == 0) {
    # CACE-TC Non-clustered

    sand_tc   <- lm_robust(form_wls, data = wlsdat_tc, weights = wgt_tc, se_type = "HC0")

    cace_tc_res    <- summary(sand_tc)
    imp_tc         <- sand_tc$coefficients[2]
    sand_se_all_tc <- sand_tc$std.error[2]
    tstat_tc       <- sand_tc$statistic[2]
    p_tc           <- sand_tc$p.value[2]

    # CACE-T

    sand_t   <- lm_robust(form_wls, data = wlsdat_t, weights = wgt_t, se_type = "HC0")

    cace_t_res    <- summary(sand_t)
    imp_t         <- sand_t$coefficients[2]
    sand_se_all_t <- sand_t$std.error[2]
    tstat_t       <- sand_t$statistic[2]
    p_t           <- sand_t$p.value[2]

    # CACE-11 Non-clustered

    sand_11   <- lm_robust(form_wls, data = wlsdat_11, weights = wgt_11, se_type = "HC0")

    cace_11_res    <- summary(sand_11)
    imp_11         <- sand_11$coefficients[2]
    sand_se_all_11 <- sand_11$std.error[2]
    tstat_11       <- sand_11$statistic[2]
    p_11           <- sand_11$p.value[2]

    # ITT Non-clustered

    sand_itt  <- lm_robust(form_wls, data = wlsdat, weights = ww, se_type = "HC0")

    itt_res     <- summary(sand_itt)
    imp_itt     <- sand_itt$coefficients[2]
    sand_se_iv  <- sand_itt$std.error[2]
    tstat_iv    <- sand_itt$statistic[2]
    p_iv        <- sand_itt$p.value[2]

    sand_se_all_iv_tau <- sand_se_iv / prop_tc

    } else if (clus_var != 0) {

      # CACE-TC Clustered

      sand_tc   <- lm_robust(form_wls, data = wlsdat_tc, weights = wgt_tc, clusters = clusv,
                           se_type = "stata")

      kk <- 0

      if (hasx == 1) {
        kk         <- length(sand_tc$coefficients)        # Number of coefficients
      }

      cace_tc_res    <- summary(sand_tc)
      imp_tc         <- sand_tc$coefficients[2]
      sand_se_all_tc <- sand_tc$std.error[2]
      tstat_tc       <- sand_tc$statistic[2]
      p_tc           <- sand_tc$p.value[2]
      nobs_tc        <- sand_tc$nobs
      mclus_tc       <- sand_tc$nclusters

      seadj_tc       <- sqrt(((mclus_tc-1)/mclus_tc)*((nobs_tc-kk)/(nobs_tc-1)))

      sand_se_all_tc <- sand_se_all_tc*seadj_tc

      # CACE-T Clustered

      sand_t    <- lm_robust(form_wls, data = wlsdat_t,  weights = wgt_t, clusters = clusv,
                           se_type = "stata")

      cace_t_res    <- summary(sand_t)
      imp_t         <- sand_t$coefficients[2]
      sand_se_all_t <- sand_t$std.error[2]
      tstat_t       <- sand_t$statistic[2]
      p_t           <- sand_t$p.value[2]
      nobs_t        <- sand_t$nobs
      mclus_t       <- sand_t$nclusters

      seadj_t       <- sqrt(((mclus_t-1)/mclus_t)*((nobs_t-kk)/(nobs_t-1)))

      sand_se_all_t <- sand_se_all_t*seadj_t

      # CACE-11 Clustered

      sand_11   <- lm_robust(form_wls, data = wlsdat_11, weights = wgt_11, clusters = clusv,
                             se_type = "stata")
      cace_11_res    <- summary(sand_11)
      imp_11         <- sand_11$coefficients[2]
      sand_se_all_11 <- sand_11$std.error[2]
      tstat_11       <- sand_11$statistic[2]
      p_11           <- sand_11$p.value[2]
      nobs_11        <- sand_11$nobs
      mclus_11       <- sand_11$nclusters

      seadj_11       <- sqrt(((mclus_11-1)/mclus_11)*((nobs_11-kk)/(nobs_11-1)))

      sand_se_all_11 <- sand_se_all_11*seadj_11

      # ITT Clustered

      sand_itt  <- lm_robust(form_wls, data = wlsdat, weights = ww, clusters = clusv,
                             se_type = "stata")

      itt_res     <- summary(sand_itt)
      imp_itt     <- sand_itt$coefficients[2]
      sand_se_iv  <- sand_itt$std.error[2]
      tstat_iv    <- sand_itt$statistic[2]
      p_iv        <- sand_itt$p.value[2]
      nobs_itt    <- sand_itt$nobs
      mclus_itt   <- sand_itt$nclusters

      seadj_iv       <- sqrt(((mclus_itt-1)/mclus_itt)*((nobs_itt-kk)/(nobs_itt-1)))

      sand_se_all_iv_tau <- sand_se_iv*seadj_iv / prop_tc

      } # end if (clus_var != 0)

  # Define IMP_IV and Use the treatment group values

  imp_iv_t  <- imp_itt / prop_tc_t
  imp_iv_c  <- imp_itt / prop_tc_c

  imp_iv    <- imp_iv_t

  } # end if (se_ipw == 0)

###
# DO THE REGRESSIONS IF WANT ADJ IPW STANDARD ERRORS
###

if (se_ipw == 1) {

# CACE_TC

# Run wls_tc regression
cace_tc <- lm(form_wls, data = wlsdat_tc, weights = wgt_tc)

cace_tc_res <- summary(cace_tc)
imp_tc      <- cace_tc$coefficients[2]
resid_tc    <- cace_tc$residuals

# CACE_T

# Run wls_t regression
cace_t <- lm(form_wls, data=wlsdat_t, weights = wgt_t)

cace_t_res <- summary(cace_t)
imp_t      <- cace_t$coefficients[2]
resid_t    <- cace_t$residuals

# CACE_11

# Run wls_11 regression

cace_11 <- lm(form_wls, data=wlsdat_11, weights = wgt_11)

cace_11_res <- summary(cace_11)
imp_11      <- cace_11$coefficients[2]
resid_11    <- cace_11$residuals

# CACE_IV uses ITT

# Run wls_itt regression
itt <- lm(form_wls, data=wlsdat, weights = ww)

itt_res   <- summary(itt)
imp_itt   <- itt$coefficients[2]
resid_itt <- itt$residuals

# Calculate x_mean and x_mean times beta_hat for IV estimator by T

nx <- 0
if (hasx == 1) {
  ncoeff <- length(itt$coefficients)     # Number of coefficients
  betax  <- itt$coefficients[3:ncoeff]   # Pull off x coefficients
  nx     <- NROW(betax)                  # Number of covariates
  nobs   <- length(resid_itt)            # Number of observations

  # Make xdat a matrix since need to perform a matrix operation
  xdatm <- matrix(unlist(xdat), ncol = nx, nrow = nobs)

  xbeta <- xdatm %*% betax               # Calculate Xb
  wlsdat$xbeta <- xbeta                  # Add Xb to wlsdat

  # Calculate weighted means of Xb for treatments and controls
  wlsdat_trt <- wlsdat[which(wlsdat$ttt == 1),]
  wlsdat_cnt <- wlsdat[which(wlsdat$ttt == 0),]

  xbt <- weighted.mean(wlsdat_trt$xbeta,wlsdat_trt$ww)
  xbc <- weighted.mean(wlsdat_cnt$xbeta,wlsdat_cnt$ww)

  # Now calculate weighted means for each X by treatment and control
  mxt <- betax
  mxc <- betax
  for (i in 1:nx) {
    vv <- x_wls[i]                                           # Pull off the x variable
    mxt[i] <- weighted.mean(wlsdat_trt[[vv]],wlsdat_trt$ww)  # Weighted mean for Ts (notice double brackets!!)
    mxc[i] <- weighted.mean(wlsdat_cnt[[vv]],wlsdat_cnt$ww)  # Weighted mean for Cs (notice double brackets!!)
    }
  }  # if hasx == 1

# Define IMP_IV and Use the treatment group values

imp_iv_t  <- imp_itt / prop_tc_t
imp_iv_c  <- imp_itt / prop_tc_c

imp_iv    <- imp_iv_t

###
# VARIANCE CALCULATIONS FOR ADJ ESTIMATORS
###

mydat$ccc <- 1 - mydat$ttt             # Create control group indicator
if (hasx == 1) {
  xdat1 <- cbind(mydat['clusv'],xdat)   # Add clus to xdat
  }

xtdat1  <- cbind(mydat['clusv'],xtdat)  # Add clus to xtdat
xcdat1  <- cbind(mydat['clusv'],xcdat)  # Add clus to xcdat

# Create u vectors with WLS residuals and clusv and same for weights
u_tc    <- cbind(mydat['clusv'],resid_tc)
wgt_tc  <- cbind(mydat['clusv'],mydat['wgt_tc'])

u_t     <- cbind(mydat['clusv'],resid_t)
wgt_t   <- cbind(mydat['clusv'],mydat['wgt_t'])

u_11    <- cbind(mydat['clusv'],resid_11)
wgt_11  <- cbind(mydat['clusv'],mydat['wgt_11'])

u_iv    <- cbind(mydat['clusv'],resid_itt)
wgt_ww  <- cbind(mydat['clusv'],mydat['ww'])

k1 <- NROW(xt_logit)              # Number of covariates in T logit model
k0 <- NROW(xc_logit)              # Number of covariates in C logit model

mclus <- length(unique(mydat$clusv))  # Number of clusters

if (nx>0) {
  xb_diff  <- xbt - xbc
  mx_diff <- mxt - mxc
  } else {
    xb_diff <- 0
    }

m1 <- 0
m0 <- 0

# CACE-TC INITIALIZE VARIABLES FOR VARIANCES FOR CACE-TC IPW ESTIMATOR

# CACE-TC Variables to calculate WLS GEE

phij   <- matrix(0,(nx+2),1)
gammaj <- matrix(0,(nx+2),(nx+2))

delta <- gammaj
gamma <- gammaj

# CACE-TC Variables to calculate Logit Treatment and Control GEE

phij_l1   <- matrix(0,(k1+1),1)
gammaj_l1 <- matrix(0,(k1+1),(k1+1))

delta_l1 <- gammaj_l1
gamma_l1 <- gammaj_l1

phij_l0   <- matrix(0,(k0+1),1)
gammaj_l0 <- matrix(0,(k0+1),(k0+1))

delta_l0 <- gammaj_l0
gamma_l0 <- gammaj_l0

gammaj_l <- matrix(0,(k1+k0+2),(k1+k0+2))

# CACE-TC Variables to calculate Off-Diagonal WLS and Logit GEE - the B Matrix in Paper

gammaj_b <- matrix(0,(nx+2),(k1+k0+2))
gamma_b  <- gammaj_b

# CACE-TC Variables for all pieces

nv <- nx + k1 + k0 + 4

phij_all <- matrix(0,nv,1)
phi_all  <- phij_all

gammaj_all <- matrix(0,nv,nv)
gamma_all  <- gammaj_all

delta_all <- gammaj_all

# CACE-T Variables to calculate WLS GEE

phij_t   <- matrix(0,(nx+2),1)
gammaj_t <- matrix(0,(nx+2),(nx+2))

delta_t <- gammaj_t
gamma_t <- gammaj_t

# CACE-T Logit

gammaj_l_t <- matrix(0,(k1+1),(k1+1))

# CACE-T Variables to calculate Off-Diagonal WLS and Logit GEE - the B Matrix in Paper

gammaj_b_t <- matrix(0,(nx+2),(k1+1))
gamma_b_t  <- gammaj_b_t

# CACE-T Variables for all pieces

nv_t <- nx + k1 + 3

phij_all_t <- matrix(0,nv_t,1)
phi_all_t  <- phij_all_t

gammaj_all_t <- matrix(0,nv_t,nv_t)
gamma_all_t  <- gammaj_all_t

delta_all_t <- gammaj_all_t

# CACE-11 INITIALIZE VARIABLES FOR VARIANCES FOR CACE-11 IPW ESTIMATOR

#  CACE-11 Variables to calculate WLS GEE

phij_11   <- matrix(0,(nx+2),1)
gammaj_11 <- matrix(0,(nx+2),(nx+2))

delta_11 <- gammaj_11
gamma_11 <- gammaj_11

# CACE-11 Variables to calculate Logit Treatment and Control GEE

gammaj_l_11 <- matrix(0,(k1+k0+2),(k1+k0+2))

# CACE-11 Variables to calculate Off-Diagonal WLS and Logit GEE - the B Matrix in Paper

gammaj_b_11 <- matrix(0,(nx+2),(k1+k0+2))
gamma_b_11  <- gammaj_b_11

# CACE-11 Variables for all pieces

nv <- nx + k1 + k0 + 4

phij_all_11 <- matrix(0,nv,1)
phi_all_11  <- phij_all_11

gammaj_all_11 <- matrix(0,nv,nv)
gamma_all_11  <- gammaj_all_11

delta_all_11 <- gammaj_all_11

# CACE-IV INITIALIZE VARIABLES FOR VARIANCES FOR CACE-IV IPW ESTIMATOR

# CACE-IV Variables to calculate WLS GEE

phij_iv   <- matrix(0,(nx+2),1)
gammaj_iv <- matrix(0,(nx+2),(nx+2))

delta_iv <- gammaj_iv
gamma_iv <- gammaj_iv

# CACE-IV Variables to calculate Off-Diagonal WLS and Logit GEE - the B Matrix in Paper

gammaj_b_iv <- matrix(0,(k0+3),(nx+2))
gamma_b_iv  <- gammaj_b_iv

# CACE-IV Variables for all pieces

nv_iv <- nx + k0 + 5

phij_all_iv <- matrix(0,nv_iv,1)
phi_all_iv  <- phij_all_iv

gammaj_all_iv <- matrix(0,nv_iv,nv_iv)
gamma_all_iv  <- gammaj_all_iv

delta_all_iv <- gammaj_all_iv

###
# DO OVER ALL CLUSTERS
###

for (j in 1:mclus) {
  # CACE_TC CALCULATE GEE VARIANCE PIECES FOR CACE_TC IPW ESTIMATOR

  # CACE-TC WLS variance

  tj <- mydat$ttt[mydat$clusv == j]
  cj <- mydat$ccc[mydat$clusv == j]

  onesj <- cj
  onesj <- 1
  rj <- mydat$serv_var[mydat$clusv == j]

  if (nx > 0) {
    xj  <- xdat1[xdat1$clusv == j,2:(nx+1)]
    xj  <- as.matrix(xj)
    }

  x1j <- xtdat1[xtdat1$clusv == j,2:(k1+1)]
  x0j <- xcdat1[xcdat1$clusv == j,2:(k0+1)]

  x1j <- as.matrix(x1j)
  x0j <- as.matrix(x0j)

  uj  <- u_tc[u_tc$clusv == j,2]

  wj  <- wgt_tc[wgt_tc$clusv == j,2]
  wwj <- wgt_ww[wgt_ww$clusv == j,2]

  wuj <- wj*uj
  wtj <- wj*tj
  wcj <- wj*cj

  if (nx > 0) {
    wxj <- sqrt(wj)*xj
    }

  phij[1,1] <- t(tj) %*% wuj
  phij[2,1] <- t(cj) %*% wuj

  if (nx>0) {
    phij[3:(2+nx),1] <- t(xj) %*% wuj
    }

  deltaj <- phij %*% t(phij)
  delta  <- delta + deltaj

  gammaj[1,1]        <- t(tj) %*% wj
  gammaj[1,2]        <- 0
  gammaj[2,1]        <- 0
  gammaj[2,2]        <- t(cj) %*% wj

  if (nx>0) {
    gammaj[1,3:(2+nx)] <- t(wtj) %*% xj
    gammaj[2,3:(2+nx)] <- t(wcj) %*% xj
    gammaj[3:(2+nx),1] <- t(xj) %*% wtj
    gammaj[3:(2+nx),2] <- t(xj) %*% wcj
    gammaj[3:(2+nx),3:(2+nx)] <- t(wxj) %*% wxj
    }

  gamma <- gamma + gammaj

  # CACE-TC Treatment logit variance

  treat <- tj[1]

  if (treat == 1) {
    m1 <- m1 + 1
    } else {
      m0 <- m0 + 1
      }

  ej1    <- mydat$predt[mydat$clusv == j]
  ej0    <- mydat$predc[mydat$clusv == j]
  etaj1  <- mydat$resid_lt[which(mydat$clusv == j),]   # Not sure why need to use this format
  etaj0  <- mydat$resid_lc[which(mydat$clusv == j),]

  wwetaj1 <- wwj*etaj1
  wwetaj0 <- wwj*etaj0
  wwtj    <- wwj*tj
  wwcj    <- wwj*cj

  twwetaj1 <- tj*wwetaj1
  wweej1   <- wwj*ej1*(1-ej1)
  twweej1  <- tj*wwj*ej1*(1-ej1)
  twweexj1 <- sqrt(twweej1)*x1j

  phij_l1[1,1]        <- t(tj) %*% wwetaj1
  phij_l1[2:(k1+1),1] <- t(x1j) %*% twwetaj1

  deltaj_l1 <- phij_l1 %*% t(phij_l1)
  delta_l1  <- delta_l1 + deltaj_l1

  gammaj_l1[1,1]        <- t(tj) %*% wweej1
  gammaj_l1[1,2:(k1+1)] <- t(twweej1) %*% x1j
  gammaj_l1[2:(k1+1),1] <- t(x1j) %*% twweej1
  gammaj_l1[2:(k1+1),2:(k1+1)] <- t(twweexj1) %*% twweexj1

  gamma_l1 <- gamma_l1 + gammaj_l1

  # CACE-TC Control logit variance

  cwwetaj0 <- cj*wwetaj0
  wweej0   <- wwj*ej0*(1-ej0)
  cwweej0  <- cj*wwj*ej0*(1-ej0)
  cwweexj0 <- sqrt(cwweej0)*x0j

  phij_l0[1,1]        <- t(cj) %*% wwetaj0
  phij_l0[2:(k0+1),1] <- t(x0j) %*% cwwetaj0

  deltaj_l0 <- phij_l0 %*% t(phij_l0)
  delta_l0  <- delta_l0 + deltaj_l0

  gammaj_l0[1,1]        <- t(cj) %*% wweej0
  gammaj_l0[1,2:(k0+1)] <- t(cwweej0) %*% x0j
  gammaj_l0[2:(k0+1),1] <- t(x0j) %*% cwweej0
  gammaj_l0[2:(k0+1),2:(k0+1)] <- t(cwweexj0) %*% cwweexj0

  gamma_l0 <- gamma_l0 + gammaj_l0

  # CACE-TC B Matrix

  nrue0    <- (1-rj)*uj*ej0*(1-ej0)
  tnrue0   <- tj*(1-rj)*uj*ej0*(1-ej0)
  if (nx>0) {
    tnruex0  <- tnrue0*xj
    }

  nrue1    <- (1-rj)*uj*ej1*(1-ej1)
  cnrue1   <- cj*(1-rj)*uj*ej1*(1-ej1)
  if (nx>0) {
    cnruex1 <- cnrue1*xj
    }

  gammaj_b[1,(k1+2)]           <- -t(tj) %*% nrue0
  gammaj_b[1,(k1+3):(k1+k0+2)] <- -t(tnrue0) %*% x0j

  gammaj_b[2,1]        <- -t(cj) %*% nrue1
  gammaj_b[2,2:(k1+1)] <- -t(cnrue1) %*% x1j

  if (nx>0) {
    gammaj_b[3:(nx+2),1]        <- -t(xj) %*% cnrue1
    gammaj_b[3:(nx+2),2:(k1+1)] <- -t(cnruex1) %*% x1j

    gammaj_b[3:(nx+2),(k1+2)]           <- -t(xj) %*% tnrue0
    gammaj_b[3:(nx+2),(k1+3):(k1+k0+2)] <- -t(tnruex0) %*% x0j
    }

  # CACE-TC The whole shebang by putting A, B, and C matrixes together

  phij_all[1:(nx+2),1]         <- phij
  phij_all[(nx+3):(nx+k1+3),1] <- phij_l1
  phij_all[(nx+k1+4):nv,1]     <- phij_l0

  deltaj_all <- phij_all %*% t(phij_all)

  delta_all  <- delta_all + deltaj_all

  gammaj_all[1:(nx+2),1:(nx+2)] <- gammaj

  gammaj_all[1:(nx+2),(nx+3):nv]<- gammaj_b

  # CACE-TC Put logits together

  gammaj_l[1:(k1+1),1:(k1+1)]                 <- gammaj_l1
  gammaj_l[(k1+2):(k0+k1+2),(k1+2):(k0+k1+2)] <- gammaj_l0

  gammaj_all[(nx+3):nv,(nx+3):nv] <- gammaj_l

  gamma_all <- gamma_all + gammaj_all

  # CACE-T CALCULATE GEE VARIANCE PIECES FOR CACE_T IPW ESTIMATOR
  # CACE_T WLS variance

  uj_t  <- u_t[u_t$clusv == j,2]

  wj_t  <- wgt_t[wgt_t$clusv == j,2]

  wuj_t <- wj_t*uj_t
  wtj_t <- wj_t*tj
  wcj_t <- wj_t*cj

  if (nx > 0) {
    wxj_t <- sqrt(wj_t)*xj
    }

  phij_t[1,1] <- t(tj) %*% wuj_t
  phij_t[2,1] <- t(cj) %*% wuj_t

  if (nx>0) {
    phij_t[3:(2+nx),1] <- t(xj) %*% wuj_t
    }

  deltaj_t <- phij_t %*% t(phij_t)
  delta_t  <- delta_t + deltaj_t

  gammaj_t[1,1]        <- t(tj) %*% wj_t
  gammaj_t[1,2]        <- 0
  gammaj_t[2,1]        <- 0
  gammaj_t[2,2]        <- t(cj) %*% wj_t

  if (nx>0) {
    gammaj_t[1,3:(2+nx)] <- t(wtj_t) %*% xj
    gammaj_t[2,3:(2+nx)] <- t(wcj_t) %*% xj
    gammaj_t[3:(2+nx),1] <- t(xj) %*% wtj_t
    gammaj_t[3:(2+nx),2] <- t(xj) %*% wcj_t
    gammaj_t[3:(2+nx),3:(2+nx)] <- t(wxj_t) %*% wxj_t
    }

  gamma_t <- gamma_t + gammaj_t

  # CACE_T B Matrix

  nrue1_t    <- uj_t*ej1*(1-ej1)
  cnrue1_t   <- cj*uj_t*ej1*(1-ej1)
  if (nx>0) {
    cnruex1_t  <- cnrue1_t*xj
    }

  gammaj_b_t[2,1]        <- -t(cj) %*% nrue1_t
  gammaj_b_t[2,2:(k1+1)] <- -t(cnrue1_t) %*% x1j

  if (nx>0) {
    gammaj_b_t[3:(nx+2),1]        <- -t(xj) %*% cnrue1_t
    gammaj_b_t[3:(nx+2),2:(k1+1)] <- -t(cnruex1_t) %*% x1j
    }

  # CACE-T The whole shebang by putting A, B, and C matrixes together

  phij_all_t[1:(nx+2),1]         <- phij_t
  phij_all_t[(nx+3):(nx+k1+3),1] <- phij_l1

  deltaj_all_t <- phij_all_t %*% t(phij_all_t)

  delta_all_t  <- delta_all_t + deltaj_all_t

  gammaj_all_t[1:(nx+2),1:(nx+2)] <- gammaj_t

  gammaj_all_t[1:(nx+2),(nx+3):nv_t] <- gammaj_b_t

  gammaj_l_t[1:(k1+1),1:(k1+1)] <- gammaj_l1

  gammaj_all_t[(nx+3):nv_t,(nx+3):nv_t] <- gammaj_l_t

  gamma_all_t <- gamma_all_t + gammaj_all_t

  # CACE-11 CALCULATE GEE VARIANCE PIECES FOR CACE_11 IPW ESTIMATOR
  # CACE-11 WLS variance

  uj_11  <- u_11[u_11$clusv == j,2]

  wj_11  <- wgt_11[wgt_11$clusv == j,2]

  wuj_11 <- wj_11*uj_11
  wtj_11 <- wj_11*tj
  wcj_11 <- wj_11*cj

  if (nx > 0) {
    wxj_11 <- sqrt(wj_11)*xj
    }

  phij_11[1,1] <- t(tj) %*% wuj_11
  phij_11[2,1] <- t(cj) %*% wuj_11

  if (nx>0) {
    phij_11[3:(2+nx),1] <- t(xj) %*% wuj_11
    }

  deltaj_11 <- phij_11 %*% t(phij_11)
  delta_11  <- delta_11 + deltaj_11

  gammaj_11[1,1]        <- t(tj) %*% wj_11
  gammaj_11[1,2]        <- 0
  gammaj_11[2,1]        <- 0
  gammaj_11[2,2]        <- t(cj) %*% wj_11

  if (nx>0) {
    gammaj_11[1,3:(2+nx)] <- t(wtj_11) %*% xj
    gammaj_11[2,3:(2+nx)] <- t(wcj_11) %*% xj
    gammaj_11[3:(2+nx),1] <- t(xj) %*% wtj_11
    gammaj_11[3:(2+nx),2] <- t(xj) %*% wcj_11
    gammaj_11[3:(2+nx),3:(2+nx)] <- t(wxj_11) %*% wxj_11
    }

  gamma_11 <- gamma_11 + gammaj_11

  # CACE-11 B Matrix

  nrue0_11    <- rj*uj_11*ej0*(1-ej0)
  tnrue0_11   <- tj*rj*uj_11*ej0*(1-ej0)
  if (nx>0) {
    tnruex0_11  <- tnrue0_11*xj
    }

  nrue1_11    <- rj*uj_11*ej1*(1-ej1)
  cnrue1_11   <- cj*rj*uj_11*ej1*(1-ej1)
  if (nx>0) {
    cnruex1_11 <- cnrue1_11*xj
    }

  gammaj_b_11[1,(k1+2)]           <- -t(tj) %*% nrue0_11
  gammaj_b_11[1,(k1+3):(k1+k0+2)] <- -t(tnrue0_11) %*% x0j

  gammaj_b_11[2,1]        <- -t(cj) %*% nrue1_11
  gammaj_b_11[2,2:(k1+1)] <- -t(cnrue1_11) %*% x1j

  if (nx>0) {
    gammaj_b_11[3:(nx+2),1]        <- -t(xj) %*% cnrue1_11
    gammaj_b_11[3:(nx+2),2:(k1+1)] <- -t(cnruex1_11) %*% x1j

    gammaj_b_11[3:(nx+2),(k1+2)]           <- -t(xj) %*% tnrue0_11
    gammaj_b_11[3:(nx+2),(k1+3):(k1+k0+2)] <- -t(tnruex0_11) %*% x0j
    }

  # CACE-11 The whole shebang by putting A, B, and C matrixes together */

  phij_all_11[1:(nx+2),1]         <- phij_11
  phij_all_11[(nx+3):(nx+k1+3),1] <- phij_l1
  phij_all_11[(nx+k1+4):nv,1]     <- phij_l0

  deltaj_all_11 <- phij_all_11 %*% t(phij_all_11)

  delta_all_11  <- delta_all_11 + deltaj_all_11

  gammaj_all_11[1:(nx+2),1:(nx+2)] <- gammaj_11

  gammaj_all_11[1:(nx+2),(nx+3):nv] <- gammaj_b_11

  # CACE-11 Put logits together

  gammaj_l_11[1:(k1+1),1:(k1+1)]                 <- gammaj_l1
  gammaj_l_11[(k1+2):(k0+k1+2),(k1+2):(k0+k1+2)] <- gammaj_l0

  gammaj_all_11[(nx+3):nv,(nx+3):nv] <- gammaj_l_11

  gamma_all_11 <- gamma_all_11 + gammaj_all_11

  # CACE-IV WLS variance

  uj_iv  <- u_iv[u_iv$clusv == j,2]

  wwj_iv <- wgt_ww[wgt_ww$clusv == j,2]   # Sample weight
  wj_iv  <- wgt_tc[wgt_tc$clusv == j,2]   # TC weight for pi denominator

  wwuj_iv <- wwj_iv*uj_iv
  wwtj_iv <- wwj_iv*tj
  wwcj_iv <- wwj_iv*cj

  if (nx > 0) {
    wwxj_iv <- sqrt(wwj_iv)*xj
    }

  phij_iv[1,1] <- t(tj) %*% wwuj_iv
  phij_iv[2,1] <- t(cj) %*% wwuj_iv

  if (nx>0) {
    phij_iv[3:(2+nx),1] <- t(xj) %*% wwuj_iv
    }

  deltaj_iv <- phij_iv %*% t(phij_iv)
  delta_iv  <- delta_iv + deltaj_iv

  gammaj_iv[1,1]        <- t(tj) %*% wwj_iv
  gammaj_iv[1,2]        <- 0
  gammaj_iv[2,1]        <- 0
  gammaj_iv[2,2]        <- t(cj) %*% wwj_iv

  if (nx>0) {
    gammaj_iv[1,3:(2+nx)] <- t(wwtj_iv) %*% xj
    gammaj_iv[2,3:(2+nx)] <- t(wwcj_iv) %*% xj
    gammaj_iv[3:(2+nx),1] <- t(xj) %*% wwtj_iv
    gammaj_iv[3:(2+nx),2] <- t(xj) %*% wwcj_iv
    gammaj_iv[3:(2+nx),3:(2+nx)] <- t(wwxj_iv) %*% wwxj_iv
    }

  gamma_iv <- gamma_iv + gammaj_iv

  # CACE-IV B Matrix

  ntja <- tj*wwj_iv           # Sum of TC weights
  ntj <- t(ntja) %*% ntja
  gammaj_b_iv[(k0+3),1] <- -1
  gammaj_b_iv[(k0+3),2] <-  1

  if (nx>0) {
    gammaj_b_iv[(k0+3),3:(2+nx)] <- 0  # Used to be mx_diff
    }

  # CACE-IV The whole shebang by putting A, B, and C matrixes together

  phij_all_iv[1:(nx+2),1]         <- phij_iv
  phij_all_iv[(nx+3):(nx+k0+3),1] <- phij_l0
  phij_all_iv[(nx+k0+4),1]        <- (ntj*prop_tc) - (t(tj) %*% wj_iv)  # Uses T weights for pi
  phij_all_iv[(nx+k0+5),1]  <- (imp_itt-xb_diff) - (imp_iv*prop_tc)     # Uses T weights for pi

  deltaj_all_iv <- phij_all_iv %*% t(phij_all_iv)

  delta_all_iv  <- delta_all_iv + deltaj_all_iv

  gammaj_all_iv[1:(nx+2),1:(nx+2)]     <- gammaj_iv       # A matrix

  gammaj_all_iv[(nx+3):nv_iv,1:(nx+2)] <- gammaj_b_iv     # B matrix

  # CACE-IV C matrix

  cwweej0  <- cj*wwj*ej0*(1-ej0)

  eej0_iv   <- (1-rj)*wwj_iv*ej0*(1-ej0)
  teej0_iv  <- tj*wwj_iv*(1-rj)*ej0*(1-ej0)

  gammaj_all_iv[(nx+3):(nx+3+k0),(nx+3):(nx+3+k0)] <- gammaj_l0
  gammaj_all_iv[(nv_iv-1),(nx+3):(nx+3)]           <- t(tj) %*% eej0_iv
  gammaj_all_iv[(nv_iv-1),(nx+4):(nx+3+k0)]        <- t(teej0_iv) %*% x0j
  gammaj_all_iv[(nv_iv-1),(nv_iv-1)]               <- -ntj
  gammaj_all_iv[nv_iv,(nv_iv-1)]                   <- imp_iv
  gammaj_all_iv[nv_iv,nv_iv]                       <- prop_tc

  gamma_all_iv <- gamma_all_iv + gammaj_all_iv

  }  # For j

###
# CALCULATE FULL GEE SEs
###

# CACE-TC
# CACE-TC WLS

delta_hat <- delta/mclus
gamma_hat <- gamma/mclus

ig <- solve(gamma_hat)     # Solve is the inverse function in R

# CACE-TC WLS All variables

sand_vara <- ig %*% delta_hat %*% t(ig)
sand_vara <- sand_vara/mclus
sand_sea  <- sqrt(diag(sand_vara))

# CACE-TC WLS Impacts by differencing T and C mean estimates

lamda <- matrix(0,1,(nx+2))
lamda[1,1] <-  1
lamda[1,2] <- -1

sand_var <- lamda %*% sand_vara %*% t(lamda)
sand_se_tc  <- sqrt(sand_var)

# Treatment Logits - All Variables

delta_hat_l1 <- delta_l1/m1
gamma_hat_l1 <- gamma_l1/m1

ig_l1 <- solve(gamma_hat_l1)

sand_vara_l1 <- ig_l1 %*% delta_hat_l1 %*% t(ig_l1)
sand_vara_l1 <- sand_vara_l1/m1
sand_sea_l1  <- sqrt(diag(sand_vara_l1))

# Control Logits - All Variables

delta_hat_l0 <- delta_l0/m0
gamma_hat_l0 <- gamma_l0/m0

ig_l0 <- solve(gamma_hat_l0)

sand_vara_l0 <- ig_l0 %*% delta_hat_l0 %*% t(ig_l0)
sand_vara_l0 <- sand_vara_l0/m0
sand_sea_l0  <- sqrt(diag(sand_vara_l0))

# CACE-TC Whole shebang

delta_hat_all <- delta_all/mclus
gamma_hat_all <- gamma_all/mclus

ig_all = solve(gamma_hat_all)

# CACE-TC Shebang All variables

sand_vara_all <- ig_all %*% delta_hat_all %*% t(ig_all)
sand_vara_all <- sand_vara_all/mclus
sand_sea_all  <- sqrt(diag(sand_vara_all))

# CACE-TC Shebang Impacts by differencing T and C mean estimates

lamda_all = matrix(0,1,nv)
lamda_all[1,1] <-  1
lamda_all[1,2] <- -1

sand_var_all <- lamda_all %*% sand_vara_all %*% t(lamda_all)
sand_se_all_tc  <- sqrt(sand_var_all)

# CACE-T
# CACE-T WLS

delta_hat_t <- delta_t/mclus
gamma_hat_t <- gamma_t/mclus

ig_t <- solve(gamma_hat_t)

# CACE-T WLS All variables

sand_vara_t <- ig_t %*% delta_hat_t %*% t(ig_t)
sand_vara_t <- sand_vara_t/mclus
sand_sea_t  <- sqrt(diag(sand_vara_t))

# CACE-T WLS Impacts by differencing T and C mean estimates

lamda <- matrix(0,1,(nx+2))
lamda[1,1] <-  1
lamda[1,2] <- -1

sand_var_t <- lamda %*% sand_vara_t %*% t(lamda)
sand_se_t  <- sqrt(sand_var_t)

# CACE-T Whole shebang

delta_hat_all_t <- delta_all_t/mclus
gamma_hat_all_t <- gamma_all_t/mclus

ig_all_t <- solve(gamma_hat_all_t)

# CACE-T Shebang All variables

sand_vara_all_t <- ig_all_t %*% delta_hat_all_t %*% t(ig_all_t)
sand_vara_all_t <- sand_vara_all_t/mclus
sand_sea_all_t  <- sqrt(diag(sand_vara_all_t))

# CACE-T Shebang Impacts by differencing T and C mean estimates

lamda_all_t = matrix(0,1,nv_t)
lamda_all_t[1,1] <-  1
lamda_all_t[1,2] <- -1

sand_var_all_t <- lamda_all_t %*% sand_vara_all_t %*% t(lamda_all_t)
sand_se_all_t  <- sqrt(sand_var_all_t)

# CACE-11
# CACE-11 WLS

delta_hat_11 <- delta_11/mclus
gamma_hat_11 <- gamma_11/mclus

ig_11 <- solve(gamma_hat_11)

# CACE-11 WLS All variables

sand_vara_11 <- ig_11 %*% delta_hat_11 %*% t(ig_11)
sand_vara_11 <- sand_vara_11/mclus
sand_sea_11  <- sqrt(diag(sand_vara_11))

# CACE-11 WLS Impacts by differencing T and C mean estimates

lamda = matrix(0,1,(nx+2))
lamda[1,1] <-  1
lamda[1,2] <- -1

sand_var_11 <- lamda %*% sand_vara_11 %*% t(lamda)
sand_se_11  <- sqrt(sand_var_11)

# CACE-11 Whole shebang

delta_hat_all_11 <- delta_all_11/mclus
gamma_hat_all_11 <- gamma_all_11/mclus

ig_all_11 <- solve(gamma_hat_all_11)

# CACE-11 Shebang All variables

sand_vara_all_11 <- ig_all_11 %*% delta_hat_all_11 %*% t(ig_all_11)
sand_vara_all_11 <- sand_vara_all_11/mclus
sand_sea_all_11  <- sqrt(diag(sand_vara_all_11))

# CACE-11 Shebang Impacts by differencing T and C mean estimates */

lamda_all <- matrix(0,1,nv)
lamda_all[1,1] <-  1
lamda_all[1,2] <- -1

sand_var_all_11 <- lamda_all %*% sand_vara_all_11 %*% t(lamda_all)
sand_se_all_11  <- sqrt(sand_var_all_11)

# CACE-IV
# CACE-IV WLS

delta_hat_iv <- delta_iv/mclus
gamma_hat_iv <- gamma_iv/mclus

ig_iv <- solve(gamma_hat_iv)

# CACE-IV WLS All variables

sand_vara_iv <- ig_iv %*% delta_hat_iv %*% t(ig_iv)
sand_vara_iv <- sand_vara_iv/mclus
sand_sea_iv  <- sqrt(diag(sand_vara_iv))

# CACE-IV WLS Impacts by differencing T and C mean estimates

lamda <- matrix(0,1,(nx+2))
lamda[1,1] <-  1
lamda[1,2] <- -1

sand_var_iv <- lamda %*% sand_vara_iv %*% t(lamda)
sand_se_iv  <- sqrt(sand_var_iv)

# CACE-IV Whole shebang

delta_hat_all_iv <- delta_all_iv/mclus
gamma_hat_all_iv <- gamma_all_iv/mclus

ig_all_iv <- solve(gamma_hat_all_iv)

# CACE-IV Shebang All variables

sand_vara_all_iv <- ig_all_iv %*% delta_hat_all_iv %*% t(ig_all_iv)

sand_vara_all_iv <- sand_vara_all_iv/mclus
sand_sea_all_iv  <- sqrt(diag(sand_vara_all_iv))

# CACE-IV Shebang Impacts by differencing T and C mean estimates

lamda_all_iv <- matrix(0,1,nv_iv)
lamda_all_iv[1,1] <-  1
lamda_all_iv[1,2] <- -1

sand_var_all_iv <- lamda_all_iv %*% sand_vara_all_iv %*% t(lamda_all_iv)
sand_se_all_iv  <- sqrt(sand_var_all_iv)

sand_var_all_iv_tau <- sand_vara_all_iv[nv_iv,nv_iv]
sand_se_all_iv_tau  <- sqrt(sand_var_all_iv_tau)

# Conduct significance testing

# Df Updated since get degrees of freedom that are too small for clustered designs

if (df_wls == 0) {
  df_itt <- mclus - nx - 2
  df_tc  <- mclus - nx - k1 - k0 - 4
  df_t   <- mclus - nx - k1 - 3
  df_11  <- mclus - nx - k1 - k0 - 4
  df_iv  <- mclus - nx - k0 - 5
} else if (df_wls == 1) {
  df_itt <- mclus - nx - 2
  df_tc  <- mclus - nx - 2
  df_t   <- mclus - nx - 2
  df_11  <- mclus - nx - 2
  df_iv  <- mclus - nx - 2
}

tstat_itt  <- abs(imp_itt)/sand_se_iv

tstat_tc  <- abs(imp_tc)/sand_se_all_tc
tstat_tca <- abs(imp_tc)/sand_se_tc

tstat_t  <- abs(imp_t)/sand_se_all_t
tstat_ta <- abs(imp_t)/sand_se_t

tstat_11  <- abs(imp_11)/sand_se_all_11
tstat_11a <- abs(imp_11)/sand_se_11

tstat_iv  <- abs(imp_iv)/sand_se_all_iv_tau

sand_se_direct <- sand_se_iv / prop_tc
tstat_iva <- abs(imp_iv)/sand_se_direct

p_itt <- 2*pt(tstat_itt, df=df_itt, lower.tail=FALSE)
p_tc  <- 2*pt(tstat_tc,  df=df_tc,  lower.tail=FALSE)
p_t   <- 2*pt(tstat_t,   df=df_t,   lower.tail=FALSE)
p_11  <- 2*pt(tstat_11,  df=df_11,  lower.tail=FALSE)
p_iv  <- 2*pt(tstat_iv,  df=df_iv,  lower.tail=FALSE)

p_tca <- 2*pt(tstat_tca, df=df_tc,  lower.tail=FALSE)
p_ta  <- 2*pt(tstat_ta,  df=df_t,   lower.tail=FALSE)
p_11a <- 2*pt(tstat_11a, df=df_11,  lower.tail=FALSE)
p_iva <- 2*pt(tstat_iva, df=df_iv,  lower.tail=FALSE)

} # END IF (se_ipw == 1)


Sum_data  <- sum_dat1
Persons   <- nobs_dat1
Clusters  <- mobs_dat1

orig_samp  <- cbind(Persons,Clusters)
#orig_samp  <- format(orig_samp, justify = "left", width = 10)
orig_samp1 <- as.data.frame(orig_samp)

Persons  <- nobs_dat
Clusters <- mobs_dat

anal_samp  <- cbind(Persons,Clusters)
#anal_samp  <- format(anal_samp, justify = "left", width = 10)
anal_samp1 <- as.data.frame(anal_samp)


# Print Data Summary and Logit and WLS Regression Results

tits1  <- sprintf("SUMMARY OF INPUT CSV DATA FILE")
tits2  <- sprintf("ORIGINAL SAMPLE SIZES")
tits3  <- sprintf("ANALYSIS SAMPLE SIZES FOR OUTCOME %s AFTER REMOVING MISSING DATA",yname)

titl1  <- sprintf("LOGIT RESULTS FOR THE TREATMENT GROUP FOR OUTCOME %s",yname)
titl1a <- sprintf("WITH WEIGHTS %s",wgt)

titl2 <- sprintf("LOGIT RESULTS FOR THE CONTROL GROUP FOR OUTCOME %s",yname)
titw1 <- sprintf("WLS RESULTS FOR THE ITT ESTIMATOR FOR OUTCOME %s",yname)
titw2 <- sprintf("WLS RESULTS FOR THE CACE_T ESTIMATOR FOR OUTCOME %s",yname)
titw3 <- sprintf("WLS RESULTS FOR THE CACE_TC2 ESTIMATOR FOR OUTCOME %s",yname)
titw4 <- sprintf("WLS RESULTS FOR THE CACE_11 ESTIMATOR FOR OUTCOME %s",yname)

blank <- c("")

# PRINT DATA SUMMARY AND LOGIT AND REGRESSION RESULTS

if (iy == 1) {
  d_app <- c('FALSE')      # Do not append the output files
  } else {
    d_app <- c('TRUE')     # Append the output files with new output
    }

if (iy == 1) {
  #cat(blank,tits1,blank, sep="\n")
  #print(Sum_data)
  cat(blank,tits2,blank, sep="\n")
  print(orig_samp1, row.names = F)
}

cat(blank,tits3,blank, sep="\n")
print(anal_samp1, row.names = F)

sink(out_regr,d_app)       # Sink writes out results

if (iy == 1) {
  cat(blank,tits1,blank, sep="\n")
  print(Sum_data)
  cat(blank,tits2,blank, sep="\n")
  print(orig_samp1, row.names = F)
  }

cat(blank,tits3,blank, sep="\n")
print(anal_samp1, row.names = F)

if (wgt == 0) {
  cat(blank,titl1,blank, sep="\n")
  print(logt_res)

  cat(blank,titl2,blank, sep="\n")
  print(logc_res)

  cat(blank,titw1,blank, sep="\n")
  print(itt_res)

  cat(blank,titw2,blank, sep="\n")
  print(cace_tc_res)

  cat(blank,titw3,blank, sep="\n")
  print(cace_t_res)

  cat(blank,titw4,blank, sep="\n")
  print(cace_11_res)

} else {
  cat(blank,titl1,titl1a,blank, sep="\n")
  print(logt_res)

  cat(blank,titl2,titl1a,blank, sep="\n")
  print(logc_res)

  cat(blank,titw1,titl1a,blank, sep="\n")
  print(itt_res)

  cat(blank,titw2,titl1a,blank, sep="\n")
  print(cace_tc_res)

  cat(blank,titw3,titl1a,blank, sep="\n")
  print(cace_t_res)

  cat(blank,titw4,titl1a, sep="\n")
  print(cace_11_res)
  }

sink()     # Turns off sink

# Create results for printing

if (se_ipw == 1) {
  # Merge the impact, se, tstat, and pvalue columns
  res_itt <- do.call("cbind", list(imp_itt, sand_se_iv, tstat_itt, p_itt,
                                   sand_se_iv, tstat_itt, p_itt))

  res_t  <- do.call("cbind", list(imp_t, sand_se_all_t, tstat_t, p_t,
                                  sand_se_t, tstat_ta, p_ta))

  res_iv <- do.call("cbind", list(imp_iv, sand_se_all_iv_tau, tstat_iv, p_iv,
                                  sand_se_direct, tstat_iva, p_iva))

  res_tc <- do.call("cbind", list(imp_tc, sand_se_all_tc, tstat_tc, p_tc,
                                  sand_se_tc, tstat_tca, p_tca))

  res_11 <- do.call("cbind", list(imp_11, sand_se_all_11, tstat_11, p_11,
                                  sand_se_11, tstat_11a, p_11a))

  } else if (se_ipw == 0) {
    res_itt <- do.call("cbind", list(imp_itt, sand_se_iv, tstat_iv, p_iv))

    res_t  <- do.call("cbind", list(imp_t, sand_se_all_t, tstat_t, p_t))

    res_iv <- do.call("cbind", list(imp_iv, sand_se_all_iv_tau, tstat_iv, p_iv))

    res_tc <- do.call("cbind", list(imp_tc, sand_se_all_tc, tstat_tc, p_tc))

    res_11 <- do.call("cbind", list(imp_11, sand_se_all_11, tstat_11, p_11))
    }

# Stack the results across estimators
resg <- rbind(res_itt,res_t)
resg <- rbind(resg,res_iv)
resg <- rbind(resg,res_tc)
resg <- rbind(resg,res_11)

# Format
Impact  <- sprintf("%.4f",resg[,1])
Std_err <- sprintf("%.4f",resg[,2])
t_value <- sprintf("%.4f",resg[,3])
p_value <- sprintf("%.4f",resg[,4])

# Create Estimator name variable
Estimator <- c("ITT","CACE_T","CACE_TC1","CACE_TC2","CACE_11")
Estimator <- format(Estimator, justify = "left", width = 9)

# Merge again
resg1 <- cbind(Estimator,Impact)
resg1 <- cbind(resg1,Std_err)
resg1 <- cbind(resg1,t_value)
resg1 <- cbind(resg1,p_value)

# Make it a data frame
resg2 <- as.data.frame(resg1)

# Add significance stars
resg2$signif <- c("")
resg2$signif[which(resg2$p_value <= .01)] <- c('***')
resg2$signif[which((resg2$p_value <= .05) & (resg2$p_value > .01))] <- c('**')
resg2$signif[which((resg2$p_value <= .10) & (resg2$p_value > .05))] <- c('*')

resg2$signif <- format(resg2$signif, justify = "left", width = 5)

# Add blank row for printing
blank <- c("")
resg2 <- rbind(blank,resg2)

# PRINT CACE ESTIMATES

if (marg == 1) {
  titc1 <- sprintf("Table %s. CACE GEE ESTIMATION RESULTS FOR THE MARGINAL LOGIT MODEL",iy)
  } else if (marg == 0) {
    titc1 <- sprintf("Table %s. CACE GEE ESTIMATION RESULTS FOR THE RANDOM EFFECTS LOGIT MODEL",iy)
    }

if (clus_var == 0) {
  titc2 <- sprintf("FOR OUTCOME VARIABLE %s: NON-CLUSTERED DESIGN",yname)
  } else if (clus_var != 0) {
    titc2 <- sprintf("FOR OUTCOME VARIABLE %s: CLUSTERED DESIGN USING %s",
                     yname,clus_var)
    }

if (se_ipw == 1) {
  titc3 <- c("(Standard errors adjust for estimation error in the IPW weights)")
  } else if (se_ipw == 0) {
    titc3 <- c("(Standard errors ignore estimation error in the IPW weights)")
    }

titc4a <- c("Notes. See Schochet (Journal of Causal inference, 2022) for details")
titc4b <- c("on the estimands and inverse probability weighting (IPW) estimators.")
titc4c <- c("ITT = intention-to-treat, CACE = complier average treatment effect.")

if (wgt == 0) {
  cat(blank,titc1,titc2,titc3,blank, sep="\n")
  print(resg2, row.names = F)
  cat(blank,titc4a,titc4b,titc4c,blank, sep="\n")
} else {
  cat(blank,titc1,titc2,titl1a,titc3,blank, sep="\n")
  print(resg2, row.names = F)
  cat(blank,titc4a,titc4b,titc4c,blank, sep="\n")
}

# Print population shares

# Format
prop_tc_t <- sprintf("%.4f",prop_tc_t)
prop_tc_c <- sprintf("%.4f",prop_tc_c)
prop_11_t <- sprintf("%.4f",prop_11_t)
prop_11_c <- sprintf("%.4f",prop_11_c)
prop_t_t  <- sprintf("%.4f",prop_t_t)
prop_t_c  <- sprintf("%.4f",prop_t_c)

# Merge T and C results
pp_t  <- cbind(prop_t_t,prop_t_c)
pp_tc <- cbind(prop_tc_t,prop_tc_c)
pp_11 <- cbind(prop_11_t,prop_11_c)

# Define population label
Population <- c("CACE_T","CACE_TC","CACE_11")
Population <- format(Population, justify = "left", width = 10)

# Merge everything together
pp_all <- rbind(pp_t,pp_tc)
pp_all <- rbind(pp_all,pp_11)

pp_all <- cbind(Population,pp_all)

# Make it a data frame
pp_alld <- as.data.frame(pp_all)

pp_alld <- rbind(blank,pp_alld)

# Rename columns and format
colnames(pp_alld)[2] <- c("Treatments")
colnames(pp_alld)[3] <- c("Controls")

pp_alld$Treatments <- format(pp_alld$Treatments, justify = "centre", width = 11)
pp_alld$Controls   <- format(pp_alld$Controls, justify = "centre", width = 8.5)

# Title
titpp1 <- sprintf("Table %sa. ESTIMATED CACE POPULATION SHARES, BY RESEARCH GROUP",iy)

cat(blank,titpp1,blank, sep="\n")
print(pp_alld, row.names = F)

# Print out results to file using sink()

sink(out_est,d_app)

if (wgt == 0) {
  cat(blank,titc1,titc2,titc3,blank, sep="\n")
  print(resg2, row.names = F)
  cat(blank,titc4a,titc4b,titc4c,blank, sep="\n")
} else {
  cat(blank,titc1,titc2,titl1a,titc3,blank, sep="\n")
  print(resg2, row.names = F)
  cat(blank,titc4a,titc4b,titc4c,blank, sep="\n")
}

cat(blank,titpp1,blank, sep="\n")
print(pp_alld, row.names = F)

sink()

}  # End Y Loop

} # End of function

