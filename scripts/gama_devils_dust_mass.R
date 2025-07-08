library(magicaxis)
library(data.table)
library(foreach)
library(celestial)
library(doParallel)
library(ProSpect)
library(dplyr)
library(scales)
library(hyper.fit)
library(dftools)
library(Highlander)
library(matrixStats)
library(rhdf5)
library(stringr)

catalogueDir = "/Users/22252335/Documents/GAMA-DEVILS-SFR-AGN/data/"

gama_AGN = data.frame(fread("~/Documents/GAMA-DEVILS-SFR-AGN/data/GAMAMasterCat_BestFitParams.csv"))
gama_noAGN = fread(paste0(catalogueDir, 'gkvProSpectV02.csv'))
devilsd10_AGN = fread("~/Documents/GAMA-DEVILS-SFR-AGN/data/AGNTotalCat_MasterCat4.csv")
devilsd10_noAGN = readRDS(paste0(catalogueDir, 'DEVILS_D10ProSpectCat_02_02_2021_v0.3.rds')) #catalogue that Jess done without AGN contribution
devilsd10_noAGN = devilsd10_noAGN$cat
devilsd10_noAGN$area = 1.5
devilsd10_AGN$area = 1.5

LSS = fread("~/Documents/CSFRD-Compendium/GAMA-DEVILS/LSS.csv")
LSS_func = approxfun(LSS$lbt, LSS$super, yleft = 1.2137843, yright = 0.6403993)

idx_AGN = which(gama_AGN$CATAID %in% gama_noAGN$CATAID)
idx_noAGN = which(gama_noAGN$CATAID %in% gama_AGN$CATAID)

#get matching sets
gama_match_AGN = gama_AGN[idx_AGN, ]
gama_match_noAGN = gama_noAGN[idx_noAGN, ]
#sort by catalogue ID
sort_match_AGN = gama_match_AGN[order(gama_match_AGN$CATAID), ]
sort_match_noAGN = gama_match_noAGN[order(gama_match_noAGN$CATAID), ]
#make sure redshift column aligns with CATAID
sort_match_noAGN$z = sort_match_AGN$z
sort_match_noAGN$area = 217.54

zvec = seq(0, 30, 0.01)
lbtvec = cosdistTravelTime(z = zvec, ref = "737")
lbt2z = approxfun(
  lbtvec, zvec
)

lbt_bins = seq(0, 12, 0.75)
lbt_mids = lbt_bins[-length(lbt_bins)] + diff(lbt_bins) / 2

zbins = lbt2z(lbt_bins)
zmids = zbins[-length(zbins)] + diff(zbins) / 2

sm_bins = seq(3.6, 15.6, 0.6)
sm_mids = sm_bins[-length(sm_bins)] + diff(sm_bins) / 2

double_schechter = function(x, p){
  
  ## best starting  parm = c(11, -1, -1, -2, -3), 
  
  Mstar = p[1]
  alpha1 = p[2]
  alpha2 = p[3]
  phi1 = 10^p[4]
  phi2 = 10^p[5]
  
  mu = 10^x / 10^Mstar
  
  ret = log(10) * exp(-1*mu) * (phi1 * mu^(alpha1 + 1) + phi2 * mu^(alpha2 + 1))
  return( ret )
}

LL = function(p, Data){
  like = sum(dnorm(
    x = Data$yy,
    mean = Data$func(Data$xx, p),
    sd = Data$yyerr, 
    log = TRUE
  )) + Data$prior(p)
  
  if(is.infinite(like)){
    like = -99999
  }
  
  return( like )
}

mdustvec = seq(0.5, 15.5, 0.1)
compute_mass_function = function(zlo, zhi, z, x, x_err, areas, sm_bins, errFloor, fit_fun, vmax_bins = NULL, vmax = NULL, vmaxErr = NULL, meanxErr = NULL, do_optim = FALSE, do_fit = FALSE, do_fit_quantiles = TRUE, Niters = 1000, do_plot = FALSE, add = TRUE, pt.col = "black", ln.alpha = 1.0){
  
  sm_mids = sm_bins[-length(sm_bins)] + diff(sm_bins) / 2
  ddmm = abs(sm_bins[1] - sm_bins[2])
  vol = 4*pi/3 * (cosdistCoDist(z = zhi, ref = 'Planck18')^3 - cosdistCoDist(z = zlo, ref = 'Planck18')^3)
  
  zidx = z >= zlo & z < zhi & x > 0
  Mdust = log10(x[zidx])
  MdustErr = x_err[zidx] / (log(10) * x[zidx])
  vol = rep(vol, sum(zidx)) * areas[zidx] / (4*pi*(180/pi)^2)

  if(is.null(vmax) & is.null(vmaxErr) & is.null(vmax_bins) & is.null(meanxErr)){
    vvmax = foreach(j = 1:length(sm_mids), .combine = "c") %do% {
      midx = Mdust >= sm_bins[j] & Mdust < sm_bins[j+1]
      sum(vol[midx]^-1)
    }/ddmm
    vvmaxErr = foreach(j = 1:length(sm_mids), .combine = "c") %do% {
      midx = Mdust >= sm_bins[j] & Mdust < sm_bins[j+1]
      sqrt(sum(vol[midx]^-2))
    }/ddmm
    mmeanxErr = foreach(j = 1:length(sm_mids), .combine = "c") %do% {
      midx = Mdust >= sm_bins[j] & Mdust < sm_bins[j+1]
      if(sum(midx) == 0){
        0
      }else{
        mean(MdustErr[midx])
      }
    }
  }else{
    vvmax = vmax
    vvmaxErr = vmaxErr
    sm_mids = vmax_bins
    mmeanxErr = meanxErr
  }

  vvmaxErr = sqrt( vvmaxErr^2 + (errFloor * vvmax)^2 )

  # Mdust_lim = pmax( , 5.5 )
  Mdust_lim = sm_mids[which.max(vvmax)] + 0.3
  Mdust_uplim = rev(sm_mids)[max(which( sapply(rev(sm_mids), function(x) sum(vvmax[sm_mids <= x & sm_mids > Mdust_lim]==0) >= 2 ) ))]
  
  if(do_fit){

    if(do_optim){
      opt = optim(
        par = c(8, -1, -1, -2, -3), 
        fn = LL, 
        # method = "L-BFGS-B",
        control = list("fnscale" = -1), 
        Data = list(
          xx = sm_mids[vvmax > 0 & sm_mids > Mdust_lim], 
          yy = vvmax[vvmax > 0 & sm_mids > Mdust_lim], 
          yyerr = vvmaxErr[vvmax > 0 & sm_mids > Mdust_lim], 
          func = double_schechter, 
          prior = function(p){
            sum(
              # dnorm(p[1], x = 8.0, sd = 5, log = TRUE),
              # dnorm(p[2], x = -1.0, sd = 8, log = TRUE),
              # dnorm(p[3], x = -1.0, sd = 8, log = TRUE),
              # dnorm(p[4], x = -2.0, sd = 2, log = TRUE),
              # dnorm(p[5], x = -3.0, sd = 2, log = TRUE),
              0
            )
          }
        )
      )
      
      q16_fit = double_schechter(mdustvec, p = opt$par)
      q50_fit = double_schechter(mdustvec, p = opt$par)
      q84_fit = double_schechter(mdustvec, p = opt$par)
      highout = opt
      
    }else{
      highout = Highlander(
        parm = c(8, -1, -1, -2, -3), 
        Data = list(
          xx = sm_mids[vvmax > 0 & sm_mids > Mdust_lim], 
          yy = vvmax[vvmax > 0 & sm_mids > Mdust_lim], 
          yyerr = vvmaxErr[vvmax > 0 & sm_mids > Mdust_lim], 
          func = double_schechter, 
          prior = function(p){
            sum(
              # dnorm(p[1], x = 8.0, sd = 5, log = TRUE),
              # dnorm(p[2], x = -1.0, sd = 8, log = TRUE),
              # dnorm(p[3], x = -1.0, sd = 8, log = TRUE),
              # dnorm(p[4], x = -2.0, sd = 2, log = TRUE),
              # dnorm(p[5], x = -3.0, sd = 2, log = TRUE),
              0
            )
          }
        ), 
        likefunc = LL, 
        liketype = "max", 
        Niters = c(Niters, Niters), 
        NfinalMCMC = Niters, 
        parm.names = c("Mstar", "alpha1", "alpha2", "phi1", "phi2")
      )
      
      fit_samples = foreach(k = 1:dim(highout$LD_last$Posterior1)[1], .combine = "rbind") %do% {
        double_schechter(mdustvec, p = highout$LD_last$Posterior1[k, ])
      }
      q16_fit = colQuantiles(
        as.matrix(fit_samples), probs = 0.16
      )
      q50_fit = colQuantiles(
        as.matrix(fit_samples), probs = 0.50
      )
      q84_fit = colQuantiles(
        as.matrix(fit_samples), probs = 0.84
      )
    }

    df_fit = data.frame(
      "x" = 10^mdustvec,
      "Q50" = q50_fit,
      "Q16" = q16_fit,
      "Q84" = q84_fit
    )
    df_fit$ERR = 0.5*(df_fit$Q84 - df_fit$Q16)
    df_cosmic = data.frame(
      "Q50" = trapz(
        x = mdustvec,
        y = 10^mdustvec * q50_fit
      ),
      "Q16" = trapz(
        x = mdustvec,
        y = 10^mdustvec * q16_fit
      ),
      "Q84" = trapz(
        x = mdustvec,
        y = 10^mdustvec * q84_fit
      ),
      "Q50_lim" = trapz(
        x = mdustvec[mdustvec >= Mdust_lim & mdustvec <= Mdust_uplim],
        y = 10^mdustvec[mdustvec >= Mdust_lim & mdustvec <= Mdust_uplim] * q50_fit[mdustvec >= Mdust_lim & mdustvec <= Mdust_uplim]
      )
    )
    df_cosmic$ERR = 0.5*(df_cosmic$Q84 - df_cosmic$Q16)
  }else{
    df_fit = list(NULL)
    df_cosmic = list(NULL)
    highout = list(NULL)
  }
  
  df_vmax = data.frame(
    "x" = sm_mids,
    "vmax" = vvmax,
    "vmaxErr" = vvmaxErr,
    "meanxErr" = mmeanxErr,
    "mlim" = as.numeric(sm_mids >= Mdust_lim)
  )

  if(do_plot){
    if(add){
      points(
        10^sm_mids, vvmax,
        pch = 1, 
        cex = 1.5,
        col = alpha(pt.col, ln.alpha)
      )
    }else{
      magplot(
        10^sm_mids, vvmax, log = "xy", pch = 1, cex = 1.5, xlim = 10^c(2.5, 13.5), ylim = c(1e-10, 1), col = alpha(pt.col, ln.alpha),
        xlab = "Dust Mass [Msun]",
        ylab = "Phi [Mpc^-3 dex^-1]"
      )
    }

    magerr(
      10^sm_mids, vvmax, ylo = vvmaxErr, xlo = sqrt( (mmeanxErr * log(10))^2 * (10^sm_mids)^2 ), col = alpha(pt.col, ln.alpha)
    )
    if(do_fit){

      lines(
        10^mdustvec, q50_fit, lw = 2, col = alpha(pt.col, ln.alpha)
      )
      
      if(do_fit_quantiles){
        lines(
          10^mdustvec, double_schechter(mdustvec, p = highout$parm), col = alpha(pt.col, ln.alpha)
        )
        lines(
          10^mdustvec, q16_fit, lw = 2, lty = 2, col = alpha(pt.col, ln.alpha)
        )
        lines(
          10^mdustvec, q84_fit, lw = 2, lty = 2, col = alpha(pt.col, ln.alpha)
        )
      }
    }
    # abline(v = 10^Mdust_lim, col = pt.col)
    # abline(v = 10^Mdust_uplim, col = pt.col)
    
    legend(
      x = "topright", 
      legend = paste0(round(zlo,3), "< z <", round(zhi,3))
    )
  }
  
  return(
    list(
      "vmax" = df_vmax,
      "fit" = df_fit,
      "highout" = highout,
      "cosmic" = df_cosmic
    )
  )
}
compute_mass_functionV2 = function(zlo, zhi, z, x, x_err, areas, sm_bins, errFloor, df=5, vmax_bins = NULL, vmax = NULL, vmaxErr = NULL, meanxErr = NULL, do_optim = FALSE, do_fit = FALSE, do_fit_quantiles = TRUE, Niters = 1000, do_plot = FALSE, add = TRUE, pt.col = "black", ln.alpha = 1.0){
  
  sm_mids = sm_bins[-length(sm_bins)] + diff(sm_bins) / 2
  ddmm = abs(sm_bins[1] - sm_bins[2])
  vol = 4*pi/3 * (cosdistCoDist(z = zhi, ref = 'Planck18')^3 - cosdistCoDist(z = zlo, ref = 'Planck18')^3)
  
  zidx = z >= zlo & z < zhi & x > 0
  Mdust = log10(x[zidx])
  MdustErr = x_err[zidx] / (log(10) * x[zidx])
  vol = rep(vol, sum(zidx)) * areas[zidx] / (4*pi*(180/pi)^2)
  
  if(is.null(vmax) & is.null(vmaxErr) & is.null(vmax_bins) & is.null(meanxErr)){
    vvmax = foreach(j = 1:length(sm_mids), .combine = "c") %do% {
      midx = Mdust >= sm_bins[j] & Mdust < sm_bins[j+1]
      sum(vol[midx]^-1)
    }/ddmm
    vvmaxErr = foreach(j = 1:length(sm_mids), .combine = "c") %do% {
      midx = Mdust >= sm_bins[j] & Mdust < sm_bins[j+1]
      sqrt(sum(vol[midx]^-2))
    }/ddmm
    mmeanxErr = foreach(j = 1:length(sm_mids), .combine = "c") %do% {
      midx = Mdust >= sm_bins[j] & Mdust < sm_bins[j+1]
      if(sum(midx) == 0){
        0
      }else{
        median(MdustErr[midx])
      }
    }
  }else{
    vvmax = vmax
    vvmaxErr = vmaxErr
    sm_mids = vmax_bins
    mmeanxErr = meanxErr
  }
  
  vvmaxErr = sqrt(vvmaxErr^2 + (errFloor * vvmax)^2)
  
  # Mdust_lim = pmax( , 5.5 )
  Mdust_lim = sm_mids[which.max(vvmax)] + 0.3
  Mdust_uplim = rev(sm_mids)[max(which( sapply(rev(sm_mids), function(x) sum(vvmax[sm_mids <= x & sm_mids > Mdust_lim]==0) >= 2 ) ))]
  
  sm_mids_fit = sm_mids[vvmax > 0 & sm_mids > Mdust_lim]
  sm_mids_fit_err = mmeanxErr[vvmax > 0 & sm_mids > Mdust_lim]
  vvmax_fit = log10( vvmax[vvmax > 0 & sm_mids > Mdust_lim] )
  vvmax_fit_err = vvmaxErr[vvmax > 0 & sm_mids > Mdust_lim]/( vvmax[vvmax > 0 & sm_mids > Mdust_lim] * log(10))
  
  # print(sm_mids_fit)
  # print(sm_mids_fit_err)
  # print(vvmax_fit)
  # print(vvmax_fit_err)
  
  if(do_fit){
    
    if(do_optim){
      
      sm_spl = smooth.spline(
        c( sm_mids_fit, Mdust_uplim ), 
        c( vvmax_fit, -15 ),
        w = c( vvmax_fit_err^-2, max(vvmax_fit_err^-2) ),
        df = df
      )
      out = 10^predict(sm_spl, mdustvec)$y
      
      q16_fit = out
      q50_fit = out
      q84_fit = out
    }else{
      
      fit_samples = foreach(ii = 1:Niters, .combine = rbind, .errorhandling = "remove") %do% {
       
        xx_sample_jostle = rnorm(n = length(sm_mids_fit), mean = sm_mids_fit, sd = sm_mids_fit_err)
        yy_sample_jostle = rnorm(n = length(vvmax_fit), mean = vvmax_fit, sd = vvmax_fit_err)
        
        Mdust_uplim_jostle = rev(xx_sample_jostle)[max(which( sapply(rev(xx_sample_jostle), function(x) sum(vvmax[sm_mids <= x & sm_mids > Mdust_lim]==0) >= 2 ) ))]
        
        sm_spl = smooth.spline(
          c( sm_mids_fit, Mdust_uplim ), 
          c( yy_sample_jostle, -15 ),
          w = c( vvmax_fit_err^-2, 100 ),
          df = df
        )
        out = 10^predict(sm_spl, mdustvec)$y
      }
      
      q16_fit = colQuantiles(
        as.matrix(fit_samples), probs = 0.16
      )
      q50_fit = colQuantiles(
        as.matrix(fit_samples), probs = 0.50
      )
      q84_fit = colQuantiles(
        as.matrix(fit_samples), probs = 0.84
      )
    }
    
    df_fit = data.frame(
      "x" = 10^mdustvec,
      "Q50" = q50_fit,
      "Q16" = q16_fit,
      "Q84" = q84_fit
    )
    df_fit$ERR = 0.5*(df_fit$Q84 - df_fit$Q16)
    df_cosmic = data.frame(
      "Q50" = trapz(
        x = mdustvec,
        y = 10^mdustvec * q50_fit
      ),
      "Q16" = trapz(
        x = mdustvec,
        y = 10^mdustvec * q16_fit
      ),
      "Q84" = trapz(
        x = mdustvec,
        y = 10^mdustvec * q84_fit
      ),
      "Q50_lim" = trapz(
        x = mdustvec[mdustvec >= Mdust_lim & mdustvec <= Mdust_uplim],
        y = 10^mdustvec[mdustvec >= Mdust_lim & mdustvec <= Mdust_uplim] * q50_fit[mdustvec >= Mdust_lim & mdustvec <= Mdust_uplim]
      )
    )
    df_cosmic$ERR = 0.5*(df_cosmic$Q84 - df_cosmic$Q16)
  }else{
    df_fit = list(NULL)
    df_cosmic = list(NULL)
    highout = list(NULL)
  }
  
  df_vmax = data.frame(
    "x" = sm_mids,
    "vmax" = vvmax,
    "vmaxErr" = vvmaxErr,
    "meanxErr" = mmeanxErr,
    "mlim" = as.numeric(sm_mids >= Mdust_lim)
  )
  
  if(do_plot){
    if(add){
      points(
        10^sm_mids, vvmax,
        pch = 1, 
        cex = 1.5,
        col = alpha(pt.col, ln.alpha)
      )
    }else{
      magplot(
        10^sm_mids, vvmax, log = "xy", pch = 1, cex = 1.5, xlim = 10^c(2.5, 13.5), ylim = c(1e-10, 1), col = alpha(pt.col, ln.alpha),
        xlab = "Dust Mass [Msun]",
        ylab = "Phi [Mpc^-3 dex^-1]"
      )
    }
    
    magerr(
      10^sm_mids, vvmax, ylo = vvmaxErr, xlo = sqrt( (mmeanxErr * log(10))^2 * (10^sm_mids)^2 ), col = alpha(pt.col, ln.alpha)
    )
    if(do_fit){
      
      lines(
        10^mdustvec, q50_fit, lw = 2, col = alpha(pt.col, ln.alpha)
      )
      
      if(do_fit_quantiles){
        lines(
          10^mdustvec, q16_fit, lw = 2, lty = 2, col = alpha(pt.col, ln.alpha)
        )
        lines(
          10^mdustvec, q84_fit, lw = 2, lty = 2, col = alpha(pt.col, ln.alpha)
        )
      }
    }
    # abline(v = 10^Mdust_lim, col = pt.col)
    # abline(v = 10^Mdust_uplim, col = pt.col)
    
    legend(
      x = "topright", 
      legend = paste0(round(zlo,3), "< z <", round(zhi,3))
    )
  }
  
  return(
    list(
      "vmax" = df_vmax,
      "fit" = df_fit,
      "cosmic" = df_cosmic
    )
  )
}

stellar_mass_density = foreach(i = 1:(length(zbins)-1)) %do% {
  
  GAMA = compute_mass_function(
    zlo = zbins[i],
    zhi = zbins[i+1],
    z = sort_match_noAGN$z,
    x = sort_match_noAGN$StellarMass_50,
    x_err = 0.5 * (sort_match_noAGN$StellarMass_84 - sort_match_noAGN$StellarMass_16),
    areas = rep(217.54, length(sort_match_noAGN$z)),
    sm_bins = sm_bins,
    errFloor = 0.1,
    do_fit = FALSE,
    do_plot = TRUE,
    add = FALSE,
    pt.col = "purple"
  )
  DEVILS = compute_mass_function(
    zlo = zbins[i],
    zhi = zbins[i+1],
    z = devilsd10_noAGN$z,
    x = devilsd10_noAGN$StellarMass,
    x_err = 0.5 * (devilsd10_noAGN$StellarMass_UB - devilsd10_noAGN$StellarMass_LB),
    areas = rep(1.5, length(devilsd10_noAGN$z)),
    sm_bins = sm_bins,
    errFloor = 0.1,
    do_fit = FALSE,
    do_plot = TRUE,
    add = TRUE,
    pt.col = "cornflowerblue"
  )
  
  vmax_combine = foreach(j = 1:(length(sm_bins)-1) , .combine = rbind) %do% {
    
    if( GAMA$vmax$vmaxErr[j] <= ifelse(DEVILS$vmax$vmaxErr[j]==0, 999, DEVILS$vmax$vmaxErr[j]) ){
      if( !GAMA$vmax$mlim[j] ){
        c(DEVILS$vmax$x[j], DEVILS$vmax$vmax[j], DEVILS$vmax$vmaxErr[j], DEVILS$vmax$meanxErr[j])
      }else{
        if( (DEVILS$vmax$vmax[j] - GAMA$vmax$vmax[j])/(GAMA$vmax$vmaxErr[j] + 1e-323) > 5 ){
          c(DEVILS$vmax$x[j], DEVILS$vmax$vmax[j], DEVILS$vmax$vmaxErr[j], DEVILS$vmax$meanxErr[j])
        }else{
          c(GAMA$vmax$x[j], GAMA$vmax$vmax[j], GAMA$vmax$vmaxErr[j], DEVILS$vmax$meanxErr[j])
        }
      }
    }else{
      c(DEVILS$vmax$x[j], DEVILS$vmax$vmax[j], DEVILS$vmax$vmaxErr[j], DEVILS$vmax$meanxErr[j])
    }
    
  }
  
  fit = compute_mass_function(
    zlo = zbins[i],
    zhi = zbins[i+1],
    
    ## Not gonna use these 
    z = devilsd10_noAGN$z,
    x = devilsd10_noAGN$dustmass.total,
    x_err = 0.5 * (devilsd10_noAGN$StellarMass_UB - devilsd10_noAGN$StellarMass_LB),
    areas = devilsd10_noAGN$area,
    
    sm_bins = sm_bins,
    vmax_bins = vmax_combine[,1],
    vmax = vmax_combine[,2],
    vmaxErr = vmax_combine[,3],
    meanxErr = vmax_combine[,4],
    errFloor = 0,
    do_fit = TRUE,
    do_plot = TRUE,
    add = TRUE
  )
  
  fit$vmax$GAMA = GAMA$vmax$vmax
  fit$vmax$DEVILS = DEVILS$vmax$vmax
  
  points(
    10^fit$vmax$x[fit$vmax$mlim == 1],
    fit$vmax$vmax[fit$vmax$mlim == 1], 
    pch = 16, col = "red"
  )
  return(fit)
}
dust_mass_density = foreach(i = 1:(length(zbins)-1)) %do% {
  
  GAMA = compute_mass_function(
    zlo = zbins[i],
    zhi = zbins[i+1],
    z = sort_match_noAGN$z,
    x = sort_match_noAGN$DustMass_50,
    x_err = 0.5 * (sort_match_noAGN$DustMass_84 - sort_match_noAGN$DustMass_16),
    areas = rep(217.54, length(sort_match_noAGN$z)),
    sm_bins = sm_bins,
    errFloor = 0.1,
    do_fit = FALSE,
    do_plot = TRUE,
    add = FALSE,
    pt.col = "purple"
  )
  DEVILS = compute_mass_function(
    zlo = zbins[i],
    zhi = zbins[i+1],
    z = devilsd10_noAGN$z,
    x = devilsd10_noAGN$dustmass.total,
    x_err = 0.5 * (devilsd10_noAGN$dustmass.total_UB - devilsd10_noAGN$dustmass.total_LB),
    areas = rep(1.5, length(devilsd10_noAGN$z)),
    sm_bins = sm_bins,
    errFloor = 0.1,
    do_fit = FALSE,
    do_plot = TRUE,
    add = TRUE,
    pt.col = "cornflowerblue"
  )
  
  vmax_combine = foreach(j = 1:(length(sm_bins)-1) , .combine = rbind) %do% {
    if( GAMA$vmax$vmaxErr[j] <= ifelse(DEVILS$vmax$vmaxErr[j]==0, 999, DEVILS$vmax$vmaxErr[j]) ){
      if( !GAMA$vmax$mlim[j] ){
        c(DEVILS$vmax$x[j], DEVILS$vmax$vmax[j], DEVILS$vmax$vmaxErr[j], DEVILS$vmax$meanxErr[j])
      }else{
        if( (DEVILS$vmax$vmax[j] - GAMA$vmax$vmax[j])/(GAMA$vmax$vmaxErr[j] + 1e-323) > 5 ){
          c(DEVILS$vmax$x[j], DEVILS$vmax$vmax[j], DEVILS$vmax$vmaxErr[j], DEVILS$vmax$meanxErr[j])
        }else{
          c(GAMA$vmax$x[j], GAMA$vmax$vmax[j], GAMA$vmax$vmaxErr[j], GAMA$vmax$meanxErr[j])
        }
      }
    }else{
      c(DEVILS$vmax$x[j], DEVILS$vmax$vmax[j], DEVILS$vmax$vmaxErr[j], DEVILS$vmax$meanxErr[j])
    }
  }
  
  fit = compute_mass_function(
    zlo = zbins[i],
    zhi = zbins[i+1],
    
    ## Not gonna use these 
    z = devilsd10_noAGN$z,
    x = devilsd10_noAGN$dustmass.total,
    x_err = 0.5 * (devilsd10_noAGN$dustlum.total_UB - devilsd10_noAGN$dustlum.total_LB),
    areas = devilsd10_noAGN$area,
    
    sm_bins = sm_bins,
    vmax_bins = vmax_combine[,1],
    vmax = vmax_combine[,2],
    vmaxErr = vmax_combine[,3],
    meanxErr = vmax_combine[,4],
    errFloor = 0.0,
    do_fit = TRUE,
    do_plot = TRUE,
    add = TRUE,
    pt.col = "black"
  )
  fit$vmax$GAMA = GAMA$vmax$vmax
  fit$vmax$DEVILS = DEVILS$vmax$vmax
  points(
    10^fit$vmax$x[fit$vmax$mlim == 1],
    fit$vmax$vmax[fit$vmax$mlim == 1], 
    pch = 16, col = "red"
  )
  return(fit)
}
dust_mass_density_wAGN = foreach(i = 1:(length(zbins)-1)) %do% {
  
  GAMA = compute_mass_function(
    zlo = zbins[i],
    zhi = zbins[i+1],
    z = sort_match_AGN$z,
    x = sort_match_AGN$dustmass.total,
    x_err = 0.5 * (sort_match_AGN$dustmass.total_UB - sort_match_AGN$dustmass.total_LB),
    areas = rep(217.54, length(sort_match_AGN$z)),
    sm_bins = sm_bins,
    errFloor = 0.1,
    do_fit = FALSE,
    do_plot = TRUE,
    add = FALSE,
    pt.col = "purple"
  )
  DEVILS = compute_mass_function(
    zlo = zbins[i],
    zhi = zbins[i+1],
    z = devilsd10_AGN$z,
    x = devilsd10_AGN$dustmass.total,
    x_err = 0.5 * (devilsd10_AGN$dustmass.total_UB - devilsd10_AGN$dustmass.birth_LB),
    areas = rep(1.5, length(devilsd10_AGN$z)),
    sm_bins = sm_bins,
    errFloor = 0.1,
    do_fit = FALSE,
    do_plot = TRUE,
    add = TRUE,
    pt.col = "cornflowerblue"
  )
  
  vmax_combine = foreach(j = 1:(length(sm_bins)-1) , .combine = rbind) %do% {
    
    if( GAMA$vmax$vmaxErr[j] <= ifelse(DEVILS$vmax$vmaxErr[j]==0, 999, DEVILS$vmax$vmaxErr[j]) ){
      if( !GAMA$vmax$mlim[j] ){
        c(DEVILS$vmax$x[j], DEVILS$vmax$vmax[j], DEVILS$vmax$vmaxErr[j], DEVILS$vmax$meanxErr[j])
      }else{
        if( (DEVILS$vmax$vmax[j] - GAMA$vmax$vmax[j])/(GAMA$vmax$vmaxErr[j] + 1e-323) > 5 ){
          c(DEVILS$vmax$x[j], DEVILS$vmax$vmax[j], DEVILS$vmax$vmaxErr[j], DEVILS$vmax$meanxErr[j])
        }else{
          c(GAMA$vmax$x[j], GAMA$vmax$vmax[j], GAMA$vmax$vmaxErr[j], GAMA$vmax$meanxErr[j])
        }
      }
    }else{
      c(DEVILS$vmax$x[j], DEVILS$vmax$vmax[j], DEVILS$vmax$vmaxErr[j], DEVILS$vmax$meanxErr[j])
    }
    
  }
  
  fit = compute_mass_function(
    zlo = zbins[i],
    zhi = zbins[i+1],
    
    ## Not gonna use these 
    z = devilsd10_noAGN$z,
    x = devilsd10_noAGN$dustmass.total,
    x_err = 0.5 * (devilsd10_AGN$dustmass.total_UB - devilsd10_AGN$dustmass.birth_LB),
    areas = devilsd10_noAGN$area,
    
    sm_bins = sm_bins,
    vmax_bins = vmax_combine[,1],
    vmax = vmax_combine[,2],
    vmaxErr = vmax_combine[,3],
    meanxErr = vmax_combine[,4],
    errFloor = 0,
    do_fit = TRUE,
    do_plot = TRUE,
    add = TRUE
  )
  fit$vmax$GAMA = GAMA$vmax$vmax
  fit$vmax$DEVILS = DEVILS$vmax$vmax
  points(
    10^fit$vmax$x[fit$vmax$mlim == 1],
    fit$vmax$vmax[fit$vmax$mlim == 1], 
    pch = 16, col = "red"
  )
  return(fit)
}

# gc()
# csmh_edd_corr = foreach(i = 1:(length(zbins)-1), .combine = c, .errorhandling = "pass") %do% {
# 
#   sink("~/Documents/DustMassDensity/scripts/temp.txt")
#   message(lbt_mids[i])
#   png(paste0("~/Documents/DustMassDensity/plots/edd_bias_smf//lbt_", lbt_mids[i], ".png"), width = 7, height = 5, units = "in", res = 240)
#   GAMA = compute_mass_function(
#     zlo = zbins[i],
#     zhi = zbins[i+1],
#     z = sort_match_noAGN$z,
#     x = sort_match_noAGN$StellarMass_50,
#     x_err = 0.5 * (sort_match_noAGN$StellarMass_84 - sort_match_noAGN$StellarMass_16),
#     areas = rep(217.54, length(sort_match_noAGN$z)),
#     sm_bins = sm_bins,
#     errFloor = 0.0,
#     do_fit = FALSE,
#     do_plot = TRUE,
#     add = FALSE,
#     pt.col = "purple"
#   )
#   DEVILS = compute_mass_function(
#     zlo = zbins[i],
#     zhi = zbins[i+1],
#     z = devilsd10_noAGN$z,
#     x = devilsd10_noAGN$StellarMass,
#     x_err = 0.5 * (devilsd10_noAGN$StellarMass_UB - devilsd10_noAGN$StellarMass_LB),
#     areas = rep(1.5, length(devilsd10_noAGN$z)),
#     sm_bins = sm_bins,
#     errFloor = 0.0,
#     do_fit = FALSE,
#     do_plot = TRUE,
#     add = TRUE,
#     pt.col = "cornflowerblue"
#   )
# 
#   vmax_combine = foreach(j = 1:(length(sm_bins)-1) , .combine = rbind) %do% {
# 
#     if( GAMA$vmax$vmaxErr[j] <= ifelse(DEVILS$vmax$vmaxErr[j]==0, 999, DEVILS$vmax$vmaxErr[j]) ){
#       if( !GAMA$vmax$mlim[j] ){
#         c(DEVILS$vmax$x[j], DEVILS$vmax$vmax[j], DEVILS$vmax$vmaxErr[j], DEVILS$vmax$meanxErr[j])
#       }else{
#         if( (DEVILS$vmax$vmax[j] - GAMA$vmax$vmax[j])/(GAMA$vmax$vmaxErr[j] + 1e-323) > 5 ){
#           c(DEVILS$vmax$x[j], DEVILS$vmax$vmax[j], DEVILS$vmax$vmaxErr[j], DEVILS$vmax$meanxErr[j])
#         }else{
#           c(GAMA$vmax$x[j], GAMA$vmax$vmax[j], GAMA$vmax$vmaxErr[j], GAMA$vmax$meanxErr[j])
#         }
#       }
#     }else{
#       c(DEVILS$vmax$x[j], DEVILS$vmax$vmax[j], DEVILS$vmax$vmaxErr[j], DEVILS$vmax$meanxErr[j])
#     }
# 
#   }
# 
#   samples = foreach(j = 1:200) %do% {
#     vmax_bins_sample =  rnorm(length(vmax_combine[,1]), mean = vmax_combine[,1], sd = vmax_combine[,4])
#     vmax_bins_sample[is.nan(vmax_bins_sample)] = vmax_combine[,1][is.nan(vmax_bins_sample)]
# 
#     vmax_sample = rnorm(length(vmax_combine[,2]), mean = vmax_combine[,2], sd = vmax_combine[,3])
#     vmax_sample[is.nan(vmax_sample)] = vmax_combine[,2][is.nan(vmax_sample)]
#     cbind(vmax_bins_sample, vmax_sample)
#   }
# 
#   temp = foreach(j = 1:length(samples), .combine = c, .errorhandling = "remove") %do% {
# 
#     fit = suppressMessages(compute_mass_function(
#       zlo = zbins[i],
#       zhi = zbins[i+1],
# 
#       ## Not gonna use these
#       z = devilsd10_noAGN$z,
#       x = devilsd10_noAGN$StellarMass,
#       x_err = 0.5 * (devilsd10_noAGN$StellarMass_UB - devilsd10_noAGN$StellarMass_LB),
#       areas = devilsd10_noAGN$area,
# 
#       sm_bins = sm_bins,
#       vmax_bins = samples[[j]][,1],
#       vmax = vmax_combine[,2],
#       vmaxErr = vmax_combine[,3],
#       meanxErr = vmax_combine[,4],
#       errFloor = 0,
#       do_fit = TRUE,
#       do_optim = FALSE,
#       do_fit_quantiles = FALSE,
#       Niters = 100,
#       do_plot = TRUE,
#       add = TRUE,
#       ln.alpha = 0.05
#     ))
# 
#     fit$cosmic$Q50
#   }
#   dev.off()
#   sink()
#   file.remove("~/Documents/DustMassDensity/scripts/temp.txt")
# 
#   # return( colMedians(as.matrix(temp), na.rm = TRUE) )
#   return( median(temp, na.rm = TRUE) )
# 
# }
# cdmh_edd_corr = foreach(i = 1:(length(zbins)-1), .combine = c, .errorhandling = "pass") %do% {
#   
#   sink("~/Documents/DustMassDensity/scripts/temp.txt")
#   message(lbt_mids[i])
#   png(paste0("~/Documents/DustMassDensity/plots/edd_bias/lbt_", lbt_mids[i], ".png"), width = 7, height = 5, units = "in", res = 240)
#   GAMA = compute_mass_function(
#     zlo = zbins[i],
#     zhi = zbins[i+1],
#     z = sort_match_noAGN$z,
#     x = sort_match_noAGN$DustMass_50,
#     x_err = 0.5 * (sort_match_noAGN$DustMass_84 - sort_match_noAGN$DustMass_16),
#     areas = rep(217.54, length(sort_match_noAGN$z)),
#     sm_bins = sm_bins,
#     errFloor = 0.0,
#     do_fit = FALSE,
#     do_plot = TRUE,
#     add = FALSE,
#     pt.col = "purple"
#   )
#   DEVILS = compute_mass_function(
#     zlo = zbins[i],
#     zhi = zbins[i+1],
#     z = devilsd10_noAGN$z,
#     x = devilsd10_noAGN$dustmass.total,
#     x_err = 0.5 * (devilsd10_noAGN$dustmass.total_UB - devilsd10_noAGN$dustmass.total_LB),
#     areas = rep(1.5, length(devilsd10_noAGN$z)),
#     sm_bins = sm_bins,
#     errFloor = 0.0,
#     do_fit = FALSE,
#     do_plot = TRUE,
#     add = TRUE,
#     pt.col = "cornflowerblue"
#   )
#   
#   vmax_combine = foreach(j = 1:(length(sm_bins)-1) , .combine = rbind) %do% {
#     
#     if( GAMA$vmax$vmaxErr[j] <= ifelse(DEVILS$vmax$vmaxErr[j]==0, 999, DEVILS$vmax$vmaxErr[j]) ){
#       if( !GAMA$vmax$mlim[j] ){
#         c(DEVILS$vmax$x[j], DEVILS$vmax$vmax[j], DEVILS$vmax$vmaxErr[j], DEVILS$vmax$meanxErr[j])
#       }else{
#         if( (DEVILS$vmax$vmax[j] - GAMA$vmax$vmax[j])/(GAMA$vmax$vmaxErr[j] + 1e-323) > 5 ){
#           c(DEVILS$vmax$x[j], DEVILS$vmax$vmax[j], DEVILS$vmax$vmaxErr[j], DEVILS$vmax$meanxErr[j])
#         }else{
#           c(GAMA$vmax$x[j], GAMA$vmax$vmax[j], GAMA$vmax$vmaxErr[j], GAMA$vmax$meanxErr[j])
#         }
#       }
#     }else{
#       c(DEVILS$vmax$x[j], DEVILS$vmax$vmax[j], DEVILS$vmax$vmaxErr[j], DEVILS$vmax$meanxErr[j])
#     }
#     
#   }
#   
#   samples = foreach(j = 1:200) %do% {
#     vmax_bins_sample =  rnorm(length(vmax_combine[,1]), mean = vmax_combine[,1], sd = vmax_combine[,4])
#     vmax_bins_sample[is.nan(vmax_bins_sample)] = vmax_combine[,1][is.nan(vmax_bins_sample)]
#     
#     vmax_sample = rnorm(length(vmax_combine[,2]), mean = vmax_combine[,2], sd = vmax_combine[,3])
#     vmax_sample[is.nan(vmax_sample)] = vmax_combine[,2][is.nan(vmax_sample)]
#     cbind(vmax_bins_sample, vmax_sample)
#   }
#   
#   temp = foreach(j = 1:length(samples), .combine = c, .errorhandling = "remove") %do% {
#     
#     fit = suppressMessages(compute_mass_function(
#       zlo = zbins[i],
#       zhi = zbins[i+1],
#       
#       ## Not gonna use these
#       z = devilsd10_noAGN$z,
#       x = devilsd10_noAGN$dustmass.total,
#       x_err = 0.5 * (devilsd10_noAGN$dustmass.total_UB - devilsd10_noAGN$dustmass.total_LB),
#       areas = devilsd10_noAGN$area,
#       
#       sm_bins = sm_bins,
#       vmax_bins = samples[[j]][,1],
#       vmax = vmax_combine[,2],
#       vmaxErr = vmax_combine[,3],
#       meanxErr = vmax_combine[,4],
#       errFloor = 0,
#       do_fit = TRUE,
#       do_optim = FALSE,
#       do_fit_quantiles = FALSE,
#       Niters = 100,
#       do_plot = TRUE,
#       add = TRUE,
#       ln.alpha = 0.05
#     ))
#     
#     fit$cosmic$Q50
#   }
#   
#   foobar=compute_mass_function(
#     zlo = zbins[i],
#     zhi = zbins[i+1],
#     
#     ## Not gonna use these
#     z = devilsd10_noAGN$z,
#     x = devilsd10_noAGN$dustmass.total,
#     x_err = 0.5 * (devilsd10_noAGN$dustmass.total_UB - devilsd10_noAGN$dustmass.total_LB),
#     areas = devilsd10_noAGN$area,
#     
#     sm_bins = sm_bins,
#     vmax_bins = vmax_combine[,1],
#     vmax = vmax_combine[,2],
#     vmaxErr = vmax_combine[,3],
#     meanxErr = vmax_combine[,4],
#     errFloor = 0,
#     do_fit = TRUE,
#     do_fit_quantiles = FALSE,
#     Niters = 100,
#     do_plot = TRUE,
#     add = TRUE,
#     ln.alpha = 1.0,
#     pt.col = "red"
#   )
#   foobar$vmax$GAMA = GAMA$vmax$vmax
#   foobar$vmax$DEVILS = DEVILS$vmax$vmax
#   points(
#     10^foobar$vmax$x[fit$vmax$mlim == 1],
#     foobar$vmax$vmax[fit$vmax$mlim == 1],
#     pch = 16, col = "red"
#   )
#   dev.off()
#   
#   sink()
#   file.remove("~/Documents/DustMassDensity/scripts/temp.txt")
#   
#   # return( colMedians(as.matrix(temp), na.rm = TRUE) )
#   return( median(temp, na.rm = TRUE) )
# }

gc()
csmh_edd_corrV2 = foreach(i = 1:(length(zbins)-1), .errorhandling = "pass") %do% {
  
  sink("~/Documents/DustMassDensity/scripts/temp.txt")
  message(lbt_mids[i])
  png(paste0("~/Documents/DustMassDensity/plots/edd_bias_smf/lbt_", lbt_mids[i], ".png"), width = 7, height = 5, units = "in", res = 240)
  
  gama_x = sort_match_noAGN$StellarMass_50
  gama_err = 0.5 * (sort_match_noAGN$StellarMass_84 - sort_match_noAGN$StellarMass_16)
  # gama_err = gama_err/(log(10) * gama_x)
  # gama_x = log10(gama_x)
  
  devils_x = devilsd10_noAGN$StellarMass
  devils_err = 0.5 * (devilsd10_noAGN$StellarMass_UB - devilsd10_noAGN$StellarMass_LB)
  # devils_err = devils_err/(log(10) * devils_x)
  # devils_x = log10(devils_x)
  
  temp = c()
  # temp = foreach(j = 1:10, .combine = c, .errorhandling = "remove") %do% {
  
  for(j in 1:100){
    
    set.seed(j)
    gama_samples_ = rnorm(n = length(gama_x), mean = gama_x, sd = gama_err)
    
    GAMA = compute_mass_function(
      zlo = zbins[i],
      zhi = zbins[i+1],
      z = sort_match_noAGN$z,
      x = gama_samples_,
      x_err = gama_err,
      areas = rep(217.54, length(sort_match_noAGN$z)),
      sm_bins = sm_bins,
      errFloor = 0.0,
      do_fit = FALSE,
      do_plot = TRUE,
      add = j!=1,
      pt.col = "purple",
      ln.alpha = 0.1
    )
    
    devils_sample_ = rnorm(n = length(devils_x), mean = devils_x, sd = devils_err)
    DEVILS = compute_mass_function(
      zlo = zbins[i],
      zhi = zbins[i+1],
      z = devilsd10_noAGN$z,
      x = devils_sample_,
      x_err = devils_err,
      areas = rep(1.5, length(devilsd10_noAGN$z)),
      sm_bins = sm_bins,
      errFloor = 0.0,
      do_fit = FALSE,
      do_plot = TRUE,
      add = TRUE,
      pt.col = "cornflowerblue",
      ln.alpha = 0.1
    )
    
    vmax_combine = foreach(k = 1:(length(sm_bins)-1) , .combine = rbind) %do% {
      
      if( GAMA$vmax$vmaxErr[k] <= ifelse(DEVILS$vmax$vmaxErr[k]==0, 999, DEVILS$vmax$vmaxErr[k]) ){
        if( !GAMA$vmax$mlim[k] ){
          c(DEVILS$vmax$x[k], DEVILS$vmax$vmax[k], DEVILS$vmax$vmaxErr[k], DEVILS$vmax$meanxErr[k])
        }else{
          if( (DEVILS$vmax$vmax[k] - GAMA$vmax$vmax[k])/(GAMA$vmax$vmaxErr[k] + 1e-323) > 5 ){
            c(DEVILS$vmax$x[k], DEVILS$vmax$vmax[k], DEVILS$vmax$vmaxErr[k], DEVILS$vmax$meanxErr[k])
          }else{
            c(GAMA$vmax$x[k], GAMA$vmax$vmax[k], GAMA$vmax$vmaxErr[k], GAMA$vmax$meanxErr[k])
          }
        }
      }else{
        c(DEVILS$vmax$x[k], DEVILS$vmax$vmax[k], DEVILS$vmax$vmaxErr[k], DEVILS$vmax$meanxErr[k])
      }
    }
    
    fit = suppressMessages(compute_mass_function(
      zlo = zbins[i],
      zhi = zbins[i+1],
      
      ## Not gonna use these
      z = devilsd10_noAGN$z,
      x = devils_x,
      x_err = devils_err,
      areas = devilsd10_noAGN$area,
      sm_bins = sm_bins,
      
      vmax_bins = vmax_combine[,1],
      vmax = vmax_combine[,2],
      vmaxErr = vmax_combine[,3],
      meanxErr = vmax_combine[,4],
      
      errFloor = 0,
      do_fit = TRUE,
      do_optim = FALSE,
      do_fit_quantiles = FALSE,
      Niters = 100,
      do_plot = TRUE,
      add = TRUE,
      ln.alpha = 0.05
    ))
    
    temp = c(temp, fit$cosmic$Q50)
    # return( fit$cosmic$Q50 )
    # return( NULL )
  }
  
  lines(
    10^mdustvec,
    stellar_mass_density[[i]]$fit$Q50, 
    col = "red"
  )
  points(
    10^stellar_mass_density[[i]]$vmax$x[fit$vmax$mlim == 1],
    stellar_mass_density[[i]]$vmax$vmax[fit$vmax$mlim == 1],
    pch = 16, col = "red"
  )
  dev.off()
  
  sink()
  file.remove("~/Documents/DustMassDensity/scripts/temp.txt")
  
  # return( colMedians(as.matrix(temp), na.rm = TRUE) )
  return( temp )
}
cdmh_edd_corrV2 = foreach(i = 1:(length(zbins)-1), .errorhandling = "pass") %do% {
  
  sink("~/Documents/DustMassDensity/scripts/temp.txt")
  message(lbt_mids[i])
  png(paste0("~/Documents/DustMassDensity/plots/edd_bias/lbt_", lbt_mids[i], ".png"), width = 7, height = 5, units = "in", res = 240)
  
  gama_x = sort_match_noAGN$DustMass_50
  gama_err = 0.5 * (sort_match_noAGN$DustMass_84 - sort_match_noAGN$DustMass_16)
  # gama_err = gama_err/(log(10) * gama_x)
  # gama_x = log10(gama_x)
  
  devils_x = devilsd10_noAGN$dustmass.total
  devils_err = 0.5 * (devilsd10_noAGN$dustmass.total_UB - devilsd10_noAGN$dustmass.total_LB)
  # devils_err = devils_err/(log(10) * devils_x)
  # devils_x = log10(devils_x)

  temp = c()
  # temp = foreach(j = 1:10, .combine = c, .errorhandling = "remove") %do% {

  for(j in 1:100){
    
    set.seed(j)
    gama_samples_ = rnorm(n = length(gama_x), mean = gama_x, sd = gama_err)
    
    GAMA = compute_mass_function(
      zlo = zbins[i],
      zhi = zbins[i+1],
      z = sort_match_noAGN$z,
      x = gama_samples_,
      x_err = gama_err,
      areas = rep(217.54, length(sort_match_noAGN$z)),
      sm_bins = sm_bins,
      errFloor = 0.0,
      do_fit = FALSE,
      do_plot = TRUE,
      add = j!=1,
      pt.col = "purple",
      ln.alpha = 0.1
    )
    
    devils_sample_ = rnorm(n = length(devils_x), mean = devils_x, sd = devils_err)
    DEVILS = compute_mass_function(
      zlo = zbins[i],
      zhi = zbins[i+1],
      z = devilsd10_noAGN$z,
      x = devils_sample_,
      x_err = devils_err,
      areas = rep(1.5, length(devilsd10_noAGN$z)),
      sm_bins = sm_bins,
      errFloor = 0.0,
      do_fit = FALSE,
      do_plot = TRUE,
      add = TRUE,
      pt.col = "cornflowerblue",
      ln.alpha = 0.1
    )
    
    vmax_combine = foreach(k = 1:(length(sm_bins)-1) , .combine = rbind) %do% {
      
      if( GAMA$vmax$vmaxErr[k] <= ifelse(DEVILS$vmax$vmaxErr[k]==0, 999, DEVILS$vmax$vmaxErr[k]) ){
        if( !GAMA$vmax$mlim[k] ){
          c(DEVILS$vmax$x[k], DEVILS$vmax$vmax[k], DEVILS$vmax$vmaxErr[k], DEVILS$vmax$meanxErr[k])
        }else{
          if( (DEVILS$vmax$vmax[k] - GAMA$vmax$vmax[k])/(GAMA$vmax$vmaxErr[k] + 1e-323) > 5 ){
            c(DEVILS$vmax$x[k], DEVILS$vmax$vmax[k], DEVILS$vmax$vmaxErr[k], DEVILS$vmax$meanxErr[k])
          }else{
            c(GAMA$vmax$x[k], GAMA$vmax$vmax[k], GAMA$vmax$vmaxErr[k], GAMA$vmax$meanxErr[k])
          }
        }
      }else{
        c(DEVILS$vmax$x[k], DEVILS$vmax$vmax[k], DEVILS$vmax$vmaxErr[k], DEVILS$vmax$meanxErr[k])
      }
    }
    
    fit = suppressMessages(compute_mass_function(
      zlo = zbins[i],
      zhi = zbins[i+1],

      ## Not gonna use these
      z = devilsd10_noAGN$z,
      x = devilsd10_noAGN$dustmass.total,
      x_err = 0.5 * (devilsd10_noAGN$dustmass.total_UB - devilsd10_noAGN$dustmass.total_LB),
      areas = devilsd10_noAGN$area,
      sm_bins = sm_bins,

      vmax_bins = vmax_combine[,1],
      vmax = vmax_combine[,2],
      vmaxErr = vmax_combine[,3],
      meanxErr = vmax_combine[,4],

      errFloor = 0,
      do_fit = TRUE,
      do_optim = FALSE,
      do_fit_quantiles = FALSE,
      Niters = 100,
      do_plot = TRUE,
      add = TRUE,
      ln.alpha = 0.05
    ))

    temp = c(temp, fit$cosmic$Q50)
    # return( fit$cosmic$Q50 )
    # return( NULL )
  }
  
  lines(
    10^mdustvec,
    dust_mass_density[[i]]$fit$Q50, 
    col = "red"
  )
  points(
    10^dust_mass_density[[i]]$vmax$x[fit$vmax$mlim == 1],
    dust_mass_density[[i]]$vmax$vmax[fit$vmax$mlim == 1],
    pch = 16, col = "red"
  )
  dev.off()
  
  sink()
  file.remove("~/Documents/DustMassDensity/scripts/temp.txt")
  
  # return( colMedians(as.matrix(temp), na.rm = TRUE) )
  return( temp )
}

gc()
csmh = data.frame(foreach(i = 1:length(stellar_mass_density), .combine = bind_rows) %do% {
  stellar_mass_density[[i]]$cosmic
})
cdmh = data.frame(foreach(i = 1:length(dust_mass_density), .combine = bind_rows) %do% {
  dust_mass_density[[i]]$cosmic
})
cdmh_wAGN = data.frame(foreach(i = 1:length(dust_mass_density_wAGN), .combine = bind_rows) %do% {
  dust_mass_density_wAGN[[i]]$cosmic
})
# cdmh_devils = data.frame(foreach(i = 1:length(dust_mass_density_devils_only), .combine = bind_rows) %do% {
#   dust_mass_density_devils_only[[i]]$cosmic
# })
# cdmh_gama = data.frame(foreach(i = 1:length(dust_mass_density_gama_only), .combine = bind_rows) %do% {
#   dust_mass_density_gama_only[[i]]$cosmic
# })

print(
  log10(sapply(csmh_edd_corrV2, sd)) 
)
print(
  log10(cdmh_edd_corrV2) - log10(cdmh$Q50)
)


# cdlh = data.frame(foreach(i = 1:length(luminosity_density), .combine = bind_rows) %do% {
#   luminosity_density[[i]]$cosmic
# })
# cdlh_wAGN = data.frame(foreach(i = 1:length(luminosity_density_wAGN), .combine = bind_rows) %do% {
#   luminosity_density_wAGN[[i]]$cosmic
# })

dsilva25 = data.frame(fread("~/Documents/DustMassDensity/data/literature_evo/csfh/DSilva25_CSFH_CAGNH_fit.csv"))
dsilva_csmh_func = approxfun(
  cosdistTravelTime(dsilva25$z, ref = "Planck18"), 
  10^dsilva25$CSMHQ50
)
LL_csmh = function(p){
  
  -1*sum(dnorm(
    x = log10( csmh$Q50 ),
    mean = p[1] + log10(dsilva_csmh_func(lbt_mids)),
    sd = csmh$ERR/(log(10) * csmh$Q50), 
    log = TRUE
  ))
  
}
norm_csmh = optimise(f = LL_csmh, interval = c(-2,1))

LSS_corr = dsilva_csmh_func(lbt_mids) / csmh$Q50 * 10^norm_csmh$minimum

h5file = '~/Documents/DustMassDensity/data/all_data.h5'

h5createFile(h5file)

h5delete(h5file, "zmids")
h5delete(h5file, "zbins")
h5delete(h5file, "lbtmids")
h5delete(h5file, "lbtbins")
h5delete(h5file, "LSSCorrection")
h5delete(h5file, "EddingtonCorrectionSMF")
h5delete(h5file, "EddingtonCorrection")

h5write(obj = zmids, h5file, name = "zmids")
h5write(obj = zbins, h5file, name = "zbins")
h5write(obj = lbt_mids, h5file, name = "lbtmids")
h5write(obj = lbt_bins, h5file, name = "lbtbins")
h5write(obj = LSS_corr, h5file, name = "LSSCorrection")
h5write(obj = csmh_edd_corrV2/csmh$Q50, h5file, name = "EddingtonCorrectionSMF")
h5write(obj = cdmh_edd_corrV2/cdmh$Q50, h5file, name = "EddingtonCorrection")


h5delete(h5file, "cosmic")
h5createGroup(h5file, "cosmic")
h5write(obj = csmh, file = h5file, name = "cosmic/Mstar")
h5write(obj = cdmh, file = h5file, name = "cosmic/Mdust")
h5write(obj = cdmh, file = h5file, name = "cosmic/Mdust_gama")
h5write(obj = cdmh, file = h5file, name = "cosmic/Mdust_devils")
h5write(obj = cdmh_wAGN, file = h5file, name = "cosmic/MdustwAGN")
# h5write(obj = cdlh, file = h5file, name = "cosmic/Ldust")
# h5write(obj = cdlh_wAGN, file = h5file, name = "cosmic/LdustwAGN")

h5delete(h5file, "vmax")
h5createGroup(h5file, "vmax")
h5createGroup(h5file, "vmax/Mstar")
for(i in 1:length(stellar_mass_density)){
  h5write(
    obj = stellar_mass_density[[i]]$vmax, 
    file = h5file,
    name = paste0("vmax/Mstar/zb",i)
  )
}
h5createGroup(h5file, "vmax/Mdust")
for(i in 1:length(dust_mass_density)){
  h5write(
    obj = dust_mass_density[[i]]$vmax, 
    file = h5file,
    name = paste0("vmax/Mdust/zb",i)
  )
}
h5createGroup(h5file, "vmax/MdustwAGN")
for(i in 1:length(dust_mass_density_wAGN)){
  h5write(
    obj = dust_mass_density_wAGN[[i]]$vmax, 
    file = h5file,
    name = paste0("vmax/MdustwAGN/zb",i)
  )
}
# h5createGroup(h5file, "vmax/Ldust")
# for(i in 1:length(luminosity_density)){
#   h5write(
#     obj = luminosity_density[[i]]$vmax, 
#     file = h5file,
#     name = paste0("vmax/Ldust/zb",i)
#   )
# }
# h5createGroup(h5file, "vmax/LdustwAGN")
# for(i in 1:length(luminosity_density_wAGN)){
#   h5write(
#     obj = luminosity_density_wAGN[[i]]$vmax, 
#     file = h5file,
#     name = paste0("vmax/LdustwAGN/zb",i)
#   )
# }

h5delete(h5file, "fit")
h5createGroup(h5file, "fit")
h5createGroup(h5file, "fit/Mstar")
for(i in 1:length(stellar_mass_density)){
  h5write(
    obj = stellar_mass_density[[i]]$fit, 
    file = h5file,
    name = paste0("fit/Mstar/zb",i)
  )
}
h5createGroup(h5file, "fit/Mdust")
for(i in 1:length(dust_mass_density)){
  h5write(
    obj = dust_mass_density[[i]]$fit, 
    file = h5file,
    name = paste0("fit/Mdust/zb",i)
  )
}
h5createGroup(h5file, "fit/MdustwAGN")
for(i in 1:length(dust_mass_density_wAGN)){
  h5write(
    obj = dust_mass_density_wAGN[[i]]$fit, 
    file = h5file,
    name = paste0("fit/MdustwAGN/zb",i)
  )
}
# h5createGroup(h5file, "fit/Ldust")
# for(i in 1:length(luminosity_density)){
#   h5write(
#     obj = luminosity_density[[i]]$fit, 
#     file = h5file,
#     name = paste0("fit/Ldust/zb",i)
#   )
# }
# h5createGroup(h5file, "fit/LdustwAGN")
# for(i in 1:length(luminosity_density_wAGN)){
#   h5write(
#     obj = luminosity_density_wAGN[[i]]$fit, 
#     file = h5file,
#     name = paste0("fit/LdustwAGN/zb",i)
#   )
# }


## Dust mass density
driver18_cdmh_raw = data.frame(fread("~/Documents/DustMassDensity/data/literature_evo/cdmh/driver18_raw.csv", skip = 3))
driver18_cdmh = data.frame(
    "lbt" = driver18_cdmh_raw$V1,
    "cdmh" = driver18_cdmh_raw$V3,
    "err_pois" = as.numeric(str_remove(driver18_cdmh_raw$V5, "±")),
    "err_cv" = as.numeric(str_remove(driver18_cdmh_raw$V7, "±")),
    "err_AGN" = as.numeric(str_remove(driver18_cdmh_raw$V8, "±"))
)
driver18_cdmh$err_pois = sqrt( (10^driver18_cdmh$cdmh * log(10) * driver18_cdmh$err_pois)^2 )
driver18_cdmh$err_cv = sqrt( (10^driver18_cdmh$cdmh * log(10) * driver18_cdmh$err_cv)^2 )
driver18_cdmh$err_AGN = sqrt( (10^driver18_cdmh$cdmh * log(10) * driver18_cdmh$err_AGN)^2 )
driver18_cdmh$err = sqrt( driver18_cdmh$err_pois^2 + driver18_cdmh$err_cv^2 + driver18_cdmh$err_AGN^2 )
driver18_cdmh$cdmh = 10^driver18_cdmh$cdmh
fwrite(driver18_cdmh, "~/Documents/DustMassDensity/data/literature_evo/cdmh/driver18.csv")

# thorne21 = fread("~/Documents/JWSTCATALOGUE/Data/literature_evo/csmh/thorne21.csv")
# driver18_csmh = fread("~/Documents/JWSTCATALOGUE/Data/literature_evo/csmh/driver2018.csv")
# 
# driver18_cdmh = data.frame(
#   "lbt" = c(0.85, 1.52, 2.16, 2.90, 3.65, 4.35, 5.11, 5.86, 6.59, 7.36, 8.11, 8.82, 9.50),
#   "cdmh" = 10^c(5.17, 5.25, 5.24, 5.24, 5.27, 5.18, 5.20, 5.32, 5.39, 5.65, 5.54, 5.55, 5.34)
# )
# magplot(
#   lbt_mids, 
#   cdmh$Q50, 
#   log = "y", 
#   pch = 16, 
#   xlab = "Lookback time [Gyr]",
#   ylab = "Cosmic Mdust density [Msun yr^-1 Mpc^-3]", 
#   ylim = c(5e2, 1e9),
#   xlim = c(0, 12.5), 
#   col = "blue"
# )
# magerr(
#   lbt_mids, 
#   cdmh$Q50, 
#   ylo = cdmh$Q50 - cdmh$Q16,
#   yhi = cdmh$Q84 - cdmh$Q50,
#   col = "blue"
# )
# points(
#   lbt_mids, 
#   cdmh_wAGN$Q50,
#   pch = 18, col = "purple"
# )
# magerr(
#   lbt_mids, 
#   cdmh_wAGN$Q50, 
#   ylo = cdmh_wAGN$Q50 - cdmh_wAGN$Q16,
#   yhi = cdmh_wAGN$Q84 - cdmh_wAGN$Q50,
#   col = "purple"
# )
# md14 = log10( 0.015 * ( (1+zvec)^2.7/(1 + ((1+zvec)/2.9)^4.6) ) ) + log10(0.63)
# lines(
#   lbtvec, 10^md14 * 2e7
# )
# points(
#   cosdistTravelTime(driver18_csmh$z, ref = "737"), 
#   10^driver18_csmh$csmd, 
#   pch = 8, 
#   col = "cornflowerblue"
# )
# points(
#   lbt_mids, 
#   csmh$Q50, 
#   pch = 16
# )
# magerr(
#   lbt_mids, 
#   csmh$Q50, 
#   ylo = csmh$Q50 - csmh$Q16,
#   yhi = csmh$Q84 - csmh$Q50,
# )
# points(
#   driver18_cdmh$lbt,
#   driver18_cdmh$cdmh, 
#   pch = 15, 
#   col = "cornflowerblue"
# )
# legend(
#   x = "bottomleft", 
#   pch = c(16, 15, NA, 8),
#   lty = c(NA, NA, 1, NA),
#   col = c("black", "cornflowerblue", "black", "black"), 
#   legend = c("GAMA+DEVILS", "Dust mass density \n Driver+18", "Madau&Dickinson (2014) \n scaled up", "Stellar mass density \n Driver+18")
# )
# 
# magplot(
#   lbt_mids, smf_evol[,5]
# )




# tot_cat = data.frame(
#   "z" = c(sort_match_noAGN$z, devilsd10_noAGN$z),
#   "Mstar" = c(sort_match_noAGN$StellarMass_50, devilsd10_noAGN$StellarMass), 
#   "MstarErr" = c(0.5*(sort_match_noAGN$StellarMass_84 - sort_match_noAGN$StellarMass_16), 
#                  0.5*(devilsd10_noAGN$StellarMass_UB - devilsd10_noAGN$StellarMass_LB)),
#   "Mdust" = c(sort_match_noAGN$DustMass_50, devilsd10_noAGN$dustmass.total),
#   "MdustErr" = c(0.5*(sort_match_noAGN$DustMass_84 - sort_match_noAGN$DustMass_16), 
#                  0.5*(devilsd10_noAGN$dustmass.total_UB - devilsd10_noAGN$dustmass.total_LB)), 
#   "surveyID" = c(rep("GAMA", dim(sort_match_noAGN)[1]), rep("DEVILS", dim(devilsd10_noAGN)[1]))
# )

# tot_cat$areas = rep(1.5, dim(tot_cat)[1])
# tot_cat$areas[tot_cat$surveyID == "GAMA"] = 217.54
# 
# maghist(
#   log(1+sort_match_AGN$z), log = "y", col = "purple"
# )
# maghist(
#   log(1+devilsd10_noAGN$cat$z), log = "y", add =TRUE, col = "cornflowerblue"
# )
# 
# dmf1 = compute_mass_function(
#   zlo = 0.02,
#   zhi = 0.08,
#   z = sort_match_noAGN$z,
#   x = sort_match_noAGN$DustMass_50,
#   areas = sort_match_noAGN$area,
#   sm_bins = seq(6.5, 12.5, 0.5), 
#   errFloor = 0.1,
#   vmax = NULL,
#   vmaxErr = NULL,
#   do_plot = TRUE,
#   do_fit = TRUE,
#   add = FALSE
# )
# dmf1 = compute_mass_function(
#   zlo = 0.02,
#   zhi = 0.08,
#   z = devilsd10_noAGN$z,
#   x = devilsd10_noAGN$dustmass.total,
#   areas = devilsd10_noAGN$area,
#   sm_bins = seq(6.5, 12.5, 0.5),
#   errFloor = 0.1,
#   vmax = NULL,
#   vmaxErr = NULL,
#   do_plot = TRUE,
#   do_fit = TRUE,
#   add = TRUE
# )
# dmf1 = compute_mass_function(
#   zlo = 0.02,
#   zhi = 0.08,
#   z = tot_cat$z,
#   x = tot_cat$Mdust,
#   areas = tot_cat$areas,
#   sm_bins = seq(6.5, 12.5, 0.5),
#   vmax = NULL,
#   vmaxErr = NULL,
#   do_plot = TRUE,
#   do_fit = FALSE,
#   add = TRUE
# )
# gama_dust_mass_lim = foreach(i = 1:(length(zbins)-1), .combine = "c") %do% {
#   zidx = sort_match_noAGN$z >= zbins[i] & sort_match_noAGN$z < zbins[i+1]
#   if(sum(zidx) > 0){
#     quantile(sort_match_noAGN$DustMass_Lim[zidx], probs = 0.95)
#   }else{
#     1e15
#   }
# }
# devils_dust_mass_lim = foreach(i = 1:(length(zbins)-1), .combine = "c") %do% {
#   zidx = sort_match_noAGN$z >= zbins[i] & sort_match_noAGN$z < zbins[i+1]
#   if(sum(zidx) > 0){
#     quantile(devilsd10_noAGN$cat$DustMass_Lim[zidx], probs = 0.95)
#   }else{
#     1e15
#   }
# }
# tot_mass_lim = foreach(i = 1:(length(zbins)-1), .combine = "c") %do% {
#   zidx = tot_cat$z >= zbins[i] & tot_cat$z < zbins[i+1]
#   if(sum(zidx) > 0){
#     quantile(tot_cat$MdustLim[zidx], probs = 0.95)
#   }else{
#     1e15
#   }
# }
# 
# gama_Mdust_completeness_spl = smooth.spline(
#   x = zmids[gama_dust_mass_lim<1e15], 
#   log10(gama_dust_mass_lim[gama_dust_mass_lim<1e15]), 
#   df = 5
# )
# gama_Mdust_completeness_predict = predict(gama_Mdust_completeness_spl, zvec)
# gama_Mdust_completeness_func = approxfun(gama_Mdust_completeness_predict$x, 10^gama_Mdust_completeness_predict$y)
# 
# devils_Mdust_completeness_spl = smooth.spline(
#   x = zmids[devils_dust_mass_lim<1e15], 
#   log10(devils_dust_mass_lim[devils_dust_mass_lim<1e15]), 
#   df = 5
# )
# devils_Mdust_completeness_predict = predict(devils_Mdust_completeness_spl, zvec)
# devils_Mdust_completeness_func = approxfun(devils_Mdust_completeness_predict$x, 10^devils_Mdust_completeness_predict$y)
# 
# tot_Mdust_completeness_spl = smooth.spline(
#   x = zmids[tot_mass_lim<1e15], 
#   log10(tot_mass_lim[tot_mass_lim<1e15]), 
#   df = 5
# )
# tot_Mdust_completeness_predict = predict(tot_Mdust_completeness_spl, zvec)
# tot_Mdust_completeness_func = approxfun(tot_Mdust_completeness_predict$x, 10^tot_Mdust_completeness_predict$y)
# 
# foobar = foreach(i = 1:(length(zbins)-1), .combine = "c") %do% {
#   
#   vol = 4*pi/3 * (cosdistCoDist(z = zbins[i+1], ref = '737')^3 - cosdistCoDist(z = zbins[i], ref = '737')^3)
#   
#   gama_lim = gama_Mdust_completeness_func(zbins[i+1])
#   devils_lim = devils_Mdust_completeness_func(zbins[i+1])
#   
#   gama_Mdust_trim = log10( sort_match_noAGN$DustMass_50[sort_match_noAGN$DustMass_50 >= gama_lim & sort_match_noAGN$z >= zbins[i] & sort_match_noAGN$z < zbins[i+1]] )
#   devils_Mdust_trim = log10( devilsd10_noAGN$cat$dustmass.total[devilsd10_noAGN$cat$dustmass.total >= devils_lim & devilsd10_noAGN$cat$z >= zbins[i] & devilsd10_noAGN$cat$z < zbins[i+1]] )
#   # 
#   gama_Mdust = log10( sort_match_noAGN$DustMass_50[sort_match_noAGN$z >= zbins[i] & sort_match_noAGN$z < zbins[i+1]] )
#   devils_Mdust = log10( devilsd10_noAGN$cat$dustmass.total[devilsd10_noAGN$cat$z >= zbins[i] & devilsd10_noAGN$cat$z < zbins[i+1]] )
#   
#   gama_vol = (217.54) / (4*pi*(180/pi)^2) * vol
#   gama_vol = rep(gama_vol, length(gama_Mdust))
#   devils_vol = (1.5) / (4*pi*(180/pi)^2) * vol
#   devils_vol = rep(devils_vol, length(devils_Mdust))
#   
#   # gama_vol_trim = (217.54) / (4*pi*(180/pi)^2) * vol
#   # gama_vol_trim = rep(gama_vol_trim, length(gama_Mdust_trim))
#   # devils_vol_trim = (1.5) / (4*pi*(180/pi)^2) * vol
#   # devils_vol_trim = rep(devils_vol_trim, length(devils_Mdust_trim))
#   # 
#   # Mdust = c(gama_Mdust, devils_Mdust)
#   # vols = c(gama_vol, devils_vol)
#   # 
#   # Mdust_trim = c(gama_Mdust_trim, devils_Mdust_trim)
#   # vols_trim = c(gama_vol_trim, devils_vol_trim)
#   
#   gama_vmax = foreach(j = 1:length(sm_mids), .combine = "c") %do% {
#     midx = gama_Mdust >= sm_bins[j] & gama_Mdust < sm_bins[j+1]
#     sum(gama_vol[midx]^-1)
#   }
#   devils_vmax = foreach(j = 1:length(sm_mids), .combine = "c") %do% {
#     midx = devils_Mdust >= sm_bins[j] & devils_Mdust < sm_bins[j+1]
#     sum(devils_vol[midx]^-1)
#   }
#   
#   gama_vmaxErr = foreach(j = 1:length(sm_mids), .combine = "c") %do% {
#     midx = gama_Mdust >= sm_bins[j] & gama_Mdust < sm_bins[j+1]
#     sqrt(sum(gama_vol[midx]^-2))
#   }
#   devils_vmaxErr = foreach(j = 1:length(sm_mids), .combine = "c") %do% {
#     midx = devils_Mdust >= sm_bins[j] & devils_Mdust < sm_bins[j+1]
#     sqrt(sum(devils_vol[midx]^-2))
#   }
#   
#   magplot(
#     10^sm_mids, gama_vmax, log = "xy", pch = 1, xlim = 10^c(2.5, 15.5), ylim = c(1e-10, 1), col = "purple"
#   )
#   magerr(
#     10^sm_mids, gama_vmax, ylo = gama_vmaxErr, col = "purple"
#   )
#   abline(
#     v = gama_lim
#   )
#   points(
#     10^sm_mids, devils_vmax, pch = 16, col = "cornflowerblue"
#   )
#   abline(
#     v = devils_lim
#   )
# }
# luminosity_density = foreach(i = 1:(length(zbins)-1)) %do% {
#   
#   GAMA = compute_mass_function(
#     zlo = zbins[i],
#     zhi = zbins[i+1],
#     z = sort_match_noAGN$z,
#     x = sort_match_noAGN$DustLum_50,
#     areas = rep(217.54, length(sort_match_noAGN$z)),
#     sm_bins = sm_bins,
#     errFloor = 0.1,
#     do_fit = FALSE,
#     do_plot = TRUE,
#     add = FALSE,
#     pt.col = "purple"
#   )
#   DEVILS = compute_mass_function(
#     zlo = zbins[i],
#     zhi = zbins[i+1],
#     z = devilsd10_noAGN$z,
#     x = devilsd10_noAGN$dustlum.total,
#     areas = rep(1.5, length(devilsd10_noAGN$z)),
#     sm_bins = sm_bins,
#     errFloor = 0.1,
#     do_fit = FALSE,
#     do_plot = TRUE,
#     add = TRUE,
#     pt.col = "cornflowerblue"
#   )
#   
#   vmax_combine = foreach(j = 1:(length(sm_bins)-1) , .combine = rbind) %do% {
#     
#     if( GAMA$vmax$vmaxErr[j] <= ifelse(DEVILS$vmax$vmaxErr[j]==0, 999, DEVILS$vmax$vmaxErr[j]) ){
#       if( !GAMA$vmax$mlim[j] ){
#         c(DEVILS$vmax$x[j], DEVILS$vmax$vmax[j], DEVILS$vmax$vmaxErr[j])
#       }else{
#         if( (DEVILS$vmax$vmax[j] - GAMA$vmax$vmax[j])/(GAMA$vmax$vmaxErr[j] + 1e-323) > 5 ){
#           c(DEVILS$vmax$x[j], DEVILS$vmax$vmax[j], DEVILS$vmax$vmaxErr[j])
#         }else{
#           c(GAMA$vmax$x[j], GAMA$vmax$vmax[j], GAMA$vmax$vmaxErr[j])
#         }
#       }
#     }else{
#       c(DEVILS$vmax$x[j], DEVILS$vmax$vmax[j], DEVILS$vmax$vmaxErr[j])
#     }
#     
#   }
#   
#   fit = compute_mass_function(
#     zlo = zbins[i],
#     zhi = zbins[i+1],
#     
#     ## Not gonna use these 
#     z = devilsd10_noAGN$z,
#     x = devilsd10_noAGN$dustmass.total,
#     areas = devilsd10_noAGN$area,
#     
#     sm_bins = sm_bins,
#     vmax_bins = vmax_combine[,1],
#     vmax = vmax_combine[,2],
#     vmaxErr = vmax_combine[,3],
#     errFloor = 0,
#     do_fit = TRUE,
#     do_plot = TRUE,
#     add = TRUE
#   )
#   fit$vmax$GAMA = GAMA$vmax$vmax
#   fit$vmax$DEVILS = DEVILS$vmax$vmax
#   points(
#     10^fit$vmax$x[fit$vmax$mlim == 1],
#     fit$vmax$vmax[fit$vmax$mlim == 1], 
#     pch = 16, col = "red"
#   )
#   return(fit)
# }
# luminosity_density_wAGN = foreach(i = 1:(length(zbins)-1)) %do% {
#   
#   GAMA = compute_mass_function(
#     zlo = zbins[i],
#     zhi = zbins[i+1],
#     z = sort_match_AGN$z,
#     x = sort_match_AGN$dustlum.total,
#     areas = rep(217.54, length(sort_match_AGN$z)),
#     sm_bins = sm_bins,
#     errFloor = 0.1,
#     do_fit = FALSE,
#     do_plot = TRUE,
#     add = FALSE,
#     pt.col = "purple"
#   )
#   DEVILS = compute_mass_function(
#     zlo = zbins[i],
#     zhi = zbins[i+1],
#     z = devilsd10_AGN$z,
#     x = devilsd10_AGN$dustlum.total,
#     areas = rep(1.5, length(devilsd10_AGN$z)),
#     sm_bins = sm_bins,
#     errFloor = 0.1,
#     do_fit = FALSE,
#     do_plot = TRUE,
#     add = TRUE,
#     pt.col = "cornflowerblue"
#   )
#   
#   vmax_combine = foreach(j = 1:(length(sm_bins)-1) , .combine = rbind) %do% {
#     
#     if( GAMA$vmax$vmaxErr[j] <= ifelse(DEVILS$vmax$vmaxErr[j]==0, 999, DEVILS$vmax$vmaxErr[j]) ){
#       if( !GAMA$vmax$mlim[j] ){
#         c(DEVILS$vmax$x[j], DEVILS$vmax$vmax[j], DEVILS$vmax$vmaxErr[j])
#       }else{
#         if( (DEVILS$vmax$vmax[j] - GAMA$vmax$vmax[j])/(GAMA$vmax$vmaxErr[j] + 1e-323) > 5 ){
#           c(DEVILS$vmax$x[j], DEVILS$vmax$vmax[j], DEVILS$vmax$vmaxErr[j])
#         }else{
#           c(GAMA$vmax$x[j], GAMA$vmax$vmax[j], GAMA$vmax$vmaxErr[j])
#         }
#       }
#     }else{
#       c(DEVILS$vmax$x[j], DEVILS$vmax$vmax[j], DEVILS$vmax$vmaxErr[j])
#     }
#     
#   }
#   
#   fit = compute_mass_function(
#     zlo = zbins[i],
#     zhi = zbins[i+1],
#     
#     ## Not gonna use these 
#     z = devilsd10_noAGN$z,
#     x = devilsd10_noAGN$dustmass.total,
#     areas = devilsd10_noAGN$area,
#     
#     sm_bins = sm_bins,
#     vmax_bins = vmax_combine[,1],
#     vmax = vmax_combine[,2],
#     vmaxErr = vmax_combine[,3],
#     errFloor = 0,
#     do_fit = TRUE,
#     do_plot = TRUE,
#     add = TRUE
#   )
#   fit$vmax$GAMA = GAMA$vmax$vmax
#   fit$vmax$DEVILS = DEVILS$vmax$vmax
#   points(
#     10^fit$vmax$x[fit$vmax$mlim == 1],
#     fit$vmax$vmax[fit$vmax$mlim == 1], 
#     pch = 16, col = "red"
#   )
#   return(fit)
# }

# dust_mass_density_devils_only = foreach(i = 1:(length(zbins)-1), .errorhandling = "remove") %do% {
#   
#   DEVILS = compute_mass_functionV2(
#     zlo = zbins[i],
#     zhi = zbins[i+1],
#     z = devilsd10_noAGN$z,
#     x = devilsd10_noAGN$dustmass.total,
#     x_err = 0.5 * (devilsd10_noAGN$dustmass.total_UB - devilsd10_noAGN$dustmass.total_LB),
#     areas = rep(1.5, length(devilsd10_noAGN$z)),
#     sm_bins = sm_bins,
#     errFloor = 0.0,
#     do_fit = TRUE,
#     do_plot = TRUE,
#     add = FALSE,
#     pt.col = "cornflowerblue"
#   )
#   
#   return(DEVILS)
# }
# dust_mass_density_gama_only = foreach(i = 1:(length(zbins)-1), .errorhandling = "remove") %do% {
#   
#   GAMA = compute_mass_functionV2(
#     zlo = zbins[i],
#     zhi = zbins[i+1],
#     z = sort_match_noAGN$z,
#     x = sort_match_noAGN$DustMass_50,
#     x_err = 0.5 * (sort_match_noAGN$DustMass_84 - sort_match_noAGN$DustMass_16),
#     areas = rep(217.54, length(sort_match_noAGN$z)),
#     sm_bins = sm_bins,
#     errFloor = 0.0,
#     do_fit = TRUE,
#     do_plot = TRUE,
#     add = FALSE,
#     pt.col = "purple"
#   )
#   
#   return(GAMA)
# }

stellar_wave = BC03hr$Wave
stellar_spectrum = BC03hr$Zspec[[1]][1,]
foobar = Dale_interp(Dale = Dale_NormTot)

magplot(
  foobar$Wave, foobar$Aspec, log = "xy"
)
lines(
  stellar_wave, stellar_spectrum  
)

fluxBC03=Lum2Flux(BC03lr$Wave, BC03lr$Zspec[[5]][161,]*1e10)
birthBC03_atten=CF_birth_atten(fluxBC03[,1],fluxBC03[,2], tau = 1e-3)
screenBC03_atten=CF_screen_atten(fluxBC03[,1],fluxBC03[,2], tau = 1e-3)
totalatten=birthBC03_atten$total_atten+screenBC03_atten$total_atten

# Here we show 3 pure greybodies, and then 3 Dale templates normalised by SFR with
# differing amounts of AGN contribution. Note even the pure SFR template is not quite
# greybody (there is an excess in the MIR and very-FIR).

plot(fluxBC03, log='xy', xlab=BC03lr$Labels$Wavelab, ylab='Flux / erg/s/cm^2/Angstrom',
     type='l', col='red', xlim=c(1e2,1e7))
lines(fluxBC03[,1], (birthBC03_atten$flux+screenBC03_atten$flux)/2, col='black')
lines(10^seq(5,7,by=0.01), greybody_norm(10^seq(5,7,by=0.01), beta=1, norm=totalatten),
      col='brown', lty=1)

lines(10^seq(5,7,by=0.01), greybody_norm(10^seq(5,7,by=0.01), beta=1.5, norm=totalatten),
      col='brown', lty=2)
lines(10^seq(5,7,by=0.01), greybody_norm(10^seq(5,7,by=0.01), beta=2, norm=totalatten),
      col='brown', lty=3)
lines(Dale_NormSFR$Wave, Dale_NormSFR$Aspec[[1]][24,]*totalatten, col='purple')
lines(Dale_NormSFR$Wave, Dale_NormSFR$Aspec[[11]][24,]*totalatten, col='orange')
lines(Dale_NormSFR$Wave, Dale_NormSFR$Aspec[[20]][24,]*totalatten, col='green')


stellar_wave = BC03lr$Wave
fluxBC03=BC03lr$Zspec[[5]][161,]
birthBC03_atten=CF_birth_atten(stellar_wave, fluxBC03)
screenBC03_atten=CF_screen_atten(stellar_wave, fluxBC03)
totalatten=birthBC03_atten$total_atten+screenBC03_atten$total_atten
plot(stellar_wave, fluxBC03, log='xy', xlab=BC03lr$Labels$Wavelab, ylab='Flux / erg/s/cm^2/Angstrom',
     type='l', col='red', xlim=c(1e2,1e7))
lines(stellar_wave, (birthBC03_atten$flux+screenBC03_atten$flux)/2, col='black')
lines(10^seq(5,7,by=0.01), greybody_norm(10^seq(5,7,by=0.01), beta=1, norm=totalatten),
      col='brown', lty=1)

lines(10^seq(5,7,by=0.01), greybody_norm(10^seq(5,7,by=0.01), beta=1.5, norm=totalatten),
      col='brown', lty=2)
lines(10^seq(5,7,by=0.01), greybody_norm(10^seq(5,7,by=0.01), beta=2, norm=totalatten),
      col='brown', lty=3)
lines(Dale_NormSFR$Wave, Dale_NormSFR$Aspec[[1]][24,]*totalatten, col='purple')
lines(Dale_NormSFR$Wave, Dale_NormSFR$Aspec[[11]][24,]*totalatten, col='orange')
lines(Dale_NormSFR$Wave, Dale_NormSFR$Aspec[[20]][24,]*totalatten, col='green')

dale_DM = dustmass(
  wave_star = stellar_wave, 
  lum_star_nodust = fluxBC03, 
  lum_star_dust = (birthBC03_atten$flux+screenBC03_atten$flux)/2, 
  wave_dust = Dale_NormSFR$Wave,
  Dale_NormSFR$Aspec[[11]][24,]
)

testSED=ProSpectSED(
  AGNlum=1e43, 
  speclib=BC03lr, 
  Dale=Dale_NormTot,
  AGN=Fritz, 
  Dale_M2L_func=Dale_M2L_func
)
plot(testSED)

testSED$dustmass


AGN = AGNinterp