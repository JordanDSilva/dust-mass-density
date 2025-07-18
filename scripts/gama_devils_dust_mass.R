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
lbtvec = cosdistTravelTime(z = zvec, ref = "Planck18")
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

mdustvec = seq(4.5, 15.5, 0.1)
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
              dnorm(p[1], x = 8.0, sd = 5, log = TRUE),
              dnorm(p[2], x = -1.0, sd = 8, log = TRUE),
              dnorm(p[3], x = -1.0, sd = 8, log = TRUE),
              dnorm(p[4], x = -2.0, sd = 2, log = TRUE),
              dnorm(p[5], x = -3.0, sd = 2, log = TRUE),
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
              dnorm(p[1], x = 8.0, sd = 5, log = TRUE),
              dnorm(p[2], x = -1.0, sd = 8, log = TRUE),
              dnorm(p[3], x = -1.0, sd = 8, log = TRUE),
              dnorm(p[4], x = -2.0, sd = 2, log = TRUE),
              dnorm(p[5], x = -3.0, sd = 2, log = TRUE),
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

gc()
smf_mc_err = foreach(i = 1:(length(zbins)-1), .errorhandling = "pass") %do% {
  
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
  
  temp = foreach(j = 1:100, .combine = rbind, .errorhandling = "remove") %do% {
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
    
    return( vmax_combine[,2] )
  }
  dev.off()

  return( colSds(as.matrix(temp), na.rm = TRUE) )
}
dmf_mc_err = foreach(i = 1:(length(zbins)-1), .errorhandling = "pass") %do% {
  
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
  
  temp = foreach(j = 1:100, .combine = rbind, .errorhandling = "remove") %do% {
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
    
    return( vmax_combine[,2] )
  }
  dev.off()
  
  return( colSds(as.matrix(temp), na.rm = TRUE) )
}
dmf_wAGN_mc_err = foreach(i = 1:(length(zbins)-1), .errorhandling = "pass") %do% {
  
  message(lbt_mids[i])
  png(paste0("~/Documents/DustMassDensity/plots/edd_bias_wAGN/lbt_", lbt_mids[i], ".png"), width = 7, height = 5, units = "in", res = 240)
  
  gama_x = sort_match_AGN$dustmass.total
  gama_err = 0.5 * (sort_match_AGN$dustmass.total_UB - sort_match_AGN$dustmass.total_LB)
  # gama_err = gama_err/(log(10) * gama_x)
  # gama_x = log10(gama_x)
  
  devils_x = devilsd10_AGN$dustmass.total
  devils_err = 0.5 * (devilsd10_AGN$dustmass.total_UB - devilsd10_AGN$dustmass.total_LB)
  # devils_err = devils_err/(log(10) * devils_x)
  # devils_x = log10(devils_x)
  
  temp = foreach(j = 1:100, .combine = rbind, .errorhandling = "remove") %do% {
    set.seed(j)
    gama_samples_ = rnorm(n = length(gama_x), mean = gama_x, sd = gama_err)
    
    GAMA = compute_mass_function(
      zlo = zbins[i],
      zhi = zbins[i+1],
      z = sort_match_AGN$z,
      x = gama_samples_,
      x_err = gama_err,
      areas = rep(217.54, length(sort_match_AGN$z)),
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
      z = devilsd10_AGN$z,
      x = devils_sample_,
      x_err = devils_err,
      areas = rep(1.5, length(devilsd10_AGN$z)),
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
    return( vmax_combine[,2] )
  }
  dev.off()
  
  return( colSds(as.matrix(temp), na.rm = TRUE) )
}

stellar_mass_density = foreach(i = 1:(length(zbins)-1)) %do% {
  
  png(paste0("~/Documents/DustMassDensity/plots/smf/lbt_",lbt_mids[i],".png"))
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
    vmaxErr = sqrt( vmax_combine[,3]^2 + smf_mc_err[[i]]^2 ),
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
  dev.off()
  return(fit)
}
dust_mass_density = foreach(i = 1:(length(zbins)-1)) %do% {
  
  png(paste0("~/Documents/DustMassDensity/plots/dmf/lbt_",lbt_mids[i],".png"))
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
    vmaxErr = sqrt( vmax_combine[,3]^2 + dmf_mc_err[[i]]^2),
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
  dev.off()
  return(fit)
}
dust_mass_density_wAGN = foreach(i = 1:(length(zbins)-1)) %do% {
  
  png(paste0("~/Documents/DustMassDensity/plots/dmf-AGN/lbt_",lbt_mids[i],".png"))
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
    vmaxErr = sqrt( vmax_combine[,3]^2 + dmf_wAGN_mc_err[[i]]^2 ),
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
  dev.off()
  return(fit)
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

h5write(obj = zmids, h5file, name = "zmids")
h5write(obj = zbins, h5file, name = "zbins")
h5write(obj = lbt_mids, h5file, name = "lbtmids")
h5write(obj = lbt_bins, h5file, name = "lbtbins")
h5write(obj = LSS_corr, h5file, name = "LSSCorrection")

h5delete(h5file, "cosmic")
h5createGroup(h5file, "cosmic")
h5write(obj = csmh, file = h5file, name = "cosmic/Mstar")
h5write(obj = cdmh, file = h5file, name = "cosmic/Mdust")
h5write(obj = cdmh, file = h5file, name = "cosmic/Mdust_gama")
h5write(obj = cdmh, file = h5file, name = "cosmic/Mdust_devils")
h5write(obj = cdmh_wAGN, file = h5file, name = "cosmic/MdustwAGN")

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

## Dust mass density
driver18_cdmh_raw = data.frame(fread("~/Documents/DustMassDensity/data/literature_evo/cdmh/driver18_raw.csv", skip = 3))
driver18_cdmh = data.frame(
    "lbt" = driver18_cdmh_raw$V1,
    "z" = lbt2z(v = unlist(driver18_cdmh_raw$V1)),
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
