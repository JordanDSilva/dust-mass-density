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

fir_input = devilsd10_AGN$FIRInput
fir_input[is.na(fir_input)] = 0

sum( fir_input == 1 & devilsd10_AGN$z > 0) /  sum(devilsd10_AGN$z > 0)

RR14_BPL = function(Z, par, doDTG = FALSE){
  
  ## Remy Ruyer+14 using metallicity dependent XCO
  ## x = 12 + log(O/H)
  ## xSol = 12 + log(O/H)Sol
  
  ## Z/Zsol = (O/H)/(O/HSol)
  ## log(Z) - log(Zsol) = log(O/H) - log(O/HSol)
  ## log(Z)+12 - log(Zsol)-12 = log(O/H)+12 - log(O/Hsol)-12
  ## log(Z) - log(Zsol) = x - xSol
  ## log(Z/0.014) + xSol = log(O/H) + 12
  
  a = 2.21
  alphaH = 1.00
  b = 0.96
  # alphaL = 3.10
  # xt = 8.10
  
  # a = par[1]
  # alphaH = par[2]
  # b = par[3]
  alphaL = par[1]
  xt = par[2]
  
  xSol = 8.69
  ZOH = log10(Z / 0.014) + xSol
  
  GTD = ifelse(
    ZOH > xt,
    a + alphaH*(xSol - ZOH),
    b + alphaL*(xSol - ZOH)
  )
  DTG = (10^GTD)^-1
  
  ##Mgas = muGal * Mhydrogen
  ## DTG = Mdust / Mgas
  ## DTH = Mdust / Mhydrogen = Mdust / (Mgas / muGal) = DTG / muGal = DTG / (1 / (1 - Ysol - Zgal))
  muGal = 1 / (1 - 0.270 - Z)
  # muGal = (1-Z) / (1 + (4/3))
  DTH = DTG/muGal
  if(doDTG){
    return(DTG)
  }else{
    return(DTH)
  }
}

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

mdustvec = seq(3.5, 15.5, 0.1)
compute_mass_function = function(zlo, zhi, z, x, x_err, areas, sm_bins, errFloor, fit_fun, vmax_bins = NULL, vmax = NULL, vmaxErr = NULL, meanxErr = NULL, do_optim = FALSE, do_fit = FALSE, do_fit_quantiles = TRUE, Niters = 2000, do_plot = FALSE, add = TRUE, pt.col = "black", ln.alpha = 1.0){
  
  sm_mids = sm_bins[-length(sm_bins)] + diff(sm_bins) / 2
  ddmm = abs(sm_bins[1] - sm_bins[2])
  vol = 4*pi/3 * (cosdistCoDist(z = zhi, ref = 'Planck18')^3 - cosdistCoDist(z = zlo, ref = 'Planck18')^3)
  
  zidx = z >= zlo & z < zhi & x > 0
  Mdust = log10(x[zidx])
  MdustErr = x_err[zidx] / (log(10) * x[zidx])
  vol = rep(vol, sum(zidx)) * areas[zidx] / (4*pi*(180/pi)^2)

  if(is.null(vmax) & is.null(vmaxErr) & is.null(vmax_bins) & is.null(meanxErr)){
    # vvmax = foreach(j = 1:length(sm_mids), .combine = "c") %do% {
    #   midx = Mdust >= sm_bins[j] & Mdust < sm_bins[j+1]
    #   sum(vol[midx]^-1)
    # }/ddmm
    # vvmaxErr = foreach(j = 1:length(sm_mids), .combine = "c") %do% {
    #   midx = Mdust >= sm_bins[j] & Mdust < sm_bins[j+1]
    #   sqrt(sum(vol[midx]^-2))
    # }/ddmm
    
    hh = maghist(Mdust, breaks = sm_bins, plot = FALSE)
    vvmax = hh$counts / (unique(vol) * ddmm)
    vvmaxErr = sqrt(hh$counts) / (unique(vol) * ddmm)
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
  Mdust_lim = sm_mids[which.max(vvmax)]
  Mdust_uplim = rev(sm_mids)[max(which( sapply(rev(sm_mids), function(x) sum(vvmax[sm_mids <= x & sm_mids > Mdust_lim]==0) >= 2 ) ))]
  
  if(do_fit){

    if(do_optim){
      opt = optim(
        par = c(8, -1, -1, -2, -3), 
        fn = LL, 
        # method = "L-BFGS-B",
        control = list("fnscale" = -1), 
        Data = list(
          xx = sm_mids[vvmax > 0 & sm_mids >= Mdust_lim], 
          yy = vvmax[vvmax > 0 & sm_mids >= Mdust_lim], 
          yyerr = vvmaxErr[vvmax > 0 & sm_mids >= Mdust_lim], 
          func = double_schechter, 
          prior = function(p){
            sum(
              dnorm(p[1], x = 8.0, sd = 4.0, log = TRUE),
              dnorm(p[2], x = -1.0, sd = 1.0, log = TRUE),
              dnorm(p[3], x = -1.0, sd = 1.0, log = TRUE),
              dnorm(p[4], x = -2.0, sd = 4.0, log = TRUE),
              dnorm(p[5], x = -3.0, sd = 4.0, log = TRUE),
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
      pinit = c(10, 0.0, 0.0, -4, -4)
      highout = Highlander(
        parm = pinit, 
        lower = c(5, -2.5, -2.5, -8, -8),
        upper = c(15, 2.5, 2.5, 1, 1),
        Data = list(
          xx = sm_mids[vvmax > 0 & sm_mids >= Mdust_lim], 
          yy = vvmax[vvmax > 0 & sm_mids >= Mdust_lim], 
          yyerr = vvmaxErr[vvmax > 0 & sm_mids >= Mdust_lim], 
          func = double_schechter, 
          prior = function(p){
            sum(
              dnorm(p[1], x = 8.0, sd = 4.0, log = TRUE),
              dnorm(p[2], x = -1.0, sd = 1.0, log = TRUE),
              dnorm(p[3], x = -1.0, sd = 1.0, log = TRUE),
              dnorm(p[4], x = -2.0, sd = 4.0, log = TRUE),
              dnorm(p[5], x = -3.0, sd = 4.0, log = TRUE),
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
      10^sm_mids, vvmax, ylo = vvmaxErr, 
      col = alpha(pt.col, ln.alpha)
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
## Fold in also RR relation error and metallicity error?
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
  
  temp = foreach(j = 1:500, .combine = rbind, .errorhandling = "remove") %do% {
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
  gama_Metallicity = sort_match_noAGN$Zgas_50
  gama_Metallicity_err = 0.5 * (sort_match_noAGN$Zgas_84 - sort_match_noAGN$Zgas_16)
  # gama_err = gama_err/(log(10) * gama_x)
  # gama_x = log10(gama_x)
  
  devils_x = devilsd10_noAGN$dustmass.total
  devils_err = 0.5 * (devilsd10_noAGN$dustmass.total_UB - devilsd10_noAGN$dustmass.total_LB)
  devils_Metallicity = devilsd10_noAGN$Zfinal
  devils_Metallicity_err = 0.5 * (devilsd10_noAGN$Zfinal_UB - devilsd10_noAGN$Zfinal_LB)
  # devils_err = devils_err/(log(10) * devils_x)
  # devils_x = log10(devils_x)
  
  # RR_samples = cbind(
  #   rnorm(n = 500, mean = 3.10, sd = 1.33),
  #   rnorm(n = 500, mean = 8.10, sd = 0.43)
  # )
  
  temp = foreach(j = 1:500, .combine = rbind, .errorhandling = "remove") %do% {
    set.seed(j)
    
    RR_temp = function(Z){RR14_BPL(Z = Z, par = c(3.10, 8.10))}
    
    gama_samples_ = rnorm(n = length(gama_x), mean = gama_x, sd = gama_err)
    # gama_metals_samples_ = rnorm(n = length(gama_Metallicity), mean = gama_Metallicity, sd = gama_Metallicity_err)
    
    GAMA = compute_mass_function(
      zlo = zbins[i],
      zhi = zbins[i+1],
      z = sort_match_noAGN$z,
      x = gama_samples_ * 1 / (0.0073/RR_temp(10^gama_Metallicity)),
      x_err = gama_err * 1 / (0.0073/RR_temp(10^gama_Metallicity)),
      areas = rep(217.54, length(sort_match_noAGN$z)),
      sm_bins = sm_bins,
      errFloor = 0.0,
      do_fit = FALSE,
      do_plot = TRUE,
      add = j!=1,
      pt.col = "purple",
      ln.alpha = 0.1
    )
    
    devils_samples_ = rnorm(n = length(devils_x), mean = devils_x, sd = devils_err)
    # devils_metals_samples_ = rnorm(n = length(devils_Metallicity), mean = devils_Metallicity, sd = devils_Metallicity_err)
    DEVILS = compute_mass_function(
      zlo = zbins[i],
      zhi = zbins[i+1],
      z = devilsd10_noAGN$z,
      x = devils_samples_ * 1 / (0.0073/RR_temp(10^devils_Metallicity)),
      x_err = devils_err * 1 / (0.0073/RR_temp(10^devils_Metallicity)),
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
  gama_Metallicity = sort_match_AGN$Zfinal
  gama_Metallicity_err = 0.5 * (sort_match_AGN$Zfinal_UB - sort_match_AGN$Zfinal_LB)
  # gama_err = gama_err/(log(10) * gama_x)
  # gama_x = log10(gama_x)
  
  devils_x = devilsd10_AGN$dustmass.total
  devils_err = 0.5 * (devilsd10_AGN$dustmass.total_UB - devilsd10_AGN$dustmass.total_LB)
  devils_Metallicity = devilsd10_AGN$Zfinal
  devils_Metallicity_err = 0.5 * (devilsd10_AGN$Zfinal_UB - devilsd10_AGN$Zfinal_LB)
  # devils_err = devils_err/(log(10) * devils_x)
  # devils_x = log10(devils_x)
  
  RR_samples = cbind(
    rnorm(n = 500, mean = 3.10, sd = 1.33),
    rnorm(n = 500, mean = 8.10, sd = 0.43)
  )
  
  temp = foreach(j = 1:500, .combine = rbind, .errorhandling = "remove") %do% {
    set.seed(j)
    
    RR_temp = function(Z){RR14_BPL(Z = Z, par = c(3.10, 8.10))}
    
    gama_samples_ = rnorm(n = length(gama_x), mean = gama_x, sd = gama_err)
    # gama_metals_samples_ = rnorm(n = length(gama_Metallicity), mean = gama_Metallicity, sd = gama_Metallicity_err)
    
    GAMA = compute_mass_function(
      zlo = zbins[i],
      zhi = zbins[i+1],
      z = sort_match_AGN$z,
      x = gama_samples_ * 1 / (0.0073/RR_temp(10^gama_Metallicity)),
      x_err = gama_err * 1 / (0.0073/RR_temp(10^gama_Metallicity)),
      areas = rep(217.54, length(sort_match_AGN$z)),
      sm_bins = sm_bins,
      errFloor = 0.0,
      do_fit = FALSE,
      do_plot = TRUE,
      add = j!=1,
      pt.col = "purple",
      ln.alpha = 0.1
    )
    
    devils_samples_ = rnorm(n = length(devils_x), mean = devils_x, sd = devils_err)
    # devils_metals_samples_ = rnorm(n = length(devils_Metallicity), mean = devils_Metallicity, sd = devils_Metallicity_err)
    
    DEVILS = compute_mass_function(
      zlo = zbins[i],
      zhi = zbins[i+1],
      z = devilsd10_AGN$z,
      x = devils_samples_ * 1 / (0.0073/RR_temp(10^devils_Metallicity)),
      x_err = devils_err * 1 / (0.0073/RR_temp(10^devils_Metallicity)),
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
dmf_mc_no_corr_err = foreach(i = 1:(length(zbins)-1), .errorhandling = "pass") %do% {
  
  message(lbt_mids[i])
  # png(paste0("~/Documents/DustMassDensity/plots/edd_bias/lbt_", lbt_mids[i], ".png"), width = 7, height = 5, units = "in", res = 240)
  
  gama_x = sort_match_noAGN$DustMass_50
  gama_err = 0.5 * (sort_match_noAGN$DustMass_84 - sort_match_noAGN$DustMass_16)

  devils_x = devilsd10_noAGN$dustmass.total
  devils_err = 0.5 * (devilsd10_noAGN$dustmass.total_UB - devilsd10_noAGN$dustmass.total_LB)

  temp = foreach(j = 1:500, .combine = rbind, .errorhandling = "remove") %do% {
    set.seed(j)
    
    RR_temp = function(Z){RR14_BPL(Z = Z, par = c(3.10, 8.10))}
    
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
    
    devils_samples_ = rnorm(n = length(devils_x), mean = devils_x, sd = devils_err)
    DEVILS = compute_mass_function(
      zlo = zbins[i],
      zhi = zbins[i+1],
      z = devilsd10_noAGN$z,
      x = devils_samples_,
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
  # dev.off()
  
  return( colSds(as.matrix(temp), na.rm = TRUE) )
}
dmf_wAGN_mc_no_corr_err = foreach(i = 1:(length(zbins)-1), .errorhandling = "pass") %do% {
  
  message(lbt_mids[i])
  # png(paste0("~/Documents/DustMassDensity/plots/edd_bias_wAGN/lbt_", lbt_mids[i], ".png"), width = 7, height = 5, units = "in", res = 240)
  
  gama_x = sort_match_AGN$dustmass.total
  gama_err = 0.5 * (sort_match_AGN$dustmass.total_UB - sort_match_AGN$dustmass.total_LB)

  devils_x = devilsd10_AGN$dustmass.total
  devils_err = 0.5 * (devilsd10_AGN$dustmass.total_UB - devilsd10_AGN$dustmass.total_LB)

  temp = foreach(j = 1:500, .combine = rbind, .errorhandling = "remove") %do% {
    set.seed(j)
    
    RR_temp = function(Z){RR14_BPL(Z = Z, par = c(3.10, 8.10))}
    
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
    
    devils_samples_ = rnorm(n = length(devils_x), mean = devils_x, sd = devils_err)
    DEVILS = compute_mass_function(
      zlo = zbins[i],
      zhi = zbins[i+1],
      z = devilsd10_AGN$z,
      x = devils_samples_,
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
  # dev.off()
  return( colSds(as.matrix(temp), na.rm = TRUE) )
}
gmf_mc_err = foreach(i = 1:(length(zbins)-1), .errorhandling = "pass") %do% {
  
  message(lbt_mids[i])
  png(paste0("~/Documents/DustMassDensity/plots/gedd_bias/lbt_", lbt_mids[i], ".png"), width = 7, height = 5, units = "in", res = 240)
  
  gama_x = sort_match_noAGN$DustMass_50
  gama_err = 0.5 * (sort_match_noAGN$DustMass_84 - sort_match_noAGN$DustMass_16)
  gama_Metallicity = sort_match_noAGN$Zgas_50
  gama_Metallicity_err = 0.5 * (sort_match_noAGN$Zgas_84 - sort_match_noAGN$Zgas_16)
  # gama_err = gama_err/(log(10) * gama_x)
  # gama_x = log10(gama_x)
  
  devils_x = devilsd10_noAGN$dustmass.total
  devils_err = 0.5 * (devilsd10_noAGN$dustmass.total_UB - devilsd10_noAGN$dustmass.total_LB)
  devils_Metallicity = devilsd10_noAGN$Zfinal
  devils_Metallicity_err = 0.5 * (devilsd10_noAGN$Zfinal_UB - devilsd10_noAGN$Zfinal_LB)
  # devils_err = devils_err/(log(10) * devils_x)
  # devils_x = log10(devils_x)
  
  ## Way too much error so dont use this
  # RR_samples = cbind(
  #   rnorm(n = 500, mean = 3.10, sd = 1.33),
  #   rnorm(n = 500, mean = 8.10, sd = 0.43)
  # )
  RR_temp = function(Z){RR14_BPL(Z = Z, par = c(3.10, 8.10), doDTG = TRUE)}
  
  temp = foreach(j = 1:500, .combine = rbind, .errorhandling = "remove") %do% {
    set.seed(j)
    
    gama_samples_ = rnorm(n = length(gama_x), mean = gama_x, sd = gama_err)
    # gama_metals_samples_ = rnorm(n = length(gama_Metallicity), mean = gama_Metallicity, sd = gama_Metallicity_err)
    GAMA = compute_mass_function(
      zlo = zbins[i],
      zhi = zbins[i+1],
      z = sort_match_noAGN$z,
      x = gama_samples_ * 1 / (0.0073/RR_temp(10^gama_Metallicity))/RR_temp(10^gama_Metallicity),
      x_err = gama_err * 1 / (0.0073/RR_temp(10^gama_Metallicity))/RR_temp(10^gama_Metallicity),
      areas = rep(217.54, length(sort_match_noAGN$z)),
      sm_bins = sm_bins,
      errFloor = 0.0,
      do_fit = FALSE,
      do_plot = TRUE,
      add = j!=1,
      pt.col = "purple",
      ln.alpha = 0.1
    )
    
    devils_samples_ = rnorm(n = length(devils_x), mean = devils_x, sd = devils_err)
    # devils_metals_samples_ = rnorm(n = length(devils_Metallicity), mean = devils_Metallicity, sd = devils_Metallicity_err)
    DEVILS = compute_mass_function(
      zlo = zbins[i],
      zhi = zbins[i+1],
      z = devilsd10_noAGN$z,
      x = devils_samples_ * 1 / (0.0073/RR_temp(10^devils_Metallicity))/RR_temp(10^devils_Metallicity),
      x_err = devils_err * 1 / (0.0073/RR_temp(10^devils_Metallicity))/RR_temp(10^devils_Metallicity),
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
gmf_wAGN_mc_err = foreach(i = 1:(length(zbins)-1), .errorhandling = "pass") %do% {
  
  message(lbt_mids[i])
  png(paste0("~/Documents/DustMassDensity/plots/gedd_bias_wAGN/lbt_", lbt_mids[i], ".png"), width = 7, height = 5, units = "in", res = 240)
  
  gama_x = sort_match_AGN$dustmass.total
  gama_err = 0.5 * (sort_match_AGN$dustmass.total_UB - sort_match_AGN$dustmass.total_LB)
  gama_Metallicity = sort_match_AGN$Zfinal
  gama_Metallicity_err = 0.5 * (sort_match_AGN$Zfinal_UB - sort_match_AGN$Zfinal)
  # gama_err = gama_err/(log(10) * gama_x)
  # gama_x = log10(gama_x)
  
  devils_x = devilsd10_AGN$dustmass.total
  devils_err = 0.5 * (devilsd10_AGN$dustmass.total_UB - devilsd10_AGN$dustmass.total_LB)
  devils_Metallicity = devilsd10_AGN$Zfinal
  devils_Metallicity_err = 0.5 * (devilsd10_AGN$Zfinal_UB - devilsd10_AGN$Zfinal)
  # devils_err = devils_err/(log(10) * devils_x)
  # devils_x = log10(devils_x)
  
  ## Way too much error so dont use this
  # RR_samples = cbind(
  #   rnorm(n = 500, mean = 3.10, sd = 1.33),
  #   rnorm(n = 500, mean = 8.10, sd = 0.43)
  # )
  RR_temp = function(Z){RR14_BPL(Z = Z, par = c(3.10, 8.10), doDTG = TRUE)}
  
  temp = foreach(j = 1:500, .combine = rbind, .errorhandling = "remove") %do% {
    set.seed(j)
    
    
    gama_samples_ = rnorm(n = length(gama_x), mean = gama_x, sd = gama_err)
    # gama_metals_samples_ = rnorm(n = length(gama_Metallicity), mean = gama_Metallicity, sd = gama_Metallicity_err)
    GAMA = compute_mass_function(
      zlo = zbins[i],
      zhi = zbins[i+1],
      z = sort_match_AGN$z,
      x = gama_samples_ * 1 / (0.0073/RR_temp(10^gama_Metallicity))/RR_temp(10^gama_Metallicity),
      x_err = gama_err * 1 / (0.0073/RR_temp(10^gama_Metallicity))/RR_temp(10^gama_Metallicity),
      areas = rep(217.54, length(sort_match_AGN$z)),
      sm_bins = sm_bins,
      errFloor = 0.0,
      do_fit = FALSE,
      do_plot = TRUE,
      add = j!=1,
      pt.col = "purple",
      ln.alpha = 0.1
    )
    
    devils_samples_ = rnorm(n = length(devils_x), mean = devils_x, sd = devils_err)
    # devils_metals_samples_ = rnorm(n = length(devils_Metallicity), mean = devils_Metallicity, sd = devils_Metallicity_err)
    DEVILS = compute_mass_function(
      zlo = zbins[i],
      zhi = zbins[i+1],
      z = devilsd10_AGN$z,
      x = devils_samples_ * 1 / (0.0073/RR_temp(10^devils_Metallicity))/RR_temp(10^devils_Metallicity),
      x_err = devils_err * 1 / (0.0073/RR_temp(10^devils_Metallicity))/RR_temp(10^devils_Metallicity),
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
    errFloor = 0.0,
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
    errFloor = 0.0,
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
    errFloor = 0.1,
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
    x = sort_match_noAGN$DustMass_50 * 1 / (0.0073/RR14_BPL(10^sort_match_noAGN$Zgas_50, par = c(3.10, 8.10))),
    x_err = 0.5 * (sort_match_noAGN$DustMass_84 - sort_match_noAGN$DustMass_16) * 1 / (0.0073/RR14_BPL(10^sort_match_noAGN$Zgas_50, par = c(3.10, 8.10))),
    areas = rep(217.54, length(sort_match_noAGN$z)),
    sm_bins = sm_bins,
    errFloor = 0.0,
    do_fit = FALSE,
    do_plot = TRUE,
    add = FALSE,
    pt.col = "purple"
  )
  DEVILS = compute_mass_function(
    zlo = zbins[i],
    zhi = zbins[i+1],
    z = devilsd10_noAGN$z,
    x = devilsd10_noAGN$dustmass.total * 1 / (0.0073/RR14_BPL(10^devilsd10_noAGN$Zfinal, par = c(3.10, 8.10))),
    x_err = 0.5 * (devilsd10_noAGN$dustmass.total_UB - devilsd10_noAGN$dustmass.total_LB) * 1 / (0.0073/RR14_BPL(10^devilsd10_noAGN$Zfinal, par = c(3.10, 8.10))),
    areas = rep(1.5, length(devilsd10_noAGN$z)),
    sm_bins = sm_bins,
    errFloor = 0.0,
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
    errFloor = 0.1,
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
dust_mass_density_no_corr = foreach(i = 1:(length(zbins)-1)) %do% {
  
  # png(paste0("~/Documents/DustMassDensity/plots/dmf/lbt_",lbt_mids[i],".png"))
  GAMA = compute_mass_function(
    zlo = zbins[i],
    zhi = zbins[i+1],
    z = sort_match_noAGN$z,
    x = sort_match_noAGN$DustMass_50,
    x_err = 0.5 * (sort_match_noAGN$DustMass_84 - sort_match_noAGN$DustMass_16),
    areas = rep(217.54, length(sort_match_noAGN$z)),
    sm_bins = sm_bins,
    errFloor = 0.0,
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
    errFloor = 0.0,
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
    vmaxErr = sqrt( vmax_combine[,3]^2 + dmf_mc_no_corr_err[[i]]^2),
    meanxErr = vmax_combine[,4],
    errFloor = 0.1,
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
    x = sort_match_AGN$dustmass.total * 1 / (0.0073/RR14_BPL(10^sort_match_AGN$Zfinal, par = c(3.10, 8.10))),
    x_err = 0.5 * (sort_match_AGN$dustmass.total_UB - sort_match_AGN$dustmass.total_LB) * 1 / (0.0073/RR14_BPL(10^sort_match_AGN$Zfinal, par = c(3.10, 8.10))),
    areas = rep(217.54, length(sort_match_AGN$z)),
    sm_bins = sm_bins,
    errFloor = 0.0,
    do_fit = FALSE,
    do_plot = TRUE,
    add = FALSE,
    pt.col = "purple"
  )
  DEVILS = compute_mass_function(
    zlo = zbins[i],
    zhi = zbins[i+1],
    z = devilsd10_AGN$z,
    x = devilsd10_AGN$dustmass.total * 1 / (0.0073/RR14_BPL(10^devilsd10_AGN$Zfinal, par = c(3.10, 8.10))),
    x_err = 0.5 * (devilsd10_AGN$dustmass.total_UB - devilsd10_AGN$dustmass.birth_LB) * 1 / (0.0073/RR14_BPL(10^devilsd10_AGN$Zfinal, par = c(3.10, 8.10))),
    areas = rep(1.5, length(devilsd10_AGN$z)),
    sm_bins = sm_bins,
    errFloor = 0.0,
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
    errFloor = 0.1,
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
dust_mass_density_wAGN_no_corr = foreach(i = 1:(length(zbins)-1)) %do% {
  
  # png(paste0("~/Documents/DustMassDensity/plots/dmf-AGN/lbt_",lbt_mids[i],".png"))
  GAMA = compute_mass_function(
    zlo = zbins[i],
    zhi = zbins[i+1],
    z = sort_match_AGN$z,
    x = sort_match_AGN$dustmass.total,
    x_err = 0.5 * (sort_match_AGN$dustmass.total_UB - sort_match_AGN$dustmass.total_LB),
    areas = rep(217.54, length(sort_match_AGN$z)),
    sm_bins = sm_bins,
    errFloor = 0.0,
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
    errFloor = 0.0,
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
    vmaxErr = sqrt( vmax_combine[,3]^2 + dmf_wAGN_mc_no_corr_err[[i]]^2 ),
    meanxErr = vmax_combine[,4],
    errFloor = 0.1,
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

gas_mass_density = foreach(i = 1:(length(zbins)-1)) %do% {
  
  png(paste0("~/Documents/DustMassDensity/plots/gmf/lbt_",lbt_mids[i],".png"))
  GAMA = compute_mass_function(
    zlo = zbins[i],
    zhi = zbins[i+1],
    z = sort_match_noAGN$z,
    x = sort_match_noAGN$DustMass_50 * 1/(0.0073/RR14_BPL(10^sort_match_noAGN$Zgas_50, par = c(3.10, 8.10), doDTG = TRUE))/RR14_BPL(10^sort_match_noAGN$Zgas_50, par = c(3.10, 8.10), doDTG = TRUE),
    x_err = 0.5 * (sort_match_noAGN$DustMass_84 - sort_match_noAGN$DustMass_16) * 1/(0.0073/RR14_BPL(10^sort_match_noAGN$Zgas_50, par = c(3.10, 8.10), doDTG = TRUE))/RR14_BPL(10^sort_match_noAGN$Zgas_50, par = c(3.10, 8.10), doDTG = TRUE),
    areas = rep(217.54, length(sort_match_noAGN$z)),
    sm_bins = sm_bins,
    errFloor = 0.0,
    do_fit = FALSE,
    do_plot = TRUE,
    add = FALSE,
    pt.col = "purple"
  )
  DEVILS = compute_mass_function(
    zlo = zbins[i],
    zhi = zbins[i+1],
    z = devilsd10_noAGN$z,
    x = devilsd10_noAGN$dustmass.total * 1/(0.0073/RR14_BPL(10^devilsd10_noAGN$Zfinal, par = c(3.10, 8.10), doDTG = TRUE))/RR14_BPL(10^devilsd10_noAGN$Zfinal, par = c(3.10, 8.10), doDTG = TRUE),
    x_err = 0.5 * (devilsd10_noAGN$dustmass.total_UB - devilsd10_noAGN$dustmass.total_LB) * 1/(0.0073/RR14_BPL(10^devilsd10_noAGN$Zfinal, par = c(3.10, 8.10), doDTG = TRUE))/RR14_BPL(10^devilsd10_noAGN$Zfinal, par = c(3.10, 8.10), doDTG = TRUE),
    areas = rep(1.5, length(devilsd10_noAGN$z)),
    sm_bins = sm_bins,
    errFloor = 0.0,
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
    vmaxErr = sqrt( vmax_combine[,3]^2 + gmf_mc_err[[i]]^2),
    meanxErr = vmax_combine[,4],
    errFloor = 0.1,
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
gas_mass_density_wAGN = foreach(i = 1:(length(zbins)-1)) %do% {
  
  png(paste0("~/Documents/DustMassDensity/plots/gmf-AGN/lbt_",lbt_mids[i],".png"))
  GAMA = compute_mass_function(
    zlo = zbins[i],
    zhi = zbins[i+1],
    z = sort_match_AGN$z,
    x = sort_match_AGN$dustmass.total * 1 / (0.0073/RR14_BPL(10^sort_match_AGN$Zfinal, par = c(3.10, 8.10), doDTG = TRUE))/RR14_BPL(10^sort_match_AGN$Zfinal, par = c(3.10, 8.10), doDTG = TRUE),
    x_err = 0.5 * (sort_match_AGN$dustmass.total_UB - sort_match_AGN$dustmass.total_LB) * 1 / (0.0073/RR14_BPL(10^sort_match_AGN$Zfinal, par = c(3.10, 8.10), doDTG = TRUE))/RR14_BPL(10^sort_match_AGN$Zfinal, par = c(3.10, 8.10), doDTG = TRUE),
    areas = rep(217.54, length(sort_match_AGN$z)),
    sm_bins = sm_bins,
    errFloor = 0.0,
    do_fit = FALSE,
    do_plot = TRUE,
    add = FALSE,
    pt.col = "purple"
  )
  DEVILS = compute_mass_function(
    zlo = zbins[i],
    zhi = zbins[i+1],
    z = devilsd10_AGN$z,
    x = devilsd10_AGN$dustmass.total * 1 / (0.0073/RR14_BPL(10^devilsd10_AGN$Zfinal, par = c(3.10, 8.10), doDTG = TRUE))/RR14_BPL(10^devilsd10_AGN$Zfinal, par = c(3.10, 8.10), doDTG = TRUE),
    x_err = 0.5 * (devilsd10_AGN$dustmass.total_UB - devilsd10_AGN$dustmass.birth_LB) * 1 / (0.0073/RR14_BPL(10^devilsd10_AGN$Zfinal, par = c(3.10, 8.10), doDTG = TRUE))/RR14_BPL(10^devilsd10_AGN$Zfinal, par = c(3.10, 8.10), doDTG = TRUE),
    areas = rep(1.5, length(devilsd10_AGN$z)),
    sm_bins = sm_bins,
    errFloor = 0.0,
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
    vmaxErr = sqrt( vmax_combine[,3]^2 + gmf_wAGN_mc_err[[i]]^2 ),
    meanxErr = vmax_combine[,4],
    errFloor = 0.1,
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

for(i in 1:10000){
  tryCatch(
    dev.off(),
    error = function(e){NULL}
  )
}

gc()
csmh = data.frame(foreach(i = 1:length(stellar_mass_density), .combine = bind_rows) %do% {
  stellar_mass_density[[i]]$cosmic
})
cdmh = data.frame(foreach(i = 1:length(dust_mass_density), .combine = bind_rows) %do% {
  dust_mass_density[[i]]$cosmic
})
cdmh_no_corr = data.frame(foreach(i = 1:length(dust_mass_density_no_corr), .combine = bind_rows) %do% {
  dust_mass_density_no_corr[[i]]$cosmic
})
cdmh_wAGN = data.frame(foreach(i = 1:length(dust_mass_density_wAGN), .combine = bind_rows) %do% {
  dust_mass_density_wAGN[[i]]$cosmic
})
cdmh_wAGN_no_corr = data.frame(foreach(i = 1:length(dust_mass_density_wAGN_no_corr), .combine = bind_rows) %do% {
  dust_mass_density_wAGN_no_corr[[i]]$cosmic
})
# cdmh_wAGN_no_corr_noERR = data.frame(foreach(i = 1:length(dust_mass_density_wAGN_no_corr_noERR), .combine = bind_rows) %do% {
#   dust_mass_density_wAGN_no_corr_noERR[[i]]$cosmic
# })
# cdmh_birth_wAGN_no_corr_noERR = data.frame(foreach(i = 1:length(dust_birth_mass_density_wAGN_no_corr_noERR), .combine = bind_rows) %do% {
#   dust_birth_mass_density_wAGN_no_corr_noERR[[i]]$cosmic
# })
# cdmh_screen_wAGN_no_corr_noERR = data.frame(foreach(i = 1:length(dust_screen_mass_density_wAGN_no_corr_noERR), .combine = bind_rows) %do% {
#   dust_screen_mass_density_wAGN_no_corr_noERR[[i]]$cosmic
# })
cgmh = data.frame(foreach(i = 1:length(gas_mass_density), .combine = bind_rows) %do% {
  gas_mass_density[[i]]$cosmic
})
cgmh_wAGN = data.frame(foreach(i = 1:length(gas_mass_density_wAGN), .combine = bind_rows) %do% {
  gas_mass_density_wAGN[[i]]$cosmic
})

# magplot(lbt_mids, (cdmh_wAGN_no_corr_noERR$Q50), ylim = 10**c(3, 7), log = "y")
# points(lbt_mids, (cdmh_birth_wAGN_no_corr_noERR$Q50), col = "red")
# magerr(lbt_mids, (cdmh_birth_wAGN_no_corr_noERR$Q50), ylo=cdmh_birth_wAGN_no_corr_noERR$ERR, col = "red")
# points(lbt_mids, (cdmh_screen_wAGN_no_corr_noERR$Q50), col = "blue")
# magerr(lbt_mids, (cdmh_screen_wAGN_no_corr_noERR$Q50), ylo=cdmh_screen_wAGN_no_corr_noERR$ERR, col = "blue")

names_par = c("M", "alpha", "beta", "phi1", "phi2")
names_par_err = paste0(names_par, "Err")
dmf_par = data.frame(foreach(i = 1:length(dust_mass_density), .combine = rbind) %do% {
  temp = colQuantiles(dust_mass_density[[i]]$highout$LD_last$Posterior1, probs = c(0.5, 0.16, 0.84))
  cbind(
    t(temp[,1])
  )
})
dmf_par_err = data.frame(foreach(i = 1:length(dust_mass_density), .combine = rbind) %do% {
  temp = colQuantiles(dust_mass_density[[i]]$highout$LD_last$Posterior1, probs = c(0.5, 0.16, 0.84))
  cbind(
    t(0.5 * (temp[,3] - temp[,2]))
  )
})
dmf_wAGN_par = data.frame(foreach(i = 1:length(dust_mass_density_wAGN), .combine = rbind) %do% {
  temp = colQuantiles(dust_mass_density_wAGN[[i]]$highout$LD_last$Posterior1, probs = c(0.5, 0.16, 0.84))
  cbind(
    t(temp[,1])
  )
})
dmf_wAGN_par_err = data.frame(foreach(i = 1:length(dust_mass_density_wAGN), .combine = rbind) %do% {
  temp = colQuantiles(dust_mass_density_wAGN[[i]]$highout$LD_last$Posterior1, probs = c(0.5, 0.16, 0.84))
  cbind(
    t(0.5 * (temp[,3] - temp[,2]))
  )
})
gmf_par = data.frame(foreach(i = 1:length(gas_mass_density), .combine = rbind) %do% {
  temp = colQuantiles(gas_mass_density[[i]]$highout$LD_last$Posterior1, probs = c(0.5, 0.16, 0.84))
  cbind(
    t(temp[,1])
  )
})
gmf_par_err = data.frame(foreach(i = 1:length(gas_mass_density), .combine = rbind) %do% {
  temp = colQuantiles(gas_mass_density[[i]]$highout$LD_last$Posterior1, probs = c(0.5, 0.16, 0.84))
  cbind(
    t(0.5 * (temp[,3] - temp[,2]))
  )
})
gmf_wAGN_par = data.frame(foreach(i = 1:length(gas_mass_density_wAGN), .combine = rbind) %do% {
  temp = colQuantiles(gas_mass_density_wAGN[[i]]$highout$LD_last$Posterior1, probs = c(0.5, 0.16, 0.84))
  cbind(
    t(temp[,1])
  )
})
gmf_wAGN_par_err = data.frame(foreach(i = 1:length(gas_mass_density_wAGN), .combine = rbind) %do% {
  temp = colQuantiles(gas_mass_density_wAGN[[i]]$highout$LD_last$Posterior1, probs = c(0.5, 0.16, 0.84))
  cbind(
    t(0.5 * (temp[,3] - temp[,2]))
  )
})
smf_par = data.frame(foreach(i = 1:length(stellar_mass_density), .combine = rbind) %do% {
  temp = colQuantiles(stellar_mass_density[[i]]$highout$LD_last$Posterior1, probs = c(0.5, 0.16, 0.84))
  cbind(
    t(temp[,1])
  )
})
smf_par_err = data.frame(foreach(i = 1:length(stellar_mass_density), .combine = rbind) %do% {
  temp = colQuantiles(stellar_mass_density[[i]]$highout$LD_last$Posterior1, probs = c(0.5, 0.16, 0.84))
  cbind(
    t(0.5 * (temp[,3] - temp[,2]))
  )
})
names(dmf_par) = names_par
names(dmf_wAGN_par) = names_par
names(smf_par) = names_par
names(gmf_par) = names_par
names(gmf_wAGN_par) = names_par
names(dmf_par_err) = names_par_err
names(dmf_wAGN_par_err) = names_par_err
names(smf_par_err) = names_par_err
names(gmf_par_err) = names_par_err
names(gmf_wAGN_par_err) = names_par_err

dmf_par = data.frame(cbind(dmf_par, dmf_par_err))
dmf_wAGN_par = data.frame(cbind(dmf_wAGN_par, dmf_wAGN_par_err))
gmf_par = data.frame(cbind(gmf_par, gmf_par_err))
gmf_wAGN_par = data.frame(cbind(gmf_wAGN_par, gmf_wAGN_par_err))
smf_par = data.frame(cbind(smf_par, smf_par_err))

dsilva25 = data.frame(fread("~/Documents/DustMassDensity/data/literature_evo/csfh/DSilva25_CSFH_CAGNH_fit.csv"))
dsilva_csmh_func = approxfun(
  dsilva25$z, 
  10^dsilva25$CSMHQ50
)
LL_csmh = function(p){
  
  -1*sum(dnorm(
    x = log10( csmh$Q50 ),
    mean = p[1] + log10(dsilva_csmh_func(zmids)),
    sd = csmh$ERR/(log(10) * csmh$Q50), 
    log = TRUE
  ))
  
}

norm_csmh = optimise(f = LL_csmh, interval = c(-2,1))

LSS_corr =dsilva_csmh_func(zmids) / csmh$Q50 * 10^norm_csmh$minimum
LSS_corr[lbt_mids >= 8] = 1.0

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
h5write(obj = cgmh, file = h5file, name = "cosmic/Mgas")
h5write(obj = cdmh_no_corr, file = h5file, name = "cosmic/MdustNoCorr")
h5write(obj = cdmh_wAGN, file = h5file, name = "cosmic/MdustwAGN")
h5write(obj = cgmh_wAGN, file = h5file, name = "cosmic/MgaswAGN")
h5write(obj = cdmh_wAGN_no_corr, file = h5file, name = "cosmic/MdustwAGNNoCorr")

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
h5createGroup(h5file, "vmax/Mgas")
for(i in 1:length(gas_mass_density)){
  h5write(
    obj = gas_mass_density[[i]]$vmax, 
    file = h5file,
    name = paste0("vmax/Mgas/zb",i)
  )
}
h5createGroup(h5file, "vmax/MgaswAGN")
for(i in 1:length(gas_mass_density_wAGN)){
  h5write(
    obj = gas_mass_density_wAGN[[i]]$vmax, 
    file = h5file,
    name = paste0("vmax/MgaswAGN/zb",i)
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
h5createGroup(h5file, "fit/Mgas")
for(i in 1:length(dust_mass_density)){
  h5write(
    obj = gas_mass_density[[i]]$fit, 
    file = h5file,
    name = paste0("fit/Mgas/zb",i)
  )
}
h5createGroup(h5file, "fit/MgaswAGN")
for(i in 1:length(gas_mass_density_wAGN)){
  h5write(
    obj = gas_mass_density_wAGN[[i]]$fit, 
    file = h5file,
    name = paste0("fit/MgaswAGN/zb",i)
  )
}

h5delete(h5file, "par")
h5createGroup(h5file, "par")
h5write(
  obj = dmf_par, 
  file = h5file,
  name = "par/DMF"
)
h5write(
  obj = dmf_wAGN_par, 
  file = h5file,
  name = "par/DMFwAGN"
)
h5write(
  obj = gmf_par, 
  file = h5file,
  name = "par/GMF"
)
h5write(
  obj = gmf_wAGN_par, 
  file = h5file,
  name = "par/GMFwAGN"
)
h5write(
  obj = smf_par, 
  file = h5file,
  name = "par/SMF"
)

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

save.image(
  file = "~/Documents/DustMassDensity/data/gama_devils_dust_mass.Rdata"
)

# dust_mass_density_wAGN_no_corr_noERR = foreach(i = 1:(length(zbins)-1)) %do% {
#   
#   # png(paste0("~/Documents/DustMassDensity/plots/dmf-AGN/lbt_",lbt_mids[i],".png"))
#   GAMA = compute_mass_function(
#     zlo = zbins[i],
#     zhi = zbins[i+1],
#     z = sort_match_AGN$z,
#     x = sort_match_AGN$dustmass.total,
#     x_err = 0.5 * (sort_match_AGN$dustmass.total_UB - sort_match_AGN$dustmass.total_LB),
#     areas = rep(217.54, length(sort_match_AGN$z)),
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
#     z = devilsd10_AGN$z,
#     x = devilsd10_AGN$dustmass.total,
#     x_err = 0.5 * (devilsd10_AGN$dustmass.total_UB - devilsd10_AGN$dustmass.birth_LB),
#     areas = rep(1.5, length(devilsd10_AGN$z)),
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
#   fit = compute_mass_function(
#     zlo = zbins[i],
#     zhi = zbins[i+1],
#     
#     ## Not gonna use these 
#     z = devilsd10_noAGN$z,
#     x = devilsd10_noAGN$dustmass.total,
#     x_err = 0.5 * (devilsd10_AGN$dustmass.total_UB - devilsd10_AGN$dustmass.birth_LB),
#     areas = devilsd10_noAGN$area,
#     
#     sm_bins = sm_bins,
#     vmax_bins = vmax_combine[,1],
#     vmax = vmax_combine[,2],
#     vmaxErr = sqrt( vmax_combine[,3]^2 ),
#     meanxErr = vmax_combine[,4],
#     errFloor = 0.1,
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
#   dev.off()
#   return(fit)
# }
# dust_birth_mass_density_wAGN_no_corr_noERR = foreach(i = 1:(length(zbins)-1)) %do% {
#   
#   # png(paste0("~/Documents/DustMassDensity/plots/dmf-AGN/lbt_",lbt_mids[i],".png"))
#   GAMA = compute_mass_function(
#     zlo = zbins[i],
#     zhi = zbins[i+1],
#     z = sort_match_AGN$z,
#     x = sort_match_AGN$dustmass.birth,
#     x_err = 0.5 * (sort_match_AGN$dustmass.birth_UB - sort_match_AGN$dustmass.birth_LB),
#     areas = rep(217.54, length(sort_match_AGN$z)),
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
#     z = devilsd10_AGN$z,
#     x = devilsd10_AGN$dustmass.birth,
#     x_err = 0.5 * (devilsd10_AGN$dustmass.birth_UB - devilsd10_AGN$dustmass.birth_LB),
#     areas = rep(1.5, length(devilsd10_AGN$z)),
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
#   fit = compute_mass_function(
#     zlo = zbins[i],
#     zhi = zbins[i+1],
#     
#     ## Not gonna use these 
#     z = devilsd10_noAGN$z,
#     x = devilsd10_noAGN$dustmass.birth,
#     x_err = 0.5 * (devilsd10_AGN$dustmass.birth_UB - devilsd10_AGN$dustmass.birth_LB),
#     areas = devilsd10_noAGN$area,
#     
#     sm_bins = sm_bins,
#     vmax_bins = vmax_combine[,1],
#     vmax = vmax_combine[,2],
#     vmaxErr = sqrt( vmax_combine[,3]^2 ),
#     meanxErr = vmax_combine[,4],
#     errFloor = 0.1,
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
#   dev.off()
#   return(fit)
# }
# dust_screen_mass_density_wAGN_no_corr_noERR = foreach(i = 1:(length(zbins)-1)) %do% {
#   
#   # png(paste0("~/Documents/DustMassDensity/plots/dmf-AGN/lbt_",lbt_mids[i],".png"))
#   GAMA = compute_mass_function(
#     zlo = zbins[i],
#     zhi = zbins[i+1],
#     z = sort_match_AGN$z,
#     x = sort_match_AGN$dustmass.screen,
#     x_err = 0.5 * (sort_match_AGN$dustmass.screen_UB - sort_match_AGN$dustmass.screen_LB),
#     areas = rep(217.54, length(sort_match_AGN$z)),
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
#     z = devilsd10_AGN$z,
#     x = devilsd10_AGN$dustmass.screen,
#     x_err = 0.5 * (devilsd10_AGN$dustmass.screen_UB - devilsd10_AGN$dustmass.screen_LB),
#     areas = rep(1.5, length(devilsd10_AGN$z)),
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
#   fit = compute_mass_function(
#     zlo = zbins[i],
#     zhi = zbins[i+1],
#     
#     ## Not gonna use these 
#     z = devilsd10_noAGN$z,
#     x = devilsd10_noAGN$dustmass.screen,
#     x_err = 0.5 * (devilsd10_AGN$dustmass.screen_UB - devilsd10_AGN$dustmass.screen_LB),
#     areas = devilsd10_noAGN$area,
#     
#     sm_bins = sm_bins,
#     vmax_bins = vmax_combine[,1],
#     vmax = vmax_combine[,2],
#     vmaxErr = sqrt( vmax_combine[,3]^2 ),
#     meanxErr = vmax_combine[,4],
#     errFloor = 0.1,
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
#   dev.off()
#   return(fit)
# }
