library(Rfits)
library(celestial)
library(magicaxis)
library(scales)
library(rhdf5)
library(matrixStats)
library(data.table)
library(foreach)
library(ProSpect)
library(doParallel)
library(Highlander)

list2env(readRDS("~/Documents/DustMassDensity/data/new_M2L_data.rds"), envir = .GlobalEnv)

magphys_cat = Rfits_read_table("~/Documents/DustMassDensity/data/GAMA/MagPhysv06.fits")
prospect_cat = Rfits_read_table("~/Documents/DustMassDensity/data/GAMA/ProSpectv03.fits")

lambdar_phot = Rfits_read_table("~/Documents/DustMassDensity/data/GAMA/LambdarCatv01.fits")
profoundgkv_phot_raw = Rfits_read_table("~/Documents/DustMassDensity/data/GAMA/gkvScienceCatv02.fits")
profoundIR_phot_raw = Rfits_read_table("~/Documents/DustMassDensity/data/GAMA/gkvFarIRv03.fits")
profoundIR_phot = profoundIR_phot_raw[profoundIR_phot_raw$uberID %in% profoundgkv_phot_raw$uberID, ]
profoundgkv_phot = profoundgkv_phot_raw[profoundgkv_phot_raw$uberID %in% profoundIR_phot$uberID, ]
message(all(profoundIR_phot$uberID == profoundgkv_phot$uberID)) ## Check that they are matched
profound_phot = cbind(profoundgkv_phot, profoundIR_phot)

dim(profound_phot)
dim(prospect_cat)

dim(lambdar_phot)
dim(magphys_cat)

## Only get photometry for those that were SED fitted
lambdar_phot_match = lambdar_phot[lambdar_phot$CATAID %in% magphys_cat$CATAID, ]
message(all(lambdar_phot_match$CATAID == magphys_cat$CATAID)) ## Check that they are matched
lambdar_magphys = cbind(lambdar_phot_match, magphys_cat)

profound_phot_match = profound_phot[profound_phot$uberID %in% prospect_cat$uberID, ]
message(all(profound_phot_match$uberID == prospect_cat$uberID)) ## Check that they are matched
message(all(profound_phot_match$CATAID == prospect_cat$CATAID)) ## Check that they are matched
profound_prospect = cbind(profound_phot_match, prospect_cat)

lambdar_idx = lambdar_magphys$P100_flux/lambdar_magphys$P100_fluxerr > 5 & 
  lambdar_magphys$P160_flux/lambdar_magphys$P160_fluxerr > 5 & 
  lambdar_magphys$S250_flux/lambdar_magphys$S250_fluxerr > 5 &
  lambdar_magphys$S350_flux/lambdar_magphys$S350_fluxerr > 5 &
  lambdar_magphys$S500_flux/lambdar_magphys$S500_fluxerr > 5

profound_idx = profound_prospect$flux_PSF_p100/profound_prospect$flux_PSF_err_p100 > 5 & 
  profound_prospect$flux_PSF_p160/profound_prospect$flux_PSF_err_p160 > 5 & 
  profound_prospect$flux_PSF_s250/profound_prospect$flux_PSF_err_s250 > 5 &
  profound_prospect$flux_PSF_s350/profound_prospect$flux_PSF_err_s350 > 5 &
  profound_prospect$flux_PSF_s500/profound_prospect$flux_PSF_err_s500 > 5

lambdar_magphys_trim = lambdar_magphys[lambdar_idx,]
profound_prospect_trim = profound_prospect[profound_idx,]

idx_coordmatch = coordmatch(
  coordref = profound_prospect_trim[, c("RAmax", "Decmax")],
  coordcompare = lambdar_magphys_trim[, c("RA", "DEC")]
)

lambdar_magphys_match = data.frame(lambdar_magphys_trim[idx_coordmatch$bestmatch$compareID, ])
profound_prospect_match = data.frame(profound_prospect_trim[idx_coordmatch$bestmatch$refID, ])

df = data.frame(
  "zProSpect" = profound_prospect_match$Z,
  "MdustProSpect" = profound_prospect_match$DustMass_50,
  "MdustErrProSpect" = 0.5*(profound_prospect_match$DustMass_84 - profound_prospect_match$DustMass_16),
  "LdustProSpect" = profound_prospect_match$DustLum_50,
  "LdustErrProSpect" = 0.5*(profound_prospect_match$DustLum_84 - profound_prospect_match$DustLum_16),
  "RAProSpect" = profound_prospect_match$RAmax,
  "DECProSpect" = profound_prospect_match$Decmax,
  
  "zMagphys" = lambdar_magphys_match$Z,
  "MdustMagphys" = 10^lambdar_magphys_match$mass_dust_percentile50,
  "MdustErrMagphys" = sqrt( (10^lambdar_magphys_match$mass_dust_percentile50 * log(10) * 0.5 * (lambdar_magphys_match$mass_dust_percentile84 - lambdar_magphys_match$mass_dust_percentile16))^2 ),
  "RAMagphys" = lambdar_magphys_match$RA,
  "DECMagphys" = lambdar_magphys_match$DEC,
  data.frame(lambdar_magphys_match)[,grep("^(mass_dust|L_dust|f_mu_IR|xi|T).*best_fit", names(lambdar_magphys_match), value = TRUE, perl = TRUE)]
)

magplot(
  log10(profound_prospect_match$DustLum_50),
  lambdar_magphys_match$L_dust_percentile50
)
abline(0,1)

wavelength = c(
  1539, 2316, 3528, 4760, 6326, 7599, 6654, 10229, 12556, 16499, 21571, 
  3.4e4, 4.65e4, 12.8e4, 22.4e4, 98.9e4, 156e4, 249e4, 350e4, 504e4
)

profound_fluxes = profound_prospect_match[, grep("^flux_(?!err)(?!.*_err)(?!.*l$).*", names(profound_prospect_match), perl = TRUE, value = TRUE)]
profound_flux_errs = profound_prospect_match[, grep("err", names(profound_prospect_match), perl = TRUE, value = TRUE)]
profound_fluxes[profound_fluxes < 0] = 0
profound_flux_errs = sqrt( profound_flux_errs^2 + (0.1 * profound_fluxes)^2 )

lambdar_fluxes = lambdar_magphys_match[, grep("^(?!z_flux$).*flux(?!err)", names(lambdar_magphys_match), perl = TRUE, value = TRUE)]
lambdar_flux_errs = lambdar_magphys_match[, grep("^(?!z_fluxerr$).*err", names(lambdar_magphys_match), perl = TRUE, value = TRUE)]
lambdar_fluxes[lambdar_fluxes < 0] = 0
lambdar_flux_errs = sqrt( lambdar_flux_errs^2 + (0.1 * lambdar_fluxes)^2 )

profound_weighted_stack = colSums(as.matrix(profound_fluxes) * as.matrix(profound_flux_errs^-2), na.rm = TRUE) / colSums(as.matrix(profound_flux_errs^-2), na.rm = TRUE)
lambdar_weighted_stack = colSums(as.matrix(lambdar_fluxes) * as.matrix(lambdar_flux_errs^-2), na.rm = TRUE) / colSums(as.matrix(lambdar_flux_errs^-2), na.rm = TRUE)
profound_weighted_err = sqrt(1/colSums(as.matrix(profound_flux_errs^-2), na.rm = TRUE))
lambdar_weighted_err = sqrt(1/colSums(as.matrix(lambdar_flux_errs^-2), na.rm = TRUE))

df_spec = data.frame(
  "wavelength" = wavelength,
  "ProFoundStack" = profound_weighted_stack,
  "ProFoundErr" = profound_weighted_err,
  "LambdarStack" = lambdar_weighted_stack,
  "LambdarErr" = lambdar_weighted_err
)

# magplot(
#   lambdar_magphys_match$mass_dust_percentile50, 
#   log10(profound_prospect_match$DustMass_50), 
#   xlab = "Magphys Dust mass",
#   ylab = "ProSpect Dust mass"
# )
# magplot(
#   lambdar_magphys_match$L_dust_percentile50, 
#   log10(profound_prospect_match$DustLum_50), 
#   xlab = "Magphys Dust mass",
#   ylab = "ProSpect Dust mass"
# )
# abline(0,1)
# 
# magplot(
#   lambdar_magphys_match$RA,
#   lambdar_magphys_match$DEC
# )
# points(
#   profound_prospect_match$RAmax,
#   profound_prospect_match$Decmax,
#   pch = 16, cex = 0.5
# )

## Fit with single component grey body?
filtout = list(
  approxfun(getfilt("FUV"), yleft = 0, yright = 0),
  approxfun(getfilt("NUV"), yleft = 0, yright = 0),
  approxfun(getfilt("u_VST"), yleft = 0, yright = 0),
  approxfun(getfilt("g_VST"), yleft = 0, yright = 0),
  approxfun(getfilt("r_VST"), yleft = 0, yright = 0),
  approxfun(getfilt("i_VST"), yleft = 0, yright = 0),
  approxfun(getfilt("Z_VISTA"), yleft = 0, yright = 0),
  approxfun(getfilt("Y_VISTA"), yleft = 0, yright = 0),
  approxfun(getfilt("J_VISTA"), yleft = 0, yright = 0),
  approxfun(getfilt("H_VISTA"), yleft = 0, yright = 0),
  approxfun(getfilt("K_VISTA"), yleft = 0, yright = 0),
  approxfun(getfilt("W1_WISE"), yleft = 0, yright = 0),
  approxfun(getfilt("W2_WISE"), yleft = 0, yright = 0),
  approxfun(getfilt("W3_WISE"), yleft = 0, yright = 0),
  approxfun(getfilt("W4_WISE"), yleft = 0, yright = 0),
  approxfun(getfilt("P100_Herschel"), yleft = 0, yright = 0),
  approxfun(getfilt("P160_Herschel"), yleft = 0, yright = 0),
  approxfun(getfilt("S250_Herschel"), yleft = 0, yright = 0),
  approxfun(getfilt("S350_Herschel"), yleft = 0, yright = 0),
  approxfun(getfilt("S500_Herschel"), yleft = 0, yright = 0)
)
wave_grey_body = 10^seq(3, 8, 0.01)
grey_body_fits = foreach(kk = 1:dim(profound_prospect_match)[1], .combine = rbind) %do% {
  if(kk %% 10 == 0){message(kk)}
  
  # Ldust = 10^magphys_final_sample$L_dust_percentile50[kk]
  # Ldust = prospect_final_sample$DustLum_50[kk]
  fit_idx = 16:20
  
  # ff = lambdar_fluxes[kk, ]
  # fferr = lambdar_flux_errs[kk, ]
  ff = profound_fluxes[kk, ]
  fferr = profound_flux_errs[kk, ]
  
  LL = function(parm){
    
    Temp = parm[1]
    Ldust = parm[2]
    
    rest_frame = greybody_norm(
      wave = wave_grey_body, 
      Temp = Temp, 
      beta = 2.0, 
      z = profound_prospect_match$Z[kk],
      norm = 10^Ldust
    )
    
    observed_frame = photom_lum(
      wave = wave_grey_body,
      lum = rest_frame,
      outtype = "Jy",
      filters = filtout[fit_idx],
      z = profound_prospect_match$Z[kk],
      ref = "Planck15"
    )
    
    likelihood = sum(dnorm(
      x = unlist(ff[fit_idx]),
      mean = observed_frame,
      sd = unlist(unlist(fferr[fit_idx])),
      log = TRUE
    ), na.rm = TRUE)
    return(likelihood)
  }
  
  opt = optim(
    par = c(15, 9),
    fn = LL,
    method = "L-BFGS-B",
    lower = c(5, 4),
    upper = c(150, 15),
    control = list(fnscale = -1),
    hessian = TRUE
  )
  
  rest_frame = greybody_norm(
    wave = wave_grey_body, 
    Temp = opt$par[1], 
    beta = 2.0, 
    z = profound_prospect_match$Z[kk],
    norm = 10^opt$par[2]
  )
  observed_frame = photom_lum(
    wave = wave_grey_body,
    lum = rest_frame,
    outtype = "Jy",
    filters = filtout[fit_idx],
    z = profound_prospect_match$Z[kk],
    ref = "Planck15"
  )
  observed_lambdaflux = Lum2Flux(
    wave = wave_grey_body,
    lum = rest_frame,
    z = profound_prospect_match$Z[kk], 
    ref = "Planck15"
  )
  observed_freqflux = convert_wave2freq(
    flux_wave = observed_lambdaflux$flux,
    wave = observed_lambdaflux$wave
  )
  
  png(paste0("~/Documents/DustMassDensity/plots/fit_greybody/", kk, ".png"))
  magplot(
    wavelength,
    unlist(ff),
    ylim = c(1e-5, 10),
    log = "xy",
    pch = 16,
    col = "red",
  )
  magerr(
    wavelength,
    unlist(ff),
    ylo = unlist(fferr),
    col = "red"
  )
  points(
    wavelength[fit_idx],
    observed_frame,
    pch = 1, 
    cex = 1.5
  )
  lines(
    observed_lambdaflux$wave,
    CGS2Jansky(observed_freqflux)
  )
  dev.off()
  
  T_best = opt$par[1]
  L_best = opt$par[2]
  err = sqrt(diag(-solve(opt$hessian)))
  T_err = err[1]
  L_err = err[2]
  
  grey_body_best = greybody(
    wave = wave_grey_body, 
    Temp = T_best, 
    beta = 2.0, 
    k850 = 0.077
  )
  Mdust_best = 10^L_best / sum(c(0,diff(wave_grey_body))*grey_body_best)
  # Mdust_best = Ldust / trapz(x = wave_grey_body, grey_body_best)
  
  # trapz(x = wave_grey_body, grey_body_best)
  
  Mdust_samples = sapply(1:1000, function(x){
    T_sample = rnorm(n = 1, mean = T_best, sd = T_err)
    L_sample = rnorm(n = 1, mean = L_best, sd = L_err)
    grey_body_sample = greybody(
      wave = wave_grey_body, 
      Temp = T_sample, 
      beta = 2.0, 
      k850 = 0.077
    )
    Mdust_sample = 10^L_sample / sum(c(0,diff(wave_grey_body))*grey_body_sample)
    return(Mdust_sample)
  })
  
  MdustQ16 = quantile(Mdust_samples, 0.16)
  MdustQ84 = quantile(Mdust_samples, 0.84)
  MdustErr = 0.5 * (MdustQ84 - MdustQ16)

  ret_ = data.frame(
    "Ldust" = L_best,
    "LdustErr" = L_err,
    "Mdust" = Mdust_best,
    "MdustErr" = MdustErr,
    "Tdust" = T_best,
    "TdustErr" = T_err
  )
  return(ret_)
}
black_body_fits = foreach(kk = 1:dim(profound_prospect_match)[1], .combine = rbind) %do% {
  if(kk %% 10 == 0){message(kk)}
  
  # Ldust = 10^magphys_final_sample$L_dust_percentile50[kk]
  # Ldust = prospect_final_sample$DustLum_50[kk]
  fit_idx = 16:20
  
  # ff = lambdar_fluxes[kk, ]
  # fferr = lambdar_flux_errs[kk, ]
  ff = profound_fluxes[kk, ]
  fferr = profound_flux_errs[kk, ]
  
  LL = function(parm){
    
    Temp = parm[1]
    Ldust = parm[2]
    
    rest_frame = blackbody_norm(
      wave = wave_grey_body, 
      Temp = Temp, 
      z = profound_prospect_match$Z[kk],
      norm = 10^Ldust
    )
    
    observed_frame = photom_lum(
      wave = wave_grey_body,
      lum = rest_frame,
      outtype = "Jy",
      filters = filtout[fit_idx],
      z = profound_prospect_match$Z[kk],
      ref = "Planck15"
    )
    
    likelihood = sum(dnorm(
      x = unlist(ff[fit_idx]),
      mean = observed_frame,
      sd = unlist(unlist(fferr[fit_idx])),
      log = TRUE
    ), na.rm = TRUE)
    return(likelihood)
  }
  
  opt = optim(
    par = c(15, 9),
    fn = LL,
    method = "L-BFGS-B",
    lower = c(5, 4),
    upper = c(150, 15),
    control = list(fnscale = -1),
    hessian = TRUE
  )
  
  rest_frame = blackbody_norm(
    wave = wave_grey_body, 
    Temp = opt$par[1], 
    z = profound_prospect_match$Z[kk],
    norm = 10^opt$par[2]
  )
  observed_frame = photom_lum(
    wave = wave_grey_body,
    lum = rest_frame,
    outtype = "Jy",
    filters = filtout[fit_idx],
    z = profound_prospect_match$Z[kk],
    ref = "Planck15"
  )
  observed_lambdaflux = Lum2Flux(
    wave = wave_grey_body,
    lum = rest_frame,
    z = profound_prospect_match$Z[kk], 
    ref = "Planck15"
  )
  observed_freqflux = convert_wave2freq(
    flux_wave = observed_lambdaflux$flux,
    wave = observed_lambdaflux$wave
  )
  
  png(paste0("~/Documents/DustMassDensity/plots/fit_blackbody/", kk, ".png"))
  magplot(
    wavelength,
    unlist(ff),
    ylim = c(1e-5, 10),
    log = "xy",
    pch = 16,
    col = "red",
  )
  magerr(
    wavelength,
    unlist(ff),
    ylo = unlist(fferr),
    col = "red"
  )
  points(
    wavelength[fit_idx],
    observed_frame,
    pch = 1, 
    cex = 1.5
  )
  lines(
    observed_lambdaflux$wave,
    CGS2Jansky(observed_freqflux)
  )
  dev.off()
  
  T_best = opt$par[1]
  L_best = opt$par[2]
  err = sqrt(diag(-solve(opt$hessian)))
  T_err = err[1]
  L_err = err[2]
  
  grey_body_best = blackbody(
    wave = wave_grey_body, 
    Temp = T_best, 
    k850 = 0.077
  )
  Mdust_best = 10^L_best / sum(c(0,diff(wave_grey_body))*grey_body_best)
  # Mdust_best = Ldust / trapz(x = wave_grey_body, grey_body_best)
  
  # trapz(x = wave_grey_body, grey_body_best)
  
  Mdust_samples = sapply(1:1000, function(x){
    T_sample = rnorm(n = 1, mean = T_best, sd = T_err)
    L_sample = rnorm(n = 1, mean = L_best, sd = L_err)
    grey_body_sample = blackbody(
      wave = wave_grey_body, 
      Temp = T_sample, 
      k850 = 0.077
    )
    Mdust_sample = 10^L_sample / sum(c(0,diff(wave_grey_body))*grey_body_sample)
    return(Mdust_sample)
  })
  
  MdustQ16 = quantile(Mdust_samples, 0.16)
  MdustQ84 = quantile(Mdust_samples, 0.84)
  MdustErr = 0.5 * (MdustQ84 - MdustQ16)
  
  ret_ = data.frame(
    "Ldust" = L_best,
    "LdustErr" = L_err,
    "Mdust" = Mdust_best,
    "MdustErr" = MdustErr,
    "Tdust" = T_best,
    "TdustErr" = T_err
  )
  return(ret_)
}

par(mfrow = c(1,2))
magplot(
  10^grey_body_fits$Ldust,
  10^lambdar_magphys_match$L_dust_percentile50,
  xlab = "Grey body Dust luminosity [Lsun]",
  ylab = "Magphys/ProSpect luminosity [Lsun]",
  col = alpha("black", 0.5),
  log = "xy",
  xlim = c(1e7,1e13),
  ylim = c(1e7, 1e13)
)
points(
  10^grey_body_fits$Ldust,
  profound_prospect_match$DustLum_50,
  pch = 16, 
  col = alpha("red", 0.5)
)
legend(
  x = "bottomright",
  pch = c(1,16),
  col = c("black", "red"), 
  legend = c("Magphys", "ProSpect")
)
abline(0,1)
maghist(
  log10(10^grey_body_fits$Ldust/profound_prospect_match$DustLum_50), 
  xlab = "Grey body - ProSpect [dex]"
)

magplot(
  grey_body_fits$Mdust, 
  10^lambdar_magphys_match$mass_dust_percentile50, 
  xlab = "Grey body Dust Mass [Msun]",
  ylab = "Magphys/ProSpect Dust Mass [Msun]",
  col = alpha("black", 0.5),
  log = "xy",
  xlim = c(1e5,1e10),
  ylim = c(1e5, 1e10)
)
points(
  grey_body_fits$Mdust,
  profound_prospect_match$DustMass_50, 
  pch = 16, 
  col = alpha("red", 0.5)
)
legend(
  x = "bottomright",
  pch = c(1,16),
  col = c("black", "red"), 
  legend = c("Magphys", "ProSpect")
)
abline(0,1)
maghist(
  log10((grey_body_fits$Mdust)/profound_prospect_match$DustMass_50), 
  xlab = "Grey body - ProSpect [dex]", 
  breaks = seq(-2.5, 1.0, 0.1), 
  ylim = c(0, 80),
  col = alpha("black", 0.4), 
)
maghist(
  log10((grey_body_fits$Mdust)/10^lambdar_magphys_match$mass_dust_percentile50), 
  add = TRUE, 
  breaks = seq(-2.5, 1.0, 0.1),
  col = alpha("blue", 0.4)
)

# magplot(
#   NA,
#   log = "xy",
#   ylim = c(1e-5, 0.5),
#   xlim = c(1e3, 1e7),
#   pch = 16
# )
# for(kk in 1:dim(profound_fluxes)[1]){
#   points(wavelength, unlist(profound_fluxes[kk,]), col = alpha("darkorange",0.1), pch = 16)
#   magerr(wavelength, unlist(profound_fluxes[kk,]), ylo = unlist(unlist(profound_flux_errs[kk,])), col = alpha("darkorange",0.1))
#   points(wavelength, unlist(lambdar_fluxes[kk,]), col = alpha("darkred",0.1), pch = 16)
#   magerr(wavelength, unlist(lambdar_fluxes[kk,]), ylo = unlist(unlist(lambdar_flux_errs[kk,])), col = alpha("darkred",0.1))
# }
# points(
#   wavelength, profound_weighted_stack, pch = 16, col = "darkorange"
# )
# magerr(
#   wavelength, profound_weighted_stack, ylo = profound_weighted_err, col = "darkorange"
# )
# points(
#   wavelength, lambdar_weighted_stack, pch =16, col = "darkred"
# )
# magerr(
#   wavelength, lambdar_weighted_stack, ylo = lambdar_weighted_err, col = "darkred"
# )

## Refit stuff
profound_Data = list(
  "flux" = profound_fluxes,
  "flux_err" = profound_flux_errs,
  "redshifts" = profound_prospect_match$Z,
  "uberID" = profound_prospect_match$uberID,
  "catID" = profound_prospect_match$CATAID,
  "RA" = profound_prospect_match$RAmax,
  "Dec" = profound_prospect_match$Decmax
  # "filtout" = filtout
)
saveRDS(profound_Data, "~/Documents/DustMassDensity/Pawsey/save_highSN_sample.rds")

prospect_fnames = list.files("/Users/22252335/Documents/DustMassDensity/Pawsey/out/", full.names = TRUE)
prospect_files = lapply(prospect_fnames, function(x){
  foo = readRDS(x)
  # plot(foo$bestfit$SEDout)
  bar = c(colQuantiles(as.matrix(foo$parm_sample), probs = c(0.16, 0.5, 0.84)), foo$RA, foo$Dec, colQuantiles(as.matrix(foo$highlander$LD_last$Posterior1[, c("alpha_SF_screen", "alpha_SF_birth", "Zfinal")]), probs = c(0.16, 0.5, 0.84)))
  names(bar) = c(paste0(names(foo$parm_sample),"Q16"), paste0(names(foo$parm_sample),"Q50"), paste0(names(foo$parm_sample),"Q84"), "RA", "Dec", paste0(c("alpha_screen", "alpha_birth", "Zfinal"),"Q16"), paste0(c("alpha_screen", "alpha_birth", "Zfinal"),"Q50"), paste0(c("alpha_screen", "alpha_birth", "Zfinal"),"Q84"))
  return(bar)
})
prospect_refit = data.frame(do.call(rbind, prospect_files))

prospect_coord_match = coordmatch(
  coordref = profound_prospect_match[, c("RAmax", "Decmax")], 
  coordcompare = prospect_refit[, c("RA", "Dec")]
)
prospect_refit$stub = prospect_fnames
prospect_refit_match = prospect_refit[prospect_coord_match$bestmatch$compareID, ]

magplot(
  prospect_refit_match$DustLumQ50,
  profound_prospect_match$DustLum_50,
  log = "xy"
)
magerr(
  prospect_refit_match$DustLumQ50,
  profound_prospect_match$DustLum_50,
  xlo = (prospect_refit_match$DustLumQ84 - prospect_refit_match$DustLumQ16) * 0.5,
  ylo = (profound_prospect_match$DustLum_84 - profound_prospect_match$DustLum_16) * 0.5,
)
abline(0,1)

magplot(
  prospect_refit_match$DustMassQ50,
  profound_prospect_match$DustMass_50,
  log = "xy"
)
magerr(
  prospect_refit_match$DustMassQ50,
  profound_prospect_match$DustMass_50,
  xlo = (prospect_refit_match$DustMassQ84 - prospect_refit_match$DustMassQ16) * 0.5,
  ylo = (profound_prospect_match$DustMass_84 - profound_prospect_match$DustMass_16) * 0.5,
)
abline(0,1)

# RR14_BPL = readRDS("~/Documents/DustMassDensity/data/RR14_BPL.rds")
new_dust_masses_solar = lapply(prospect_refit_match$stub, function(x){
  message(x)
  foo = readRDS(x)
  # plot(foo$bestfit$SEDout)
  # bar = c(colQuantiles(as.matrix(foo$parm_sample), probs = c(0.16, 0.5, 0.84)), foo$RA, foo$Dec, colQuantiles(as.matrix(foo$highlander$LD_last$Posterior1[, c("alpha_SF_screen", "alpha_SF_birth", "Zfinal")]), probs = c(0.16, 0.5, 0.84)))
  # names(bar) = c(paste0(names(foo$parm_sample),"Q16"), paste0(names(foo$parm_sample),"Q50"), paste0(names(foo$parm_sample),"Q84"), "RA", "Dec", paste0(c("alpha_screen", "alpha_birth", "Zfinal"),"Q16"), paste0(c("alpha_screen", "alpha_birth", "Zfinal"),"Q50"), paste0(c("alpha_screen", "alpha_birth", "Zfinal"),"Q84"))
  
  alpha_screen = foo$highlander$LD_last$Posterior1[, "alpha_SF_screen"]
  alpha_birth = foo$highlander$LD_last$Posterior1[, "alpha_SF_birth"]
  
  lum_screen = foo$parm_sample$DustLumScreen
  lum_birth = foo$parm_sample$DustLumBirth
  
  # Zfinal = 10^foo$highlander$LD_last$Posterior1[, "Zfinal"]
  # RR14_DTH = RR14_BPL(Z = Zfinal, doDTG = FALSE)
  
  new_dust_mass_solar = (lum_screen/Dale_vM2L_func(alpha_screen) + lum_birth/Dale_vM2L_func(alpha_birth))
  
  bar = c(
    "NewDustMassQ50" = median(new_dust_mass_solar),
    "NewDustMassQ16" = quantile(new_dust_mass_solar, 0.16),
    "NewDustMassQ84" = quantile(new_dust_mass_solar, 0.84),
    "NewDustMassErr" = 0.5 * (quantile(new_dust_mass_solar, 0.84) - quantile(new_dust_mass_solar, 0.16))
  )
  
  return(bar)
})
new_dust_masses_solar_DF = data.frame(do.call(rbind, new_dust_masses_solar))
names(new_dust_masses_solar_DF) = c("NewDustMassQ50","NewDustMassQ16","NewDustMassQ84","NewDustMassErr")

print(
  median(prospect_refit_match$DustMassQ50/new_dust_masses_solar_DF$NewDustMassQ50)
)
print(
  median(profound_prospect_match$DustMass_50/new_dust_masses_solar_DF$NewDustMassQ50)
)
print(
  log10(median(10^lambdar_magphys_match$mass_dust_percentile50/new_dust_masses_solar_DF$NewDustMassQ50))
)
Md_corr = h5read('~/Documents/DustMassDensity/data/all_data.h5', "Md_corr")
print(
  Md_corr
)

magplot(
  new_dust_masses_solar_DF$NewDustMassQ50,
  profound_prospect_match$DustMass_50,
  log = "xy",
  xlim = c(1e5,1e10),
  ylim = c(1e5, 1e10)
)
magerr(
  new_dust_masses_solar_DF$NewDustMassQ50,
  profound_prospect_match$DustMass_50,
  xlo = new_dust_masses_solar_DF$NewDustMassErr,
  ylo = (profound_prospect_match$DustMass_84-profound_prospect_match$DustMass_16)*0.5
)
points(
  10^lambdar_magphys_match$mass_dust_percentile50,
  profound_prospect_match$DustMass_50, 
  col = alpha("red", 0.5),
  pch = 16
)
points(
  grey_body_fits$Mdust*1.1*1.6,
  profound_prospect_match$DustMass_50, 
  pch = 16, 
  col = alpha("blue", 0.5)
)
abline(0,1)


# maghist(
#   log10(
#     grey_body_fits$Mdust*1.1*1.6 / new_dust_masses_solar_DF$NewDustMassQ50
#   ), 
#   xlab = "Grey body dust mass \n scaled up for total luminosity and total mass in dust",
# )
# maghist(
#   log10(
#     10^lambdar_magphys_match$mass_dust_percentile50 / new_dust_masses_solar_DF$NewDustMassQ50
#   ), 
#   xlab = "Grey body dust mass \n scaled up for total luminosity and total mass in dust",
# )

magplot(
  grey_body_fits$Mdust*1.14*1.6,
  new_dust_masses_solar_DF$NewDustMassQ50,
  ylab = "New dust masses from ProSpect \n using the new variable DTH",
  xlab = "Grey body dust mass \n scaled up for total luminosity and total mass in dust",
  log = "xy",
  xlim = c(1e5,1e10),
  ylim = c(1e5, 1e10)
)
abline(0,1)
legend(
  x = "bottomright",
  legend = c(
   paste0("Median difference (x-y) = ", round(log10(median(grey_body_fits$Mdust*1.14*1.6 / new_dust_masses_solar_DF$NewDustMassQ50)),3), " dex")
  )
)

grey_body_fits_df = grey_body_fits
names(grey_body_fits_df) = paste0("greybody", names(grey_body_fits_df))

h5file = '~/Documents/DustMassDensity/data/all_data.h5'
h5delete(h5file, "Photometry")
h5createGroup(h5file, "Photometry")
h5write(
  obj = cbind(df, new_dust_masses_solar_DF, grey_body_fits_df), file = h5file, name = "Photometry/FinalSample"
)
h5write(
  obj = df_spec, file = h5file, name = "Photometry/WeightedSpec"
)

save.image("~/Documents/DustMassDensity/data/compare_prospect_magphys.Rdata")

