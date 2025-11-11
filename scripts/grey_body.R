library(ProSpect)
library(foreach)
library(doParallel)
library(magicaxis)
library(rhdf5)
library(data.table)

# load("~/Documents/DustMassDensity/data/compare_prospect_magphys.Rdata")
foo = readRDS("~/Documents/DustMassDensity/Pawsey/save_highSN_sample.rds")
list2env(foo, envir = .GlobalEnv)
## Look more closely at grey body fits of 218 FIR galaxies

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
wavelength = sapply( 
  list(
    getfilt("FUV"),
    getfilt("NUV"),
    getfilt("u_VST"),
    getfilt("g_VST"),
    getfilt("r_VST"),
    getfilt("i_VST"),
    getfilt("Z_VISTA"),
    getfilt("Y_VISTA"),
    getfilt("J_VISTA"),
    getfilt("H_VISTA"),
    getfilt("K_VISTA"),
    getfilt("W1_WISE"),
    getfilt("W2_WISE"),
    getfilt("W3_WISE"),
    getfilt("W4_WISE"),
    getfilt("P100_Herschel"), 
    getfilt("P160_Herschel"),
    getfilt("S250_Herschel"),
    getfilt("S350_Herschel"),
    getfilt("S500_Herschel")
  ), 
  cenwavefunc
)

wave_grey_body = 10^seq(3, 8, 0.01)

cores = 6
fit_grey_body = function(beta = 2, kappa = 0.077){
  registerDoParallel(cores = cores)
  grey_body_fits = foreach(kk = 1:dim(flux)[1], .combine = rbind) %dopar% {
    if(kk %% 10 == 0){message(kk)}
    
    # Ldust = 10^magphys_final_sample$L_dust_percentile50[kk]
    # Ldust = prospect_final_sample$DustLum_50[kk]
    fit_idx = 16:20
    
    # ff = lambdar_fluxes[kk, ]
    # fferr = lambdar_flux_errs[kk, ]
    ff = flux[kk, ]
    fferr = flux_err[kk, ]
    z = redshifts[kk]
    
    LL = function(parm){
      
      Temp = parm[1]
      Ldust = parm[2]
      
      rest_frame = greybody_norm(
        wave = wave_grey_body, 
        Temp = Temp, 
        beta = beta, 
        z = redshifts[kk],
        norm = 10^Ldust
      )
      
      observed_frame = photom_lum(
        wave = wave_grey_body,
        lum = rest_frame,
        outtype = "Jy",
        filters = filtout[fit_idx],
        z = z,
        ref = "Planck18"
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
      beta = beta, 
      z = z,
      norm = 10^opt$par[2]
    )
    observed_frame = photom_lum(
      wave = wave_grey_body,
      lum = rest_frame,
      outtype = "Jy",
      filters = filtout[fit_idx],
      z = z,
      ref = "Planck18"
    )
    observed_lambdaflux = Lum2Flux(
      wave = wave_grey_body,
      lum = rest_frame,
      z = z, 
      ref = "Planck18"
    )
    observed_freqflux = convert_wave2freq(
      flux_wave = observed_lambdaflux$flux,
      wave = observed_lambdaflux$wave
    )
    
    stub = paste0("~/Documents/DustMassDensity/plots/more_greybodies/", beta, "_", kappa, "/")
    dir.create(stub, recursive = TRUE)
    png(paste0(stub, kk, ".png"))
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
      beta = beta, 
      k850 = kappa
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
        beta = beta, 
        k850 = kappa
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
  stopImplicitCluster()
  return(grey_body_fits)
}

beta2.0_kappa0.077 = fit_grey_body()
beta3.0_kappa0.077 = fit_grey_body(beta = 3.0)
beta1.5_kappa0.077 = fit_grey_body(beta = 1.5)
beta1.0_kappa0.077 = fit_grey_body(beta = 1.0)
beta0.5_kappa0.077 = fit_grey_body(beta = 0.5)
beta2.0_kappa0.77 = fit_grey_body(kappa = 0.77)
beta2.0_kappa0.03 = fit_grey_body(kappa = 0.077/3.0)


maghist(
  beta3.0_kappa0.077$Ldust/beta2.0_kappa0.077$Ldust
)
maghist(
  beta1.5_kappa0.077$Ldust/beta2.0_kappa0.077$Ldust
)
maghist(
  beta1.0_kappa0.077$Ldust/beta2.0_kappa0.077$Ldust
)
maghist(
  beta0.5_kappa0.077$Ldust/beta3.0_kappa0.077$Ldust
)

maghist(
  beta3.0_kappa0.077$Mdust/beta2.0_kappa0.077$Mdust
)

maghist(
  beta1.5_kappa0.077$Mdust/beta2.0_kappa0.077$Mdust
)
maghist(
  beta1.0_kappa0.077$Mdust/beta2.0_kappa0.077$Mdust
)
maghist(
  beta0.5_kappa0.077$Mdust/beta2.0_kappa0.077$Mdust
)

maghist(
  beta0.5_kappa0.077$Mdust/beta3.0_kappa0.077$Mdust
)

maghist(
  beta3.0_kappa0.077$Tdust/beta2.0_kappa0.077$Tdust
)
maghist(
  beta1.5_kappa0.077$Tdust/beta2.0_kappa0.077$Tdust
)
maghist(
  beta1.0_kappa0.077$Tdust/beta2.0_kappa0.077$Tdust
)
maghist(
  beta0.5_kappa0.077$Tdust/beta3.0_kappa0.077$Tdust
)

maghist(
  beta2.0_kappa0.03$Mdust/beta2.0_kappa0.077$Mdust
)

df = data.frame(
  "beta20" = beta2.0_kappa0.077$Mdust,
  "beta30" = beta3.0_kappa0.077$Mdust,
  "beta15" = beta1.5_kappa0.077$Mdust,
  "beta10" = beta1.0_kappa0.077$Mdust,
  "beta05" = beta0.5_kappa0.077$Mdust, 
  "RA" = RA,
  "Dec" = Dec,
  "z" = redshifts,
  "uberID" = uberID
)
h5file = '~/Documents/DustMassDensity/data/all_data.h5'
h5createGroup(
  file = h5file,
  group = "greybody"
)
h5write(
  file = h5file, 
  obj = df,
  name = "greybody/Mdust"
)

kappa_lambda = (fread("~/Documents/DustMassDensity/data/literature_evo/james2002_kappa_lambda.csv"))
kappa_lambda_func = approxfun(
  log10(kappa_lambda$lambda),
  log10(kappa_lambda$kappa)
)

10^kappa_lambda_func(log10(500)) / 10^kappa_lambda_func(log10(850))



