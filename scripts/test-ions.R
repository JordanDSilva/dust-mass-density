library(Highlander)
library(ProSpect)
library(data.table)
library(magicaxis)
library(celestial)
library(foreach)
library(doParallel)
library(dplyr)
library(Cairo)
library(scales)
library(LaplacesDemon)
library(matrixStats)
library(stringr)
library(rhdf5)

set.seed(666)

speedOfLight = 299792458 ##m/s
planckConst = 6.62607015e-34 ##J/s
photCrossSectionV0 = 6.3e-18 ##cm^-2
freqV0 = 911.75 ## Angstrom
thomsonCrossSection = 6.6524587051e-29 ##m^-2

mProton = 1.67262192e-27 ##kg
mElectron = 9.1093837e-31 ##kg

# recombinationCoeff = 4.2e-13 ## cm^3 s^-1 case A
recombinationCoeff = 2.6e-13 ## cm^3 s^-1
recombinationCoeffMpc = recombinationCoeff * (3.24078e-25)^3 ## Mpc^3 s^-1

km2m = 1000
Mpc2m = 3.08567758e+22
print(Mpc2m/km2m)
Msol2kg = 1.9891e+30

OmegaBaryon = 0.0224 / 0.7^2 ##Planck18

ClumpFactor = 3 ## default value
YHe = 0.24 ## Helium fraction
XH = 1 - YHe ## Hydrogen fraction

AgeOfUniverse = cosdistUniAgeAtz(z = 0, ref = "Planck18")
zgrid = c(seq(0.01, 15.51, 0.25)) ## Mason+19 calculate tau CMB at z=15
UniAge = cosdistUniAgeAtz(z = zgrid, ref = "Planck18")
LumDist = cosdistLumDist(z = zgrid, ref = "Planck18")
cosmo_grid = list(
  "zgrid" = zgrid,
  "UniAge" = UniAge,
  "LumDist" = LumDist
)

compute_ion = function(wavelum) {
  
  integral(
    f = function(x){
      sel = wavelum$wave <= (freqV0 + 1000) ## Only Lyman continuum
      photon_energies = (planckConst * speedOfLight) / (1e-10 * wavelum$wave[sel]) ## Joules
      appx = approxfun(
        x = wavelum$wave[sel], 
        y = (wavelum$lum[sel] * 3.86e26) / (photon_energies) / 1e51, ## Nphotons per second per angstrom
        yleft = 0, yright = 0
      )
      yy = appx(x)
      return(yy)
    },
    xmin = 0, xmax = freqV0, method = "Kronrod", 
  ) * 1e51
}

compute_photo_ion = function(wavelum) {
  
  integral(
    f = function(x){
      sel = wavelum$wave <= (freqV0 + 1000) ## Only Lyman continuum
      photon_energies = (planckConst * speedOfLight) / (1e-10 * wavelum$wave[sel]) ## Joules
      appx = approxfun(
        x = wavelum$wave[sel], 
        y = (wavelum$lum[sel] * 3.86e26) / (photon_energies) / 1e51 / (4 * pi * (8 * 2.561e22)^2) * (freqV0/wavelum$wave[sel])^-3 * 6.3e-18, ## Nphotons per second per angstrom
        yleft = 0, yright = 0
      )
      yy = appx(x)
      return(yy)
    },
    xmin = 0, xmax = freqV0, method = "Kronrod", 
  ) * 1e51
}

## How does viewing angle of AGN manifest as AGN escape fraction in Fritz+06 model
complete_Unobscured = Fritz_interp(
  an = 90, ta = 1e-3, lum = 1e44, 
  
  ## same set up for fitting
  rm = 60,
  be = -0.5,
  al = 4.0,
  ct = 100
)
UnobsNions = compute_ion(wavelum = data.frame("wave"=complete_Unobscured$wave, "lum" = complete_Unobscured$lum/3.826e33))
anGrid = seq(0,90,1.0)
AGNfesc_Fritz = sapply(anGrid, function(an){
  test = Fritz_interp(
    an = an, ta = 1e-3,  lum = 1e44,
    rm = 60,
    be = -0.5,
    al = 4.0,
    ct = 100
  )
  testNions = compute_ion(wavelum = data.frame("wave"=test$wave, "lum" = test$lum/3.826e33))
  testNions / UnobsNions
  
})
FritzFesc2an = approxfun(
  AGNfesc_Fritz, anGrid, rule = 2
)

UV_bandpass = function(lambda){
  ifelse(lambda < 1450 | lambda > 1550, 0, 1)
}

h5file = "~/Documents/JWSTCATALOGUE/Data/all_data.h5"
evol_fits = h5read(
  file = h5file,
  name = "cosmic_evolution/history_fits"
)

cosmic_SEDs = function(cosmo_grid, params, cores, evo = evol_fits, do_quantile = FALSE){
  start = Sys.time()
  
  ## Spit out cosmic SEDs with metallicity Zfinal, for fitting purposes
  zgrid = cosmo_grid$zgrid
  UniAge = cosmo_grid$UniAge
  LumDist = cosmo_grid$LumDist
  
  tau_screen = params$tau_screen
  escape_frac = params$escape_frac
  AGN_escape_frac = params$AGN_escape_frac
  
  if(is.function(escape_frac)){
    escape_frac_Vec = escape_frac(zgrid)
  }else{
    escape_frac_Vec = escape_frac
  }
  
  if(is.function(AGN_escape_frac)){
    AGN_escape_frac_Vec = AGN_escape_frac(zgrid)
  }else{
    AGN_escape_frac_Vec = AGN_escape_frac
  }
  
  tau_birth = params$tau_birth
  tau_AGN = params$tau_AGN
  
  Zstart = params$Zstart
  Zfinal = params$Zfinal
  # Zvec = seq(Zstart, (Zfinal), length.out = length(zgrid))
  
  CSFH_CAGNH_approx_vals = lapply(c(50, 16, 84), function(qq){
    AGNlum = 10^approx(
      x = evo$z,
      y = evo[[paste0("CAGNHQ", qq)]], 
      xout = zgrid
    )$y
    mSFR = 10^approx(
      x = evo$z,
      y = evo[[paste0("CSFHQ", qq)]], 
      xout = zgrid
    )$y
    
    SFHappx_all = approxfun(
      x = 1e9*(AgeOfUniverse - UniAge), ## lookback time
      y = mSFR,
      yleft = 0, yright = 0
    )
    tempZ = Zfunc_massmap_box(
      age = 1e9*(AgeOfUniverse - UniAge),
      Zstart = Zstart,
      Zfinal = Zfinal,
      Zagemax = AgeOfUniverse,
      massfunc = SFHappx_all
    )
    
    data.frame(
      "mSFR" = mSFR,
      "AGNlum" = AGNlum,
      "tempZ" = tempZ
    )
  })
  names(CSFH_CAGNH_approx_vals) = c(paste0("vals_", c(50,16,84)))
  
  registerDoParallel(cores = cores)
  df = foreach(i = 1:length(zgrid), .combine = rbind) %dopar% {
    
    z = zgrid[i]
    message(z)
    
    if(length(escape_frac_Vec)>1){
      escape_frac = escape_frac_Vec[i]
    }
    
    if(length(AGN_escape_frac_Vec)>1){
      AGN_escape_frac = AGN_escape_frac_Vec[i]
    }
    
    agnan = FritzFesc2an(AGN_escape_frac)
    
    dff_func = function(qq){
      
      AGNlum = approx(
        x = zgrid, 
        y = CSFH_CAGNH_approx_vals[[paste0("vals_", qq)]]$AGNlum, 
        xout = z
      )$y
      
      # Truncate the SFH and ZH at each redshift slice
      SFHappx = approxfun(
        x = 1e9*(UniAge[i] - c(UniAge[zgrid>=z],0.09850496)), ## Magemax z = 30, Universe age = 0.0985 Gyr. Need at least 2 values to interpolate in highest z bin / max of zgrid
        y = c(CSFH_CAGNH_approx_vals[[paste0("vals_", qq)]]$mSFR[zgrid >= z],0),
        yleft = 0, yright = 0
      )
      tempZHappx = function(x,...){
        approxfun(
          x = 1e9*(UniAge[i] - c(UniAge[zgrid>=z],0.09850496)), ## Magemax z = 30, Universe age = 0.0985 Gyr. Need at least 2 values to interpolate in highest z bin / max of zgrid
          y = c(CSFH_CAGNH_approx_vals[[paste0("vals_", qq)]]$tempZ[zgrid >= z], Zstart),
          yleft = Zfinal, yright = Zstart)(x)
      }
      
      fiducial_SED = ProSpectSED(
        
        ref = "737",
        H0 = 70,
        OmegaM = 0.3,
        OmegaL = 0.7,
        z = z,
        LumDist_Mpc = LumDist[i],
        agemax = UniAge[i]*1e9,
        speclib = BC03hr,
        
        IGMabsorb = 0.0,
        escape_frac = escape_frac,
        emission = FALSE,
        
        SFH = SFHfunc,
        massfunc = SFHappx,
        
        AGN = Fritz,
        
        ## Same set up for fitting
        AGNrm = 60,
        AGNbe = -0.5,
        AGNal = 4.0,
        AGNct = 100,
        
        Dale = Dale_NormTot,
        Dale_M2L_func = Dale_M2L_func,
        
        tau_screen = tau_screen,
        tau_birth = tau_birth,
        
        AGNlum = ifelse(AGN_escape_frac == 0, 0, AGNlum), 
        AGNta = tau_AGN,
        AGNan = agnan,
        
        Z = tempZHappx,
        waveout = seq(1, 9.35, by = 0.01),
      )
      # par(bg = "grey")
      # magplot(fiducial_SED$FinalLum, xlim = c(10, 1e4), ylim = 10^c(1, 7), log = "xy", type = "l", lwd = 5)
      # lines(fiducial_SED$StarsUnAtten, col = "blue", lty = 2, lwd = 3)
      # lines(fiducial_SED$StarsAtten, col = "darkgreen", lwd = 3)
      # lines(fiducial_SED$AGN, col = "purple", lwd = 3)
      
      names_dff = c(
        paste0("TotSED", qq),
        paste0("StarsSED",qq),
        paste0("AGNSED",qq)
      )
      if(AGN_escape_frac == 0){
        fiducial_SED$AGN = list(
          "wave" = fiducial_SED$StarsAtten$wave,
          "lum" = rep(0, length(fiducial_SED$StarsAtten$wave))
        )
        fiducial_SED$FinalLum = fiducial_SED$StarsAtten
      }
      
      if(escape_frac == 0){
        fiducial_SED$StarsAtten = list(
          "wave" = fiducial_SED$AGN$lum,
          "lum" = rep(0, length(fiducial_SED$AGN$wave))
        )
        fiducial_SED$FinalLum = fiducial_SED$AGN
      }
      yy = list(
        fiducial_SED$FinalLum,
        fiducial_SED$StarsAtten,
        fiducial_SED$AGN
      )
      names(yy) = names_dff
      yy
    }
    
    med_ = dff_func(50)
    if(do_quantile){
      q16_ = dff_func(16)
      q84_ = dff_func(84)
      return( list(med_, q16_, q84_) )
    }else{
      return( list(med_) )
    }
  }
  
  stopImplicitCluster()
  
  TotSEDQ50 = do.call("cbind", lapply(df[,1], function(x){x$TotSED50$lum}))
  StarsSEDQ50 = do.call("cbind", lapply(df[,1], function(x){x$StarsSED50$lum}))
  AGNSEDQ50 = do.call("cbind", lapply(df[,1], function(x){x$AGNSED50$lum}))
  
  TotWavelength = df[[1]]$TotSED50$wave  
  StarsWavelength = df[[1]]$StarsSED50$wave
  AGNWavelength = df[[1]]$AGNSED50$wave
  
  ## Lyc 
  tot_nion = foreach(i = 1:length(zgrid), .combine = 'c') %do% {
    compute_ion(
      wavelum = list(
        "wave" = TotWavelength,
        "lum" = TotSEDQ50[,i]
      )
    )
  }
  stars_nion = foreach(i = 1:length(zgrid), .combine = 'c') %do% {
    compute_ion(
      wavelum = list(
        "wave" = StarsWavelength,
        "lum" = StarsSEDQ50[,i]
      )
    )
  }
  agn_nion = foreach(i = 1:length(zgrid), .combine = 'c') %do% {
    compute_ion(
      wavelum = list(
        "wave" = AGNWavelength,
        "lum" = AGNSEDQ50[,i]
      )
    )
  }
  
  ## UV at 1500Ang
  tot_UV = foreach(i = 1:length(zgrid), .combine = 'c') %do% {
    photom_lum(wave = TotWavelength, 
               lum = TotSEDQ50[,i], 
               outtype = "CGS", 
               filters = list("UV1500"=UV_bandpass), 
               z = 0, ref = "737") *  4*pi*(10 * 3.086e18)^2
  }
  stars_UV = foreach(i = 1:length(zgrid), .combine = 'c') %do% {
    photom_lum(wave = StarsWavelength, 
               lum = StarsSEDQ50[,i], 
               outtype = "CGS", 
               filters = list("UV1500"=UV_bandpass), 
               z = 0, ref = "737") *  4*pi*(10 * 3.086e18)^2
  }
  agn_UV = foreach(i = 1:length(zgrid), .combine = 'c') %do% {
    photom_lum(wave = AGNWavelength, 
               lum = AGNSEDQ50[,i], 
               outtype = "CGS", 
               filters = list("UV1500"=UV_bandpass), 
               z = 0, ref = "737") *  4*pi*(10 * 3.086e18)^2
  }
  return_ = list(
    "TotWavelength" =TotWavelength,
    "StarsWavelength" = StarsWavelength,
    "AGNWavelength" = AGNWavelength,
    
    "TotSEDQ50" = TotSEDQ50,
    "AGNSEDQ50" = AGNSEDQ50,
    "StarsSEDQ50" = StarsSEDQ50,
    
    "NionsTotalQ50" = tot_nion,
    "NionsStellarQ50" = stars_nion,
    "NionsAGNQ50" = agn_nion,
    
    "UVTotalQ50" = tot_UV,
    "UVStellarQ50" = stars_UV,
    "UVAGNQ50" = agn_UV
  )
  
  if(do_quantile){
    TotSEDQ16 = do.call("cbind", lapply(df[,2], function(x){x$TotSED16$lum}))
    StarsSEDQ16 = do.call("cbind", lapply(df[,2], function(x){x$StarsSED16$lum}))
    AGNSEDQ16 = do.call("cbind", lapply(df[,2], function(x){x$AGNSED16$lum}))
    TotSEDQ84 = do.call("cbind", lapply(df[,3], function(x){x$TotSED84$lum}))
    StarsSEDQ84 = do.call("cbind", lapply(df[,3], function(x){x$StarsSED84$lum}))
    AGNSEDQ84 = do.call("cbind", lapply(df[,3], function(x){x$AGNSED84$lum}))
    
    tot_nionQ16 = foreach(i = 1:length(zgrid), .combine = 'c') %do% {
      compute_ion(
        wavelum = list(
          "wave" = TotWavelength,
          "lum" = TotSEDQ16[,i]
        )
      )
    }
    stars_nionQ16 = foreach(i = 1:length(zgrid), .combine = 'c') %do% {
      compute_ion(
        wavelum = list(
          "wave" = StarsWavelength,
          "lum" = StarsSEDQ16[,i]
        )
      )
    }
    agn_nionQ16 = foreach(i = 1:length(zgrid), .combine = 'c') %do% {
      compute_ion(
        wavelum = list(
          "wave" = AGNWavelength,
          "lum" = AGNSEDQ16[,i]
        )
      )
    }
    
    tot_nionQ84 = foreach(i = 1:length(zgrid), .combine = 'c') %do% {
      compute_ion(
        wavelum = list(
          "wave" = TotWavelength,
          "lum" = TotSEDQ84[,i]
        )
      )
    }
    stars_nionQ84 = foreach(i = 1:length(zgrid), .combine = 'c') %do% {
      compute_ion(
        wavelum = list(
          "wave" = StarsWavelength,
          "lum" = StarsSEDQ84[,i]
        )
      )
    }
    agn_nionQ84 = foreach(i = 1:length(zgrid), .combine = 'c') %do% {
      compute_ion(
        wavelum = list(
          "wave" = AGNWavelength,
          "lum" = AGNSEDQ84[,i]
        )
      )
    }
    
    tot_UV16 = foreach(i = 1:length(zgrid), .combine = 'c') %do% {
      photom_lum(wave = TotWavelength, 
                 lum = TotSEDQ16[,i], 
                 outtype = "CGS", 
                 filters = list("UV1500"=UV_bandpass), 
                 z = 0, ref = "737") *  4*pi*(10 * 3.086e18)^2
    }
    stars_UV16 = foreach(i = 1:length(zgrid), .combine = 'c') %do% {
      photom_lum(wave = StarsWavelength, 
                 lum = StarsSEDQ16[,i], 
                 outtype = "CGS", 
                 filters = list("UV1500"=UV_bandpass), 
                 z = 0, ref = "737") *  4*pi*(10 * 3.086e18)^2
    }
    agn_UV16 = foreach(i = 1:length(zgrid), .combine = 'c') %do% {
      photom_lum(wave = AGNWavelength, 
                 lum = AGNSEDQ16[,i], 
                 outtype = "CGS", 
                 filters = list("UV1500"=UV_bandpass), 
                 z = 0, ref = "737") *  4*pi*(10 * 3.086e18)^2
    }
    
    tot_UV84 = foreach(i = 1:length(zgrid), .combine = 'c') %do% {
      photom_lum(wave = TotWavelength, 
                 lum = TotSEDQ84[,i], 
                 outtype = "CGS", 
                 filters = list("UV1500"=UV_bandpass), 
                 z = 0, ref = "737") *  4*pi*(10 * 3.086e18)^2
    }
    stars_UV84 = foreach(i = 1:length(zgrid), .combine = 'c') %do% {
      photom_lum(wave = StarsWavelength, 
                 lum = StarsSEDQ84[,i], 
                 outtype = "CGS", 
                 filters = list("UV1500"=UV_bandpass), 
                 z = 0, ref = "737") *  4*pi*(10 * 3.086e18)^2
    }
    agn_UV84 = foreach(i = 1:length(zgrid), .combine = 'c') %do% {
      photom_lum(wave = AGNWavelength, 
                 lum = AGNSEDQ84[,i], 
                 outtype = "CGS", 
                 filters = list("UV1500"=UV_bandpass), 
                 z = 0, ref = "737") *  4*pi*(10 * 3.086e18)^2
    }
    
    return_$TotSEDQ16 = TotSEDQ16
    return_$StarsSEDQ16 = StarsSEDQ16
    return_$AGNSEDQ16 = AGNSEDQ16
    
    return_$TotSED84 = TotSEDQ84
    return_$StarsSEDQ84 = StarsSEDQ84
    return_$AGNSEDQ84 = AGNSEDQ84
    
    return_$NionsTotalQ16 = tot_nionQ16
    return_$NionsStellarQ16 = stars_nionQ16
    return_$NionsAGNQ16 = agn_nionQ16
    
    return_$NionsTotalQ84 = tot_nionQ84
    return_$NionsStellarQ84 = stars_nionQ84
    return_$NionsAGNQ84 = agn_nionQ84
    
    return_$UVTotalQ16 = tot_UV16
    return_$UVStarsQ16 = stars_UV16
    return_$UVAGNQ16 = agn_UV16
    
    return_$UVTotalQ84 = tot_UV84
    return_$UVStarsQ84 = stars_UV84
    return_$UVAGNQ84 = agn_UV84
  }
  
  print(
    Sys.time() - start
  )
  
  return( 
    return_
  ) 
}

## Compute maximum and minimum contributions
AverageLyc = cosmic_SEDs(
  cosmo_grid = cosmo_grid,
  params=list(
    "tau_screen" = 1e-3,
    "tau_AGN" = 1e-3,
    "tau_birth" = 1e-3,
    "AGN_escape_frac" = 1.0,
    "escape_frac" = 1.0,
    'Zstart' = 1e-4,
    'Zfinal' = 0.02
  ), 
  evo = evol_fits,
  cores = 1, 
  do_quantile = TRUE
)

## Mass in the Stromgren sphere
n = 100
Q = AverageLyc$NionsTotalQ50

Rs = ((3 * Q) / (4 * pi * recombinationCoeff * n^2))^(1/3)

Q * mProton/Msol2kg / (recombinationCoeff * n)

magplot(
  zgrid, 
  ( (4 * pi / 3) * Rs^3 * mProton * n / Msol2kg ), 
  log = "y",
  xlim = c(0, 4),
  ylim = c(1e3, 1e6),
  xlab = "Redshift", 
  ylab = "Mass of Stromgren sphere [Msun Mpc^-3]",
  lty = 1, 
  pch = 16
)

df = data.frame(
  "z" = zgrid, 
  "HIIQ50" = AverageLyc$NionsTotalQ50 * mProton/Msol2kg / (recombinationCoeff * n),
  "HIIQ16" = AverageLyc$NionsTotalQ16 * mProton/Msol2kg / (recombinationCoeff * n),
  "HIIQ84" = AverageLyc$NionsTotalQ84 * mProton/Msol2kg / (recombinationCoeff * n)
)
fwrite(
  df,
  "~/Documents/DustMassDensity/data/literature_evo/HIIDensity.csv"
)
