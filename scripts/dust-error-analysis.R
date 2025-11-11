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

load("~/Documents/DustMassDensity/data/compare_prospect_magphys.Rdata") ## From compare_prospect_magphys.R

## df_spec Stack of 218 FIR S/N > 5 confirmed galaxies

ff = unlist(unname(df_spec$ProFoundStack))
fferr = sqrt( unlist(unname(df_spec$ProFoundErr))^2 + (0.1 * ff)^2 )
wavelength = df_spec$wavelength

ff[ff == -999] = NA

FIR_cutoff_wavelength = 224000 ## W4, to be safe (next one is P100 and the FIR Herscel bands)

do_fit = function(cutoff_wave, pfunc){
  
  flux_data = list()
  flux_data$filter = c("FUV", "NUV", "u_VST", "g_VST" ,"r_VST" ,"i_VST" ,"Z_VISTA" ,"Y_VISTA" ,"J_VISTA" ,"H_VISTA" ,"K_VISTA" ,"W1_WISE" ,"W2_WISE" ,"W3_WISE" ,"W4_WISE" ,"P100_Herschel" ,"P160_Herschel" ,"S250_Herschel" ,"S350_Herschel" ,"S500_Herschel")
  flux_data$cenwave = c(1539, 2316, 3528, 4760, 6326, 7599, 6654, 10229, 12556, 16499, 21571, 3.4e4, 4.65e4, 12.8e4, 22.4e4, 98.9e4, 156e4, 249e4, 350e4, 504e4)
  flux_data$flux = ff
  flux_data$fluxerr = fferr
  flux_data = data.frame(flux_data)
  
  filtout = list(
    approxfun(getfilt("FUV")),
    approxfun(getfilt("NUV")),
    approxfun(getfilt("u_VST")),
    approxfun(getfilt("g_VST")),
    approxfun(getfilt("r_VST")),
    approxfun(getfilt("i_VST")),
    approxfun(getfilt("Z_VISTA")),
    approxfun(getfilt("Y_VISTA")),
    approxfun(getfilt("J_VISTA")),
    approxfun(getfilt("H_VISTA")),
    approxfun(getfilt("K_VISTA")),
    approxfun(getfilt("W1_WISE")),
    approxfun(getfilt("W2_WISE")),
    approxfun(getfilt("W3_WISE")),
    approxfun(getfilt("W4_WISE")),
    approxfun(getfilt("P100_Herschel")),
    approxfun(getfilt("P160_Herschel")),
    approxfun(getfilt("S250_Herschel")),
    approxfun(getfilt("S350_Herschel")),
    approxfun(getfilt("S500_Herschel"))
  )
  
  filtout = filtout[flux_data$cenwave < cutoff_wave]
  flux_data = flux_data[flux_data$cenwave < cutoff_wave, ]

  
  ## Take mean redshift of sample 
  z = median(df$zProSpect)
  LookbackTime = cosdistTravelTime(z, ref = 'Planck15')*1e9
  agemax = 13.38e9 - LookbackTime
  upperlimit = (agemax/1e9) 
  lowerlimit = (0 - LookbackTime)/1e9
  
  magemax = agemax/1e9
  Zagemax = agemax/1e9
  LumDist_Mpc = cosdistLumDist(z = z, ref = "Planck15")
  
  pro_data=list(
    flux=flux_data,
    arglist=list(
      z = z, 
      agemax = agemax, 
      LumDist_Mpc = LumDist_Mpc,
      
      magemax = magemax, 
      Zagemax = Zagemax, 
      
      massfunc=massfunc_snorm_trunc,
      Z=Zfunc_massmap_lin,
      Zstart = 1e-4,
      
      ref="Planck15",
      H0 = 67.8,
      OmegaM = 0.308,
      
      SMstar = TRUE
    ),
    filtout=filtout,
    
    speclib=BC03lr,
    SFH=SFHfunc,
    
    Dale=Dale_NormTot,
    Dale_M2L_func = Dale_M2L_func,
    
    AGN = NULL,
    
    parm.names=c('mSFR','mpeak','mperiod','mskew','tau_birth','tau_screen','alpha_SF_birth','alpha_SF_screen','Zfinal'), 
    mon.names=c("LP","masstot","dustmass.birth", "dustmass.screen", "dustmass.total", "dustlum.birth", "dustlum.screen", 
                "dustlum.total", "SFRburst", paste("flux.",flux_data$filter,sep='')), 
    logged=c(T,F,T,F,T,T,F,F,T),
    intervals=list(lo=c(-3,-2,log10(0.3),-0.5,-2.5,-5,0,0,-4), 
                   hi=c(4,agemax/1e9,2,1,1.5,1,4,4,-1.3)), 
    fit = 'LD',
    N=dim(flux_data)[1],
    prior=pfunc,
    verbose=FALSE
  )
  
  pro_data$flux$cenwave = flux_data$cenwave
  
  startpoint = (pro_data$intervals$lo+pro_data$intervals$hi)/2
  names(startpoint) = pro_data$parm.names
  
  highout = Highlander(
    startpoint, pro_data, ProSpectSEDlike, 
    Niters=c(2000, 2000),  NfinalMCMC = 2000, 
    lower=pro_data$intervals$lo, upper=pro_data$intervals$hi,
    seed=666, optim_iters = 2, likefunctype = "LD"
  )
  pro_data$fit = 'check'
  bestfit = ProSpectSEDlike(highout$par, Data=pro_data)
  
  return(
    list(
      "bestfit" = bestfit, 
      "highout" = highout
    )
  )
}

prior_func = function(parm){
  sum(
    c(
      100 * erf(parm['mperiod'] + 2) - 100,
      dnorm(parm['tau_birth'],mean=0.2,sd=0.5,log=TRUE),
      -20 * erf(parm['tau_screen'] - 2),
      dnorm(parm['alpha_SF_birth'],mean=2,sd=1,log=TRUE),
      dnorm(parm['alpha_SF_screen'],mean=2,sd=1,log=TRUE)
    )
  )
}
reduced_prior_func = function(parm){
  sum(
    c(
      100 * erf(parm['mperiod'] + 2) - 100,
      dnorm(parm['tau_birth'],mean=0.2,sd=0.5,log=TRUE),
      -20 * erf(parm['tau_screen'] - 2)
      # dnorm(parm['alpha_SF_birth'],mean=2,sd=1,log=TRUE),
      # dnorm(parm['alpha_SF_screen'],mean=2,sd=1,log=TRUE)
    )
  )
}
default_fit = do_fit(cutoff_wave = Inf, pfunc = prior_func)
fir_no_prior_fit = do_fit(cutoff_wave = Inf, pfunc = reduced_prior_func)
better_fit = do_fit(cutoff_wave = FIR_cutoff_wavelength, pfunc = prior_func)
worse_fit = do_fit(cutoff_wave = FIR_cutoff_wavelength, pfunc = reduced_prior_func)

plot(
  density(
    log10(worse_fit$highout$LD_last$Monitor[,'dustmass.total'])
  ), xlim = c(4, 9), ylim = c(1e-3, 7), xlab = "log Dust mass", ylab = "Posterior", lwd = 3, 
  col = "cornflowerblue"
)
lines(
  density(
    log10(better_fit$highout$LD_last$Monitor[,'dustmass.total'])
  ), 
  col = "orange", lwd = 3
)
lines(
  density(
    log10(fir_no_prior_fit$highout$LD_last$Monitor[,'dustmass.total'])
  ), 
  col = "forestgreen", lwd = 3
)
lines(
  density(
    log10(default_fit$highout$LD_last$Monitor[,'dustmass.total'])
  ), 
  col = "black", lwd = 3
)
lines(
  density(
    log10(df$MdustProSpect), weights = (df$MdustErrProSpect/(log(10)*df$MdustProSpect))^-2
  ), col = "red"
)
lines(
  density(
    log10(df$MdustMagphys)
  ), col = "red"
)


weighted.mean(
  df$MdustProSpect, w = df$MdustErrProSpect^-2
)
1/sqrt(sum(df$MdustErrProSpect^-2))

plot(
  default_fit$bestfit
)
plot(
  worse_fit$bestfit
)

magplot(
  default_fit$bestfit$SEDout$DustEmit, log = "xy", type = "l", xlim = c(1e5, 1e7)
)
lines(
  worse_fit$bestfit$SEDout$DustEmit
)
legend(
  x = "topleft", 
  legend = c(
    "FIR+Prior",
    "FIR+no Prior",
    "no FIR+Prior", 
    "no Fir+ no Prior"
  ),
  col = c(
    "black",
    "forestgreen",
    "orange", 
    "cornflowerblue"
  ),
  lwd = 3
)

lambdaMax_default = default_fit$bestfit$SEDout$DustEmit$wave[default_fit$bestfit$SEDout$DustEmit$wave>2e5][which.max(default_fit$bestfit$SEDout$DustEmit$lum[default_fit$bestfit$SEDout$DustEmit$wave>2e5])]
lambdaMax_worse = worse_fit$bestfit$SEDout$DustEmit$wave[worse_fit$bestfit$SEDout$DustEmit$wave>2e5][which.max(worse_fit$bestfit$SEDout$DustEmit$lum[worse_fit$bestfit$SEDout$DustEmit$wave>2e5])]
abline(v = lambdaMax_default)
abline(v = lambdaMax_worse)

2.9e-3/(lambdaMax_default * 1e-10)
2.9e-3/(lambdaMax_worse * 1e-10)

median(
  log10(default_fit$highout$LD_last$Monitor[,'dustmass.total'])
) - median(
  log10(worse_fit$highout$LD_last$Monitor[,'dustmass.total'])
)


df_posteriors = data.frame(
  "FirPrior" = log10(default_fit$highout$LD_last$Monitor[,'dustmass.total']),
  "FirnoPrior" = log10(fir_no_prior_fit$highout$LD_last$Monitor[,'dustmass.total']), 
  "noFiIRPrior" = log10(better_fit$highout$LD_last$Monitor[,'dustmass.total']),
  "noFIRnoPrior" = log10(worse_fit$highout$LD_last$Monitor[,'dustmass.total'])
)




median(
  log10(worse_fit$highout$LD_last$Monitor[,'dustmass.total'])
)
median(
  log10(default_fit$highout$LD_last$Monitor[,'dustmass.total'])
)

legend(
  x = "topleft", 
  legend = c("No FIR + No prior", "No FIR + Prior", "FIR+Prior"), 
  lwd = 3, 
  col = c("cornflowerblue", "orange", "black")
)


catalogueDir = "/Users/22252335/Documents/GAMA-DEVILS-SFR-AGN/data/"
devilsd10_AGN = fread("~/Documents/GAMA-DEVILS-SFR-AGN/data/AGNTotalCat_MasterCat4.csv")
devilsd10_noAGN = readRDS(paste0(catalogueDir, 'DEVILS_D10ProSpectCat_02_02_2021_v0.3.rds')) #catalogue that Jess done without AGN contribution
devilsd10_noAGN = devilsd10_noAGN$cat
devilsd10_noAGN$area = 1.47
devilsd10_AGN$area = 1.47
devilsd10_noAGN = data.frame(devilsd10_noAGN[order(devilsd10_noAGN$UID),])
devilsd10_AGN = data.frame(devilsd10_AGN[order(devilsd10_AGN$UID),])
devilsd10_AGN$FIRInput = ifelse( is.na(devilsd10_AGN$FIRInput), 0L, devilsd10_AGN$FIRInput)
devilsd10_noAGN$FIRInput = devilsd10_AGN$FIRInput


## make a hybrid catalogue
devilsd10_col_names = c(
  "UID", "z", 
  "StellarMass", "dustmass.birth","dustmass.screen","dustmass.total","dustlum.birth","dustlum.screen","dustlum.total","Zfinal",
  "StellarMass_LB", "dustmass.birth_LB","dustmass.screen_LB","dustmass.total_LB","dustlum.birth_LB","dustlum.screen_LB","dustlum.total_LB","Zfinal_LB",
  "StellarMass_UB", "dustmass.birth_UB","dustmass.screen_UB","dustmass.total_UB","dustlum.birth_UB","dustlum.screen_UB","dustlum.total_UB","Zfinal_UB",
  "FIRInput"
)
devilsd10_hybrid = devilsd10_noAGN[, devilsd10_col_names]
devilsd10_AGN_idx = devilsd10_AGN$AGNfrac >= 0.1 & devilsd10_AGN$LP > devilsd10_noAGN$LP
message("AGN preferred fraction: ", sum(devilsd10_AGN_idx)/dim(devilsd10_hybrid)[1])
devilsd10_hybrid[devilsd10_AGN_idx, devilsd10_col_names] = devilsd10_AGN[devilsd10_AGN_idx, devilsd10_col_names]

median(devilsd10_hybrid$z[devilsd10_hybrid$FIRInput == 1])
zvec = seq(0, 30, 0.01)
lbtvec = cosdistTravelTime(z = zvec, ref = "Planck18")
lbt2z = approxfun(
  lbtvec, zvec
)

lbt_bins = seq(0, 12, 0.75)
lbt_mids = lbt_bins[-length(lbt_bins)] + diff(lbt_bins) / 2

zbins = lbt2z(lbt_bins)
zmids = zbins[-length(zbins)] + diff(zbins) / 2

devils_z_hist = foreach(ii = 1:(length(zbins)-1), .combine = c) %do% {
  zidx = devilsd10_hybrid$z >= zbins[ii] & devilsd10_hybrid$z < zbins[ii + 1]
  sum(zidx)
}
FIR_input_hist = foreach(ii = 1:(length(zbins)-1), .combine = c) %do% {
  zidx = devilsd10_hybrid$z >= zbins[ii] & devilsd10_hybrid$z < zbins[ii + 1]
  sum(devilsd10_hybrid$FIRInput[zidx])
}
magplot(
  zmids, FIR_input_hist/devils_z_hist, pch = 16, type = "s", log = "y"
)

devils_FIR_stats = data.frame(
  "mids" = zmids,
  "stats" = FIR_input_hist/devils_z_hist
)

h5file = '~/Documents/DustMassDensity/data/all_data.h5'
# h5delete(h5file, "cosmic")
h5createGroup(h5file, "FIRDustMassErrs")
h5write(obj = devils_FIR_stats, file = h5file, name = "FIRDustMassErrs/FIRStats")
h5write(obj = df_posteriors, file = h5file, name = "FIRDustMassErrs/DustMassPosteriors")
