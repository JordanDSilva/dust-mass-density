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
library(Rfits)
library(fields)

# 06 March 2026 11:36AM, Add the variable M2L thing for the _corr results. 

h5file = '~/Documents/DustMassDensity/data/all_data.h5'
catalogueDir = "~/Documents/DustMassDensity/data/gama_devils_catalogues/"

gama_noAGN = data.frame(Rfits_read_table("~/Documents/DustMassDensity/data/GAMA/ProSpectv03.fits"))
gama_noAGN$area = 217.54
gama_AGN = data.frame(Rfits_read_table("~/Documents/DustMassDensity/data/GAMA/ProSpectAGNv02.fits"))
gama_AGN$area = 217.54
gama_science = data.frame(Rfits_read_table("~/Documents/DustMassDensity/data/GAMA/gkvScienceCatv02.fits"))
gama_trim = gama_science[gama_science$SC >= 4 & gama_science$NQ > 2,]
gama_noAGN$z = gama_trim$Z
gama_noAGN_match = gama_noAGN[gama_noAGN$CATAID %in% gama_AGN$CATAID, ]
gama_AGN$z = gama_noAGN_match$z
gama_noAGN = gama_noAGN_match
rm(gama_science)
rm(gama_trim)
rm(gama_noAGN_match)

devilsd10_AGN = fread("~/Documents/DustMassDensity/data/gama_devils_catalogues/AGNTotalCat_MasterCat4.csv")
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
  "StellarMass", "dustmass.birth","dustmass.screen","dustmass.total","dustlum.birth","dustlum.screen","dustlum.total","Zfinal", "alpha_SF_birth", "alpha_SF_screen",
  "StellarMass_LB", "dustmass.birth_LB","dustmass.screen_LB","dustmass.total_LB","dustlum.birth_LB","dustlum.screen_LB","dustlum.total_LB","Zfinal_LB", "alpha_SF_birth_LB", "alpha_SF_screen_LB",
  "StellarMass_UB", "dustmass.birth_UB","dustmass.screen_UB","dustmass.total_UB","dustlum.birth_UB","dustlum.screen_UB","dustlum.total_UB","Zfinal_UB", "alpha_SF_birth_UB", "alpha_SF_screen_UB",
  "FIRInput"
)
devilsd10_hybrid = devilsd10_noAGN[, devilsd10_col_names]
devilsd10_AGN_idx = devilsd10_AGN$AGNfrac >= 0.1 & devilsd10_AGN$LP > devilsd10_noAGN$LP
message("AGN preferred fraction: ", sum(devilsd10_AGN_idx)/dim(devilsd10_hybrid)[1])
devilsd10_hybrid[devilsd10_AGN_idx, devilsd10_col_names] = devilsd10_AGN[devilsd10_AGN_idx, devilsd10_col_names]

gama_col_names = c(
  "uberID", "CATAID", "z",
  "StellarMass_50", "StellarMass_16", "StellarMass_84",  
  "DustMass_50", "DustMass_16", "DustMass_84",
  "DustLum_50", "DustLum_16", "DustLum_84",
  "Zgas_50", "Zgas_16", "Zgas_84"
)
gama_hybrid = gama_noAGN[,gama_col_names]
gama_AGN_idx = gama_AGN$fAGN_bestfit >= 0.1 & gama_AGN$LP > (gama_noAGN$LP)
# gama_AGN_idx = gama_AGN$fAGN_bestfit >= 0.1 & -2*gama_AGN$LP+12*log(24) <-2*gama_noAGN$LP+9*log(24)
message("AGN preferred fraction: ", sum(gama_AGN_idx)/dim(gama_hybrid)[1])
gama_hybrid[gama_AGN_idx, gama_col_names] = gama_noAGN[gama_AGN_idx, gama_col_names]

RR14_BPL = function(Z, doDTG = FALSE){
  
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
  alphaL = 3.10
  xt = 8.10
  
  # a = par[1]
  # alphaH = par[2]
  # b = par[3]
  # alphaL = par[1]
  # xt = par[2]
  
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
  ## DTH = Mdust / Mhydrogen = Mdust / (Mgas / muGal) = DTG * muGal = DTG * (1 / (1 - Ysol - Zgal))
  muGal = 1 / (1 - 0.270 - Z)
  # muGal = (1-Z) / (1 + (4/3))
  DTH = DTG*muGal
  if(doDTG){
    return(DTG)
  }else{
    return(DTH)
  }
}
saveRDS(RR14_BPL, "~/Documents/DustMassDensity/data/RR14_BPL.rds")
deVis19_PL = function(Z, doDTG = FALSE){
  #2.45 * (12+log(O/H)) - 23.30
  
  xSol = 8.69
  ZOH = log10(Z / 0.014) + xSol
  
  a = 2.45
  b = -23.30
  
  DTG = 10^(a * ZOH + b)
  
  fz = (Z/0.014) * 0.014
  xiGal = (1 - (0.2485 + fz*1.41) - fz)^-1
  # muGal = (1-Z) / (1 + (4/3))
  DTH = DTG*xiGal
  if(doDTG){
    return(DTG)
  }else{
    return(DTH)
  }
}

## Define qPAH early on 
qPAH_VSG = 0.035
message("Using qPAH_VSG = ", qPAH_VSG)

shivaei24 = fread("~/Documents/DustMassDensity/data/shivaeqPAHZ.csv")
shivaei24_fit = approxfun(
  shivaei24$OH,
  shivaei24$q, 
  rule = 2
)
shivaei24_qPAHZ = function(Z){
  xSol = 8.69
  ZOH = log10(Z / 0.014) + xSol
  
  ff = c(shivaei24_fit(ZOH)) * 0.01
  ff[ZOH >= 8.4] =  0.035
  ff[ZOH < 8.1] =  0.01
  return(ff)
}
metallicity_vec = 10^seq(-4, log10(0.06), 0.001)
qPAHZcorr = sapply(metallicity_vec, function(z){Dale_M2L_variableDTH_func(Dale_M2L$alpha_SF, qPAH_VSG = shivaei24_qPAHZ(z)) / Dale_M2L_variableDTH_func(Dale_M2L$alpha_SF, qPAH_VSG = qPAH_VSG)})

qPAHZ_corr_func = function(Z, alpha) {
  alpha_clamped = pmax(min(Dale_M2L$alpha_SF), pmin(max(Dale_M2L$alpha_SF), alpha))
  Z_clamped     = pmax(min(metallicity_vec),     pmin(max(metallicity_vec),     Z))
  
  interp.surface(
    list(x = Dale_M2L$alpha_SF, y = metallicity_vec, z = qPAHZcorr),
    cbind(alpha_clamped,Z_clamped)
  )
}
save(shivaei24, shivaei24_fit, shivaei24_qPAHZ, metallicity_vec, qPAHZcorr, qPAHZ_corr_func, file = "~/Documents/DustMassDensity/data/qPAHZ.rda")

# Laod in constant correction factor
Md_corr = as.numeric(h5read(
  h5file,
  "Md_corr"
))

# Invert Dale_M2L
Dale_L2M = approxfun(
  x = Dale_M2L_func(seq(0,4,0.01)),
  y = seq(0,4,0.01), 
  rule = 2
)

zvec = seq(0, 30, 0.01)
lbtvec = cosdistTravelTime(z = zvec, ref = "Planck18")
lbt2z = approxfun(
  lbtvec, zvec
)

lbt_bins = seq(0, 12, 0.75)
lbt_mids = lbt_bins[-length(lbt_bins)] + diff(lbt_bins) / 2

zbins = lbt2z(lbt_bins)
zmids = zbins[-length(zbins)] + diff(zbins) / 2
vol_mids = 4*pi/3 * (cosdistCoDist(z = zmids + diff(zbins) / 2, ref = 'Planck18')^3 - cosdistCoDist(z = zmids - diff(zbins) / 2, ref = 'Planck18')^3)

sm_bins = seq(3.6, 15.6, 0.6)
sm_mids = sm_bins[-length(sm_bins)] + diff(sm_bins) / 2

mdustvec = seq(3.5, 15.5, 0.1)
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

## Get the monte carlo error on the DMF
Nsamples_MCERR = 1000

mc_err_gama_noAGN_dmf = foreach(i = 1:length(zmids)) %do% {
  message(i)
  vol = vol_mids[i] * 217.54 / (4*pi*(180/pi)^2)
  zidx = gama_noAGN$z >= zbins[i] & gama_noAGN$z < zbins[i+1]
  Mdust = gama_noAGN$DustMass_50[zidx]
  MdustErr = 0.5 * (gama_noAGN$DustMass_84[zidx] - gama_noAGN$DustMass_16[zidx])
  
  samples = matrix(rnorm(n = Nsamples_MCERR*sum(zidx), mean = Mdust, sd = MdustErr), nrow = Nsamples_MCERR, ncol = sum(zidx))
  
  hist_samples = foreach(j = 1:Nsamples_MCERR, .combine = "rbind") %do% {
    m_sample = samples[j,]
    log_m_sample = log10(m_sample[m_sample > 0])
    hh = maghist(x = log_m_sample, breaks = sm_bins, plot = FALSE, verbose = FALSE)
    bin_dmf = hh$counts
    return(bin_dmf)
  }
  hist_quant = colQuantiles(hist_samples, probs = c(0.5, 0.16, 0.84))
}
mc_err_gama_AGN_dmf = foreach(i = 1:length(zmids)) %do% {
  message(i)
  vol = vol_mids[i] * 217.54 / (4*pi*(180/pi)^2)
  zidx = gama_AGN$z >= zbins[i] & gama_AGN$z < zbins[i+1]
  Mdust = gama_AGN$DustMass_50[zidx]
  MdustErr = 0.5 * (gama_AGN$DustMass_84[zidx] - gama_AGN$DustMass_16[zidx])
  
  samples = matrix(rnorm(n = Nsamples_MCERR*sum(zidx), mean = Mdust, sd = MdustErr), nrow = Nsamples_MCERR, ncol = sum(zidx))
  
  hist_samples = foreach(j = 1:Nsamples_MCERR, .combine = "rbind") %do% {
    m_sample = samples[j,]
    log_m_sample = log10(m_sample[m_sample > 0])
    hh = maghist(x = log_m_sample, breaks = sm_bins, plot = FALSE, verbose = FALSE)
    bin_dmf = hh$counts
    return(bin_dmf)
  }
  hist_quant = colQuantiles(hist_samples, probs = c(0.5, 0.16, 0.84))
}
mc_err_gama_hybrid_dmf = foreach(i = 1:length(zmids)) %do% {
  message(i)
  vol = vol_mids[i] * 217.54 / (4*pi*(180/pi)^2)
  zidx = gama_hybrid$z >= zbins[i] & gama_hybrid$z < zbins[i+1]
  
  Mdust = gama_hybrid$DustMass_50[zidx]
  MdustErr = 0.5 * (gama_hybrid$DustMass_84[zidx] - gama_hybrid$DustMass_16[zidx])
  
  samples = matrix(rnorm(n = Nsamples_MCERR*sum(zidx), mean = Mdust, sd = MdustErr), nrow = Nsamples_MCERR, ncol = sum(zidx))
  
  hist_samples = foreach(j = 1:Nsamples_MCERR, .combine = "rbind") %do% {
    m_sample = samples[j,]
    log_m_sample = log10(m_sample[m_sample > 0])
    hh = maghist(x = log_m_sample, breaks = sm_bins, plot = FALSE, verbose = FALSE)
    bin_dmf = hh$counts
    return(bin_dmf)
  }
  hist_quant = colQuantiles(hist_samples, probs = c(0.5, 0.16, 0.84))
}
mc_err_devilsd10_noAGN_dmf = foreach(i = 1:length(zmids)) %do% {
  message(i)
  vol = vol_mids[i] * 1.47 / (4*pi*(180/pi)^2)
  zidx = devilsd10_noAGN$z >= zbins[i] & devilsd10_noAGN$z < zbins[i+1]
  Mdust = devilsd10_noAGN$dustmass.total[zidx]
  MdustErr = 0.5 * (devilsd10_noAGN$dustmass.total_UB[zidx] - devilsd10_noAGN$dustmass.total_LB[zidx])
  
  samples = matrix(rnorm(n = Nsamples_MCERR*sum(zidx), mean = Mdust, sd = MdustErr), nrow = Nsamples_MCERR, ncol = sum(zidx))
  
  hist_samples = foreach(j = 1:Nsamples_MCERR, .combine = "rbind") %do% {
    m_sample = samples[j,]
    log_m_sample = log10(m_sample[m_sample > 0])
    hh = maghist(x = log_m_sample, breaks = sm_bins, plot = FALSE, verbose = FALSE)
    bin_dmf = hh$counts
    return(bin_dmf)
  }
  hist_quant = colQuantiles(hist_samples, probs = c(0.5, 0.16, 0.84))
}
mc_err_devilsd10_AGN_dmf = foreach(i = 1:length(zmids)) %do% {
  message(i)
  vol = vol_mids[i] * 1.47 / (4*pi*(180/pi)^2)
  zidx = devilsd10_AGN$z >= zbins[i] & devilsd10_AGN$z < zbins[i+1]
  Mdust = devilsd10_AGN$dustmass.total[zidx]
  MdustErr = 0.5 * (devilsd10_AGN$dustmass.total_UB[zidx] - devilsd10_AGN$dustmass.total_LB[zidx])
  
  samples = matrix(rnorm(n = Nsamples_MCERR*sum(zidx), mean = Mdust, sd = MdustErr), nrow = Nsamples_MCERR, ncol = sum(zidx))
  
  hist_samples = foreach(j = 1:Nsamples_MCERR, .combine = "rbind") %do% {
    m_sample = samples[j,]
    log_m_sample = log10(m_sample[m_sample > 0])
    hh = maghist(x = log_m_sample, breaks = sm_bins, plot = FALSE, verbose = FALSE)
    bin_dmf = hh$counts
    return(bin_dmf)
  }
  hist_quant = colQuantiles(hist_samples, probs = c(0.5, 0.16, 0.84))
}
mc_err_devilsd10_hybrid_dmf = foreach(i = 1:length(zmids)) %do% {
  message(i)
  vol = vol_mids[i] * 1.47 / (4*pi*(180/pi)^2)
  zidx = devilsd10_hybrid$z >= zbins[i] & devilsd10_hybrid$z < zbins[i+1]
  
  Mdust = devilsd10_hybrid$dustmass.total[zidx]
  MdustErr = 0.5 * (devilsd10_hybrid$dustmass.total_UB[zidx] - devilsd10_hybrid$dustmass.total_LB[zidx])
  
  samples = matrix(rnorm(n = Nsamples_MCERR*sum(zidx), mean = Mdust, sd = MdustErr), nrow = Nsamples_MCERR, ncol = sum(zidx))
  
  hist_samples = foreach(j = 1:Nsamples_MCERR, .combine = "rbind") %do% {
    m_sample = samples[j,]
    log_m_sample = log10(m_sample[m_sample > 0])
    hh = maghist(x = log_m_sample, breaks = sm_bins, plot = FALSE, verbose = FALSE)
    bin_dmf = hh$counts
    return(bin_dmf)
  }
  hist_quant = colQuantiles(hist_samples, probs = c(0.5, 0.16, 0.84))
}

mc_err_gama_noAGN_dmf_corr = foreach(i = 1:length(zmids)) %do% {
  message(i)
  vol = vol_mids[i] * 217.54 / (4*pi*(180/pi)^2)
  zidx = gama_noAGN$z >= zbins[i] & gama_noAGN$z < zbins[i+1]
  
  Mdust = gama_noAGN$DustMass_50[zidx]
  MdustErr = 0.5 * (gama_noAGN$DustMass_84[zidx] - gama_noAGN$DustMass_16[zidx])
  
  Ldust = gama_noAGN$DustLum_50[zidx]
  LdustErr = 0.5 * (gama_noAGN$DustLum_84[zidx] - gama_noAGN$DustLum_16[zidx])
  
  Z = gama_noAGN$Zgas_50[zidx]
  Zerr = 0.5 * (gama_noAGN$Zgas_84[zidx] - gama_noAGN$Zgas_16[zidx])
  
  Mdust_samples = matrix(rnorm(n = Nsamples_MCERR*sum(zidx), mean = Mdust, sd = MdustErr), nrow = Nsamples_MCERR, ncol = sum(zidx))
  Ldust_samples = matrix(rnorm(n = Nsamples_MCERR*sum(zidx), mean = Ldust, sd = LdustErr), nrow = Nsamples_MCERR, ncol = sum(zidx))
  Z_samples = matrix(rnorm(n = Nsamples_MCERR*sum(zidx), mean = log10(Z), sd = Zerr/(Z*log(10))), nrow = Nsamples_MCERR, ncol = sum(zidx))
  
  hist_samples = foreach(j = 1:Nsamples_MCERR, .combine = "rbind") %do% {
    alpha_sample = Dale_L2M(Ldust_samples[j,]/Mdust_samples[j,])
    new_dust_sample = Ldust_samples[j,] / (Dale_M2L_variableDTH_func(alpha_sample, qPAH_VSG = qPAH_VSG) * qPAHZ_corr_func(10^Z_samples[j,], alpha_sample))
    # new_dust_sample = Ldust_samples[j,] / Dale_M2L_variableDTH_func(alpha_sample, qPAH_VSG = qPAH_VSG)
    m_sample = new_dust_sample * RR14_BPL(Z = 10^Z_samples[j,])/0.0073
    log_m_sample = log10(m_sample[m_sample > 0])
    hh = maghist(x = log_m_sample, breaks = sm_bins, plot = FALSE, verbose = FALSE)
    bin_dmf = hh$counts
    return(bin_dmf)
  }
  hist_quant = colQuantiles(hist_samples, probs = c(0.5, 0.16, 0.84))
}
mc_err_gama_AGN_dmf_corr = foreach(i = 1:length(zmids)) %do% {
  message(i)
  vol = vol_mids[i] * 217.54 / (4*pi*(180/pi)^2)
  zidx = gama_AGN$z >= zbins[i] & gama_AGN$z < zbins[i+1]
  
  Mdust = gama_AGN$DustMass_50[zidx]
  MdustErr = 0.5 * (gama_AGN$DustMass_84[zidx] - gama_AGN$DustMass_16[zidx])
  
  Ldust = gama_AGN$DustLum_50[zidx]
  LdustErr = 0.5 * (gama_AGN$DustLum_84[zidx] - gama_AGN$DustLum_16[zidx])
  
  Z = gama_AGN$Zgas_50[zidx]
  Zerr = 0.5 * (gama_AGN$Zgas_84[zidx] - gama_AGN$Zgas_16[zidx])
  
  Mdust_samples = matrix(rnorm(n = Nsamples_MCERR*sum(zidx), mean = Mdust, sd = MdustErr), nrow = Nsamples_MCERR, ncol = sum(zidx))
  Ldust_samples = matrix(rnorm(n = Nsamples_MCERR*sum(zidx), mean = Ldust, sd = LdustErr), nrow = Nsamples_MCERR, ncol = sum(zidx))
  Z_samples = matrix(rnorm(n = Nsamples_MCERR*sum(zidx), mean = log10(Z), sd = Zerr/(Z*log(10))), nrow = Nsamples_MCERR, ncol = sum(zidx))
  
  hist_samples = foreach(j = 1:Nsamples_MCERR, .combine = "rbind") %do% {
    alpha_sample = Dale_L2M(Ldust_samples[j,]/Mdust_samples[j,])
    new_dust_sample = Ldust_samples[j,] / (Dale_M2L_variableDTH_func(alpha_sample, qPAH_VSG = qPAH_VSG) * qPAHZ_corr_func(10^Z_samples[j,], alpha_sample))
    m_sample = new_dust_sample * RR14_BPL(Z = 10^Z_samples[j,])/0.0073
    log_m_sample = log10(m_sample[m_sample > 0])
    hh = maghist(x = log_m_sample, breaks = sm_bins, plot = FALSE, verbose = FALSE)
    bin_dmf = hh$counts
    return(bin_dmf)
  }
  hist_quant = colQuantiles(hist_samples, probs = c(0.5, 0.16, 0.84))
}
mc_err_gama_hybrid_dmf_corr = foreach(i = 1:length(zmids)) %do% {
  message(i)
  vol = vol_mids[i] * 217.54 / (4*pi*(180/pi)^2)
  zidx = gama_hybrid$z >= zbins[i] & gama_hybrid$z < zbins[i+1]

  Mdust = gama_hybrid$DustMass_50[zidx]
  MdustErr = 0.5 * (gama_hybrid$DustMass_84[zidx] - gama_hybrid$DustMass_16[zidx])
  
  Ldust = gama_hybrid$DustLum_50[zidx]
  LdustErr = 0.5 * (gama_hybrid$DustLum_84[zidx] - gama_hybrid$DustLum_16[zidx])

  Z = gama_hybrid$Zgas_50[zidx]
  Zerr = 0.5 * (gama_hybrid$Zgas_84[zidx] - gama_hybrid$Zgas_16[zidx])
  
  Mdust_samples = matrix(rnorm(n = Nsamples_MCERR*sum(zidx), mean = Mdust, sd = MdustErr), nrow = Nsamples_MCERR, ncol = sum(zidx))
  Ldust_samples = matrix(rnorm(n = Nsamples_MCERR*sum(zidx), mean = Ldust, sd = LdustErr), nrow = Nsamples_MCERR, ncol = sum(zidx))
  Z_samples = matrix(rnorm(n = Nsamples_MCERR*sum(zidx), mean = log10(Z), sd = Zerr/(Z*log(10))), nrow = Nsamples_MCERR, ncol = sum(zidx))

  hist_samples = foreach(j = 1:Nsamples_MCERR, .combine = "rbind") %do% {
    alpha_sample = Dale_L2M(Ldust_samples[j,]/Mdust_samples[j,])
    new_dust_sample = Ldust_samples[j,] / (Dale_M2L_variableDTH_func(alpha_sample, qPAH_VSG = qPAH_VSG) * qPAHZ_corr_func(10^Z_samples[j,], alpha_sample))
    m_sample = new_dust_sample * RR14_BPL(Z = 10^Z_samples[j,])/0.0073
    log_m_sample = log10(m_sample[m_sample > 0])
    hh = maghist(x = log_m_sample, breaks = sm_bins, plot = FALSE, verbose = FALSE)
    bin_dmf = hh$counts
    return(bin_dmf)
  }
  hist_quant = colQuantiles(hist_samples, probs = c(0.5, 0.16, 0.84))
}
mc_err_gama_hybrid_dmf_corr_Zsol = foreach(i = 1:length(zmids)) %do% {
  message(i)
  vol = vol_mids[i] * 217.54 / (4*pi*(180/pi)^2)
  zidx = gama_hybrid$z >= zbins[i] & gama_hybrid$z < zbins[i+1]
  
  Mdust = gama_hybrid$DustMass_50[zidx]
  MdustErr = 0.5 * (gama_hybrid$DustMass_84[zidx] - gama_hybrid$DustMass_16[zidx])
  
  Ldust = gama_hybrid$DustLum_50[zidx]
  LdustErr = 0.5 * (gama_hybrid$DustLum_84[zidx] - gama_hybrid$DustLum_16[zidx])
  
  # Z = gama_hybrid$Zgas_50[zidx]
  # Zerr = 0.5 * (gama_hybrid$Zgas_84[zidx] - gama_hybrid$Zgas_16[zidx])
  
  Mdust_samples = matrix(rnorm(n = Nsamples_MCERR*sum(zidx), mean = Mdust, sd = MdustErr), nrow = Nsamples_MCERR, ncol = sum(zidx))
  Ldust_samples = matrix(rnorm(n = Nsamples_MCERR*sum(zidx), mean = Ldust, sd = LdustErr), nrow = Nsamples_MCERR, ncol = sum(zidx))
  # Z_samples = matrix(rnorm(n = Nsamples_MCERR*sum(zidx), mean = log10(Z), sd = Zerr/(Z*log(10))), nrow = Nsamples_MCERR, ncol = sum(zidx))
  
  hist_samples = foreach(j = 1:Nsamples_MCERR, .combine = "rbind") %do% {
    alpha_sample = Dale_L2M(Ldust_samples[j,]/Mdust_samples[j,])
    new_dust_sample = Ldust_samples[j,] / Dale_M2L_variableDTH_func(alpha_sample, qPAH_VSG = shivaei24_qPAHZ(0.014))
    m_sample = new_dust_sample * RR14_BPL(Z = 0.014)/0.0073
    log_m_sample = log10(m_sample[m_sample > 0])
    hh = maghist(x = log_m_sample, breaks = sm_bins, plot = FALSE, verbose = FALSE)
    bin_dmf = hh$counts
    return(bin_dmf)
  }
  hist_quant = colQuantiles(hist_samples, probs = c(0.5, 0.16, 0.84))
}
mc_err_devilsd10_noAGN_dmf_corr = foreach(i = 1:length(zmids)) %do% {
  message(i)
  vol = vol_mids[i] * 1.47 / (4*pi*(180/pi)^2)
  zidx = devilsd10_noAGN$z >= zbins[i] & devilsd10_noAGN$z < zbins[i+1]
  
  Mdust = devilsd10_noAGN$dustmass.total[zidx]
  MdustErr = 0.5 * (devilsd10_noAGN$dustmass.total_UB[zidx] - devilsd10_noAGN$dustmass.total_LB[zidx])
  
  Ldust = devilsd10_noAGN$dustlum.total[zidx]
  LdustErr = 0.5 * (devilsd10_noAGN$dustlum.total_UB[zidx] - devilsd10_noAGN$dustlum.total_LB[zidx])
  
  Z = devilsd10_noAGN$Zfinal[zidx]
  Zerr = 0.5 * (devilsd10_noAGN$Zfinal_UB[zidx] - devilsd10_noAGN$Zfinal_LB[zidx])
  
  Mdust_samples = matrix(rnorm(n = Nsamples_MCERR*sum(zidx), mean = Mdust, sd = MdustErr), nrow = Nsamples_MCERR, ncol = sum(zidx))
  Ldust_samples = matrix(rnorm(n = Nsamples_MCERR*sum(zidx), mean = Ldust, sd = LdustErr), nrow = Nsamples_MCERR, ncol = sum(zidx))
  Z_samples = matrix(rnorm(n = Nsamples_MCERR*sum(zidx), mean = Z, sd = Zerr), nrow = Nsamples_MCERR, ncol = sum(zidx))
  
  hist_samples = foreach(j = 1:Nsamples_MCERR, .combine = "rbind") %do% {
    alpha_sample = Dale_L2M(Ldust_samples[j,]/Mdust_samples[j,])
    new_dust_sample = Ldust_samples[j,] / (Dale_M2L_variableDTH_func(alpha_sample, qPAH_VSG = qPAH_VSG) * qPAHZ_corr_func(10^Z_samples[j,], alpha_sample))
    m_sample = new_dust_sample * RR14_BPL(Z = 10^Z_samples[j,])/0.0073
    log_m_sample = log10(m_sample[m_sample > 0])
    hh = maghist(x = log_m_sample, breaks = sm_bins, plot = FALSE, verbose = FALSE)
    bin_dmf = hh$counts
    return(bin_dmf)
  }
  hist_quant = colQuantiles(hist_samples, probs = c(0.5, 0.16, 0.84))
}
mc_err_devilsd10_AGN_dmf_corr = foreach(i = 1:length(zmids)) %do% {
  message(i)
  vol = vol_mids[i] * 1.47 / (4*pi*(180/pi)^2)
  zidx = devilsd10_AGN$z >= zbins[i] & devilsd10_AGN$z < zbins[i+1]
  
  Mdust = devilsd10_AGN$dustmass.total[zidx]
  MdustErr = 0.5 * (devilsd10_AGN$dustmass.total_UB[zidx] - devilsd10_AGN$dustmass.total_LB[zidx])
  
  Ldust = devilsd10_AGN$dustlum.total[zidx]
  LdustErr = 0.5 * (devilsd10_AGN$dustlum.total_UB[zidx] - devilsd10_AGN$dustlum.total_LB[zidx])
  
  Z = devilsd10_AGN$Zfinal[zidx]
  Zerr = 0.5 * (devilsd10_AGN$Zfinal_UB[zidx] - devilsd10_AGN$Zfinal_LB[zidx])
  
  Mdust_samples = matrix(rnorm(n = Nsamples_MCERR*sum(zidx), mean = Mdust, sd = MdustErr), nrow = Nsamples_MCERR, ncol = sum(zidx))
  Ldust_samples = matrix(rnorm(n = Nsamples_MCERR*sum(zidx), mean = Ldust, sd = LdustErr), nrow = Nsamples_MCERR, ncol = sum(zidx))
  Z_samples = matrix(rnorm(n = Nsamples_MCERR*sum(zidx), mean = Z, sd = Zerr), nrow = Nsamples_MCERR, ncol = sum(zidx))
  
  hist_samples = foreach(j = 1:Nsamples_MCERR, .combine = "rbind") %do% {
    alpha_sample = Dale_L2M(Ldust_samples[j,]/Mdust_samples[j,])
    new_dust_sample = Ldust_samples[j,] / (Dale_M2L_variableDTH_func(alpha_sample, qPAH_VSG = qPAH_VSG) * qPAHZ_corr_func(10^Z_samples[j,], alpha_sample))
    m_sample = new_dust_sample * RR14_BPL(Z = 10^Z_samples[j,])/0.0073
    log_m_sample = log10(m_sample[m_sample > 0])
    hh = maghist(x = log_m_sample, breaks = sm_bins, plot = FALSE, verbose = FALSE)
    bin_dmf = hh$counts
    return(bin_dmf)
  }
  hist_quant = colQuantiles(hist_samples, probs = c(0.5, 0.16, 0.84))
}
mc_err_devilsd10_hybrid_dmf_corr = foreach(i = 1:length(zmids)) %do% {
  message(i)
  vol = vol_mids[i] * 1.47 / (4*pi*(180/pi)^2)
  zidx = devilsd10_hybrid$z >= zbins[i] & devilsd10_hybrid$z < zbins[i+1]

  Mdust = devilsd10_hybrid$dustmass.total[zidx]
  MdustErr = 0.5 * (devilsd10_hybrid$dustmass.total_UB[zidx] - devilsd10_hybrid$dustmass.total_LB[zidx])
  
  Ldust = devilsd10_hybrid$dustlum.total[zidx]
  LdustErr = 0.5 * (devilsd10_hybrid$dustlum.total_UB[zidx] - devilsd10_hybrid$dustlum.total_LB[zidx])
  
  Z = devilsd10_hybrid$Zfinal[zidx]
  Zerr = 0.5 * (devilsd10_hybrid$Zfinal_UB[zidx] - devilsd10_hybrid$Zfinal_LB[zidx])
  
  Mdust_samples = matrix(rnorm(n = Nsamples_MCERR*sum(zidx), mean = Mdust, sd = MdustErr), nrow = Nsamples_MCERR, ncol = sum(zidx))
  Ldust_samples = matrix(rnorm(n = Nsamples_MCERR*sum(zidx), mean = Ldust, sd = LdustErr), nrow = Nsamples_MCERR, ncol = sum(zidx))
  Z_samples = matrix(rnorm(n = Nsamples_MCERR*sum(zidx), mean = Z, sd = Zerr), nrow = Nsamples_MCERR, ncol = sum(zidx))

  hist_samples = foreach(j = 1:Nsamples_MCERR, .combine = "rbind") %do% {
    alpha_sample = Dale_L2M(Ldust_samples[j,]/Mdust_samples[j,])
    new_dust_sample = Ldust_samples[j,] / (Dale_M2L_variableDTH_func(alpha_sample, qPAH_VSG = qPAH_VSG) * qPAHZ_corr_func(10^Z_samples[j,], alpha_sample))
    # new_dust_sample = Ldust_samples[j,] / Dale_M2L_variableDTH_func(alpha_sample, qPAH_VSG = qPAH_VSG)
    m_sample = new_dust_sample * RR14_BPL(Z = 10^Z_samples[j,])/0.0073
    log_m_sample = log10(m_sample[m_sample > 0])
    hh = maghist(x = log_m_sample, breaks = sm_bins, plot = FALSE, verbose = FALSE)
    bin_dmf = hh$counts
    return(bin_dmf)
  }
  hist_quant = colQuantiles(hist_samples, probs = c(0.5, 0.16, 0.84))
}
mc_err_devilsd10_hybrid_dmf_corr_Zsol = foreach(i = 1:length(zmids)) %do% {
  message(i)
  vol = vol_mids[i] * 1.47 / (4*pi*(180/pi)^2)
  zidx = devilsd10_hybrid$z >= zbins[i] & devilsd10_hybrid$z < zbins[i+1]
  
  Mdust = devilsd10_hybrid$dustmass.total[zidx]
  MdustErr = 0.5 * (devilsd10_hybrid$dustmass.total_UB[zidx] - devilsd10_hybrid$dustmass.total_LB[zidx])
  
  Ldust = devilsd10_hybrid$dustlum.total[zidx]
  LdustErr = 0.5 * (devilsd10_hybrid$dustlum.total_UB[zidx] - devilsd10_hybrid$dustlum.total_LB[zidx])
  
  # Z = devilsd10_hybrid$Zfinal[zidx]
  # Zerr = 0.5 * (devilsd10_hybrid$Zfinal_UB[zidx] - devilsd10_hybrid$Zfinal_LB[zidx])
  
  Mdust_samples = matrix(rnorm(n = Nsamples_MCERR*sum(zidx), mean = Mdust, sd = MdustErr), nrow = Nsamples_MCERR, ncol = sum(zidx))
  Ldust_samples = matrix(rnorm(n = Nsamples_MCERR*sum(zidx), mean = Ldust, sd = LdustErr), nrow = Nsamples_MCERR, ncol = sum(zidx))
  # Z_samples = matrix(rnorm(n = Nsamples_MCERR*sum(zidx), mean = Z, sd = Zerr), nrow = Nsamples_MCERR, ncol = sum(zidx))
  
  hist_samples = foreach(j = 1:Nsamples_MCERR, .combine = "rbind") %do% {
    alpha_sample = Dale_L2M(Ldust_samples[j,]/Mdust_samples[j,])
    new_dust_sample = Ldust_samples[j,] / Dale_M2L_variableDTH_func(alpha_sample, qPAH_VSG = shivaei24_qPAHZ(0.014))
    m_sample = new_dust_sample * RR14_BPL(Z = 0.014)/0.0073
    log_m_sample = log10(m_sample[m_sample > 0])
    hh = maghist(x = log_m_sample, breaks = sm_bins, plot = FALSE, verbose = FALSE)
    bin_dmf = hh$counts
    return(bin_dmf)
  }
  hist_quant = colQuantiles(hist_samples, probs = c(0.5, 0.16, 0.84))
}

mc_err_gama_noAGN_gmf_corr = foreach(i = 1:length(zmids)) %do% {
  message(i)
  vol = vol_mids[i] * 217.54 / (4*pi*(180/pi)^2)
  zidx = gama_noAGN$z >= zbins[i] & gama_noAGN$z < zbins[i+1]
  
  Mdust = gama_noAGN$DustMass_50[zidx]
  MdustErr = 0.5 * (gama_noAGN$DustMass_84[zidx] - gama_noAGN$DustMass_16[zidx])
  
  Ldust = gama_noAGN$DustLum_50[zidx]
  LdustErr = 0.5 * (gama_noAGN$DustLum_84[zidx] - gama_noAGN$DustLum_16[zidx])
  
  Z = gama_noAGN$Zgas_50[zidx]
  Zerr = 0.5 * (gama_noAGN$Zgas_84[zidx] - gama_noAGN$Zgas_16[zidx])
  
  Mdust_samples = matrix(rnorm(n = Nsamples_MCERR*sum(zidx), mean = Mdust, sd = MdustErr), nrow = Nsamples_MCERR, ncol = sum(zidx))
  Ldust_samples = matrix(rnorm(n = Nsamples_MCERR*sum(zidx), mean = Ldust, sd = LdustErr), nrow = Nsamples_MCERR, ncol = sum(zidx))
  Z_samples = matrix(rnorm(n = Nsamples_MCERR*sum(zidx), mean = log10(Z), sd = Zerr/(Z*log(10))), nrow = Nsamples_MCERR, ncol = sum(zidx))
  
  hist_samples = foreach(j = 1:Nsamples_MCERR, .combine = "rbind") %do% {
    alpha_sample = Dale_L2M(Ldust_samples[j,]/Mdust_samples[j,])
    new_dust_sample = Ldust_samples[j,] / (Dale_M2L_variableDTH_func(alpha_sample, qPAH_VSG = qPAH_VSG) * qPAHZ_corr_func(10^Z_samples[j,], alpha_sample))
    m_sample = new_dust_sample * RR14_BPL(Z = 10^Z_samples[j,])/0.0073/RR14_BPL(Z = 10^Z_samples[j,], doDTG = TRUE)
    log_m_sample = log10(m_sample[m_sample > 0])
    hh = maghist(x = log_m_sample, breaks = sm_bins, plot = FALSE, verbose = FALSE)
    bin_dmf = hh$counts
    return(bin_dmf)
  }
  hist_quant = colQuantiles(hist_samples, probs = c(0.5, 0.16, 0.84))
}
mc_err_gama_AGN_gmf_corr = foreach(i = 1:length(zmids)) %do% {
  message(i)
  vol = vol_mids[i] * 217.54 / (4*pi*(180/pi)^2)
  zidx = gama_AGN$z >= zbins[i] & gama_AGN$z < zbins[i+1]
  
  Mdust = gama_AGN$DustMass_50[zidx]
  MdustErr = 0.5 * (gama_AGN$DustMass_84[zidx] - gama_AGN$DustMass_16[zidx])
  
  Ldust = gama_AGN$DustLum_50[zidx]
  LdustErr = 0.5 * (gama_AGN$DustLum_84[zidx] - gama_AGN$DustLum_16[zidx])
  
  Z = gama_AGN$Zgas_50[zidx]
  Zerr = 0.5 * (gama_AGN$Zgas_84[zidx] - gama_AGN$Zgas_16[zidx])
  
  Mdust_samples = matrix(rnorm(n = Nsamples_MCERR*sum(zidx), mean = Mdust, sd = MdustErr), nrow = Nsamples_MCERR, ncol = sum(zidx))
  Ldust_samples = matrix(rnorm(n = Nsamples_MCERR*sum(zidx), mean = Ldust, sd = LdustErr), nrow = Nsamples_MCERR, ncol = sum(zidx))
  Z_samples = matrix(rnorm(n = Nsamples_MCERR*sum(zidx), mean = log10(Z), sd = Zerr/(Z*log(10))), nrow = Nsamples_MCERR, ncol = sum(zidx))
  
  hist_samples = foreach(j = 1:Nsamples_MCERR, .combine = "rbind") %do% {
    alpha_sample = Dale_L2M(Ldust_samples[j,]/Mdust_samples[j,])
    new_dust_sample = Ldust_samples[j,] / (Dale_M2L_variableDTH_func(alpha_sample, qPAH_VSG = qPAH_VSG) * qPAHZ_corr_func(10^Z_samples[j,], alpha_sample))
    m_sample = new_dust_sample * RR14_BPL(Z = 10^Z_samples[j,])/0.0073/RR14_BPL(Z = 10^Z_samples[j,], doDTG = TRUE)
    log_m_sample = log10(m_sample[m_sample > 0])
    hh = maghist(x = log_m_sample, breaks = sm_bins, plot = FALSE, verbose = FALSE)
    bin_dmf = hh$counts
    return(bin_dmf)
  }
  hist_quant = colQuantiles(hist_samples, probs = c(0.5, 0.16, 0.84))
}
mc_err_gama_hybrid_gmf_corr = foreach(i = 1:length(zmids)) %do% {
  message(i)
  vol = vol_mids[i] * 217.54 / (4*pi*(180/pi)^2)
  zidx = gama_hybrid$z >= zbins[i] & gama_hybrid$z < zbins[i+1]
  
  Mdust = gama_hybrid$DustMass_50[zidx]
  MdustErr = 0.5 * (gama_hybrid$DustMass_84[zidx] - gama_hybrid$DustMass_16[zidx])
  
  Ldust = gama_hybrid$DustLum_50[zidx]
  LdustErr = 0.5 * (gama_hybrid$DustLum_84[zidx] - gama_hybrid$DustLum_16[zidx])
  
  Z = gama_hybrid$Zgas_50[zidx]
  Zerr = 0.5 * (gama_hybrid$Zgas_84[zidx] - gama_hybrid$Zgas_16[zidx])
  
  Mdust_samples = matrix(rnorm(n = Nsamples_MCERR*sum(zidx), mean = Mdust, sd = MdustErr), nrow = Nsamples_MCERR, ncol = sum(zidx))
  Ldust_samples = matrix(rnorm(n = Nsamples_MCERR*sum(zidx), mean = Ldust, sd = LdustErr), nrow = Nsamples_MCERR, ncol = sum(zidx))
  Z_samples = matrix(rnorm(n = Nsamples_MCERR*sum(zidx), mean = log10(Z), sd = Zerr/(Z*log(10))), nrow = Nsamples_MCERR, ncol = sum(zidx))
  
  hist_samples = foreach(j = 1:Nsamples_MCERR, .combine = "rbind") %do% {
    alpha_sample = Dale_L2M(Ldust_samples[j,]/Mdust_samples[j,])
    new_dust_sample = Ldust_samples[j,] / (Dale_M2L_variableDTH_func(alpha_sample, qPAH_VSG = qPAH_VSG) * qPAHZ_corr_func(10^Z_samples[j,], alpha_sample))
    m_sample = new_dust_sample * RR14_BPL(Z = 10^Z_samples[j,])/0.0073/RR14_BPL(Z = 10^Z_samples[j,], doDTG = TRUE)
    log_m_sample = log10(m_sample[m_sample > 0])
    hh = maghist(x = log_m_sample, breaks = sm_bins, plot = FALSE, verbose = FALSE)
    bin_dmf = hh$counts
    return(bin_dmf)
  }
  hist_quant = colQuantiles(hist_samples, probs = c(0.5, 0.16, 0.84))
}
mc_err_devilsd10_noAGN_gmf_corr = foreach(i = 1:length(zmids)) %do% {
  message(i)
  vol = vol_mids[i] * 1.47 / (4*pi*(180/pi)^2)
  zidx = devilsd10_noAGN$z >= zbins[i] & devilsd10_noAGN$z < zbins[i+1]
  
  Mdust = devilsd10_noAGN$dustmass.total[zidx]
  MdustErr = 0.5 * (devilsd10_noAGN$dustmass.total_UB[zidx] - devilsd10_noAGN$dustmass.total_LB[zidx])
  
  Ldust = devilsd10_noAGN$dustlum.total[zidx]
  LdustErr = 0.5 * (devilsd10_noAGN$dustlum.total_UB[zidx] - devilsd10_noAGN$dustlum.total_LB[zidx])
  
  Z = devilsd10_noAGN$Zfinal[zidx]
  Zerr = 0.5 * (devilsd10_noAGN$Zfinal_UB[zidx] - devilsd10_noAGN$Zfinal_LB[zidx])
  
  Mdust_samples = matrix(rnorm(n = Nsamples_MCERR*sum(zidx), mean = Mdust, sd = MdustErr), nrow = Nsamples_MCERR, ncol = sum(zidx))
  Ldust_samples = matrix(rnorm(n = Nsamples_MCERR*sum(zidx), mean = Ldust, sd = LdustErr), nrow = Nsamples_MCERR, ncol = sum(zidx))
  Z_samples = matrix(rnorm(n = Nsamples_MCERR*sum(zidx), mean = Z, sd = Zerr), nrow = Nsamples_MCERR, ncol = sum(zidx))
  
  hist_samples = foreach(j = 1:Nsamples_MCERR, .combine = "rbind") %do% {
    alpha_sample = Dale_L2M(Ldust_samples[j,]/Mdust_samples[j,])
    new_dust_sample = Ldust_samples[j,] / (Dale_M2L_variableDTH_func(alpha_sample, qPAH_VSG = qPAH_VSG) * qPAHZ_corr_func(10^Z_samples[j,], alpha_sample))
    m_sample = new_dust_sample * RR14_BPL(Z = 10^Z_samples[j,])/0.0073/RR14_BPL(Z = 10^Z_samples[j,], doDTG = TRUE)
    log_m_sample = log10(m_sample[m_sample > 0])
    hh = maghist(x = log_m_sample, breaks = sm_bins, plot = FALSE, verbose = FALSE)
    bin_dmf = hh$counts
    return(bin_dmf)
  }
  hist_quant = colQuantiles(hist_samples, probs = c(0.5, 0.16, 0.84))
}
mc_err_devilsd10_AGN_gmf_corr = foreach(i = 1:length(zmids)) %do% {
  message(i)
  vol = vol_mids[i] * 1.47 / (4*pi*(180/pi)^2)
  zidx = devilsd10_AGN$z >= zbins[i] & devilsd10_AGN$z < zbins[i+1]
  
  Mdust = devilsd10_AGN$dustmass.total[zidx]
  MdustErr = 0.5 * (devilsd10_AGN$dustmass.total_UB[zidx] - devilsd10_AGN$dustmass.total_LB[zidx])
  
  Ldust = devilsd10_AGN$dustlum.total[zidx]
  LdustErr = 0.5 * (devilsd10_AGN$dustlum.total_UB[zidx] - devilsd10_AGN$dustlum.total_LB[zidx])
  
  Z = devilsd10_AGN$Zfinal[zidx]
  Zerr = 0.5 * (devilsd10_AGN$Zfinal_UB[zidx] - devilsd10_AGN$Zfinal_LB[zidx])
  
  Mdust_samples = matrix(rnorm(n = Nsamples_MCERR*sum(zidx), mean = Mdust, sd = MdustErr), nrow = Nsamples_MCERR, ncol = sum(zidx))
  Ldust_samples = matrix(rnorm(n = Nsamples_MCERR*sum(zidx), mean = Ldust, sd = LdustErr), nrow = Nsamples_MCERR, ncol = sum(zidx))
  Z_samples = matrix(rnorm(n = Nsamples_MCERR*sum(zidx), mean = Z, sd = Zerr), nrow = Nsamples_MCERR, ncol = sum(zidx))
  
  hist_samples = foreach(j = 1:Nsamples_MCERR, .combine = "rbind") %do% {
    alpha_sample = Dale_L2M(Ldust_samples[j,]/Mdust_samples[j,])
    new_dust_sample = Ldust_samples[j,] / (Dale_M2L_variableDTH_func(alpha_sample, qPAH_VSG = qPAH_VSG) * qPAHZ_corr_func(10^Z_samples[j,], alpha_sample))
    m_sample = new_dust_sample * RR14_BPL(Z = 10^Z_samples[j,])/0.0073/RR14_BPL(Z = 10^Z_samples[j,], doDTG = TRUE)
    log_m_sample = log10(m_sample[m_sample > 0])
    hh = maghist(x = log_m_sample, breaks = sm_bins, plot = FALSE, verbose = FALSE)
    bin_dmf = hh$counts
    return(bin_dmf)
  }
  hist_quant = colQuantiles(hist_samples, probs = c(0.5, 0.16, 0.84))
}
mc_err_devilsd10_hybrid_gmf_corr = foreach(i = 1:length(zmids)) %do% {
  message(i)
  vol = vol_mids[i] * 1.47 / (4*pi*(180/pi)^2)
  zidx = devilsd10_hybrid$z >= zbins[i] & devilsd10_hybrid$z < zbins[i+1]
  
  Mdust = devilsd10_hybrid$dustmass.total[zidx]
  MdustErr = 0.5 * (devilsd10_hybrid$dustmass.total_UB[zidx] - devilsd10_hybrid$dustmass.total_LB[zidx])
  
  Ldust = devilsd10_hybrid$dustlum.total[zidx]
  LdustErr = 0.5 * (devilsd10_hybrid$dustlum.total_UB[zidx] - devilsd10_hybrid$dustlum.total_LB[zidx])
  
  Z = devilsd10_hybrid$Zfinal[zidx]
  Zerr = 0.5 * (devilsd10_hybrid$Zfinal_UB[zidx] - devilsd10_hybrid$Zfinal_LB[zidx])
  
  Mdust_samples = matrix(rnorm(n = Nsamples_MCERR*sum(zidx), mean = Mdust, sd = MdustErr), nrow = Nsamples_MCERR, ncol = sum(zidx))
  Ldust_samples = matrix(rnorm(n = Nsamples_MCERR*sum(zidx), mean = Ldust, sd = LdustErr), nrow = Nsamples_MCERR, ncol = sum(zidx))
  Z_samples = matrix(rnorm(n = Nsamples_MCERR*sum(zidx), mean = Z, sd = Zerr), nrow = Nsamples_MCERR, ncol = sum(zidx))
  
  hist_samples = foreach(j = 1:Nsamples_MCERR, .combine = "rbind") %do% {
    alpha_sample = Dale_L2M(Ldust_samples[j,]/Mdust_samples[j,])
    new_dust_sample = Ldust_samples[j,] / (Dale_M2L_variableDTH_func(alpha_sample, qPAH_VSG = qPAH_VSG) * qPAHZ_corr_func(10^Z_samples[j,], alpha_sample))
    m_sample = new_dust_sample * RR14_BPL(Z = 10^Z_samples[j,])/0.0073/RR14_BPL(Z = 10^Z_samples[j,], doDTG = TRUE)
    log_m_sample = log10(m_sample[m_sample > 0])
    hh = maghist(x = log_m_sample, breaks = sm_bins, plot = FALSE, verbose = FALSE)
    bin_dmf = hh$counts
    return(bin_dmf)
  }
  hist_quant = colQuantiles(hist_samples, probs = c(0.5, 0.16, 0.84))
}

mc_err_gama_noAGN_smf = foreach(i = 1:length(zmids)) %do% {
  message(i)
  vol = vol_mids[i] * 217.54 / (4*pi*(180/pi)^2)
  zidx = gama_noAGN$z >= zbins[i] & gama_noAGN$z < zbins[i+1]
  
  Mdust = gama_noAGN$StellarMass_50[zidx]
  MdustErr = 0.5 * (gama_noAGN$StellarMass_84[zidx] - gama_noAGN$StellarMass_84[zidx])
  
  samples = matrix(rnorm(n = Nsamples_MCERR*sum(zidx), mean = Mdust, sd = MdustErr), nrow = Nsamples_MCERR, ncol = sum(zidx))
  
  hist_samples = foreach(j = 1:Nsamples_MCERR, .combine = "rbind") %do% {
    m_sample = samples[j,]
    log_m_sample = log10(m_sample[m_sample > 0])
    hh = maghist(x = log_m_sample, breaks = sm_bins, plot = FALSE, verbose = FALSE)
    bin_dmf = hh$counts
    return(bin_dmf)
  }
  hist_quant = colQuantiles(hist_samples, probs = c(0.5, 0.16, 0.84))
}
mc_err_gama_AGN_smf = foreach(i = 1:length(zmids)) %do% {
  message(i)
  vol = vol_mids[i] * 217.54 / (4*pi*(180/pi)^2)
  zidx = gama_AGN$z >= zbins[i] & gama_AGN$z < zbins[i+1]
  
  Mdust = gama_AGN$StellarMass_50[zidx]
  MdustErr = 0.5 * (gama_AGN$StellarMass_84[zidx] - gama_AGN$StellarMass_84[zidx])
  
  samples = matrix(rnorm(n = Nsamples_MCERR*sum(zidx), mean = Mdust, sd = MdustErr), nrow = Nsamples_MCERR, ncol = sum(zidx))
  
  hist_samples = foreach(j = 1:Nsamples_MCERR, .combine = "rbind") %do% {
    m_sample = samples[j,]
    log_m_sample = log10(m_sample[m_sample > 0])
    hh = maghist(x = log_m_sample, breaks = sm_bins, plot = FALSE, verbose = FALSE)
    bin_dmf = hh$counts
    return(bin_dmf)
  }
  hist_quant = colQuantiles(hist_samples, probs = c(0.5, 0.16, 0.84))
}
mc_err_gama_hybrid_smf = foreach(i = 1:length(zmids)) %do% {
  message(i)
  vol = vol_mids[i] * 217.54 / (4*pi*(180/pi)^2)
  zidx = gama_hybrid$z >= zbins[i] & gama_hybrid$z < zbins[i+1]
  
  Mdust = gama_hybrid$StellarMass_50[zidx]
  MdustErr = 0.5 * (gama_hybrid$StellarMass_84[zidx] - gama_hybrid$StellarMass_84[zidx])
  
  samples = matrix(rnorm(n = Nsamples_MCERR*sum(zidx), mean = Mdust, sd = MdustErr), nrow = Nsamples_MCERR, ncol = sum(zidx))
  
  hist_samples = foreach(j = 1:Nsamples_MCERR, .combine = "rbind") %do% {
    m_sample = samples[j,]
    log_m_sample = log10(m_sample[m_sample > 0])
    hh = maghist(x = log_m_sample, breaks = sm_bins, plot = FALSE, verbose = FALSE)
    bin_dmf = hh$counts
    return(bin_dmf)
  }
  hist_quant = colQuantiles(hist_samples, probs = c(0.5, 0.16, 0.84))
}
mc_err_devilsd10_noAGN_smf = foreach(i = 1:length(zmids)) %do% {
  message(i)
  vol = vol_mids[i] * 1.47 / (4*pi*(180/pi)^2)
  zidx = devilsd10_noAGN$z >= zbins[i] & devilsd10_noAGN$z < zbins[i+1]
  
  Mdust = devilsd10_noAGN$StellarMass[zidx]
  MdustErr = 0.5 * (devilsd10_noAGN$StellarMass_UB[zidx] - devilsd10_noAGN$StellarMass_LB[zidx])
  
  samples = matrix(rnorm(n = Nsamples_MCERR*sum(zidx), mean = Mdust, sd = MdustErr), nrow = Nsamples_MCERR, ncol = sum(zidx))
  
  hist_samples = foreach(j = 1:Nsamples_MCERR, .combine = "rbind") %do% {
    m_sample = samples[j,]
    log_m_sample = log10(m_sample[m_sample > 0])
    hh = maghist(x = log_m_sample, breaks = sm_bins, plot = FALSE, verbose = FALSE)
    bin_dmf = hh$counts
    return(bin_dmf)
  }
  hist_quant = colQuantiles(hist_samples, probs = c(0.5, 0.16, 0.84))
}
mc_err_devilsd10_AGN_smf = foreach(i = 1:length(zmids)) %do% {
  message(i)
  vol = vol_mids[i] * 1.47 / (4*pi*(180/pi)^2)
  zidx = devilsd10_AGN$z >= zbins[i] & devilsd10_AGN$z < zbins[i+1]
  
  Mdust = devilsd10_AGN$StellarMass[zidx]
  MdustErr = 0.5 * (devilsd10_AGN$StellarMass_UB[zidx] - devilsd10_AGN$StellarMass_LB[zidx])
  
  samples = matrix(rnorm(n = Nsamples_MCERR*sum(zidx), mean = Mdust, sd = MdustErr), nrow = Nsamples_MCERR, ncol = sum(zidx))
  
  hist_samples = foreach(j = 1:Nsamples_MCERR, .combine = "rbind") %do% {
    m_sample = samples[j,]
    log_m_sample = log10(m_sample[m_sample > 0])
    hh = maghist(x = log_m_sample, breaks = sm_bins, plot = FALSE, verbose = FALSE)
    bin_dmf = hh$counts
    return(bin_dmf)
  }
  hist_quant = colQuantiles(hist_samples, probs = c(0.5, 0.16, 0.84))
}
mc_err_devilsd10_hybrid_smf = foreach(i = 1:length(zmids)) %do% {
  message(i)
  vol = vol_mids[i] * 1.47 / (4*pi*(180/pi)^2)
  zidx = devilsd10_hybrid$z >= zbins[i] & devilsd10_hybrid$z < zbins[i+1]
  
  Mdust = devilsd10_hybrid$StellarMass[zidx]
  MdustErr = 0.5 * (devilsd10_hybrid$StellarMass_UB[zidx] - devilsd10_hybrid$StellarMass_LB[zidx])
  
  samples = matrix(rnorm(n = Nsamples_MCERR*sum(zidx), mean = Mdust, sd = MdustErr), nrow = Nsamples_MCERR, ncol = sum(zidx))
  
  hist_samples = foreach(j = 1:Nsamples_MCERR, .combine = "rbind") %do% {
    m_sample = samples[j,]
    log_m_sample = log10(m_sample[m_sample > 0])
    hh = maghist(x = log_m_sample, breaks = sm_bins, plot = FALSE, verbose = FALSE)
    bin_dmf = hh$counts
    return(bin_dmf)
  }
  hist_quant = colQuantiles(hist_samples, probs = c(0.5, 0.16, 0.84))
}

## Save the mc_errs
save(
  mc_err_gama_noAGN_dmf, mc_err_gama_AGN_dmf, mc_err_gama_hybrid_dmf,
  mc_err_devilsd10_noAGN_dmf, mc_err_devilsd10_AGN_dmf, mc_err_devilsd10_hybrid_dmf, 
  mc_err_gama_noAGN_dmf_corr, mc_err_gama_AGN_dmf_corr, mc_err_gama_hybrid_dmf_corr, 
  mc_err_devilsd10_noAGN_dmf_corr, mc_err_devilsd10_AGN_dmf_corr, mc_err_devilsd10_hybrid_dmf_corr, 
  mc_err_gama_noAGN_gmf_corr, mc_err_gama_AGN_gmf_corr, mc_err_gama_hybrid_gmf_corr, 
  mc_err_devilsd10_noAGN_gmf_corr, mc_err_devilsd10_AGN_gmf_corr, mc_err_devilsd10_hybrid_gmf_corr, 
  mc_err_gama_noAGN_smf, mc_err_gama_AGN_smf, mc_err_gama_hybrid_smf, 
  mc_err_devilsd10_noAGN_smf, mc_err_devilsd10_AGN_smf, mc_err_devilsd10_hybrid_smf,
  file = "~/Documents/DustMassDensity/data/dmf_mc_errs.Rdata"
)
load("~/Documents/DustMassDensity/data/dmf_mc_errs.Rdata")

# Use constant corr factor for GAMA and alpha corrections for DEVILS
err_floor = 0.1
do_MC_err = 1.0 # make 1 to do 
gama_noAGN_dmf = foreach(i = 1:length(zmids)) %do% {
  message(i)
  vol = vol_mids[i] * 217.54 / (4*pi*(180/pi)^2)
  zidx = gama_noAGN$z >= zbins[i] & gama_noAGN$z < zbins[i+1]
  Mdust = gama_noAGN$DustMass_50[zidx]
  mc_err = 0.5*(mc_err_gama_noAGN_dmf[[i]][,3] - mc_err_gama_noAGN_dmf[[i]][,2])

  log_m = log10(Mdust[Mdust > 0])
  hh = maghist(x = log_m, breaks = sm_bins, plot = FALSE, verbose = FALSE)
  bin_dmf = hh$counts/(vol * diff(sm_bins))
  pois_dmf = sqrt(hh$counts)/(vol * diff(sm_bins))
  
  err_dmf = sqrt( sqrt(hh$counts)^2 + (err_floor*hh$counts)^2 + do_MC_err*mc_err^2 )/(vol * diff(sm_bins))
  df = data.frame(cbind(bin_dmf, err_dmf))
  names(df) = c("phi", "err")
  df$x = hh$mids
  return(df)
}
gama_AGN_dmf = foreach(i = 1:length(zmids)) %do% {
  message(i)
  vol = vol_mids[i] * 217.54 / (4*pi*(180/pi)^2)
  zidx = gama_AGN$z >= zbins[i] & gama_AGN$z < zbins[i+1]
  Mdust = gama_AGN$DustMass_50[zidx]
  mc_err = 0.5*(mc_err_gama_AGN_dmf[[i]][,3] - mc_err_gama_AGN_dmf[[i]][,2])
  
  log_m = log10(Mdust[Mdust > 0])
  hh = maghist(x = log_m, breaks = sm_bins, plot = FALSE, verbose = FALSE)
  bin_dmf = hh$counts/(vol * diff(sm_bins))
  pois_dmf = sqrt(hh$counts)/(vol * diff(sm_bins))
  
  err_dmf = sqrt( sqrt(hh$counts)^2 + (err_floor*hh$counts)^2 + do_MC_err*mc_err^2 )/(vol * diff(sm_bins))
  df = data.frame(cbind(bin_dmf, err_dmf))
  names(df) = c("phi", "err")
  df$x = hh$mids
  return(df)
}
devilsd10_noAGN_dmf = foreach(i = 1:length(zmids)) %do% {
  message(i)
  vol = vol_mids[i] * 1.47 / (4*pi*(180/pi)^2)
  zidx = devilsd10_noAGN$z >= zbins[i] & devilsd10_noAGN$z < zbins[i+1]
  Mdust = devilsd10_noAGN$dustmass.total[zidx]
  mc_err = 0.5*(mc_err_devilsd10_noAGN_dmf[[i]][,3] - mc_err_devilsd10_noAGN_dmf[[i]][,2])
  
  log_m = log10(Mdust[Mdust > 0])
  hh = maghist(x = log_m, breaks = sm_bins, plot = FALSE, verbose = FALSE)
  bin_dmf = hh$counts/(vol * diff(sm_bins))
  pois_dmf = sqrt(hh$counts)/(vol * diff(sm_bins))
  
  err_dmf = sqrt( sqrt(hh$counts)^2 + (err_floor*hh$counts)^2 + do_MC_err*mc_err^2 )/(vol * diff(sm_bins))
  df = data.frame(cbind(bin_dmf, err_dmf))
  names(df) = c("phi", "err")
  df$x = hh$mids
  return(df)
}
devilsd10_AGN_dmf = foreach(i = 1:length(zmids)) %do% {
  message(i)
  vol = vol_mids[i] * 1.47 / (4*pi*(180/pi)^2)
  zidx = devilsd10_AGN$z >= zbins[i] & devilsd10_AGN$z < zbins[i+1]
  Mdust = devilsd10_AGN$dustmass.total[zidx]
  mc_err = 0.5*(mc_err_devilsd10_AGN_dmf[[i]][,3] - mc_err_devilsd10_AGN_dmf[[i]][,2])
  
  log_m = log10(Mdust[Mdust > 0])
  hh = maghist(x = log_m, breaks = sm_bins, plot = FALSE, verbose = FALSE)
  bin_dmf = hh$counts/(vol * diff(sm_bins))
  pois_dmf = sqrt(hh$counts)/(vol * diff(sm_bins))
  
  err_dmf = sqrt( sqrt(hh$counts)^2 + (err_floor*hh$counts)^2 + do_MC_err*mc_err^2 )/(vol * diff(sm_bins))
  df = data.frame(cbind(bin_dmf, err_dmf))
  names(df) = c("phi", "err")
  df$x = hh$mids
  return(df)
}
gama_hybrid_dmf = foreach(i = 1:length(zmids)) %do% {
  message(i)
  vol = vol_mids[i] * 217.54 / (4*pi*(180/pi)^2)
  zidx = gama_hybrid$z >= zbins[i] & gama_hybrid$z < zbins[i+1]
  Mdust = gama_hybrid$DustMass_50[zidx]
  mc_err = 0.5*(mc_err_gama_hybrid_dmf[[i]][,3] - mc_err_gama_hybrid_dmf[[i]][,2])
  
  log_m = log10(Mdust[Mdust > 0])
  hh = maghist(x = log_m, breaks = sm_bins, plot = FALSE, verbose = FALSE)
  bin_dmf = hh$counts/(vol * diff(sm_bins))
  pois_dmf = sqrt(hh$counts)/(vol * diff(sm_bins))
  
  err_dmf = sqrt( sqrt(hh$counts)^2 + (err_floor*hh$counts)^2 + do_MC_err*mc_err^2 )/(vol * diff(sm_bins))
  df = data.frame(cbind(bin_dmf, err_dmf))
  names(df) = c("phi", "err")
  df$x = hh$mids
  return(df)
}
devilsd10_hybrid_dmf = foreach(i = 1:length(zmids)) %do% {
  message(i)
  vol = vol_mids[i] * 1.47 / (4*pi*(180/pi)^2)
  zidx = devilsd10_hybrid$z >= zbins[i] & devilsd10_hybrid$z < zbins[i+1]
  Mdust = devilsd10_hybrid$dustmass.total[zidx]
  mc_err = 0.5*(mc_err_devilsd10_hybrid_dmf[[i]][,3] - mc_err_devilsd10_hybrid_dmf[[i]][,2])
  
  log_m = log10(Mdust[Mdust > 0])
  hh = maghist(x = log_m, breaks = sm_bins, plot = FALSE, verbose = FALSE)
  bin_dmf = hh$counts/(vol * diff(sm_bins))
  pois_dmf = sqrt(hh$counts)/(vol * diff(sm_bins))
  
  err_dmf = sqrt( sqrt(hh$counts)^2 + (err_floor*hh$counts)^2 + do_MC_err*mc_err^2 )/(vol * diff(sm_bins))
  df = data.frame(cbind(bin_dmf, err_dmf))
  names(df) = c("phi", "err")
  df$x = hh$mids
  return(df)
}

gama_noAGN_dmf_corr = foreach(i = 1:length(zmids)) %do% {
  message(i)
  vol = vol_mids[i] * 217.54 / (4*pi*(180/pi)^2)
  zidx = gama_noAGN$z >= zbins[i] & gama_noAGN$z < zbins[i+1]
  
  gama_alpha = Dale_L2M(gama_noAGN$DustLum_50[zidx]/gama_noAGN$DustMass_50[zidx])
  Mdust_new = gama_noAGN$DustLum_50[zidx] / ( Dale_M2L_variableDTH_func(gama_alpha, qPAH_VSG = qPAH_VSG) * qPAHZ_corr_func(gama_noAGN$Zgas_50[zidx], alpha = gama_alpha) )
  
  Mdust = Mdust_new * RR14_BPL(Z = gama_noAGN$Zgas_50[zidx], doDTG = FALSE)/0.0073
  mc_err = 0.5*(mc_err_gama_noAGN_dmf_corr[[i]][,3] - mc_err_gama_noAGN_dmf_corr[[i]][,2])
  
  log_m = log10(Mdust[Mdust > 0])
  hh = maghist(x = log_m, breaks = sm_bins, plot = FALSE, verbose = FALSE)
  bin_dmf = hh$counts/(vol * diff(sm_bins))
  pois_dmf = sqrt(hh$counts)/(vol * diff(sm_bins))
  
  err_dmf = sqrt( sqrt(hh$counts)^2 + (err_floor*hh$counts)^2 + do_MC_err*mc_err^2 )/(vol * diff(sm_bins))
  df = data.frame(cbind(bin_dmf, err_dmf))
  names(df) = c("phi", "err")
  df$x = hh$mids
  return(df)
}
gama_AGN_dmf_corr = foreach(i = 1:length(zmids)) %do% {
  message(i)
  vol = vol_mids[i] * 217.54 / (4*pi*(180/pi)^2)
  zidx = gama_AGN$z >= zbins[i] & gama_AGN$z < zbins[i+1]
  
  gama_alpha = Dale_L2M(gama_AGN$DustLum_50[zidx]/gama_AGN$DustMass_50[zidx])
  Mdust_new = gama_AGN$DustLum_50[zidx] / ( Dale_M2L_variableDTH_func(gama_alpha, qPAH_VSG = qPAH_VSG) * qPAHZ_corr_func(gama_AGN$Zgas_50[zidx], alpha = gama_alpha) )
  
  Mdust = Mdust_new * RR14_BPL(Z = gama_AGN$Zgas_50[zidx], doDTG = FALSE)/0.0073
  mc_err = 0.5*(mc_err_gama_AGN_dmf_corr[[i]][,3] - mc_err_gama_AGN_dmf_corr[[i]][,2])
  
  log_m = log10(Mdust[Mdust > 0])
  hh = maghist(x = log_m, breaks = sm_bins, plot = FALSE, verbose = FALSE)
  bin_dmf = hh$counts/(vol * diff(sm_bins))
  pois_dmf = sqrt(hh$counts)/(vol * diff(sm_bins))
  
  err_dmf = sqrt( sqrt(hh$counts)^2 + (err_floor*hh$counts)^2 + do_MC_err*mc_err^2 )/(vol * diff(sm_bins))
  df = data.frame(cbind(bin_dmf, err_dmf))
  names(df) = c("phi", "err")
  df$x = hh$mids
  return(df)
}
gama_hybrid_dmf_corr = foreach(i = 1:length(zmids)) %do% {
  message(i)
  vol = vol_mids[i] * 217.54 / (4*pi*(180/pi)^2)
  zidx = gama_hybrid$z >= zbins[i] & gama_hybrid$z < zbins[i+1]
  
  gama_alpha = Dale_L2M(gama_hybrid$DustLum_50[zidx]/gama_hybrid$DustMass_50[zidx])
  Mdust_new = gama_hybrid$DustLum_50[zidx]/( Dale_M2L_variableDTH_func(gama_alpha, qPAH_VSG = qPAH_VSG) * qPAHZ_corr_func(gama_hybrid$Zgas_50[zidx], alpha = gama_alpha) )
  
  Mdust = Mdust_new * RR14_BPL(Z = gama_hybrid$Zgas_50[zidx], doDTG = FALSE)/0.0073
  mc_err = 0.5*(mc_err_gama_hybrid_dmf_corr[[i]][,3] - mc_err_gama_hybrid_dmf_corr[[i]][,2])
  
  log_m = log10(Mdust[Mdust > 0])
  hh = maghist(x = log_m, breaks = sm_bins, plot = FALSE, verbose = FALSE)
  bin_dmf = hh$counts/(vol * diff(sm_bins))
  pois_dmf = sqrt(hh$counts)/(vol * diff(sm_bins))
  
  err_dmf = sqrt( sqrt(hh$counts)^2 + (err_floor*hh$counts)^2 + do_MC_err*mc_err^2 )/(vol * diff(sm_bins))
  df = data.frame(cbind(bin_dmf, err_dmf))
  names(df) = c("phi", "err")
  df$x = hh$mids
  return(df)
}
gama_hybrid_dmf_corr_Zsol = foreach(i = 1:length(zmids)) %do% {
  message(i)
  vol = vol_mids[i] * 217.54 / (4*pi*(180/pi)^2)
  zidx = gama_hybrid$z >= zbins[i] & gama_hybrid$z < zbins[i+1]
  
  gama_alpha = Dale_L2M(gama_hybrid$DustLum_50[zidx]/gama_hybrid$DustMass_50[zidx])
  Mdust_new = gama_hybrid$DustLum_50[zidx]/Dale_M2L_variableDTH_func(gama_alpha, qPAH_VSG = shivaei24_qPAHZ(0.014))
  
  Mdust = Mdust_new * RR14_BPL(Z = 0.014, doDTG = FALSE)/0.0073
  mc_err = 0.5*(mc_err_gama_hybrid_dmf_corr_Zsol[[i]][,3] - mc_err_gama_hybrid_dmf_corr_Zsol[[i]][,2])
  
  log_m = log10(Mdust[Mdust > 0])
  hh = maghist(x = log_m, breaks = sm_bins, plot = FALSE, verbose = FALSE)
  bin_dmf = hh$counts/(vol * diff(sm_bins))
  pois_dmf = sqrt(hh$counts)/(vol * diff(sm_bins))
  
  err_dmf = sqrt( sqrt(hh$counts)^2 + (err_floor*hh$counts)^2 + do_MC_err*mc_err^2 )/(vol * diff(sm_bins))
  df = data.frame(cbind(bin_dmf, err_dmf))
  names(df) = c("phi", "err")
  df$x = hh$mids
  return(df)
}
devilsd10_noAGN_dmf_corr = foreach(i = 1:length(zmids)) %do% {
  message(i)
  vol = vol_mids[i] * 1.47 / (4*pi*(180/pi)^2)
  zidx = devilsd10_noAGN$z >= zbins[i] & devilsd10_noAGN$z < zbins[i+1]

  devils_alpha = Dale_L2M(devilsd10_noAGN$dustlum.total[zidx]/devilsd10_noAGN$dustmass.total[zidx])
  Mdust_new = devilsd10_noAGN$dustlum.total[zidx]/( Dale_M2L_variableDTH_func(devils_alpha, qPAH_VSG = qPAH_VSG) * qPAHZ_corr_func(10^devilsd10_noAGN$Zfinal[zidx], alpha = devils_alpha) )
  
  Mdust = Mdust_new * RR14_BPL(Z = 10^devilsd10_noAGN$Zfinal[zidx], doDTG = FALSE)/0.0073
  mc_err = 0.5*(mc_err_devilsd10_noAGN_dmf_corr[[i]][,3] - mc_err_devilsd10_noAGN_dmf_corr[[i]][,2])
  
  log_m = log10(Mdust[Mdust > 0])
  hh = maghist(x = log_m, breaks = sm_bins, plot = FALSE, verbose = FALSE)
  bin_dmf = hh$counts/(vol * diff(sm_bins))
  pois_dmf = sqrt(hh$counts)/(vol * diff(sm_bins))
  
  err_dmf = sqrt( sqrt(hh$counts)^2 + (err_floor*hh$counts)^2 + do_MC_err*mc_err^2 )/(vol * diff(sm_bins))
  df = data.frame(cbind(bin_dmf, err_dmf))
  names(df) = c("phi", "err")
  df$x = hh$mids
  return(df)
}
devilsd10_AGN_dmf_corr = foreach(i = 1:length(zmids)) %do% {
  message(i)
  vol = vol_mids[i] * 1.47 / (4*pi*(180/pi)^2)
  zidx = devilsd10_AGN$z >= zbins[i] & devilsd10_AGN$z < zbins[i+1]
  
  devils_alpha = Dale_L2M(devilsd10_AGN$dustlum.total[zidx]/devilsd10_AGN$dustmass.total[zidx])
  Mdust_new = devilsd10_AGN$dustlum.total[zidx]/( Dale_M2L_variableDTH_func(devils_alpha, qPAH_VSG = qPAH_VSG) * qPAHZ_corr_func(10^devilsd10_AGN$Zfinal[zidx], alpha = devils_alpha) )
  
  Mdust = Mdust_new * RR14_BPL(Z = 10^devilsd10_AGN$Zfinal[zidx], doDTG = FALSE)/0.0073
  mc_err = 0.5*(mc_err_devilsd10_AGN_dmf_corr[[i]][,3] - mc_err_devilsd10_AGN_dmf_corr[[i]][,2])
  
  log_m = log10(Mdust[Mdust > 0])
  hh = maghist(x = log_m, breaks = sm_bins, plot = FALSE, verbose = FALSE)
  bin_dmf = hh$counts/(vol * diff(sm_bins))
  pois_dmf = sqrt(hh$counts)/(vol * diff(sm_bins))
  
  err_dmf = sqrt( sqrt(hh$counts)^2 + (err_floor*hh$counts)^2 + do_MC_err*mc_err^2 )/(vol * diff(sm_bins))
  df = data.frame(cbind(bin_dmf, err_dmf))
  names(df) = c("phi", "err")
  df$x = hh$mids
  return(df)
}
devilsd10_hybrid_dmf_corr = foreach(i = 1:length(zmids)) %do% {
  message(i)
  vol = vol_mids[i] * 1.47 / (4*pi*(180/pi)^2)
  zidx = devilsd10_hybrid$z >= zbins[i] & devilsd10_hybrid$z < zbins[i+1]
  
  devils_alpha = Dale_L2M(devilsd10_hybrid$dustlum.total[zidx]/devilsd10_hybrid$dustmass.total[zidx])
  Mdust_new = devilsd10_hybrid$dustlum.total[zidx]/( Dale_M2L_variableDTH_func(devils_alpha, qPAH_VSG = qPAH_VSG) * qPAHZ_corr_func(10^devilsd10_hybrid$Zfinal[zidx], alpha = devils_alpha) )
  
  # Mdust_new_component = devilsd10_hybrid$dustlum.birth[zidx]/Dale_M2L_variableDTH_func(devilsd10_hybrid$alpha_SF_birth[zidx], qPAH_VSG = qPAH_VSG) + devilsd10_hybrid$dustlum.screen[zidx]/Dale_M2L_variableDTH_func(devilsd10_hybrid$alpha_SF_screen[zidx], qPAH_VSG = qPAH_VSG)
  ## I think Valid to assume this
  ## Compare using individual components, simplify things for GAMA
  # magplot(
  #   Mdust_new, Mdust_new_component, log = "xy"
  # )
  # maghist(
  #   Mdust_new/Mdust_new
  # )
  
  Mdust = Mdust_new * RR14_BPL(Z = 10^devilsd10_hybrid$Zfinal[zidx], doDTG = FALSE)/0.0073
  mc_err = 0.5*(mc_err_devilsd10_hybrid_dmf_corr[[i]][,3] - mc_err_devilsd10_hybrid_dmf_corr[[i]][,2])
  
  log_m = log10(Mdust[Mdust > 0])
  hh = maghist(x = log_m, breaks = sm_bins, plot = FALSE, verbose = FALSE)
  bin_dmf = hh$counts/(vol * diff(sm_bins))
  pois_dmf = sqrt(hh$counts)/(vol * diff(sm_bins))
  
  err_dmf = sqrt( sqrt(hh$counts)^2 + (err_floor*hh$counts)^2 + do_MC_err*mc_err^2 )/(vol * diff(sm_bins))
  df = data.frame(cbind(bin_dmf, err_dmf))
  names(df) = c("phi", "err")
  df$x = hh$mids
  return(df)
}
devilsd10_hybrid_dmf_corr_Zsol = foreach(i = 1:length(zmids)) %do% {
  message(i)
  vol = vol_mids[i] * 1.47 / (4*pi*(180/pi)^2)
  zidx = devilsd10_hybrid$z >= zbins[i] & devilsd10_hybrid$z < zbins[i+1]
  
  # Mdust_new_component = devilsd10_hybrid$dustlum.birth[zidx]/Dale_M2L_variableDTH_func(devilsd10_hybrid$alpha_SF_birth[zidx]) + devilsd10_hybrid$dustlum.screen[zidx]/Dale_M2L_variableDTH_func(devilsd10_hybrid$alpha_SF_screen[zidx])
  devils_alpha = Dale_L2M(devilsd10_hybrid$dustlum.total[zidx]/devilsd10_hybrid$dustmass.total[zidx])
  Mdust_new = devilsd10_hybrid$dustlum.total[zidx]/Dale_M2L_variableDTH_func(devils_alpha, qPAH_VSG = shivaei24_qPAHZ(0.014))
  
  ## I think Valid to assume this
  # magplot(
  #   Mdust_new, Mdust_new_component, log = "xy"
  # )
  # maghist(
  #   Mdust_new/Mdust_new_component 
  # )
  
  Mdust = Mdust_new * RR14_BPL(Z = 0.014, doDTG = FALSE)/0.0073
  mc_err = 0.5*(mc_err_devilsd10_hybrid_dmf_corr_Zsol[[i]][,3] - mc_err_devilsd10_hybrid_dmf_corr_Zsol[[i]][,2])
  
  log_m = log10(Mdust[Mdust > 0])
  hh = maghist(x = log_m, breaks = sm_bins, plot = FALSE, verbose = FALSE)
  bin_dmf = hh$counts/(vol * diff(sm_bins))
  pois_dmf = sqrt(hh$counts)/(vol * diff(sm_bins))
  
  err_dmf = sqrt( sqrt(hh$counts)^2 + (err_floor*hh$counts)^2 + do_MC_err*mc_err^2 )/(vol * diff(sm_bins))
  df = data.frame(cbind(bin_dmf, err_dmf))
  names(df) = c("phi", "err")
  df$x = hh$mids
  return(df)
}

gama_noAGN_gmf_corr = foreach(i = 1:length(zmids)) %do% {
  message(i)
  vol = vol_mids[i] * 217.54 / (4*pi*(180/pi)^2)
  zidx = gama_noAGN$z >= zbins[i] & gama_noAGN$z < zbins[i+1]
  
  gama_alpha = Dale_L2M(gama_noAGN$DustLum_50[zidx]/gama_noAGN$DustMass_50[zidx])
  Mdust_new = gama_noAGN$DustLum_50[zidx] / ( Dale_M2L_variableDTH_func(gama_alpha, qPAH_VSG = qPAH_VSG) * qPAHZ_corr_func(gama_noAGN$Zgas_50[zidx], alpha = gama_alpha) )
  
  Mdust = Mdust_new * RR14_BPL(Z = gama_noAGN$Zgas_50[zidx], doDTG = FALSE)/0.0073/RR14_BPL(Z = gama_noAGN$Zgas_50[zidx], doDTG = TRUE)
  mc_err = 0.5*(mc_err_gama_noAGN_gmf_corr[[i]][,3] - mc_err_gama_noAGN_gmf_corr[[i]][,2]) ## already factors in the error on the dust and the Z so that when we scale the dust the error is unchanged
  
  log_m = log10(Mdust[Mdust > 0])
  hh = maghist(x = log_m, breaks = sm_bins, plot = FALSE, verbose = FALSE)
  bin_dmf = hh$counts/(vol * diff(sm_bins))
  pois_dmf = sqrt(hh$counts)/(vol * diff(sm_bins))
  
  err_dmf = sqrt( sqrt(hh$counts)^2 + (err_floor*hh$counts)^2 + do_MC_err*mc_err^2 )/(vol * diff(sm_bins))
  df = data.frame(cbind(bin_dmf, err_dmf))
  names(df) = c("phi", "err")
  df$x = hh$mids
  return(df)
}
gama_AGN_gmf_corr = foreach(i = 1:length(zmids)) %do% {
  message(i)
  vol = vol_mids[i] * 217.54 / (4*pi*(180/pi)^2)
  zidx = gama_AGN$z >= zbins[i] & gama_AGN$z < zbins[i+1]
  
  gama_alpha = Dale_L2M(gama_AGN$DustLum_50[zidx]/gama_AGN$DustMass_50[zidx])
  Mdust_new = gama_AGN$DustLum_50[zidx]/ ( Dale_M2L_variableDTH_func(gama_alpha, qPAH_VSG = qPAH_VSG) * qPAHZ_corr_func(gama_AGN$Zgas_50[zidx], alpha = gama_alpha) )

  
  Mdust = Mdust_new * RR14_BPL(Z = gama_AGN$Zgas_50[zidx], doDTG = FALSE)/0.0073/RR14_BPL(Z = gama_AGN$Zgas_50[zidx], doDTG = TRUE)
  mc_err = 0.5*(mc_err_gama_AGN_gmf_corr[[i]][,3] - mc_err_gama_AGN_gmf_corr[[i]][,2]) ## already factors in the error on the dust and the Z so that when we scale the dust the error is unchanged
  
  log_m = log10(Mdust[Mdust > 0])
  hh = maghist(x = log_m, breaks = sm_bins, plot = FALSE, verbose = FALSE)
  bin_dmf = hh$counts/(vol * diff(sm_bins))
  pois_dmf = sqrt(hh$counts)/(vol * diff(sm_bins))
  
  err_dmf = sqrt( sqrt(hh$counts)^2 + (err_floor*hh$counts)^2 + do_MC_err*mc_err^2 )/(vol * diff(sm_bins))
  df = data.frame(cbind(bin_dmf, err_dmf))
  names(df) = c("phi", "err")
  df$x = hh$mids
  return(df)
}
gama_hybrid_gmf_corr = foreach(i = 1:length(zmids)) %do% {
  message(i)
  vol = vol_mids[i] * 217.54 / (4*pi*(180/pi)^2)
  zidx = gama_hybrid$z >= zbins[i] & gama_hybrid$z < zbins[i+1]
  
  gama_alpha = Dale_L2M(gama_hybrid$DustLum_50[zidx]/gama_hybrid$DustMass_50[zidx])
  Mdust_new = gama_hybrid$DustLum_50[zidx] / ( Dale_M2L_variableDTH_func(gama_alpha, qPAH_VSG = qPAH_VSG) * qPAHZ_corr_func(gama_hybrid$Zgas_50[zidx], alpha = gama_alpha) )
  
  Mdust = Mdust_new * RR14_BPL(Z = gama_hybrid$Zgas_50[zidx], doDTG = FALSE)/0.0073/RR14_BPL(Z = gama_hybrid$Zgas_50[zidx], doDTG = TRUE)
  mc_err = 0.5*(mc_err_gama_hybrid_gmf_corr[[i]][,3] - mc_err_gama_hybrid_gmf_corr[[i]][,2]) ## already factors in the error on the dust and the Z so that when we scale the dust the error is unchanged
  
  log_m = log10(Mdust[Mdust > 0])
  hh = maghist(x = log_m, breaks = sm_bins, plot = FALSE, verbose = FALSE)
  bin_dmf = hh$counts/(vol * diff(sm_bins))
  pois_dmf = sqrt(hh$counts)/(vol * diff(sm_bins))
  
  err_dmf = sqrt( sqrt(hh$counts)^2 + (err_floor*hh$counts)^2 + do_MC_err*mc_err^2 )/(vol * diff(sm_bins))
  df = data.frame(cbind(bin_dmf, err_dmf))
  names(df) = c("phi", "err")
  df$x = hh$mids
  return(df)
}
devilsd10_noAGN_gmf_corr = foreach(i = 1:length(zmids)) %do% {
  message(i)
  vol = vol_mids[i] * 1.47 / (4*pi*(180/pi)^2)
  zidx = devilsd10_noAGN$z >= zbins[i] & devilsd10_noAGN$z < zbins[i+1]
  
  devils_alpha = Dale_L2M(devilsd10_noAGN$dustlum.total[zidx]/devilsd10_noAGN$dustmass.total[zidx])
  Mdust_new = devilsd10_noAGN$dustlum.total[zidx] / ( Dale_M2L_variableDTH_func(gama_alpha, qPAH_VSG = qPAH_VSG) * qPAHZ_corr_func(10^devilsd10_noAGN$Zfinal[zidx], alpha = devils_alpha) )
  
  Mdust = Mdust_new * RR14_BPL(Z = 10^devilsd10_noAGN$Zfinal[zidx], doDTG = FALSE)/0.0073/RR14_BPL(Z = 10^devilsd10_noAGN$Zfinal[zidx], doDTG = TRUE)
  
  mc_err = 0.5*(mc_err_devilsd10_noAGN_gmf_corr[[i]][,3] - mc_err_devilsd10_noAGN_gmf_corr[[i]][,2])
  
  log_m = log10(Mdust[Mdust > 0])
  hh = maghist(x = log_m, breaks = sm_bins, plot = FALSE, verbose = FALSE)
  bin_dmf = hh$counts/(vol * diff(sm_bins))
  pois_dmf = sqrt(hh$counts)/(vol * diff(sm_bins))
  
  err_dmf = sqrt( sqrt(hh$counts)^2 + (err_floor*hh$counts)^2 + do_MC_err*mc_err^2 )/(vol * diff(sm_bins))
  df = data.frame(cbind(bin_dmf, err_dmf))
  names(df) = c("phi", "err")
  df$x = hh$mids
  return(df)
}
devilsd10_AGN_gmf_corr = foreach(i = 1:length(zmids)) %do% {
  message(i)
  vol = vol_mids[i] * 1.47 / (4*pi*(180/pi)^2)
  zidx = devilsd10_AGN$z >= zbins[i] & devilsd10_AGN$z < zbins[i+1]
  
  devils_alpha = Dale_L2M(devilsd10_AGN$dustlum.total[zidx]/devilsd10_AGN$dustmass.total[zidx])
  Mdust_new = devilsd10_AGN$dustlum.total[zidx] / ( Dale_M2L_variableDTH_func(gama_alpha, qPAH_VSG = qPAH_VSG) * qPAHZ_corr_func(10^devilsd10_AGN$Zfinal[zidx], alpha = devils_alpha) )
  
  Mdust = Mdust_new * RR14_BPL(Z = 10^devilsd10_AGN$Zfinal[zidx], doDTG = FALSE)/0.0073/RR14_BPL(Z = 10^devilsd10_AGN$Zfinal[zidx], doDTG = TRUE)
  
  mc_err = 0.5*(mc_err_devilsd10_AGN_gmf_corr[[i]][,3] - mc_err_devilsd10_AGN_gmf_corr[[i]][,2])
  
  log_m = log10(Mdust[Mdust > 0])
  hh = maghist(x = log_m, breaks = sm_bins, plot = FALSE, verbose = FALSE)
  bin_dmf = hh$counts/(vol * diff(sm_bins))
  pois_dmf = sqrt(hh$counts)/(vol * diff(sm_bins))
  
  err_dmf = sqrt( sqrt(hh$counts)^2 + (err_floor*hh$counts)^2 + do_MC_err*mc_err^2 )/(vol * diff(sm_bins))
  df = data.frame(cbind(bin_dmf, err_dmf))
  names(df) = c("phi", "err")
  df$x = hh$mids
  return(df)
}
devilsd10_hybrid_gmf_corr = foreach(i = 1:length(zmids)) %do% {
  message(i)
  vol = vol_mids[i] * 1.47 / (4*pi*(180/pi)^2)
  zidx = devilsd10_hybrid$z >= zbins[i] & devilsd10_hybrid$z < zbins[i+1]
  
  devils_alpha = Dale_L2M(devilsd10_hybrid$dustlum.total[zidx]/devilsd10_hybrid$dustmass.total[zidx])
  Mdust_new = devilsd10_hybrid$dustlum.total[zidx] / ( Dale_M2L_variableDTH_func(gama_alpha, qPAH_VSG = qPAH_VSG) * qPAHZ_corr_func(10^devilsd10_hybrid$Zfinal[zidx], alpha = devils_alpha) )

  Mdust = Mdust_new * RR14_BPL(Z = 10^devilsd10_hybrid$Zfinal[zidx], doDTG = FALSE)/0.0073/RR14_BPL(Z = 10^devilsd10_hybrid$Zfinal[zidx], doDTG = TRUE)
  
  mc_err = 0.5*(mc_err_devilsd10_hybrid_gmf_corr[[i]][,3] - mc_err_devilsd10_hybrid_gmf_corr[[i]][,2])
  
  log_m = log10(Mdust[Mdust > 0])
  hh = maghist(x = log_m, breaks = sm_bins, plot = FALSE, verbose = FALSE)
  bin_dmf = hh$counts/(vol * diff(sm_bins))
  pois_dmf = sqrt(hh$counts)/(vol * diff(sm_bins))
  
  err_dmf = sqrt( sqrt(hh$counts)^2 + (err_floor*hh$counts)^2 + do_MC_err*mc_err^2 )/(vol * diff(sm_bins))
  df = data.frame(cbind(bin_dmf, err_dmf))
  names(df) = c("phi", "err")
  df$x = hh$mids
  return(df)
}

gama_noAGN_smf = foreach(i = 1:length(zmids)) %do% {
  message(i)
  vol = vol_mids[i] * 217.54 / (4*pi*(180/pi)^2)
  zidx = gama_noAGN$z >= zbins[i] & gama_noAGN$z < zbins[i+1]
  Mdust = gama_noAGN$StellarMass_50[zidx]
  mc_err = 0.5*(mc_err_gama_noAGN_smf[[i]][,3] - mc_err_gama_noAGN_smf[[i]][,2])
  
  log_m = log10(Mdust[Mdust > 0])
  hh = maghist(x = log_m, breaks = sm_bins, plot = FALSE, verbose = FALSE)
  bin_dmf = hh$counts/(vol * diff(sm_bins))
  pois_dmf = sqrt(hh$counts)/(vol * diff(sm_bins))
  
  err_dmf = sqrt( sqrt(hh$counts)^2 + (err_floor*hh$counts)^2 + do_MC_err*mc_err^2 )/(vol * diff(sm_bins))
  df = data.frame(cbind(bin_dmf, err_dmf))
  names(df) = c("phi", "err")
  df$x = hh$mids
  return(df)
}
gama_AGN_smf = foreach(i = 1:length(zmids)) %do% {
  message(i)
  vol = vol_mids[i] * 217.54 / (4*pi*(180/pi)^2)
  zidx = gama_AGN$z >= zbins[i] & gama_AGN$z < zbins[i+1]
  Mdust = gama_AGN$StellarMass_50[zidx]
  mc_err = 0.5*(mc_err_gama_AGN_smf[[i]][,3] - mc_err_gama_AGN_smf[[i]][,2])
  
  log_m = log10(Mdust[Mdust > 0])
  hh = maghist(x = log_m, breaks = sm_bins, plot = FALSE, verbose = FALSE)
  bin_dmf = hh$counts/(vol * diff(sm_bins))
  pois_dmf = sqrt(hh$counts)/(vol * diff(sm_bins))
  
  err_dmf = sqrt( sqrt(hh$counts)^2 + (err_floor*hh$counts)^2 + do_MC_err*mc_err^2 )/(vol * diff(sm_bins))
  df = data.frame(cbind(bin_dmf, err_dmf))
  names(df) = c("phi", "err")
  df$x = hh$mids
  return(df)
}
gama_hybrid_smf = foreach(i = 1:length(zmids)) %do% {
  message(i)
  vol = vol_mids[i] * 217.54 / (4*pi*(180/pi)^2)
  zidx = gama_hybrid$z >= zbins[i] & gama_hybrid$z < zbins[i+1]
  Mdust = gama_hybrid$StellarMass_50[zidx]
  mc_err = 0.5*(mc_err_gama_hybrid_smf[[i]][,3] - mc_err_gama_hybrid_smf[[i]][,2])
  
  log_m = log10(Mdust[Mdust > 0])
  hh = maghist(x = log_m, breaks = sm_bins, plot = FALSE, verbose = FALSE)
  bin_dmf = hh$counts/(vol * diff(sm_bins))
  pois_dmf = sqrt(hh$counts)/(vol * diff(sm_bins))
  
  err_dmf = sqrt( sqrt(hh$counts)^2 + (err_floor*hh$counts)^2 + do_MC_err*mc_err^2 )/(vol * diff(sm_bins))
  df = data.frame(cbind(bin_dmf, err_dmf))
  names(df) = c("phi", "err")
  df$x = hh$mids
  return(df)
}
devilsd10_noAGN_smf = foreach(i = 1:length(zmids)) %do% {
  message(i)
  vol = vol_mids[i] * 1.47 / (4*pi*(180/pi)^2)
  zidx = devilsd10_noAGN$z >= zbins[i] & devilsd10_noAGN$z < zbins[i+1]
  Mdust = devilsd10_noAGN$StellarMass[zidx]
  mc_err = 0.5*(mc_err_devilsd10_noAGN_smf[[i]][,3] - mc_err_devilsd10_noAGN_smf[[i]][,2])
  
  log_m = log10(Mdust[Mdust > 0])
  hh = maghist(x = log_m, breaks = sm_bins, plot = FALSE, verbose = FALSE)
  bin_dmf = hh$counts/(vol * diff(sm_bins))
  pois_dmf = sqrt(hh$counts)/(vol * diff(sm_bins))
  
  err_dmf = sqrt( sqrt(hh$counts)^2 + (err_floor*hh$counts)^2 + do_MC_err*mc_err^2 )/(vol * diff(sm_bins))
  df = data.frame(cbind(bin_dmf, err_dmf))
  names(df) = c("phi", "err")
  df$x = hh$mids
  return(df)
}
devilsd10_AGN_smf = foreach(i = 1:length(zmids)) %do% {
  message(i)
  vol = vol_mids[i] * 1.47 / (4*pi*(180/pi)^2)
  zidx = devilsd10_AGN$z >= zbins[i] & devilsd10_AGN$z < zbins[i+1]
  Mdust = devilsd10_AGN$StellarMass[zidx]
  mc_err = 0.5*(mc_err_devilsd10_AGN_smf[[i]][,3] - mc_err_devilsd10_AGN_smf[[i]][,2])
  
  log_m = log10(Mdust[Mdust > 0])
  hh = maghist(x = log_m, breaks = sm_bins, plot = FALSE, verbose = FALSE)
  bin_dmf = hh$counts/(vol * diff(sm_bins))
  pois_dmf = sqrt(hh$counts)/(vol * diff(sm_bins))
  
  err_dmf = sqrt( sqrt(hh$counts)^2 + (err_floor*hh$counts)^2 + do_MC_err*mc_err^2 )/(vol * diff(sm_bins))
  df = data.frame(cbind(bin_dmf, err_dmf))
  names(df) = c("phi", "err")
  df$x = hh$mids
  return(df)
}
devilsd10_hybrid_smf = foreach(i = 1:length(zmids)) %do% {
  message(i)
  vol = vol_mids[i] * 1.47 / (4*pi*(180/pi)^2)
  zidx = devilsd10_hybrid$z >= zbins[i] & devilsd10_hybrid$z < zbins[i+1]
  Mdust = devilsd10_hybrid$StellarMass[zidx]
  mc_err = 0.5*(mc_err_devilsd10_hybrid_smf[[i]][,3] - mc_err_devilsd10_hybrid_smf[[i]][,2])
  
  log_m = log10(Mdust[Mdust > 0])
  hh = maghist(x = log_m, breaks = sm_bins, plot = FALSE, verbose = FALSE)
  bin_dmf = hh$counts/(vol * diff(sm_bins))
  pois_dmf = sqrt(hh$counts)/(vol * diff(sm_bins))
  
  err_dmf = sqrt( sqrt(hh$counts)^2 + (err_floor*hh$counts)^2 + do_MC_err*mc_err^2 )/(vol * diff(sm_bins))
  df = data.frame(cbind(bin_dmf, err_dmf))
  names(df) = c("phi", "err")
  df$x = hh$mids
  return(df)
}

## Fit functions 
## WHICH = 0L neither, 1L = DEVILS, 2L = GAMA
Nsamples = 10000
combine_noAGN_dmf = foreach(i = 1:length(zmids)) %do% {
  
  gama_dmf = gama_noAGN_dmf[[i]]
  gama_err = gama_dmf[,2]
  
  devilsd10_dmf = devilsd10_noAGN_dmf[[i]]
  devilsd10_err = devilsd10_dmf[,2]
  
  gama_xlim = sm_mids[which.max(gama_dmf[,1])]
  devilsd10_xlim = sm_mids[which.max(devilsd10_dmf[,1])]
  gama_SN = gama_dmf[,1]/gama_err
  devilsd10_SN = devilsd10_dmf[,1]/devilsd10_err
  
  gama_idx = sm_mids >= gama_xlim
  devilsd10_idx = sm_mids >= devilsd10_xlim
  gama_dmf[gama_idx, ]
  
  combine_dmf = foreach(j = 1:length(sm_mids), .combine = rbind) %do% {
    if( gama_err[j] <= ifelse(devilsd10_err[j]==0, 999, devilsd10_err[j]) ){ ## check if GAMA or DEVILS has better S/N
      if( !gama_idx[j] ){ ## If below GAMA completeness, use devils 
        c(devilsd10_dmf[j,1], devilsd10_err[j], 1L)
      }else{ 
        if( (devilsd10_dmf[j,1] - gama_dmf[j,1])/(gama_err[j] + 1e-323) > 5 ){ 
          c(devilsd10_dmf[j,1], devilsd10_err[j], 1L) ## If above GAMA completeness and devils is more than 5sigma above then use devils
        }else{
          c(gama_dmf[j,1], gama_err[j], 2L)
        }
      }
    }else{
      c(devilsd10_dmf[j,1], devilsd10_err[j], 1L)
    }
  }
  combine_dmf[sm_mids < devilsd10_xlim, ] = 0

  pinit = c(10, 0.0, 0.0, -4, -4)
  highout = Highlander(
    parm = pinit, 
    lower = c(5, -2.0, -1.1, -8, -8),
    upper = c(15, -0.8, 1.5, 1, 1),
    Data = list(
      xx = sm_mids[combine_dmf[,1] > 0], 
      yy = combine_dmf[,1][combine_dmf[,1] > 0], 
      yyerr = combine_dmf[,2][combine_dmf[,1] > 0], 
      func = double_schechter, 
      prior = function(p){
        sum(
          # dnorm(p[1], x = 8.0, sd = 4.0, log = TRUE),
          # dnorm(p[2], x = -1.0, sd = 1.0, log = TRUE),
          # dnorm(p[3], x = -1.0, sd = 1.0, log = TRUE),
          # dnorm(p[4], x = -2.0, sd = 4.0, log = TRUE),
          # dnorm(p[5], x = -3.0, sd = 4.0, log = TRUE),
          0
        )
      }
    ), 
    likefunc = LL, 
    liketype = "max", 
    Niters = c(Nsamples, Nsamples), 
    NfinalMCMC = Nsamples, 
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
  
  df_fit = data.frame(
    "x" = 10^mdustvec,
    "Q50" = q50_fit,
    "Q16" = q16_fit,
    "Q84" = q84_fit
  )
  df_fit$ERR = 0.5*(df_fit$Q84 - df_fit$Q16)
  
  df_cosmic_samples = foreach(k = 1:dim(highout$LD_last$Posterior1)[1], .combine = c) %do% {
    trapz(
      x = mdustvec,
      y = 10^mdustvec * fit_samples[k,]
    )
  }
  df_cosmic = data.frame(
    "Q50" = quantile(df_cosmic_samples, probs = 0.5),
    "Q16" = quantile(df_cosmic_samples, probs = 0.16),
    "Q84" = quantile(df_cosmic_samples, probs = 0.84)
  )
  df_cosmic$ERR = 0.5 * (df_cosmic$Q84 - df_cosmic$Q16)
  
  combine_dmf = data.frame(combine_dmf)
  names(combine_dmf) = c("phi", "err", "WHICH")
  combine_dmf$x = sm_mids
  ret_ = list(
    "dmf" = combine_dmf,
    "fit" = df_fit, 
    "cosmic" = df_cosmic,
    "highout" = highout
  )
  
  png(paste0("~/Documents/DustMassDensity/plots/dmf_noAGN/lbt", round(lbt_mids[i],3), ".png"))
  magplot(
    NA,
    log = "y",
    xlim = c(3.5, 13.5),
    ylim = c(1e-8, 10),
    xlab = "Mdust",
    ylab = "Phi",
    pch = 1,
    col = "red"
  )
  lines(
    mdustvec, q50_fit, col = "red"
  )
  lines(
    mdustvec, q16_fit, col = "red", lty = 2
  )
  lines(
    mdustvec, q84_fit, col = "red", lty = 2
  )
  points(
    sm_mids,
    combine_dmf[,1],
    pch = 16,
    cex = 2,
    col = alpha("red", 0.6)
  )
  magerr(
    sm_mids,
    combine_dmf[,1],
    ylo = combine_dmf[,2],
    col = "red"
  )
  points(
    sm_mids,
    gama_dmf[,1],
    pch = "x",
    col = "purple"
  )
  magerr(
    sm_mids,
    gama_dmf[,1],
    ylo = gama_err,
    col = "purple"
  )
  points(
    sm_mids,
    devilsd10_dmf[,1],
    pch = "x",
    col = "cornflowerblue"
  )
  magerr(
    sm_mids,
    devilsd10_dmf[,1],
    ylo = devilsd10_err,
    col = "cornflowerblue"
  )
  dev.off()
  return(ret_)
}
combine_AGN_dmf = foreach(i = 1:length(zmids)) %do% {
  
  gama_dmf = gama_AGN_dmf[[i]]
  gama_err = gama_dmf[,2]
  
  devilsd10_dmf = devilsd10_AGN_dmf[[i]]
  devilsd10_err = devilsd10_dmf[,2]
  
  gama_xlim = sm_mids[which.max(gama_dmf[,1])]
  devilsd10_xlim = sm_mids[which.max(devilsd10_dmf[,1])]
  gama_SN = gama_dmf[,1]/gama_err
  devilsd10_SN = devilsd10_dmf[,1]/devilsd10_err
  
  gama_idx = sm_mids >= gama_xlim
  devilsd10_idx = sm_mids >= devilsd10_xlim
  gama_dmf[gama_idx, ]
  
  combine_dmf = foreach(j = 1:length(sm_mids), .combine = rbind) %do% {
    if( gama_err[j] <= ifelse(devilsd10_err[j]==0, 999, devilsd10_err[j]) ){ ## check if GAMA or DEVILS has better S/N
      if( !gama_idx[j] ){ ## If below GAMA completeness, use devils 
        c(devilsd10_dmf[j,1], devilsd10_err[j])
      }else{ 
        if( (devilsd10_dmf[j,1] - gama_dmf[j,1])/(gama_err[j] + 1e-323) > 5 ){ 
          c(devilsd10_dmf[j,1], devilsd10_err[j]) ## If above GAMA completeness and devils is more than 5sigma above then use devils
        }else{
          c(gama_dmf[j,1], gama_err[j])
        }
      }
    }else{
      c(devilsd10_dmf[j,1], devilsd10_err[j])
    }
  }
  combine_dmf[sm_mids < devilsd10_xlim, ] = 0

  pinit = c(10, 0.0, 0.0, -4, -4)
  highout = Highlander(
    parm = pinit, 
    lower = c(5, -2.0, -1.1, -8, -8),
    upper = c(15, -0.8, 1.5, 1, 1),
    Data = list(
      xx = sm_mids[combine_dmf[,1] > 0], 
      yy = combine_dmf[,1][combine_dmf[,1] > 0], 
      yyerr = combine_dmf[,2][combine_dmf[,1] > 0], 
      func = double_schechter, 
      prior = function(p){
        sum(
          # dnorm(p[1], x = 8.0, sd = 4.0, log = TRUE),
          # dnorm(p[2], x = -1.0, sd = 1.0, log = TRUE),
          # dnorm(p[3], x = -1.0, sd = 1.0, log = TRUE),
          # dnorm(p[4], x = -2.0, sd = 4.0, log = TRUE),
          # dnorm(p[5], x = -3.0, sd = 4.0, log = TRUE),
          0
        )
      }
    ), 
    likefunc = LL, 
    liketype = "max", 
    Niters = c(Nsamples, Nsamples), 
    NfinalMCMC = Nsamples, 
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
  
  df_fit = data.frame(
    "x" = 10^mdustvec,
    "Q50" = q50_fit,
    "Q16" = q16_fit,
    "Q84" = q84_fit
  )
  df_fit$ERR = 0.5*(df_fit$Q84 - df_fit$Q16)
   
  df_cosmic_samples = foreach(k = 1:dim(highout$LD_last$Posterior1)[1], .combine = c) %do% {
    trapz(
      x = mdustvec,
      y = 10^mdustvec * fit_samples[k,]
    )
  }
  df_cosmic = data.frame(
    "Q50" = quantile(df_cosmic_samples, probs = 0.5),
    "Q16" = quantile(df_cosmic_samples, probs = 0.16),
    "Q84" = quantile(df_cosmic_samples, probs = 0.84)
  )
  df_cosmic$ERR = 0.5 * (df_cosmic$Q84 - df_cosmic$Q16)
  
  combine_dmf = data.frame(combine_dmf)
  names(combine_dmf) = c("phi", "err")
  combine_dmf$x = sm_mids
  
  ret_ = list(
    "dmf" = combine_dmf,
    "fit" = df_fit, 
    "cosmic" = df_cosmic,
    "highout" = highout
  )
  png(paste0("~/Documents/DustMassDensity/plots/dmf_AGN/lbt", round(lbt_mids[i],3), ".png"))
  magplot(
    NA,
    log = "y",
    xlim = c(3.5, 13.5),
    ylim = c(1e-8, 10),
    xlab = "Mdust",
    ylab = "Phi",
    pch = 1,
    col = "red"
  )
  lines(
    mdustvec, q50_fit, col = "red"
  )
  lines(
    mdustvec, q16_fit, col = "red", lty = 2
  )
  lines(
    mdustvec, q84_fit, col = "red", lty = 2
  )
  points(
    sm_mids,
    combine_dmf[,1],
    pch = 16,
    cex = 2,
    col = alpha("red", 0.6)
  )
  magerr(
    sm_mids,
    combine_dmf[,1],
    ylo = combine_dmf[,2],
    col = "red"
  )
  points(
    sm_mids,
    gama_dmf[,1],
    pch = "x",
    col = "purple"
  )
  magerr(
    sm_mids,
    gama_dmf[,1],
    ylo = gama_err,
    col = "purple"
  )
  points(
    sm_mids,
    devilsd10_dmf[,1],
    pch = "x",
    col = "cornflowerblue"
  )
  magerr(
    sm_mids,
    devilsd10_dmf[,1],
    ylo = devilsd10_err,
    col = "cornflowerblue"
  )
  dev.off()
  return(ret_)
}
combine_hybrid_dmf = foreach(i = 1:length(zmids)) %do% {
  
  gama_dmf = gama_hybrid_dmf[[i]]
  gama_err = gama_dmf[,2]
  
  devilsd10_dmf = devilsd10_hybrid_dmf[[i]]
  devilsd10_err = devilsd10_dmf[,2]
  
  gama_xlim = sm_mids[which.max(gama_dmf[,1])]
  devilsd10_xlim = sm_mids[which.max(devilsd10_dmf[,1])]
  gama_SN = gama_dmf[,1]/gama_err
  devilsd10_SN = devilsd10_dmf[,1]/devilsd10_err
  
  gama_idx = sm_mids >= gama_xlim
  devilsd10_idx = sm_mids >= devilsd10_xlim
  gama_dmf[gama_idx, ]
  
  combine_dmf = foreach(j = 1:length(sm_mids), .combine = rbind) %do% {
    if( gama_err[j] <= ifelse(devilsd10_err[j]==0, 999, devilsd10_err[j]) ){ ## check if GAMA or DEVILS has better S/N
      if( !gama_idx[j] ){ ## If below GAMA completeness, use devils 
        c(devilsd10_dmf[j,1], devilsd10_err[j], 1L)
      }else{ 
        if( (devilsd10_dmf[j,1] - gama_dmf[j,1])/(gama_err[j] + 1e-323) > 5 ){ 
          c(devilsd10_dmf[j,1], devilsd10_err[j], 1L) ## If above GAMA completeness and devils is more than 5sigma above then use devils
        }else{
          c(gama_dmf[j,1], gama_err[j], 2L)
        }
      }
    }else{
      c(devilsd10_dmf[j,1], devilsd10_err[j], 1L)
    }
  }
  combine_dmf[sm_mids < devilsd10_xlim, ] = 0

  pinit = c(10, 0.0, 0.0, -4, -4)
  highout = Highlander(
    parm = pinit, 
    lower = c(5, -2.0, -1.1, -8, -8),
    upper = c(15, -0.8, 1.5, 1, 1),
    Data = list(
      xx = sm_mids[combine_dmf[,1] > 0], 
      yy = combine_dmf[,1][combine_dmf[,1] > 0], 
      yyerr = combine_dmf[,2][combine_dmf[,1] > 0], 
      func = double_schechter, 
      prior = function(p){
        sum(
          # dnorm(p[1], x = 8.0, sd = 4.0, log = TRUE),
          # dnorm(p[2], x = -1.0, sd = 1.0, log = TRUE),
          # dnorm(p[3], x = -1.0, sd = 1.0, log = TRUE),
          # dnorm(p[4], x = -2.0, sd = 4.0, log = TRUE),
          # dnorm(p[5], x = -3.0, sd = 4.0, log = TRUE),
          0
        )
      }
    ), 
    likefunc = LL, 
    liketype = "max", 
    Niters = c(Nsamples, Nsamples), 
    NfinalMCMC = Nsamples, 
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
  
  df_fit = data.frame(
    "x" = 10^mdustvec,
    "Q50" = q50_fit,
    "Q16" = q16_fit,
    "Q84" = q84_fit
  )
  df_fit$ERR = 0.5*(df_fit$Q84 - df_fit$Q16)
  df_cosmic_samples = foreach(k = 1:dim(highout$LD_last$Posterior1)[1], .combine = c) %do% {
    trapz(
      x = mdustvec,
      y = 10^mdustvec * fit_samples[k,]
    )
  }
  df_cosmic = data.frame(
    "Q50" = quantile(df_cosmic_samples, probs = 0.5),
    "Q16" = quantile(df_cosmic_samples, probs = 0.16),
    "Q84" = quantile(df_cosmic_samples, probs = 0.84)
  )
  df_cosmic$ERR = 0.5 * (df_cosmic$Q84 - df_cosmic$Q16)
  
  combine_dmf = data.frame(combine_dmf)
  names(combine_dmf) = c("phi", "err", "WHICH")
  combine_dmf$x = sm_mids
  
  ret_ = list(
    "dmf" = combine_dmf,
    "fit" = df_fit, 
    "cosmic" = df_cosmic,
    "highout" = highout
  )
  
  png(paste0("~/Documents/DustMassDensity/plots/dmf_hybrid/lbt_", round(lbt_mids[i],3), ".png"))
  magplot(
    NA,
    log = "y",
    xlim = c(3.5, 13.5),
    ylim = c(1e-8, 10),
    xlab = "Mdust",
    ylab = "Phi",
    pch = 1,
    col = "red"
  )
  lines(
    mdustvec, q50_fit, col = "red"
  )
  lines(
    mdustvec, q16_fit, col = "red", lty = 2
  )
  lines(
    mdustvec, q84_fit, col = "red", lty = 2
  )
  points(
    sm_mids,
    combine_dmf[,1],
    pch = 16,
    cex = 2,
    col = alpha("red", 0.6)
  )
  magerr(
    sm_mids,
    combine_dmf[,1],
    ylo = combine_dmf[,2],
    col = "red"
  )
  points(
    sm_mids,
    gama_dmf[,1],
    pch = "x",
    col = "purple"
  )
  magerr(
    sm_mids,
    gama_dmf[,1],
    ylo = gama_err,
    col = "purple"
  )
  points(
    sm_mids,
    devilsd10_dmf[,1],
    pch = "x",
    col = "cornflowerblue"
  )
  magerr(
    sm_mids,
    devilsd10_dmf[,1],
    ylo = devilsd10_err,
    col = "cornflowerblue"
  )
  dev.off()
  return(ret_)
}

combine_noAGN_dmf_corr = foreach(i = 1:length(zmids)) %do% {
  
  gama_dmf = gama_noAGN_dmf_corr[[i]]
  gama_err = gama_dmf[,2]
  
  devilsd10_dmf = devilsd10_noAGN_dmf_corr[[i]]
  devilsd10_err = devilsd10_dmf[,2]
  
  gama_xlim = sm_mids[which.max(gama_dmf[,1])]
  devilsd10_xlim = sm_mids[which.max(devilsd10_dmf[,1])]
  gama_SN = gama_dmf[,1]/gama_err
  
  gama_idx = sm_mids >= gama_xlim
  devilsd10_idx = sm_mids >= devilsd10_xlim
  gama_dmf[gama_idx, ]
  
  combine_dmf = foreach(j = 1:length(sm_mids), .combine = rbind) %do% {
    if( gama_err[j] <= ifelse(devilsd10_err[j]==0, 999, devilsd10_err[j]) ){ ## check if GAMA or DEVILS has better S/N
      if( !gama_idx[j] ){ ## If below GAMA completeness, use devils 
        c(devilsd10_dmf[j,1], devilsd10_err[j], 1L)
      }else{ 
        if( (devilsd10_dmf[j,1] - gama_dmf[j,1])/(gama_err[j] + 1e-323) > 5 ){ 
          c(devilsd10_dmf[j,1], devilsd10_err[j], 1L) ## If above GAMA completeness and devils is more than 5sigma above then use devils
        }else{
          c(gama_dmf[j,1], gama_err[j], 2L)
        }
      }
    }else{
      c(devilsd10_dmf[j,1], devilsd10_err[j], 1L)
    }
  }
  combine_dmf[sm_mids < devilsd10_xlim, ] = 0
  
  pinit = c(6, 0.0, 0.0, -4, -4)
  highout = Highlander(
    parm = pinit, 
    lower = c(5, -2.0, -1.1, -8, -8),
    upper = c(15, -0.8, 1.5, 1, 1),
    Data = list(
      xx = sm_mids[combine_dmf[,1] > 0], 
      yy = combine_dmf[,1][combine_dmf[,1] > 0], 
      yyerr = combine_dmf[,2][combine_dmf[,1] > 0], 
      func = double_schechter, 
      prior = function(p){
        sum(
          # dnorm(p[1], x = 8.0, sd = 4.0, log = TRUE),
          # dnorm(p[2], x = -1.0, sd = 1.0, log = TRUE),
          # dnorm(p[3], x = -1.0, sd = 1.0, log = TRUE),
          # dnorm(p[4], x = -2.0, sd = 4.0, log = TRUE),
          # dnorm(p[5], x = -3.0, sd = 4.0, log = TRUE),
          0
        )
      }
    ), 
    likefunc = LL, 
    liketype = "max", 
    Niters = c(Nsamples, Nsamples), 
    NfinalMCMC = Nsamples, 
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
  
  df_fit = data.frame(
    "x" = 10^mdustvec,
    "Q50" = q50_fit,
    "Q16" = q16_fit,
    "Q84" = q84_fit
  )
  df_fit$ERR = 0.5*(df_fit$Q84 - df_fit$Q16)
  
  df_cosmic_samples = foreach(k = 1:dim(highout$LD_last$Posterior1)[1], .combine = c) %do% {
    trapz(
      x = mdustvec,
      y = 10^mdustvec * fit_samples[k,]
    )
  }
  df_cosmic = data.frame(
    "Q50" = quantile(df_cosmic_samples, probs = 0.5),
    "Q16" = quantile(df_cosmic_samples, probs = 0.16),
    "Q84" = quantile(df_cosmic_samples, probs = 0.84)
  )
  df_cosmic$ERR = 0.5*(df_cosmic$Q84 - df_cosmic$Q16)
  
  combine_dmf = data.frame(combine_dmf)
  names(combine_dmf) = c("phi", "err", "WHICH")
  combine_dmf$x = sm_mids
  
  ret_ = list(
    "dmf" = combine_dmf,
    "fit" = df_fit, 
    "cosmic" = df_cosmic,
    "highout" = highout
  )
  
  png(paste0("~/Documents/DustMassDensity/plots/dmf_noAGN_corr/lbt_", round(lbt_mids[i],3), ".png"))
  magplot(
    NA,
    log = "y",
    xlim = c(3.5, 13.5),
    ylim = c(1e-8, 10),
    xlab = "Mdust",
    ylab = "Phi",
    pch = 1,
    col = "red"
  )
  lines(
    mdustvec, q50_fit, col = "red"
  )
  lines(
    mdustvec, q16_fit, col = "red", lty = 2
  )
  lines(
    mdustvec, q84_fit, col = "red", lty = 2
  )
  points(
    sm_mids,
    combine_dmf[,1],
    pch = 16,
    cex = 2,
    col = alpha("red", 0.6)
  )
  magerr(
    sm_mids,
    combine_dmf[,1],
    ylo = combine_dmf[,2],
    col = "red"
  )
  points(
    sm_mids,
    gama_dmf[,1],
    pch = "x",
    col = "purple"
  )
  magerr(
    sm_mids,
    gama_dmf[,1],
    ylo = gama_err,
    col = "purple"
  )
  points(
    sm_mids,
    devilsd10_dmf[,1],
    pch = "x",
    col = "cornflowerblue"
  )
  magerr(
    sm_mids,
    devilsd10_dmf[,1],
    ylo = devilsd10_err,
    col = "cornflowerblue"
  )
  dev.off()
  return(ret_)
}
combine_AGN_dmf_corr = foreach(i = 1:length(zmids)) %do% {
  
  gama_dmf = gama_AGN_dmf_corr[[i]]
  gama_err = gama_dmf[,2]
  
  devilsd10_dmf = devilsd10_AGN_dmf_corr[[i]]
  devilsd10_err = devilsd10_dmf[,2]
  
  gama_xlim = sm_mids[which.max(gama_dmf[,1])]
  devilsd10_xlim = sm_mids[which.max(devilsd10_dmf[,1])]
  gama_SN = gama_dmf[,1]/gama_err
  
  gama_idx = sm_mids >= gama_xlim
  devilsd10_idx = sm_mids >= devilsd10_xlim
  gama_dmf[gama_idx, ]
  
  combine_dmf = foreach(j = 1:length(sm_mids), .combine = rbind) %do% {
    if( gama_err[j] <= ifelse(devilsd10_err[j]==0, 999, devilsd10_err[j]) ){ ## check if GAMA or DEVILS has better S/N
      if( !gama_idx[j] ){ ## If below GAMA completeness, use devils 
        c(devilsd10_dmf[j,1], devilsd10_err[j], 1L)
      }else{ 
        if( (devilsd10_dmf[j,1] - gama_dmf[j,1])/(gama_err[j] + 1e-323) > 5 ){ 
          c(devilsd10_dmf[j,1], devilsd10_err[j], 1L) ## If above GAMA completeness and devils is more than 5sigma above then use devils
        }else{
          c(gama_dmf[j,1], gama_err[j], 2L)
        }
      }
    }else{
      c(devilsd10_dmf[j,1], devilsd10_err[j], 1L)
    }
  }
  combine_dmf[sm_mids < devilsd10_xlim, ] = 0
  
  pinit = c(6, 0.0, 0.0, -4, -4)
  highout = Highlander(
    parm = pinit, 
    lower = c(5, -2.0, -1.1, -8, -8),
    upper = c(15, -0.8, 1.5, 1, 1),
    Data = list(
      xx = sm_mids[combine_dmf[,1] > 0], 
      yy = combine_dmf[,1][combine_dmf[,1] > 0], 
      yyerr = combine_dmf[,2][combine_dmf[,1] > 0], 
      func = double_schechter, 
      prior = function(p){
        sum(
          # dnorm(p[1], x = 8.0, sd = 4.0, log = TRUE),
          # dnorm(p[2], x = -1.0, sd = 1.0, log = TRUE),
          # dnorm(p[3], x = -1.0, sd = 1.0, log = TRUE),
          # dnorm(p[4], x = -2.0, sd = 4.0, log = TRUE),
          # dnorm(p[5], x = -3.0, sd = 4.0, log = TRUE),
          0
        )
      }
    ), 
    likefunc = LL, 
    liketype = "max", 
    Niters = c(Nsamples, Nsamples), 
    NfinalMCMC = Nsamples, 
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
  
  df_fit = data.frame(
    "x" = 10^mdustvec,
    "Q50" = q50_fit,
    "Q16" = q16_fit,
    "Q84" = q84_fit
  )
  df_fit$ERR = 0.5*(df_fit$Q84 - df_fit$Q16)
  
  df_cosmic_samples = foreach(k = 1:dim(highout$LD_last$Posterior1)[1], .combine = c) %do% {
    trapz(
      x = mdustvec,
      y = 10^mdustvec * fit_samples[k,]
    )
  }
  df_cosmic = data.frame(
    "Q50" = quantile(df_cosmic_samples, probs = 0.5),
    "Q16" = quantile(df_cosmic_samples, probs = 0.16),
    "Q84" = quantile(df_cosmic_samples, probs = 0.84)
  )
  df_cosmic$ERR = 0.5*(df_cosmic$Q84 - df_cosmic$Q16)
  
  combine_dmf = data.frame(combine_dmf)
  names(combine_dmf) = c("phi", "err", "WHICH")
  combine_dmf$x = sm_mids
  
  ret_ = list(
    "dmf" = combine_dmf,
    "fit" = df_fit, 
    "cosmic" = df_cosmic,
    "highout" = highout
  )
  
  png(paste0("~/Documents/DustMassDensity/plots/dmf_AGN_corr/lbt_", round(lbt_mids[i],3), ".png"))
  magplot(
    NA,
    log = "y",
    xlim = c(3.5, 13.5),
    ylim = c(1e-8, 10),
    xlab = "Mdust",
    ylab = "Phi",
    pch = 1,
    col = "red"
  )
  lines(
    mdustvec, q50_fit, col = "red"
  )
  lines(
    mdustvec, q16_fit, col = "red", lty = 2
  )
  lines(
    mdustvec, q84_fit, col = "red", lty = 2
  )
  points(
    sm_mids,
    combine_dmf[,1],
    pch = 16,
    cex = 2,
    col = alpha("red", 0.6)
  )
  magerr(
    sm_mids,
    combine_dmf[,1],
    ylo = combine_dmf[,2],
    col = "red"
  )
  points(
    sm_mids,
    gama_dmf[,1],
    pch = "x",
    col = "purple"
  )
  magerr(
    sm_mids,
    gama_dmf[,1],
    ylo = gama_err,
    col = "purple"
  )
  points(
    sm_mids,
    devilsd10_dmf[,1],
    pch = "x",
    col = "cornflowerblue"
  )
  magerr(
    sm_mids,
    devilsd10_dmf[,1],
    ylo = devilsd10_err,
    col = "cornflowerblue"
  )
  dev.off()
  return(ret_)
}
combine_hybrid_dmf_corr = foreach(i = 1:length(zmids)) %do% {
  
  gama_dmf = gama_hybrid_dmf_corr[[i]]
  gama_err = gama_dmf[,2]
  
  devilsd10_dmf = devilsd10_hybrid_dmf_corr[[i]]
  devilsd10_err = devilsd10_dmf[,2]
  
  gama_xlim = sm_mids[which.max(gama_dmf[,1])]
  devilsd10_xlim = sm_mids[which.max(devilsd10_dmf[,1])]
  gama_SN = gama_dmf[,1]/gama_err
  
  gama_idx = sm_mids >= gama_xlim 
  devilsd10_idx = sm_mids >= devilsd10_xlim 
  gama_dmf[gama_idx, ]
  
  combine_dmf = foreach(j = 1:length(sm_mids), .combine = rbind) %do% {
    if( gama_err[j] <= ifelse(devilsd10_err[j]==0, 999, devilsd10_err[j]) ){ ## check if GAMA or DEVILS has better S/N
      if( !gama_idx[j] ){ ## If below GAMA completeness, use devils 
        c(devilsd10_dmf[j,1], devilsd10_err[j], 1L) ## include check on what data set
      }else{ 
        if( (devilsd10_dmf[j,1] - gama_dmf[j,1])/(gama_err[j] + 1e-323) > 5 ){ 
          c(devilsd10_dmf[j,1], devilsd10_err[j], 1L) ## If above GAMA completeness and devils is more than 5sigma above then use devils
        }else{
          c(gama_dmf[j,1], gama_err[j], 2L)
        }
      }
    }else{
      c(devilsd10_dmf[j,1], devilsd10_err[j], 1L)
    }
  }
  combine_dmf[sm_mids < devilsd10_xlim, ] = 0

  pinit = c(10, 0.0, 0.0, -4, -4)
  highout = Highlander(
    parm = pinit, 
    lower = c(5, -2.0, -1.1, -8, -8),
    upper = c(15, -0.8, 1.5, 1, 1),
    Data = list(
      xx = sm_mids[combine_dmf[,1] > 0], 
      yy = combine_dmf[,1][combine_dmf[,1] > 0], 
      yyerr = combine_dmf[,2][combine_dmf[,1] > 0], 
      func = double_schechter, 
      prior = function(p){
        sum(
          # dnorm(p[1], x = 8.0, sd = 4.0, log = TRUE),
          # dnorm(p[2], x = -1.0, sd = 1.0, log = TRUE),
          # dnorm(p[3], x = -1.0, sd = 1.0, log = TRUE),
          # dnorm(p[4], x = -2.0, sd = 4.0, log = TRUE),
          # dnorm(p[5], x = -3.0, sd = 4.0, log = TRUE),
          0
        )
      }
    ), 
    likefunc = LL, 
    liketype = "max", 
    Niters = c(Nsamples, Nsamples), 
    NfinalMCMC = Nsamples, 
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
  
  df_fit = data.frame(
    "x" = 10^mdustvec,
    "Q50" = q50_fit,
    "Q16" = q16_fit,
    "Q84" = q84_fit
  )
  df_fit$ERR = 0.5*(df_fit$Q84 - df_fit$Q16)
  
  df_cosmic_samples = foreach(k = 1:dim(highout$LD_last$Posterior1)[1], .combine = c) %do% {
    trapz(
      x = mdustvec,
      y = 10^mdustvec * fit_samples[k,]
    )
  }
  df_cosmic = data.frame(
    "Q50" = quantile(df_cosmic_samples, probs = 0.5),
    "Q16" = quantile(df_cosmic_samples, probs = 0.16),
    "Q84" = quantile(df_cosmic_samples, probs = 0.84)
  )
  df_cosmic$ERR = 0.5*(df_cosmic$Q84 - df_cosmic$Q16)
  
  combine_dmf = data.frame(combine_dmf)
  names(combine_dmf) = c("phi", "err", "WHICH")
  combine_dmf$x = sm_mids
  
  ret_ = list(
    "dmf" = combine_dmf,
    "fit" = df_fit, 
    "cosmic" = df_cosmic,
    "highout" = highout
  )
  
  png(paste0("~/Documents/DustMassDensity/plots/dmf_hybrid_corr/lbt_", round(lbt_mids[i],3), ".png"))
  magplot(
    NA,
    log = "y",
    xlim = c(3.5, 13.5),
    ylim = c(1e-8, 10),
    xlab = "Mdust",
    ylab = "Phi",
    pch = 1,
    col = "red"
  )
  lines(
    mdustvec, q50_fit, col = "red"
  )
  lines(
    mdustvec, q16_fit, col = "red", lty = 2
  )
  lines(
    mdustvec, q84_fit, col = "red", lty = 2
  )
  points(
    sm_mids,
    combine_dmf[,1],
    pch = 16,
    cex = 2,
    col = alpha("red", 0.6)
  )
  magerr(
    sm_mids,
    combine_dmf[,1],
    ylo = combine_dmf[,2],
    col = "red"
  )
  points(
    sm_mids,
    gama_dmf[,1],
    pch = "x",
    col = "purple"
  )
  magerr(
    sm_mids,
    gama_dmf[,1],
    ylo = gama_err,
    col = "purple"
  )
  points(
    sm_mids,
    devilsd10_dmf[,1],
    pch = "x",
    col = "cornflowerblue"
  )
  magerr(
    sm_mids,
    devilsd10_dmf[,1],
    ylo = devilsd10_err,
    col = "cornflowerblue"
  )
  dev.off()
  return(ret_)
}
combine_hybrid_dmf_corr_Zsol = foreach(i = 1:length(zmids)) %do% {
  
  gama_dmf = gama_hybrid_dmf_corr_Zsol[[i]]
  gama_err = gama_dmf[,2]
  
  devilsd10_dmf = devilsd10_hybrid_dmf_corr_Zsol[[i]]
  devilsd10_err = devilsd10_dmf[,2]
  
  gama_xlim = sm_mids[which.max(gama_dmf[,1])]
  devilsd10_xlim = sm_mids[which.max(devilsd10_dmf[,1])]
  gama_SN = gama_dmf[,1]/gama_err
  
  gama_idx = sm_mids >= gama_xlim
  devilsd10_idx = sm_mids >= devilsd10_xlim
  gama_dmf[gama_idx, ]
  
  combine_dmf = foreach(j = 1:length(sm_mids), .combine = rbind) %do% {
    if( gama_err[j] <= ifelse(devilsd10_err[j]==0, 999, devilsd10_err[j]) ){ ## check if GAMA or DEVILS has better S/N
      if( !gama_idx[j] ){ ## If below GAMA completeness, use devils 
        c(devilsd10_dmf[j,1], devilsd10_err[j], 1L) ## include check on what data set
      }else{ 
        if( (devilsd10_dmf[j,1] - gama_dmf[j,1])/(gama_err[j] + 1e-323) > 5 ){ 
          c(devilsd10_dmf[j,1], devilsd10_err[j], 1L) ## If above GAMA completeness and devils is more than 5sigma above then use devils
        }else{
          c(gama_dmf[j,1], gama_err[j], 2L)
        }
      }
    }else{
      c(devilsd10_dmf[j,1], devilsd10_err[j], 1L)
    }
  }
  combine_dmf[sm_mids < devilsd10_xlim, ] = 0
  
  pinit = c(6, 0.0, 0.0, -4, -4)
  highout = Highlander(
    parm = pinit, 
    lower = c(5, -2.0, -1.1, -8, -8),
    upper = c(15, -0.8, 1.5, 1, 1),
    Data = list(
      xx = sm_mids[combine_dmf[,1] > 0], 
      yy = combine_dmf[,1][combine_dmf[,1] > 0], 
      yyerr = combine_dmf[,2][combine_dmf[,1] > 0], 
      func = double_schechter, 
      prior = function(p){
        sum(
          # dnorm(p[1], x = 8.0, sd = 4.0, log = TRUE),
          # dnorm(p[2], x = -1.0, sd = 1.0, log = TRUE),
          # dnorm(p[3], x = -1.0, sd = 1.0, log = TRUE),
          # dnorm(p[4], x = -2.0, sd = 4.0, log = TRUE),
          # dnorm(p[5], x = -3.0, sd = 4.0, log = TRUE),
          0
        )
      }
    ), 
    likefunc = LL, 
    liketype = "max", 
    Niters = c(Nsamples, Nsamples), 
    NfinalMCMC = Nsamples, 
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
  
  df_fit = data.frame(
    "x" = 10^mdustvec,
    "Q50" = q50_fit,
    "Q16" = q16_fit,
    "Q84" = q84_fit
  )
  df_fit$ERR = 0.5*(df_fit$Q84 - df_fit$Q16)
  
  df_cosmic_samples = foreach(k = 1:dim(highout$LD_last$Posterior1)[1], .combine = c) %do% {
    trapz(
      x = mdustvec,
      y = 10^mdustvec * fit_samples[k,]
    )
  }
  df_cosmic = data.frame(
    "Q50" = quantile(df_cosmic_samples, probs = 0.5),
    "Q16" = quantile(df_cosmic_samples, probs = 0.16),
    "Q84" = quantile(df_cosmic_samples, probs = 0.84)
  )
  df_cosmic$ERR = 0.5*(df_cosmic$Q84 - df_cosmic$Q16)
  
  combine_dmf = data.frame(combine_dmf)
  names(combine_dmf) = c("phi", "err", "WHICH")
  combine_dmf$x = sm_mids
  
  ret_ = list(
    "dmf" = combine_dmf,
    "fit" = df_fit, 
    "cosmic" = df_cosmic,
    "highout" = highout
  )
  
  png(paste0("~/Documents/DustMassDensity/plots/dmf_hybrid_corr_Zsol/lbt_", round(lbt_mids[i],3), ".png"))
  magplot(
    NA,
    log = "y",
    xlim = c(3.5, 13.5),
    ylim = c(1e-8, 10),
    xlab = "Mdust",
    ylab = "Phi",
    pch = 1,
    col = "red"
  )
  lines(
    mdustvec, q50_fit, col = "red"
  )
  lines(
    mdustvec, q16_fit, col = "red", lty = 2
  )
  lines(
    mdustvec, q84_fit, col = "red", lty = 2
  )
  points(
    sm_mids,
    combine_dmf[,1],
    pch = 16,
    cex = 2,
    col = alpha("red", 0.6)
  )
  magerr(
    sm_mids,
    combine_dmf[,1],
    ylo = combine_dmf[,2],
    col = "red"
  )
  points(
    sm_mids,
    gama_dmf[,1],
    pch = "x",
    col = "purple"
  )
  magerr(
    sm_mids,
    gama_dmf[,1],
    ylo = gama_err,
    col = "purple"
  )
  points(
    sm_mids,
    devilsd10_dmf[,1],
    pch = "x",
    col = "cornflowerblue"
  )
  magerr(
    sm_mids,
    devilsd10_dmf[,1],
    ylo = devilsd10_err,
    col = "cornflowerblue"
  )
  dev.off()
  return(ret_)
}

combine_noAGN_gmf_corr = foreach(i = 1:length(zmids)) %do% {
  
  gama_dmf = gama_noAGN_gmf_corr[[i]]
  gama_err = gama_dmf[,2]
  
  devilsd10_dmf = devilsd10_noAGN_gmf_corr[[i]]
  devilsd10_err = devilsd10_dmf[,2]
  
  gama_xlim = sm_mids[which.max(gama_dmf[,1])]
  devilsd10_xlim = sm_mids[which.max(devilsd10_dmf[,1])]
  gama_SN = gama_dmf[,1]/gama_err
  devilsd10_SN = devilsd10_dmf[,1]/devilsd10_err
  
  gama_idx = sm_mids >= gama_xlim
  devilsd10_idx = sm_mids >= devilsd10_xlim
  gama_dmf[gama_idx, ]
  
  combine_dmf = foreach(j = 1:length(sm_mids), .combine = rbind) %do% {
    if( gama_err[j] <= ifelse(devilsd10_err[j]==0, 999, devilsd10_err[j]) ){ ## check if GAMA or DEVILS has better S/N
      if( !gama_idx[j] ){ ## If below GAMA completeness, use devils 
        c(devilsd10_dmf[j,1], devilsd10_err[j], 1L)
      }else{ 
        if( (devilsd10_dmf[j,1] - gama_dmf[j,1])/(gama_err[j] + 1e-323) > 5 ){ 
          c(devilsd10_dmf[j,1], devilsd10_err[j], 1L) ## If above GAMA completeness and devils is more than 5sigma above then use devils
        }else{
          c(gama_dmf[j,1], gama_err[j], 2L)
        }
      }
    }else{
      c(devilsd10_dmf[j,1], devilsd10_err[j], 1L)
    }
  }
  combine_dmf[sm_mids < devilsd10_xlim, ] = 0
  
  pinit = c(10, 0.0, 0.0, -4, -4)
  highout = Highlander(
    parm = pinit, 
    lower = c(5, -2.0, -1.1, -8, -8),
    upper = c(15, -0.8, 1.5, 1, 1),
    Data = list(
      xx = sm_mids[combine_dmf[,1] > 0], 
      yy = combine_dmf[,1][combine_dmf[,1] > 0], 
      yyerr = combine_dmf[,2][combine_dmf[,1] > 0], 
      func = double_schechter, 
      prior = function(p){
        sum(
          # dnorm(p[1], x = 8.0, sd = 4.0, log = TRUE),
          # dnorm(p[2], x = -1.0, sd = 1.0, log = TRUE),
          # dnorm(p[3], x = -1.0, sd = 1.0, log = TRUE),
          # dnorm(p[4], x = -2.0, sd = 4.0, log = TRUE),
          # dnorm(p[5], x = -3.0, sd = 4.0, log = TRUE),
          0
        )
      }
    ), 
    likefunc = LL, 
    liketype = "max", 
    Niters = c(Nsamples, Nsamples), 
    NfinalMCMC = Nsamples, 
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
  
  df_fit = data.frame(
    "x" = 10^mdustvec,
    "Q50" = q50_fit,
    "Q16" = q16_fit,
    "Q84" = q84_fit
  )
  df_fit$ERR = 0.5*(df_fit$Q84 - df_fit$Q16)
  df_cosmic_samples = foreach(k = 1:dim(highout$LD_last$Posterior1)[1], .combine = c) %do% {
    trapz(
      x = mdustvec,
      y = 10^mdustvec * fit_samples[k,]
    )
  }
  df_cosmic = data.frame(
    "Q50" = quantile(df_cosmic_samples, probs = 0.5),
    "Q16" = quantile(df_cosmic_samples, probs = 0.16),
    "Q84" = quantile(df_cosmic_samples, probs = 0.84)
  )
  df_cosmic$ERR = 0.5 * (df_cosmic$Q84 - df_cosmic$Q16)
  
  combine_dmf = data.frame(combine_dmf)
  names(combine_dmf) = c("phi", "err", "WHICH")
  combine_dmf$x = sm_mids
  
  ret_ = list(
    "dmf" = combine_dmf,
    "fit" = df_fit, 
    "cosmic" = df_cosmic,
    "highout" = highout
  )
  
  png(paste0("~/Documents/DustMassDensity/plots/gmf_noAGN_corr/lbt_", round(lbt_mids[i],3), ".png"))
  magplot(
    NA,
    log = "y",
    xlim = c(3.5, 13.5),
    ylim = c(1e-8, 10),
    xlab = "Mgas",
    ylab = "Phi",
    pch = 1,
    col = "red"
  )
  lines(
    mdustvec, q50_fit, col = "red"
  )
  lines(
    mdustvec, q16_fit, col = "red", lty = 2
  )
  lines(
    mdustvec, q84_fit, col = "red", lty = 2
  )
  points(
    sm_mids,
    combine_dmf[,1],
    pch = 16,
    cex = 2,
    col = alpha("red", 0.6)
  )
  magerr(
    sm_mids,
    combine_dmf[,1],
    ylo = combine_dmf[,2],
    col = "red"
  )
  points(
    sm_mids,
    gama_dmf[,1],
    pch = "x",
    col = "purple"
  )
  magerr(
    sm_mids,
    gama_dmf[,1],
    ylo = gama_err,
    col = "purple"
  )
  points(
    sm_mids,
    devilsd10_dmf[,1],
    pch = "x",
    col = "cornflowerblue"
  )
  magerr(
    sm_mids,
    devilsd10_dmf[,1],
    ylo = devilsd10_err,
    col = "cornflowerblue"
  )
  dev.off()
  return(ret_)
}
combine_AGN_gmf_corr = foreach(i = 1:length(zmids)) %do% {
  
  gama_dmf = gama_AGN_gmf_corr[[i]]
  gama_err = gama_dmf[,2]
  
  devilsd10_dmf = devilsd10_AGN_gmf_corr[[i]]
  devilsd10_err = devilsd10_dmf[,2]
  
  gama_xlim = sm_mids[which.max(gama_dmf[,1])]
  devilsd10_xlim = sm_mids[which.max(devilsd10_dmf[,1])]
  gama_SN = gama_dmf[,1]/gama_err
  devilsd10_SN = devilsd10_dmf[,1]/devilsd10_err
  
  gama_idx = sm_mids >= gama_xlim
  devilsd10_idx = sm_mids >= devilsd10_xlim
  gama_dmf[gama_idx, ]
  
  combine_dmf = foreach(j = 1:length(sm_mids), .combine = rbind) %do% {
    if( gama_err[j] <= ifelse(devilsd10_err[j]==0, 999, devilsd10_err[j]) ){ ## check if GAMA or DEVILS has better S/N
      if( !gama_idx[j] ){ ## If below GAMA completeness, use devils 
        c(devilsd10_dmf[j,1], devilsd10_err[j], 1L)
      }else{ 
        if( (devilsd10_dmf[j,1] - gama_dmf[j,1])/(gama_err[j] + 1e-323) > 5 ){ 
          c(devilsd10_dmf[j,1], devilsd10_err[j], 1L) ## If above GAMA completeness and devils is more than 5sigma above then use devils
        }else{
          c(gama_dmf[j,1], gama_err[j], 2L)
        }
      }
    }else{
      c(devilsd10_dmf[j,1], devilsd10_err[j], 1L)
    }
  }
  combine_dmf[sm_mids < devilsd10_xlim, ] = 0
  
  pinit = c(10, 0.0, 0.0, -4, -4)
  highout = Highlander(
    parm = pinit, 
    lower = c(5, -2.0, -1.1, -8, -8),
    upper = c(15, -0.8, 1.5, 1, 1),
    Data = list(
      xx = sm_mids[combine_dmf[,1] > 0], 
      yy = combine_dmf[,1][combine_dmf[,1] > 0], 
      yyerr = combine_dmf[,2][combine_dmf[,1] > 0], 
      func = double_schechter, 
      prior = function(p){
        sum(
          # dnorm(p[1], x = 8.0, sd = 4.0, log = TRUE),
          # dnorm(p[2], x = -1.0, sd = 1.0, log = TRUE),
          # dnorm(p[3], x = -1.0, sd = 1.0, log = TRUE),
          # dnorm(p[4], x = -2.0, sd = 4.0, log = TRUE),
          # dnorm(p[5], x = -3.0, sd = 4.0, log = TRUE),
          0
        )
      }
    ), 
    likefunc = LL, 
    liketype = "max", 
    Niters = c(Nsamples, Nsamples), 
    NfinalMCMC = Nsamples, 
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
  
  df_fit = data.frame(
    "x" = 10^mdustvec,
    "Q50" = q50_fit,
    "Q16" = q16_fit,
    "Q84" = q84_fit
  )
  df_fit$ERR = 0.5*(df_fit$Q84 - df_fit$Q16)
  df_cosmic_samples = foreach(k = 1:dim(highout$LD_last$Posterior1)[1], .combine = c) %do% {
    trapz(
      x = mdustvec,
      y = 10^mdustvec * fit_samples[k,]
    )
  }
  df_cosmic = data.frame(
    "Q50" = quantile(df_cosmic_samples, probs = 0.5),
    "Q16" = quantile(df_cosmic_samples, probs = 0.16),
    "Q84" = quantile(df_cosmic_samples, probs = 0.84)
  )
  df_cosmic$ERR = 0.5 * (df_cosmic$Q84 - df_cosmic$Q16)
  
  combine_dmf = data.frame(combine_dmf)
  names(combine_dmf) = c("phi", "err")
  combine_dmf$x = sm_mids
  
  ret_ = list(
    "dmf" = combine_dmf,
    "fit" = df_fit, 
    "cosmic" = df_cosmic,
    "highout" = highout
  )
  
  png(paste0("~/Documents/DustMassDensity/plots/gmf_AGN_corr/lbt_", round(lbt_mids[i],3), ".png"))
  magplot(
    NA,
    log = "y",
    xlim = c(3.5, 13.5),
    ylim = c(1e-8, 10),
    xlab = "Mgas",
    ylab = "Phi",
    pch = 1,
    col = "red"
  )
  lines(
    mdustvec, q50_fit, col = "red"
  )
  lines(
    mdustvec, q16_fit, col = "red", lty = 2
  )
  lines(
    mdustvec, q84_fit, col = "red", lty = 2
  )
  points(
    sm_mids,
    combine_dmf[,1],
    pch = 16,
    cex = 2,
    col = alpha("red", 0.6)
  )
  magerr(
    sm_mids,
    combine_dmf[,1],
    ylo = combine_dmf[,2],
    col = "red"
  )
  points(
    sm_mids,
    gama_dmf[,1],
    pch = "x",
    col = "purple"
  )
  magerr(
    sm_mids,
    gama_dmf[,1],
    ylo = gama_err,
    col = "purple"
  )
  points(
    sm_mids,
    devilsd10_dmf[,1],
    pch = "x",
    col = "cornflowerblue"
  )
  magerr(
    sm_mids,
    devilsd10_dmf[,1],
    ylo = devilsd10_err,
    col = "cornflowerblue"
  )
  dev.off()
  return(ret_)
}
combine_hybrid_gmf_corr = foreach(i = 1:length(zmids)) %do% {
  
  gama_dmf = gama_hybrid_gmf_corr[[i]]
  gama_err = gama_dmf[,2]
  
  devilsd10_dmf = devilsd10_hybrid_gmf_corr[[i]]
  devilsd10_err = devilsd10_dmf[,2]
  
  gama_xlim = sm_mids[which.max(gama_dmf[,1])]
  devilsd10_xlim = sm_mids[which.max(devilsd10_dmf[,1])]
  gama_SN = gama_dmf[,1]/gama_err
  devilsd10_SN = devilsd10_dmf[,1]/devilsd10_err
  
  gama_idx = sm_mids >= gama_xlim
  devilsd10_idx = sm_mids >= devilsd10_xlim
  gama_dmf[gama_idx, ]
  
  combine_dmf = foreach(j = 1:length(sm_mids), .combine = rbind) %do% {
    if( gama_err[j] <= ifelse(devilsd10_err[j]==0, 999, devilsd10_err[j]) ){ ## check if GAMA or DEVILS has better S/N
      if( !gama_idx[j] ){ ## If below GAMA completeness, use devils 
        c(devilsd10_dmf[j,1], devilsd10_err[j], 1L)
      }else{ 
        if( (devilsd10_dmf[j,1] - gama_dmf[j,1])/(gama_err[j] + 1e-323) > 5 ){ 
          c(devilsd10_dmf[j,1], devilsd10_err[j], 1L) ## If above GAMA completeness and devils is more than 5sigma above then use devils
        }else{
          c(gama_dmf[j,1], gama_err[j], 2L)
        }
      }
    }else{
      c(devilsd10_dmf[j,1], devilsd10_err[j], 1L)
    }
  }
  combine_dmf[sm_mids < devilsd10_xlim, ] = 0
  
  pinit = c(10, 0.0, 0.0, -4, -4)
  highout = Highlander(
    parm = pinit, 
    lower = c(5, -2.0, -1.1, -8, -8),
    upper = c(15, -0.8, 1.5, 1, 1),
    Data = list(
      xx = sm_mids[combine_dmf[,1] > 0], 
      yy = combine_dmf[,1][combine_dmf[,1] > 0], 
      yyerr = combine_dmf[,2][combine_dmf[,1] > 0], 
      func = double_schechter, 
      prior = function(p){
        sum(
          # dnorm(p[1], x = 8.0, sd = 4.0, log = TRUE),
          # dnorm(p[2], x = -1.0, sd = 1.0, log = TRUE),
          # dnorm(p[3], x = -1.0, sd = 1.0, log = TRUE),
          # dnorm(p[4], x = -2.0, sd = 4.0, log = TRUE),
          # dnorm(p[5], x = -3.0, sd = 4.0, log = TRUE),
          0
        )
      }
    ), 
    likefunc = LL, 
    liketype = "max", 
    Niters = c(Nsamples, Nsamples), 
    NfinalMCMC = Nsamples, 
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
  
  df_fit = data.frame(
    "x" = 10^mdustvec,
    "Q50" = q50_fit,
    "Q16" = q16_fit,
    "Q84" = q84_fit
  )
  df_fit$ERR = 0.5*(df_fit$Q84 - df_fit$Q16)
  df_cosmic_samples = foreach(k = 1:dim(highout$LD_last$Posterior1)[1], .combine = c) %do% {
    trapz(
      x = mdustvec,
      y = 10^mdustvec * fit_samples[k,]
    )
  }
  df_cosmic = data.frame(
    "Q50" = quantile(df_cosmic_samples, probs = 0.5),
    "Q16" = quantile(df_cosmic_samples, probs = 0.16),
    "Q84" = quantile(df_cosmic_samples, probs = 0.84)
  )
  df_cosmic$ERR = 0.5 * (df_cosmic$Q84 - df_cosmic$Q16)
  
  combine_dmf = data.frame(combine_dmf)
  names(combine_dmf) = c("phi", "err", "WHICH")
  combine_dmf$x = sm_mids
  
  ret_ = list(
    "dmf" = combine_dmf,
    "fit" = df_fit, 
    "cosmic" = df_cosmic,
    "highout" = highout
  )
  
  png(paste0("~/Documents/DustMassDensity/plots/gmf_hybrid_corr/lbt_", round(lbt_mids[i],3), ".png"))
  magplot(
    NA,
    log = "y",
    xlim = c(3.5, 13.5),
    ylim = c(1e-8, 10),
    xlab = "Mgas",
    ylab = "Phi",
    pch = 1,
    col = "red"
  )
  lines(
    mdustvec, q50_fit, col = "red"
  )
  lines(
    mdustvec, q16_fit, col = "red", lty = 2
  )
  lines(
    mdustvec, q84_fit, col = "red", lty = 2
  )
  points(
    sm_mids,
    combine_dmf[,1],
    pch = 16,
    cex = 2,
    col = alpha("red", 0.6)
  )
  magerr(
    sm_mids,
    combine_dmf[,1],
    ylo = combine_dmf[,2],
    col = "red"
  )
  points(
    sm_mids,
    gama_dmf[,1],
    pch = "x",
    col = "purple"
  )
  magerr(
    sm_mids,
    gama_dmf[,1],
    ylo = gama_err,
    col = "purple"
  )
  points(
    sm_mids,
    devilsd10_dmf[,1],
    pch = "x",
    col = "cornflowerblue"
  )
  magerr(
    sm_mids,
    devilsd10_dmf[,1],
    ylo = devilsd10_err,
    col = "cornflowerblue"
  )
  dev.off()
  return(ret_)
}

## SMF
thorne21_evol = function(lbt){
  list(
    Mstar = 0.0006*lbt + 10.7792,
    logphi1 = -7e-8*lbt^7 - 2.5825,
    logphi2 = 0.0843*lbt - 2.6863
  )
}
combine_noAGN_smf = foreach(i = 1:length(zmids)) %do% {
  message("Stellar mass reduced integration range, Mstar=[8,12] and Thorne+21 priors")
  gama_dmf = gama_noAGN_smf[[i]]
  gama_err = gama_dmf[,2]
  
  devilsd10_dmf = devilsd10_noAGN_smf[[i]]
  devilsd10_err = devilsd10_dmf[,2]
  
  gama_xlim = sm_mids[which.max(gama_dmf[,1])]
  devilsd10_xlim = sm_mids[which.max(devilsd10_dmf[,1])]
  gama_SN = gama_dmf[,1]/gama_err
  devilsd10_SN = devilsd10_dmf[,1]/devilsd10_err
  
  gama_idx = sm_mids >= gama_xlim
  devilsd10_idx = sm_mids >= devilsd10_xlim
  gama_dmf[gama_idx, ]
  
  combine_dmf = foreach(j = 1:length(sm_mids), .combine = rbind) %do% {
    if( gama_err[j] <= ifelse(devilsd10_err[j]==0, 999, devilsd10_err[j]) ){ ## check if GAMA or DEVILS has better S/N
      if( !gama_idx[j] ){ ## If below GAMA completeness, use devils 
        c(devilsd10_dmf[j,1], devilsd10_err[j], 1L)
      }else{ 
        if( (devilsd10_dmf[j,1] - gama_dmf[j,1])/(gama_err[j] + 1e-323) > 5 ){ 
          c(devilsd10_dmf[j,1], devilsd10_err[j], 1L) ## If above GAMA completeness and devils is more than 5sigma above then use devils
        }else{
          c(gama_dmf[j,1], gama_err[j], 2L)
        }
      }
    }else{
      c(devilsd10_dmf[j,1], devilsd10_err[j], 1L)
    }
  }
  combine_dmf[sm_mids < devilsd10_xlim, ] = 0
  
  thorne21 = thorne21_evol(lbt = lbt_mids[i])
  
  pinit = c(10, 0.0, 0.0, -4, -4)
  highout = Highlander(
    parm = pinit, 
    lower = c(5, -2.0, -1.1, -8, -8),
    upper = c(15, -0.8, 1.5, 1, 1),
    Data = list(
      xx = sm_mids[combine_dmf[,1] > 0], 
      yy = combine_dmf[,1][combine_dmf[,1] > 0], 
      yyerr = combine_dmf[,2][combine_dmf[,1] > 0], 
      func = double_schechter, 
      prior = function(p){
        sum(
          dnorm(p[1], x = thorne21$Mstar, sd = 0.5, log = TRUE),
          dnorm(p[2], x = -1.5, sd = 0.1, log = TRUE),
          dnorm(p[3], x = -0.5, sd = 0.1, log = TRUE),
          dnorm(p[4], x = thorne21$logphi1, sd = 0.5, log = TRUE),
          dnorm(p[5], x = thorne21$logphi2, sd = 0.5, log = TRUE),
          0
        )
      }
    ), 
    likefunc = LL, 
    liketype = "max", 
    Niters = c(Nsamples, Nsamples), 
    NfinalMCMC = Nsamples, 
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
  
  df_fit = data.frame(
    "x" = 10^mdustvec,
    "Q50" = q50_fit,
    "Q16" = q16_fit,
    "Q84" = q84_fit
  )
  df_fit$ERR = 0.5*(df_fit$Q84 - df_fit$Q16)

  df_cosmic_samples = foreach(k = 1:dim(highout$LD_last$Posterior1)[1], .combine = c) %do% {
    trapz(
      x = mdustvec[mdustvec >= 8 & mdustvec <= 12],
      y = 10^mdustvec[mdustvec >= 8 & mdustvec <= 12] * fit_samples[k,mdustvec >= 8 & mdustvec <= 12]
    )
  }
  df_cosmic = data.frame(
    "Q50" = quantile(df_cosmic_samples, probs = 0.5),
    "Q16" = quantile(df_cosmic_samples, probs = 0.16),
    "Q84" = quantile(df_cosmic_samples, probs = 0.84)
  )
  df_cosmic$ERR = 0.5 * (df_cosmic$Q84 - df_cosmic$Q16)
  
  combine_dmf = data.frame(combine_dmf)
  names(combine_dmf) = c("phi", "err", "WHICH")
  combine_dmf$x = sm_mids
  
  ret_ = list(
    "dmf" = combine_dmf,
    "fit" = df_fit, 
    "cosmic" = df_cosmic,
    "highout" = highout
  )
  
  png(paste0("~/Documents/DustMassDensity/plots/smf_noAGN/lbt_", round(lbt_mids[i],3), ".png"))
  par(mfcol = c(1,2), oma = rep(1.5,4), mar = rep(1.5,4))
  magplot(
    NA,
    log = "y",
    xlim = c(3.5, 13.5),
    ylim = c(1e-8, 10),
    xlab = "Mstar",
    ylab = "Phi",
    pch = 1,
    col = "red"
  )
  lines(
    mdustvec, q50_fit, col = "red"
  )
  lines(
    mdustvec, q16_fit, col = "red", lty = 2
  )
  lines(
    mdustvec, q84_fit, col = "red", lty = 2
  )
  points(
    sm_mids,
    combine_dmf[,1],
    pch = 16,
    cex = 2,
    col = alpha("red", 0.6)
  )
  magerr(
    sm_mids,
    combine_dmf[,1],
    ylo = combine_dmf[,2],
    col = "red"
  )
  points(
    sm_mids,
    gama_dmf[,1],
    pch = "x",
    col = "purple"
  )
  magerr(
    sm_mids,
    gama_dmf[,1],
    ylo = gama_err,
    col = "purple"
  )
  points(
    sm_mids,
    devilsd10_dmf[,1],
    pch = "x",
    col = "cornflowerblue"
  )
  magerr(
    sm_mids,
    devilsd10_dmf[,1],
    ylo = devilsd10_err,
    col = "cornflowerblue"
  )
  magplot(
    NA,
    log = "y",
    xlim = c(3.5, 13.5),
    ylim = c(10**3.5, 10**9.5),
    xlab = "Mstar",
    ylab = "Phi",
    pch = 1,
    col = "red"
  )
  lines(
    mdustvec, q50_fit * 10^mdustvec, col = "red"
  )
  lines(
    mdustvec, q16_fit * 10^mdustvec, col = "red", lty = 2
  )
  lines(
    mdustvec, q84_fit * 10^mdustvec, col = "red", lty = 2
  )
  dev.off()
  return(ret_)
}
combine_AGN_smf = foreach(i = 1:length(zmids)) %do% {
  message("Stellar mass reduced integration range, Mstar=[8,12] and Thorne+21 priors")
  gama_dmf = gama_AGN_smf[[i]]
  gama_err = gama_dmf[,2]
  
  devilsd10_dmf = devilsd10_AGN_smf[[i]]
  devilsd10_err = devilsd10_dmf[,2]
  
  gama_xlim = sm_mids[which.max(gama_dmf[,1])]
  devilsd10_xlim = sm_mids[which.max(devilsd10_dmf[,1])]
  gama_SN = gama_dmf[,1]/gama_err
  devilsd10_SN = devilsd10_dmf[,1]/devilsd10_err
  
  gama_idx = sm_mids >= gama_xlim
  devilsd10_idx = sm_mids >= devilsd10_xlim
  gama_dmf[gama_idx, ]
  
  combine_dmf = foreach(j = 1:length(sm_mids), .combine = rbind) %do% {
    if( gama_err[j] <= ifelse(devilsd10_err[j]==0, 999, devilsd10_err[j]) ){ ## check if GAMA or DEVILS has better S/N
      if( !gama_idx[j] ){ ## If below GAMA completeness, use devils 
        c(devilsd10_dmf[j,1], devilsd10_err[j], 1L)
      }else{ 
        if( (devilsd10_dmf[j,1] - gama_dmf[j,1])/(gama_err[j] + 1e-323) > 5 ){ 
          c(devilsd10_dmf[j,1], devilsd10_err[j], 1L) ## If above GAMA completeness and devils is more than 5sigma above then use devils
        }else{
          c(gama_dmf[j,1], gama_err[j], 2L)
        }
      }
    }else{
      c(devilsd10_dmf[j,1], devilsd10_err[j], 1L)
    }
  }
  combine_dmf[sm_mids < devilsd10_xlim, ] = 0
  
  thorne21 = thorne21_evol(lbt = lbt_mids[i])
  
  pinit = c(10, 0.0, 0.0, -4, -4)
  highout = Highlander(
    parm = pinit, 
    lower = c(5, -2.0, -1.1, -8, -8),
    upper = c(15, -0.8, 1.5, 1, 1),
    Data = list(
      xx = sm_mids[combine_dmf[,1] > 0], 
      yy = combine_dmf[,1][combine_dmf[,1] > 0], 
      yyerr = combine_dmf[,2][combine_dmf[,1] > 0], 
      func = double_schechter, 
      prior = function(p){
        sum(
          dnorm(p[1], x = thorne21$Mstar, sd = 0.5, log = TRUE),
          dnorm(p[2], x = -1.5, sd = 0.1, log = TRUE),
          dnorm(p[3], x = -0.5, sd = 0.1, log = TRUE),
          dnorm(p[4], x = thorne21$logphi1, sd = 0.5, log = TRUE),
          dnorm(p[5], x = thorne21$logphi2, sd = 0.5, log = TRUE),
          0
        )
      }
    ), 
    likefunc = LL, 
    liketype = "max", 
    Niters = c(Nsamples, Nsamples), 
    NfinalMCMC = Nsamples, 
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
  
  df_fit = data.frame(
    "x" = 10^mdustvec,
    "Q50" = q50_fit,
    "Q16" = q16_fit,
    "Q84" = q84_fit
  )
  df_fit$ERR = 0.5*(df_fit$Q84 - df_fit$Q16)
  
  df_cosmic_samples = foreach(k = 1:dim(highout$LD_last$Posterior1)[1], .combine = c) %do% {
    trapz(
      x = mdustvec[mdustvec >= 8 & mdustvec <= 12],
      y = 10^mdustvec[mdustvec >= 8 & mdustvec <= 12] * fit_samples[k,mdustvec >= 8 & mdustvec <= 12]
    )
  }
  df_cosmic = data.frame(
    "Q50" = quantile(df_cosmic_samples, probs = 0.5),
    "Q16" = quantile(df_cosmic_samples, probs = 0.16),
    "Q84" = quantile(df_cosmic_samples, probs = 0.84)
  )
  df_cosmic$ERR = 0.5 * (df_cosmic$Q84 - df_cosmic$Q16)
  
  combine_dmf = data.frame(combine_dmf)
  names(combine_dmf) = c("phi", "err", "WHICH")
  combine_dmf$x = sm_mids
  
  ret_ = list(
    "dmf" = combine_dmf,
    "fit" = df_fit, 
    "cosmic" = df_cosmic,
    "highout" = highout
  )
  
  png(paste0("~/Documents/DustMassDensity/plots/smf_AGN/lbt_", round(lbt_mids[i],3), ".png"))
  par(mfcol = c(1,2), oma = rep(1.5,4), mar = rep(1.5,4))
  magplot(
    NA,
    log = "y",
    xlim = c(3.5, 13.5),
    ylim = c(1e-8, 10),
    xlab = "Mstar",
    ylab = "Phi",
    pch = 1,
    col = "red"
  )
  lines(
    mdustvec, q50_fit, col = "red"
  )
  lines(
    mdustvec, q16_fit, col = "red", lty = 2
  )
  lines(
    mdustvec, q84_fit, col = "red", lty = 2
  )
  points(
    sm_mids,
    combine_dmf[,1],
    pch = 16,
    cex = 2,
    col = alpha("red", 0.6)
  )
  magerr(
    sm_mids,
    combine_dmf[,1],
    ylo = combine_dmf[,2],
    col = "red"
  )
  points(
    sm_mids,
    gama_dmf[,1],
    pch = "x",
    col = "purple"
  )
  magerr(
    sm_mids,
    gama_dmf[,1],
    ylo = gama_err,
    col = "purple"
  )
  points(
    sm_mids,
    devilsd10_dmf[,1],
    pch = "x",
    col = "cornflowerblue"
  )
  magerr(
    sm_mids,
    devilsd10_dmf[,1],
    ylo = devilsd10_err,
    col = "cornflowerblue"
  )
  magplot(
    NA,
    log = "y",
    xlim = c(3.5, 13.5),
    ylim = c(10**3.5, 10**9.5),
    xlab = "Mstar",
    ylab = "Phi",
    pch = 1,
    col = "red"
  )
  lines(
    mdustvec, q50_fit * 10^mdustvec, col = "red"
  )
  lines(
    mdustvec, q16_fit * 10^mdustvec, col = "red", lty = 2
  )
  lines(
    mdustvec, q84_fit * 10^mdustvec, col = "red", lty = 2
  )
  dev.off()
  return(ret_)
}
combine_hybrid_smf = foreach(i = 1:length(zmids)) %do% {
  message("Stellar mass reduced integration range, Mstar=[8,12] and Thorne+21 priors")
  gama_dmf = gama_hybrid_smf[[i]]
  gama_err = gama_dmf[,2]
  
  devilsd10_dmf = devilsd10_hybrid_smf[[i]]
  devilsd10_err = devilsd10_dmf[,2]
  
  gama_xlim = sm_mids[which.max(gama_dmf[,1])]
  devilsd10_xlim = sm_mids[which.max(devilsd10_dmf[,1])]
  gama_SN = gama_dmf[,1]/gama_err
  devilsd10_SN = devilsd10_dmf[,1]/devilsd10_err
  
  gama_idx = sm_mids >= gama_xlim
  devilsd10_idx = sm_mids >= devilsd10_xlim
  gama_dmf[gama_idx, ]
  
  combine_dmf = foreach(j = 1:length(sm_mids), .combine = rbind) %do% {
    if( gama_err[j] <= ifelse(devilsd10_err[j]==0, 999, devilsd10_err[j]) ){ ## check if GAMA or DEVILS has better S/N
      if( !gama_idx[j] ){ ## If below GAMA completeness, use devils 
        c(devilsd10_dmf[j,1], devilsd10_err[j], 1L)
      }else{ 
        if( (devilsd10_dmf[j,1] - gama_dmf[j,1])/(gama_err[j] + 1e-323) > 5 ){ 
          c(devilsd10_dmf[j,1], devilsd10_err[j], 1L) ## If above GAMA completeness and devils is more than 5sigma above then use devils
        }else{
          c(gama_dmf[j,1], gama_err[j], 2L)
        }
      }
    }else{
      c(devilsd10_dmf[j,1], devilsd10_err[j], 1L)
    }
  }
  combine_dmf[sm_mids < devilsd10_xlim, ] = 0
  
  thorne21 = thorne21_evol(lbt = lbt_mids[i])
  
  pinit = c(10, 0.0, 0.0, -4, -4)
  highout = Highlander(
    parm = pinit, 
    lower = c(5, -2.0, -1.1, -8, -8),
    upper = c(15, -0.8, 1.5, 1, 1),
    Data = list(
      xx = sm_mids[combine_dmf[,1] > 0], 
      yy = combine_dmf[,1][combine_dmf[,1] > 0], 
      yyerr = combine_dmf[,2][combine_dmf[,1] > 0], 
      func = double_schechter, 
      prior = function(p){
        sum(
          dnorm(p[1], x = thorne21$Mstar, sd = 0.5, log = TRUE),
          dnorm(p[2], x = -1.5, sd = 0.1, log = TRUE),
          dnorm(p[3], x = -0.5, sd = 0.1, log = TRUE),
          dnorm(p[4], x = thorne21$logphi1, sd = 0.5, log = TRUE),
          dnorm(p[5], x = thorne21$logphi2, sd = 0.5, log = TRUE),
          0
        )
      }
    ), 
    likefunc = LL, 
    liketype = "max", 
    Niters = c(Nsamples, Nsamples), 
    NfinalMCMC = Nsamples, 
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
  
  df_fit = data.frame(
    "x" = 10^mdustvec,
    "Q50" = q50_fit,
    "Q16" = q16_fit,
    "Q84" = q84_fit
  )
  df_fit$ERR = 0.5*(df_fit$Q84 - df_fit$Q16)
  
  df_cosmic_samples = foreach(k = 1:dim(highout$LD_last$Posterior1)[1], .combine = c) %do% {
    trapz(
      x = mdustvec[mdustvec >= 8 & mdustvec <= 12],
      y = 10^mdustvec[mdustvec >= 8 & mdustvec <= 12] * fit_samples[k,mdustvec >= 8 & mdustvec <= 12]
    )
  }
  df_cosmic = data.frame(
    "Q50" = quantile(df_cosmic_samples, probs = 0.5),
    "Q16" = quantile(df_cosmic_samples, probs = 0.16),
    "Q84" = quantile(df_cosmic_samples, probs = 0.84)
  )
  df_cosmic$ERR = 0.5 * (df_cosmic$Q84 - df_cosmic$Q16)
  
  combine_dmf = data.frame(combine_dmf)
  names(combine_dmf) = c("phi", "err", "WHICH")
  combine_dmf$x = sm_mids
  
  ret_ = list(
    "dmf" = combine_dmf,
    "fit" = df_fit, 
    "cosmic" = df_cosmic,
    "highout" = highout
  )
  
  png(paste0("~/Documents/DustMassDensity/plots/smf_hybrid/lbt_", round(lbt_mids[i],3), ".png"))
  par(mfcol = c(1,2), oma = rep(1.5,4), mar = rep(1.5,4))
  magplot(
    NA,
    log = "y",
    xlim = c(3.5, 13.5),
    ylim = c(1e-8, 10),
    xlab = "Mstar",
    ylab = "Phi",
    pch = 1,
    col = "red"
  )
  lines(
    mdustvec, q50_fit, col = "red"
  )
  lines(
    mdustvec, q16_fit, col = "red", lty = 2
  )
  lines(
    mdustvec, q84_fit, col = "red", lty = 2
  )
  points(
    sm_mids,
    combine_dmf[,1],
    pch = 16,
    cex = 2,
    col = alpha("red", 0.6)
  )
  magerr(
    sm_mids,
    combine_dmf[,1],
    ylo = combine_dmf[,2],
    col = "red"
  )
  points(
    sm_mids,
    gama_dmf[,1],
    pch = "x",
    col = "purple"
  )
  magerr(
    sm_mids,
    gama_dmf[,1],
    ylo = gama_err,
    col = "purple"
  )
  points(
    sm_mids,
    devilsd10_dmf[,1],
    pch = "x",
    col = "cornflowerblue"
  )
  magerr(
    sm_mids,
    devilsd10_dmf[,1],
    ylo = devilsd10_err,
    col = "cornflowerblue"
  )
  magplot(
    NA,
    log = "y",
    xlim = c(3.5, 13.5),
    ylim = c(10**3.5, 10**9.5),
    xlab = "Mstar",
    ylab = "Phi",
    pch = 1,
    col = "red"
  )
  lines(
    mdustvec, q50_fit * 10^mdustvec, col = "red"
  )
  lines(
    mdustvec, q16_fit * 10^mdustvec, col = "red", lty = 2
  )
  lines(
    mdustvec, q84_fit * 10^mdustvec, col = "red", lty = 2
  )
  dev.off()
  return(ret_)
}

## Only FIR DEVILS sources 
devilsd10_hybrid_cdmh_FIR = (foreach(i = 1:length(zmids)) %do% {
  message(i)
  vol = vol_mids[i] * 1.47 / (4*pi*(180/pi)^2)
  zidx = devilsd10_hybrid$z >= zbins[i] & devilsd10_hybrid$z < zbins[i+1] & devilsd10_hybrid$FIRInput == 1
  Mdust = devilsd10_hybrid$dustmass.total[zidx] * RR14_BPL(Z = 10^devilsd10_hybrid$Zfinal[zidx], doDTG = FALSE)/0.0073
  
  df = sum(Mdust) / vol
  # mc_err = 0.5*(mc_err_devilsd10_hybrid_dmf_corr[[i]][,3] - mc_err_devilsd10_hybrid_dmf_corr[[i]][,2])
  # 
  # log_m = log10(Mdust[Mdust > 0])
  # hh = maghist(x = log_m, breaks = sm_bins, plot = FALSE, verbose = FALSE)
  # bin_dmf = hh$counts/(vol * diff(sm_bins))
  # pois_dmf = sqrt(hh$counts)/(vol * diff(sm_bins))
  # 
  # err_dmf = sqrt( pois_dmf^2 + (err_floor*bin_dmf)^2 + do_MC_err*mc_err^2 )
  # df = data.frame(cbind(bin_dmf, err_dmf))
  # names(df) = c("phi", "err")
  # df$x = hh$mids
  return(df)
})
devilsd10_hybrid_cdmh_FIR = data.frame("CDMH_FIR" = unlist(devilsd10_hybrid_cdmh_FIR))k
## Make cosmic
cdmh_noAGN = data.frame(foreach(i = 1:length(combine_noAGN_dmf), .combine = bind_rows) %do% {
  combine_noAGN_dmf[[i]]$cosmic
})
cdmh_AGN = data.frame(foreach(i = 1:length(combine_AGN_dmf), .combine = bind_rows) %do% {
  combine_AGN_dmf[[i]]$cosmic
})
cdmh_hybrid = data.frame(foreach(i = 1:length(combine_hybrid_dmf), .combine = bind_rows) %do% {
  combine_hybrid_dmf[[i]]$cosmic
})

cdmh_noAGN_corr = data.frame(foreach(i = 1:length(combine_noAGN_dmf_corr), .combine = bind_rows) %do% {
  combine_noAGN_dmf_corr[[i]]$cosmic
})
cdmh_AGN_corr = data.frame(foreach(i = 1:length(combine_AGN_dmf_corr), .combine = bind_rows) %do% {
  combine_AGN_dmf_corr[[i]]$cosmic
})
cdmh_hybrid_corr = data.frame(foreach(i = 1:length(combine_hybrid_dmf_corr), .combine = bind_rows) %do% {
  combine_hybrid_dmf_corr[[i]]$cosmic
})
cdmh_hybrid_corr_Zsol = data.frame(foreach(i = 1:length(combine_hybrid_dmf_corr), .combine = bind_rows) %do% {
  combine_hybrid_dmf_corr_Zsol[[i]]$cosmic
})

cgmh_noAGN_corr = data.frame(foreach(i = 1:length(combine_noAGN_gmf_corr), .combine = bind_rows) %do% {
  combine_noAGN_gmf_corr[[i]]$cosmic
})
cgmh_AGN_corr = data.frame(foreach(i = 1:length(combine_AGN_gmf_corr), .combine = bind_rows) %do% {
  combine_AGN_gmf_corr[[i]]$cosmic
})
cgmh_hybrid_corr = data.frame(foreach(i = 1:length(combine_hybrid_gmf_corr), .combine = bind_rows) %do% {
  combine_hybrid_gmf_corr[[i]]$cosmic
})

csmh_noAGN = data.frame(foreach(i = 1:length(combine_noAGN_smf), .combine = bind_rows) %do% {
  combine_noAGN_smf[[i]]$cosmic
})
csmh_AGN = data.frame(foreach(i = 1:length(combine_AGN_smf), .combine = bind_rows) %do% {
  combine_AGN_smf[[i]]$cosmic
})
csmh_hybrid = data.frame(foreach(i = 1:length(combine_hybrid_smf), .combine = bind_rows) %do% {
  combine_hybrid_smf[[i]]$cosmic
})

# driver18 = fread("~/Documents/DustMassDensity/data/literature_evo/cdmh/driver18.csv")
# magplot(
#   driver18$lbt, driver18$cdmh, log = "y", ylim = c(1e4, 1e7)
# )
# lines(
#   lbt_mids, cdmh_hybrid_corr$Q50
# )

names_par = c("M", "alpha", "beta", "phi1", "phi2")
names_par_err = paste0(names_par, "Err")
combine_noAGN_dmf_par = data.frame(foreach(i = 1:length(combine_noAGN_dmf), .combine = rbind) %do% {
  temp = colQuantiles(combine_noAGN_dmf[[i]]$highout$LD_last$Posterior1, probs = c(0.5, 0.16, 0.84))
  cbind(
    t(temp[,1])
  )
})
combine_noAGN_dmf_par_err = data.frame(foreach(i = 1:length(combine_noAGN_dmf), .combine = rbind) %do% {
  temp = colQuantiles(combine_noAGN_dmf[[i]]$highout$LD_last$Posterior1, probs = c(0.5, 0.16, 0.84))
  cbind(
    t(0.5 * (temp[,3] - temp[,2]))
  )
})
combine_AGN_dmf_par = data.frame(foreach(i = 1:length(combine_AGN_dmf), .combine = rbind) %do% {
  temp = colQuantiles(combine_AGN_dmf[[i]]$highout$LD_last$Posterior1, probs = c(0.5, 0.16, 0.84))
  cbind(
    t(temp[,1])
  )
})
combine_AGN_dmf_par_err = data.frame(foreach(i = 1:length(combine_AGN_dmf), .combine = rbind) %do% {
  temp = colQuantiles(combine_AGN_dmf[[i]]$highout$LD_last$Posterior1, probs = c(0.5, 0.16, 0.84))
  cbind(
    t(0.5 * (temp[,3] - temp[,2]))
  )
})
combine_hybrid_dmf_par = data.frame(foreach(i = 1:length(combine_hybrid_dmf), .combine = rbind) %do% {
  temp = colQuantiles(combine_hybrid_dmf[[i]]$highout$LD_last$Posterior1, probs = c(0.5, 0.16, 0.84))
  cbind(
    t(temp[,1])
  )
})
combine_hybrid_dmf_par_err = data.frame(foreach(i = 1:length(combine_hybrid_dmf), .combine = rbind) %do% {
  temp = colQuantiles(combine_hybrid_dmf[[i]]$highout$LD_last$Posterior1, probs = c(0.5, 0.16, 0.84))
  cbind(
    t(0.5 * (temp[,3] - temp[,2]))
  )
})
combine_hybrid_dmf_corr_par = data.frame(foreach(i = 1:length(combine_hybrid_dmf_corr), .combine = rbind) %do% {
  temp = colQuantiles(combine_hybrid_dmf_corr[[i]]$highout$LD_last$Posterior1, probs = c(0.5, 0.16, 0.84))
  cbind(
    t(temp[,1])
  )
})
combine_hybrid_dmf_corr_par_err = data.frame(foreach(i = 1:length(combine_hybrid_dmf_corr), .combine = rbind) %do% {
  temp = colQuantiles(combine_hybrid_dmf_corr[[i]]$highout$LD_last$Posterior1, probs = c(0.5, 0.16, 0.84))
  cbind(
    t(0.5 * (temp[,3] - temp[,2]))
  )
})
combine_hybrid_gmf_corr_par = data.frame(foreach(i = 1:length(combine_hybrid_gmf_corr), .combine = rbind) %do% {
  temp = colQuantiles(combine_hybrid_gmf_corr[[i]]$highout$LD_last$Posterior1, probs = c(0.5, 0.16, 0.84))
  cbind(
    t(temp[,1])
  )
})
combine_hybrid_gmf_corr_par_err = data.frame(foreach(i = 1:length(combine_hybrid_gmf_corr), .combine = rbind) %do% {
  temp = colQuantiles(combine_hybrid_gmf_corr[[i]]$highout$LD_last$Posterior1, probs = c(0.5, 0.16, 0.84))
  cbind(
    t(0.5 * (temp[,3] - temp[,2]))
  )
})
combine_hybrid_smf_par = data.frame(foreach(i = 1:length(combine_hybrid_smf), .combine = rbind) %do% {
  temp = colQuantiles(combine_hybrid_smf[[i]]$highout$LD_last$Posterior1, probs = c(0.5, 0.16, 0.84))
  cbind(
    t(temp[,1])
  )
})
combine_hybrid_smf_par_err = data.frame(foreach(i = 1:length(combine_hybrid_smf), .combine = rbind) %do% {
  temp = colQuantiles(combine_hybrid_smf[[i]]$highout$LD_last$Posterior1, probs = c(0.5, 0.16, 0.84))
  cbind(
    t(0.5 * (temp[,3] - temp[,2]))
  )
})
names(combine_noAGN_dmf_par) = names_par
names(combine_AGN_dmf_par) = names_par
names(combine_hybrid_dmf_par) = names_par
names(combine_hybrid_dmf_corr_par) = names_par
names(combine_hybrid_gmf_corr_par) = names_par
names(combine_hybrid_smf_par) = names_par
names(combine_noAGN_dmf_par_err) = names_par_err
names(combine_AGN_dmf_par_err) = names_par_err
names(combine_hybrid_dmf_par_err) = names_par_err
names(combine_hybrid_dmf_corr_par_err) = names_par_err
names(combine_hybrid_gmf_corr_par_err) = names_par_err
names(combine_hybrid_smf_par_err) = names_par_err

combine_noAGN_dmf_par = data.frame(cbind(combine_noAGN_dmf_par, combine_noAGN_dmf_par_err))
combine_AGN_dmf_par = data.frame(cbind(combine_AGN_dmf_par, combine_AGN_dmf_par_err))
combine_hybrid_dmf_par = data.frame(cbind(combine_hybrid_dmf_par, combine_hybrid_dmf_par_err))
combine_hybrid_dmf_corr_par = data.frame(cbind(combine_hybrid_dmf_corr_par, combine_hybrid_dmf_corr_par_err))
combine_hybrid_gmf_corr_par = data.frame(cbind(combine_hybrid_gmf_corr_par, combine_hybrid_gmf_corr_par_err))
combine_hybrid_smf_par = data.frame(cbind(combine_hybrid_smf_par, combine_hybrid_smf_par_err))

# LSS correction
LL_SMF_evol =  function(p, Data){
  sum(dnorm(
    x = Data$x,
    mean = Data$fitfunc(p, Data$z),
    sd = Data$err,
    log = TRUE
  ))
}

Mstar_Evol = optim(
  par = c(11.0, 0.5, -0.1),
  fn = LL_SMF_evol,
  Data = list(
    x = combine_hybrid_smf_par$M,
    err = combine_hybrid_smf_par$MErr,
    lbt = lbt_mids,
    z = zmids,
    fitfunc = function(p, t){p[1] + p[2]*t + p[3]*t^2}
  ),
  control = list(fnscale = -1),
  hessian = TRUE
)
alpha_Evol = optim(
  par = c(-0.5, 0.0, 0.0),
  fn = LL_SMF_evol,
  Data = list(
    x = combine_hybrid_smf_par$alpha,
    err = combine_hybrid_smf_par$alphaErr,
    lbt = lbt_mids,
    z = zmids,
    fitfunc = function(p, t){p[1] + p[2]*t + p[3]*t^2}
  ),
  control = list(fnscale = -1),
  hessian = TRUE
)
beta_Evol = optim(
  par =  c(-0.5, 0.0, 0.0),
  fn = LL_SMF_evol,
  Data = list(
    x = combine_hybrid_smf_par$beta,
    err = combine_hybrid_smf_par$betaErr,
    lbt = lbt_mids,
    z = zmids,
    fitfunc = function(p, t){p[1] + p[2]*t + p[3]*t^2}
  ),
  control = list(fnscale = -1),
  hessian = TRUE
)
phi1_Evol = optim(
  par = c(-2.0, 0.0, 0.0),
  fn = LL_SMF_evol,
  Data = list(
    x = combine_hybrid_smf_par$phi1,
    err = combine_hybrid_smf_par$phi1Err,
    lbt = lbt_mids,
    z = zmids,
    fitfunc = function(p, t){p[1] + p[2]*t + p[3]*t^2}
  ),
  control = list(fnscale = -1),
  hessian = TRUE
)
phi2_Evol = optim(
  par = c(-3.0, 0.0, 0.0),
  fn = LL_SMF_evol,
  Data = list(
    x = combine_hybrid_smf_par$phi2,
    err = combine_hybrid_smf_par$phi2Err,
    lbt = lbt_mids,
    z = zmids,
    fitfunc = function(p, t){p[1] + p[2]*t + p[3]*t^2}
  ),
  control = list(fnscale = -1),
  hessian = TRUE
)

# par(mfcol = c(5,1), oma = rep(0.5,4), mar = c(3.5, 3.5, 0.0, 0.0))
# magplot(
#   zmids, combine_hybrid_smf_par$M, ylim = c(10,11.5), ylab = "M*"
# )
# magerr(
#   zmids, combine_hybrid_smf_par$M, ylo = combine_hybrid_smf_par$MErr
# )
# lines(
#   zmids, Mstar_Evol$par[1] + zmids*Mstar_Evol$par[2] + zmids^2*Mstar_Evol$par[3]
# )
# magplot(
#   zmids, combine_hybrid_smf_par$phi1, ylim = c(-4, -2.5), ylab = "Phi1"
# )
# magerr(
#   zmids, combine_hybrid_smf_par$phi1, ylo = combine_hybrid_smf_par$phi1Err
# )
# lines(
#   zmids, phi1_Evol$par[1] + zmids*phi1_Evol$par[2] + zmids^2*phi1_Evol$par[3]
# )
# magplot(
#   zmids, combine_hybrid_smf_par$phi2, ylim = c(-5, -2), ylab = "Phi2"
# )
# magerr(
#   zmids, combine_hybrid_smf_par$phi2, ylo = combine_hybrid_smf_par$phi2Err
# )
# lines(
#   zmids, phi2_Evol$par[1] + zmids*phi2_Evol$par[2] + zmids^2*phi2_Evol$par[3]
# )
# magplot(
#   zmids, combine_hybrid_smf_par$alpha, ylim = c(-3.5,-0.5), ylab = "alpha"
# )
# magerr(
#   zmids, combine_hybrid_smf_par$alpha, ylo = combine_hybrid_smf_par$alphaErr
# )
# lines(
#   zmids, alpha_Evol$par[1] + zmids*alpha_Evol$par[2] + zmids^2*alpha_Evol$par[3]
# )
# magplot(
#   zmids, combine_hybrid_smf_par$beta, ylim = c(-1.0,0.0), ylab = "beta", xlab = "Lookback time [Gyr]"
# )
# magerr(
#   zmids, combine_hybrid_smf_par$beta, ylo = combine_hybrid_smf_par$betaErr
# )
# lines(
#   zmids, beta_Evol$par[1] + zmids*beta_Evol$par[2] + zmids^2*beta_Evol$par[3]
# )

Mstar_Evol_Samples = mvrnorm(
  n = Nsamples,
  mu = Mstar_Evol$par,
  Sigma = -1*solve(Mstar_Evol$hessian)
)
alpha_Evol_Samples = mvrnorm(
  n = Nsamples,
  mu = alpha_Evol$par,
  Sigma = -1*solve(alpha_Evol$hessian)
)
beta_Evol_Samples = mvrnorm(
  n = Nsamples,
  mu = beta_Evol$par,
  Sigma = -1*solve(beta_Evol$hessian)
)
phi1_Evol_Samples = mvrnorm(
  n = Nsamples,
  mu = phi1_Evol$par,
  Sigma = -1*solve(phi1_Evol$hessian)
)
phi2_Evol_Samples = mvrnorm(
  n = Nsamples,
  mu = phi2_Evol$par,
  Sigma = -1*solve(phi2_Evol$hessian)
)

smf_regressed_fits = c(
  Mstar_Evol$par,
  sqrt(diag(-1*solve(Mstar_Evol$hessian))),
  alpha_Evol$par,
  sqrt(diag(-1*solve(alpha_Evol$hessian))),
  beta_Evol$par,
  sqrt(diag(-1*solve(beta_Evol$hessian))),
  phi1_Evol$par,
  sqrt(diag(-1*solve(phi1_Evol$hessian))),
  phi2_Evol$par,
  sqrt(diag(-1*solve(phi2_Evol$hessian)))
)
names(smf_regressed_fits) = c("MstarT1", "MstarT2", "MstarT3", "MstarT1Err", "MstarT2Err", "MstarT3Err", "alphaT1", "alphaT2", "alphaT3", "alphaT1Err", "alphaT2Err", "alphaT3Err", "betaT1", "betaT2", "betaT3", "betaT1Err", "betaT2Err", "betaT3Err", "phi1T1", "phi1T2", "phi1T3", "phi1T1Err", "phi1T2Err", "phi1T3Err", "phi2T1", "phi2T2", "phi2T3", "phi2T1Err", "phi2T2Err", "phi2T3Err")
smf_regressed_fits = data.frame(t(smf_regressed_fits))

smf_regressed_evol = data.frame(
  "MstarQ50" = colQuantiles(as.matrix(foreach(i = 1:Nsamples, .combine = rbind)%do%{Mstar_Evol_Samples[i,1] + Mstar_Evol_Samples[i,2]*zmids + Mstar_Evol_Samples[i,3]*zmids^2}), probs = 0.50),
  "MstarQ16" = colQuantiles(as.matrix(foreach(i = 1:Nsamples, .combine = rbind)%do%{Mstar_Evol_Samples[i,1] + Mstar_Evol_Samples[i,2]*zmids + Mstar_Evol_Samples[i,3]*zmids^2}), probs = 0.16),
  "MstarQ84" = colQuantiles(as.matrix(foreach(i = 1:Nsamples, .combine = rbind)%do%{Mstar_Evol_Samples[i,1] + Mstar_Evol_Samples[i,2]*zmids + Mstar_Evol_Samples[i,3]*zmids^2}), probs = 0.84),
  
  "alphaQ50" = colQuantiles(as.matrix(foreach(i = 1:Nsamples, .combine = rbind)%do%{alpha_Evol_Samples[i,1] + alpha_Evol_Samples[i,2]*zmids + alpha_Evol_Samples[i,3]*zmids^2}), probs = 0.50),
  "alphaQ16" = colQuantiles(as.matrix(foreach(i = 1:Nsamples, .combine = rbind)%do%{alpha_Evol_Samples[i,1] + alpha_Evol_Samples[i,2]*zmids + alpha_Evol_Samples[i,3]*zmids^2}), probs = 0.16),
  "alphaQ84" = colQuantiles(as.matrix(foreach(i = 1:Nsamples, .combine = rbind)%do%{alpha_Evol_Samples[i,1] + alpha_Evol_Samples[i,2]*zmids + alpha_Evol_Samples[i,3]*zmids^2}), probs = 0.84),
  
  "betaQ50" = colQuantiles(as.matrix(foreach(i = 1:Nsamples, .combine = rbind)%do%{beta_Evol_Samples[i,1] + beta_Evol_Samples[i,2]*zmids + beta_Evol_Samples[i,3]*zmids^2}), probs = 0.50),
  "betaQ16" = colQuantiles(as.matrix(foreach(i = 1:Nsamples, .combine = rbind)%do%{beta_Evol_Samples[i,1] + beta_Evol_Samples[i,2]*zmids + beta_Evol_Samples[i,3]*zmids^2}), probs = 0.16),
  "betaQ84" = colQuantiles(as.matrix(foreach(i = 1:Nsamples, .combine = rbind)%do%{beta_Evol_Samples[i,1] + beta_Evol_Samples[i,2]*zmids + beta_Evol_Samples[i,3]*zmids^2}), probs = 0.84),
  
  "phi1Q50" = colQuantiles(as.matrix(foreach(i = 1:Nsamples, .combine = rbind)%do%{phi1_Evol_Samples[i,1] + phi1_Evol_Samples[i,2]*zmids + phi1_Evol_Samples[i,3]*zmids^2}), probs = 0.50),
  "phi1Q16" = colQuantiles(as.matrix(foreach(i = 1:Nsamples, .combine = rbind)%do%{phi1_Evol_Samples[i,1] + phi1_Evol_Samples[i,2]*zmids + phi1_Evol_Samples[i,3]*zmids^2}), probs = 0.16),
  "phi1Q84" = colQuantiles(as.matrix(foreach(i = 1:Nsamples, .combine = rbind)%do%{phi1_Evol_Samples[i,1] + phi1_Evol_Samples[i,2]*zmids + phi1_Evol_Samples[i,3]*zmids^2}), probs = 0.84),
  
  "phi2Q50" = colQuantiles(as.matrix(foreach(i = 1:Nsamples, .combine = rbind)%do%{phi2_Evol_Samples[i,1] + phi2_Evol_Samples[i,2]*zmids + phi2_Evol_Samples[i,3]*zmids^2}), probs = 0.50),
  "phi2Q16" = colQuantiles(as.matrix(foreach(i = 1:Nsamples, .combine = rbind)%do%{phi2_Evol_Samples[i,1] + phi2_Evol_Samples[i,2]*zmids + phi2_Evol_Samples[i,3]*zmids^2}), probs = 0.16),
  "phi2Q84" = colQuantiles(as.matrix(foreach(i = 1:Nsamples, .combine = rbind)%do%{phi2_Evol_Samples[i,1] + phi2_Evol_Samples[i,2]*zmids + phi2_Evol_Samples[i,3]*zmids^2}), probs = 0.84)
)
smf_regressed_evol$MstarErr = 0.5 * (smf_regressed_evol$MstarQ84 - smf_regressed_evol$MstarQ16)
smf_regressed_evol$alphaErr = 0.5 * (smf_regressed_evol$alphaQ84 - smf_regressed_evol$alphaQ16)
smf_regressed_evol$betaErr = 0.5 * (smf_regressed_evol$betaQ84 - smf_regressed_evol$betaQ16)
smf_regressed_evol$phi1Err = 0.5 * (smf_regressed_evol$phi1Q84 - smf_regressed_evol$phi1Q16)
smf_regressed_evol$phi2Err = 0.5 * (smf_regressed_evol$phi2Q84 - smf_regressed_evol$phi2Q16)

smf_regressed_samples = foreach(i = 1:length(lbt_mids)) %do% {
  samples_ = foreach(j = 1:Nsamples, .combine = rbind) %do% {
    M = Mstar_Evol_Samples[j,1] +     Mstar_Evol_Samples[j,2] * zmids[i] + Mstar_Evol_Samples[j,3] * zmids[i]^2
    phi1 = phi1_Evol_Samples[j,1] +   phi1_Evol_Samples[j,2]  * zmids[i] + phi1_Evol_Samples[j,3]  * zmids[i]^2
    phi2 = phi2_Evol_Samples[j,1] +   phi2_Evol_Samples[j,2]  * zmids[i] + phi2_Evol_Samples[j,3]  * zmids[i]^2
    alpha = alpha_Evol_Samples[j,1] + alpha_Evol_Samples[j,2] * zmids[i] + alpha_Evol_Samples[j,3] * zmids[i]^2
    beta = beta_Evol_Samples[j,1] +   beta_Evol_Samples[j,2]  * zmids[i] + beta_Evol_Samples[j,3]  * zmids[i]^2
    
    regressed_Schechter = double_schechter(
      x = mdustvec,
      p = c(M, alpha, beta, phi1, phi2)
    )
  }
  df = colQuantiles(
    as.matrix(samples_), probs = c(0.5, 0.16, 0.84)
  )
  df = data.frame(df)
  names(df) = c("Q50", "Q16", "Q84")
  df$ERR = 0.5 * (df$Q84 - df$Q16)
  df$x = 10^mdustvec
  return(df)
}
csmh_regressed = data.frame(foreach(i = 1:length(lbt_mids), .combine = rbind) %do% {
  csmh_samples = sapply(1:Nsamples, function(j){
    # M = Mstar_Evol_Samples[j,1] +     Mstar_Evol_Samples[j,2] * lbt_mids[i]
    # phi1 = phi1_Evol_Samples[j,1] +   phi1_Evol_Samples[j,2]  * lbt_mids[i]^7
    # phi2 = phi2_Evol_Samples[j,1] +   phi2_Evol_Samples[j,2]  * lbt_mids[i]
    # alpha = alpha_Evol_Samples[j,1] + alpha_Evol_Samples[j,2] * lbt_mids[i]
    # beta = beta_Evol_Samples[j,1] +   beta_Evol_Samples[j,2]  * lbt_mids[i]
    M = Mstar_Evol_Samples[j,1] +     Mstar_Evol_Samples[j,2] * zmids[i] + Mstar_Evol_Samples[j,3] * zmids[i]^2
    phi1 = phi1_Evol_Samples[j,1] +   phi1_Evol_Samples[j,2]  * zmids[i] + phi1_Evol_Samples[j,3]  * zmids[i]^2
    phi2 = phi2_Evol_Samples[j,1] +   phi2_Evol_Samples[j,2]  * zmids[i] + phi2_Evol_Samples[j,3]  * zmids[i]^2
    alpha = alpha_Evol_Samples[j,1] + alpha_Evol_Samples[j,2] * zmids[i] + alpha_Evol_Samples[j,3] * zmids[i]^2
    beta = beta_Evol_Samples[j,1] +   beta_Evol_Samples[j,2]  * zmids[i] + beta_Evol_Samples[j,3]  * zmids[i]^2

    regressed_Schechter = double_schechter(
      x = mdustvec,
      p = c(M, alpha, beta, phi1, phi2)
    )
    trapz(
      mdustvec,
      10^mdustvec * regressed_Schechter
    )
  })
  quantile(csmh_samples, c(0.50, 0.16, 0.84))
})
names(csmh_regressed) = c("Q50", "Q16", "Q84")
csmh_regressed$ERR = 0.5 * (csmh_regressed$Q84 - csmh_regressed$Q16)

LL_csmh = function(p){
  -1*sum(dnorm(
    x = log10(csmh_hybrid$Q50),
    mean = p[1] + log10(csmh_regressed$Q50),
    sd = csmh_hybrid$ERR/(log(10) * csmh_hybrid$Q50),
    log = TRUE
  ))
}
norm_csmh = optimise(f = LL_csmh, interval = c(-2,1))

# driver22_csmh = data.frame(
#   "z" = 0.0/2.0,
#   "lbt" = cosdistTravelTime(0.0/2.0, ref = "Planck18"),
#   "rho" = 2.97e8
# ) ## propagated to z=0
# magplot(
#   lbt_mids,
#   csmh_hybrid$Q50,
#   log = "y",
#   xlim = c(-0.5, 12),
#   ylim = c(1e7, 5e8),
#   xlab = "Lookback time",
#   ylab = "Stellar mass history"
# )
# magerr(
#   lbt_mids,
#   csmh_hybrid$Q50,
#   ylo = csmh_hybrid$ERR
# )
# lines(
#   lbt_mids, csmh_regressed$Q50, lty = 2
# )
# lines(
#   lbt_mids, csmh_regressed$Q50 * 10^norm_csmh$minimum,
#   lty = 1
# )
# points(
#   driver22_csmh$z,
#   driver22_csmh$rho,
#   pch = 16, col = "blue", cex = 1.5
# )
# legend(
#   x = "topright",
#   pch = c(16, 1, NA, NA),
#   col = c("blue", "black", "black", "black"),
#   lty = c(NA, NA, 2, 1),
#   legend = c("Driver+22", "ProHybrid", "Regressed fit", "Regressed fit scaled")
# )

LSS_corr = csmh_regressed$Q50 / csmh_hybrid$Q50 * 10^norm_csmh$minimum
LSS_corrErr = abs(LSS_corr) * sqrt((csmh_regressed$ERR/csmh_regressed$Q50)^2 + (csmh_hybrid$ERR/csmh_hybrid$Q50)^2) 
# LSS_corr[lbt_mids >= 7] = 1.0
# LSS_corrErr[lbt_mids >= 7] = 0.0

## Also show integrated M&D14 curve 
zvec = seq(0, 14, 0.01)
Hz = approxfun(
  zvec, 
  cosgrowH(z = zvec, ref = "Planck18")
)
psi = 0.015 * ((1+zvec)^2.7 / (1 + ((1+zvec)/2.9)^5.6))
make_mass = rev(
  cumtrapz(
    max(zvec) - rev(zvec), 
    rev(psi * ((Hz(zvec)*(1/3.09e19) / (1/3.15e7))*(1+zvec))^-1)
  ) 
) 
returnFrac = 0.41
ret = (1-returnFrac) * make_mass
md14_csmh = approxfun(
  zvec, 
  ret, 
  rule = 2
) 
norm_csmh_md14 = optimise(
  f = function(p){
    -1*sum(dnorm(
      x = log10(csmh_hybrid$Q50),
      mean = p[1] + log10(md14_csmh(zmids)),
      sd = csmh_hybrid$ERR/(log(10) * csmh_hybrid$Q50),
      log = TRUE
    ))
  },
  interval = c(-2,1)
)
LSS_md14 = md14_csmh(zmids) / csmh_hybrid$Q50 * 10^norm_csmh_md14$minimum

LSS = data.frame(
  "LSS" = LSS_corr,
  "ERR" = LSS_corrErr,
  "LSS_md14" = LSS_md14
)

# Md_corr = as.numeric(h5read(h5file, "Md_corr"))
cdmh_hybrid_corr$CORR = LSS_corr * cdmh_hybrid_corr$Q50
cdmh_hybrid_corr$CORR_ERR = (LSS_corr * cdmh_hybrid_corr$Q50) * sqrt( (cdmh_hybrid_corr$ERR/cdmh_hybrid_corr$Q50)^2 + (LSS_corrErr/LSS_corr)^2 )
cdmh_hybrid_corr$CORRQ16 = cdmh_hybrid_corr$CORR - cdmh_hybrid_corr$CORR_ERR
cdmh_hybrid_corr$CORRQ84 = cdmh_hybrid_corr$CORR + cdmh_hybrid_corr$CORR_ERR
cdmh_noAGN_corr$CORR = LSS_corr * cdmh_noAGN_corr$Q50
cdmh_noAGN_corr$CORR_ERR = (LSS_corr * cdmh_noAGN_corr$Q50) * sqrt( (cdmh_noAGN_corr$ERR/cdmh_noAGN_corr$Q50)^2 + (LSS_corrErr/LSS_corr)^2 )
cdmh_noAGN_corr$CORRQ16 = cdmh_noAGN_corr$CORR - cdmh_noAGN_corr$CORR_ERR
cdmh_noAGN_corr$CORRQ84 = cdmh_noAGN_corr$CORR + cdmh_noAGN_corr$CORR_ERR
cdmh_AGN_corr$CORR = LSS_corr * cdmh_AGN_corr$Q50
cdmh_AGN_corr$CORR_ERR = (LSS_corr * cdmh_AGN_corr$Q50) * sqrt( (cdmh_AGN_corr$ERR/cdmh_AGN_corr$Q50)^2 + (LSS_corrErr/LSS_corr)^2 )
cdmh_AGN_corr$CORRQ16 = cdmh_AGN_corr$CORR - cdmh_AGN_corr$CORR_ERR
cdmh_AGN_corr$CORRQ84 = cdmh_AGN_corr$CORR + cdmh_AGN_corr$CORR_ERR

cgmh_noAGN_corr$CORR = LSS_corr * cgmh_noAGN_corr$Q50
cgmh_noAGN_corr$CORR_ERR = (LSS_corr * cgmh_noAGN_corr$Q50) * sqrt( (cgmh_noAGN_corr$ERR/cgmh_noAGN_corr$Q50)^2 + (LSS_corrErr/LSS_corr)^2 )
cgmh_noAGN_corr$CORRQ16 = cgmh_noAGN_corr$CORR - cgmh_noAGN_corr$CORR_ERR
cgmh_noAGN_corr$CORRQ84 = cgmh_noAGN_corr$CORR + cgmh_noAGN_corr$CORR_ERR
cgmh_AGN_corr$CORR = LSS_corr * cgmh_AGN_corr$Q50
cgmh_AGN_corr$CORR_ERR = (LSS_corr * cgmh_AGN_corr$Q50) * sqrt( (cgmh_AGN_corr$ERR/cgmh_AGN_corr$Q50)^2 + (LSS_corrErr/LSS_corr)^2 )
cgmh_AGN_corr$CORRQ16 = cgmh_AGN_corr$CORR - cgmh_AGN_corr$CORR_ERR
cgmh_AGN_corr$CORRQ84 = cgmh_AGN_corr$CORR + cgmh_AGN_corr$CORR_ERR
cgmh_hybrid_corr$CORR = LSS_corr * cgmh_hybrid_corr$Q50
cgmh_hybrid_corr$CORR_ERR = (LSS_corr * cgmh_hybrid_corr$Q50) * sqrt( (cgmh_hybrid_corr$ERR/cgmh_hybrid_corr$Q50)^2 + (LSS_corrErr/LSS_corr)^2 )
cgmh_hybrid_corr$CORRQ16 = cgmh_hybrid_corr$CORR - cgmh_hybrid_corr$CORR_ERR
cgmh_hybrid_corr$CORRQ84 = cgmh_hybrid_corr$CORR + cgmh_hybrid_corr$CORR_ERR

csmh_hybrid$CORR = csmh_hybrid$Q50 * LSS_corr
csmh_hybrid$CORR_ERR = (LSS_corr * csmh_hybrid$Q50) * sqrt( (csmh_hybrid$ERR/csmh_hybrid$Q50)^2 + (LSS_corrErr/LSS_corr)^2 )
csmh_hybrid$CORRQ16 = csmh_hybrid$CORR - csmh_hybrid$CORR_ERR
csmh_hybrid$CORRQ84 = csmh_hybrid$CORR + csmh_hybrid$CORR_ERR

csmh_regressed$Q50_scaled = csmh_regressed$Q50 * 10^norm_csmh$minimum
csmh_regressed$Q16_scaled = csmh_regressed$Q16 * 10^norm_csmh$minimum
csmh_regressed$Q84_scaled = csmh_regressed$Q84 * 10^norm_csmh$minimum
csmh_regressed$ERR_scaled = csmh_regressed$ERR * 10^norm_csmh$minimum

# Compute difference on the cdmh between AGN and noAGN, in dex IDGAF
csmh_noAGN_AGN_difference = foreach(i = 1:dim(cdmh_AGN_corr)[1], .combine = bind_rows) %do% {
  noAGN_samples = rnorm(n = 1000, mean = csmh_noAGN$Q50[i], sd = csmh_noAGN$ERR[i])
  AGN_samples = rnorm(n = 1000, mean = csmh_AGN$Q50[i], sd = csmh_AGN$ERR[i])
  
  difference = log10(noAGN_samples/AGN_samples)
  list(
    "Q50" = median(difference, na.rm = TRUE), 
    "Q16" = quantile(difference, prob = 0.16, na.rm = TRUE),
    "Q84" = quantile(difference, prob = 0.84, na.rm = TRUE),
    "ERR" = 0.5 * (quantile(difference, prob = 0.84, na.rm = TRUE) - quantile(difference, prob = 0.16, na.rm = TRUE))
  )
}
cdmh_noAGN_AGN_difference = foreach(i = 1:dim(cdmh_AGN_corr)[1], .combine = bind_rows) %do% {
  noAGN_samples = rnorm(n = 1000, mean = cdmh_noAGN_corr$CORR[i], sd = cdmh_noAGN_corr$CORR_ERR[i])
  AGN_samples = rnorm(n = 1000, mean = cdmh_AGN_corr$CORR[i], sd = cdmh_AGN_corr$CORR_ERR[i])
  
  difference = log10(noAGN_samples/AGN_samples)
  list(
    "Q50" = median(difference, na.rm = TRUE), 
    "Q16" = quantile(difference, prob = 0.16, na.rm = TRUE),
    "Q84" = quantile(difference, prob = 0.84, na.rm = TRUE),
    "ERR" = 0.5 * (quantile(difference, prob = 0.84, na.rm = TRUE) - quantile(difference, prob = 0.16, na.rm = TRUE))
  )
}
cgmh_noAGN_AGN_difference = foreach(i = 1:dim(cdmh_AGN_corr)[1], .combine = bind_rows) %do% {
  noAGN_samples = rnorm(n = 1000, mean = cgmh_noAGN_corr$CORR[i], sd = cgmh_noAGN_corr$CORR_ERR[i])
  AGN_samples = rnorm(n = 1000, mean = cgmh_AGN_corr$CORR[i], sd = cgmh_AGN_corr$CORR_ERR[i])
  
  difference = log10(noAGN_samples/AGN_samples)
  list(
    "Q50" = median(difference, na.rm = TRUE), 
    "Q16" = quantile(difference, prob = 0.16, na.rm = TRUE),
    "Q84" = quantile(difference, prob = 0.84, na.rm = TRUE),
    "ERR" = 0.5 * (quantile(difference, prob = 0.84, na.rm = TRUE) - quantile(difference, prob = 0.16, na.rm = TRUE))
  )
}
# magplot(
#   lbt_mids,
#   cdmh_noAGN_AGN_difference$Q50, 
#   ylim = c(-0.1, 0.5)
# )
# magerr(
#   lbt_mids,
#   cdmh_noAGN_AGN_difference$Q50, 
#   ylo = cdmh_noAGN_AGN_difference$ERR
# )
# points(
#   lbt_mids,
#   cgmh_noAGN_AGN_difference$Q50, 
#   col = "red"
# )
# magerr(
#   lbt_mids,
#   cgmh_noAGN_AGN_difference$Q50, 
#   ylo = cgmh_noAGN_AGN_difference$ERR, 
#   col = "red"
# )
# 
# maghist(
#   log10(devilsd10_AGN$dustlum.total/devilsd10_noAGN$dustlum.total)
# )

h5ls(h5file)
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
h5write(obj = LSS, h5file, name = "LSSCorrection")

h5delete(h5file, "cosmic")
h5createGroup(h5file, "cosmic")
h5write(obj = cdmh_noAGN, file = h5file, name = "cosmic/MdustnoAGN")
h5write(obj = cdmh_AGN, file = h5file, name = "cosmic/MdustAGN")
h5write(obj = cdmh_hybrid, file = h5file, name = "cosmic/MdustHybrid")
h5write(obj = devilsd10_hybrid_cdmh_FIR, file = h5file, name = "cosmic/MdustDEVILSFIR")
h5write(obj = cdmh_noAGN_corr, file = h5file, name = "cosmic/MdustnoAGNCorr")
h5write(obj = cdmh_AGN_corr, file = h5file, name = "cosmic/MdustAGNCorr")
# h5delete(h5file, "cosmic/MdustHybridCorr")
h5write(obj = cdmh_hybrid_corr, file = h5file, name = "cosmic/MdustHybridCorr")
h5write(obj = cgmh_noAGN_corr, file = h5file, name = "cosmic/MgasnoAGNCorr")
h5write(obj = cgmh_AGN_corr, file = h5file, name = "cosmic/MgasAGNCorr")
h5write(obj = cgmh_hybrid_corr, file = h5file, name = "cosmic/MgasHybridCorr")
h5write(obj = csmh_hybrid, file = h5file, name = "cosmic/MstarHybrid")
h5write(obj = csmh_regressed, file = h5file, name = "cosmic/MstarRegressed")

h5write(obj = csmh_noAGN_AGN_difference, file = h5file, name = "cosmic/csmhAGNnoAGNDiff")
h5write(obj = cdmh_noAGN_AGN_difference, file = h5file, name = "cosmic/cdmhAGNnoAGNDiff")
h5write(obj = cgmh_noAGN_AGN_difference, file = h5file, name = "cosmic/cgmhAGNnoAGNDiff")

h5delete(h5file, "binDF")
h5createGroup(h5file, "binDF")
h5createGroup(h5file, "binDF/MdustnoAGN")
for(i in 1:length(combine_noAGN_dmf)){
  h5write(
    obj = combine_noAGN_dmf[[i]]$dmf, 
    file = h5file,
    name = paste0("binDF/MdustnoAGN/zb",i)
  )
}
h5createGroup(h5file, "binDF/MdustAGN")
for(i in 1:length(combine_noAGN_dmf)){
  h5write(
    obj = combine_AGN_dmf[[i]]$dmf, 
    file = h5file,
    name = paste0("binDF/MdustAGN/zb",i)
  )
}
h5createGroup(h5file, "binDF/MdustHybrid")
for(i in 1:length(combine_hybrid_dmf)){
  h5write(
    obj = combine_hybrid_dmf[[i]]$dmf, 
    file = h5file,
    name = paste0("binDF/MdustHybrid/zb",i)
  )
}
h5createGroup(h5file, "binDF/MdustHybridCorr")
for(i in 1:length(combine_hybrid_dmf_corr)){
  h5write(
    obj = combine_hybrid_dmf_corr[[i]]$dmf, 
    file = h5file,
    name = paste0("binDF/MdustHybridCorr/zb",i)
  )
}
h5createGroup(h5file, "binDF/MgasHybridCorr")
for(i in 1:length(combine_hybrid_gmf_corr)){
  h5write(
    obj = combine_hybrid_gmf_corr[[i]]$dmf, 
    file = h5file,
    name = paste0("binDF/MgasHybridCorr/zb",i)
  )
}
h5createGroup(h5file, "binDF/MstarHybrid")
for(i in 1:length(combine_hybrid_smf)){
  h5write(
    obj = combine_hybrid_smf[[i]]$dmf, 
    file = h5file,
    name = paste0("binDF/MstarHybrid/zb",i)
  )
}

## Gama and devils
h5createGroup(h5file, "binDF/MdustHybridGAMA")
for(i in 1:length(gama_hybrid_dmf_corr)){
  h5write(
    obj = gama_hybrid_dmf_corr[[i]], 
    file = h5file,
    name = paste0("binDF/MdustHybridGAMA/zb",i)
  )
}
h5createGroup(h5file, "binDF/MdustHybridDEVILSD10")
for(i in 1:length(devilsd10_hybrid_dmf_corr)){
  h5write(
    obj = devilsd10_hybrid_dmf_corr[[i]], 
    file = h5file,
    name = paste0("binDF/MdustHybridDEVILSD10/zb",i)
  )
}
h5createGroup(h5file, "binDF/MgasHybridGAMA")
for(i in 1:length(gama_hybrid_gmf_corr)){
  h5write(
    obj = gama_hybrid_gmf_corr[[i]], 
    file = h5file,
    name = paste0("binDF/MgasHybridGAMA/zb",i)
  )
}
h5createGroup(h5file, "binDF/MgasHybridDEVILSD10")
for(i in 1:length(devilsd10_hybrid_gmf_corr)){
  h5write(
    obj = devilsd10_hybrid_gmf_corr[[i]], 
    file = h5file,
    name = paste0("binDF/MgasHybridDEVILSD10/zb",i)
  )
}

h5delete(h5file, "fitDF")
h5createGroup(h5file, "fitDF")
h5createGroup(h5file, "fitDF/MdustnoAGN")
for(i in 1:length(combine_noAGN_dmf)){
  h5write(
    obj = combine_noAGN_dmf[[i]]$fit, 
    file = h5file,
    name = paste0("fitDF/MdustnoAGN/zb",i)
  )
}
h5createGroup(h5file, "fitDF/MdustAGN")
for(i in 1:length(combine_AGN_dmf)){
  h5write(
    obj = combine_AGN_dmf[[i]]$fit, 
    file = h5file,
    name = paste0("fitDF/MdustAGN/zb",i)
  )
}
h5createGroup(h5file, "fitDF/MdustHybrid")
for(i in 1:length(combine_hybrid_dmf)){
  h5write(
    obj = combine_hybrid_dmf[[i]]$fit, 
    file = h5file,
    name = paste0("fitDF/MdustHybrid/zb",i)
  )
}
h5createGroup(h5file, "fitDF/MdustHybridCorr")
for(i in 1:length(combine_hybrid_dmf_corr)){
  h5write(
    obj = combine_hybrid_dmf_corr[[i]]$fit, 
    file = h5file,
    name = paste0("fitDF/MdustHybridCorr/zb",i)
  )
}
h5createGroup(h5file, "fitDF/MgasHybridCorr")
for(i in 1:length(combine_hybrid_gmf_corr)){
  h5write(
    obj = combine_hybrid_gmf_corr[[i]]$fit, 
    file = h5file,
    name = paste0("fitDF/MgasHybridCorr/zb",i)
  )
}
h5createGroup(h5file, "fitDF/MstarHybrid")
for(i in 1:length(combine_hybrid_smf)){
  h5write(
    obj = combine_hybrid_smf[[i]]$fit, 
    file = h5file,
    name = paste0("fitDF/MstarHybrid/zb",i)
  )
}
h5createGroup(h5file, "fitDF/MstarHybridRegressed")
for(i in 1:length(smf_regressed_samples)){
  h5write(
    obj = smf_regressed_samples[[i]], 
    file = h5file,
    name = paste0("fitDF/MstarHybridRegressed/zb",i)
  )
}

h5delete(h5file, "par")
h5createGroup(h5file, "par")
h5write(
  obj = combine_noAGN_dmf_par, 
  file = h5file,
  name = "par/noAGNDMF"
)
h5write(
  obj = combine_AGN_dmf_par, 
  file = h5file,
  name = "par/AGNDMF"
)
h5write(
  obj = combine_hybrid_dmf_par, 
  file = h5file,
  name = "par/HybridDMF"
)
h5write(
  obj = combine_hybrid_dmf_corr_par, 
  file = h5file,
  name = "par/HybridCorrDMF"
)
h5write(
  obj = combine_hybrid_gmf_corr_par, 
  file = h5file,
  name = "par/HybridCorrGMF"
)
h5write(
  obj = combine_hybrid_smf_par, 
  file = h5file,
  name = "par/HybridSMF"
)
h5write(
  obj = smf_regressed_evol, 
  file = h5file,
  name = "par/HybridSMFEvol"
)
h5write(
  obj = smf_regressed_fits, 
  file = h5file,
  name = "par/HybridSMFEvolPar"
)

save.image("~/Documents/DustMassDensity/data/dmf.Rdata")