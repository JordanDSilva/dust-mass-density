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

catalogueDir = "/Users/22252335/Documents/GAMA-DEVILS-SFR-AGN/data/"

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

devilsd10_AGN = fread("~/Documents/GAMA-DEVILS-SFR-AGN/data/AGNTotalCat_MasterCat4.csv")
devilsd10_noAGN = readRDS(paste0(catalogueDir, 'DEVILS_D10ProSpectCat_02_02_2021_v0.3.rds')) #catalogue that Jess done without AGN contribution
devilsd10_noAGN = devilsd10_noAGN$cat
devilsd10_noAGN$area = 1.47
devilsd10_AGN$area = 1.47
devilsd10_noAGN = data.frame(devilsd10_noAGN[order(devilsd10_noAGN$UID),])
devilsd10_AGN = data.frame(devilsd10_AGN[order(devilsd10_AGN$UID),])

## make a hybrid catalogue
devilsd10_col_names = c(
  "UID", "z", 
  "StellarMass", "dustmass.birth","dustmass.screen","dustmass.total","dustlum.birth","dustlum.screen","dustlum.total","Zfinal",
  "StellarMass_LB", "dustmass.birth_LB","dustmass.screen_LB","dustmass.total_LB","dustlum.birth_LB","dustlum.screen_LB","dustlum.total_LB","Zfinal_LB",
  "StellarMass_UB", "dustmass.birth_UB","dustmass.screen_UB","dustmass.total_UB","dustlum.birth_UB","dustlum.screen_UB","dustlum.total_UB","Zfinal_UB"
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
gama_AGN_idx = gama_AGN$fAGN_bestfit >= 0.1 & gama_AGN$LP > gama_noAGN$LP
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
Nsamples = 1000

mc_err_gama_noAGN_dmf = foreach(i = 1:length(zmids)) %do% {
  message(i)
  vol = vol_mids[i] * 217.54 / (4*pi*(180/pi)^2)
  zidx = gama_noAGN$z >= zbins[i] & gama_noAGN$z < zbins[i+1]
  Mdust = gama_noAGN$DustMass_50[zidx]
  MdustErr = 0.5 * (gama_noAGN$DustMass_84[zidx] - gama_noAGN$DustMass_16[zidx])
  
  samples = matrix(rnorm(n = Nsamples*sum(zidx), mean = Mdust, sd = MdustErr), nrow = Nsamples, ncol = sum(zidx))
  
  hist_samples = foreach(j = 1:Nsamples, .combine = "rbind") %do% {
    m_sample = samples[j,]
    log_m_sample = log10(m_sample[m_sample > 0])
    hh = maghist(x = log_m_sample, breaks = sm_bins, plot = FALSE, verbose = FALSE)
    bin_dmf = hh$counts/(vol * diff(sm_bins))
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
  
  samples = matrix(rnorm(n = Nsamples*sum(zidx), mean = Mdust, sd = MdustErr), nrow = Nsamples, ncol = sum(zidx))
  
  hist_samples = foreach(j = 1:Nsamples, .combine = "rbind") %do% {
    m_sample = samples[j,]
    log_m_sample = log10(m_sample[m_sample > 0])
    hh = maghist(x = log_m_sample, breaks = sm_bins, plot = FALSE, verbose = FALSE)
    bin_dmf = hh$counts/(vol * diff(sm_bins))
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
  
  samples = matrix(rnorm(n = Nsamples*sum(zidx), mean = Mdust, sd = MdustErr), nrow = Nsamples, ncol = sum(zidx))
  
  hist_samples = foreach(j = 1:Nsamples, .combine = "rbind") %do% {
    m_sample = samples[j,]
    log_m_sample = log10(m_sample[m_sample > 0])
    hh = maghist(x = log_m_sample, breaks = sm_bins, plot = FALSE, verbose = FALSE)
    bin_dmf = hh$counts/(vol * diff(sm_bins))
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
  
  samples = matrix(rnorm(n = Nsamples*sum(zidx), mean = Mdust, sd = MdustErr), nrow = Nsamples, ncol = sum(zidx))
  
  hist_samples = foreach(j = 1:Nsamples, .combine = "rbind") %do% {
    m_sample = samples[j,]
    log_m_sample = log10(m_sample[m_sample > 0])
    hh = maghist(x = log_m_sample, breaks = sm_bins, plot = FALSE, verbose = FALSE)
    bin_dmf = hh$counts/(vol * diff(sm_bins))
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
  
  samples = matrix(rnorm(n = Nsamples*sum(zidx), mean = Mdust, sd = MdustErr), nrow = Nsamples, ncol = sum(zidx))
  
  hist_samples = foreach(j = 1:Nsamples, .combine = "rbind") %do% {
    m_sample = samples[j,]
    log_m_sample = log10(m_sample[m_sample > 0])
    hh = maghist(x = log_m_sample, breaks = sm_bins, plot = FALSE, verbose = FALSE)
    bin_dmf = hh$counts/(vol * diff(sm_bins))
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
  
  samples = matrix(rnorm(n = Nsamples*sum(zidx), mean = Mdust, sd = MdustErr), nrow = Nsamples, ncol = sum(zidx))
  
  hist_samples = foreach(j = 1:Nsamples, .combine = "rbind") %do% {
    m_sample = samples[j,]
    log_m_sample = log10(m_sample[m_sample > 0])
    hh = maghist(x = log_m_sample, breaks = sm_bins, plot = FALSE, verbose = FALSE)
    bin_dmf = hh$counts/(vol * diff(sm_bins))
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
  
  Z = gama_hybrid$Zgas_50[zidx]
  Zerr = 0.5 * (gama_hybrid$Zgas_84[zidx] - gama_hybrid$Zgas_16[zidx])
  
  samples = matrix(rnorm(n = Nsamples*sum(zidx), mean = Mdust, sd = MdustErr), nrow = Nsamples, ncol = sum(zidx))
  Z_samples = matrix(rnorm(n = Nsamples*sum(zidx), mean = log10(Z), sd = Zerr/(Z*log(10))), nrow = Nsamples, ncol = sum(zidx))
  
  hist_samples = foreach(j = 1:Nsamples, .combine = "rbind") %do% {
    m_sample = samples[j,] * RR14_BPL(Z = 10^Z_samples[j,])/0.0073
    log_m_sample = log10(m_sample[m_sample > 0])
    hh = maghist(x = log_m_sample, breaks = sm_bins, plot = FALSE, verbose = FALSE)
    bin_dmf = hh$counts/(vol * diff(sm_bins))
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
  Z = devilsd10_hybrid$Zfinal[zidx]
  Zerr = 0.5 * (devilsd10_hybrid$Zfinal_UB[zidx] - devilsd10_hybrid$Zfinal_LB[zidx])
  
  samples = matrix(rnorm(n = Nsamples*sum(zidx), mean = Mdust, sd = MdustErr), nrow = Nsamples, ncol = sum(zidx))
  Z_samples = matrix(rnorm(n = Nsamples*sum(zidx), mean = Z, sd = Zerr), nrow = Nsamples, ncol = sum(zidx))
  
  hist_samples = foreach(j = 1:Nsamples, .combine = "rbind") %do% {
    m_sample = samples[j,] * RR14_BPL(Z = 10^Z_samples[j,])/0.0073
    log_m_sample = log10(m_sample[m_sample > 0])
    hh = maghist(x = log_m_sample, breaks = sm_bins, plot = FALSE, verbose = FALSE)
    bin_dmf = hh$counts/(vol * diff(sm_bins))
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
  
  samples = matrix(rnorm(n = Nsamples*sum(zidx), mean = Mdust, sd = MdustErr), nrow = Nsamples, ncol = sum(zidx))
  
  hist_samples = foreach(j = 1:Nsamples, .combine = "rbind") %do% {
    m_sample = samples[j,]
    log_m_sample = log10(m_sample[m_sample > 0])
    hh = maghist(x = log_m_sample, breaks = sm_bins, plot = FALSE, verbose = FALSE)
    bin_dmf = hh$counts/(vol * diff(sm_bins))
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
  
  samples = matrix(rnorm(n = Nsamples*sum(zidx), mean = Mdust, sd = MdustErr), nrow = Nsamples, ncol = sum(zidx))
  
  hist_samples = foreach(j = 1:Nsamples, .combine = "rbind") %do% {
    m_sample = samples[j,]
    log_m_sample = log10(m_sample[m_sample > 0])
    hh = maghist(x = log_m_sample, breaks = sm_bins, plot = FALSE, verbose = FALSE)
    bin_dmf = hh$counts/(vol * diff(sm_bins))
    return(bin_dmf)
  }
  hist_quant = colQuantiles(hist_samples, probs = c(0.5, 0.16, 0.84))
}

err_floor = 0.1
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
  
  err_dmf = sqrt( pois_dmf^2 + (err_floor*bin_dmf)^2 + mc_err^2 )
  df = data.frame(cbind(bin_dmf, err_dmf))
  names(df) = c("phi", "err")
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
  
  err_dmf = sqrt( pois_dmf^2 + (err_floor*bin_dmf)^2 + mc_err^2 )
  df = data.frame(cbind(bin_dmf, err_dmf))
  names(df) = c("phi", "err")
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
  
  err_dmf = sqrt( pois_dmf^2 + (err_floor*bin_dmf)^2 + mc_err^2 )
  df = data.frame(cbind(bin_dmf, err_dmf))
  names(df) = c("phi", "err")
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
  
  err_dmf = sqrt( pois_dmf^2 + (err_floor*bin_dmf)^2 + mc_err^2 )
  df = data.frame(cbind(bin_dmf, err_dmf))
  names(df) = c("phi", "err")
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
  
  err_dmf = sqrt( pois_dmf^2 + (err_floor*bin_dmf)^2 + mc_err^2 )
  df = data.frame(cbind(bin_dmf, err_dmf))
  names(df) = c("phi", "err")
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
  
  err_dmf = sqrt( pois_dmf^2 + (err_floor*bin_dmf)^2 + mc_err^2 )
  df = data.frame(cbind(bin_dmf, err_dmf))
  names(df) = c("phi", "err")
  return(df)
}
gama_hybrid_dmf_corr = foreach(i = 1:length(zmids)) %do% {
  message(i)
  vol = vol_mids[i] * 217.54 / (4*pi*(180/pi)^2)
  zidx = gama_hybrid$z >= zbins[i] & gama_hybrid$z < zbins[i+1]
  Mdust = gama_hybrid$DustMass_50[zidx] * RR14_BPL(Z = gama_hybrid$Zgas_50[zidx], doDTG = FALSE)/0.0073
  mc_err = 0.5*(mc_err_gama_hybrid_dmf_corr[[i]][,3] - mc_err_gama_hybrid_dmf_corr[[i]][,2])
  
  log_m = log10(Mdust[Mdust > 0])
  hh = maghist(x = log_m, breaks = sm_bins, plot = FALSE, verbose = FALSE)
  bin_dmf = hh$counts/(vol * diff(sm_bins))
  pois_dmf = sqrt(hh$counts)/(vol * diff(sm_bins))
  
  err_dmf = sqrt( pois_dmf^2 + (err_floor*bin_dmf)^2 + mc_err^2 )
  df = data.frame(cbind(bin_dmf, err_dmf))
  names(df) = c("phi", "err")
  return(df)
}
devilsd10_hybrid_dmf_corr = foreach(i = 1:length(zmids)) %do% {
  message(i)
  vol = vol_mids[i] * 1.47 / (4*pi*(180/pi)^2)
  zidx = devilsd10_hybrid$z >= zbins[i] & devilsd10_hybrid$z < zbins[i+1]
  Mdust = devilsd10_hybrid$dustmass.total[zidx] * RR14_BPL(Z = 10^devilsd10_hybrid$Zfinal[zidx], doDTG = FALSE)/0.0073
  mc_err = 0.5*(mc_err_devilsd10_hybrid_dmf_corr[[i]][,3] - mc_err_devilsd10_hybrid_dmf_corr[[i]][,2])
  
  log_m = log10(Mdust[Mdust > 0])
  hh = maghist(x = log_m, breaks = sm_bins, plot = FALSE, verbose = FALSE)
  bin_dmf = hh$counts/(vol * diff(sm_bins))
  pois_dmf = sqrt(hh$counts)/(vol * diff(sm_bins))
  
  err_dmf = sqrt( pois_dmf^2 + (err_floor*bin_dmf)^2 + mc_err^2 )
  df = data.frame(cbind(bin_dmf, err_dmf))
  names(df) = c("phi", "err")
  return(df)
}
gama_hybrid_gmf_corr = foreach(i = 1:length(zmids)) %do% {
  message(i)
  vol = vol_mids[i] * 217.54 / (4*pi*(180/pi)^2)
  zidx = gama_hybrid$z >= zbins[i] & gama_hybrid$z < zbins[i+1]
  Mdust = gama_hybrid$DustMass_50[zidx] * RR14_BPL(Z = gama_hybrid$Zgas_50[zidx], doDTG = FALSE)/0.0073/RR14_BPL(Z = gama_hybrid$Zgas_50[zidx], doDTG = TRUE)
  mc_err = 0.5*(mc_err_gama_hybrid_dmf_corr[[i]][,3] - mc_err_gama_hybrid_dmf_corr[[i]][,2])
  
  log_m = log10(Mdust[Mdust > 0])
  hh = maghist(x = log_m, breaks = sm_bins, plot = FALSE, verbose = FALSE)
  bin_dmf = hh$counts/(vol * diff(sm_bins))
  pois_dmf = sqrt(hh$counts)/(vol * diff(sm_bins))
  
  err_dmf = sqrt( pois_dmf^2 + (err_floor*bin_dmf)^2 + mc_err^2 )
  df = data.frame(cbind(bin_dmf, err_dmf))
  names(df) = c("phi", "err")
  return(df)
}
devilsd10_hybrid_gmf_corr = foreach(i = 1:length(zmids)) %do% {
  message(i)
  vol = vol_mids[i] * 1.47 / (4*pi*(180/pi)^2)
  zidx = devilsd10_hybrid$z >= zbins[i] & devilsd10_hybrid$z < zbins[i+1]
  Mdust = devilsd10_hybrid$dustmass.total[zidx] * RR14_BPL(Z = 10^devilsd10_hybrid$Zfinal[zidx], doDTG = FALSE)/0.0073/RR14_BPL(Z = 10^devilsd10_hybrid$Zfinal[zidx], doDTG = TRUE)
  mc_err = 0.5*(mc_err_devilsd10_hybrid_dmf_corr[[i]][,3] - mc_err_devilsd10_hybrid_dmf_corr[[i]][,2])
  
  log_m = log10(Mdust[Mdust > 0])
  hh = maghist(x = log_m, breaks = sm_bins, plot = FALSE, verbose = FALSE)
  bin_dmf = hh$counts/(vol * diff(sm_bins))
  pois_dmf = sqrt(hh$counts)/(vol * diff(sm_bins))
  
  err_dmf = sqrt( pois_dmf^2 + (err_floor*bin_dmf)^2 + mc_err^2 )
  df = data.frame(cbind(bin_dmf, err_dmf))
  names(df) = c("phi", "err")
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
  
  err_dmf = sqrt( pois_dmf^2 + (err_floor*bin_dmf)^2 + mc_err^2 )
  df = data.frame(cbind(bin_dmf, err_dmf))
  names(df) = c("phi", "err")
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
  
  err_dmf = sqrt( pois_dmf^2 + (err_floor*bin_dmf)^2 + mc_err^2 )
  df = data.frame(cbind(bin_dmf, err_dmf))
  names(df) = c("phi", "err")
  return(df)
}

## Make cosmic 
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
    lower = c(5, -2.5, -2.5, -8, -8),
    upper = c(15, 2.5, 2.5, 1, 1),
    Data = list(
      xx = sm_mids[combine_dmf[,1] > 0], 
      yy = combine_dmf[,1][combine_dmf[,1] > 0], 
      yyerr = combine_dmf[,2][combine_dmf[,1] > 0], 
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
    )
  )
  df_cosmic$ERR = 0.5*(df_cosmic$Q84 - df_cosmic$Q16)
  
  combine_dmf = data.frame(combine_dmf)
  names(combine_dmf) = c("phi", "err")
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
    lower = c(5, -2.5, -2.5, -8, -8),
    upper = c(15, 2.5, 2.5, 1, 1),
    Data = list(
      xx = sm_mids[combine_dmf[,1] > 0], 
      yy = combine_dmf[,1][combine_dmf[,1] > 0], 
      yyerr = combine_dmf[,2][combine_dmf[,1] > 0], 
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
    )
  )
  df_cosmic$ERR = 0.5*(df_cosmic$Q84 - df_cosmic$Q16)
  
  combine_dmf = data.frame(combine_dmf)
  names(combine_dmf) = c("phi", "err")
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
    lower = c(5, -2.5, -2.5, -8, -8),
    upper = c(15, 2.5, 2.5, 1, 1),
    Data = list(
      xx = sm_mids[combine_dmf[,1] > 0], 
      yy = combine_dmf[,1][combine_dmf[,1] > 0], 
      yyerr = combine_dmf[,2][combine_dmf[,1] > 0], 
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
    )
  )
  df_cosmic$ERR = 0.5*(df_cosmic$Q84 - df_cosmic$Q16)
  
  combine_dmf = data.frame(combine_dmf)
  names(combine_dmf) = c("phi", "err")
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
combine_hybrid_dmf_corr = foreach(i = 1:length(zmids)) %do% {
  
  gama_dmf = gama_hybrid_dmf_corr[[i]]
  gama_err = gama_dmf[,2]
  
  devilsd10_dmf = devilsd10_hybrid_dmf_corr[[i]]
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
    lower = c(5, -2.5, -2.5, -8, -8),
    upper = c(15, 2.5, 2.5, 1, 1),
    Data = list(
      xx = sm_mids[combine_dmf[,1] > 0], 
      yy = combine_dmf[,1][combine_dmf[,1] > 0], 
      yyerr = combine_dmf[,2][combine_dmf[,1] > 0], 
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
    )
  )
  df_cosmic$ERR = 0.5*(df_cosmic$Q84 - df_cosmic$Q16)
  
  combine_dmf = data.frame(combine_dmf)
  names(combine_dmf) = c("phi", "err")
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
    lower = c(5, -2.5, -2.5, -8, -8),
    upper = c(15, 2.5, 2.5, 1, 1),
    Data = list(
      xx = sm_mids[combine_dmf[,1] > 0], 
      yy = combine_dmf[,1][combine_dmf[,1] > 0], 
      yyerr = combine_dmf[,2][combine_dmf[,1] > 0], 
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
    )
  )
  df_cosmic$ERR = 0.5*(df_cosmic$Q84 - df_cosmic$Q16)
  
  combine_dmf = data.frame(combine_dmf)
  names(combine_dmf) = c("phi", "err")
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
combine_hybrid_smf = foreach(i = 1:length(zmids)) %do% {
  
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
    lower = c(5, -2.5, -2.5, -8, -8),
    upper = c(15, 2.5, 2.5, 1, 1),
    Data = list(
      xx = sm_mids[combine_dmf[,1] > 0], 
      yy = combine_dmf[,1][combine_dmf[,1] > 0], 
      yyerr = combine_dmf[,2][combine_dmf[,1] > 0], 
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
    )
  )
  df_cosmic$ERR = 0.5*(df_cosmic$Q84 - df_cosmic$Q16)
  
  combine_dmf = data.frame(combine_dmf)
  names(combine_dmf) = c("phi", "err")
  ret_ = list(
    "dmf" = combine_dmf,
    "fit" = df_fit, 
    "cosmic" = df_cosmic,
    "highout" = highout
  )
  
  png(paste0("~/Documents/DustMassDensity/plots/smf_hybrid/lbt_", round(lbt_mids[i],3), ".png"))
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
  dev.off()
  return(ret_)
}

cdmh_noAGN = data.frame(foreach(i = 1:length(combine_noAGN_dmf), .combine = bind_rows) %do% {
  combine_noAGN_dmf[[i]]$cosmic
})
cdmh_AGN = data.frame(foreach(i = 1:length(combine_AGN_dmf), .combine = bind_rows) %do% {
  combine_AGN_dmf[[i]]$cosmic
})
cdmh_hybrid = data.frame(foreach(i = 1:length(combine_hybrid_dmf), .combine = bind_rows) %do% {
  combine_hybrid_dmf[[i]]$cosmic
})
cdmh_hybrid_corr = data.frame(foreach(i = 1:length(combine_hybrid_dmf_corr), .combine = bind_rows) %do% {
  combine_hybrid_dmf_corr[[i]]$cosmic
})
cgmh_hybrid_corr = data.frame(foreach(i = 1:length(combine_hybrid_gmf_corr), .combine = bind_rows) %do% {
  combine_hybrid_gmf_corr[[i]]$cosmic
})
csmh_hybrid = data.frame(foreach(i = 1:length(combine_hybrid_smf), .combine = bind_rows) %do% {
  combine_hybrid_smf[[i]]$cosmic
})

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

dsilva25 = data.frame(fread("~/Documents/DustMassDensity/data/literature_evo/csfh/DSilva25_CSFH_CAGNH_fit.csv"))
dsilva_csmh_func = approxfun(
  dsilva25$z, 
  10^dsilva25$CSMHQ50
)
LL_csmh = function(p){
  
  -1*sum(dnorm(
    x = log10( csmh_hybrid$Q50 ),
    mean = p[1] + log10(dsilva_csmh_func(zmids)),
    sd = csmh_hybrid$ERR/(log(10) * csmh_hybrid$Q50), 
    log = TRUE
  ))
  
}
norm_csmh = optimise(f = LL_csmh, interval = c(-2,1))
LSS_corr =dsilva_csmh_func(zmids) / csmh_hybrid$Q50 * 10^norm_csmh$minimum
LSS_corr[lbt_mids >= 8] = 1.0

h5file = '~/Documents/DustMassDensity/data/all_data.h5'
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
h5write(obj = LSS_corr, h5file, name = "LSSCorrection")

h5delete(h5file, "cosmic")
h5createGroup(h5file, "cosmic")
h5write(obj = cdmh_noAGN, file = h5file, name = "cosmic/MdustnoAGN")
h5write(obj = cdmh_AGN, file = h5file, name = "cosmic/MdustAGN")
h5write(obj = cdmh_hybrid, file = h5file, name = "cosmic/MdustHybrid")
h5write(obj = cdmh_hybrid_corr, file = h5file, name = "cosmic/MdustHybridCorr")
h5write(obj = cgmh_hybrid_corr, file = h5file, name = "cosmic/MgasHybridCorr")
h5write(obj = csmh_hybrid, file = h5file, name = "cosmic/MstarHybrid")

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

# save.image("~/Documents/DustMassDensity/data/dmf.Rdata")
