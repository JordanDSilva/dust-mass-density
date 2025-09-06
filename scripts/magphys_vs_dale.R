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
library(RColorBrewer)

# install.packages("extrafont")
# library(extrafont)
# font_import("~/Documents/dejavu-sans/")

#loadfonts()

set.seed(666)
## Try to implement the Magphys dust model

waveout = 10^seq(1, 9.35, 0.001)
totSED = ProSpectSED(z = 0, ref = "Planck18", waveout = log10(waveout), Dale = Dale_NormTot, Dale_M2L_func = Dale_M2L_func)
totSED_prior = ProSpectSED(z = 0, ref = "Planck18", waveout = log10(waveout), Dale = Dale_NormTot, Dale_M2L_func = Dale_M2L_func, alpha_SF_birth = 2, alpha_SF_screen = 2, tau_birth = 10^-0.2, tau_screen = 10^-2.3)
# plot(totSED)

## PAH emission template from Madden+06 of a photodissociation region in the MW
m17_pdr_pah = read.table("~/Documents/DustMassDensity/data/M17_PDR_PAH.csv", header = TRUE, sep = ",")
names(m17_pdr_pah) = c("mu", "flux")
m17_approx = cubicspline(x = c(m17_pdr_pah$mu*1e4), y = c(m17_pdr_pah$flux), xi = waveout, endp2nd = TRUE, der = c(1.721362e-06, -6.885449e-06))
m17_approx[m17_approx < 0] = 0
# m17_approx = approx(x = c(m17_pdr_pah$mu*1e4), y = c(m17_pdr_pah$flux), xout = waveout, yleft = 0, yright = 0)$y
Florentzian = (1 + ((waveout/1e8)^-1 - 3039.1)^2/(19.4)^2)^-1 * approx(x = waveout/1e4, y = m17_approx, xout = 11.3)$y * 0.1

# magplot(
#   waveout/1e4,
#   m17_approx,
#   log = "xy",
#   type = "l",
#   xlim = c(1, 20),
#   ylim = c(1e-6, 1)
# )
# points(
#   m17_pdr_pah$mu, 
#   m17_pdr_pah$flux,
#   pch = 16, col = "red"
# )
# lines(
#   waveout/1e4, Florentzian, col = "blue"
# )
# lines(
#   waveout/1e4, Florentzian + m17_approx, col = "green"
# )

PAH_emission_lines = (Florentzian + m17_approx) / waveout^2 ## Flambda
PAH_emission_lines = PAH_emission_lines / trapz(x = waveout, y = PAH_emission_lines)

grey_continuum = greybody_norm(
  wave = waveout, Temp = 850, beta = 1.0
)
G44 = approx(waveout/1e4, grey_continuum, 4.4)$y ## Continuum level at 4.4 micron

PAH77 = mean(approx(waveout/1e4, PAH_emission_lines, seq(7.5-0.3,7.5+0.3,0.1))$y) * 0.11

grey_continuum_scaled = grey_continuum * PAH77/G44

PAH_total = PAH_emission_lines + grey_continuum_scaled
PAH_total = PAH_total / trapz(waveout, PAH_total)

m17_norm = approxfun(x = waveout, y = PAH_total, yleft = 0, yright = 0)
m17_norm = splinefun(x = waveout, y = PAH_total)

# magplot(
#   waveout/1e4, 
#   PAH_emission_lines * waveout, 
#   log = "xy",
#   type = "l",
#   xlim = c(1, 20),
#   col = "blue"
#   # ylim= c(1e-6, 1)
# )
# lines(
#   waveout/1e4, grey_continuum_scaled*waveout, col = "red"
# )
# lines(
#   waveout/1e4, waveout*m17_norm(waveout), col = "black"
# )

l_MIR = function(wave){
  2.000002^-1 * (greybody_norm(wave, Temp = 130, beta = 1.5) + greybody_norm(wave, Temp = 250, beta = 1.5))
}

## Note that the greybody calculation for the dust mass already has 4pi factor in it as per the ProSpect equations...
L_BC = function(wave, p){
  L_dust_tot = p['L_dust_tot']
  f_mu_frac = p['f_mu_frac']
  xi_PAH_BC = p['xi_PAH_BC']
  xi_MIR_BC = p['xi_MIR_BC']
  xi_W_BC = p['xi_W_BC']
  Tw_BC = p['Tw_BC']
  # Tw_BC = 45
  df = list(
    "PAH" = xi_PAH_BC * m17_norm(wave) * (1 - f_mu_frac) * L_dust_tot,
    "MIR" = xi_MIR_BC * l_MIR(wave) * (1 - f_mu_frac) * L_dust_tot,
    "WBC" = xi_W_BC * greybody_norm(wave = wave, Temp = Tw_BC, beta = 1.5) * (1 - f_mu_frac) * L_dust_tot,
    "Ltot" = (xi_PAH_BC * m17_norm(wave) + xi_MIR_BC * l_MIR(wave) + xi_W_BC*greybody_norm(wave = wave, Temp = Tw_BC, beta = 1.5)) * (1 - f_mu_frac) * L_dust_tot,
    "MdustWBC" = xi_W_BC * (1 - f_mu_frac) * L_dust_tot * (trapz(wave, greybody(wave = wave, Temp = Tw_BC, beta = 1.5)))^-1
  )
  return(df)
}
L_ISM = function(wave, p){
  L_dust_tot = p['L_dust_tot']
  f_mu_frac = p['f_mu_frac']
  # xi_PAH_ISM = p['xi_PAH_ISM']
  # xi_MIR_ISM = p['xi_MIR_ISM']
  # xi_W_ISM = p['xi_W_ISM']
  xi_C_ISM = p['xi_C_ISM']
  xi_PAH_ISM = 0.550 * (1 - xi_C_ISM)
  xi_MIR_ISM = 0.275 * (1 - xi_C_ISM)
  xi_W_ISM = 0.175 * (1 - xi_C_ISM)
  # Tw_ISM = p['Tw_ISM']
  Tw_ISM = 45
  Tc_ISM = p['Tc_ISM']
  df = list(
    "PAH" = xi_PAH_ISM * m17_norm(wave) * (f_mu_frac) * L_dust_tot,
    "MIR" = xi_MIR_ISM * l_MIR(wave) * (f_mu_frac) * L_dust_tot,
    "WISM" = xi_W_ISM * greybody_norm(wave = wave, Temp = Tw_ISM, beta = 1.5) * (f_mu_frac) * L_dust_tot,
    "CISM" = xi_C_ISM * greybody_norm(wave = wave, Temp = Tc_ISM, beta = 2.0) * (f_mu_frac) * L_dust_tot,
    "Ltot" = (xi_PAH_ISM*m17_norm(wave) + xi_MIR_ISM*l_MIR(wave) + xi_W_ISM*greybody_norm(wave = wave, Temp = Tw_ISM, beta = 1.5) + xi_C_ISM*greybody_norm(wave = wave, Temp = Tc_ISM, beta = 2.0)) * (f_mu_frac) * L_dust_tot,
    "MdustWISM" = xi_W_ISM * f_mu_frac * L_dust_tot * (trapz(wave, greybody(wave = wave, Temp = Tw_ISM, beta = 1.5)))^-1,
    "MdustCISM" = xi_C_ISM * f_mu_frac * L_dust_tot * (trapz(wave, greybody(wave = wave, Temp = Tc_ISM, beta = 2.0)))^-1
  )
  return(df)
}
L_dust = function(wave, p){
  
  LBC = L_BC(wave, p)
  LISM = L_ISM(wave, p)
  
  df = list(
    "wave" = wave,
    "BC" = LBC,
    "ISM" = LISM,
    "Ltot" = LBC$Ltot + LISM$Ltot, 
    "Mdust" = unname(1.1 * (LBC$MdustWBC + LISM$MdustWISM + LISM$MdustCISM))
  )
  return(df)
}

## standardd model in Magphys paper
p = c(
  "L_dust_tot" = unname(totSED$dustlum["total"]) * 1.0,
  # "L_dust_tot" = 1,
  "f_mu_frac" = 0.6,
  
  "xi_PAH_BC" = 0.05,
  "xi_MIR_BC" = 0.15,
  "xi_W_BC" = 0.8,

  "xi_C_ISM" = 0.80,
  
  "Tw_BC" = 48,
  "Tc_ISM" = 22
)
names_magphys = names(p)

## Draw parameters from fits?
# magphys_cat = Rfits_read_table("~/Documents/DustMassDensity/data/GAMA/MagPhysv06.fits")
# devilsd10_noAGN = readRDS(paste0(catalogueDir, 'DEVILS_D10ProSpectCat_02_02_2021_v0.3.rds')) #catalogue that Jess done without AGN contribution
# devilsd10_noAGN = devilsd10_noAGN$cat

sample_from_distribution = function(Niters, data, xmin, xmax){
  hh = maghist(
    x = data, 
    xlim = c(xmin, xmax), 
    freq = FALSE
  )
  
  draws = approx(
    x = cumtrapz(
      x = hh$mids, 
      y = hh$density
    ),
    y = hh$mids, 
    xout = runif(n = Niters, min = 0.0, max = 1.0), 
    rule = 2
  )$y
  
  maghist(
    draws, freq = FALSE, breaks = hh$breaks, col = alpha("red", 0.4), add = TRUE
  )
  return(draws)
}

Niters = 1000
# alpha_birth_samples = pmin(pmax(rnorm(Niters, mean = 2, sd = 1), 0) ,4)
# alpha_screen_samples = pmin(pmax(rnorm(Niters, mean = 2, sd = 1), 0) ,4)
# tau_birth_samples = pmin(pmax(rnorm(Niters, mean = -0.2, sd = 0.5), -2.5), 1.5)
# tau_screen_prior = approxfun(
#   cumtrapz(
#     x = seq(-5.0, 1.0, 0.01),
#     y = -20 * erf(seq(-5.0, 1.0, 0.01)-2)
#   )/118.9949,
#   seq(-5.0, 1.0, 0.01)
# )
# tau_screen_samples = pmin(pmax(tau_screen_prior(runif(n = Niters, min = 0.0, max = 1.0)),-5.0), 1.0)

alpha_birth_samples = runif(Niters, 0.0, 4.0)
alpha_screen_samples = runif(Niters, 0.0, 4.0)
tau_birth_samples = runif(Niters, -2.5, 1.5)
tau_screen_samples = runif(Niters, -5.0, 1.0)
prospect_dale_samples = foreach(i = 1:Niters) %do% {
  if(i %% 100 == 0){
    message(i)
  }
  testSED = ProSpectSED(
    z = 0,
    ref = "Planck18",
    waveout = log10(waveout),
    Dale = Dale_NormTot,
    Dale_M2L_func = Dale_M2L_func,
    alpha_SF_birth = alpha_birth_samples[i],
    alpha_SF_screen = alpha_screen_samples[i],
    tau_birth = 10^tau_birth_samples[i],
    tau_screen = 10^tau_screen_samples[i]
  )
  
  df = list(
    "wave" = testSED$DustEmit$wave,
    "lum" = testSED$DustEmit$lum,
    "dustlum" = testSED$dustlum,
    "dustmass" = testSED$dustmass
  )
  return(df)
}
dale_dust_lum_samples = sapply(prospect_dale_samples, function(x)unname(x$dustlum["total"]))
dale_dust_lum_birth_samples = sapply(prospect_dale_samples, function(x)unname(x$dustlum["birth"]))
dale_dust_lum_screen_samples = sapply(prospect_dale_samples, function(x)unname(x$dustlum["screen"]))
dale_dust_mass_samples = sapply(prospect_dale_samples, function(x)unname(x$dustmass['total']))
dale_dust_mass_birth_samples = sapply(prospect_dale_samples, function(x)unname(x$dustmass['birth']))
dale_dust_mass_screen_samples = sapply(prospect_dale_samples, function(x)unname(x$dustmass['screen']))
dale_FIR_samples = sapply(prospect_dale_samples, function(x)trapz(x$wave[x$wave>2e5], x$lum[x$wave>2e5]))

## Draw parameters as per the magphys paper to construct the 50k models library
magphys_parm_mat = matrix(data = 0L, nrow = length(p)-1, ncol = Niters)
for(i in 1:Niters){
  
  frac_BC = runif(1, 0, 1)
  
  xi_W_BC_sample = runif(1, 0, 1.0)
  xi_MIR_BC_sample = runif(1, 0, 1.0-xi_W_BC_sample)
  xi_PAH_BC_sample = 1.0 - xi_W_BC_sample - xi_MIR_BC_sample

  magphys_parm_mat[,i] = c(
    runif(1, 0, 1),
    xi_PAH_BC_sample,
    xi_MIR_BC_sample,
    xi_W_BC_sample,
    runif(1, 0.5, 1),
    runif(1, 30, 60),
    runif(1, 15, 25)
  )
}
magphys_dust_lum = foreach(i = 1:Niters) %do% {
  if(i %% 100 == 0){message(i)}
  magphys_parm = c(dale_dust_lum_samples[i], magphys_parm_mat[,i])
  names(magphys_parm) = names(p)
  magphys_lum = L_dust(totSED$DustEmit$wave, p = magphys_parm)
}
magphys_dust_lum_samples = sapply(magphys_dust_lum, function(x)trapz(x$wave, x$Ltot))
magphys_dust_lum_Mass_samples = sapply(magphys_dust_lum, function(x){trapz(x$wave, x$BC$WBC) + trapz(x$wave, x$ISM$WISM) + trapz(x$wave, x$ISM$CISM)})
magphys_dust_lum_MassBC_samples = sapply(magphys_dust_lum, function(x){trapz(x$wave, x$BC$WBC)})
magphys_dust_lum_BCtot_samples = sapply(magphys_dust_lum, function(x){trapz(x$wave, x$BC$Ltot)})
magphys_dust_lum_MassISM_samples = sapply(magphys_dust_lum, function(x){trapz(x$wave, x$ISM$WISM) + trapz(x$wave, x$ISM$CISM)})
magphys_dust_lum_ISMtot_samples = sapply(magphys_dust_lum, function(x){trapz(x$wave, x$ISM$Ltot)})
magphys_dust_mass_samples = sapply(magphys_dust_lum,  function(x)unname(x$Mdust))
magphys_FIR_samples = sapply(magphys_dust_lum, function(x)trapz(x$wave[x$wave>2e5], x$Ltot[x$wave>2e5]))

magphys_standard =  L_dust(totSED$DustEmit$wave, p = p)

magphys_Mass_contribution = approxfun(
  totSED$DustEmit$wave,
  (magphys_standard$BC$WBC+magphys_standard$ISM$WISM+magphys_standard$ISM$CISM)/magphys_standard$Ltot,
  rule = 2
)
Dale_M2L_new = sapply(1:64, function(i){
  Msol = 1.989e30
  mH = 1.674e-27
  Lsol = 3.828e26
  DTH = 0.0073
  
  qPAH_VSG = 0.10
  message("I am using q = ", qPAH_VSG)
  qBIG = 1 - qPAH_VSG
  
  # weight = ifelse(
  #   Dale_Orig$Wave > 2e5, 
  #   1 - qPAH_VSG,
  #   qPAH_VSG
  # )
  
  weight = magphys_Mass_contribution(Dale_Orig$Wave) * (qBIG - qPAH_VSG) + qPAH_VSG
  weight[Dale_Orig$Wave >  1e7] = qBIG

  foo = ( (Dale_Orig$Aspec[[1]][i, ] / Lsol) / (DTH*weight * mH/Msol))/Dale_Orig$Wave
  sum(
    c(0, diff(Dale_Orig$Wave)) * foo
  )
})
Dale_vM2L_func = approxfun(
  Dale_M2L$alpha_SF, Dale_M2L_new, rule = 2
)
new_M2L_data = list(
  "magphys_standard" = magphys_standard, 
  "magphys_Mass_contribution" = magphys_Mass_contribution,
  "Dale_M2L_new" = Dale_M2L_new,
  "Dale_vM2L_func" = Dale_vM2L_func
)

new_M2L_DF = data.frame(
  "wave" = Dale_Orig$Wave,
  "weight" = magphys_Mass_contribution(Dale_Orig$Wave)
)
saveRDS(new_M2L_data, "~/Documents/DustMassDensity/data/new_M2L_data.rds")
dale_mass_birth_corr = dale_dust_lum_birth_samples/Dale_vM2L_func(alpha_birth_samples)
dale_mass_screen_corr = dale_dust_lum_screen_samples/Dale_vM2L_func(alpha_screen_samples)
dale_mass = dale_mass_screen_corr + dale_mass_birth_corr

Md_corr = median(dale_dust_mass_samples/dale_mass)
message("Average factor of difference? = ", Md_corr, sep = " ")

dale_spec_samples = as.matrix(foreach(i = 1:Niters, .combine = rbind) %do% {
  prospect_dale_samples[[i]]$lum 
})
magphys_spec_samples = as.matrix(foreach(i = 1:Niters, .combine = rbind) %do% {
  magphys_dust_lum[[i]]$Ltot 
})

df_samples = data.frame(
  "DaleLum" = dale_dust_lum_samples,
  "DaleLumBC" = dale_dust_lum_birth_samples,
  "DaleLumISM" = dale_dust_lum_screen_samples,
  "DaleMdust" = dale_dust_mass_samples,
  "DaleMdustBC" = dale_dust_mass_birth_samples,
  "DaleMdustISM" = dale_dust_mass_screen_samples,
  "DaleFIR" = dale_FIR_samples,
  "MagphysLum" = magphys_dust_lum_samples,
  "MagphysLumMass" = magphys_dust_lum_Mass_samples,
  "MagphysMdust" = magphys_dust_mass_samples,
  "MagphysFIR" = magphys_FIR_samples
)

spec_samples = data.frame(
  "wave" = unname(prospect_dale_samples[[1]]$wave), 
  "MagphysQ50" = colQuantiles(magphys_spec_samples, probs = 0.5),
  "MagphysQ05" = colQuantiles(magphys_spec_samples, probs = 0.05),
  "MagphysQ16" = colQuantiles(magphys_spec_samples, probs = 0.16),
  "MagphysQ84" = colQuantiles(magphys_spec_samples, probs = 0.84),
  "MagphysQ95" = colQuantiles(magphys_spec_samples, probs = 0.95),
  "DaleQ50" = colQuantiles(dale_spec_samples, probs = 0.5),
  "DaleQ05" = colQuantiles(dale_spec_samples, probs = 0.05),
  "DaleQ16" = colQuantiles(dale_spec_samples, probs = 0.16),
  "DaleQ84" = colQuantiles(dale_spec_samples, probs = 0.84),
  "DaleQ95" = colQuantiles(dale_spec_samples, probs = 0.95)
)

spec_standard = data.frame(
  "wave" = totSED$DustEmit$wave,
  "DaleDefault" = totSED$DustEmit$lum,
  "Dale" = totSED_prior$DustEmit$lum,
  "Magphys" = magphys_standard$Ltot,
  "MagphysBC" = magphys_standard$BC$Ltot,
  "MagphysISM" = magphys_standard$ISM$Ltot,
  "MagphysBC_PAH" = magphys_standard$BC$PAH,
  "MagphysBC_MIR" = magphys_standard$BC$MIR,
  "MagphysBC_WBC" = magphys_standard$BC$WBC,
  "MagphysISM_PAH" = magphys_standard$ISM$PAH,
  "MagphysISM_MIR" = magphys_standard$ISM$MIR,
  "MagphysISM_WISM" = magphys_standard$ISM$WISM,
  "MagphysISM_CISM" = magphys_standard$ISM$CISM
)

h5file = '~/Documents/DustMassDensity/data/all_data.h5'
h5delete(h5file, "DaleMagphys")
h5createGroup(h5file, "DaleMagphys")
h5delete(h5file, "Md_corr")
h5write(obj = unname(Md_corr), file = h5file, name = 'Md_corr')
h5write(
  obj = df_samples, file = h5file, name = "DaleMagphys/Samples"
)
h5write(
  obj = spec_samples, file = h5file, name = "DaleMagphys/Spec"
)
h5write(
  obj = spec_standard, file = h5file, name = "DaleMagphys/StandardSpec"
)
# h5delete(h5file, "DaleMagphys/DTHWeight")
h5write(
  obj = new_M2L_DF, file = h5file, name = "DaleMagphys/DTHWeight"
)

simonMagphys = h5read(file = h5file, name = "Photometry/FinalSample")

simonCompare = foreach(kk = 1:dim(simonMagphys)[1], .combine = rbind) %do% {
  
  ## See equations 19-22 in daCunha+08
  fmu = simonMagphys$f_mu_IR_best_fit[kk]
  xi_C = simonMagphys$xi_C_tot_best_fit[kk]/fmu
  xi_W = (simonMagphys$xi_W_tot_best_fit[kk] - 0.175*(1 - xi_C)*fmu)/(1-fmu)
  xi_MIR = (simonMagphys$xi_MIR_tot_best_fit[kk] - 0.275*(1 - xi_C)*fmu)/(1-fmu)
  xi_PAH = (simonMagphys$xi_PAH_tot_best_fit[kk] - 0.550*(1 - xi_C)*fmu)/(1-fmu)
  
  ptest = c(
    simonMagphys$L_dust_best_fit[kk],
    fmu,
    xi_PAH,
    xi_MIR,
    xi_W,
    xi_C,
    simonMagphys$T_W_BC_best_fit[kk],
    simonMagphys$T_C_ISM_best_fit[kk]
  )
  names(ptest) = names(p)
  dustTest = L_dust(waveout, p = ptest)
  
  list(log10(simonMagphys$mass_dust_best_fit[kk]), log10(dustTest$Mdust), log10(simonMagphys$L_dust_best_fit[kk]), log10(trapz(dustTest$wave, dustTest$Ltot)))
}
## Confirmed that the magphys has been pretty much reproduced
maghist(unlist(simonCompare[,1]) - unlist(simonCompare[,2]))
maghist(unlist(simonCompare[,3]) - unlist(simonCompare[,4]))

save.image(
  file = "~/Documents/DustMassDensity/data/magphys_vs_dale.Rdata"
)


# magphys_Mass_SED = approxfun(
#   totSED$DustEmit$wave,
#   (magphys_standard$BC$WBC+magphys_standard$ISM$WISM+magphys_standard$ISM$CISM),
#   rule = 2
# )
# Dale_M2L_new_ = sapply(1:64, function(i){
#   
#   Msol = 1.989e30
#   mH = 1.674e-27
#   Lsol = 3.828e26
#   DTH = 0.0073
#   
#   Dale_ff = Dale_Orig$Aspec[[1]][i, ] / Dale_Orig$Wave
#   Dale_ff_norm = Dale_ff / sum(
#     c(0, diff(Dale_Orig$Wave)) * Dale_ff
#   )
#   
#   # wave_idx = Dale_Orig$Wave > Dale_Orig$Wave[which.max(Dale_Orig$Aspec[[1]][i,])] & Dale_Orig$Wave < 1e7
#   wave_idx = Dale_Orig$Wave > 2e5
#   
#   LL = function(parm){
#     
#     Temp = parm[1]
#     # Ldust = parm[2]
# 
#     rest_frame = greybody_norm(
#       wave = Dale_Orig$Wave, 
#       Temp = Temp, 
#       beta = 1.5,
#       z = 0,
#       norm = 1
#     )
#     
#     likelihood = sum( (((rest_frame[wave_idx]) - (Dale_ff_norm[wave_idx]))/Dale_ff_norm[wave_idx])^2 )
#     return(likelihood)
#   }
#   
#   opt = optim(
#     par = c(50),
#     fn = LL,
#     method = "L-BFGS-B",
#     lower = c(5),
#     upper = c(150),
#     hessian = TRUE
#   )
#   
#   rest_frame = greybody_norm(
#     wave = Dale_Orig$Wave, 
#     Temp = opt$par[1], 
#     beta = 1.5,
#     z = 0,
#     norm = 1
#   )
#   # magphys_mass = magphys_Mass_SED(Dale_Orig$Wave)
#   # magphys_mass_norm = magphys_mass / sum(c(0, diff(Dale_Orig$Wave)) * magphys_mass)
#   
#   qPAH_VSG = 0.14
#   # weight_magphys = pmin(pmax(magphys_Mass_contribution(Dale_Orig$Wave), qPAH_VSG), 1-qPAH_VSG)
#   # weight_magphys[Dale_Orig$Wave > 1e7] = 1-qPAH_VSG
#   
#   png(paste0("~/Documents/DustMassDensity/plots/dale_M2L_new/Dale_", i, ".png"), width = 10, height = 6, units = "in", res = 240)
#   par(mfrow = c(1,2))
#   magplot(
#     Dale_Orig$Wave,
#     Dale_ff_norm, 
#     log = "xy", type = "l", lwd = 2, 
#     xlim = c(1e4, 1e7),
#     ylim = c(1e-10, 1e-5),
#     xlab = "Wave [Ang]",
#     ylab = "Flux"
#   )
#   lines(
#     Dale_Orig$Wave,
#     rest_frame, 
#     col = "red"
#   )
#   # lines(
#   #   Dale_Orig$Wave,
#   #   magphys_mass_norm, 
#   #   col = "blue"
#   # )
#   legend(
#     x = "bottomleft", 
#     legend = paste0("alpha = ", Dale_Orig$alpha_SF[i])
#   )
#   weight = pmin(pmax(rest_frame/Dale_ff_norm, qPAH_VSG), 1-qPAH_VSG)
#   magplot(
#     Dale_Orig$Wave, weight, 
#     xlim = c(1e4, 1e7),
#     ylim = c(0,1),
#     log = "x", type = "l", lwd = 3,
#     xlab = "Wave [Ang]",
#     ylab = "DTH weight"
#   )
#   weight[Dale_Orig$Wave > Dale_Orig$Wave[which.max(weight)]] = max(weight)
#   lines(
#     Dale_Orig$Wave, weight, col = "red", lwd = 1
#   )
#   dev.off()
# 
#   foo = ( (Dale_Orig$Aspec[[1]][i, ] / Lsol) / (DTH*weight * mH/Msol))/Dale_Orig$Wave
#   sum(
#     c(0, diff(Dale_Orig$Wave)) * foo
#   )
# })
