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
# plot(totSED)

## PAH emission template from Madden+06 of a photodissociation region in the MW
m17_pdr_pah = read.table("~/Documents/DustMassDensity/data/M17_PDR_PAH.csv", header = TRUE, sep = ",")
names(m17_pdr_pah) = c("mu", "flux")
m17_approx = cubicspline(x = c(m17_pdr_pah$mu*1e4), y = c(m17_pdr_pah$flux), xi = waveout, endp2nd = TRUE, der = c(1.721362e-06, -6.885449e-06))
m17_approx[m17_approx < 0] = 0
# m17_approx = approx(x = c(m17_pdr_pah$mu*1e4), y = c(m17_pdr_pah$flux), xout = waveout, yleft = 0, yright = 0)$y
Florentzian = (1 + ((waveout/1e8)^-1 - 3039.1)^2/(19.4)^2)^-1 * approx(x = waveout/1e4, y = m17_approx, xout = 11.3)$y * 0.1

magplot(
  waveout/1e4,
  m17_approx,
  log = "xy",
  type = "l",
  xlim = c(1, 20),
  ylim = c(1e-6, 1)
)
points(
  m17_pdr_pah$mu, 
  m17_pdr_pah$flux,
  pch = 16, col = "red"
)
lines(
  waveout/1e4, Florentzian, col = "blue"
)
lines(
  waveout/1e4, Florentzian + m17_approx, col = "green"
)

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

magplot(
  waveout/1e4, 
  PAH_emission_lines * waveout, 
  log = "xy",
  type = "l",
  xlim = c(1, 20),
  col = "blue"
  # ylim= c(1e-6, 1)
)
lines(
  waveout/1e4, grey_continuum_scaled*waveout, col = "red"
)
lines(
  waveout/1e4, waveout*m17_norm(waveout), col = "black"
)

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

Niters = 1000
alpha_birth_samples = runif(Niters, 0, 4)
alpha_screen_samples = runif(Niters, 0, 4)
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
  )
  
  df = list(
    "wave" = testSED$DustEmit$wave,
    "lum" = testSED$DustEmit$lum,
    "dustlum" = testSED$dustlum,
    "dustmass" = testSED$dustmass
  )
  return(df)
}

magphys_parm_mat = matrix(data = 0L, nrow = length(p)-1, ncol = Niters)
for(i in 1:Niters){
  
  frac_BC = runif(1, 0, 1)
  xi_W_BC_sample = runif(1, 0, 1.0)
  xi_MIR_BC_sample = (1-xi_W_BC_sample) * frac_BC
  xi_PAH_BC_sample = (1-xi_W_BC_sample) * (1-frac_BC)
  magphys_parm_mat[,i] = c(
    runif(1, 0, 1),
    xi_PAH_BC_sample,
    xi_MIR_BC_sample,
    xi_W_BC_sample,
    runif(1, 0, 1),
    runif(1, 10, 100),
    runif(1, 10, 100)
  )
  
  # frac_BC = 10^runif(1, -4, 0)
  # xi_W_BC_sample = 10^runif(1, -4, 0)
  # xi_MIR_BC_sample = (1-xi_W_BC_sample) * frac_BC
  # xi_PAH_BC_sample = (1-xi_W_BC_sample) * (1-frac_BC)
  # magphys_parm_mat[,i] = c(
  #   10^runif(1, -4, 0),
  #   xi_PAH_BC_sample,
  #   xi_MIR_BC_sample,
  #   xi_W_BC_sample,
  #   10^runif(1, -4, 0),
  #   10^runif(1, 1, 2),
  #   10^runif(1, 1, 2)
  # )
}
magphys_dust_lum = foreach(i = 1:Niters) %do% {
  if(i %% 100 == 0){message(i)}
  magphys_parm = c(unname(totSED$dustlum["total"]), magphys_parm_mat[,i])
  names(magphys_parm) = names(p)
  magphys_lum = L_dust(totSED$DustEmit$wave, p = magphys_parm)
}

par(mar = c(2.5, 2.5, 1.5, 1.5), oma = rep(1.5, 4), family = "DejaVu Sans", cex = 2.0)
magplot(
  NA,
  type = "l", 
  log = "xy",
  lwd = 3,
  xlim = c(1e3, 1e7),
  ylim = c(1e7, 3e10),
  xlab = "Wavelength [Ang]",
  ylab = "Lum [𝐿⊙ / Ang]",
  grid = FALSE, 
  cex.lab = 2.0
)
for(i in 1:length(prospect_dale_samples)){
  lines(
    prospect_dale_samples[[i]]$wave,
    prospect_dale_samples[[i]]$lum * prospect_dale_samples[[i]]$wave,
    col = alpha("darkorange", 0.1)
  )
}
for(i in 1:Niters){
  lines(
    magphys_dust_lum[[i]]$wave, 
    magphys_dust_lum[[i]]$wave * magphys_dust_lum[[i]]$Ltot, 
    col = alpha("darkred", 0.1),
    lwd = 1, lty = 1.0
  )
}
legend(
  x = "topleft",
  legend = c("Dale samples", "Magphys samples"),
  lwd = 2,
  lty = c(1,1),
  col = c("darkorange", "darkred")
)

magphys_standard =  L_dust(totSED$DustEmit$wave, p = p)
magplot(
  totSED$FinalLum$wave, totSED$FinalLum$lum * totSED$FinalLum$wave, 
  type = "l", 
  log = "xy",
  lwd = 3,
  xlim = c(1e2, 1e7),
  ylim = c(1e7, 3e10),
  xlab = "Wavelength [Ang]",
  ylab = "Lum [Lsun / Ang]"
)
lines(
  totSED$DustEmit$wave, totSED$DustEmit$lum * totSED$DustEmit$wave,
  col = alpha("darkorange", 0.8),
  lwd = 5
)
lines(
  totSED$DustEmit$wave, magphys_standard$Ltot * totSED$DustEmit$wave,
  col = alpha("darkred", 0.8),
  lwd = 5
)
lines(
  totSED$DustEmit$wave, magphys_standard$BC$Ltot * totSED$DustEmit$wave,
  col = alpha("darkred", 0.8),
  lwd = 3, lty = 2
)
lines(
  totSED$DustEmit$wave, magphys_standard$ISM$Ltot * totSED$DustEmit$wave,
  col = alpha("darkred", 0.8),
  lwd = 3, lty = 3
)
lines(
  totSED$DustEmit$wave, magphys_standard$ISM$PAH * totSED$DustEmit$wave,
  col = alpha("red", 0.8),
  lwd = 3, lty = 3
)
lines(
  totSED$DustEmit$wave, magphys_standard$ISM$WISM * totSED$DustEmit$wave,
  col = alpha("red", 0.8),
  lwd = 3, lty = 3
)
lines(
  totSED$DustEmit$wave, magphys_standard$BC$WBC * totSED$DustEmit$wave,
  col = alpha("red", 0.8),
  lwd = 3, lty = 3
)
legend(
  x = "topleft",
  legend = c("Final lum", "Dale from ProSpect", "Reproduce magphys", "Magphys BC", "Magphys ISM"),
  lwd = 2,
  lty = c(1,1,1,2,3),
  col = c("black", "darkorange", "darkred", "darkred", "darkred")
)

dale_spec_samples = as.matrix(foreach(i = 1:Niters, .combine = rbind) %do% {
  prospect_dale_samples[[i]]$lum 
})
magphys_spec_samples = as.matrix(foreach(i = 1:Niters, .combine = rbind) %do% {
  magphys_dust_lum[[i]]$Ltot 
})

dale_dust_lum_samples = sapply(prospect_dale_samples, function(x)unname(x$dustlum["total"]))
dale_dust_lum_birth_samples = sapply(prospect_dale_samples, function(x)unname(x$dustlum["birth"]))
dale_dust_lum_screen_samples = sapply(prospect_dale_samples, function(x)unname(x$dustlum["screen"]))
dale_dust_mass_samples = sapply(prospect_dale_samples, function(x)unname(x$dustmass['total']))
dale_dust_mass_birth_samples = sapply(prospect_dale_samples, function(x)unname(x$dustmass['birth']))
dale_dust_mass_screen_samples = sapply(prospect_dale_samples, function(x)unname(x$dustmass['screen']))
dale_FIR_samples = sapply(prospect_dale_samples, function(x)trapz(x$wave[x$wave>2e5], x$lum[x$wave>2e5]))

magphys_dust_lum_samples = sapply(magphys_dust_lum, function(x)trapz(x$wave, x$Ltot))
magphys_dust_lum_no_PAH_samples = sapply(magphys_dust_lum, function(x){trapz(x$wave, x$BC$WBC) + trapz(x$wave, x$ISM$WISM) + trapz(x$wave, x$ISM$CISM)})
magphys_dust_mass_samples = sapply(magphys_dust_lum,  function(x)unname(x$Mdust))
magphys_FIR_samples = sapply(magphys_dust_lum, function(x)trapz(x$wave[x$wave>2e5], x$Ltot[x$wave>2e5]))

df_samples = data.frame(
  "DaleLum" = dale_dust_lum_samples,
  "DaleLumBC" = dale_dust_lum_birth_samples,
  "DaleLumISM" = dale_dust_lum_screen_samples,
  "DaleMdust" = dale_dust_mass_samples,
  "DaleMdustBC" = dale_dust_mass_birth_samples,
  "DaleMdustISM" = dale_dust_mass_screen_samples,
  "DaleFIR" = dale_FIR_samples,
  "MagphysLum" = magphys_dust_lum_samples,
  "MagphysLumNoPAH" = magphys_dust_lum_no_PAH_samples,
  "MagphysMdust" = magphys_dust_mass_samples,
  "MagphysFIR" = magphys_FIR_samples
)

Md_corr = (median(dale_FIR_samples / magphys_FIR_samples) * median(magphys_dust_lum_samples / magphys_dust_lum_no_PAH_samples))
message("Average factor of difference? = ", Md_corr)

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
  "Dale" = totSED$DustEmit$lum,
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
h5write(obj = Md_corr, file = h5file, name = 'Md_corr')
h5write(
  obj = df_samples, file = h5file, name = "DaleMagphys/Samples"
)
h5write(
  obj = spec_samples, file = h5file, name = "DaleMagphys/Spec"
)
h5write(
  obj = spec_standard, file = h5file, name = "DaleMagphys/StandardSpec"
)


maghist(
  log10(dale_dust_mass_samples), seq(0.5, 8.5, 0.25), col = alpha("darkorange", 0.5), log = "y", ylim = c(1,500), 
  ylab = "Frequency", xlab = "log10(Dust mass/Msun)"
)
maghist(
  log10(dale_dust_mass_birth_samples), seq(3.5, 8.5, 0.25), col = alpha("green", 0.5), log = "y", ylim = c(1,150), add = TRUE,
  ylab = "Frequency", xlab = "log10(Dust mass/Msun)"
)
maghist(
  log10(dale_dust_mass_screen_samples), seq(3.5, 8.5, 0.25), col = alpha("blue", 0.5), log = "y", ylim = c(1,150), add = TRUE,
  ylab = "Frequency", xlab = "log10(Dust mass/Msun)"
)
maghist(
  log10(magphys_dust_mass_samples), seq(3.5, 8.5, 0.25), add = TRUE, col = alpha("darkred", 0.5), log = "y", ylim = c(1,150)
)
legend(
  x = "topleft",
  legend = c("Dale", "Magphys"),
  lwd = 5,
  lty = c(1,1),
  col = c(alpha("darkorange", 0.5), alpha("darkred", 0.5)),
  cex = 1.5
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
