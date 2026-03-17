library(data.table)
library(magicaxis)
library(rhdf5)
library(pracma)
library(ProSpect)

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

data_dir = "/Volumes/JordanData/GAMA-DEVILS-SFR-AGN/"
devilsd10_AGN = fread(paste0(data_dir, "/data/AGNTotalCat_MasterCat4.csv"))
devilsd10_noAGN = readRDS(paste0(data_dir, '/data/DEVILS_D10ProSpectCat_02_02_2021_v0.3.rds')) #catalogue that Jess done without AGN contribution
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
  "StellarMass", "dustmass.birth","dustmass.screen","dustmass.total","dustlum.birth","dustlum.screen","dustlum.total","Zfinal", "SFRburst", "tau_screen", "tau_birth", "alpha_SF_birth", "alpha_SF_screen",
  "StellarMass_LB", "dustmass.birth_LB","dustmass.screen_LB","dustmass.total_LB","dustlum.birth_LB","dustlum.screen_LB","dustlum.total_LB","Zfinal_LB", "SFRburst_LB", "tau_screen_LB", "tau_birth_LB", "alpha_SF_birth_LB", "alpha_SF_screen_LB",
  "StellarMass_UB", "dustmass.birth_UB","dustmass.screen_UB","dustmass.total_UB","dustlum.birth_UB","dustlum.screen_UB","dustlum.total_UB","Zfinal_UB", "SFRburst_UB", "tau_screen_UB", "tau_birth_UB", "alpha_SF_birth_UB", "alpha_SF_screen_UB",
  "FIRInput"
)
devilsd10_hybrid = devilsd10_noAGN[, devilsd10_col_names]
devilsd10_hybrid$FITTYPE = "STELLAR"
devilsd10_hybrid$AGN = 0
devilsd10_hybrid$AGNfrac = 0
devilsd10_AGN_idx = devilsd10_AGN$AGNfrac >= 0.1 & devilsd10_AGN$LP > devilsd10_noAGN$LP
message("AGN preferred fraction: ", sum(devilsd10_AGN_idx)/dim(devilsd10_hybrid)[1])
devilsd10_hybrid[devilsd10_AGN_idx, devilsd10_col_names] = devilsd10_AGN[devilsd10_AGN_idx, devilsd10_col_names]
devilsd10_hybrid$FITTYPE[devilsd10_AGN_idx] = "STELLAR+AGN"
devilsd10_hybrid$AGN[devilsd10_AGN_idx] = devilsd10_AGN$AGNlum[devilsd10_AGN_idx]
devilsd10_hybrid$AGNfrac[devilsd10_AGN_idx] = devilsd10_AGN$AGNfrac[devilsd10_AGN_idx]

devilsd10_hybrid$new_dust = devilsd10_hybrid$dustlum.birth/Dale_M2L_variableDTH_func(devilsd10_hybrid$alpha_SF_birth) + devilsd10_hybrid$dustlum.screen/Dale_M2L_variableDTH_func(devilsd10_hybrid$alpha_SF_screen)

devilsd10_hybrid$Mgas = devilsd10_hybrid$new_dust * RR14_BPL(Z = 10^devilsd10_hybrid$Zfinal, doDTG = FALSE)/0.0073/RR14_BPL(Z = 10^devilsd10_hybrid$Zfinal, doDTG = TRUE)
devilsd10_hybrid$N100_Y = devilsd10_AGN$N100_Y
devilsd10_hybrid$R50_Y = devilsd10_AGN$R50_Y

devilsd10_highGas = devilsd10_hybrid[devilsd10_hybrid$Mgas >= 10^11.5, ]
dim(devilsd10_highGas)

median(devilsd10_highGas$SFRburst)
log10(median(devilsd10_highGas$Mgas))

maghist(
  devilsd10_highGas$z
)
maghist(
  devilsd10_hybrid$SFRburst[devilsd10_hybrid$SFRburst>1e-4], 
  log = "xy"
)
maghist(
  devilsd10_highGas$SFRburst, 
  log = "xy", col = "red", add = TRUE
)
maghist(
  devilsd10_hybrid$tau_screen, 
  log = "y", freq = FALSE
)
maghist(
  devilsd10_highGas$tau_screen, 
  log = "y", 
  add = TRUE, 
  col = "red", freq = FALSE
)
maghist(
  devilsd10_hybrid$tau_birth, 
  log = "y"
)
maghist(
  devilsd10_highGas$tau_birth, 
  log = "y", 
  add = TRUE, 
  col = "red"
)

maghist(
  devilsd10_highGas$alpha_SF_screen_UB[devilsd10_highGas$FIRInput == 1L] - devilsd10_highGas$alpha_SF_screen_LB[devilsd10_highGas$FIRInput == 1L]
)
maghist(
  devilsd10_hybrid$alpha_SF_screen_UB[devilsd10_hybrid$FIRInput == 1L] - devilsd10_hybrid$alpha_SF_screen_LB[devilsd10_hybrid$FIRInput == 1L]
)

maghist(
  devilsd10_hybrid$alpha_SF_birth[devilsd10_hybrid$FIRInput == 1L], 
  log = "y"
)
maghist(
  devilsd10_highGas$alpha_SF_birth[devilsd10_highGas$FIRInput == 1L], 
  log = "y", 
  add = TRUE, 
  col = "red"
)
maghist(
  devilsd10_hybrid$alpha_SF_screen[devilsd10_highGas$FIRInput == 1L], 
  log = "y"
)
maghist(
  devilsd10_highGas$alpha_SF_screen[devilsd10_highGas$FIRInput == 1L], 
  log = "y", 
  add = TRUE, 
  col = "red"
)

magplot(
  devilsd10_highGas$dustlum.total, 
  devilsd10_highGas$dustmass.total,
  log = "xy", 
  xlim = c(1e1, 1e14),
  ylim = c(1e1, 1e14)
)
points(
  devilsd10_hybrid$dustlum.total, 
  devilsd10_hybrid$dustmass.total,
  col = "red", 
  pch = "."
)

magplot(
  devilsd10_highGas$dustlum.total, 
  devilsd10_highGas$tau_screen,
  log = "x", 
  xlim = c(1e1, 1e14),
  ylim = c(-1, 1.5)
)
points(
  devilsd10_hybrid$dustlum.total, 
  devilsd10_hybrid$tau_screen,
  col = "red", 
  pch = "."
)


maghist(
  devilsd10_highGas$dustmass.total/devilsd10_highGas$Mgas
)
maghist(
  devilsd10_hybrid$dustmass.total/devilsd10_hybrid$Mgas
)

pie(c(sum(devilsd10_highGas$FIRInput), length(devilsd10_highGas$FIRInput)-sum(devilsd10_highGas$FIRInput)),
    labels = c("FIRInput", "No FIRInput"))

maghist(
  devilsd10_highGas$StellarMass[devilsd10_highGas$FIRInput == 1L], 
  log = "xy", xlim = c(1e9, 1e13)
)
maghist(
  devilsd10_highGas$StellarMass[devilsd10_highGas$FIRInput == 0L], 
  log = "xy", xlim = c(1e9, 1e13)
)

magplot(
  devilsd10_highGas$Zfinal,
  devilsd10_highGas$Mgas/devilsd10_highGas$dustmass.total
)
magplot(
  devilsd10_highGas$Zfinal,
  1/RR14_BPL(10^devilsd10_highGas$Zfinal), 
  log = "y"
)

magplot(
  devilsd10_highGas$dustmass.total, 
  devilsd10_highGas$Mgas,
  log = "xy"
)
magplot(
  devilsd10_highGas$StellarMass, devilsd10_highGas$Mgas, 
  z = devilsd10_highGas$AGNfrac,
  log = "xy"
)
magplot(
  devilsd10_highGas$dustmass.total, 
  devilsd10_highGas$SFRburst,
  log = "xy"
)
magplot(
  devilsd10_highGas$z, 
  devilsd10_highGas$SFRburst,
  z = devilsd10_highGas$FIRInput,
  log = "y"
)


maghist(
  devilsd10_highGas$z[devilsd10_highGas$StellarMass<=1e8]
)

devilsd10_highGas[devilsd10_highGas$StellarMass<=1e8,]

magplot(
  devilsd10_highGas$StellarMass, devilsd10_highGas$dustmass.total, 
  log = "xy"
)

h5file = '~/Documents/DustMassDensity/data/all_data.h5'

rho_gas = h5read(h5file, "cosmic/MgasHybridCorr")

for(i in c(13,14,15,16)){
  dfGas = h5read(h5file, paste0("binDF/MgasHybridCorr/zb",i))
  fitGas = h5read(h5file, paste0("fitDF/MgasHybridCorr/zb",i))
  
  spl = smooth.spline(
    x = c(dfGas$x[dfGas$phi > 0]),
    y = c(log10(dfGas$phi[dfGas$phi > 0])), 
    w = c(dfGas$err[dfGas$phi > 0]/(log(10)*dfGas$phi[dfGas$phi > 0]))^-2, 
    df = 5
  )
  spl_out = predict(
    object = spl,
    x = log10(fitGas$x)
  )
  
  # magplot(
  #   10^dfGas$x, 
  #   dfGas$phi, 
  #   log = "xy",
  #   pch = 16, 
  #   xlim = c(1e4, 1e14), 
  #   ylim = c(1e-8, 1e-1)
  # )
  # magerr(
  #   10^dfGas$x, 
  #   dfGas$phi, 
  #   ylo = dfGas$err
  # )
  # lines(
  #   fitGas$x,
  #   fitGas$Q50
  # )
  # lines(
  #   10^spl_out$x,
  #   10^spl_out$y, 
  #   col = "blue"
  # )
  
  magplot(
    fitGas$x,
    fitGas$x * fitGas$Q50,
    type = "l",
    log = "xy",
    xlim = c(1e4, 1e14),
    ylim = c(1e5, 1e8)
  )
  lines(
    10^spl_out$x,
    10^spl_out$y * 10^spl_out$x,
    col = "blue"
  )
  
  spline_csmh = log10(trapz(x = spl_out$x[spl_out$x>=min(dfGas$x[dfGas$phi > 0])], y = 10^spl_out$x[spl_out$x>=min(dfGas$x[dfGas$phi > 0])] * 10^spl_out$y[spl_out$x>=min(dfGas$x[dfGas$phi > 0])]))
  
  ds_csmh = log10(trapz(x = log10(fitGas$x[log10(fitGas$x)>=min(dfGas$x[dfGas$phi > 0])]), y = fitGas$x[log10(fitGas$x)>=min(dfGas$x[dfGas$phi > 0])] * fitGas$Q50[log10(fitGas$x)>=min(dfGas$x[dfGas$phi > 0])]))
 
  print(spline_csmh)
  print(ds_csmh)
  print("\n")
}


magplot(
  zhang21$MHI[zhang21$z == 0.85],
  zhang21$phi[zhang21$z == 0.85],
  ylim = c(-5.5, 0.5),
  xlim = c(8, 12)
)
magerr(
  zhang21$MHI[zhang21$z == 0.85],
  zhang21$phi[zhang21$z == 0.85],
  ylo = zhang21$philo[zhang21$z == 0.85],
  yhi = zhang21$phihi[zhang21$z == 0.85],
)
points(
  dfGas$x, 
  log10(dfGas$phi), 
  pch = 16
)

magplot(
  zhang21$MHI[zhang21$z == 1.20],
  zhang21$phi[zhang21$z == 1.20],
  ylim = c(-5.5, 0.5),
  xlim = c(8, 12)
)
magerr(
  zhang21$MHI[zhang21$z == 1.20],
  zhang21$phi[zhang21$z == 1.20],
  ylo = zhang21$philo[zhang21$z == 1.20],
  yhi = zhang21$phihi[zhang21$z == 1.20],
)
points(
  dfGas$x, 
  log10(dfGas$phi), 
  pch = 16
)


magplot(
  seq(1,4,0.1),
  Dale_M2L_variableDTH_func(seq(1,4,0.1))/Dale_M2L_func(seq(1,4,0.1))
)
magplot(
  devilsd10_hybrid$z,
  Dale_M2L_variableDTH_func(devilsd10_hybrid$alpha_SF_birth)/Dale_M2L_func(devilsd10_hybrid$alpha_SF_birth), 
  pch = "."
)
magplot(
  devilsd10_hybrid$z,
  Dale_M2L_variableDTH_func(devilsd10_hybrid$alpha_SF_screen)/Dale_M2L_func(devilsd10_hybrid$alpha_SF_screen), 
  pch = ".", 
  col = alpha("black", 0.1)
)
