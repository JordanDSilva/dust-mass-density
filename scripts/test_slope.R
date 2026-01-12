library(ProSpect)
library(magicaxis)
library(rhdf5)
library(Highlander)
library(celestial)
library(data.table)
library(foreach)
library(matrixStats)

tvec = seq(0, 13, 0.1)

h5file = '~/Documents/DustMassDensity/data/all_data.h5'
zmids = h5read(file = h5file, name = "zmids")
lbtmids = h5read(file = h5file, name = "lbtmids")
mdust = h5read(file = h5file, name = "cosmic/MdustHybridCorr")
points(
  lbtmids, mdust$CORR, pch = 16
)

dunne2011 = data.frame(
  fread(
    "~/Documents/DustMassDensity/data/literature_evo/cdmh/dunne11.csv"
  )
)
driver2018 = data.frame(
  fread(
    "~/Documents/DustMassDensity/data/literature_evo/cdmh/driver18.csv"
  )
)
pozzi2020 = data.frame(
  fread(
    "~/Documents/DustMassDensity/data/literature_evo/cdmh/pozzi20.csv"
  )
)
beeston2024 = data.frame(
  fread(
    "~/Documents/DustMassDensity/data/literature_evo/cdmh/beeston24.csv"
  )
)
eales2024 = data.frame(
  fread(
    "~/Documents/DustMassDensity/data/literature_evo/cdmh/eales24.csv"
  )
)
berta2025 = data.frame(
  fread(
    "~/Documents/DustMassDensity/data/literature_evo/cdmh/berta25.csv"
  )
)
chiang2025 = data.frame(
  fread(
    "~/Documents/DustMassDensity/data/literature_evo/cdmh/chiang25.csv"
  )
)

# lit_combined = data.frame(
#   "z" = c(dunne2011$z, pozzi2020$z, beeston2024$z, eales2024$z, chiang2025$z, driver2018$z, zmids),
#   "lbt" = cosdistTravelTime(c(dunne2011$z, pozzi2020$z, beeston2024$z, eales2024$z, chiang2025$z, driver2018$z, zmids), ref = "Planck18"),
#   "cdmh" = c(dunne2011$cdmh, pozzi2020$cdmh, beeston2024$cdmh, eales2024$cdmh, chiang2025$cdmh, driver2018$cdmh, mdust$CORR),
#   "err" = c(dunne2011$err, pozzi2020$err, beeston2024$err, eales2024$err, 0.5*(chiang2025$errhi+chiang2025$errlo), driver2018$err, mdust$CORR_ERR)
# )
lit_combined = data.frame(
  "z" = c(dunne2011$z, pozzi2020$z, beeston2024$z, eales2024$z, berta2025$z, chiang2025$z),
  "lbt" = cosdistTravelTime(c(dunne2011$z, pozzi2020$z, beeston2024$z, eales2024$z, berta2025$z, chiang2025$z), ref = "Planck18"),
  "cdmh" = c(dunne2011$cdmh, pozzi2020$cdmh, beeston2024$cdmh, eales2024$cdmh, berta2025$cdmh, chiang2025$cdmh),
  "err" = c(dunne2011$err, pozzi2020$err, beeston2024$err, eales2024$err, 0.5*(berta2025$errhi+berta2025$errlo), 0.5*(chiang2025$errhi+chiang2025$errlo))
)

LL = function(p, Data){
  mSFR = 10^p[1]
  mpeak = p[2]
  mperiod = p[3]
  mskew = p[4]
  
  model = massfunc_snorm_trunc(
    age = Data$x*1e9,
    mSFR = mSFR,
    mpeak = mpeak,
    mperiod = mperiod,
    mskew = mskew, 
    
    magemax = 13.8,
    mtrunc = 2.0
  )
  
  loglikeihood = dnorm(
    x = log10(Data$y),
    mean = log10(model),
    sd = Data$err/(log(10)*Data$y),
    log = TRUE
  )
  logposterior = sum(loglikeihood, na.rm = TRUE) + Data$prior(p)
  if(is.infinite(logposterior)){
    logposterior = -999999
  }
  return(logposterior)
}

pfunc = function(p){
  sum(
    # dnorm(p[2], 10, 0.5, log = TRUE),
    0
  )
}

Niters = 5000
fit_me = Highlander(
  parm = c(6.5, 10, 1.0, 0.5),
  Data = list(
    "x" = lbtmids,
    "y" = mdust$CORR,
    "err" = mdust$CORR_ERR,
    prior = pfunc
  ),
  likefunc = LL,
  liketype = "max",
  seed = 42069,
  lower = c(4.5, 0, 0.1, -5.0),
  upper = c(8.5, 13.8, 10.0, 5.0),
  Niters = c(Niters,Niters), 
  optim_iters = 2
)
fit_lit_combined = Highlander(
  parm = c(6.5, 10, 1.0, 0.5),
  Data = list(
    # "x" = c(dunne2011$z, pozzi2020$z, beeston2024$z, eales2024$z),
    # "y" = c(dunne2011$cdmh, pozzi2020$cdmh, beeston2024$cdmh, eales2024$cdmh),
    # "err" = c(dunne2011$err, pozzi2020$err, beeston2024$err, eales2024$err),
    "x" = lit_combined$lbt,
    "y" = lit_combined$cdmh,
    "err" = lit_combined$err,
    prior = pfunc
  ),
  likefunc = LL,
  liketype = "max",
  seed = 42069,
  lower = c(4.5, 0, 0.1, -5.0),
  upper = c(8.5, 13.8, 10.0, 5.0),
  Niters = c(Niters,Niters), 
  optim_iters = 2
)
fit_dunne2011 = Highlander(
  parm = c(6.5, 10, 1.0, 0.5),
  Data = list(
    "x" = cosdistTravelTime(dunne2011$z, ref = "Planck18"),
    "y" = dunne2011$cdmh,
    "err" = dunne2011$err,
    prior = pfunc
  ),
  likefunc = LL,
  liketype = "max",
  seed = 42069,
  lower = c(4.5, 0, 0.1, -5.0),
  upper = c(8.5, 13.8, 10.0, 5.0),
  Niters = c(Niters,Niters), 
  optim_iters = 2
)
fit_driver2018 = Highlander(
  parm = c(6.5, 10, 1.0, 0.5),
  Data = list(
    "x" = driver2018$lbt,
    "y" = driver2018$cdmh,
    "err" = driver2018$err,
    prior = pfunc
  ),
  likefunc = LL,
  liketype = "max",
  seed = 42069,
  lower = c(4.5, 0, 0.1, -5.0),
  upper = c(8.5, 13.8, 10.0, 5.0),
  Niters = c(Niters,Niters), 
  optim_iters = 2
)
fit_pozzi2020 = Highlander(
  parm = c(6.5, 10, 1.0, 0.5),
  Data = list(
    "x" = cosdistTravelTime(pozzi2020$z, ref = "Planck18"),
    "y" = pozzi2020$cdmh,
    "err" = pozzi2020$err,
    prior = pfunc
  ),
  likefunc = LL,
  liketype = "max",
  seed = 42069,
  lower = c(4.5, 0, 0.1, -5.0),
  upper = c(8.5, 13.8, 10.0, 5.0),
  Niters = c(Niters,Niters), 
  optim_iters = 2
)
fit_beeston2024 = Highlander(
  parm = c(6.5, 10, 1.0, 0.5),
  Data = list(
    "x" = cosdistTravelTime(beeston2024$z, ref = "Planck18"),
    "y" = beeston2024$cdmh,
    "err" = beeston2024$err,
    prior = pfunc
  ),
  likefunc = LL,
  liketype = "max",
  seed = 42069,
  lower = c(4.5, 0, 0.1, -5.0),
  upper = c(8.5, 13.8, 10.0, 5.0),
  Niters = c(Niters,Niters), 
  optim_iters = 2
)
fit_eales2024 = Highlander(
  parm = c(6.5, 10, 1.0, 0.5),
  Data = list(
    "x" = cosdistTravelTime(eales2024$z, ref = "Planck18"),
    "y" = eales2024$cdmh,
    "err" = eales2024$err,
    prior = pfunc
  ),
  likefunc = LL,
  liketype = "max",
  seed = 42069,
  lower = c(4.5, 0, 0.1, -5.0),
  upper = c(8.5, 13.8, 10.0, 5.0),
  Niters = c(Niters,Niters), 
  optim_iters = 2
)
fit_chiang2025 = Highlander(
  parm = c(6.5, 10, 1.0, 0.5),
  Data = list(
    "x" = cosdistTravelTime(chiang2025$z, ref = "Planck18"),
    "y" = chiang2025$cdmh,
    "err" = 0.5*(chiang2025$errhi+chiang2025$errlo),
    prior = pfunc
  ),
  likefunc = LL,
  liketype = "max",
  seed = 42069,
  lower = c(4.5, 0, 0.1, -5.0),
  upper = c(8.5, 13.8, 10.0, 5.0),
  Niters = c(Niters,Niters), 
  optim_iters = 2
)


tvec = seq(0, 13, 0.1)
yy_me = massfunc_snorm_trunc(
  age = tvec*1e9, 
  mSFR = 10^fit_me$parm[1], 
  mpeak = fit_me$parm[2], 
  mperiod = fit_me$parm[3], 
  mskew = fit_me$parm[4], 
  magemax = 13.8,
  mtrunc = 2.0
)
yy_combined = massfunc_snorm_trunc(
  age = tvec*1e9, 
  mSFR = 10^fit_lit_combined$parm[1], 
  mpeak = fit_lit_combined$parm[2], 
  mperiod = fit_lit_combined$parm[3], 
  mskew = fit_lit_combined$parm[4], 
  magemax = 13.8,
  mtrunc = 2.0
)
yy_dunne2011 = massfunc_snorm_trunc(
  age = tvec*1e9, 
  mSFR = 10^fit_dunne2011$parm[1], 
  mpeak = fit_dunne2011$parm[2], 
  mperiod = fit_dunne2011$parm[3], 
  mskew = fit_dunne2011$parm[4], 
  magemax = 13.8,
  mtrunc = 2.0
)
yy_driver2018 = massfunc_snorm_trunc(
  age = tvec*1e9, 
  mSFR = 10^fit_driver2018$parm[1], 
  mpeak = fit_driver2018$parm[2], 
  mperiod = fit_driver2018$parm[3], 
  mskew = fit_driver2018$parm[4], 
  magemax = 13.8,
  mtrunc = 2.0
)
yy_pozzi2020 = massfunc_snorm_trunc(
  age = tvec*1e9, 
  mSFR = 10^fit_pozzi2020$parm[1], 
  mpeak = fit_pozzi2020$parm[2], 
  mperiod = fit_pozzi2020$parm[3], 
  mskew = fit_pozzi2020$parm[4], 
  magemax = 13.8,
  mtrunc = 2.0
)
yy_beeston2024 = massfunc_snorm_trunc(
  age = tvec*1e9, 
  mSFR = 10^fit_beeston2024$parm[1], 
  mpeak = fit_beeston2024$parm[2], 
  mperiod = fit_beeston2024$parm[3], 
  mskew = fit_beeston2024$parm[4], 
  magemax = 13.8,
  mtrunc = 2.0
)
yy_eales2024 = massfunc_snorm_trunc(
  age = tvec*1e9, 
  mSFR = 10^fit_eales2024$parm[1], 
  mpeak = fit_eales2024$parm[2], 
  mperiod = fit_eales2024$parm[3], 
  mskew = fit_eales2024$parm[4], 
  magemax = 13.8,
  mtrunc = 2.0
)
yy_chiang2025 = massfunc_snorm_trunc(
  age = tvec*1e9, 
  mSFR = 10^fit_chiang2025$parm[1], 
  mpeak = fit_chiang2025$parm[2], 
  mperiod = fit_chiang2025$parm[3], 
  mskew = fit_chiang2025$parm[4], 
  magemax = 13.8,
  mtrunc = 2.0
)

compute_snorm_quantile = function(posterior){
  
  yy = foreach(i = 1:dim(posterior)[1], .combine = cbind) %do% {
    massfunc_snorm_trunc(
      age = tvec*1e9, 
      mSFR = 10^posterior[i,1], 
      mpeak = posterior[i,2], 
      mperiod = posterior[i,3], 
      mskew = posterior[i,4], 
      magemax = 13.8,
      mtrunc = 2.0
    )
  }
  calc=rowQuantiles(
    as.matrix(yy),
    probs = c(0.5, 0.16, 0.84)
  )
  colnames(calc) = c("Q50", "Q16", "Q84")
  return(data.frame(calc))
}

quantiles_me = compute_snorm_quantile(fit_me$LD_last$Posterior1)
quantiles_driver2018 = compute_snorm_quantile(fit_driver2018$LD_last$Posterior1)
quantiles_lit_combined = compute_snorm_quantile(fit_lit_combined$LD_last$Posterior1)


magplot(
  tvec, quantiles_me$Q50, log = "y", ylim = c(3e4,4e5), xlim = c(0, 13.8), type = "l", lwd = 5, 
  col = "cornflowerblue", xlab = "Lookback time", ylab = "Dust density"
)
lines(
  tvec, quantiles_me$Q16, lwd = 1
)
lines(
  tvec, quantiles_me$Q84, lwd = 1
)
# points(
#   lbtmids, mdust$CORR, pch = 16, col = "blue"
# )
# magerr(
#   lbtmids, mdust$CORR, ylo = mdust$CORR_ERR, col = "blue"
# )
lines(
  tvec, quantiles_lit_combined$Q50, col = "darkgrey", lwd = 3, lty = 1
)
# lines(
#   tvec, quantiles_lit_combined$Q16, col = "grey", lwd = 1, lty = 2
# )
# lines(
#   tvec, quantiles_lit_combined$Q84, col = "grey", lwd = 1, lty = 2
# )
lines(
  tvec, quantiles_driver2018$Q50, col = "darkgreen", lwd = 3
)
points(
  lit_combined$lbt, lit_combined$cdmh, pch = 1, col = "grey"
)
magerr(
  lit_combined$lbt, lit_combined$cdmh, ylo = lit_combined$err, col = "grey"
)
points(
  lbtmids, mdust$CORR, pch = 16, col = "blue"
)
magerr(
  lbtmids, mdust$CORR, ylo = mdust$CORR_ERR, col = "blue"
)
points(
  driver2018$lbt, driver2018$cdmh, pch = 1, col = "darkgreen"
)
magerr(
  driver2018$lbt, driver2018$cdmh, ylo = driver2018$err, pch = 16, col = "darkgreen"
)
# points(
#   cosdistTravelTime(z = dunne2011$z, ref = "Planck18"), dunne2011$cdmh, pch = 16, col = "red"
# )
# lines(
#   tvec, yy_dunne2011, col = "red"
)
# points(
#   driver2018$lbt, driver2018$cdmh, pch = 16, col = "blue"
# )
lines(
  tvec, yy_driver2018, col = "blue"
)
points(
  cosdistTravelTime(z = pozzi2020$z, ref = "Planck18"), pozzi2020$cdmh, pch = 16, col = "purple"
)
lines(
  tvec, yy_pozzi2020, col = "purple"
)
points(
  cosdistTravelTime(z = beeston2024$z, ref = "Planck18"), beeston2024$cdmh, pch = 16, col = "orange"
)
lines(
  tvec, yy_beeston2024, col = "orange"
)
points(
  cosdistTravelTime(z = eales2024$z, ref = "Planck18"), eales2024$cdmh, pch = 16, col = "darkgreen"
)
lines(
  tvec, yy_eales2024, col = "darkgreen"
)
points(
  cosdistTravelTime(z = chiang2025$z, ref = "Planck18"), chiang2025$cdmh, pch = 16, col = "magenta"
)
lines(
  tvec, yy_chiang2025, col = "magenta"
)

legend(
  x = "bottomleft", 
  col = c("blue", "darkgreen", "grey"),
  legend = c("DSilva+26", "Driver+18", "Various Herschel papers \n incl. Eales+"), 
  lty = 1, 
  lwd = 3
)

magtri(
  fit_me$LD_last$Posterior1
)
magtri(
  fit_lit_combined$LD_last$Posterior1
)
magtri(
  fit_pozzi2020$LD_last$Posterior1
)

fit_me$parm

fit_pozzi2020$parm
fit_eales2024$parm

fit_me$parm
fit_me$parm


magplot(
  tvec, 
  yy_me/yy_driver2018, 
  log = "y", 
  xlim = c(0,10)
)

magplot(
  tvec, 
  yy_me/yy_beeston2024, 
  # log = "y", 
  xlim = c(0,10), ylim = c(0,2)
)

calc_slope = function(x, y, x1, x2){
  
  y1 = approx(
    x = x,
    y = y,
    xout = x1
  )$y
  
  y2 = approx(
    x = x,
    y = y,
    xout = x2
  )$y
  
  (y2-y1)/(x2-x1)
}

calc_slope(
  x = tvec,
  y = (quantiles_me$Q50),
  x1 = 0.01,
  x2 = 5
)
calc_slope(
  x = tvec,
  y = (quantiles_driver2018$Q50),
  x1 = 0.01,
  x2 = 5
)
calc_slope(
  x = tvec,
  y = (quantiles_lit_combined$Q50),
  x1 = 0.01,
  x2 = 5
)
calc_slope(
  x = tvec,
  y = (yy_chiang2025),
  x1 = 0.01,
  x2 = 5
)

