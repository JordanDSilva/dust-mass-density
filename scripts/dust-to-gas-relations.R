library(magicaxis)
library(Rfits)
library(ProSpect)
## Test different DTG

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
## Use their default PG16S metallicity calibrator 
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

draine07 = 0.0073
RR14_BPL(0.02)
deVis19_PL(0.02)

magplot(
  10^seq(-4, log10(0.05), 0.01), 
  RR14_BPL(10^seq(-4, log10(0.05), 0.01)), 
  log = "xy"
)
lines(
  10^seq(-4, log10(0.05), 0.1), 
  deVis19_PL(10^seq(-4, log10(0.05), 0.1))
)
abline(h = draine07)

magplot(
  10^seq(-4, log10(0.05), 0.01),
  deVis19_PL(10^seq(-4, log10(0.05), 0.01))/RR14_BPL(10^seq(-4, log10(0.05), 0.01)), 
  log = "x", 
)
abline(v = 0.02)
abline(v = 0.014)


Mdust = 10^9
maghist(
  Mdust/RR14_BPL(10^seq(-4, log10(0.05), 0.01)), log = "xy"
)
maghist(
  Mdust/deVis19_PL(10^seq(-4, log10(0.05), 0.01)), log = "xy"
)


prospect_cat = Rfits_read_table("~/Documents/DustMassDensity/data/GAMA/ProSpectv03.fits")
prospect_cat$D