library(celestial)
library(rhdf5)
library(magicaxis)
library(foreach)
library(data.table)
library(hyper.fit)

zvec = seq(0, 30, 0.01)
lbtvec = cosdistTravelTime(z = zvec, ref = "Planck18")
lbt2z = approxfun(
  lbtvec, zvec
)

dsilva25 = data.frame(fread("~/Documents/JWSTCATALOGUE/Archive/DSilva25_CSFH_CAGNH_fit.csv"))

## Log AGN to BHAR
agn2bhar = log10( 3.154e7 * 1e-7 * (0.1 * 3e8**2)^-1 / 2e30)
  
YHe = 0.24 ## Helium fraction
XH = 1 - YHe ## Hydrogen fraction
mProton = 1.67262192e-27 ##kg
Msol2kg = 1.9891e+30
OmegaBaryon = 0.0224 / (68.4/100)^2 ##Planck18
rhoH = OmegaBaryon * cosgrowRhoCrit(z = 0, ref = "Planck18", Dist = "Co") 


zvec = dsilva25$z
Hgrid = cosgrow(zvec, ref = "Planck18")$H
smbhQ50 = 10^(dsilva25$CAGNHQ50 + agn2bhar) / ((Hgrid * (1/3.09e19) / (1/3.15e7)) * (1 + zvec))
smbhQ50 = rev(cumtrapz(x = max(zvec)-rev(zvec), rev(smbhQ50)))
smbhQ16 = 10^(dsilva25$CAGNHQ16 + agn2bhar) / ((Hgrid * (1/3.09e19) / (1/3.15e7)) * (1 + zvec))
smbhQ16 = rev(cumtrapz(x = max(zvec)-rev(zvec), rev(smbhQ16)))
smbhQ84 = 10^(dsilva25$CAGNHQ84 + agn2bhar) / ((Hgrid * (1/3.09e19) / (1/3.15e7)) * (1 + zvec))
smbhQ84 = rev(cumtrapz(x = max(zvec)-rev(zvec), rev(smbhQ84)))

h5file = '~/Documents/DustMassDensity/data/all_data.h5'
zmids = h5read(file = h5file, name = "zmids")
mdust = h5read(file = h5file, name = "cosmic/Mdust")
mdustwAGN = h5read(file = h5file, name = "cosmic/MdustwAGN")
Md_corr = as.numeric(h5read(file = h5file, name = "Md_corr"))
# Md_corr = 1

LSS_corr = h5read(file = h5file, name = "LSSCorrection")

bellstedt20_gasZdensity = data.frame(fread("~/Documents/DustMassDensity/data/literature_evo/metallicity/bellstedt20_gazZdensity.csv"))
bellstedt20_meanZ = data.frame(fread("~/Documents/DustMassDensity/data/literature_evo/metallicity/bellstedt20_meanZgas.csv"))

peroux_ZDTG = data.frame(fread("~/Documents/DustMassDensity/data/literature_evo/omega/aa58_peroux_suppl_tables/Peroux_SupplementalTable4_Fig7_10.asc"))

# magplot(
#   peroux_ZDTG$log_Z,
#   peroux_ZDTG$log_DTG
# )
# hf_ZDTG = hyper.fit(
#   X = cbind(peroux_ZDTG$log_Z,peroux_ZDTG$log_DTG), 
#   algo.func = "LD"
# ) 
# plot(hf_ZDTG)
# 
# z_DTG_fun = approxfun(
#   
#   lbt2z(rev(bellstedt20_meanZ$X)),
#   10^( hf_ZDTG$parm[1] * rev(log10(bellstedt20_meanZ$Y)) + hf_ZDTG$parm[2] ), 
#   rule = 2
#   
# )

rhoCrit0 = cosgrowRhoCrit(z = 0, ref = "Planck18", Dist = "Co")
rhoCritz = cosgrowRhoCrit(z = zmids, ref = "Planck18", Dist = "Co")


df = data.frame(

  "z" = zmids,
  "rhoDustQ50" = (mdust$Q50*LSS_corr/Md_corr),
  "rhoDustQ16" = (mdustwAGN$Q16*LSS_corr/Md_corr),
  "rhoDustQ84" = (mdust$Q84*LSS_corr/Md_corr),
  "rhoStarQ50" = approx(x = zvec, y = 10^dsilva25$CSMHQ50, xout = zmids)$y,
  "rhoStarQ16" = approx(x = zvec, y = 10^dsilva25$CSMHQ16, xout = zmids)$y,
  "rhoStarQ84" = approx(x = zvec, y = 10^dsilva25$CSMHQ84, xout = zmids)$y,
  "rhoBHQ50" = approx(x = zvec, y = smbhQ50, xout = zmids)$y,
  "rhoBHQ16" = approx(x = zvec, y = smbhQ16, xout = zmids)$y,
  "rhoBHQ84" = approx(x = zvec, y = smbhQ84, xout = zmids)$y,
  "rhoZgas" = approx(x = lbt2z(bellstedt20_gasZdensity$X), y = bellstedt20_gasZdensity$Y, xout = zmids, rule = 2)$y,
  "rhoCritz" = rhoCritz, 
  "rhoBaryon" = cosgrowOmegaM(z=zmids, ref = "Planck18") * OmegaBaryon/cosgrowOmegaM(z=0, ref = "Planck18") * rhoCritz,
  # "rhoBaryon" = OmegaBaryon * rhoCrit0 * (1+zmids)^3,
  "rhoMatter" = cosgrowOmegaM(z=zmids, ref = "Planck18") * rhoCritz,
  "rhoLambda" = cosgrowOmegaL(z=zmids, ref = "Planck18") * rhoCritz,
  
  "OmegaDustQ50" = (mdust$Q50*LSS_corr/Md_corr)/rhoCritz,
  "OmegaDustQ16" = (mdustwAGN$Q16*LSS_corr/Md_corr)/rhoCritz,
  "OmegaDustQ84" = (mdust$Q84*LSS_corr/Md_corr)/rhoCritz,
  "OmegaStarQ50" = approx(x = zvec, y = 10^dsilva25$CSMHQ50, xout = zmids)$y/rhoCritz,
  "OmegaStarQ16" = approx(x = zvec, y = 10^dsilva25$CSMHQ16, xout = zmids)$y/rhoCritz,
  "OmegaStarQ84" = approx(x = zvec, y = 10^dsilva25$CSMHQ84, xout = zmids)$y/rhoCritz,
  "OmegaBHQ50" = approx(x = zvec, y = smbhQ50, xout = zmids)$y/rhoCritz,
  "OmegaBHQ16" = approx(x = zvec, y = smbhQ16, xout = zmids)$y/rhoCritz,
  "OmegaBHQ84" = approx(x = zvec, y = smbhQ84, xout = zmids)$y/rhoCritz,
  "OmegaZgas" = approx(x = lbt2z(bellstedt20_gasZdensity$X), y = bellstedt20_gasZdensity$Y, xout = zmids, rule = 2)$y/rhoCritz,
  "OmegaBaryon" = cosgrowOmegaM(z=zmids, ref = "Planck18") * OmegaBaryon/cosgrowOmegaM(z=0, ref = "Planck18"),
  "OmegaMatter" = cosgrowOmegaM(z=zmids, ref = "Planck18"),
  "OmegaLambda" = cosgrowOmegaL(z=zmids, ref = "Planck18"),
  
  "Omega0DustQ50" = (mdust$Q50*LSS_corr/Md_corr)/cosgrowRhoCrit(z = 0, ref = "Planck18", Dist = "Co"),
  "Omega0DustQ16" = (mdustwAGN$Q16*LSS_corr/Md_corr)/cosgrowRhoCrit(z = 0, ref = "Planck18", Dist = "Co"),
  "Omega0DustQ84" = (mdust$Q84*LSS_corr/Md_corr)/cosgrowRhoCrit(z = 0, ref = "Planck18", Dist = "Co"),
  "Omega0StarQ50" = approx(x = zvec, y = 10^dsilva25$CSMHQ50, xout = zmids)$y/cosgrowRhoCrit(z = 0, ref = "Planck18", Dist = "Co"),
  "Omega0StarQ16" = approx(x = zvec, y = 10^dsilva25$CSMHQ16, xout = zmids)$y/cosgrowRhoCrit(z = 0, ref = "Planck18", Dist = "Co"),
  "Omega0StarQ84" = approx(x = zvec, y = 10^dsilva25$CSMHQ84, xout = zmids)$y/cosgrowRhoCrit(z = 0, ref = "Planck18", Dist = "Co"),
  "Omega0BHQ50" = approx(x = zvec, y = smbhQ50, xout = zmids)$y/cosgrowRhoCrit(z = 0, ref = "Planck18", Dist = "Co"),
  "Omega0BHQ16" = approx(x = zvec, y = smbhQ16, xout = zmids)$y/cosgrowRhoCrit(z = 0, ref = "Planck18", Dist = "Co"),
  "Omega0BHQ84" = approx(x = zvec, y = smbhQ84, xout = zmids)$y/cosgrowRhoCrit(z = 0, ref = "Planck18", Dist = "Co"),
  "Omega0Zgas" = approx(x = lbt2z(bellstedt20_gasZdensity$X), y = bellstedt20_gasZdensity$Y, xout = zmids, rule = 2)$y/cosgrowRhoCrit(z = 0, ref = "Planck18", Dist = "Co")
)

h5delete(file = h5file, name = "Omega")
h5write(file = h5file, obj = df, name = "Omega")
