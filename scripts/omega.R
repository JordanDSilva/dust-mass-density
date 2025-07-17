library(celestial)
library(rhdf5)
library(magicaxis)
library(foreach)
library(data.table)

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

dsilva25$CAGNHQ50 + log10( 3.154e7 * 1e-7 * (0.1 * 3e8**2)^-1 / 2e30)

magplot(zvec, rhoH * (1+zvec)^3, log= "y", ylim = c(1e3, 1e15))
lines(zvec, 10^dsilva25$CSMHQ50, col = "blue")
lines(zvec, smbhQ50)

h5file = '~/Documents/DustMassDensity/data/all_data.h5'
zmids = h5read(file = h5file, name = "zmids")
mdust = h5read(file = h5file, name = "cosmic/Mdust")
Md_corr = as.numeric(h5read(file = h5file, name = "Md_corr"))

bellstedt20_gasZdensity = data.frame(fread("~/Documents/DustMassDensity/data/literature_evo/metallicity/bellstedt20_gazZdensity.csv"))

df = data.frame(
  
  "z" = zmids,
  "OmegaDustQ50" = (mdust$Q50/Md_corr)/cosgrowRhoCrit(z = 0, ref = "Planck18", Dist = "Co"),
  "OmegaDustQ16" = (mdust$Q16/Md_corr)/cosgrowRhoCrit(z = 0, ref = "Planck18", Dist = "Co"),
  "OmegaDustQ84" = (mdust$Q84/Md_corr)/cosgrowRhoCrit(z = 0, ref = "Planck18", Dist = "Co"), 
  
  "OmegaStarQ50" = approx(x = zvec, y = 10^dsilva25$CSMHQ50, xout = zmids)$y/cosgrowRhoCrit(z = 0, ref = "Planck18", Dist = "Co"),
  "OmegaStarQ16" = approx(x = zvec, y = 10^dsilva25$CSMHQ16, xout = zmids)$y/cosgrowRhoCrit(z = 0, ref = "Planck18", Dist = "Co"),
  "OmegaStarQ84" = approx(x = zvec, y = 10^dsilva25$CSMHQ84, xout = zmids)$y/cosgrowRhoCrit(z = 0, ref = "Planck18", Dist = "Co"),
  
  "OmegaBHQ50" = approx(x = zvec, y = smbhQ50, xout = zmids)$y/cosgrowRhoCrit(z = 0, ref = "Planck18", Dist = "Co"),
  "OmegaBHQ16" = approx(x = zvec, y = smbhQ16, xout = zmids)$y/cosgrowRhoCrit(z = 0, ref = "Planck18", Dist = "Co"),
  "OmegaBHQ84" = approx(x = zvec, y = smbhQ84, xout = zmids)$y/cosgrowRhoCrit(z = 0, ref = "Planck18", Dist = "Co"), 
  
  "OmegaZgas" = approx(x = lbt2z(bellstedt20_gasZdensity$X), y = bellstedt20_gasZdensity$Y, xout = zmids, rule = 2)$y/cosgrowRhoCrit(z = 0, ref = "Planck18", Dist = "Co")
)
h5delete(file = h5file, name = "Omega")
h5write(file = h5file, obj = df, name = "Omega")


# bellstedt20_meanZgas = data.frame(fread("~/Documents/DustMassDensity/data/literature_evo/metallicity/bellstedt20_meanZgas.csv"))
# 
# bellstedt20_meanOAbund = log10(bellstedt20_meanZgas$Y/0.02 * 10^(8.69-12)) + 12
# bellstedt20_meanOAbund = seq(7.0, 8.7, 0.1)
# park24_Z_DGR = -1.1381*(bellstedt20_meanOAbund)^2 + 19.6996*(bellstedt20_meanOAbund) - 85.4377
# magplot(bellstedt20_meanOAbund, park24_Z_DGR, type = "l", xlim = c(7, 8.6), ylim = c(-6, -1.5))
# 
# magplot(bellstedt20_meanZgas$X, bellstedt20_meanZgas$Y/0.0134)
