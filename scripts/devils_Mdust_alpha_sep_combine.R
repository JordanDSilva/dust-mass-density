library(ProSpect)
library(magicaxis)
library(data.table)
library(Rfits)

# 7 March 2026

h5file = '~/Documents/DustMassDensity/data/all_data.h5'
catalogueDir = "~/Documents/DustMassDensity/data/gama_devils_catalogues/"

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

Mdust_new = devilsd10_hybrid$dustlum.birth/Dale_M2L_variableDTH_func(devilsd10_hybrid$alpha_SF_birth) + devilsd10_hybrid$dustlum.screen/Dale_M2L_variableDTH_func(devilsd10_hybrid$alpha_SF_screen)

Dale_L2M = approxfun(
  x = Dale_M2L_func(seq(0,4,0.01)),
  y = seq(0,4,0.01), 
  rule = 2
)

devils_alpha = Dale_L2M(devilsd10_hybrid$dustlum.total/devilsd10_hybrid$dustmass.total)
Mdust_new2 = devilsd10_hybrid$dustlum.total/Dale_M2L_variableDTH_func(devils_alpha)

print( median(devilsd10_hybrid$dustmass.total/Mdust_new2) )

maghist(
  log10(devilsd10_hybrid$dustmass.total/Mdust_new2), 
  xlab = "Separate alpha Mdust - Combined alpha Mdust", 
  main = "DEVILS galaxies"
)


profoundgkv_phot_raw = Rfits_read_table("~/Documents/DustMassDensity/data/GAMA/gkvScienceCatv02.fits")
profoundIR_phot_raw = Rfits_read_table("~/Documents/DustMassDensity/data/GAMA/gkvFarIRv03.fits")
profoundIR_phot = profoundIR_phot_raw[profoundIR_phot_raw$uberID %in% profoundgkv_phot_raw$uberID, ]
profoundgkv_phot = profoundgkv_phot_raw[profoundgkv_phot_raw$uberID %in% profoundIR_phot$uberID, ]
message(all(profoundIR_phot$uberID == profoundgkv_phot$uberID)) ## Check that they are matched
profound_phot = cbind(profoundgkv_phot, profoundIR_phot)

gama_Y = -2.5*log10(profound_phot$flux_Yt)+8.9
gama_W4 = -2.5*log10(profound_phot$flux_PSF_w4)+8.9
magplot(
  gama_Y,
  gama_W4, 
  pch =".", 
  xlim = c(5, 30),
  ylim = c(5, 30)
)
maghist(
  gama_Y[gama_Y<= 21.2] - gama_W4[gama_Y<= 21.2], 
  xlim = c(-10, 10)
)

maghist(
  log10(devilsd10_AGN$StellarMass/devilsd10_noAGN$StellarMass), log = "y"
)
