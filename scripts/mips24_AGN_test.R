library(ProSpect)
library(magicaxis)
library(data.table)

# 7 March 2026

h5file = '~/Documents/DustMassDensity/data/all_data.h5'
catalogueDir = "~/Documents/DustMassDensity/data/gama_devils_catalogues/"

devilsd10_AGN = fread("~/Documents/DustMassDensity/data/gama_devils_catalogues/AGNTotalCat_MasterCat4.csv")
mag_Y = -2.5*log10(devilsd10_AGN$flux_Y)+8.9 
mag_M24 = -2.5*log10(devilsd10_AGN$flux_M24)+8.9
maghist(
  mag_Y[mag_Y <= 21.2] - mag_M24[mag_Y <= 21.2], xlim = c(-10, 10), breaks = 30, 
  xlab = "Y-Mips24 [mag]"
)
legend(
  x = "topleft", 
  legend = c(paste0("Median = ", round(median(mag_Y[mag_Y <= 21.2] - mag_M24[mag_Y <= 21.2], na.rm = TRUE), 3)))
)

## These were matched back

devilsd10_noAGN = readRDS(paste0(catalogueDir, 'DEVILS_D10ProSpectCat_02_02_2021_v0.3.rds')) #catalogue that Jess done without AGN contribution
devilsd10_noAGN = devilsd10_noAGN$cat
devilsd10_noAGN$area = 1.47
devilsd10_AGN$area = 1.47
devilsd10_noAGN = data.frame(devilsd10_noAGN[order(devilsd10_noAGN$UID),])
devilsd10_AGN = data.frame(devilsd10_AGN[order(devilsd10_AGN$UID),])
devilsd10_AGN$FIRInput = ifelse( is.na(devilsd10_AGN$FIRInput), 0L, devilsd10_AGN$FIRInput)
devilsd10_noAGN$FIRInput = devilsd10_AGN$FIRInput

sum( mag_Y >= 21.2 & devilsd10_AGN$FIRInput == 1L )

maghist(
  mag_Y[mag_Y >= 21.2 & devilsd10_AGN$FIRInput == 1L] - mag_M24[mag_Y >= 21.2 & devilsd10_AGN$FIRInput == 1L], xlim = c(-10, 10), breaks = 30, 
  xlab = "Y-Mips24 [mag]"
)



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

magplot(
  -2.5*log10(devilsd10_AGN$flux_Y[devilsd10_AGN_idx])+8.9, 
  -2.5*log10(devilsd10_AGN$flux_Y[devilsd10_AGN_idx])+8.9 - (-2.5*log10(devilsd10_AGN$flux_M24[devilsd10_AGN_idx])+8.9), 
  xlim = c(15,70), 
  ylim = c(-5,5)
)

non_na_idx = !is.na(devilsd10_AGN$flux_Y) & !is.na(devilsd10_AGN$flux_M24)

idx = devilsd10_AGN$FIRInput == 1L


devilsd10_AGN$
magbin(
  mag_Y[idx],
  mag_Y[idx] - mag_M24[idx], 
  ylim = c(-10,10)
)
abline(h = 3, col = "red")

nc_M24 = maghist(mag_M24, breaks = seq(10,100,0.5), plot = FALSE)
magplot(
  nc_M24$mids,
  nc_M24$counts/(0.5*1.47),
  log = "y", 
  xlim = c(10,100), 
  ylim = c(1, 1e6)
)

maghist(
  mag_Y[mag_Y<21.2 & mag_M24<17] - mag_M24[mag_Y<21.2 & mag_M24<17], log = "y"
)

mdix = 
devilsd10_AGN$flux_M24[]

(sum(is.na(mag_M24)) + sum(idx))/length(mag_M24)

median(
  mag_Y[mag_Y<21.2 & idx] - mag_M24[mag_Y<21.2 & idx], na.rm = TRUE
)
