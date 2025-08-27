library(Rfits)
library(celestial)
library(magicaxis)
library(scales)
library(rhdf5)
library(matrixStats)
library(data.table)

magphys_cat = Rfits_read_table("~/Documents/DustMassDensity/data/GAMA/MagPhysv06.fits")
prospect_cat = Rfits_read_table("~/Documents/DustMassDensity/data/GAMA/ProSpectv03.fits")

lambdar_phot = Rfits_read_table("~/Documents/DustMassDensity/data/GAMA/LambdarCatv01.fits")
profoundgkv_phot_raw = Rfits_read_table("~/Documents/DustMassDensity/data/GAMA/gkvScienceCatv02.fits")
profoundIR_phot_raw = Rfits_read_table("~/Documents/DustMassDensity/data/GAMA/gkvFarIRv03.fits")
profoundIR_phot = profoundIR_phot_raw[profoundIR_phot_raw$uberID %in% profoundgkv_phot_raw$uberID, ]
profoundgkv_phot = profoundgkv_phot_raw[profoundgkv_phot_raw$uberID %in% profoundIR_phot$uberID, ]
message(all(profoundIR_phot$uberID == profoundgkv_phot$uberID)) ## Check that they are matched
# profoundIR_phot$Z = profoundgkv_phot$Z
# profoundIR_phot$RA = profoundgkv_phot$RAmax
# profoundIR_phot$DEC = profoundgkv_phot$Decmax
profound_phot = cbind(profoundgkv_phot, profoundIR_phot)

lambdar_idx = lambdar_phot$P100_flux/lambdar_phot$P100_fluxerr > 5 & 
  lambdar_phot$P160_flux/lambdar_phot$P160_fluxerr > 5 & 
  lambdar_phot$S250_flux/lambdar_phot$S250_fluxerr > 5 &
  lambdar_phot$S350_flux/lambdar_phot$S350_fluxerr > 5 &
  lambdar_phot$S500_flux/lambdar_phot$S500_fluxerr > 5

profound_idx = profound_phot$flux_PSF_p100/profound_phot$flux_PSF_err_p100 > 5 & 
  profound_phot$flux_PSF_p160/profound_phot$flux_PSF_err_p160 > 5 & 
  profound_phot$flux_PSF_s250/profound_phot$flux_PSF_err_s250 > 5 &
  profound_phot$flux_PSF_s350/profound_phot$flux_PSF_err_s350 > 5 &
  profound_phot$flux_PSF_s500/profound_phot$flux_PSF_err_s500 > 5

lambdar_trim = lambdar_phot[lambdar_idx,]
profound_trim = profound_phot[profound_idx,]

idx_coordmatch = coordmatch(
  coordref = profound_trim[, c("RAmax", "Decmax")],
  coordcompare = lambdar_trim[, c("RA", "DEC")]
)

profound_match = profound_trim[idx_coordmatch$bestmatch$refID, ]
lambdar_match = lambdar_trim[idx_coordmatch$bestmatch$compareID, ]

magphys_match = magphys_cat[magphys_cat$CATAID %in% lambdar_match$CATAID, ]
lambdar_magphys_match = lambdar_match[lambdar_match$CATAID %in% magphys_match$CATAID, ]
magphys_match$RA =  lambdar_magphys_match$RA
magphys_match$DEC =  lambdar_magphys_match$DEC

prospect_match = prospect_cat[prospect_cat$uberID %in% profound_match$uberID, ]
profound_prospect_match = profound_match[profound_match$uberID %in% prospect_cat$uberID, ]
prospect_match$RA = profound_prospect_match$RAmax
prospect_match$DEC = profound_prospect_match$Decmax
prospect_match$Z = profound_prospect_match$Z

idx2_coordmatch = coordmatch(
  coordref = prospect_match[, c("RA", "DEC")],
  coordcompare = magphys_match[, c("RA", "DEC")]
)

prospect_final_sample = prospect_match[idx2_coordmatch$bestmatch$refID, ]
magphys_final_sample = magphys_match[idx2_coordmatch$bestmatch$compareID, ]

profound_final_sample = data.frame(profound_prospect_match[idx2_coordmatch$bestmatch$refID, ])
lambdar_final_sample = data.frame(lambdar_magphys_match[idx2_coordmatch$bestmatch$compareID, ])

# magplot(
#   prospect_final_sample$RA, prospect_final_sample$DEC, pch = 1, lwd=3, col = "red", cex = 2, type = "p", xlab = "RA", ylab = "DEC"
# )
# points(
#   magphys_final_sample$RA, magphys_final_sample$DEC, pch = 4, col = "darkorange", cex = 2, lwd = 3
# )
# legend(
#   x = "bottomleft",
#   legend = c("ProFound", "Lambdar"),
#   pch = c(1, 4), 
#   cex = 1.3,
#   lwd = NA,
#   col = c("red", "darkorange")
# )
# 
# magplot(
#   (prospect_final_sample$RA - magphys_final_sample$RA), 
#   prospect_final_sample$DEC - magphys_final_sample$DEC
# )
# 
# magplot(prospect_final_sample$Z, magphys_final_sample$Z)
# 
# par(mfrow = c(2,1), mar = c(2.5, 2.5, 0.5, 0.5), oma = rep(1.0,4))
# magplot(
#   log10(prospect_final_sample$DustMass_50),
#   magphys_final_sample$mass_dust_percentile50, 
#   xlab = "ProSpect Mdust [Msun]",
#   ylab = "Magphys Mdust [Msun]"
# )
# magerr(
#   log10(prospect_final_sample$DustMass_50),
#   magphys_final_sample$mass_dust_percentile50,
#   xlo = (prospect_final_sample$DustMass_50-prospect_final_sample$DustMass_16)/(log(10) * prospect_final_sample$DustMass_50),
#   xhi = (prospect_final_sample$DustMass_84-prospect_final_sample$DustMass_50)/(log(10) * prospect_final_sample$DustMass_50),
#   ylo = magphys_final_sample$mass_dust_percentile50 - magphys_final_sample$mass_dust_percentile16,
#   yhi = magphys_final_sample$mass_dust_percentile84 - magphys_final_sample$mass_dust_percentile50
# )
# abline(0,1, col = "red", lwd = 3)
# 
# maghist(
#   log10(prospect_final_sample$DustMass_50) - magphys_final_sample$mass_dust_percentile50, breaks = seq(-3,3,0.1), 
#   ylab = "Frequency", xlab = "ProSpect - Magphys [dex]"
# )
# abline(v = 0.48889, col = "red", lwd = 3)
# legend(
#   x = "topleft", 
#   col = c("red"), 
#   lwd = 3,
#   legend = paste0("Median = ", 0.48889)
# )
# 
# par(mfrow = c(2,1), mar = c(2.5, 2.5, 0.5, 0.5), oma = rep(1.0,4))
# magplot(
#   prospect_final_sample$Z,
#   magphys_final_sample$Z, 
#   xlab = "ProSpect redshift",
#   ylab = "Magphys redshift"
# )
# abline(0,1, col = "red", lwd = 3)
# maghist(
#   prospect_final_sample$Z - magphys_final_sample$Z, 
#   xlab = "ProSpect redshift - Magphys redshift",
#   ylab = "Frequency", 
#   breaks = 40
# )
# abline(v = -0.0004900, col = "red", lwd = 3)
# legend(
#   x = "topleft", 
#   col = c("red"), 
#   lwd = 3,
#   legend = paste0("Median = ", -0.0004900)
# )

df = data.frame(
  "zProSpect" = prospect_final_sample$Z,
  "MdustProSpect" = prospect_final_sample$DustMass_50,
  "MdustErrProSpect" = 0.5*(prospect_final_sample$DustMass_84 - prospect_final_sample$DustMass_16),
  "RAProSpect" = prospect_final_sample$RA,
  "DECProSpect" = prospect_final_sample$DEC,
  
  "zMagphys" = magphys_final_sample$Z,
  "MdustMagphys" = 10^magphys_final_sample$mass_dust_percentile50,
  "MdustErrMagphys" = sqrt( (10^magphys_final_sample$mass_dust_percentile50 * log(10) * 0.5 * (magphys_final_sample$mass_dust_percentile84 - magphys_final_sample$mass_dust_percentile16))^2 ),
  "RAMagphys" = magphys_final_sample$RA,
  "DECMagphys" = magphys_final_sample$DEC,
  data.frame(magphys_final_sample)[,grep("^(mass_dust|L_dust|f_mu_IR|xi|T).*best_fit", names(magphys_cat), value = TRUE, perl = TRUE)]
)

wavelength = c(
  1539, 2316, 3528, 4760, 6326, 7599, 6654, 10229, 12556, 16499, 21571, 
  3.4e4, 4.65e4, 12.8e4, 22.4e4, 98.9e4, 156e4, 249e4, 350e4, 504e4
)

profound_fluxes = profound_final_sample[, grep("^flux_(?!err)(?!.*_err)(?!.*l$).*", names(profound_final_sample), perl = TRUE, value = TRUE)]
profound_flux_errs = profound_final_sample[, grep("err", names(profound_final_sample), perl = TRUE, value = TRUE)]
profound_fluxes[profound_fluxes < 0] = 0
profound_flux_errs = sqrt( profound_flux_errs^2 + (0.1 * profound_fluxes)^2 )

lambdar_fluxes = lambdar_final_sample[, grep("^(?!z_flux$).*flux(?!err)", names(lambdar_final_sample), perl = TRUE, value = TRUE)]
lambdar_flux_errs = lambdar_final_sample[, grep("^(?!z_fluxerr$).*err", names(lambdar_final_sample), perl = TRUE, value = TRUE)]
lambdar_fluxes[lambdar_fluxes == -999] = NA
lambdar_fluxes[lambdar_fluxes < 0] = 0
lambdar_flux_errs = sqrt( lambdar_flux_errs^2 + (0.1 * lambdar_fluxes)^2 )

profound_weighted_stack = colSums(as.matrix(profound_fluxes) * as.matrix(profound_flux_errs^-2), na.rm = TRUE) / colSums(as.matrix(profound_flux_errs^-2), na.rm = TRUE)
lambdar_weighted_stack = colSums(as.matrix(lambdar_fluxes) * as.matrix(lambdar_flux_errs^-2), na.rm = TRUE) / colSums(as.matrix(lambdar_flux_errs^-2), na.rm = TRUE)
profound_weighted_err = sqrt(1/colSums(as.matrix(profound_flux_errs^-2), na.rm = TRUE))
lambdar_weighted_err = sqrt(1/colSums(as.matrix(lambdar_flux_errs^-2), na.rm = TRUE))

df_spec = data.frame(
  "wavelength" = wavelength,
  "ProFoundStack" = profound_weighted_stack,
  "ProFoundErr" = profound_weighted_err,
  "LambdarStack" = lambdar_weighted_stack,
  "LambdarErr" = lambdar_weighted_err
)

# magplot(
#   NA,
#   log = "xy", 
#   ylim = c(1e-5, 0.5),
#   xlim = c(1e3, 1e7),
#   pch = 16
# )
# for(kk in 1:dim(profound_fluxes)[1]){
#   points(wavelength, unlist(profound_fluxes[kk,]), col = alpha("darkorange",0.1), pch = 16)
#   magerr(wavelength, unlist(profound_fluxes[kk,]), ylo = unlist(unlist(profound_flux_errs[kk,])), col = alpha("darkorange",0.1))
#   points(wavelength, unlist(lambdar_fluxes[kk,]), col = alpha("darkred",0.1), pch = 16)
#   magerr(wavelength, unlist(lambdar_fluxes[kk,]), ylo = unlist(unlist(lambdar_flux_errs[kk,])), col = alpha("darkred",0.1))
# }
# points(
#   wavelength, profound_weighted_stack, pch = 16, col = "darkorange"
# )
# magerr(
#   wavelength, profound_weighted_stack, ylo = profound_weighted_err, col = "darkorange"
# )
# points(
#   wavelength, lambdar_weighted_stack, pch =16, col = "darkred"
# )
# magerr(
#   wavelength, lambdar_weighted_stack, ylo = lambdar_weighted_err, col = "darkred"
# )

h5file = '~/Documents/DustMassDensity/data/all_data.h5'
h5delete(h5file, "Photometry")
h5createGroup(h5file, "Photometry")
h5write(
  obj = df, file = h5file, name = "Photometry/FinalSample"
)
h5write(
  obj = df_spec, file = h5file, name = "Photometry/WeightedSpec"
)

foobar = h5read(
 file = h5file, name = "Photometry/FinalSample"
)

## Save data from prospect fitting
# filtout = list(
#   approxfun(getfilt("FUV"),rule = 2),
#   approxfun(getfilt("NUV"),rule = 2),
#   approxfun(getfilt("u_VST"),rule = 2),
#   approxfun(getfilt("g_VST"),rule = 2),
#   approxfun(getfilt("r_VST"),rule = 2),
#   approxfun(getfilt("i_VST"),rule = 2),
#   approxfun(getfilt("Z_VISTA"),rule = 2),
#   approxfun(getfilt("Y_VISTA"),rule = 2),
#   approxfun(getfilt("J_VISTA"),rule = 2),
#   approxfun(getfilt("H_VISTA"),rule = 2),
#   approxfun(getfilt("K_VISTA"),rule = 2),
#   approxfun(getfilt("W1_WISE"),rule = 2),
#   approxfun(getfilt("W2_WISE"),rule = 2),
#   approxfun(getfilt("W3_WISE"),rule = 2),
#   approxfun(getfilt("W4_WISE"),rule = 2),
#   approxfun(getfilt("P100_Herschel"),rule = 2),
#   approxfun(getfilt("P160_Herschel"),rule = 2),
#   approxfun(getfilt("S250_Herschel"),rule = 2),
#   approxfun(getfilt("S350_Herschel"),rule = 2),
#   approxfun(getfilt("S500_Herschel"),rule = 2)
# )
# profound_Data = list(
#   "flux" = profound_fluxes,
#   "flux_err" = profound_flux_errs,
#   "redshifts" = prospect_match$Z,
#   "magphys_redshifts" = magphys_match$Z,
#   "uberID" = profound_final_sample$uberID,
#   "RA" = profound_final_sample$RAmax,
#   "Dec" = profound_final_sample$Decmax,
#   "filtout" = filtout,
#   "cenwave" = wavelength
# )
# saveRDS(profound_Data, "~/Documents/DustMassDensity/data/save_highSN_sample.rds")
# foo = readRDS("~/Documents/DustMassDensity/data/save_highSN_sample.rds")
# list2env(foo, envir = .GlobalEnv)
# rm(foo)
# rm(cenwave)
# save.image("~/Documents/DustMassDensity/scripts/Pawsey/save_highSN_sample.rda")
# file.remove("~/Documents/DustMassDensity/data/save_highSN_sample.rds")

prospect_fnames = list.files("/Users/22252335/Documents/DustMassDensity/Pawsey/magphysLike/", full.names = TRUE)
prospect_files = lapply(prospect_fnames, function(x){
  foo = readRDS(x)
  # plot(foo$bestfit$SEDout)
  bar = c(colQuantiles(as.matrix(foo$parm_sample), probs = c(0.16, 0.5, 0.84)), foo$RA, foo$Dec, colQuantiles(as.matrix(foo$highlander$LD_last$Posterior1[, c("alpha_SF_screen", "alpha_SF_birth")]), probs = c(0.16, 0.5, 0.84)))
  names(bar) = c(paste0(names(foo$parm_sample),"Q16"), paste0(names(foo$parm_sample),"Q50"), paste0(names(foo$parm_sample),"Q84"), "RA", "Dec", paste0(c("alpha_screen", "alpha_birth"),"Q16"), paste0(c("alpha_screen", "alpha_birth"),"Q50"), paste0(c("alpha_screen", "alpha_birth"),"Q84"))
  return(bar)
})
prospect_refit = data.frame(do.call(rbind, prospect_files))

prospect_coord_match = coordmatch(
  coordref = profound_final_sample[, c("RAmax", "Decmax")], 
  coordcompare = prospect_refit[, c("RA", "Dec")]
)
prospect_refit_match = prospect_refit[prospect_coord_match$bestmatch$compareID, ]
new_dust_mass = prospect_refit_match$DustLumBirthQ50/Dale_vM2L_func(alpha_SF = prospect_refit_match$alpha_birthQ50, Zfinal = 0.02) + prospect_refit_match$DustLumScreenQ50/Dale_vM2L_func(alpha_SF = prospect_refit_match$alpha_screenQ50, Zfinal = 0.02)

## ITS JUST CHEWY THANKS FOR THE FIFTY GIFTED 
magplot(
  log10(prospect_refit_match$DustMassQ50),
  magphys_final_sample$mass_dust_percentile50 , 
  xlim = c(2, 10),
  yli = c(2,10)
)
points(
  log10(prospect_final_sample$DustMass_50), 
  magphys_final_sample$mass_dust_percentile50, 
  pch = 16, col = "blue"
)
abline(0,1)

maghist(
  magphys_final_sample$mass_dust_percentile50 - log10(prospect_final_sample$DustMass_50)
)
maghist(
  magphys_final_sample$mass_dust_percentile50 - log10(prospect_refit_match$DustMassQ50)
)
maghist(
  magphys_final_sample$mass_dust_percentile50 - log10(new_dust_mass) 
)

