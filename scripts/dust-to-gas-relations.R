library(magicaxis)
library(Rfits)
library(ProSpect)
library(data.table)
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
deVis19_PL = function(Z, doDTG = FALSE, a = 2.45, b=-23.30){
  #2.45 * (12+log(O/H)) - 23.30
  
  xSol = 8.69
  ZOH = log10(Z / 0.014) + xSol
  
  # a = a
  # b = -23.30
  
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
lines(
  10^seq(-4, log10(0.05), 0.1), 
  deVis19_PL(10^seq(-4, log10(0.05), 0.1), a = 2.0, b = -19.56)
)
lines(
  10^seq(-4, log10(0.05), 0.1), 
  deVis19_PL(10^seq(-4, log10(0.05), 0.1), a = 2.13, b = -20.93)
)
lines(
  10^seq(-4, log10(0.05), 0.1), 
  deVis19_PL(10^seq(-4, log10(0.05), 0.1), a = 2.15, b = -21.19)
)
lines(
  10^seq(-4, log10(0.05), 0.1), 
  deVis19_PL(10^seq(-4, log10(0.05), 0.1), a = 2.10, b = -20.91)
)
lines(
  10^seq(-4, log10(0.05), 0.1), 
  deVis19_PL(10^seq(-4, log10(0.05), 0.1), a = 1.78, b = -18.52)
)
lines(
  10^seq(-4, log10(0.05), 0.1), 
  deVis19_PL(10^seq(-4, log10(0.05), 0.1), a = 1.95, b = -19.96)
)
abline(h = draine07)
abline(v = 0.014)

magplot(
  10^seq(-4, log10(0.05), 0.01), 
  RR14_BPL(10^seq(-4, log10(0.05), 0.01))/RR14_BPL(10^seq(-4, log10(0.05), 0.01)), 
  log = "xy"
)
lines(
  10^seq(-4, log10(0.05), 0.1), 
  deVis19_PL(10^seq(-4, log10(0.05), 0.1))/RR14_BPL(10^seq(-4, log10(0.05), 0.1))
)
lines(
  10^seq(-4, log10(0.05), 0.1), 
  deVis19_PL(10^seq(-4, log10(0.05), 0.1), a = 2.0, b = -19.56)/RR14_BPL(10^seq(-4, log10(0.05), 0.1))
)
lines(
  10^seq(-4, log10(0.05), 0.1), 
  deVis19_PL(10^seq(-4, log10(0.05), 0.1), a = 2.13, b = -20.93)/RR14_BPL(10^seq(-4, log10(0.05), 0.1))
)
lines(
  10^seq(-4, log10(0.05), 0.1), 
  deVis19_PL(10^seq(-4, log10(0.05), 0.1), a = 2.15, b = -21.19)/RR14_BPL(10^seq(-4, log10(0.05), 0.1))
)
lines(
  10^seq(-4, log10(0.05), 0.1), 
  deVis19_PL(10^seq(-4, log10(0.05), 0.1), a = 2.10, b = -20.91)/RR14_BPL(10^seq(-4, log10(0.05), 0.1))
)
lines(
  10^seq(-4, log10(0.05), 0.1), 
  deVis19_PL(10^seq(-4, log10(0.05), 0.1), a = 1.78, b = -18.52)/RR14_BPL(10^seq(-4, log10(0.05), 0.1))
)
lines(
  10^seq(-4, log10(0.05), 0.1), 
  deVis19_PL(10^seq(-4, log10(0.05), 0.1), a = 1.95, b = -19.96)/RR14_BPL(10^seq(-4, log10(0.05), 0.1))
)


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

## RR14 metallicities 
df <- data.frame(
  Name = c("Haro11", "Haro2", "Haro3", "He2-10", "HS0017+1055", "HS0052+2536",
           "HS0822+3542", "HS1222+3741", "HS1236+3937", "HS1304+3529", "HS1319+3224",
           "HS1330+3651", "HS1442+4250", "HS2352+2733", "IZw18", "IC10", "IIZw40",
           "Mrk1089", "Mrk1450", "Mrk153", "Mrk209", "Mrk930", "NGC1140", "NGC1569",
           "NGC1705", "NGC2366", "NGC4214", "NGC4449", "NGC4861", "NGC5253", "NGC625",
           "NGC6822", "Pox186", "SBS0335-052", "SBS1159+545", "SBS1211+540",
           "SBS1249+493", "SBS1415+437", "SBS1533+574", "Tol0618-402", "Tol1214-277",
           "UGC4483", "UGCA20"),
  log_OH = c(8.36, 8.23, 8.28, 8.43, 7.63, 8.04, 7.32, 7.79, 7.72, 7.93, 7.81,
             7.98, 7.60, 8.40, 7.14, 8.17, 8.23, 8.10, 7.84, 7.86, 7.74, 8.03,
             8.38, 8.02, 8.27, 7.70, 8.26, 8.20, 7.89, 8.25, 8.22, 7.96, 7.70,
             7.25, 7.44, 7.58, 7.68, 7.55, 8.05, 8.09, 7.52, 7.46, 7.50)
)
df2 <- data.frame(
  Name = c("UM133", "UM311", "UM448", "UM461", "VIIZw403",
           "NGC0337", "NGC0584", "NGC0628", "NGC0855", "NGC0925",
           "NGC1097", "NGC1266", "NGC1291", "NGC1316", "NGC1377",
           "NGC1404", "IC0342", "NGC1482", "NGC1512", "NGC2146",
           "HoII", "DDO053", "NGC2798", "NGC2841", "NGC2915",
           "HoI", "NGC2976", "NGC3049", "NGC3077", "M81dwB",
           "NGC3190", "NGC3184", "NGC3198", "IC2574", "NGC3265",
           "NGC3351", "NGC3521", "NGC3621", "NGC3627", "NGC3773",
           "NGC3938", "NGC4236"),
  log_OH = c(7.82, 8.36, 8.32, 7.73, 7.66,
             8.18, 8.43, 8.35, 8.29, 8.25,
             8.47, 8.29, 8.52, 8.77, 8.29,
             8.54, 8.46, 8.11, 8.56, 8.68,
             7.72, 7.60, 8.34, 8.54, 7.94,
             7.61, 8.36, 8.53, 8.52, 7.84,
             8.49, 8.51, 8.34, 7.85, 8.27,
             8.60, 8.39, 8.27, 8.34, 8.43,
             8.42, 8.17)
)
df3 <- data.frame(
  Name = c("NGC4254", "NGC4321", "NGC4536", "NGC4559", "NGC4569",
           "NGC4579", "NGC4594", "NGC4625", "NGC4631", "NGC4725",
           "NGC4736", "DDO154", "NGC4826", "DDO165", "NGC5055",
           "NGC5398", "NGC5408", "NGC5457", "NGC5474", "NGC5713",
           "NGC5866", "NGC6946", "NGC7331", "NGC7793",
           "M83", "NGC1808", "NGC7552", "M82", "NGC1068",
           "NGC0891", "MGC+02-04-025", "NGC7469", "NGC5256",
           "NGC5953", "M51", "NGC3995", "NGC3994", "NGC6052",
           "NGC1222", "NGC7674", "NGC4670"),
  log_OH = c(8.45, 8.50, 8.21, 8.29, 8.58,
             8.54, 8.54, 8.35, 8.12, 8.35,
             8.31, 7.54, 8.54, 7.63, 8.40,
             8.35, 7.81, 8.42, 8.31, 8.24,
             8.47, 8.40, 8.34, 8.31,
             8.62, 9.10, 8.35, 8.51, 9.00,
             8.90, 8.52, 8.70, 8.53,
             8.73, 8.55, 8.66, 8.61, 8.65,
             8.57, 8.14, 8.30)
)
df_all <- rbind(df, df2, df3)
df_all$log_OH

xSol = 8.69
ZOH = log10(Z / 0.014) + xSol
maghist(
  10^(df$log_OH - xSol)
)

deVis19 = fread("~/Downloads/DP_metallicities_global.csv")

maghist(
  10^(df$log_OH - xSol) 
)
maghist(
  10^(deVis19$`12+logOH_PG16_S` - xSol)
)
0.05/0.014
