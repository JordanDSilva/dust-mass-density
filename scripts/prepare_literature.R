library(data.table)

## Prepare literature

## DMF
dunne11_raw = data.frame(fread("~/Documents/DustMassDensity/data/literature_evo/dmf/dunne11_raw.csv"))
dunne11_raw$z = 0.05
dunne11_raw$z[13:18] = 0.15
dunne11_raw$z[19:23] = 0.25
dunne11_raw$z[24:29] = 0.35
dunne11_raw$z[30:33] = 0.45
fwrite(
  dunne11_raw, 
  "~/Documents/DustMassDensity/data/literature_evo/dmf/dunne11.csv"
)

beeston24_raw = data.frame(fread("~/Documents/DustMassDensity/data/literature_evo/dmf/beeston24_raw.csv"))
beeston24_raw$z = 0.05
beeston24_raw$z[15:21] = 0.15
beeston24_raw$z[22:25] = 0.25
beeston24_raw$z[26:29] = 0.35
beeston24_raw$z[30:32] = 0.45
fwrite(
  beeston24_raw, 
  "~/Documents/DustMassDensity/data/literature_evo/dmf/beeston24.csv"
)

pozzi_z = c(0.175, 0.325, 0.6, 1, 1.5, 2.15)
kk = 1
pozzi_raw = lapply(
  list.files(path = "~/Documents/DustMassDensity/data/literature_evo/dmf/", pattern = glob2rx("*pozzi*raw*"), full.names = TRUE), 
  function(x){
    foo = data.frame(fread(x))
    foo$Mdust = 10^foo$Mdust
    foo$phi = 10^foo$phi
    foo$z = pozzi_z[kk]
    kk <<- kk + 1
    return(foo)
  }
)
pozzi20 = rbindlist(pozzi_raw)
fwrite(
  pozzi20, 
  "~/Documents/DustMassDensity/data/literature_evo/dmf/pozzi20.csv"
)

## CDMH
pozzi20 = data.frame(
  "z" = c(0.175, 0.325, 0.6, 1, 1.5, 2.15), 
  "cdmh" = c(1.32, 2.24, 3.62, 3.85, 2.22, 1.38) * 1e5,
  "err" = c(0.17, 0.29, 0.46, 0.48, 0.29, 0.17) * 1e5
)
fwrite(
  pozzi20, 
  "~/Documents/DustMassDensity/data/literature_evo/cdmh/pozzi20.csv"
)

beeston24 = data.frame(
  "z" = c(0.05, 0.15, 0.25, 0.35, 0.45), 
  "cdmh" = c(1.47, 2.04, 2.47, 2.97, 2.76)*1e5,
  "err" = c(0.13, 0.14, 0.04, 0.12, 0.36)*1e5
)
fwrite(
  beeston24, 
  "~/Documents/DustMassDensity/data/literature_evo/cdmh/beeston24.csv"
)

dunne11 = data.frame(
  "z" = c(0.05, 0.15, 0.25, 0.35, 0.45), 
  "cdmh" = c(0.98, 1.51, 2.08, 2.06, 2.26)*1e5,
  "err" = c(0.14, 0.16, 0.29, 0.75, 0.41)*1e5
)
fwrite(
  beeston24, 
  "~/Documents/DustMassDensity/data/literature_evo/cdmh/dunne11.csv"
)

chiangz = c(0.10, 0.22, 0.34, 0.48, 0.64, 0.81, 1.00, 1.20, 1.43, 1.68, 1.96, 2.27, 2.61, 2.98, 3.40, 3.85)
chiang25_raw = data.frame(fread("~/Documents/DustMassDensity/data/literature_evo/cdmh/chiang25_raw.csv"))
chiang25_up = chiang25_raw$cdmh[1:16]
chiang25_med = chiang25_raw$cdmh[17:32]
chiang25_down = chiang25_raw$cdmh[33:48]
plot(chiangz, chiang25_up)
points(chiangz, chiang25_med)
points(chiangz, chiang25_down)

chiang25 = data.frame(
  "z" = chiangz, 
  "cdmh" = 10^approx(x = chiang25_raw$z[17:32], y = chiang25_raw$cdmh[17:32], xout = chiangz, rule = 2)$y,
  "errhi" = (approx(x = chiang25_raw$z[1:16], y = chiang25_raw$cdmh[1:16], xout = chiangz, rule = 2)$y - approx(x = chiang25_raw$z[17:32], y = chiang25_raw$cdmh[17:32], xout = chiangz, rule = 2)$y) * 10^approx(x = chiang25_raw$z[17:32], y = chiang25_raw$cdmh[17:32], xout = chiangz, rule = 2)$y * log(10),
  "errlo" = (approx(x = chiang25_raw$z[17:32], y = chiang25_raw$cdmh[17:32], xout = chiangz, rule = 2)$y - approx(x = chiang25_raw$z[33:48], y = chiang25_raw$cdmh[33:48], xout = chiangz, rule = 2)$y) * 10^approx(x = chiang25_raw$z[17:32], y = chiang25_raw$cdmh[17:32], xout = chiangz, rule = 2)$y * log(10)
)
fwrite(
  chiang25,
  "~/Documents/DustMassDensity/data/literature_evo/cdmh/chiang25.csv"
)

eales24 <- data.frame(
  z = c(0.28, 0.82, 1.36, 1.90, 2.44, 2.98, 3.52, 4.05, 4.59, 5.13),
  zmin = c(0.0, 0.5, 1.1, 1.6, 2.2, 2.7, 3.2, 3.8, 4.3, 4.9),
  zmax = c(0.5, 1.1, 1.6, 2.2, 2.7, 3.2, 3.8, 4.3, 4.9, 5.5),
  cdmh = c(0.85, 2.27, 2.69, 2.96, 2.09, 1.68, 0.78, 0.54, 0.40, 0.67)*1e5,
  err = c(0.11, 0.13, 0.12, 0.12, 0.12, 0.11, 0.11, 0.14, 0.24, 0.41)*1e5,
  correction = c(0.10, 0.22, 0.20, 0.34, 0.32, 0.20, 0.22, 0.19, 0.15, 0.26)*1e5
)
fwrite(
  eales24,
  "~/Documents/DustMassDensity/data/literature_evo/cdmh/eales24.csv"
)

walterH2 <- data.frame(
  z = c(
    0.03, 0.185, 0.45, 0.93, 1.375, 2.56, 3.74, 0.98, 1.51, 2.51, 3.64, 4.75,
    2.40, 5.8, 2.40, 0.45, 0.8, 1.3, 1.95, 2.75, 0.35, 0.65, 0.95, 1.3, 1.75,
    2.25, 2.75, 3.5, 0.4, 0.85, 1.5
  ),
  value = c(
    0.104, 0.009, 0.11, 0.46, 0.55, 0.29, 0.24, 0.41, 0.48, 0.37, 0.15, 0.10,
    0.27, 0.047, 0.28, 0.20, 0.30, 0.51, 0.33, 0.42, 0.43, 0.47, 0.59, 0.84,
    0.90, 0.74, 0.63, 0.36, 0.23, 0.60, 0.48
  )*1e8,
  err_neg = c(
    0.009, 0.007, 0.05, 0.18, 0.15, 0.11, 0.07, 0.12, 0.11, 0.11, 0.07, 0.06,
    0.11, 0.023, 0.12, 0.14, 0.09, 0.055, 0.035, 0.045, 0.05, 0.07, 0.09, 0.13,
    0.13, 0.11, 0.09, 0.05, 0.05, 0.17, 0.07
  )*1e8,
  err_pos = c(
    0.009, 0.019, 0.10, 0.27, 0.20, 0.15, 0.09, 0.11, 0.12, 0.10, 0.07, 0.07,
    0.16, 0.034, 0.18, 0.16, 0.09, 0.055, 0.035, 0.045, 0.06, 0.08, 0.10, 0.14,
    0.16, 0.12, 0.11, 0.06, 0.06, 0.18, 0.09
  )*1e8
)

walterHI <- data.frame(
  z = c(
    0.0, 0.0, 0.0, 0.028, 0.096, 0.06, 0.1, 0.2, 0.24, 0.3, 1.265,
    1.0, 0.75, 1.27, 0.245, 0.805, 1.775, 2.35, 3.1, 4.175, 2.14, 2.97, 3.855,
    2.15, 2.45, 2.75, 3.05, 3.35, 2.3, 2.55, 2.85, 3.15, 3.35, 4.75,
    2.975, 3.615, 4.43, 4.005, 4.88
  ),
  value = c(
    0.60, 0.51, 1.03, 0.65, 0.73, 0.37, 0.53, 0.55, 1.11, 0.62, NA,
    0.99, 1.11, 0.96, 0.40, 0.33, 1.01, 1.45, 1.33, 0.82, 2.13, 0.79, 1.78,
    1.42, 1.25, 1.50, 1.58, 1.83, 0.74, 1.00, 1.00, 1.40, 1.62, 1.58,
    1.48, 1.31, 1.16, 1.87, 1.56
  )*1e8,
  err_neg = c(
    0.05, 0.045, 0.085, 0.06, 0.10, 0.055, 0.04, 0.07, 0.245, 0.045, NA,
    0.36, 0.185, 0.115, 0.19, 0.16, 0.29, 0.30, 0.30, 0.27, 0.86, 0.38, 0.42,
    0.07, 0.06, 0.07, 0.11, 0.19, 0.13, 0.11, 0.10, 0.12, 0.28, 0.30,
    0.50, 0.45, 0.42, 0.40, 0.27
  )*1e8,
  err_pos = c(
    0.05, 0.045, 0.085, 0.13, 0.13, 0.055, 0.04, 0.07, 0.245, 0.045, NA,
    0.55, 0.185, 0.115, 0.31, 0.26, 0.34, 0.34, 0.35, 0.30, 1.28, 0.63, 0.49,
    0.07, 0.06, 0.07, 0.11, 0.19, 0.13, 0.11, 0.10, 0.12, 0.28, 0.34,
    0.66, 0.51, 0.52, 0.42, 0.31
  )*1e8
)
magplot(
  walterH2$z,
  walterH2$value, 
  log = "y"
)
magerr(
  walterH2$z,
  walterH2$value,
  ylo = walterH2$err_neg,
  yhi = walterH2$err_pos
)

magplot(
  walterHI$z,
  walterHI$value, 
  log = "y", 
  ylim = c(4e6, 5e8)
)
magerr(
  walterHI$z,
  walterHI$value,
  ylo = walterHI$err_neg,
  yhi = walterHI$err_pos
)

fwrite(
  walterHI, 
  "~/Documents/DustMassDensity/data/literature_evo/cgmh/walter20HI.csv"
)
fwrite(
  walterH2, 
  "~/Documents/DustMassDensity/data/literature_evo/cgmh/walter20H2.csv"
)

## HII mass function
santoro22 = data.frame(
  "Ltot" = 10^c(40.38, 39.71, 41.24, 41.47, 40.75, 40.72 , 40.87, 40.14, 41.62, 41.75, 41.00, 40.84, 40.43, 41.40, 41.31, 41.58, 41.56, 40.69, 41.44),
  "Mstar" = 10^c(9.33, 9.49, 9.96, 9.94, 10.01, 9.87, 10.15, 10.41, 10.37, 10.71, 10.44, 10.61, 10.66, 10.73, 10.62, 10.73, 10.69, 10.85, 10.99),
  "SFR" = 10^c(0.14, 0.05, 0.70, 0.98, 0.29, 0.37, 0.37, 0.11, 1.29, 2.02, 0.42, 0.40, 0.23, 1.23, 0.80, 1.46, 0.84, 0.38, 1.29)/0.63,
  "Distance" = c(9.01, 9.84, 15.85, 18.99, 19.57, 17.22, 18.63, 18.83, 17.69, 19.4,12.22, 9.96, 11.32, 13.1, 16.99, 15.21, 15.77, 5.2, 18.72)
)
santoro22$Nions = santoro22$Ltot / 1.36e-12
santoro22$HIImass = santoro22$Nions * (1.67262192e-27 / 1.9891e+30) / (2.6e-13 * 100)
