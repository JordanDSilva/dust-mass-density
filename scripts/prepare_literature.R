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
