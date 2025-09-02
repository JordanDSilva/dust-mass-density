## spare crap

LL_SMF_ultra =  function(p, Data){
  
  Mstar_t1 = p[1]
  Mstar_t2 = p[2]
  Mstar_t3 = p[3]
  
  alpha_t1 = p[4]
  alpha_t2 = p[5]
  alpha_t3 = p[6]
  
  beta_t1 = p[7]
  beta_t2 = p[8]
  beta_t3 = p[9]
  
  phi1_t1 = p[10]
  phi1_t2 = p[11]
  phi1_t3 = p[12]
  
  phi2_t1 = p[13]
  phi2_t2 = p[14]
  phi2_t3 = p[15]
  
  norm_ = p[16]
  
  csmh_test = sapply(Data$z, function(zz){
    M = Mstar_t1 + Mstar_t2*zz + Mstar_t3*zz^2
    alpha = alpha_t1 + alpha_t2*zz + alpha_t3*zz^2
    beta = beta_t1 + beta_t2*zz + beta_t3*zz^2
    phi1 = phi1_t1 + phi1_t2*zz + phi1_t3*zz^2
    phi2 = phi2_t1 + phi2_t2*zz + phi2_t3*zz^2
    
    regressed_Schechter = double_schechter(
      x = mdustvec,
      p = c(M, alpha, beta, phi1, phi2)
    )
    
    trapz(
      mdustvec,
      10^mdustvec * regressed_Schechter
    ) 
  })
  
  likelihood = dnorm(
    x = c(
      Data$Mstar, Data$alpha, Data$beta, Data$phi1, Data$phi2, Data$csmh
    ),
    mean = c(
      M, alpha, beta, phi1, phi2, log10(csmh_test) + norm_
    ),
    sd = c(
      Data$MstarErr, Data$alphaErr, Data$betaErr, Data$phi1Err, Data$phi2Err, Data$csmh
    ),
    log = TRUE
  )
  
  # prior = dnorm(
  #   x = c(
  #     Mstar_t1, Mstar_t2, Mstar_t3,
  #     alpha_t1, alpha_t2, alpha_t3,
  #     beta_t1, beta_t2, beta_t3,
  #     phi1_t1, phi1_t2, phi1_t3,
  #     phi2_t1, phi2_t2, phi2_t3,
  #     norm_
  #   ),
  #   mean = c(
  #     10.83, 0.15, -0.03,
  #     -0.58, 0.05, 0.02,
  #     -1.49, -0.09, 0.02,
  #     -2.31, -0.66, 0.02,
  #     -3.33, -0.16, -0.002, 
  #     0.0
  #   ),
  #   sd = c(
  #     0.1, 0.1, 0.1,
  #     0.1, 0.1, 0.1,
  #     0.1, 0.1, 0.1,
  #     0.1, 0.1, 0.1, 
  #     0.1, 0.1, 0.1,
  #     0.5
  #   ), 
  #   log = TRUE
  # )
  # 
  LL = sum(likelihood, na.rm = TRUE)
  # if(is.infinite(LL)){
  #   return(-99999999)
  # }else{
  #   return(LL)
  # }
  return(LL)
}
SMF_ultra = function(z, p){
  Mstar_t1 = p[1]
  Mstar_t2 = p[2]
  Mstar_t3 = p[3]
  
  alpha_t1 = p[4]
  alpha_t2 = p[5]
  alpha_t3 = p[6]
  
  beta_t1 = p[7]
  beta_t2 = p[8]
  beta_t3 = p[9]
  
  phi1_t1 = p[10]
  phi1_t2 = p[11]
  phi1_t3 = p[12]
  
  phi2_t1 = p[13]
  phi2_t2 = p[14]
  phi2_t3 = p[15]
  
  norm_t1 = p[16]
  
  magplot(
    NA, xlim = c(6, 12), ylim = c(1e-8, 10), 
    log = "y"
  )
  colors_z = colorRampPalette(
    colors = c("blue", "purple", "red")
  )
  colors_z = colors_z(length(z))
  
  count = 1
  csmh_test = sapply(z, function(zz){
    M = Mstar_t1 + Mstar_t2*zz + Mstar_t3*zz^2
    alpha = alpha_t1 + alpha_t2*zz + alpha_t3*zz^2
    beta = beta_t1 + beta_t2*zz + beta_t3*zz^2
    phi1 = phi1_t1 + phi1_t2*zz + phi1_t3*zz^2
    phi2 = phi2_t1 + phi2_t2*zz + phi2_t3*zz^2
    
    regressed_Schechter = double_schechter(
      x = mdustvec,
      p = c(M, alpha, beta, phi1, phi2)
    )
    lines(
      mdustvec, 
      regressed_Schechter, 
      col = colors_z[count]
    )
    count <<- count + 1
    
    trapz(
      mdustvec,
      10^mdustvec * regressed_Schechter
    )
  })
  return(csmh_test)
  
}

LL_ultra_opt = optim(
  par = c(c(11.0, 0.0, 0.0), c(-1.0, 0.0, 0.0), c(-2.0, 0.0, 0.0), c(-2.0, 0.0, 0.0), c(-3.0, 0.0, 0.0), c(0)),
  fn = LL_SMF_ultra,
  Data = list(
    Mstar = combine_hybrid_smf_par$M,
    alpha = combine_hybrid_smf_par$alpha,
    beta = combine_hybrid_smf_par$beta,
    phi1 = combine_hybrid_smf_par$phi1,
    phi2 = combine_hybrid_smf_par$phi2,
    
    MstarErr = combine_hybrid_smf_par$MErr,
    alphaErr = combine_hybrid_smf_par$alphaErr,
    betaErr = combine_hybrid_smf_par$betaErr,
    phi1Err = combine_hybrid_smf_par$phi1Err,
    phi2Err = combine_hybrid_smf_par$phi2Err,
    
    csmh = log10(csmh_hybrid$Q50),
    csmhErr = csmh_hybrid$ERR/(log(10)*csmh_hybrid$Q50),
    
    z = zmids
  ),
  control = list(fnscale = -1),
  hessian = TRUE
)

csmh_foobar = SMF_ultra(
  zmids, p = LL_ultra_opt$par
)