cond_HMM <- function(){
  for(j in 1:n_site){
    z[j] ~ dcat(phi[j,1:3])  
    phi[j,1] <- 1 - psi       # Prob. of non-occupation
    phi[j,2] <- psi * (1 - R) # Prob. of presence
    phi[j,3] <- psi * R       # Prob. of occupancy and presence
    s[j,1] ~ dcat(delta[])
    for(i in 2:20){
      s[j,i] ~ dcat(Gamma[z[j],s[j,i-1],])
    }
    states[j,1] ~ dcat(Gamma[z[j],s[j,20],])
    #states[j,1] ~ dcat(Gamma[z[j],z[j],])
    x[j,1]~dpois(lambda[states[j,1]])
    for(i in 2:n_visit){
      states[j,i]~dcat(Gamma[z[j],states[j,i-1],])
      x[j,i]~dpois(lambda[states[j,i]])
    }
  }
  # Priors
  psi ~ dunif(0,1)
  R ~ dunif(0,1)
  #dirichlet prior
  for(i in 1:m){
    delta[i] <- 1/m
    v[i] <- 1
  }
  for(i in 1:m){
    tau[i]~dgamma(1,0.08)
  }
  lambda[1]<-tau[1]
  for(i in 2:m){
    lambda[i]<-lambda[i-1]+tau[i]
  }
  
  Gamma[1,1,1:2] ~ ddirch(v[1:2])
  Gamma[1,1,3]  <- 0
  Gamma[1,2,1:2] ~ ddirch(v[1:2])
  Gamma[1,2,3] <- 0
  Gamma[1,3,1] <- 0
  Gamma[1,3,2] <- 0
  Gamma[1,3,3] <- 0
  Gamma[2,1,1:2] ~ ddirch(v[1:2])
  Gamma[2,1,3]  <- 0
  Gamma[2,2,1:2] ~ ddirch(v[1:2])
  Gamma[2,2,3] <- 0
  Gamma[2,3,1] <- 0
  Gamma[2,3,2] <- 0
  Gamma[2,3,3] <- 0
  Gamma[3,1,1:m] ~ ddirch(v[])
  Gamma[3,2,1:m] ~ ddirch(v[])
  Gamma[3,3,1:m] ~ ddirch(v[])
}

x = as.matrix(HMM_df[,-1])
n_visit = ncol(x)
n_site = nrow(x)
n = 20
m = 3
m_bio = 3

inits <- function(){
  list("z"=rep(1,n_site))
}
mod = jags(data = list("x", "n", "m", "m_bio", "n_site", "n_visit" ), inits = inits, 
           parameters.to.save = c("z","lambda","Gamma", "psi", "R"), 
           model.file = cond_HMM, n.burnin = 5000, n.iter = 10000, n.chains = 1)


out_summary <- mod$BUGSoutput$summary
state_rows <- grepl(rownames(out_summary), pattern = "^z\\[")
probs <- out_summary[state_rows, "50%"]

TS_adjusted <- c(3L, 0L, 0L, 0L, 3L, 0L, 0L, 2L, 0L, 3L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 4L, 4L, 4L, 4L, 0L, 0L, 1L, 0L, 2L, 0L, 0L, 
                 4L, 4L, 4L, 4L, 4L, 4L, 4L, 4L, 0L, 0L, 3L, 0L, 3L, 3L, 3L, 3L, 3L, 1L, 1L, 3L, 2L, 2L, 1L, 3L, 1L, 1L, 2L, 1L, 0L, 4L, 4L, 4L, 4L, 4L, 4L, 4L)
TS_adjusted[TS_adjusted<1] <- NA
TS_adjusted[TS_adjusted>3] <- 3

result <- confusionMatrix(factor(TS_adjusted[which(!is.na(TS_adjusted))], levels = c(1,2,3)), 
                          factor(round(probs[which(!is.na(TS_adjusted))]), levels = c(1,2,3)),
                          mode = "prec_recall")