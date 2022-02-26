MS <- function(){
  # Priors
  psi ~ dunif(0,1)
  R ~ dunif(0,1)
  # p2 ~ dunif(0,1)
  # p3 ~ dunif(0,1)
  delta ~ dunif(0,1)
  
  mu_p2[1] ~ dnorm(0,0.368)
  mu_p3[1] ~ dnorm(0,0.368)
  beta ~ dnorm(0,0.368)
  
  
  
  # Define state vector
  for (s in 1:n_site){
    logit(p2_2[s]) <- mu_p2[1] + beta * num_weeks[s]
    logit(p3_2[s]) <- mu_p3[1] + beta * num_weeks[s]
    for (k in 1:n_visit) {
      # Define observation matrix
      p2[s,1,k,1] <- 1
      p2[s,1,k,2] <- 0
      p2[s,1,k,3] <- 0
      p2[s,2,k,1] <- 1 - p2_2[s]
      p2[s,2,k,2] <- p2_2[s]
      p2[s,2,k,3] <- 0
      p2[s,3,k,1] <- 1 - p3_2[s]
      p2[s,3,k,2] <- p3_2[s] * (1 - delta)
      p2[s,3,k,3] <- p3_2[s] * delta
    } #k
    
    z[s] ~ dcat(phi[s,1:3])
    #phi[site,TRUE]
    phi[s,1] <- 1 - psi            # Prob. of non-occupation
    phi[s,2] <- psi * (1 - R)        # Prob. of singlet occupancy
    phi[s,3] <- psi * R  # Prob. of occupancy and repro
    
    for (k in f[s]:l[s]) {
      x[s, k] ~ dcat( p2[s, z[s], k, 1:3] ) # encounter data
    } #k
  } #s
  
  # derived parameters
  logit(dp[1]) <- mu_p2[1]
  logit(dp[2]) <- mu_p3[1]
}


thresholds <- expand.grid(c(0.80, 0.85, 0.90, 0.95),c(1:3),c(5:15))
thresholds$accuracy <- rep(NA, nrow(thresholds))
thresholds$Kappa <- rep(NA, nrow(thresholds))

get.first <- function(x) min(which(x!=0))
get.last <- function(x) max(which(x!=0))
TS_adjusted <- c(3L, 0L, 0L, 0L, 3L, 0L, 0L, 2L, 0L, 3L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 4L, 4L, 4L, 4L, 0L, 0L, 1L, 0L, 2L, 0L, 0L, 
                 4L, 4L, 4L, 4L, 4L, 4L, 4L, 4L, 0L, 0L, 3L, 0L, 3L, 3L, 3L, 3L, 3L, 1L, 1L, 3L, 2L, 2L, 1L, 3L, 1L, 1L, 2L, 1L, 0L, 4L, 4L, 4L, 4L, 4L, 4L, 4L)
TS_adjusted[TS_adjusted<1] <- NA
TS_adjusted[TS_adjusted>3] <- 3


for(i in 1:nrow(thresholds)){
  print(paste0("Iteration: ",i))
  tbl_plot <- all_result %>%
    dplyr::select(site, date, Offset, BRMA) %>% 
    mutate(week = week(date)
    ) %>% 
    group_by(site,week) %>% 
    summarise(n = sum(BRMA>=thresholds[i,1]))
  HMM_df <- tbl_plot %>% 
    drop_na(week) %>% 
    mutate(n_week = as.numeric(week)) %>% 
    pivot_wider(site,names_from = n_week, values_from = n, names_sort = TRUE) #%>% 
  x <- as.matrix(HMM_df[,-1])
  x[x<thresholds[i,2]] <- 1
  x[x>=thresholds[i,2] & x<thresholds[i,3]] <- 2
  x[x>=thresholds[i,3]] <- 3
  print(table(apply(x, 1, max, na.rm=TRUE)))
  sum_weeks <- function(i){
    i[i>1] <- 1
    sum(i, na.rm=TRUE)
  }
  num_weeks <- apply(x, 1, sum_weeks)
  # Create vector with occasion of marking
  f <- apply(x, 1, get.first)
  # Create vector with occasion of marking
  l <- apply(x, 1, get.last)
  n_visit = ncol(x)
  n_site = nrow(x)
  inits <- function(){
    list(
      "z"=rep(3,n_site)
    )
  }
  mod <- jags(data = list("x", "n_site", "n_visit", "f", "l", "num_weeks" ), inits = inits,
             parameters.to.save = c("z", "psi", "R", "dp"),#, "beta"
             model.file = MS, n.burnin = 5000, n.iter = 10000, n.chains = 3)
  out_summary <- mod$BUGSoutput$summary
  state_rows <- grepl(rownames(out_summary), pattern = "^z\\[")
  probs <- out_summary[state_rows, "50%"]
  result <- confusionMatrix( 
                            factor(round(probs[which(!is.na(TS_adjusted))]), levels = c(1,2,3)),
                            factor(TS_adjusted[which(!is.na(TS_adjusted))], levels = c(1,2,3)),
                            mode = "prec_recall")
  thresholds$accuracy[i] <- result$overall['Accuracy']
  thresholds$Kappa[i] <- result$overall['Kappa']
  print(thresholds[i,])
  
}
#write.csv(thresholds, "./data/clean/MS_threshold_assessment.csv")
