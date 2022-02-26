
occ <- function(){
  # Likelihood
  for(s in 1:n_site) { 
    # biological model
    z[s] ~ dbern(psi)
    for(t in f[s]:l[s]) {
      x[s,t] ~ dbern(p*z[s])
    } #t
  } #s
  
  # Priors
  p ~ dunif(0, 1) # Uninformative prior
  psi ~ dunif(0, 1)
}

thresholds <- expand.grid(c(0.80, 0.85, 0.90, 0.95),c(1:5))
thresholds$accuracy <- rep(NA, nrow(thresholds))
thresholds$Kappa <- rep(NA, nrow(thresholds))

get.first <- function(x) min(which(!is.na(x)))
get.last <- function(x) max(which(!is.na(x)))

TS_adjusted <- c(3L, 0L, 0L, 0L, 3L, 0L, 0L, 2L, 0L, 3L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 4L, 4L, 4L, 4L, 0L, 0L, 1L, 0L, 2L, 0L, 0L, 
                 4L, 4L, 4L, 4L, 4L, 4L, 4L, 4L, 0L, 0L, 3L, 0L, 3L, 3L, 3L, 3L, 3L, 1L, 1L, 3L, 2L, 2L, 1L, 3L, 1L, 1L, 2L, 1L, 0L, 4L, 4L, 4L, 4L, 4L, 4L, 4L)
TS_adjusted[TS_adjusted<1] <- NA
TS_adjusted[TS_adjusted<3] <- 0
TS_adjusted[TS_adjusted>2] <- 1


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
 
  
  x[x<thresholds[i,2]] <- 0
  x[x>=thresholds[i,2]] <- 1
  # Create vector with occasion of marking
  f <- apply(x, 1, get.first)
  # Create vector with occasion of marking
  l <- apply(x, 1, get.last)
  n_visit = ncol(x)
  n_site = nrow(x)
  inits <- function(){
    list(
      "z"=rep(1,n_site)
    )
  }
  mod <- jags(data = list("x", "n_site", "n_visit", "f", "l"), inits = inits,
             parameters.to.save = c("z", "psi", "p"),
             model.file = occ, n.burnin = 5000, n.iter = 10000, n.chains = 3)
  out_summary <- mod$BUGSoutput$summary
  state_rows <- grepl(rownames(out_summary), pattern = "^z\\[")
  probs <- out_summary[state_rows, "50%"]
  
  result <- confusionMatrix(
                            factor(round(probs[which(!is.na(TS_adjusted))]), levels = c(0,1)),
                            factor(TS_adjusted[which(!is.na(TS_adjusted))], levels = c(0,1)), 
                            mode = "prec_recall")
  thresholds$accuracy[i] <- result$overall['Accuracy']
  thresholds$Kappa[i] <- result$overall['Kappa']
  print(thresholds[i,])
  
}
#write.csv(thresholds, "./data/clean/occ_threshold_assessment.csv")