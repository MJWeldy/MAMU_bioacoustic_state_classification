HMM <- function(){
  for(i in 1:m){
    delta[i] <- 1/m
    v[i] <- 1
  }
  s[1] ~ dcat(delta[])
  for(i in 2:n_visit){
    s[i] ~ dcat(Gamma[s[i-1],])
  }
  
  for(j in 1:n_site){
    
    states[j,1] ~ dcat(Gamma[s[n_visit],])
    x[j,1]~dpois(lambda[states[j,1]])
    for(i in 2:n_visit){
      states[j,i]~dcat(Gamma[states[j,i-1],])
      x[j,i]~dpois(lambda[states[j,i]])
    }
  }
  for(i in 1:m){
    tau[i]~dgamma(1,0.08)
    Gamma[i,1:m]~ddirch(v[])
  }
  lambda[1]<-tau[1]
  for(i in 2:m){
    lambda[i]<-lambda[i-1]+tau[i]
  }
}

x <- as.matrix(HMM_df[,-1])
n_visit = ncol(x)
n_site = nrow(x)
n = 20
m = 4


mod = jags(data = list("x", "n", "m", "n_site", "n_visit" ), inits = NULL, 
           parameters.to.save = c("lambda","Gamma","states"), 
           model.file = HMM, n.iter = 10000, n.chains = 3)

out_summary <- mod$BUGSoutput$summary
state_rows <- grepl(rownames(out_summary), pattern = "^states\\[")
probs <- out_summary[state_rows, "50%"]

decoded_states <- as.data.frame(
  cbind(HMM_df$site, matrix(probs, nrow= n_site, byrow=FALSE))
)
names(decoded_states) <- c("site", paste0(13:37))

long_decoded_df <- decoded_states %>%
  pivot_longer(!site, names_to = "week", values_to = "state")

long_HMM_df <- HMM_df %>%
  pivot_longer(!site, names_to = "week", values_to = "n")

decoded_state_plot <- left_join(long_HMM_df,long_decoded_df,by = c("site", "week"))

lambda <- grepl(rownames(out_summary), pattern = "^lambda\\[")
lambdas <- out_summary[lambda, "mean"]
decoded_state_plot %>% 
  ggplot(aes(x = week, y = n, group = site, col = as.factor(state)))+
  geom_point()+
  facet_wrap(~site)+
  geom_hline(yintercept=lambdas[1],alpha = 0.5)+
  geom_hline(yintercept=lambdas[2],alpha = 0.5)+
  geom_hline(yintercept=lambdas[3],alpha = 0.5)+
  geom_hline(yintercept=lambdas[4],alpha = 0.5)+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90))


TS_adjusted <- c(3L, 0L, 0L, 0L, 3L, 0L, 0L, 2L, 0L, 3L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 4L, 4L, 4L, 4L, 0L, 0L, 1L, 0L, 2L, 0L, 0L, 
                 4L, 4L, 4L, 4L, 4L, 4L, 4L, 4L, 0L, 0L, 3L, 0L, 3L, 3L, 3L, 3L, 3L, 1L, 1L, 3L, 2L, 2L, 1L, 3L, 1L, 1L, 2L, 1L, 0L, 4L, 4L, 4L, 4L, 4L, 4L, 4L)
TS_adjusted[TS_adjusted<1] <- NA
TS_adjusted[TS_adjusted>3] <- 3

max_state <- as.numeric(apply(decoded_states[,-1], 1, max, na.rm=TRUE))
max_state[max_state>3] <- 3
confusionMatrix(
                factor(max_state[!is.na(TS_adjusted)] , levels = c(1,2,3)),
                factor(TS_adjusted[!is.na(TS_adjusted)], levels = c(1,2,3)),
                mode = "prec_recall")


