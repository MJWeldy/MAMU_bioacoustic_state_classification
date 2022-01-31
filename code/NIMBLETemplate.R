# Nimble Template
library(nimble)

Nimble_Code <- nimble <- nimbleCode({

})

n_params <- c("beta0")
n_constants <- list(
  n.ind = nrow(neos),
  n.years = 3
)
n_data <- list(
  y = y,
  on = c(0, 0, 0, 0, 0)
)
n_inits <- list(
  beta0 = runif(1),
  beta1 = runif(1)
)

Nimble_Model <- nimbleModel(
  code = Nimble_Code,
  constants = n.constants,
  data = n.data,
  inits = n.inits
)
Nimble_Model$initializeInfo() # Check for NA values we need to initialize
# Nimble_Model$simulate("beta.year")
# Nimble_Model$calculate("phi")

MCMC_Model <- configureMCMC(Nimble_Model, monitors = n.params, print = T, enableWAIC = T)
Model_MCMC <- buildMCMC(MCMC_Model)
Comp_Model <- compileNimble(Model_MCMC, project = Nimble_Model)
Comp_Model$run(niter = 1000, nburnin = 500) # Testing MCMC
tail(as.mcmc(as.matrix(Comp_Model$mvSamples)), n = 5)

#Longer Run
niter <- 1000
set.seed(1)
Model_samples <- runMCMC(Comp_Model, niter = niter)
# WAIC
Comp_Model$calculateWAIC()

#Change beta switches
Nimble_Model$on <- c(1,0,0,0,1)