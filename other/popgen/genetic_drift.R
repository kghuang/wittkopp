probA <- c()
probA[1] <- 0.3
total_pop <- 100
gen <- 1
total_gen <- 100

for (i in 1:(total_gen-1)) { 
  popA <- 0
  for (j in 1:total_pop){
    random <- runif(1)
    if (random < probA[gen]){
      popA <- popA + 1
      }
  }
  probA[gen+1] <- popA/total_pop
  gen <- gen+1
}

plot(1:total_gen, probA, xlab = "Generations", ylab = "Proportion Gene A", type = "l")