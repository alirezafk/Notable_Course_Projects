n0 <- 10
size <- 1000
Dn <- matrix(0, size-n0+1)
for(n in c(n0:size)) {
  X <- rnorm(n)
  Fn <- ecdf(X)
  Y <- sort(X)
  Xn <- pnorm(Y)
  D <- abs(Xn-Fn(Y))
  Dn[n-n0] <- max(D)
}

plot(Dn, type = 'l', col = 'blue')

