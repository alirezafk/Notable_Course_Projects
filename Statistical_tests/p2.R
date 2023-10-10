n0 <- 10
size <- 1000
Dn <- matrix(0, size-n0+1)
Dn_kol <- matrix(0, size-n0+1)
for(n in c(n0:size)) {
  X <- rnorm(n)
  Fn <- ecdf(X)
  Y <- sort(X)
  Xn <- pnorm(Y)
  D <- abs(Y-qnorm(Fn(Y)))
  Dn[n-n0] <- max(D)
  D_kol <- abs(Xn-Fn(Y))
  Dn_kol[n-n0] <- max(D_kol)
}

plot(Dn_kol)