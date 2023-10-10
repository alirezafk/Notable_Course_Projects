a <- 0.05
sig <- 10
n <- 1000
z <- qnorm(a/2)
D <- seq(from = 0, to = 4, by = 0.005)
N <- c(5:1000)
pow <- matrix(0, length(D))
pow_n <- matrix(0, length(N))

i <- 1
for (d in D) {
  power <- 1-pnorm(z - (d/sig)*sqrt(n/2), lower.tail = FALSE)+pnorm(-z - (d?sig)*sqrt(n/2), lower.tail = FALSE)
  pow[i] <- power
  i <- i+1
}
# Plotting power versus delta
plot(D, pow, type="l", col = 'blue')

d <- 2
i <- 1
for (n in N) {
  power <- 1-pnorm(z - (d/sig)*sqrt(n/2), lower.tail = FALSE)+pnorm(-z - (d/sig)*sqrt(n/2), ?ower.tail = FALSE)
  pow_n[i] <- power
  i <- i+1
}

# Plotting power versus delta
plot(N, pow_n, type="l", col = 'red')

z1 <- outer(N, D, function(N, D){
  power <- 1-pnorm(z - (D/sig)*sqrt(N/2), lower.tail = FALSE)+pnorm(-z - (D/sig)*sqrt(N/2), lower.ta?l = FALSE)
})

filled.contour(N, D, z1, nlevels = 10, main = "color plot", pty="s")
