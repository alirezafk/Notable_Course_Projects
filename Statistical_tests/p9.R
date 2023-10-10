#sample from beta dist
size <- 50
a <- 2
b <- 5
m0 <- 0.6
M <- matrix(0, m2-m1+1)
alpha <- 0.05

for (m0 in c(m1:m2)) {
  X <- rbeta(size, a, b)
  m <- sum(X > m0)
  p_value <- pbinom(m, size, 0.5, lower.tail = FALSE)
  
  # Calculating p-value using wilcox.test
  wilcox.test(X, mu = m0, alternative = 'greater')
  
  # Calculating power of the test
  p <- 1 - pbeta(m0, a, b)
  decision <- qbinom(1 - alpha, size, 0.5)
  power <- pbinom(decision, size, p, lower.tail = FALSE)
}
