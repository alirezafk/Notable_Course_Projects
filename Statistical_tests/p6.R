#Generating two normal samples
N <- 10000
X <- rnorm(N, mean = 0, sd = sqrt(3))
Y <- rnorm(N, mean = 0, sd = sqrt(10))

#F test
var.test(X, Y, ratio = 1, alternative = c("two.sided"),
         conf.level = 0.95)
