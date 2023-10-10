# example 1
N <- 70
Y1 <- rnorm(N, 2, 3)
qqnorm(Y1)
qqline(Y1, col = "red")

# example 2
Y2 <- rnorm(N, 10, 4)
qqnorm(Y2)
qqline(Y2, col = "red")