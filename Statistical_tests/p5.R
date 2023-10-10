e <- matrix(0, 13)
n <- 6115
for(i in c(1:13)){
  e[i] <- n*dbinom(i-1, size = 12, prob = 0.5)
}

f <- c(7, 45, 181, 478, 829, 1112,1343, 1033, 670, 286, 104, 24, 3)
f[2] <- f[1]+f[2]
e[2] <- e[1]+e[2]
f[12] <- f[12]+f[13]
e[12] <- e[12]+e[13]
T_stat <- sum((e[2:12]-f[2:12])^2/e[2:12])

#Chi square dist with df=9
df <- 9
N <- 10000
smpl <- rchisq(N, df)
hist(smpl)

