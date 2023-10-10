f <- file.choose()

#Importing csv data
heights <- read.csv(file = f)
heights <- as.matrix(heights)
fathers <- heights[,2]
n1 = 5
n2 = 500
pow <- matrix(0, n2-n1+1)

smpl_size <- 10
pop_mean <- mean(fathers)
pop_sd <- sd(fathers)

smpl <- sample(fathers, smpl_size, replace = FALSE, prob = NULL)

d_f <- smpl_size - 1

test <- t.test(smpl, mu = 60, alternative = "two.sided")

confidence <- test$conf.int
u <- confidence[2]
l <- confidence[1]

# Exact Power using normal estimation
power1 <- 0
u1 <- (u-pop_mean)/pop_sd
l1 <- (l-pop_mean)/pop_sd
power1 <- pt(u1, df = d_f, lower.tail = FALSE)
power1 <- power1 + pt(l1, df = d_f)

# Power using normal estimation
power <- 0
power <- pnorm(u, pop_mean, pop_sd, lower.tail = FALSE)
power <- power + pnorm(l, pop_mean, pop_sd)
#pow[n-n1+1] <- power1

#plot(c(n1:n2), pow)
