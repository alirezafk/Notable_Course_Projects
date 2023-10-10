f <- file.choose()

#Importing csv data
heights <- read.csv(file = f)
heights <- as.matrix(heights)
childs <- heights[,3]

size <- 1e4
u <- matrix(0, 1, size)
l <- matrix(0, 1, size)
a <- 1e-1/2
m <- mean(childs)
smpl_size <- 10
d_f <- smpl_size - 1
output <- 0

#Without replacement
for (i in c(1:size)) {
  smpl <- sample(childs, smpl_size, replace = FALSE, prob = NULL)
  test <- t.test(smpl, alternative = "two.sided", conf.level = 0.9)
  confidence <- test$conf.int
  l[i] <- confidence[1]
  if (u[i]>m && m>l[i]) {
    output <- output + 1
  }
}

require(plotrix)
num <- 20
plotCI(c(1:num), y = matrix(m, 1, num), ui = u[1:num], li = l[1:num], main = "Without replacement")
