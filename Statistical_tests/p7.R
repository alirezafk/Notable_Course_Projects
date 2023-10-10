f <- file.choose()

#Importing csv data
grades <- read.csv(file = f, head = FALSE)
grades <- as.matrix(grades)
X <- grades[,1]
Y <- grades[,2]

#Parametric t test
t.test(X,Y)


wilcox.test(X,Y)

