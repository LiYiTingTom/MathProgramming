# drop Na
if (all(is.na(CO2))) {
  df <- na.omit(CO2)
} else {
  df <- CO2
}

# check dim
N <- dim(df)[1]

# splite target and features
X <- as.matrix(df[1:N, 4:5])
y <- as.vector(df[1:N, 2])

# Max Min Normalization
concMax <- as.numeric(max(X[,2]))
concMin <- as.numeric(min(X[,2]))
for (i in c(1:N)) {
  X[i,2] <- (as.numeric(X[i,2])-concMin)/(concMax-concMin)
}

# encoding
y[which(y == "Quebec")] <- -1
y[which(y == "Mississippi")] <- 1
y <- as.numeric(y)
plot(X, col = y+3, pch = 20)



# solve QP prob.
P1 <- 1/2 * diag(y) %*% (X %*% t(X)) %*% diag(y)
q1 <- -rep(1, N)
s1 <- 0

C1 <- rbind(y, -y, diag(1, N), -diag(1, N))
lambda <- 0.1
d1 <- c(0,0, rep(0, N), rep(-1 / (2*N*lambda), N))

# ADMM derived update formula of x, z
p1 <- nrow(P1)
p2 <- nrow(C1)
x <- x0 <- rep(0, p1)
z <- z0 <- rep(1, p2)
u <- u0 <- rep(0, p2)

rho <- 1

M1 <- -solve(1/rho * P1 + rho * t(C1) %*% C1)

# Barrier Method
Mu <- (1/2) ^ (1:10)
Ta <- 0.001 * (1/(1:10))

for (i2 in 1:length(Mu))
{
  mu <- Mu[i2]
  # ADMM
  for (i1 in 1:1000)
  {
    v1 <- -z - d1 + u
    x <- M1 %*% (q1 / rho + t(C1) %*% v1)

    v2 <- C1 %*% x -d1 + u
    z <- (1/2) * (v2 + sqrt(v2^2 + 4 * mu / rho))

    if (sum(abs(C1 %*% x - z - d1)) < Ta[i2]) break
    u <- u + (C1 %*% x -z -d1)
  }
}

x.star <- x
#which(C1 %*% x.star - d1 >= 0)

w1 = colMeans( (c(x.star) * y) %o% rep(1, ncol(X)) * X)

# Estimation of b
b1 = mean(X %*% w1 - y)

# Binary classifier
h = function(x) {
  x = matrix(x, ncol=1)
  return(sign(w1 %*% x - b1))
}
apply(X, 1, h)

# Plot classification boundary
a1 = b1/w1[2]
b1 = -w1[1]/w1[2]
plot(X, col = y+3, pch = 20)
abline(a = a1, b = b1)
