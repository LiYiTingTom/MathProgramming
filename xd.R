# Section 2. Precursors(Dua1 ascend method and method of multipliers) #

# Convex linear constrained optimization problems
# min f(x)
#subject to Ax = b                               # Linear equa1ity constraints
# f(x): objective function; Ax = b: constaint

# Example
# min (x1^2 + 2*x2^2)
# x1 + x2 = 1

f = function(x){
 x1 = x[1]
 x2 = x[2]
 out= (1*x1^2 + 2*x2^2)
 return(out)
}

c1 = function(x){
 A = matrix(c(1,1), 1, 2)
 b = 1
 return(A %*% x - b)
}

# Graph of the problem
x2 <- x1 <- seq(-2, 2, 0.1)
x12 <- expand.grid(x1, x2)
Z1 <- matrix(apply(x12, 1, f), length(x1), length(x2))
Z2 <- matrix(apply(x12, 1, c1), length(x1), length(x2))
#
image(x1, x2, Z1, main = "z = f(x1, x2)")
contour(x1, x2, Z1, levels = c(0.1, 0.2, 0.3), nlevels = 20, add = TRUE)
contour(x1, x2, Z2, levels = 0, lwd = 2, col = 3, add = TRUE)
abline(a = 1, #a = intercept, b = -1, #b = slope, lwd = 2, col = 3)
abline(h = seq(-2, 2, 0.5), v = seq(-2, 2, 0.5), lty = 3)

#===========================================================================  3/6 Example
# min (x1^2 - 1*x1*x2 + x2^2)
# 2 * x1 + x2 = 1.5
f1 = function(x){
 x1 = x[1]
 x2 = x[2]
 out = (x1^2 - x1*x2 + x2^2)
 return(out)
}

c2 = function(x){
 A = matrix(c(2,1), 1, 2)
 b = 1.5
 return(A %*% x - b)
}

# Graph of the problem
x2 <- x1 <- seq(-3, 3, 0.1)
x12 <- expand.grid(x1, x2)
Z1 <- matrix(apply(x12, 1, f1), length(x1), length(x2))
Z2 <- matrix(apply(x12, 1, c2), length(x1), length(x2))
image(x1, x2, Z1, main = "z = f1(x1, x2)")
contour(x1, x2, Z1, nlevels = 20, add = TRUE)
abline(a = 1.5, b = -2, lwd = 2, col = 3)
abline(h = seq(-3, 3, 0.5), v = seq(-3, 3, 0.5), lty = 3)




#################################
### 1. Lagrange dua1ity method: # 
#################################
# 1. Formulate the Lagreangian of prima1 optima1ity problem:
#  L(x, y) = f(x) + y^T(Ax - b)
 L <- function(x, y){
 z <- c(f(x) + matrix(y, nrow = 1) %*% as.matrix(c1(x)))
 return(z)
}
# x <- rep(0, 2); y <- rep(0, 1); L(x, y)

# 2. Solve numerica1ly the Lagrange dua1 function:
# g(y) = inf_x L(x, y)
g <- function(y){
 L.x <- function(x){
  return(L(x, y))
 }
 a1 <- optim(rep(0, 2), L.x)
 return(a1$val)
}
# g(0)
yy <- matrix(seq(-3, 3, 0.1), ncol = 1)
zz <- apply(yy, 1, g)
plot(yy, zz, type = "l",
     main = "Lagrang dua1 function of y")
# See the Lagrang dua1 function of y in p7

### 2. Define a constraint function g(x) such that the constraint is given by g(x) = 0.
# 
# 
# 
g.1 = function(x) {
	A1 = matrix(c(1,1), 1, 2)
	b1 = 1
	return(A1 %*% x - b1)
}
g.1(c(1, 2))

Z2.1 <- matrix(apply(x12, 1, g.1), length(x1), length(x2))
image(x1, x2, Z1, main = "z = f(x1, x2)")
contour(x1, x2, Z1, nlevels = 10, add = TRUE)
contour(x1, x2, Z2.1, levels = 0, add = TRUE)

### 3.

g.2 = function(x) {
	A2 = matrix(c(1, -1), 1, 2)
	b2 = -0.5
	return(A2 %*% x - b2)
}
g.2(c(1, 2))

Z2.2 <- matrix(apply(x12, 1, g.2), length(x1), length(x2))
image(x1, x2, Z1, main = "z = f(x1, x2)")
contour(x1, x2, Z1, nlevels = 10, add = TRUE)
contour(x1, x2, Z2.1, levels = 0, col = "red", lwd = 2, add = TRUE)
contour(x1, x2, Z2.2, levels = 0, col = "blue", lwd = 2,add = TRUE)

### 4.# 
A = matrix(c(1, 1, 1, -1), 2, 2, byrow = T)
b = matrix(c(1, -0.5), 2, 1)
solve(A, b)


#===============================================================1. Lagrange dua1ity (LD) method

### 1. Formulate the L(x, y) of the prima1 optima1ity problem

# L(x, y) = f(x) + y^T %*% (Ax -b)
L <- function(x, y){
 z <- f(x) + matrix(y, nrow = 1) %*% (A %*% x -b)
  return (z)
}

# Single :
A <- matrix(c(1, 1), 1, 2)
b <- 1
x <- rep(0, 2)
y <- rep(0, 1)
L(x, y)

### 2. Solve numerica1ly the LD(x, y)
# LD(y) = inf_x L(x, y)
LD = function(y) {
 L.x = function(x) {
  return (L(x, y))
 }
 a1 = optim(rep(0, 2), L.x)
 return (a1$val)
}

LD(0)
yy = matrix(seq(-3, 3, 0.1), ncol = 1)
zz = apply(yy, 1, LD)
plot(yy, zz, type = "l", main = "Lagrang dua1 function of y")

### 3. Solve the dual problem: y* = argmax_y g(y)
a1 = optim(rep(0, 1), LD, control = list(fnsca1e = -1))
y.star = a1$par

### 4. Recover x* =argmin_x L(x, y*) = the solution of the problem
L.y.star = function(x){
 return (L(x, y.star))
}
a1 = optim(rep(0, 2), L.y.star)
x.star = a1$par
x.star
#
A %*% x.star - b
image(x1, x2, Z1, main = "z = f(x1, x2)")
contour(x1, x2, Z1, add = TRUE)
abline(a = b/A[, 2], b = -A[, 1]/A[, 2], lwd = 2, col = 3)
points(x.star[1], x.star[2], col = 4, pch = 20, cex = 2)

###===================================================================== Exercise
# Find the min value of 
# f(x, y, z) = x^2 + y^2 + z^2
# ont he line that is the intersection of the plane
# x + y = 2
# and the plane
# y + z = 1

### 1.
f = function(x) {
  x1 = c(1, 2, 3) * x
  return (t(x1) %*% x1)
}

g = function(x) {
  A = matrix(c(1, 1, 0, 0, 1, 1), 2, 3, byrow = TRUE)
  b = c(2, 1)
  return (A %*% x - b)
}

L = function(x, y) {
  return (f(x) + y %*% g(x))
}

### 2.
LD = function(y) {
  L.x = function(x) {
    return (L(x, y))
  }
  a1 = optim(rep(0, 3), L.x)
  return (a1$val)
}

### 3.
a1 = optim(rep(0, 2), LD, control = list(fnscale = -1))
y.star = a1$par
round(y.star, 4)

### 4.
L.y.star = function(x) {
  return (L(x, y.star))
}
a1 = optim(rep(0, 3), L.y.star)
x.star = a1$par
x.star
#
g(x.star)

