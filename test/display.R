inputs <- t(read.table('data/inputs'))
targets <- t(read.table('data/targets'))
means <- t(read.table('data/means'))
vars <- t(read.table('data/variances'))
stds <- sqrt(vars)
sigma <- 0.5
sigma2 <- sigma * sigma

pred_vars <- vars + sigma2
pred_stds <- sqrt(pred_vars)

plot(inputs, targets, pch=20, col='blue')

f = function (x) sin(x) / x

# asdf
points(inputs, f(inputs), pch=20, col='darkgreen')
points(inputs, f(inputs) + sigma, pch=20, col='lightgreen')
points(inputs, f(inputs) - sigma, pch=20, col='lightgreen')

points(inputs, means, pch=20, col='darkred')
points(inputs, means + sigma, pch=20, col='red')
points(inputs, means - sigma, pch=20, col='red')

points(inputs, means + pred_stds, pch=20, col='pink')
points(inputs, means - pred_stds, pch=20, col='pink')

points(inputs, means + stds, pch=20, col='yellow')
points(inputs, means - stds, pch=20, col='yellow')
