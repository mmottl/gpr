inputs <- t(read.table('data2/inputs'))
targets <- t(read.table('data2/targets'))
means <- read.table('data2/means')
vars <- t(read.table('data2/variances'))
stds <- sqrt(vars)
sigma2 <- (read.table('data2/sigma2'))[,1]
sigma <- sqrt(sigma2)
sample1 <- read.table('data2/sample1')
sample2 <- (read.table('data2/sample2'))
sample3 <- (read.table('data2/sample3'))
sample4 <- (read.table('data2/sample4'))
sample5 <- (read.table('data2/sample5'))
inducing <- t(read.table('data2/inducing_inputs'))
inducing_means <- t(read.table('data2/inducing_means'))

pred_vars <- vars + sigma2
pred_stds <- sqrt(pred_vars)

plot(inputs, targets, type='p', lwd=1, pch=20, col='blue')

f = function (x) sin(x) / x

true = data.frame(inputs, f(inputs))[order(inputs),]
stdsu = data.frame(inputs, means + stds)[order(inputs),]
stdsl = data.frame(inputs, means - stds)[order(inputs),]
means = data.frame(inputs, means)[order(inputs),]
sample1 = data.frame(inputs, sample1)[order(inputs),]
sample2 = data.frame(inputs, sample2)[order(inputs),]
sample3 = data.frame(inputs, sample3)[order(inputs),]
sample4 = data.frame(inputs, sample4)[order(inputs),]
sample5 = data.frame(inputs, sample5)[order(inputs),]

lines(true, lty=1, lwd=4, col='red')
lines(means, lty=1, lwd=5, col='green')
lines(stdsu, lty=1, lwd=3, col='lightgreen')
lines(stdsl, lty=1, lwd=3, col='lightgreen')
lines(sample1, lty=2, lwd=2, col='black')
lines(sample2, lty=2, lwd=2, col='black')
lines(sample3, lty=2, lwd=2, col='black')
lines(sample4, lty=2, lwd=2, col='black')
lines(sample5, lty=2, lwd=2, col='black')

points(inducing, inducing_means, lwd=10, col='yellow')

#points(inputs, f(inputs2) + sigma, pch=20, col='darkgreen')
#points(inputs, f(inputs2) - sigma, pch=20, col='darkgreen')

#points(inputs, means, pch=20, col='red')
#points(inputs, means + sigma, pch=20, col='darkred')
#points(inputs, means - sigma, pch=20, col='darkred')

#points(inputs, means + pred_stds, pch=20, col='pink')
#points(inputs, means - pred_stds, pch=20, col='pink')

#points(inputs, means + stds, pch=20, col='yellow')
#points(inputs, means - stds, pch=20, col='yellow')
