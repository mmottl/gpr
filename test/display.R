inputs <- t(read.table('data/inputs'))
targets <- t(read.table('data/targets'))
inducing <- t(read.table('data/inducing_inputs'))
inducing_means <- t(read.table('data/inducing_means'))
means <- read.table('data/means')
vars <- t(read.table('data/variances'))
stds <- sqrt(vars)
sigma2 <- (read.table('data/sigma2'))[,1]
sigma <- sqrt(sigma2)
sample1 <- read.table('data/sample1')
sample2 <- read.table('data/sample2')
sample3 <- read.table('data/sample3')
fic_sample1 <- read.table('data/fic_sample1')
fic_sample2 <- read.table('data/fic_sample2')
fic_sample3 <- read.table('data/fic_sample3')

pred_vars <- vars + sigma2
pred_stds <- sqrt(pred_vars)

t = qnorm(1 - 0.05/2, 0, 1)

f = function (x) sin(3 * x) / x + (x - 3) / (x * x + 1)

mk_frame <- function(v) { data.frame(inputs, v)[order(inputs),] }

plot_means <- function (means, lty, lwd, col) {
  lines(mk_frame(means), lty=lty, lwd=lwd, col=col)
}

plot_means_stds <- function (means, mlty, stds, slty, col) {
  plot_means(means, mlty, 4, col)
  plot_std <- function (v) { lines(mk_frame(v), lty=slty, lwd=2, col=col) }
  plot_std(means + t * stds)
  plot_std(means - t * stds)
}

plot(inputs, targets, type='p', lwd=1, pch=20, col='blue')

plot_means_stds(f(inputs), 1, sigma, 3, 'green')
plot_means_stds(means, 1, stds, 2, 'red')
plot_means_stds(means, 1, pred_stds, 3, 'red')

plot_sample <- function (sample) {
  plot_means(sample, 1, 1, 'black')
}

plot_fic_sample <- function (sample) {
  plot_means(sample, 1, 1, 'violet')
}

lapply(list(sample1, sample2, sample3), plot_sample)
lapply(list(fic_sample1, fic_sample2, fic_sample3), plot_fic_sample)

points(inducing, inducing_means, lwd=10, col='orange')

legend(
  -5, 5,
  c(
    "targets",
    "true mean", "true 95% confidence interval",
    "mean", "95% confidence interval (mean)", "95% confidence interval",
    "FITC mean samples", "FIC mean samples", "Inducing inputs"
    ),
  col=c("blue", "green", "green", "red", "red", "red", "black", "violet", "orange"),
  pch=c(20,32,32,32,32,32,32,32),
  lty=c(0,1,3,1,3,3,1,1),
  lwd=c(1,4,2,4,3,2,1,1,10))
