load data/inputs
load data/targets
load data/inducing_inputs
load data/sigma2

sigma = sqrt(sigma2);
[dim, N] = size(inputs);
[dim, M] = size(inducing_inputs);

function res = eval_rbf2(r2)
  res = exp(-0.5 * r2);
end

function res = k(x, y)
  [dim, n1] = size(x);
  n2 = size(y, 2);
  repmat(sum(y' .* y', 2), 1, n1);
  r2 = repmat(sum(x' .* x', 2), 1, n2) - 2 * x' * y + repmat(sum(y' .* y', 2)', n1, 1);
  res = eval_rbf2(r2);
end

km = k(inducing_inputs, inducing_inputs);
kmn = k(inducing_inputs, inputs);
kn = k(inputs, inputs);

km_chol = chol(km);

qn = kmn' * inv(km) * kmn;
lam = diag(diag(kn - qn));
lam_sigma2 = lam + sigma2 * eye(N);
inv_lam_sigma2 = inv(lam_sigma2);
inv_lam_sigma = sqrt(inv_lam_sigma2);
kmn_ = kmn * inv_lam_sigma;
y = targets;
y_ = inv_lam_sigma * y;
log_det_lam_sigma2 = log(det(lam_sigma2));
b = km + kmn * inv_lam_sigma2 * kmn';
b_chol = chol(b);
kmn_y_ = kmn_ * y_;

l1_2 = log(det(b)) - log(det(km)) + log(det(lam_sigma2))
l2_2 = y' * inv(qn + lam_sigma2) * y

neg_log_likelihood = (l1_2 + l2_2 + N * log(2*pi)) / 2

% Trained
km_chol
kmn_y_
b_chol
neg_log_likelihood

% Ed's stuff
hyp = [0; 0; log(sigma2)];
w = [reshape(inducing_inputs', M*dim, 1); hyp];
his_neg_log_likelihood = spgp_lik(w, y, inputs', M)
