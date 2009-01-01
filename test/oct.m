load data/inputs
load data/targets
load data/inducing_inputs
load data/sigma2

sigma = sqrt(sigma2);
[dim, N] = size(inputs);
[dim, M] = size(inducing_inputs);

function res = eval_rbf2(r2)
  res = exp(-0.5 - 1.0 * r2);
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

jitter = 10e-9;
km = km + jitter*eye(M);
kn = kn + jitter*eye(N);

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

log_det_b = log(det(b));
log_det_km = log(det(km));
log_det_lam_sigma2 = log(det(lam_sigma2));

model_nll = (log_det_b - log_det_km + log_det_lam_sigma2 + N * log(2*pi)) / 2;
trained_nll = (y' * inv(qn + lam_sigma2) * y) / 2;

nll = (model_nll + trained_nll);
evidence = - nll;

model_nll_dsigma2 = trace(inv(qn + lam_sigma2)) / 2;
model_evidence_dsigma2 = - model_nll_dsigma2;

% Trained
evidence

% Ed's stuff
hyp = [-log(0.5); -1 / 2; log(sigma2)];
w = [reshape(inducing_inputs', M*dim, 1); hyp];
[eds_neg_log_likelihood, dfw] = spgp_lik(w, y, inputs', M);
eds_evidence = -eds_neg_log_likelihood
eds_dsigma2 = -dfw(end) / sigma2
